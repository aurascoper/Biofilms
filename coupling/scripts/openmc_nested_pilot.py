#!/usr/bin/env python
"""Phase 2B: a small nested OpenMC pilot that buys a budget before spending one.

    micromamba run -n openmc-biofilms python coupling/scripts/openmc_nested_pilot.py \
        --snapshot snap.h5 --config config/reference_synthetic_biofilm_e2e.toml \
        --outdir artifacts/pilot

IT ADAPTS RATHER THAN EXECUTING A MATRIX. Scenarios run in a declared order:
zero-effect first, because a false positive there invalidates everything after
it; then the clearly detectable case, which shows the machinery spans the
decision; and only then the uncertainty-dominated case, whose expected verdict
cannot honestly be declared in advance because it depends on the spread this
pilot is measuring. It stops on budget exhaustion and returns an indeterminate
verdict rather than a cheap one.

THE GATE STATISTIC IS A DEBIASED SQUARED EFFECT, estimated across ordered
replicate pairs r != s:

    S_hat = (1 / R(R-1)) sum_{r != s} d_r^T W d_s,     d_r = D1_r - D0_r

Because replicates are independent, E[d_r^T W d_s] = E[d_r]^T W E[d_s], so this
estimates ||E[d]||^2_W with NO noise term. It is normalised by the same
cross-replicate construction on the baseline — E[||D0_r||^2_W] carries the
baseline's own noise power, so a naive squared-norm denominator is biased high —
and compared against delta^2. Squared throughout: taking a root near zero
reintroduces exactly the bias this removes. S_hat MAY BE NEGATIVE, and must be,
about half the time when nothing changed.

WHY, AND WHAT IT REPLACED. The former primary metric was the mass-weighted
relative L2 norm E = sqrt(sum_b m_b (D1-D0)^2) / sqrt(sum_b m_b D0^2). A norm is
UNSIGNED, so with zero true effect and residual Monte Carlo noise it returns
|noise|, not 0 — and returns it *stably*, so replicate scatter never reveals the
bias. Measured here, that bias was 0.0545 at the finest mesh, larger than four of
the six material levers of interest. It is still computed and still reported,
because it remains the right evidence for the N^-1/2 history scaling and for what
common random numbers buy; it is no longer what any verdict rests on.

Omega_b requires biological occupancy, positive mass, and a declared dose floor —
three PHYSICAL criteria, declared before the run. It deliberately does not
require acceptable transport uncertainty: a statistical criterion on a region
definition made the region move with the history count and made the mesh factors
cover different fractions of the domain. See `omega_b`.

The `noise_floor` scenario runs identical material with DECORRELATED seeds. Under
the raw norm it reports that norm's bias; under the debiased statistic it should
be indistinguishable from zero, and a run where it is not has a defect the gate
cannot see.

EVERY VARIANCE SOURCE IS MEASURED SEPARATELY, because the A0 lesson is that
Monte Carlo noise falls with histories while resolution loss is a floor that
histories cannot remove. A single combined number cannot say what to buy.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "calibration"))

from physical_contract import git_provenance

from biofilm_calibration.joint_uncertainty import CorrelationGroup, sample_joint

from biofilm_openmc.config import BIOFILM_CYLINDER, load_dose_rate_config
from biofilm_openmc.dose import dose_from_per_source, per_source_from_statepoint
from biofilm_openmc.feedback_uq import (debiased_squared_effect,
                                        relative_l2_effect)
from biofilm_openmc.fingerprint import transport_state_hash
from biofilm_openmc.lineage import aggregate_by_label
from biofilm_openmc.materials import (mesh_material_masses_kg,
                                      mesh_material_volumes, voxel_mass_kg)
from biofilm_openmc.mesh import resolve_mesh_dimension, upsample_field
from biofilm_openmc.snapshot import load_snapshot
from biofilm_openmc.synthetic_gate import (S0Verdict, ThresholdPolicy,
                                           VarianceBudget, decide)

# THE GATE DECIDES ON THE DEBIASED SQUARED EFFECT, not on the raw norm. The
# same 0.10 and 0.02 magnitudes are declared, squared, because the statistic is
# squared — comparing a threshold denominated in E against draws denominated in
# E^2 is a silent factor-of-ten error, which is why `metric_id` is carried and
# `decide()` refuses a mismatch.
POLICY = ThresholdPolicy(metric_id="debiased_relative_l2_squared",
                         effect_threshold=0.10 ** 2,
                         practical_importance_floor=0.02 ** 2)
# The raw unsigned norm remains REPORTED, as the diagnostic it is: it is what
# shows the N^-1/2 history scaling and what common random numbers buy.
RAW_METRIC_ID = "relative_l2"
# Declared BEFORE the run: a bin below this carries no usable dose information
# and would make a relative metric a measure of the floor.
DOSE_FLOOR_FRACTION = 1e-3
# DIAGNOSTIC CEILING ONLY. This is deliberately no longer an Omega_b membership
# criterion — see `omega_b` for why a statistical cut on a region definition
# made the mesh factors incomparable.
MAX_REL_ERR = 0.25

# The density lever moves rho rather than the composition, so it cannot be
# expressed as a feedback element map like every other scenario.
DENSITY_SCALE = {"lever_density_x1.35": 1.35}


def _nuclear_data_id() -> str:
    """Identify the cross-section library ACTUALLY IN USE, not a literal.

    `synthetic_e2e.py` hardcodes the string "endfb-viii.0". That is a claim
    about the environment made by the source file rather than read from it, so
    it stays correct only as long as nobody points OPENMC_CROSS_SECTIONS
    somewhere else. Reading the env var reports what the run really used, and
    the recorded SHA-256 (written alongside the library per docs/openmc_stack.md)
    distinguishes two libraries that share a filename.
    """
    import os
    xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
    if not xs:
        return "OPENMC_CROSS_SECTIONS unset"
    path = Path(xs)
    library = path.parent
    # The recorded checksum is the one for the downloaded ARCHIVE, not for
    # cross_sections.xml — `docs/openmc_stack.md` promises the latter and only
    # the former exists. The archive digest is the stronger identity anyway
    # (it pins every data file, not just the index), so look for it beside the
    # extracted directory before giving up.
    candidates = [path.with_name(path.name + ".sha256")]
    if library.parent.exists():
        candidates += sorted(library.parent.glob(f"{library.name}*.sha256"))
    for candidate in candidates:
        if candidate.exists():
            digest = candidate.read_text().split()[0]
            return (f"{library.name}/{path.name} "
                    f"sha256={digest[:16]} ({candidate.name})")
    return f"{library.name}/{path.name} sha256=unrecorded"


def _seed(draw: int, rep: int, state: str, paired: bool) -> int:
    """Common random numbers across the two material states — except where the
    scenario exists to measure the noise floor.

    Matched seeds are legitimate variance reduction: the two states transport
    through nearly the same geometry, so shared histories cancel common noise.
    But the old `zero_effect` control combined matched seeds with an UNCHANGED
    material, so the two models were bit-identical and the metric was forced to
    ~1e-16 by construction — at any energy, for any physics, however broken the
    transport. That is a determinism check, not a false-positive control.

    It also left the metric's noise floor unmeasured, which matters more than
    it looks: `relative_l2_effect` is an UNSIGNED norm, so residual Monte Carlo
    noise enters as positive bias that never cancels, and it does so *stably*,
    meaning replicate spread cannot reveal it either.

    `paired=False` decorrelates the two states on identical material, so the
    effect that scenario reports IS the noise floor.
    """
    return 1000 * (draw + 1) + rep + (0 if paired or state == "baseline"
                                      else 500_000)


# A coarse bin counts as biological when at least this much of its VOLUME is
# biomass. Declared, because the alternative is an arbitrary representative.
MIN_BIOMASS_VOLUME_FRACTION = 0.5


def omega_b(dose0, mass, biomass_fraction) -> np.ndarray:
    """Biological occupancy AND positive mass AND above the dose floor.

    OCCUPANCY COMES FROM THE EXACT CSG VOLUMES, not from subsampling the lattice.
    Taking `occ[::factor]` picks one representative voxel per coarse bin, so a
    bin that is 12% biomass and 88% medium counts as fully biological — and the
    two mesh factors then measure DIFFERENT regions, which shows up as
    discretisation error that is really a definitional artifact. Asking what
    fraction of the bin's volume is actually biomass makes the factors
    comparable, and the volumes are already computed for the mass denominator.

    MEMBERSHIP MUST NOT DEPEND ON THE ESTIMATOR. This previously also required
    `rel_err < MAX_REL_ERR`, which is a STATISTICAL criterion, and it bit the
    two mesh factors unequally: at factor 1 it removed bins (983/949/968 of
    8000 = 12.1%, median rel_err 0.1837 against a 0.25 cut), while at factor 2
    it never bound (a constant 215/1000 = 21.5%, median 0.0837). The factors
    were therefore compared over regions differing ~2x in fractional coverage,
    and part of the 22% factor-1-to-2 movement previously attributed to
    discretisation was that. Worse, the region would have SHIFTED on buying
    more histories, so the metric's own domain moved with its precision.

    `MAX_REL_ERR` survives as a reported diagnostic, not as a gate: the
    rel_err distribution over Omega_b is recorded per record so a
    poorly-resolved run is visible rather than silently trimmed.
    """
    floor = DOSE_FLOOR_FRACTION * float(np.nanmax(dose0)) if np.isfinite(
        np.nanmax(dose0)) else 0.0
    return ((biomass_fraction >= MIN_BIOMASS_VOLUME_FRACTION)
            & (mass > 0) & (dose0 > floor))


def biomass_volume_fraction(volumes: dict, dimension) -> np.ndarray:
    """Per-bin biomass share of the total material volume, logical (x,y,z)."""
    dim = tuple(int(d) for d in dimension)
    total = np.zeros(int(np.prod(dim)), dtype=float)
    bio = np.zeros_like(total)
    for name, per_element in volumes.items():
        arr = np.asarray(per_element, dtype=float)
        total += arr
        if name == "baseline_biomass":
            bio += arr
    frac = np.divide(bio, total, out=np.zeros_like(bio), where=total > 0)
    return frac.reshape(dim[::-1]).transpose(2, 1, 0)


def run_state(model, mesh, dim, outdir, name, mass):
    import openmc

    cwd = outdir / name
    cwd.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    sp_path = model.run(cwd=str(cwd), output=False)
    wall = time.perf_counter() - t0
    with openmc.StatePoint(sp_path) as sp:
        per_source = per_source_from_statepoint(sp, mass)
    size = Path(sp_path).stat().st_size
    return per_source, wall, size


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--outer-draws", type=int, default=4)
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--particles", type=int, default=5000,
                    help="reduced for the pilot; the point is the variance and "
                         "the cost, not a converged answer")
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--volume-samples", type=int, default=500_000)
    ap.add_argument("--mesh-factors", default="1,2,4,5",
                    help="all must divide the lattice N evenly. More than two "
                         "because the numerics bound over exactly two factors "
                         "is a single difference, which cannot show whether "
                         "the field has begun to converge")
    ap.add_argument("--levers", action="store_true",
                    help="re-measure the material lever ranking as real "
                         "scenarios that write sample rows. The ranking used "
                         "to be six hardcoded constants in the budget artifact "
                         "with no run behind them")
    ap.add_argument("--budget-seconds", type=float, default=1800.0)
    ap.add_argument("--seed", type=int, default=20260816)
    ap.add_argument("--include-near-threshold", action="store_true",
                    help="add the 20%% Gd case, whose effect sits just above "
                         "the declared threshold; its verdict is a finding")
    args = ap.parse_args(argv)

    from biofilm_openmc.model import build_biofilm_cylinder_model

    args.outdir.mkdir(parents=True, exist_ok=True)
    base = load_dose_rate_config(args.config.with_suffix(".dosimetry.toml"),
                                 kind=BIOFILM_CYLINDER)
    snapshot = load_snapshot(args.snapshot)
    n = snapshot.cell_id.shape[0]
    factors = [int(f) for f in args.mesh_factors.split(",")]

    # Outer draws: correlated calibration inputs, sampled jointly and ONCE, so
    # the same physical realisation is used for both material states.
    joint = sample_joint(
        [CorrelationGroup("gravimetry_coupon",
                          {"density_g_cm3": base.materials["baseline_biomass"].density_g_cm3},
                          shared_relative_sigma=0.02,
                          independent_relative_sigma={"density_g_cm3": 0.005})],
        {}, seed=args.seed, draws=args.outer_draws)
    densities = joint.flat()["density_g_cm3"]

    # THE FEEDBACK LEVER IS COMPOSITION, NOT DENSITY, because dose is energy
    # PER UNIT MASS: raising density raises the heating and the mass
    # denominator together, so the two nearly cancel and only second-order
    # self-shielding survives. Composition is not automatically stronger
    # either — at Cs-137 energies Compton dominates and mu/rho tracks Z/A,
    # which is ~0.5 for almost every element. What moves the dose is a
    # DEPARTURE from Z/A = 0.5: hydrogen (Z/A = 1.0), and high-Z material
    # where the photoelectric term still contributes.
    #
    # A ranked lever table used to be asserted here as six hardcoded constants
    # and copied verbatim into the budget artifact. Those numbers were
    # transcribed from exploratory runs that wrote no sample rows, no seeds and
    # no replicate counts, and one entry (Gd 20% at 0.1241) disagreed with this
    # pilot's own registry (0.1384) for what was nominally the same scenario. A
    # number with no run behind it does not belong in a provenance artifact, so
    # they are gone; `--levers` re-measures them as real scenarios that write
    # real rows.
    #
    # Gd 40% is a deliberately EXTREME synthetic material with no biological
    # pretension: its only job is to prove the gate can return a pass. 20% Gd
    # sits just above the threshold and is therefore the honest candidate for
    # the near-boundary case, whose verdict is the pilot's finding rather than
    # a declaration.
    #
    # THE THIRD FIELD IS SEED PAIRING, and it is why `noise_floor` replaced the
    # old `zero_effect`. See the `--levers` help and `_seed` below.
    scenarios = [
        ("noise_floor", None, False),
        ("clearly_detectable", {"H": 0.111894, "O": 0.488106, "Gd": 0.40}, True),
    ]
    if args.levers:
        # The lever ranking, re-measured. Density is a `None` composition with
        # a scaled density, applied at the state branch below.
        #
        # EACH LEVER RUNS BOTH PAIRED AND UNPAIRED, because the levers are
        # exactly the effects the noise floor calls into question. The floor
        # measured on IDENTICAL material is an upper bound on the metric's
        # noise: it is what survives when nothing cancels. Common random
        # numbers cancel some of it, and how much depends on how far the
        # material change decorrelates the histories — which is small for a
        # small lever and large for Gd 40%. Running each lever both ways
        # measures that retention at the lever's OWN effect size instead of
        # assuming it, so a lever below the unpaired floor can still be shown
        # resolved if pairing demonstrably removes the noise.
        for name, elements in (
                ("lever_density_x1.35", None),
                ("lever_Fe_5pct", {"H": 0.111894, "O": 0.838106, "Fe": 0.05}),
                ("lever_Gd_5pct", {"H": 0.111894, "O": 0.838106, "Gd": 0.05}),
                ("lever_dehydration_H_0.06", {"H": 0.06, "O": 0.94})):
            scenarios.append((name, elements, True))
            scenarios.append((name + "_unpaired", elements, False))
    if args.include_near_threshold:
        # RUNS THIRD AND ONLY ON REQUEST, because its verdict depends on the
        # spread this pilot measures and therefore cannot be declared in
        # advance. Adding it before the bracketing cases were known correct
        # would have made a surprising verdict uninterpretable.
        scenarios.append(
            ("near_threshold", {"H": 0.111894, "O": 0.688106, "Gd": 0.20}, True))
    started = time.perf_counter()
    runs = 0
    total_histories = 0
    statepoint_bytes = 0
    records: list[dict] = []
    debiased: list[dict] = []
    stopped_early = None

    for factor in factors:
        cfg_f = replace(base, mesh_coarsening_factor=factor,
                        batches=args.batches, particles=args.particles)
        dim = resolve_mesh_dimension((n, n, n), factor)
        # Raytrace the CSG volumes ONCE per geometry: the counterfactual varies
        # density, not geometry, so re-raytracing would inject sampling noise
        # into a comparison whose whole point is a small difference.
        probe = build_biofilm_cylinder_model(snapshot, cfg_f)
        mesh = probe.tallies[0].filters[0].mesh
        t0 = time.perf_counter()
        volumes = mesh_material_volumes(mesh, probe, n_samples=args.volume_samples)
        raytrace_wall = time.perf_counter() - t0
        bio_fraction = biomass_volume_fraction(volumes, dim)

        for scenario, feedback_elements, paired_seeds in scenarios:
            for draw, rho in enumerate(densities):
                if time.perf_counter() - started > args.budget_seconds:
                    stopped_early = ("budget exhausted before "
                                     f"{scenario}/draw{draw}/factor{factor}")
                    break
                per_state = {}
                for state in ("baseline", "feedback"):
                    mats = dict(cfg_f.materials)
                    bm = mats["baseline_biomass"]
                    elements = bm.elements
                    rho_state = float(rho)
                    if state == "feedback":
                        if feedback_elements is not None:
                            elements = tuple(sorted(feedback_elements.items()))
                        rho_state *= DENSITY_SCALE.get(
                            scenario.removesuffix("_unpaired"), 1.0)
                    mats["baseline_biomass"] = type(bm)(
                        bm.name, rho_state, elements)
                    cfg_s = replace(cfg_f, materials=mats)
                    model = build_biofilm_cylinder_model(snapshot, cfg_s)
                    m = mesh_material_masses_kg(mesh, model, dim,
                                                volumes=volumes)
                    fields = []
                    for rep in range(args.replicates):
                        cfg_r = replace(
                            cfg_s, seed=_seed(draw, rep, state, paired_seeds))
                        mdl = build_biofilm_cylinder_model(snapshot, cfg_r)
                        ps, wall, size = run_state(
                            mdl, mesh, dim, args.outdir,
                            f"f{factor}_{scenario}_d{draw}_{state}_r{rep}", m)
                        fields.append(ps)
                        runs += 1
                        total_histories += cfg_r.batches * cfg_r.particles
                        statepoint_bytes += size
                    per_state[state] = (fields, m)

                (b_fields, mass), (f_fields, _) = (per_state["baseline"],
                                                   per_state["feedback"])
                mask = omega_b(b_fields[0].field, mass, bio_fraction)

                # THE GATE STATISTIC, one per outer draw. The U-statistic runs
                # over ORDERED REPLICATE PAIRS r != s, so it consumes the
                # replicate dimension: there is no per-replicate debiased value
                # to record, and its spread comes from the jackknife rather
                # than from replicate scatter.
                deb = debiased_squared_effect([b.field for b in b_fields],
                                              [f.field for f in f_fields],
                                              mass, mask)
                debiased.append({
                    "scenario": scenario, "mesh_factor": factor,
                    "outer_draw": draw, "density_g_cm3": float(rho),
                    "paired_seeds": bool(paired_seeds),
                    "omega_b_bins": int(mask.sum()), **deb.as_dict(),
                })

                for i, (b, f) in enumerate(zip(b_fields, f_fields)):
                    err = b.rel_err[mask]
                    records.append({
                        "scenario": scenario, "mesh_factor": factor,
                        "outer_draw": draw, "replicate": i,
                        "density_g_cm3": float(rho),
                        "paired_seeds": bool(paired_seeds),
                        "effect_rel_l2": relative_l2_effect(b.field, f.field,
                                                            mass, mask),
                        "omega_b_bins": int(mask.sum()),
                        # DIAGNOSTIC, not a membership criterion — see omega_b.
                        "reported_rel_err_median": float(
                            np.median(err)) if mask.any() else float("nan"),
                        "rel_err_frac_over_ceiling": float(
                            np.mean(err >= MAX_REL_ERR)) if mask.any()
                            else float("nan"),
                    })
            if stopped_early:
                break
        if stopped_early:
            break

    wall_total = time.perf_counter() - started
    return _report(records, debiased, args, runs, total_histories,
                   statepoint_bytes,
                   wall_total, raytrace_wall, stopped_early, base, snapshot)


def _report(records, debiased, args, runs, histories, sp_bytes, wall,
            raytrace_wall,
            stopped_early, base, snapshot) -> int:
    import openmc

    def eff(pred):
        return np.array([r["effect_rel_l2"] for r in records if pred(r)
                         and np.isfinite(r["effect_rel_l2"])], dtype=float)

    # VARIANCE PER SCENARIO, never a global maximum. Judging one scenario
    # against another's spread makes a verdict depend on an unrelated case —
    # and a global max in particular hands every small-effect scenario the
    # largest-effect scenario's discretisation error.
    scenario_names = sorted({r["scenario"] for r in records})
    factors = sorted({r["mesh_factor"] for r in records})
    budgets, verdicts = {}, {}
    for s in scenario_names:
        # EVERY BUDGET TERM IS IN E^2, the units the gate decides in. Mixing a
        # raw-norm variance into a squared-effect verdict would repeat, in a
        # subtler place, exactly the units error `metric_id` now refuses.
        within, per_draw = [], []
        for f in factors:
            for d in sorted({r["outer_draw"] for r in records}):
                rows = [x for x in debiased
                        if x["mesh_factor"] == f and x["scenario"] == s
                        and x["outer_draw"] == d
                        and x["e_squared"] is not None]
                # TRANSPORT SPREAD COMES FROM THE JACKKNIFE, not from replicate
                # scatter: the U-statistic runs over replicate PAIRS, so it
                # consumes that dimension and leaves one value per outer draw.
                within += [x["jackknife_var"] for x in rows
                           if x["jackknife_var"] is not None]
                per_draw += [(f, x["e_squared"]) for x in rows]
        transport_var = float(np.mean(within)) if within else 0.0
        calib_var = 0.0
        for f in factors:
            vals = [m[1] for m in per_draw if m[0] == f]
            if len(vals) > 1:
                calib_var = max(calib_var, float(np.var(vals, ddof=1)))
        by_factor = {f: float(np.mean([m[1] for m in per_draw if m[0] == f]))
                     for f in factors if any(m[0] == f for m in per_draw)}

        # NUMERICS IS A BOUND, NOT A VARIANCE, and this is the project's own
        # declared rule rather than a preference: `mesh_coarsening_factor`
        # carries sampling_role = convergence_axis in
        # data/uncertainty/feedback_parameter_distributions.csv, which forbids
        # drawing it — "a mesh factor is an engineering decision, and giving it
        # a probability mixes that decision into the epistemic answer", and the
        # row instructs carrying a residual discretisation BOUND instead.
        #
        # This previously computed np.var(by_factor, ddof=1) over exactly two
        # factors. A two-sample variance is algebraically (d1 - d0)^2 / 2, so
        # the published 3.8186e-04 was a squared difference wearing a
        # variance's name — and it was then divided by a genuine multi-replicate
        # transport variance and reported as a "70-115x" ratio, which compared
        # two different kinds of quantity.
        #
        # The bound is the worst squared deviation of any coarser mesh from the
        # FINEST one available. Squared to stay on the scale the budget mixes,
        # referenced to the finest rather than to the mean because coarsening
        # is monotone error, not scatter about a centre.
        if len(by_factor) > 1:
            finest = min(by_factor)
            ref = by_factor[finest]
            numerics_bound = max((v - ref) ** 2
                                 for f, v in by_factor.items() if f != finest)
        else:
            numerics_bound = 0.0

        b = VarianceBudget(transport=transport_var, numerics=numerics_bound,
                           calibration=calib_var, model_form=0.0)
        budgets[s] = {
            **b.as_dict(),
            "units": POLICY.metric_id,
            "transport_basis": "mean_delete_one_replicate_jackknife_variance",
            "numerics_basis": "max_squared_deviation_from_finest_mesh",
            "numerics_n_factors": len(by_factor),
            "numerics_per_factor_mean": {str(f): v for f, v in by_factor.items()},
        }
        draws = np.array([x["e_squared"] for x in debiased
                          if x["scenario"] == s and x["e_squared"] is not None],
                         dtype=float)
        v = (decide(draws, b, POLICY, metric_id=POLICY.metric_id)
             if draws.size else
             S0Verdict("NOT_EVALUATED", "no finite effects recorded"))
        verdicts[s] = {
            **v.as_dict(),
            "metric_id": POLICY.metric_id,
            "per_factor_mean": {str(f): v for f, v in by_factor.items()},
            # The raw unsigned norm alongside, as a DIAGNOSTIC, so the bias the
            # debiased statistic removes stays visible rather than merely
            # asserted.
            "raw_relative_l2_median": {
                "metric_id": RAW_METRIC_ID,
                **{str(f): float(np.median(
                    eff(lambda r, f=f, s=s: r["mesh_factor"] == f
                        and r["scenario"] == s)))
                   for f in factors
                   if eff(lambda r, f=f, s=s: r["mesh_factor"] == f
                          and r["scenario"] == s).size}},
        }
    budget = VarianceBudget(**{k: max(b[k] for b in budgets.values())
                               for k in ("transport", "numerics",
                                         "calibration", "model_form")})

    # THE RAW METRIC'S NOISE FLOOR, still measured and still reported — but now
    # as a DIAGNOSTIC rather than as the significance test. `noise_floor` runs
    # identical material with DECORRELATED seeds, so what it reports under the
    # raw unsigned norm is that norm's positive bias.
    nf = eff(lambda r: r["scenario"] == "noise_floor")
    floor = float(np.quantile(nf, 0.99)) if nf.size else float("nan")

    # SIGNIFICANCE NOW COMES FROM THE ESTIMATOR'S OWN STANDARD ERROR, not from
    # a comparison against another scenario. The debiased statistic is centred
    # on zero when nothing changed, so "is this distinguishable from no effect"
    # is a question about its jackknife spread — and unlike the 3x-floor rule it
    # replaces, this needs no second scenario and acquires no status as an
    # undeclared second threshold.
    for s, v in verdicts.items():
        rows = [x for x in debiased if x["scenario"] == s
                and x["e_squared"] is not None
                and x["jackknife_var"] is not None]
        v["raw_noise_floor_q99"] = floor if math.isfinite(floor) else None
        if rows:
            e2 = float(np.mean([x["e_squared"] for x in rows]))
            se = float(np.sqrt(np.mean([x["jackknife_var"] for x in rows])))
            v["e_squared_mean"] = e2
            v["e_squared_standard_error"] = se
            v["z_versus_zero"] = e2 / se if se > 0 else None
            v["distinguishable_from_zero"] = bool(se > 0 and e2 > 3.0 * se)
        else:
            v["e_squared_mean"] = v["e_squared_standard_error"] = None
            v["z_versus_zero"] = v["distinguishable_from_zero"] = None

    prov = git_provenance(str(REPO))
    budget_doc = {
        "schema_version": 1, "tier": "S0", "target_calibration": False,
        "openmc_runs": runs, "histories": histories,
        "replicate_count": args.replicates, "outer_draws": args.outer_draws,
        "mesh_factors": [int(f) for f in args.mesh_factors.split(",")],
        "raytrace_samples": args.volume_samples,
        "raytrace_wall_seconds": round(raytrace_wall, 3),
        "wall_seconds": round(wall, 3),
        "statepoint_bytes_total": int(sp_bytes),
        "statepoint_bytes_mean": int(sp_bytes / runs) if runs else 0,
        "histories_per_second": round(histories / wall, 1) if wall else None,
        "seconds_per_run": round(wall / runs, 3) if runs else None,
        "variance_by_source_worst_case": budget.as_dict(),
        "variance_by_source_per_scenario": budgets,
        "stopped_early": stopped_early,
        # `material_lever_sensitivity_rel_l2` used to be emitted here as six
        # hardcoded constants, unconditionally, whatever the run actually did.
        # No sample row, seed, replicate count or mesh factor existed for any
        # of them, and one entry disagreed with this pilot's own registry. Run
        # with `--levers` to measure the ranking into `scenarios` below, where
        # every number has a row behind it.
        "noise_floor_basis": "decorrelated_seeds_identical_material",
        "hardware_note": ("wall times are MACHINE-SPECIFIC and not portable; "
                          "the histories count is. OpenMC here is CPU-only "
                          "(MPI/OpenMP), so no GPU throughput is reported "
                          "because none can be measured."),
        "openmc_version": openmc.__version__,
        "lattice_n": int(snapshot.cell_id.shape[0]),
        "provenance": prov,
    }
    verdict_doc = {
        "schema_version": 1, "tier": "S0",
        "threshold_policy_id": POLICY.threshold_policy_id,
        "target_calibration": False,
        "claim_class": "gate_logic_validation",
        "reference_d_verdict": "NOT_EVALUATED",
        "biological_calibration": "NOT_EVALUATED",
        "scenarios": verdicts,
        "provenance": prov,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "openmc_nested_pilot_budget.json").write_text(
        json.dumps(budget_doc, indent=2, allow_nan=False) + "\n")
    (args.outdir / "openmc_nested_pilot_verdict.json").write_text(
        json.dumps(verdict_doc, indent=2, allow_nan=False) + "\n")

    if not records:
        print("no records — nothing to write")
        return 1

    # PROVENANCE TRAVELS WITH THE COMMITTED TABLE. Previously it lived only in
    # the two sibling JSONs, which sit under artifacts/ and are gitignored — so
    # the committed CSV named no commit, no OpenMC version and no nuclear data,
    # and could not be tied to the run that produced it.
    #
    # An artifact CANNOT record the commit that will contain it: it is written
    # before that commit exists, so `git_commit` names the PARENT and carries a
    # `-dirty` marker whenever the tree has uncommitted changes, which during a
    # development run it always does. The honest identity is that pair plus the
    # config and snapshot hashes, and it is written as such rather than implied.
    prov_lines = [
        f"# provenance: git_commit={prov.get('git_commit')}",
        f"#   (names the PARENT commit; an artifact cannot name the commit that contains it)",
        f"# provenance: openmc_version={openmc.__version__}",
        f"# provenance: nuclear_data={_nuclear_data_id()}",
        f"# provenance: particles={args.particles} batches={args.batches}"
        f" replicates={args.replicates} outer_draws={args.outer_draws}",
        f"# provenance: mesh_factors={args.mesh_factors}"
        f" volume_samples={args.volume_samples}",
        f"# provenance: threshold_policy_id={POLICY.threshold_policy_id}"
        f" metric_id={POLICY.metric_id}",
        f"# provenance: snapshot={args.snapshot.name} config={args.config.name}",
        f"# provenance: argv={' '.join(sys.argv[1:])}",
    ]
    import csv as _csv

    def _write(path, header, rows):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fh.write(header)
            fh.write("#\n" + "\n".join(prov_lines) + "\n")
            w = _csv.DictWriter(fh, fieldnames=list(rows[0]),
                                lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    data_dir = REPO / "data" / "calibration"
    _write(data_dir / "openmc_effect_samples.csv",
           "# openmc_effect_samples — one row per transport replicate.\n"
           "# tier S0, target_calibration = false.\n"
           "#\n"
           "# THIS TABLE IS A DIAGNOSTIC, NOT THE GATE INPUT. effect_rel_l2 is\n"
           "# an UNSIGNED norm, so it is positively biased under the null and\n"
           "# cannot return zero when nothing changed. The gate decides on the\n"
           "# debiased squared effect in openmc_debiased_effects.csv. This\n"
           "# table remains the right evidence for the N^-1/2 history scaling\n"
           "# and for what common random numbers buy.\n"
           "#\n"
           "# paired_seeds = false marks the noise-floor scenario, which runs\n"
           "# identical material with decorrelated seeds so the effect it\n"
           "# reports IS this unsigned metric's bias.\n"
           "# reported_rel_err_* are DIAGNOSTICS: rel_err is deliberately not\n"
           "# an Omega_b membership criterion, because a statistical cut made\n"
           "# the mesh factors cover different regions.\n",
           records)

    if debiased:
        _write(data_dir / "openmc_debiased_effects.csv",
               "# openmc_debiased_effects — THE GATE INPUT. One row per\n"
               "# (scenario, mesh_factor, outer_draw); there is deliberately no\n"
               "# replicate column, because the U-statistic runs over ordered\n"
               "# replicate PAIRS r != s and consumes that dimension.\n"
               "#\n"
               "# s_hat estimates ||E[d]||^2_W with no noise term, because for\n"
               "# independent replicates E[d_r^T W d_s] = E[d_r]^T W E[d_s].\n"
               "# IT MAY BE NEGATIVE, and that is the point: an unbiased\n"
               "# estimate of a non-negative quantity must fall below zero\n"
               "# about half the time when the true effect is zero. A negative\n"
               "# row is evidence of no effect, never a defect.\n"
               "#\n"
               "# e_squared = s_hat / (debiased ||D0||^2_W) and is compared\n"
               "# against delta^2. Do NOT take its square root before gating:\n"
               "# a root near zero reintroduces the positive bias this exists\n"
               "# to remove.\n"
               "# jackknife_var is the delete-one-replicate variance OF THE\n"
               "# RATIO, since the ratio is what the gate consumes.\n",
               debiased)

    print(f"\nruns={runs} histories={histories:,} wall={wall:.1f}s "
          f"({budget_doc['histories_per_second']:,} hist/s)")
    print(f"variance ({POLICY.metric_id}): {budget.as_dict()}"
          f"  dominant={budget.dominant}")
    print(f"raw unsigned-norm bias (q99 of noise_floor): {floor:.4g}"
          "   [diagnostic]")
    print(f"{'scenario':26} {'E^2':>11} {'+/- se':>10} {'z':>7}  verdict")
    for s, v in verdicts.items():
        e2, se, z = (v["e_squared_mean"], v["e_squared_standard_error"],
                     v["z_versus_zero"])
        flag = {True: "", False: "  [NOT DISTINGUISHABLE FROM ZERO]",
                None: "  [no estimate]"}[v["distinguishable_from_zero"]]
        cells = ("%11.3e %10.1e %7.1f" % (e2, se, z)
                 if e2 is not None and se is not None and z is not None
                 else f"{'-':>11} {'-':>10} {'-':>7}")
        print(f"  {s:24} {cells}  {v['verdict']}{flag}")
    if stopped_early:
        print(f"STOPPED EARLY: {stopped_early}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
