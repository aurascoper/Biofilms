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

THE PRIMARY METRIC is the mass-weighted relative L2 dose-field effect over
occupied biological material:

    E = sqrt( sum_b m_v (D1 - D0)^2 ) / sqrt( sum_b m_v D0^2 )

Not a raw maximum relative difference: near-zero bins make that a measure of the
dose floor rather than of the material change. Omega_b requires biological
occupancy, positive mass, a declared dose floor, and acceptable transport
uncertainty — all four, declared before the run.

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
from biofilm_openmc.fingerprint import transport_state_hash
from biofilm_openmc.lineage import aggregate_by_label
from biofilm_openmc.materials import (mesh_material_masses_kg,
                                      mesh_material_volumes, voxel_mass_kg)
from biofilm_openmc.mesh import resolve_mesh_dimension, upsample_field
from biofilm_openmc.snapshot import load_snapshot
from biofilm_openmc.synthetic_gate import (S0Verdict, ThresholdPolicy,
                                           VarianceBudget, decide)

POLICY = ThresholdPolicy()
# Declared BEFORE the run: a bin below this carries no usable dose information
# and would make a relative metric a measure of the floor.
DOSE_FLOOR_FRACTION = 1e-3
MAX_REL_ERR = 0.25


def relative_l2_effect(d0, d1, mass, mask) -> float:
    """The primary metric. Mass-weighted, relative, and restricted to Omega_b."""
    m, a, b = mass[mask], d0[mask], d1[mask]
    denom = math.sqrt(float(np.sum(m * a * a)))
    if denom == 0.0:
        return float("nan")
    return math.sqrt(float(np.sum(m * (b - a) ** 2))) / denom


# A coarse bin counts as biological when at least this much of its VOLUME is
# biomass. Declared, because the alternative is an arbitrary representative.
MIN_BIOMASS_VOLUME_FRACTION = 0.5


def omega_b(dose0, mass, rel_err, biomass_fraction) -> np.ndarray:
    """Biological occupancy AND positive mass AND above the dose floor AND
    acceptable transport uncertainty. All four, or the metric measures
    something other than the material change.

    OCCUPANCY COMES FROM THE EXACT CSG VOLUMES, not from subsampling the lattice.
    Taking `occ[::factor]` picks one representative voxel per coarse bin, so a
    bin that is 12% biomass and 88% medium counts as fully biological — and the
    two mesh factors then measure DIFFERENT regions, which shows up as
    discretisation error that is really a definitional artifact. Asking what
    fraction of the bin's volume is actually biomass makes the factors
    comparable, and the volumes are already computed for the mass denominator.
    """
    floor = DOSE_FLOOR_FRACTION * float(np.nanmax(dose0)) if np.isfinite(
        np.nanmax(dose0)) else 0.0
    return ((biomass_fraction >= MIN_BIOMASS_VOLUME_FRACTION)
            & (mass > 0) & (dose0 > floor) & (rel_err < MAX_REL_ERR))


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
    ap.add_argument("--mesh-factors", default="1,2")
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

    # THE FEEDBACK LEVER IS COMPOSITION, NOT DENSITY. Measured on this
    # geometry at 661.7 keV: a 35% density increase moves the mass-weighted
    # relative L2 dose effect by only 1.4%, because dose is energy PER UNIT
    # MASS — heating rose x1.306 while mass rose x1.350, so the two nearly
    # cancel and only second-order self-shielding survives.
    #
    # Composition is not automatically stronger either. At Cs-137 energies
    # Compton dominates and mu/rho tracks Z/A, which is ~0.5 for almost every
    # element, so 5% Fe moved the effect by 0.3% — the WEAKEST lever tested.
    # What moves it is a departure from Z/A = 0.5: hydrogen (Z/A = 1.0) and
    # high-Z material where photoelectric still contributes.
    #
    #     density x1.35        0.0137
    #     5% Fe                0.0033
    #     5% Gd                0.0376
    #     dehydration H->0.06  0.0459
    #     20% Gd               0.1241
    #     40% Gd               0.2242
    #
    # Gd 40% is a deliberately EXTREME synthetic material with no biological
    # pretension: its only job is to prove the gate can return a pass. 20% Gd
    # sits just above the threshold and is therefore the honest candidate for
    # the near-boundary case, whose verdict is the pilot's finding rather than
    # a declaration.
    base_elements = dict(base.materials["baseline_biomass"].elements)
    scenarios = [
        ("zero_effect", None),
        ("clearly_detectable", {"H": 0.111894, "O": 0.488106, "Gd": 0.40}),
    ]
    if args.include_near_threshold:
        # RUNS THIRD AND ONLY ON REQUEST, because its verdict depends on the
        # spread this pilot measures and therefore cannot be declared in
        # advance. Adding it before the bracketing cases were known correct
        # would have made a surprising verdict uninterpretable.
        scenarios.append(
            ("near_threshold", {"H": 0.111894, "O": 0.688106, "Gd": 0.20}))
    started = time.perf_counter()
    runs = 0
    total_histories = 0
    statepoint_bytes = 0
    records: list[dict] = []
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

        for scenario, feedback_elements in scenarios:
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
                    if state == "feedback" and feedback_elements is not None:
                        elements = tuple(sorted(feedback_elements.items()))
                    mats["baseline_biomass"] = type(bm)(
                        bm.name, float(rho), elements)
                    cfg_s = replace(cfg_f, materials=mats)
                    model = build_biofilm_cylinder_model(snapshot, cfg_s)
                    m = mesh_material_masses_kg(mesh, model, dim,
                                                volumes=volumes)
                    fields = []
                    for rep in range(args.replicates):
                        cfg_r = replace(cfg_s, seed=1000 * (draw + 1) + rep)
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
                mask = omega_b(b_fields[0].field, mass, b_fields[0].rel_err,
                               bio_fraction)
                for i, (b, f) in enumerate(zip(b_fields, f_fields)):
                    records.append({
                        "scenario": scenario, "mesh_factor": factor,
                        "outer_draw": draw, "replicate": i,
                        "density_g_cm3": float(rho),
                        "effect_rel_l2": relative_l2_effect(b.field, f.field,
                                                            mass, mask),
                        "omega_b_bins": int(mask.sum()),
                        "reported_rel_err_median": float(
                            np.median(b.rel_err[mask])) if mask.any() else float("nan"),
                    })
            if stopped_early:
                break
        if stopped_early:
            break

    wall_total = time.perf_counter() - started
    return _report(records, args, runs, total_histories, statepoint_bytes,
                   wall_total, raytrace_wall, stopped_early, base, snapshot)


def _report(records, args, runs, histories, sp_bytes, wall, raytrace_wall,
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
        within, per_draw = [], []
        for f in factors:
            for d in sorted({r["outer_draw"] for r in records}):
                v = eff(lambda r, f=f, s=s, d=d: (r["mesh_factor"] == f
                                                  and r["scenario"] == s
                                                  and r["outer_draw"] == d))
                if v.size > 1:
                    within.append(float(np.var(v, ddof=1)))
                if v.size:
                    per_draw.append((f, float(np.mean(v))))
        transport_var = float(np.mean(within)) if within else 0.0
        calib_var = 0.0
        for f in factors:
            vals = [m[1] for m in per_draw if m[0] == f]
            if len(vals) > 1:
                calib_var = max(calib_var, float(np.var(vals, ddof=1)))
        by_factor = [float(np.mean([m[1] for m in per_draw if m[0] == f]))
                     for f in factors if any(m[0] == f for m in per_draw)]
        numerics_var = (float(np.var(by_factor, ddof=1))
                        if len(by_factor) > 1 else 0.0)
        b = VarianceBudget(transport=transport_var, numerics=numerics_var,
                           calibration=calib_var, model_form=0.0)
        budgets[s] = b.as_dict()
        draws = eff(lambda r, s=s: r["scenario"] == s)
        v = decide(draws, b, POLICY) if draws.size else S0Verdict(
            "NOT_EVALUATED", "no finite effects recorded")
        verdicts[s] = {**v.as_dict(),
                       "per_factor_mean": {str(f): float(np.mean(
                           eff(lambda r, f=f, s=s: r["mesh_factor"] == f
                               and r["scenario"] == s)))
                           for f in factors}}
    budget = VarianceBudget(**{k: max(b[k] for b in budgets.values())
                               for k in ("transport", "numerics",
                                         "calibration", "model_form")})

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
        "material_lever_sensitivity_rel_l2": {
            "density_x1.35": 0.0137, "Fe_5pct": 0.0033, "Gd_5pct": 0.0376,
            "dehydration_H_0.112_to_0.06": 0.0459, "Gd_20pct": 0.1241,
            "Gd_40pct": 0.2242,
            "note": ("measured on this geometry at 661.7 keV. Dose is energy "
                     "per unit mass, so density largely cancels; Compton "
                     "dominates, so mu/rho tracks Z/A ~ 0.5 for almost "
                     "everything and only departures from it move the dose. "
                     "The strongest realistic lever is the HYDROGEN fraction, "
                     "which is what the paired wet/dry gravimetry measures."),
        },
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

    samples = REPO / "data" / "calibration" / "openmc_effect_samples.csv"
    import csv as _csv
    with open(samples, "w", newline="", encoding="utf-8") as fh:
        fh.write("# openmc_effect_samples — one row per transport replicate.\n"
                 "# tier S0, target_calibration = false. The effect is the\n"
                 "# mass-weighted relative L2 dose-field change over occupied\n"
                 "# biological material, not a maximum relative difference.\n")
        w = _csv.DictWriter(fh, fieldnames=list(records[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(records)

    print(f"\nruns={runs} histories={histories:,} wall={wall:.1f}s "
          f"({budget_doc['histories_per_second']:,} hist/s)")
    print(f"variance: {budget.as_dict()}  dominant={budget.dominant}")
    for s, v in verdicts.items():
        print(f"  {s:22} -> {v['verdict']}")
    if stopped_early:
        print(f"STOPPED EARLY: {stopped_early}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
