#!/usr/bin/env python
"""Phase B: is the residual mesh movement TALLY discretisation?

    micromamba run -n openmc-biofilms python coupling/scripts/subvoxel_refinement.py \
        --snapshot snap.h5 --config config/reference_synthetic_biofilm_e2e.toml \
        --outdir artifacts/refinement

THE QUESTION. The pilot found the effect rising monotonically as the tally mesh
COARSENS — 0.1224 at the lattice resolution to 0.1821 at factor 5 — and could
not say whether the lattice resolution itself is converged, because factor 1 was
the finest tally the code could build. This study goes the other way: 1x, 2x and
4x subvoxel tally bins per CPM voxel, on a FIXED snapshot.

WHAT REFINEMENT DOES AND DOES NOT ADD. The material lattice stays piecewise
constant at the CPM pitch, so a finer tally adds TRANSPORT resolution and never
biological detail. It resolves how dose varies WITHIN one parcel of uniform
material. That is exactly the discretisation question and it is emphatically not
a morphology question — a finer snapshot is a separate study with a separate
confound (`V_target` is in sites, so refining N shrinks every parcel).

HOW THE RATIOS ARE COMPARED. Not at their own resolutions, which would compare
different regions over different bin counts. Heating and mass are both EXTENSIVE,
so each ratio's fields are block-summed back to the lattice grid and divided
afterwards (`mesh.coarsen_ratio`), giving the mass-weighted mean dose the fine
calculation implies for each lattice voxel. Omega_b is then defined ONCE, from
the lattice-resolution raytrace, so every ratio is evaluated over an identical
region by construction.

THE STATISTIC is the debiased squared effect, so a ratio that is noisier because
its bins are smaller does not thereby look like a larger effect.

RAY COUNT IS NOT SCALED NAIVELY. Each ray crosses `dim` bins along its axis, so
`n_samples = P * dim^2` delivers about P samples per bin — the sampling per bin
does not fall as dim^3. It also must not be raised carelessly: `material_volumes`
lays rays on a regular grid, and at some counts a ray passes within
floating-point epsilon of a lattice plane and OpenMC ABORTS THE PROCESS with
"could not be located after crossing a boundary of lattice". That abort is not
catchable and it is not a property of refinement — measured here, 20^3 bins abort
at 2.56e6 rays while 80^3 bins succeed at 5e6. The counts below are verified.
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

from biofilm_openmc.config import BIOFILM_CYLINDER, load_dose_rate_config
from biofilm_openmc.dose import per_source_from_statepoint
from biofilm_openmc.feedback_uq import (cross_replicate_inner,
                                        debiased_squared_effect,
                                        relative_l2_effect)
from biofilm_openmc.fingerprint import transport_state_hash
from biofilm_openmc.materials import mesh_material_masses_kg, mesh_material_volumes
from biofilm_openmc.viewer import Grid, Layer, write_bundle
from biofilm_openmc.mesh import (biofilm_mesh_extent_cm, coarsen_field,
                                 coarsen_ratio, resolve_mesh_dimension,
                                 upsample_field)

# Declared before the run, and identical to the pilot's so the two are comparable.
MIN_BIOMASS_VOLUME_FRACTION = 0.5

# Samples per mesh bin from the CSG raytrace. 250 is ample: the volumes are
# nearly binary at fine resolution, since a bin well inside a lattice cell holds
# exactly one material and only cylinder/membrane boundary bins are mixed.
RAYS_PER_BIN = 250
MIN_RAY_SAMPLES = 500_000

# The contrast whose mesh dependence is in question, and the null control.
SCENARIOS = {
    "near_threshold": {"H": 0.111894, "O": 0.688106, "Gd": 0.20},
    "noise_floor": None,
}


def ray_samples(dimension) -> int:
    """P samples per bin, because each ray crosses `dim` bins along its axis."""
    return max(MIN_RAY_SAMPLES, RAYS_PER_BIN * int(dimension[0]) ** 2)


def biomass_volumes(volumes: dict, dimension):
    """Per-bin (biomass, total) material VOLUMES, logical (x,y,z).

    Returned as the two extensive arrays rather than as their ratio, because a
    coarse bin's biomass fraction is sum(bio) / sum(total) over the block it
    covers — NOT the mean of its fine bins' fractions. Those differ whenever the
    fine bins hold different total material volume, which is exactly what
    happens along the cylinder and membrane boundaries.
    """
    dim = tuple(int(d) for d in dimension)
    total = np.zeros(int(np.prod(dim)), dtype=float)
    bio = np.zeros_like(total)
    for name, per_element in volumes.items():
        arr = np.asarray(per_element, dtype=float)
        total += arr
        if name == "baseline_biomass":
            bio += arr
    shape = lambda a: a.reshape(dim[::-1]).transpose(2, 1, 0)
    return shape(bio), shape(total)


def volume_fraction(bio, total) -> np.ndarray:
    return np.divide(bio, total, out=np.zeros_like(bio), where=total > 0)


def omega_b(mass, biomass_fraction) -> np.ndarray:
    """Two criteria, both GEOMETRIC — see the pilot's `omega_b`.

    A dose floor taken from one replicate makes the region a random set drawn
    from that replicate, so the cross-replicate independence the debiased
    estimator rests on fails for every ordered pair involving it.
    """
    return ((biomass_fraction >= MIN_BIOMASS_VOLUME_FRACTION) & (mass > 0))


def _seed(rep: int, state: str, paired: bool) -> int:
    """Common random numbers across states, except for the null control, which
    must be able to fail. Same rule as the pilot."""
    return 7000 + rep + (0 if paired or state == "baseline" else 500_000)


def _write_bundle(args, snapshot, grids, layers, lattice_grid, common_mask,
                  config, openmc_version) -> None:
    """Assemble the multi-grid bundle from what this run already computed.

    Costs no extra transport. The point is that the refinement study is the
    first thing in the repository to produce a genuinely NATIVE high-resolution
    dose field — before it, a viewer could only have shown the lattice grid, and
    anything finer would have been a coarse field broadcast to look detailed.
    """
    from biofilm_openmc.results import provenance_attrs

    labels = [
        Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
              snapshot.cell_id.astype(np.int32),
              note="a computational biomass PARCEL id, never an organism"),
        Layer("lineage_id", "cpm_labels", "dimensionless", "categorical",
              snapshot.lineage_id.astype(np.int32)),
        Layer("generation", "cpm_labels", "dimensionless", "categorical",
              snapshot.generation.astype(np.int32),
              note="generation 0 is a founder, not a missing value"),
        Layer("omega_b", "cpm_labels", "dimensionless", "boolean",
              np.asarray(common_mask, dtype=bool), derivation=None,
              note="the metric's region: biomass volume fraction >= 0.5, "
                   "positive mass, above the dose floor. No uncertainty cut - "
                   "a statistical criterion on a region definition makes the "
                   "region move with the history count"),
    ]
    path = args.outdir / "viewer_bundle.h5"
    doc = write_bundle(
        path, [lattice_grid] + grids, labels + layers,
        provenance={**provenance_attrs(config, repo_root=str(REPO)),
                    "openmc_version": openmc_version,
                    "study": "subvoxel_tally_refinement",
                    "particles": args.particles, "batches": args.batches})
    print(f"wrote {path} - {len(doc['grids'])} grids, {len(doc['layers'])} layers")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--ratios", default="1,2,4")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--particles", type=int, default=20000)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--budget-seconds", type=float, default=3600.0)
    ap.add_argument("--bundle", action="store_true",
                    help="also write a multi-grid viewer bundle. Costs no "
                         "extra transport: the native dose field at every "
                         "ratio is already in memory here")
    args = ap.parse_args(argv)

    import openmc

    from biofilm_openmc.model import build_biofilm_cylinder_model
    from biofilm_openmc.snapshot import load_snapshot

    args.outdir.mkdir(parents=True, exist_ok=True)
    base = load_dose_rate_config(args.config.with_suffix(".dosimetry.toml"),
                                 kind=BIOFILM_CYLINDER)
    snapshot = load_snapshot(args.snapshot)
    n = snapshot.cell_id.shape[0]
    # Sorted so ratio 1 runs first: it defines the common Omega_b that every
    # finer ratio is then evaluated on.
    ratios = sorted({int(r) for r in args.ratios.split(",")})
    if ratios[0] != 1:
        raise SystemExit("ratio 1 is the reference grid and must be included")
    # EVERY RATIO MUST DIVIDE THE FINEST, checked before any transport. The mass
    # denominator is raytraced once at the finest ratio and coarsened by
    # `finest // ratio`; integer division silently gives step 1 for a ratio like
    # 3 under a finest of 4, pairing a 4x mass array with a 3x tally. The
    # mismatch surfaces only when dose conversion combines them -- after the
    # raytrace and the histories have been paid for.
    bad = [r for r in ratios if ratios[-1] % r]
    if bad:
        raise SystemExit(
            f"ratios {bad} do not divide the finest ratio {ratios[-1]}: the "
            "shared mass denominator is coarsened by finest // ratio, so a "
            "non-divisor pairs mismatched arrays and fails only after transport")

    started = time.perf_counter()
    runs = total_histories = statepoint_bytes = 0
    rows: list[dict] = []
    common_mask = None          # Omega_b at lattice resolution, defined once
    lattice_grid = None
    bundle_grids = [] if args.bundle else None
    bundle_layers = []
    reference_remapped: dict[str, float] = {}
    reference_lattice_dose: dict = {}
    field_drift: dict = {}
    stopped_early = None

    # ONE RAYTRACE, AT THE FINEST RATIO, COARSENED DOWN FOR THE REST.
    #
    # Raytracing each ratio independently looks natural and is wrong. The mass
    # denominator is an O(1/sqrt(n)) ESTIMATE, so independent raytraces differ:
    # measured here, up to 3.0% per voxel and 0.13% in total between a 1x and a
    # coarsened-4x mass. Dose is heating over mass, so that error lands directly
    # in every dose field and makes the ratios differ for a reason that has
    # nothing to do with resolution.
    #
    # Mass is EXTENSIVE, so coarsening the finest estimate is exact and gives
    # every ratio a mutually consistent denominator by construction. It is also
    # cheaper: one raytrace instead of one per ratio.
    finest = ratios[-1]
    cfg_fine = replace(base, mesh_coarsening_factor=1,
                       mesh_refinement_factor=finest,
                       batches=args.batches, particles=args.particles)
    dim_fine = resolve_mesh_dimension((n, n, n), 1, finest)
    probe_fine = build_biofilm_cylinder_model(snapshot, cfg_fine)
    t0 = time.perf_counter()
    n_rays = ray_samples(dim_fine)
    volumes_fine = mesh_material_volumes(
        probe_fine.tallies[0].filters[0].mesh, probe_fine, n_samples=n_rays)
    raytrace_wall = time.perf_counter() - t0
    mass_fine = mesh_material_masses_kg(
        probe_fine.tallies[0].filters[0].mesh, probe_fine, dim_fine,
        volumes=volumes_fine)
    bio_fine, total_fine = biomass_volumes(volumes_fine, dim_fine)

    for ratio in ratios:
        # Coarsening is forced to 1: the two are opposite ends of one axis and
        # the resolver refuses both, so a config carrying coarsening_factor = 2
        # would otherwise silently contradict the ratio being swept.
        cfg_r = replace(base, mesh_coarsening_factor=1, mesh_refinement_factor=ratio,
                        batches=args.batches, particles=args.particles)
        dim = resolve_mesh_dimension((n, n, n), 1, ratio)
        probe = build_biofilm_cylinder_model(snapshot, cfg_r)
        mesh = probe.tallies[0].filters[0].mesh

        step = finest // ratio
        mass_ratio = coarsen_field(mass_fine, step)
        extent_cm = biofilm_mesh_extent_cm(cfg_r, n)
        bio_vol, total_vol = coarsen_field(bio_fine, step), coarsen_field(total_fine, step)
        bio_fraction = volume_fraction(bio_vol, total_vol)

        for scenario, feedback_elements in SCENARIOS.items():
            paired = feedback_elements is not None
            if time.perf_counter() - started > args.budget_seconds:
                stopped_early = f"budget exhausted before r{ratio}/{scenario}"
                break
            per_state = {}
            for state in ("baseline", "feedback"):
                mats = dict(cfg_r.materials)
                bm = mats["baseline_biomass"]
                elements = bm.elements
                if state == "feedback" and feedback_elements is not None:
                    elements = tuple(sorted(feedback_elements.items()))
                mats["baseline_biomass"] = type(bm)(bm.name, bm.density_g_cm3,
                                                    elements)
                cfg_s = replace(cfg_r, materials=mats)
                model = build_biofilm_cylinder_model(snapshot, cfg_s)
                mass = mass_ratio
                fields = []
                for rep in range(args.replicates):
                    cfg_x = replace(cfg_s, seed=_seed(rep, state, paired))
                    mdl = build_biofilm_cylinder_model(snapshot, cfg_x)
                    cwd = args.outdir / f"r{ratio}_{scenario}_{state}_{rep}"
                    cwd.mkdir(parents=True, exist_ok=True)
                    sp_path = mdl.run(cwd=str(cwd), output=False)
                    with openmc.StatePoint(sp_path) as sp:
                        fields.append(per_source_from_statepoint(sp, mass))
                    runs += 1
                    total_histories += cfg_x.batches * cfg_x.particles
                    statepoint_bytes += Path(sp_path).stat().st_size
                per_state[state] = (fields, mass)

            (b_fields, mass), (f_fields, _) = (per_state["baseline"],
                                               per_state["feedback"])

            # Heating is recovered as dose * mass, both extensive, and summed
            # back to the lattice grid. `coarsen_ratio` does the divide after
            # the sums, which is what makes the coarse value the mass-weighted
            # mean of the fine dose rather than its unweighted average.
            lattice_mass = coarsen_field(mass, ratio)
            b_lat = [coarsen_ratio(f.field * mass, mass, ratio) for f in b_fields]
            f_lat = [coarsen_ratio(f.field * mass, mass, ratio) for f in f_fields]

            if bundle_grids is not None and scenario == "near_threshold":
                # THE BASELINE FIELD AT ITS OWN RESOLUTION, kept as measured.
                # Nothing is resampled onto the lattice here: the whole point of
                # the bundle is that a 4x dose field has structure the CPM grid
                # cannot hold, and flattening it to fit would destroy exactly
                # what makes it worth showing.
                gid = f"dose_refinement_{ratio}"
                bundle_grids.append(Grid(
                    gid, tuple(int(d) for d in dim),
                    tuple(float(v) for v in base.origin_cm),
                    tuple(float(e) / int(d) for e, d in zip(extent_cm, dim)),
                    material_resolution_grid="cpm_labels",
                    note=("native OpenMC tally; the MATERIAL it integrates is "
                          "piecewise constant at the CPM pitch, so sub-voxel "
                          "structure here is transport structure, not biology")))
                bundle_layers.append(Layer(
                    f"dose_per_source_r{ratio}", gid, "Gy/source-particle",
                    "intensive", b_fields[0].field))
                bundle_layers.append(Layer(
                    f"mass_kg_r{ratio}", gid, "kg", "extensive", mass,
                    note="exact CSG raytraced volumes x density"))
                bundle_layers.append(Layer(
                    f"rel_err_r{ratio}", gid, "dimensionless", "intensive",
                    b_fields[0].rel_err,
                    note="tally relative error; a DIAGNOSTIC, never an "
                         "Omega_b membership criterion"))
                if ratio > 1:
                    # BOTH SIDES OF THE RULE, on real data, because a viewer
                    # overlaying dose on labels needs exactly these two.
                    #
                    # Reducing the fine dose onto the lattice is exact - it
                    # reproduces the native ratio-1 field to nine significant
                    # figures - so it stays quotable.
                    bundle_layers.append(Layer(
                        f"dose_on_lattice_from_r{ratio}", "cpm_labels",
                        "Gy/source-particle", "intensive", b_lat[0],
                        native=False, authoritative_for_quantitation=True,
                        source_grid_id=gid, derivation="mass_weighted_mean",
                        note="block-summed heating over block-summed mass; "
                             "exact, and equal to the native lattice field"))
                    # Broadcasting labels the other way is exact too, because
                    # the material lattice is piecewise constant. It is still
                    # marked unquotable: the direction is what the rule keys
                    # on, and a reader who learns to trust ONE upsampled layer
                    # will trust the next one, which will be a dose field.
                    bundle_layers.append(Layer(
                        f"cell_id_on_r{ratio}", gid, "dimensionless",
                        "categorical",
                        upsample_field(snapshot.cell_id.astype(np.int32), ratio),
                        native=False, authoritative_for_quantitation=False,
                        source_grid_id="cpm_labels",
                        derivation="display_resampling",
                        note="for overlay only. Exact here, since the lattice "
                             "is piecewise constant, but marked unquotable "
                             "because the DIRECTION is what the rule keys on"))

            if common_mask is None:
                # Defined ONCE, at ratio 1, and reused unchanged for every
                # finer ratio, so the comparison is over an identical region by
                # construction rather than by an argument that the regions are
                # close. The fraction is sum(bio)/sum(total) over the block.
                lattice_fraction = volume_fraction(coarsen_field(bio_vol, ratio),
                                                   coarsen_field(total_vol, ratio))
                common_mask = omega_b(lattice_mass, lattice_fraction)
                lattice_grid = Grid(
                    "cpm_labels", (n, n, n),
                    tuple(float(v) for v in base.origin_cm),
                    (base.voxel_pitch_cm,) * 3,
                    note="the CPM material lattice; every material boundary in "
                         "the model lies on this grid")

            # THE FINE REGION IS THE COARSE REGION, UPSAMPLED. Using each
            # ratio's own Omega_b would compare different physical regions, and
            # the fine one moves because per-bin statistics are worse. Taking
            # the lattice mask and expanding each voxel into its r^3 sub-bins
            # makes the two evaluations cover the identical physical volume, so
            # the only thing that differs is the resolution at which the effect
            # is computed.
            fine_mask = upsample_field(common_mask, ratio)
            native = debiased_squared_effect([f.field for f in b_fields],
                                             [f.field for f in f_fields],
                                             mass, fine_mask)
            remapped = debiased_squared_effect(b_lat, f_lat, lattice_mass,
                                               common_mask)
            # CONSERVATION CHECK, NOT A RESULT. The tally mesh is a histogram
            # of the same events: with the seed independent of the ratio, every
            # ratio transports bit-identical histories, and block-summing a
            # finer histogram recovers the coarse bin exactly. So the remapped
            # value MUST equal the ratio-1 value to floating-point. If it ever
            # does not, either the remap is not conservative or the refinement
            # perturbed the transport, and both are defects.
            # THE FIELD ITSELF, not only the statistic. An earlier version
            # checked only E^2 and read its agreement as proof that the dose
            # field was reproduced. It was not: E^2 is a ratio in which the
            # mass denominator appears in both sums and cancels to high order,
            # so it agreed to nine figures while the underlying dose field
            # differed by 1.5% from independently raytraced masses. Comparing
            # the FIELD is the claim worth making, and it only became true once
            # every ratio shared one raytrace.
            if ratio == 1:
                reference_lattice_dose[scenario] = b_lat[0].copy()
            else:
                ref_field = reference_lattice_dose.get(scenario)
                if ref_field is not None:
                    peak = float(np.abs(ref_field).max())
                    drift = (float(np.abs(b_lat[0] - ref_field).max()) / peak
                             if peak else 0.0)
                    # 1e-7, and the slack has a MECHANISM rather than being
                    # numerical slop. Measured: 0 at r=1, 1.5e-13 at r=2, and
                    # 3.0e-09 at r=4. Float accumulation would predict the 8-term
                    # and 64-term sums to differ by ~8x, not by four orders of
                    # magnitude, so it is not accumulation.
                    #
                    # It is the ORPHANED-ENERGY DISCARD, and it is
                    # resolution-dependent. `specific_energy_per_source` drops
                    # heating found in bins the raytrace gave zero volume. At
                    # r=4 a sliver sub-bin is its own bin and its energy is
                    # discarded; at r=1 that same energy sits inside a coarse
                    # voxel with plenty of mass and is kept. Localised exactly
                    # where predicted: of the voxels carrying drift above 1% of
                    # the maximum there is EXACTLY ONE, and it holds 28 of 64
                    # massless sub-bins, against 19% of the remaining voxels
                    # holding any. Mass itself conserves to 4.8e-16.
                    #
                    # So this check detects more than a bad remap: it also
                    # bounds how much the sliver tolerance costs. The defect it
                    # was written for - each ratio raytracing its own mass -
                    # showed up at 1.5e-02, five orders of magnitude above.
                    field_drift[(scenario, ratio)] = drift
                    if drift > 1e-7:
                        raise SystemExit(
                            f"reduced dose field at r={ratio} ({scenario}) "
                            f"differs from the native lattice field by "
                            f"{drift:.3e} — the remap is not conservative, or "
                            "the ratios do not share a mass denominator")

            if ratio == 1:
                reference_remapped[scenario] = remapped.e_squared
            else:
                ref = reference_remapped.get(scenario)
                if ref is not None and np.isfinite(ref) and ref != 0:
                    drift = abs(remapped.e_squared - ref) / abs(ref)
                    # 1e-6, not tighter: the fine grid sums 8x or 64x more
                    # terms, so float accumulation alone moves the near-null
                    # scenarios at the 1e-8 level. Genuine non-conservation --
                    # a non-extensive remap, or refinement perturbing the
                    # transport -- is a percent-level effect and still caught.
                    if drift > 1e-6:
                        raise SystemExit(
                            f"remap is not conservative: r={ratio} {scenario} "
                            f"gives {remapped.e_squared:.8e} against the "
                            f"ratio-1 reference {ref:.8e} (relative {drift:.2e})")

            # How much dose structure lives BELOW the lattice pitch — DEBIASED
            # the same way as everything else. A single replicate's residual
            # about its voxel mean is (real sub-voxel structure + Monte Carlo
            # noise), and at r = 4 the per-bin relative error reaches 0.74, so
            # the naive version measures almost pure noise and reports it as
            # structure. The cross-replicate inner product removes it.
            w = mass.ravel()[fine_mask.ravel()]
            resid = np.array([(f.field - upsample_field(lat, ratio)).ravel()[
                fine_mask.ravel()] for f, lat in zip(b_fields, b_lat)])
            level = np.array([f.field.ravel()[fine_mask.ravel()]
                              for f in b_fields])
            sub_power = cross_replicate_inner(resid, w)
            tot_power = cross_replicate_inner(level, w)
            subvoxel_rms = (math.copysign(math.sqrt(abs(sub_power)), sub_power)
                            if np.isfinite(sub_power) else float("nan"))
            voxel_rms = (math.sqrt(tot_power) if np.isfinite(tot_power)
                         and tot_power > 0 else float("nan"))

            rows.append({
                "scenario": scenario, "ratio": ratio,
                "subvoxel_dose_rms_over_voxel_rms": (
                    subvoxel_rms / voxel_rms if voxel_rms else float("nan")),
                "mesh_dimension": int(dim[0]), "bins": int(np.prod(dim)),
                "ray_samples": n_rays,
                "mass_basis": "single_raytrace_at_finest_ratio_coarsened",
                "lattice_field_drift_vs_ratio1": field_drift.get(
                    (scenario, ratio), 0.0),
                "raytrace_seconds": round(raytrace_wall, 2),
                "native_omega_b_bins": int(fine_mask.sum()),
                "common_omega_b_bins": int(common_mask.sum()),
                "native_e_squared": native.e_squared,
                "native_jackknife_var": native.jackknife_var,
                # Reported separately because e_squared is a RATIO of two
                # estimators, and a ratio is only approximately unbiased: the
                # approximation degrades as the denominator's own noise grows,
                # which is exactly what refinement does to it.
                "native_s_hat": native.s_hat,
                "native_denominator": native.denominator,
                "remapped_e_squared": remapped.e_squared,
                "remapped_jackknife_var": remapped.jackknife_var,
                "raw_relative_l2_lattice": float(np.median(
                    [relative_l2_effect(b, f, lattice_mass, common_mask)
                     for b, f in zip(b_lat, f_lat)])),
                "median_rel_err": float(np.median(
                    b_fields[0].rel_err[fine_mask])) if fine_mask.any()
                    else float("nan"),
                "transport_state_hash": transport_state_hash(
                    snapshot, cfg_r, "endfb-viii.0")[:16],
            })
            print(f"  r={ratio} {scenario:15} dim={dim[0]:>3}^3 "
                  f"native E2={native.e_squared:+.5e} "
                  f"remapped E2={remapped.e_squared:+.5e}", flush=True)
        if stopped_early:
            break

    if bundle_grids is not None and bundle_grids:
        _write_bundle(args, snapshot, bundle_grids, bundle_layers,
                      lattice_grid, common_mask, base, openmc.__version__)

    wall = time.perf_counter() - started
    doc = {
        "schema_version": 1, "tier": "S0", "target_calibration": False,
        "study": "subvoxel_tally_refinement",
        "question": ("does the effect keep moving as the TALLY is refined below "
                     "the CPM lattice pitch, at fixed snapshot and geometry"),
        "ratios": ratios, "replicates": args.replicates,
        "particles": args.particles, "batches": args.batches,
        "openmc_runs": runs, "histories": total_histories,
        "wall_seconds": round(wall, 1),
        "statepoint_bytes_total": int(statepoint_bytes),
        "comparison_basis": "block_sum_heating_and_mass_to_lattice_then_divide",
        "omega_b_basis": "defined_once_at_lattice_resolution",
        "metric_id": "debiased_relative_l2_squared",
        "openmc_version": openmc.__version__,
        "lattice_n": int(n),
        "stopped_early": stopped_early,
        "rows": rows,
        "provenance": git_provenance(str(REPO)),
    }
    (args.outdir / "subvoxel_refinement.json").write_text(
        json.dumps(doc, indent=2, allow_nan=False) + "\n")

    print(f"\nruns={runs} histories={total_histories:,} wall={wall:.1f}s")
    print(f"{'scenario':16}{'ratio':>6}{'dim':>6}{'remapped E^2':>15}{'+/- se':>11}")
    for r in rows:
        se = (np.sqrt(r["remapped_jackknife_var"])
              if r["remapped_jackknife_var"] is not None
              and np.isfinite(r["remapped_jackknife_var"]) else float("nan"))
        print(f"{r['scenario']:16}{r['ratio']:>6}{r['mesh_dimension']:>6}"
              f"{r['remapped_e_squared']:>15.5e}{se:>11.2e}")
    if stopped_early:
        print(f"STOPPED EARLY: {stopped_early}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
