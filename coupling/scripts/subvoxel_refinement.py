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
from biofilm_openmc.mesh import (coarsen_field, coarsen_ratio,
                                 resolve_mesh_dimension, upsample_field)

# Declared before the run, and identical to the pilot's so the two are comparable.
DOSE_FLOOR_FRACTION = 1e-3
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


def omega_b(dose0, mass, biomass_fraction) -> np.ndarray:
    """Three PHYSICAL criteria. Deliberately no uncertainty cut — see the
    pilot's `omega_b`: a statistical criterion on a region definition makes the
    region move with the history count."""
    peak = float(np.nanmax(dose0))
    floor = DOSE_FLOOR_FRACTION * peak if np.isfinite(peak) else 0.0
    return ((biomass_fraction >= MIN_BIOMASS_VOLUME_FRACTION)
            & (mass > 0) & (dose0 > floor))


def _seed(rep: int, state: str, paired: bool) -> int:
    """Common random numbers across states, except for the null control, which
    must be able to fail. Same rule as the pilot."""
    return 7000 + rep + (0 if paired or state == "baseline" else 500_000)


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

    started = time.perf_counter()
    runs = total_histories = statepoint_bytes = 0
    rows: list[dict] = []
    common_mask = None          # Omega_b at lattice resolution, defined once
    reference_remapped: dict[str, float] = {}
    stopped_early = None

    for ratio in ratios:
        # Coarsening is forced to 1: the two are opposite ends of one axis and
        # the resolver refuses both, so a config carrying coarsening_factor = 2
        # would otherwise silently contradict the ratio being swept.
        cfg_r = replace(base, mesh_coarsening_factor=1, mesh_refinement_factor=ratio,
                        batches=args.batches, particles=args.particles)
        dim = resolve_mesh_dimension((n, n, n), 1, ratio)
        probe = build_biofilm_cylinder_model(snapshot, cfg_r)
        mesh = probe.tallies[0].filters[0].mesh

        t0 = time.perf_counter()
        n_rays = ray_samples(dim)
        volumes = mesh_material_volumes(mesh, probe, n_samples=n_rays)
        raytrace_wall = time.perf_counter() - t0
        bio_vol, total_vol = biomass_volumes(volumes, dim)
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
                mass = mesh_material_masses_kg(mesh, model, dim, volumes=volumes)
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

            if common_mask is None:
                # Defined ONCE, at ratio 1, and reused unchanged for every
                # finer ratio, so the comparison is over an identical region by
                # construction rather than by an argument that the regions are
                # close. The fraction is sum(bio)/sum(total) over the block.
                lattice_fraction = volume_fraction(coarsen_field(bio_vol, ratio),
                                                   coarsen_field(total_vol, ratio))
                common_mask = omega_b(b_lat[0], lattice_mass, lattice_fraction)

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
