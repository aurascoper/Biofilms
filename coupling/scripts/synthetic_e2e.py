#!/usr/bin/env python
"""Run the synthetic one-way biofilm dose chain, and check what it proves.

    snapshot -> OpenMC biofilm-cylinder run -> heating eV/source
      -> exact CSG mass -> Gy/source -> Gy/s -> dose field for Julia
      -> cell / lineage / generation attribution -> STOP

IT STOPS THERE, deliberately. The Hamiltonian, melanin production, membrane
damage and CPM time advancement are not applied. `import_dose_field!` is called
on the Julia side with both transforms `nothing`, and `accrue_dose!` stays
blocked on `seconds_per_mcs = NaN`, which is the contract's own guard rather
than a flag added here.

The run licenses one claim — that the architecture executed — and the verdict
JSON says so in fields rather than prose, because a wording edit must not be
able to flip a scientific safety check.

    export OPENMC_CROSS_SECTIONS=$HOME/openmc-data/endfb-viii.0-hdf5/cross_sections.xml
    micromamba run -n openmc-biofilms python coupling/scripts/synthetic_e2e.py \
        --snapshot snap.h5 --config config/reference_synthetic_biofilm_e2e.toml \
        --outdir /tmp/e2e
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from biofilm_openmc.config import (BIOFILM_CYLINDER, load_dose_rate_config,
                                   load_transport_config)
from biofilm_openmc.dose import (dose_from_per_source, per_source_from_statepoint,
                                 sparsity_report)
from biofilm_openmc.fingerprint import (dose_state_hash, needs_rerun,
                                        transport_state_hash)
from biofilm_openmc.lineage import aggregate_by_label
from biofilm_openmc.materials import mesh_material_masses_kg, voxel_mass_kg
from biofilm_openmc.mesh import (biofilm_mesh_extent_cm, coarsen_field,
                                 resolve_mesh_dimension, upsample_field)
from biofilm_openmc.model import build_biofilm_cylinder_model, source_position_cm
from biofilm_openmc.results import (provenance_attrs, write_dose_attribution,
                                    write_transport_result)
from biofilm_openmc.snapshot import load_snapshot, write_dose_field

REPO = Path(__file__).resolve().parents[2]
# Group separation for the ordering invariant, in standard deviations. Declared
# rather than tuned: a Monte Carlo field is noisy and a strict per-parcel
# ordering would be a test of the seed, not of the physics.
ORDERING_K = 3.0


def sha256(path: Path) -> str:
    with open(path, "rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


class Checks:
    """Invariants, recorded as data. A check that only prints is a check that
    can be read past."""

    def __init__(self):
        self.results: list[dict] = []

    def record(self, name: str, ok: bool, detail: str) -> bool:
        self.results.append({"invariant": name, "passed": bool(ok),
                             "detail": detail})
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        return bool(ok)

    @property
    def all_passed(self) -> bool:
        return all(r["passed"] for r in self.results)


def run_transport(model, outdir: Path, name: str):
    """Same invocation `drivers._openmc_runner` uses — one run directory per
    named problem, statepoint reopened from the path the run reports."""
    import openmc

    cwd = outdir / name
    cwd.mkdir(parents=True, exist_ok=True)
    return openmc.StatePoint(model.run(cwd=str(cwd), output=False))


def _write_viewer_bundle(args, cfg, snapshot, dose, dose_lattice, sd_lattice,
                         lattice_mass, occupied, n) -> None:
    """Package what this run already computed as a multi-grid viewer bundle.

    TWO GRIDS, NOT ONE, and the direction of every resampling is declared.
    The native dose lives on the OpenMC tally grid; the CPM labels live on the
    lattice. The lattice-resolution dose field this script already builds for
    attribution is an UPSAMPLING -- `upsample_field` broadcasts one coarse
    value across every fine voxel in its bin -- so it is written
    `authoritative_for_quantitation=False` with `derivation="upsampled_coarse_dose"`,
    and `viewer.bundle_problems` refuses the alternative outright. Displayed on
    the lattice it looks like per-parcel dose; it is one coarse value repeated.

    The point is not the picture. Putting both grids in one bundle subjects
    them to `bundle_problems`' coordinate check -- same physical volume, one
    origin -- which this path never had, because dose and labels were
    registered only by index arithmetic before now.
    """
    from biofilm_openmc.mesh import resolve_mesh_dimension
    from biofilm_openmc.viewer import Grid, Layer, write_bundle

    pitch = cfg.voxel_pitch_cm
    origin = tuple(float(v) for v in cfg.origin_cm)
    extent = biofilm_mesh_extent_cm(cfg, n)
    dim = resolve_mesh_dimension((n, n, n), cfg.mesh_coarsening_factor,
                                 getattr(cfg, "mesh_refinement_factor", 1))
    dose_spacing = tuple(e / d for e, d in zip(extent, dim))

    lattice_grid = Grid("cpm_labels", (n, n, n), origin, (pitch, pitch, pitch))
    dose_grid = Grid("dose_tally", tuple(int(d) for d in dim), origin,
                     dose_spacing, material_resolution_grid="cpm_labels",
                     note="the OpenMC tally mesh; same cube as the lattice, "
                          "coarsened by mesh_coarsening_factor")

    layers = [
        Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
              snapshot.cell_id.astype(np.int32), background=0,
              note="a computational biomass PARCEL id, never an organism"),
        Layer("lineage_id", "cpm_labels", "dimensionless", "categorical",
              snapshot.lineage_id.astype(np.int32), occupancy_from="cell_id"),
        Layer("generation", "cpm_labels", "dimensionless", "categorical",
              snapshot.generation.astype(np.int32), occupancy_from="cell_id",
              note="generation 0 is a founder AND the zero-fill for empty "
                   "sites; cell_id is what says which voxels hold biomass"),
        Layer("occupied", "cpm_labels", "dimensionless", "boolean",
              np.asarray(occupied, dtype=bool), background=0),
        Layer("mass_kg", "cpm_labels", "kg", "extensive", lattice_mass),
        Layer("dose_rate_Gy_s", "dose_tally", "Gy/s", "intensive",
              dose.dose_rate_mean_Gy_s,
              note="NATIVE here. Gy per second, because a source rate was "
                   "supplied; the per-source field is a different quantity"),
        Layer("dose_rate_sd_Gy_s", "dose_tally", "Gy/s", "intensive",
              dose.dose_rate_sd_Gy_s),
        Layer("dose_rate_on_lattice", "cpm_labels", "Gy/s", "intensive",
              dose_lattice, native=False,
              authoritative_for_quantitation=False,
              source_grid_id="dose_tally",
              derivation="upsampled_coarse_dose",
              note="what attribution actually consumed. One coarse value "
                   "repeated across each bin -- drawable, never quotable"),
        Layer("dose_rate_sd_on_lattice", "cpm_labels", "Gy/s", "intensive",
              sd_lattice, native=False,
              authoritative_for_quantitation=False,
              source_grid_id="dose_tally",
              derivation="upsampled_coarse_dose"),
    ]

    path = args.outdir / "viewer_bundle.h5"
    doc = write_bundle(path, [lattice_grid, dose_grid], layers,
                       provenance={**provenance_attrs(cfg, repo_root=str(REPO)),
                                   "study": "synthetic_one_way_e2e"})
    print(f"wrote {path} - {len(doc['grids'])} grids, "
          f"{len(doc['layers'])} layers")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--dose-config", type=Path,
                    help="dosimetry-stage config; defaults to "
                         "<config>.dosimetry.toml beside the transport one")
    ap.add_argument("--volume-samples", type=int, default=10_000_000,
                    help="rays for the exact CSG mass; the denominator's own "
                         "sampling error scales as 1/sqrt of this")
    ap.add_argument("--verdict", type=Path)
    args = ap.parse_args(argv)

    import openmc

    args.outdir.mkdir(parents=True, exist_ok=True)
    dose_config_path = args.dose_config or args.config.with_suffix(".dosimetry.toml")

    cfg = load_transport_config(args.config, kind=BIOFILM_CYLINDER)
    dose_cfg = load_dose_rate_config(dose_config_path, kind=BIOFILM_CYLINDER)
    snapshot = load_snapshot(args.snapshot)
    n = snapshot.cell_id.shape[0]
    checks = Checks()

    print(f"reference system : {cfg.reference_system_id}")
    print(f"execution class  : {cfg.execution_class} "
          f"(evidence_policy={cfg.evidence_policy}, "
          f"target_calibration={cfg.target_calibration})")
    print(f"lattice          : {n}^3 at {cfg.voxel_pitch_cm} cm "
          f"-> {n * cfg.voxel_pitch_cm} cm side")
    print(f"occupied voxels  : {int((snapshot.cell_id > 0).sum())}")

    # --- transport stage: eV/source, no activity read --------------------
    model = build_biofilm_cylinder_model(snapshot, cfg)
    mesh = model.tallies[0].filters[0].mesh
    dim = resolve_mesh_dimension((n, n, n), cfg.mesh_coarsening_factor)
    print(f"\nrunning transport: {cfg.batches} x {cfg.particles} "
          f"on a {tuple(dim)} mesh")
    with run_transport(model, args.outdir, "transport") as sp:
        # --- invariant 5: the mass denominator ---------------------------
        exact_mass = mesh_material_masses_kg(mesh, model, dim,
                                             n_samples=args.volume_samples)
        full_bin = coarsen_field(voxel_mass_kg(snapshot, cfg),
                                 cfg.mesh_coarsening_factor)
        bias = full_bin.sum() / exact_mass.sum()
        checks.record(
            "5a_csg_mass_differs_from_full_bin", bias > 1.0,
            f"full-bin mass overstates the CSG mass by {bias:.4f}x "
            f"({full_bin.sum():.4g} vs {exact_mass.sum():.4g} kg) — the full-bin "
            "denominator cannot see the cylinder wall or the membrane")

        per_source = per_source_from_statepoint(sp, exact_mass)
        total_eV = float(np.asarray(sp.get_tally(name="heating").mean).sum())

    # Factor 1 at the same histories: mass-weighted total energy is
    # coarsening-invariant, or the f^3 trap is back.
    cfg_fine = type(cfg)(**{**cfg.__dict__, "mesh_coarsening_factor": 1})
    model_fine = build_biofilm_cylinder_model(snapshot, cfg_fine)
    mesh_fine = model_fine.tallies[0].filters[0].mesh
    with run_transport(model_fine, args.outdir, "transport_f1") as sp1:
        exact_fine = mesh_material_masses_kg(mesh_fine, model_fine, (n, n, n),
                                             n_samples=args.volume_samples)
        per_source_fine = per_source_from_statepoint(sp1, exact_fine)
        total_eV_fine = float(np.asarray(sp1.get_tally(name="heating").mean).sum())

    energy_ratio = total_eV_fine / total_eV if total_eV else float("inf")
    checks.record(
        "5b_total_energy_is_coarsening_invariant",
        abs(energy_ratio - 1.0) < 0.05,
        f"factor 1 vs {cfg.mesh_coarsening_factor}: total heating ratio "
        f"{energy_ratio:.5f} (same histories, different bin counts)")
    mass_ratio = exact_fine.sum() / exact_mass.sum()
    checks.record(
        "5c_mass_is_coarsening_invariant", abs(mass_ratio - 1.0) < 0.01,
        f"exact CSG mass ratio {mass_ratio:.5f} at {args.volume_samples:.0e} "
        "rays; tolerance sized from the raytrace sampling error, not asserted "
        "as equality")

    # --- invariant 4: the unit chain stays distinguishable ----------------
    dose = dose_from_per_source(per_source, dose_cfg.photons_per_second)
    t_hash = transport_state_hash(snapshot, cfg, "endfb-viii.0")
    # The SAME transport hash from the dosimetry config: the source rate is not
    # part of the transport identity, because heating is per source particle.
    t_hash_dose_cfg = transport_state_hash(snapshot, dose_cfg, "endfb-viii.0")
    d_hash = dose_state_hash(t_hash, dose_cfg.photons_per_second)
    d_hash_other = dose_state_hash(t_hash, dose_cfg.photons_per_second * 2)
    ratio = np.divide(dose.dose_rate_mean_Gy_s, per_source.field,
                      out=np.zeros_like(per_source.field),
                      where=per_source.field > 0)
    checks.record(
        "4_gy_per_source_and_gy_per_second_are_distinct",
        type(per_source) is not type(dose)
        and t_hash == t_hash_dose_cfg and d_hash != d_hash_other
        and np.allclose(ratio[per_source.field > 0], dose_cfg.photons_per_second),
        f"transport hash is activity-independent ({t_hash[:12]}), dose hash is "
        f"not; Gy/s over Gy/source is exactly the source rate "
        f"{dose_cfg.photons_per_second:.3g}")

    # --- invariant 8: labels alone must not force a rerun -----------------
    relabelled = type(snapshot)(**{**snapshot.__dict__,
                                   "label_state_hash": "different-labels",
                                   "lineage_id": snapshot.lineage_id + 100})
    checks.record(
        "8_label_only_change_reuses_the_field",
        not needs_rerun(t_hash, relabelled, cfg, "endfb-viii.0"),
        "renaming lineages leaves the transport identity unchanged — labels "
        "never become materials, so the field is reused and only attribution "
        "is recomputed")

    # --- invariant 6 + attribution ---------------------------------------
    occupied = snapshot.cell_id > 0
    dose_lattice = upsample_field(dose.dose_rate_mean_Gy_s,
                                  cfg.mesh_coarsening_factor)
    sd_lattice = upsample_field(dose.dose_rate_sd_Gy_s,
                                cfg.mesh_coarsening_factor)
    lattice_mass = voxel_mass_kg(snapshot, cfg)

    aggregates = {}
    for kind, labels in (("cell", snapshot.cell_id),
                         ("lineage", snapshot.lineage_id),
                         ("generation", snapshot.generation)):
        aggregates[kind] = aggregate_by_label(dose_lattice, labels, lattice_mass,
                                              occupied=occupied)

    voxel_energy = float((dose_lattice[occupied] * lattice_mass[occupied]).sum())
    for kind, table in aggregates.items():
        label_energy = sum(v["dose_rate_Gy_s"] * v["mass_kg"] for v in table.values())
        checks.record(
            f"6_{kind}_attribution_conserves_energy",
            math.isclose(label_energy, voxel_energy, rel_tol=1e-9),
            f"{len(table)} labels, sum(D*m) = {label_energy:.6g} vs "
            f"{voxel_energy:.6g} over occupied voxels")

    generations = sorted(aggregates["generation"])
    checks.record(
        "6d_founder_generation_zero_is_present", 0 in generations,
        f"generations attributed: {generations} — every parcel in this snapshot "
        "is a founder, so a `labels > 0` fallback would have returned nothing")

    # --- invariant 7: dose falls off with distance from the source --------
    pos = source_position_cm(cfg, *( [cfg.origin_cm[0] + n * cfg.voxel_pitch_cm / 2.0,
                                      cfg.origin_cm[1] + n * cfg.voxel_pitch_cm / 2.0,
                                      cfg.origin_cm[2],
                                      cfg.origin_cm[2] + cfg.cylinder_length_cm] ))
    idx = np.indices(snapshot.cell_id.shape).astype(float)
    centres = [(idx[a] + 0.5) * cfg.voxel_pitch_cm + cfg.origin_cm[a] for a in range(3)]
    # A POINT source organizes the field by 3-D distance to the point, not by
    # cylindrical radius: two parcels at equal radius but different z do not
    # see the same fluence.
    distance = np.sqrt(sum((c - p) ** 2 for c, p in zip(centres, pos)))

    d_occ, dist_occ = dose_lattice[occupied], distance[occupied]
    sd_occ = sd_lattice[occupied]
    median = np.median(dist_occ)
    near, far = dist_occ <= median, dist_occ > median
    dn, df = d_occ[near].mean(), d_occ[far].mean()
    # Standard error of each group mean, from the per-voxel sd.
    sn = np.sqrt((sd_occ[near] ** 2).sum()) / near.sum()
    sf = np.sqrt((sd_occ[far] ** 2).sum()) / far.sum()
    separation = (dn - df) / math.sqrt(sn ** 2 + sf ** 2)
    checks.record(
        "7_dose_falls_off_with_distance_to_the_source",
        separation > ORDERING_K,
        f"near {dn:.4g} vs far {df:.4g} Gy/s, separated by {separation:.1f} "
        f"sigma (declared k={ORDERING_K}); grouped by 3-D distance to the "
        "point source, not cylindrical radius")

    # --- invariant 3 + 9: write the field and the provenance --------------
    prov = provenance_attrs(cfg, config_sha256=sha256(args.config),
                            snapshot_sha256=sha256(args.snapshot),
                            repo_root=str(REPO))
    extra = {
        "mesh_dimension": np.array(dim, dtype=np.int64),
        "mesh_coarsening_factor": cfg.mesh_coarsening_factor,
        "histories": cfg.batches * cfg.particles,
        "dose_state_hash": d_hash,
        "label_state_hash": snapshot.label_state_hash,
        "material_volume_n_samples": args.volume_samples,
        "full_bin_mass_bias": bias,
        "stage": "dosimetry",
    }
    result_path = args.outdir / "transport_result.h5"
    write_transport_result(result_path, dose, transport_hash=t_hash,
                           openmc_version=openmc.__version__,
                           nuclear_data_id="endfb-viii.0",
                           batches=cfg.batches, particles=cfg.particles,
                           seed=cfg.seed, provenance=prov, extra=extra)

    field_path = args.outdir / "dose_field.h5"
    write_dose_field(field_path, dose_lattice, sd_lattice,
                     {**{k: v for k, v in prov.items()},
                      "transport_state_hash": t_hash,
                      "source_rate_photons_per_s": dose.source_rate})
    write_dose_attribution(args.outdir / "dose_attribution.h5", aggregates,
                           transport_hash=t_hash,
                           label_hash=snapshot.label_state_hash,
                           uncertainty_method="quadrature_independent_bins_approx")

    # THE COUPLING THIS REPOSITORY REPORTS HAD NO BUNDLE, and therefore none of
    # the co-registration checking `viewer.bundle_problems` applies: grids in
    # one bundle must describe the same physical volume and share one origin.
    # `subvoxel_refinement.py` gets that for free and this path did not, so the
    # dose and the labels were only ever registered by index arithmetic.
    #
    # Costs no transport: everything below is already computed above. Pure
    # numpy + h5py, so the pinned env is untouched by this.
    _write_viewer_bundle(args, cfg, snapshot, dose, dose_lattice, sd_lattice,
                         lattice_mass, occupied, n)

    import h5py
    with h5py.File(result_path) as f:
        written = dict(f.attrs)
    required = {"reference_system_id", "system_provenance", "evidence_policy",
                "execution_class", "target_calibration", "config_sha256",
                "snapshot_sha256", "git_commit", "transport_state_hash",
                "dose_state_hash", "label_state_hash", "openmc_version",
                "nuclear_data", "batches", "particles", "seed",
                "mesh_coarsening_factor", "material_volume_n_samples"}
    missing = sorted(required - set(written))
    checks.record("9_result_carries_its_own_provenance", not missing,
                  f"{len(required)} required attrs present"
                  if not missing else f"missing: {missing}")

    checks.record(
        "3_dose_field_is_written_at_lattice_resolution",
        dose_lattice.shape == snapshot.cell_id.shape,
        f"upsampled {tuple(dim)} -> {dose_lattice.shape}; import_dose_field! "
        "enforces shape equality with the CPM lattice")

    sparsity = sparsity_report(dose, occupied[::cfg.mesh_coarsening_factor,
                                              ::cfg.mesh_coarsening_factor,
                                              ::cfg.mesh_coarsening_factor])
    print(f"\nsparsity on occupied bins: {sparsity}")

    # --- the verdict, as fields ------------------------------------------
    verdict = {
        "execution_status": "PASSED" if checks.all_passed else "FAILED",
        "reference_system_id": cfg.reference_system_id,
        "execution_class": cfg.execution_class,
        "evidence_policy": cfg.evidence_policy,
        "system_provenance": cfg.system_provenance,
        "target_calibration": bool(cfg.target_calibration),
        "claim_class": "integration_path_validation",
        "biological_calibration": "NOT_EVALUATED",
        "production_mesh_selected": False,
        "seconds_per_mcs_calibrated": False,
        "hamiltonian_applied": False,
        "melanin_applied": False,
        "membrane_response_applied": False,
        "cpm_time_advanced": False,
        "invariants": checks.results,
        "provenance": {k: (v.tolist() if isinstance(v, np.ndarray) else
                           (bool(v) if isinstance(v, np.bool_) else
                            (v.item() if isinstance(v, np.generic) else v)))
                       for k, v in written.items()},
    }
    out = args.verdict or (args.outdir / "verdict.json")
    out.write_text(json.dumps(verdict, indent=2), encoding="utf-8")
    print(f"\nverdict: {verdict['execution_status']} -> {out}")
    print("this run validated the INTEGRATION PATH. It calibrated nothing.")
    return 0 if checks.all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
