#!/usr/bin/env python
"""Run the A0 resolution/history sweep and write data/a0_transport_convergence.csv.

Usage (inside the pinned conda env, with OPENMC_CROSS_SECTIONS set):

    micromamba run -n openmc-biofilms python coupling/scripts/a0_sweep.py \
        --config config/reference_a0_water_phantom.toml \
        --out data/a0_transport_convergence.csv

Every point is per source particle, so no source activity is involved. Fixed
history counts with explicit zero-score accounting, not tally triggers.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import os
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from biofilm_openmc.config import WATER_PHANTOM, load_transport_config
from biofilm_openmc.convergence import (CSV_COLUMNS, field_difference,
                                        in_cylinder_mask, point_metrics,
                                        resolution_loss, row_for,
                                        seed_spread)
from biofilm_openmc.dose import extract_heating
from biofilm_openmc.materials import phantom_mass_kg
from biofilm_openmc.mesh import phantom_mesh_extent_cm, resolve_mesh_dimension
from biofilm_openmc.model import build_water_phantom_model

E_SRC_EV = 661655.0


def _peak_child_rss_kb() -> int:
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss


def run_point(cfg, factor: int, batches: int, particles: int, seed: int,
              workdir: Path):
    """One OpenMC run. Returns (metrics, heating_mean, dimension, timing)."""
    import openmc

    point = dataclasses.replace(cfg, mesh_coarsening_factor=factor,
                                batches=batches, particles=particles, seed=seed)
    dimension = resolve_mesh_dimension(point.mesh_base_dimension, factor)
    model = build_water_phantom_model(point)

    workdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    sp_path = model.run(cwd=str(workdir), output=False)
    runtime = time.time() - t0
    with openmc.StatePoint(sp_path) as sp:
        mean, sd = extract_heating(sp, dimension)
    stat_bytes = os.path.getsize(sp_path)

    mass = phantom_mass_kg(point, dimension, phantom_mesh_extent_cm(point))
    mask = in_cylinder_mask(point, dimension)
    metrics = point_metrics(mean, sd, mass, mask, e_src_eV=E_SRC_EV)
    return point, dimension, mass, mask, mean, metrics, runtime, stat_bytes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--factors", default="1,2,4,8,16")
    ap.add_argument("--histories", default="250000,1000000,4000000,16000000")
    ap.add_argument("--seeds", default="1,2,3,4,5")
    ap.add_argument("--particles-per-batch", type=int, default=50000)
    # The reference must out-resolve every swept point, and must not share a
    # seed with one: an identical (factor, histories, seed) run would reproduce
    # it exactly and report a spurious zero difference.
    ap.add_argument("--reference-histories", type=int, default=0,
                    help="default: 4x the largest swept history count")
    ap.add_argument("--reference-seed", type=int, default=9999)
    args = ap.parse_args(argv)

    xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
    if not xs or not Path(xs).exists():
        print("OPENMC_CROSS_SECTIONS not available — refusing to run")
        return 2

    cfg = load_transport_config(args.config, kind=WATER_PHANTOM)
    factors = [int(f) for f in args.factors.split(",")]
    histories = [int(h) for h in args.histories.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]
    ppb = args.particles_per_batch

    rows: list[dict] = []
    ref_factor = min(factors)
    ref_hist = args.reference_histories or 4 * max(histories)

    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        print(f"reference: factor {ref_factor}, {ref_hist:,} histories, "
              f"seed {args.reference_seed}")
        (_p, ref_dim, ref_mass, ref_mask, ref_heating, ref_met, rt, _sz) = run_point(
            cfg, ref_factor, max(1, ref_hist // ppb), ppb,
            args.reference_seed, root / "reference")
        print(f"  {rt:.1f}s  zero {ref_met['zero_score_fraction_roi']:.2%} "
              f"med_rel {ref_met['median_rel_err']:.4f}")

        for factor in factors:
            for hist in histories:
                batches = max(1, hist // ppb)
                actual = batches * ppb
                for seed in seeds:
                    tag = f"f{factor}_h{actual}_s{seed}"
                    (point, dim, mass, mask, mean, metrics, rt, sz) = run_point(
                        cfg, factor, batches, ppb, seed, root / tag)
                    # compare against the high-statistics reference, aggregating
                    # its deposited energy onto this (coarser or equal) grid
                    fd = field_difference(mean, ref_heating, factor, mass, mask)
                    rl = resolution_loss(mean, mass, factor, ref_heating,
                                         ref_mass, ref_mask)
                    rows.append(row_for(point, factor, dim, seed, metrics,
                                        runtime_s=rt, statepoint_bytes=sz,
                                        peak_rss_kb=_peak_child_rss_kb(),
                                        field_diff=fd, resolution_loss_=rl))
                    print(f"  {tag}: zero {metrics['zero_score_fraction_roi']:.1%} "
                          f"med_rel {metrics['median_rel_err']:.3f} "
                          f"noise {fd:.4f} resloss {rl:.4f} {rt:.1f}s")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        fh.write(
            "# A0 resolution/history convergence sweep — see\n"
            "# docs/a0_transport_resolution_decision.md for the decision and\n"
            "# data/parameter_provenance.csv for what these rows are allowed to\n"
            "# rank (sensitivity_domain = transport_numerical only).\n"
            "#\n"
            "# Geometry, materials and spectrum are HELD FIXED throughout; only\n"
            "# the tally resolution and the history count move. Per source\n"
            "# particle: no source activity is involved.\n"
            "#\n"
            "# zero_score_fraction_roi counts bins INSIDE the cylinder only —\n"
            "# the mesh spans the circumscribing cube, so ~21.5% of bins are\n"
            "# void by construction and must score zero.\n"
            "# field_diff_vs_reference aggregates the reference's DEPOSITED\n"
            "# ENERGY onto each grid before comparing; it is not a mean of\n"
            "# per-bin dose.\n"
            "# peak_child_rss_kb_cumulative is a running high-water mark across\n"
            "# all child processes, so only the largest value is meaningful.\n")
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {len(rows)} rows to {out}")
    top = max(r["histories"] for r in rows)
    print(f"seed spread of the mass-weighted mean at {top:,} histories:")
    for factor in factors:
        vals = [r["mass_weighted_mean_dose_Gy_per_src"] for r in rows
                if r["coarsening_factor"] == factor and r["histories"] == top]
        print(f"  factor {factor}: {seed_spread(vals):.4%}  ({len(vals)} seeds)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
