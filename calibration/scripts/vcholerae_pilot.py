#!/usr/bin/env python
"""Public-data pilot: run a real 3-D biofilm stack through the spatial harness.

Reference system: `public_vcholerae_surrogate`. Vibrio cholerae is NOT one of
the seven modelled species, so a successful run licenses exactly one claim —

    "the spatial and temporal calibration machinery was validated against a
     named public 3-D biofilm reference system"

— and never "the biological pitch was calibrated". The script writes
`target_calibration = false` into its own output so the distinction survives
being copied out of context.

Usage (needs the optional dataset extra and a local file):

    pip install -e "calibration[dataset]"
    python calibration/scripts/vcholerae_pilot.py \
        --nd2 /path/to/Dataset_2_Fig4a_Timecourse.nd2 \
        --out data/calibration/spatial/pilot_vcholerae.csv

OBTAINING THE FILE. Dryad doi:10.5061/dryad.zcrjdfnph is CC0 and freely
downloadable, but the site is behind a proof-of-work anti-bot challenge that
exists specifically to stop automated clients. Fetch it with a browser rather
than defeating that check; the expected SHA-256 of each file is recorded in
data/calibration/spatial/dataset_candidates.csv so the copy can be verified.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from biofilm_calibration.spatial.field import apply_occupancy, coarse_grain
from biofilm_calibration.spatial.pitch_selection import (OBSERVABLE_TOLERANCES,
                                                         evaluate_pitch)
from biofilm_calibration.spatial.structure import correlation_length, summarise

REFERENCE_SYSTEM = "public_vcholerae_surrogate"


def spatial_pass(field, factors, *, occupancy_mapping, tau, seed, tolerances,
                 held_out=False) -> list[dict]:
    """Coarse-graining loss curve for one stack."""
    rows = []
    for f in factors:
        if any(d % f for d in field.values.shape):
            continue                       # nesting is required, not optional
        ev = evaluate_pitch(field.values, field.voxel_size_um, f,
                            sample_id=field.sample_id, held_out=held_out,
                            occupancy_mapping=occupancy_mapping,
                            tolerances=tolerances, tau=tau, seed=seed)
        row = {"reference_system_id": REFERENCE_SYSTEM,
               "target_calibration": "false",
               "sample_id": field.sample_id,
               "segmentation_basis": field.segmentation_basis,
               "held_out": str(held_out).lower(),
               "coarsening_factor": f,
               "pitch_um": round(ev.pitch_um, 6),
               "passed": str(ev.passed).lower()}
        row.update({f"err_{k}": round(v, 6) for k, v in ev.errors.items()})
        row.update({f"coarse_{k}": round(v, 6) if isinstance(v, float) else v
                    for k, v in ev.coarse.items()})
        rows.append(row)
    return rows


def temporal_pass(fields, voxel_um) -> dict:
    """Decay of the biomass-field autocorrelation — the SELECTED observable.

    Field-level on both sides, so no cell-to-parcel mapping is needed. This is
    the feasibility half of the seconds/MCS blocker: it says whether the
    statistic is computable and well-behaved on real data, not what dt_MCS is.
    """
    stack = np.stack([f.values for f in fields])          # (T, x, y, z)
    mean = stack.mean(axis=0)
    dev = stack - mean

    # C(dt) normalised by C(0) computed the SAME WAY — averaged over all
    # available pairs. Dividing instead by the variance at t=0 alone would let
    # C(0) exceed 1 whenever the field's variance grows over the series, which
    # it does whenever a biofilm is developing.
    lags, raw = [], []
    for dt in range(len(fields)):
        pairs = [(dev[t] * dev[t + dt]).mean() for t in range(len(fields) - dt)]
        lags.append(dt)
        raw.append(float(np.mean(pairs)))
    if raw[0] <= 0:
        return {"computable": False, "reason": "no temporal fluctuation"}
    corr = [c / raw[0] for c in raw]

    below = [i for i, c in enumerate(corr) if c < 1 / np.e]
    tau_frames = float(below[0]) if below else float("nan")
    return {"computable": True, "n_timepoints": len(fields),
            "lags_frames": lags, "correlation": [round(c, 5) for c in corr],
            "tau_e_folding_frames": tau_frames,
            "spatial_correlation_length_um_t0":
                round(correlation_length(fields[0].values, voxel_um), 4)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nd2", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--segmentation-basis", default="cells_only",
                    help="what the stained channel encloses. Default cells_only "
                         "because a single fluorescence channel usually is — and "
                         "that basis CANNOT produce a density (see the volume "
                         "basis contract). Declare it honestly.")
    ap.add_argument("--channel", type=int, default=0)
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--occupancy-mapping", default="threshold",
                    choices=["threshold", "mass_preserving"])
    ap.add_argument("--tau", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--factors", default="1,2,4,8")
    ap.add_argument("--tolerance", type=float, default=0.05)
    ap.add_argument("--max-timepoints", type=int, default=8)
    args = ap.parse_args(argv)

    path = Path(args.nd2)
    if not path.exists():
        print(f"{path} not found. Dryad doi:10.5061/dryad.zcrjdfnph is CC0 but "
              "is served behind a proof-of-work anti-bot challenge; fetch it "
              "with a browser and verify the SHA-256 recorded in "
              "data/calibration/spatial/dataset_candidates.csv.")
        return 2

    from biofilm_calibration.spatial.readers import read_nd2, read_nd2_timecourse

    tolerances = {k: args.tolerance for k in OBSERVABLE_TOLERANCES.values()}
    factors = [int(f) for f in args.factors.split(",")]

    try:
        fields = read_nd2_timecourse(
            path, segmentation_basis=args.segmentation_basis,
            sample_id=path.stem, source_id="DRYAD_ZCRJDFNPH",
            channel=args.channel, threshold=args.threshold)
        fields = fields[:args.max_timepoints]
    except Exception:                       # not a time course
        fields = [read_nd2(path, segmentation_basis=args.segmentation_basis,
                           sample_id=path.stem, source_id="DRYAD_ZCRJDFNPH",
                           channel=args.channel, threshold=args.threshold)]

    print(f"{len(fields)} volume(s), shape {fields[0].shape}, "
          f"voxel {fields[0].voxel_size_um} um, "
          f"basis {fields[0].segmentation_basis}")

    rows = []
    for i, f in enumerate(fields):
        # the LAST timepoint is held out, so the loss curve is validated on a
        # stack it was not chosen on
        rows += spatial_pass(f, factors, occupancy_mapping=args.occupancy_mapping,
                             tau=args.tau, seed=args.seed, tolerances=tolerances,
                             held_out=(i == len(fields) - 1))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        fh.write(f"# Public-data pilot: {REFERENCE_SYSTEM}\n"
                 "# target_calibration = false. Vibrio cholerae is NOT one of the\n"
                 "# seven modelled species, so this validates the MACHINERY against a\n"
                 "# named public reference system and calibrates nothing.\n"
                 f"# source: Dryad doi:10.5061/dryad.zcrjdfnph (CC0), file {path.name}\n")
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    print(f"wrote {len(rows)} rows to {out}")

    if len(fields) > 1:
        temporal = temporal_pass(fields, fields[0].voxel_size_um)
        print("\ntemporal observable (biomass_autocorrelation_decay):")
        print(json.dumps(temporal, indent=2))

    print(f"\nreference_system_id = {REFERENCE_SYSTEM}; target_calibration = false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
