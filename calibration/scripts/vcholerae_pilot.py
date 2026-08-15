#!/usr/bin/env python
"""Public-data pilot: run a real 3-D biofilm stack through the spatial harness.

Reference system: `public_vcholerae_surrogate`. Vibrio cholerae is NOT one of
the seven modelled species, so a successful run licenses exactly one claim —

    "the spatial calibration machinery was validated against a named public
     3-D biofilm reference system"

— and never "the biological pitch was calibrated". The script writes
`target_calibration = false` into its own output so the distinction survives
being copied out of context.

WHICH FILE. Preflight on 2026-08-15 established that Dryad
doi:10.5061/dryad.zcrjdfnph does not contain what its ledger row assumed:

    Dataset_1  T=47, Y, X, S=3    UNCALIBRATED RGB brightfield, 2-D
    Dataset_2  T=49, C=3, Y, X    calibrated, 2-D  (a time course, not a stack)
    Dataset_3  Z=94, C=3, Y, X    calibrated 3-D stack, NO time axis

No file is both 3-D and time-resolved, so the spatial and temporal passes
cannot be run on the same file, and the temporal observable cannot be run at
all — a biomass field is 3-D by contract. This pilot therefore runs the
SPATIAL pass against Dataset_3, and the temporal half stays blocked on
acquisition rather than being fabricated from a single plane.

Usage (needs the optional dataset extra and a local file):

    pip install -e "calibration[dataset]"
    python calibration/scripts/vcholerae_pilot.py --preflight \
        --nd2 ~/biofilm-pilot-data/Dataset_3_Fig3c_VPS-Bap1.nd2
    python calibration/scripts/vcholerae_pilot.py \
        --nd2 ~/biofilm-pilot-data/Dataset_3_Fig3c_VPS-Bap1.nd2 \
        --segmentation-basis cells_only --channel 0 \
        --out data/calibration/spatial/pilot_vcholerae.csv

OBTAINING THE FILE. The dataset is CC0 and freely downloadable, but the site is
behind a proof-of-work anti-bot challenge that exists specifically to stop
automated clients. Fetch it with a browser rather than defeating that check;
the expected SHA-256 of every file is in
data/calibration/spatial/dataset_candidates.sha256 and is verified here before
the file is opened.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from biofilm_calibration.acquisition import SEGMENTATION_BASIS
from biofilm_calibration.spatial.pitch_selection import (OBSERVABLE_TOLERANCES,
                                                         evaluate_pitch)
from biofilm_calibration.spatial.structure import correlation_length

REFERENCE_SYSTEM = "public_vcholerae_surrogate"
SOURCE_ID = "DRYAD_ZCRJDFNPH"
MANIFEST = (Path(__file__).resolve().parents[2] / "data" / "calibration" /
            "spatial" / "dataset_candidates.sha256")

EXIT_NO_FILE = 2
EXIT_DIGEST = 3


def expected_digest(manifest: Path, name: str) -> str | None:
    """The declared SHA-256 for one filename, from the sha256sum manifest."""
    if not manifest.exists():
        return None
    for line in manifest.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, _, fname = line.partition("  ")
        if fname.strip() == name:
            return digest.strip()
    return None


def verify_digest(path: Path, manifest: Path) -> str:
    """Refuse to open a file whose bytes are not the declared bytes."""
    want = expected_digest(manifest, path.name)
    if want is None:
        print(f"{path.name} is not declared in {manifest}. A file with no "
              "declared digest has no chain of custody; add it to the manifest "
              "or pass --no-verify-digest and accept that the output will say "
              "so.", file=sys.stderr)
        raise SystemExit(EXIT_DIGEST)
    with open(path, "rb") as fh:
        got = hashlib.file_digest(fh, "sha256").hexdigest()
    if got != want:
        print(f"DIGEST MISMATCH for {path}\n  expected {want}\n  got      {got}",
              file=sys.stderr)
        raise SystemExit(EXIT_DIGEST)
    return got


def git_provenance() -> dict:
    """The commit this ran at, and whether the tree was dirty.

    A dirty tree recorded as a clean commit is a false provenance claim, so the
    marker is part of the string rather than a separate field somebody can drop.
    """
    root = Path(__file__).resolve().parents[2]
    try:
        head = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True,
                              ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": None, "git_dirty": None}
    return {"git_commit": head + ("-dirty" if dirty else ""), "git_dirty": dirty}


def spatial_pass(field, factors, *, occupancy_mapping, tau, seed, tolerances,
                 transform, holdout_kind="none") -> list[dict]:
    """Coarse-graining loss curve for one stack."""
    rows = []
    for f in factors:
        if any(d % f for d in field.values.shape):
            continue                       # nesting is required, not optional
        # held_out=False ALWAYS. PitchEvaluation.held_out is the flag
        # select_pitch reads as "validated on an independent stack"; nothing
        # this pilot produces is one, so an evaluation from here must be
        # refused by select_pitch rather than silently satisfying it.
        ev = evaluate_pitch(field.values, field.voxel_size_um, f,
                            sample_id=field.sample_id, held_out=False,
                            occupancy_mapping=occupancy_mapping,
                            tolerances=tolerances, tau=tau, seed=seed)
        row = {"reference_system_id": REFERENCE_SYSTEM,
               "target_calibration": "false",
               "sample_id": field.sample_id,
               "segmentation_basis": field.segmentation_basis,
               "segmentation_method": transform[0],
               "temporal_calibration_eligible": "false",
               "holdout_kind": holdout_kind,
               "coarsening_factor": f,
               "pitch_um": round(ev.pitch_um, 6),
               "passed": str(ev.passed).lower()}
        row.update({f"err_{k}": round(v, 6) for k, v in ev.errors.items()})
        row.update({f"coarse_{k}": round(v, 6) if isinstance(v, float) else v
                    for k, v in ev.coarse.items()})
        rows.append(row)
    return rows


def _n(x):
    """NaN is not JSON. A number that does not exist is null."""
    return None if x is None or not np.isfinite(x) else float(x)


def _jsonable(obj):
    """Recursively replace non-finite floats with null.

    json.dumps emits bare `NaN` by default, which no strict parser accepts. The
    sidecar is written with allow_nan=False so that a future non-finite value
    fails loudly here rather than producing a file that only some readers can
    load.
    """
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float):
        return _n(obj)
    return obj


def temporal_pass(fields, voxel_um, *, transform=None) -> dict:
    """The DECLARED dynamic observable, plus the one that was substituted for it.

    time_observable.csv declares `biomass_autocorrelation_decay` as the
    correlation length of the biomass field tracked across frames — that is the
    SELECTED observable and it is the primary result here.

    A temporal autocorrelation C_B(dt) was written into that row's notes on
    2026-08-14 without the declaration changing, and the first version of this
    script implemented the note. It is reported as a secondary, explicitly
    unselected candidate (`biomass_temporal_autocorrelation_decay`) so it can
    inform a future comparison without standing in for the declaration.
    """
    method = transform[0] if transform else None
    eligible = method is not None and method != "smoke_minmax_per_frame"
    stamp = {"reference_system_id": REFERENCE_SYSTEM,
             "target_calibration": "false",
             "segmentation_method": method,
             "temporal_calibration_eligible": str(bool(eligible)).lower(),
             "seconds_per_mcs_eligible": "false",
             "n_timepoints": len(fields)}

    # --- the selected observable ------------------------------------------
    xi = [_n(correlation_length(f.values, voxel_um)) for f in fields]
    stamp["selected_observable"] = "biomass_autocorrelation_decay"
    stamp["spatial_correlation_length_um_by_t"] = [
        None if v is None else round(v, 4) for v in xi]

    # --- the substituted one, kept apart and labelled ----------------------
    stack = np.stack([f.values for f in fields])          # (T, x, y, z)
    dev = stack - stack.mean(axis=0)
    raw = [float(np.mean([(dev[t] * dev[t + dt]).mean()
                          for t in range(len(fields) - dt)]))
           for dt in range(len(fields))]
    if raw[0] <= 0:
        stamp.update({"temporal_observable_computable": False,
                      "reason": "no temporal fluctuation"})
        return stamp

    # normalised by C(0) computed the SAME WAY — averaged over all available
    # pairs. Dividing by the variance at t=0 alone lets C(0) exceed 1 whenever
    # the field's variance grows, which is what a developing biofilm does.
    corr = [c / raw[0] for c in raw]
    below = [i for i, c in enumerate(corr) if c < 1 / np.e]
    stamp.update({
        "temporal_observable_computable": True,
        "secondary_observable": "biomass_temporal_autocorrelation_decay",
        "secondary_observable_selected": "false",
        "lags_frames": list(range(len(fields))),
        "correlation": [round(c, 5) for c in corr],
        "tau_e_folding_frames": _n(float(below[0]) if below else float("nan")),
        "reason": None if eligible else
                  "per-frame rescaling is not comparable across frames",
    })
    return stamp


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nd2", required=True)
    ap.add_argument("--out")
    ap.add_argument("--preflight", action="store_true",
                    help="read metadata only, allocate no array, and print what "
                         "the file declares about itself. Run this first.")
    ap.add_argument("--segmentation-basis", choices=sorted(SEGMENTATION_BASIS),
                    help="REQUIRED for a run: what the stained channel "
                         "encloses. An external experimental judgement, never "
                         "read from the file, and no default is honest. A "
                         "single fluorescence channel is usually cells_only — "
                         "and that basis CANNOT produce a density (see the "
                         "volume basis contract). Declare it honestly.")
    ap.add_argument("--channel", type=int, default=None,
                    help="no default: a file with several channels does not "
                         "tell you which one is biomass")
    ap.add_argument("--position", type=int, default=None)
    ap.add_argument("--transform", choices=["fixed_threshold",
                                            "fixed_background_and_scale",
                                            "smoke_minmax_per_frame"],
                    default="smoke_minmax_per_frame")
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--scale-lo", type=float, default=None)
    ap.add_argument("--scale-hi", type=float, default=None)
    ap.add_argument("--occupancy-mapping", default="threshold",
                    choices=["threshold", "mass_preserving"])
    ap.add_argument("--tau", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--factors", default="1,2,4,8")
    ap.add_argument("--tolerance", type=float, default=0.05)
    ap.add_argument("--max-timepoints", type=int, default=8)
    ap.add_argument("--no-verify-digest", action="store_true",
                    help="skip the SHA-256 check. Does not buy silence: the "
                         "output records that the bytes were unverified.")
    args = ap.parse_args(argv)

    path = Path(args.nd2)
    if not path.exists():
        print(f"{path} not found. Dryad doi:10.5061/dryad.zcrjdfnph is CC0 but "
              "is served behind a proof-of-work anti-bot challenge; fetch it "
              "with a browser and verify against "
              "data/calibration/spatial/dataset_candidates.sha256.")
        return EXIT_NO_FILE

    from biofilm_calibration.spatial import readers

    if args.preflight:
        info = readers.preflight(path)
        print(json.dumps(info, indent=2, default=str))
        if not info["is_3d"]:
            print("\nNOT A 3-D STACK — a biomass field is 3-D by contract, so "
                  "this file cannot enter the spatial harness.", file=sys.stderr)
        return 0

    if not args.segmentation_basis:
        ap.error("--segmentation-basis is required for a run "
                 f"(one of {sorted(SEGMENTATION_BASIS)})")
    if not args.out:
        ap.error("--out is required for a run")

    sha = None
    if not args.no_verify_digest:
        sha = verify_digest(path, MANIFEST)      # before the reader is used

    if args.transform == "fixed_threshold":
        if args.threshold is None:
            ap.error("--transform fixed_threshold needs --threshold")
        transform = ("fixed_threshold", args.threshold, float("nan"))
    elif args.transform == "fixed_background_and_scale":
        if args.scale_lo is None or args.scale_hi is None:
            ap.error("--transform fixed_background_and_scale needs --scale-lo "
                     "and --scale-hi")
        transform = (args.transform, args.scale_lo, args.scale_hi)
    else:
        transform = readers.SMOKE

    tolerances = {k: args.tolerance for k in OBSERVABLE_TOLERANCES.values()}
    factors = [int(f) for f in args.factors.split(",")]

    read_kw = dict(segmentation_basis=args.segmentation_basis,
                   source_id=SOURCE_ID, channel=args.channel,
                   position=args.position, transform=transform)
    try:
        fields = readers.read_nd2_timecourse(
            path, sample_id=path.stem, max_timepoints=args.max_timepoints,
            **read_kw)
    except readers.NoTimeAxisError:
        # the ONLY fallback: this file is a single volume. Every other failure
        # — a bad channel, an absent calibration, a corrupt file — propagates
        # with its own diagnosis instead of being retried as a static image.
        fields = [readers.read_nd2(path, sample_id=path.stem, **read_kw)]

    print(f"{len(fields)} volume(s), shape {fields[0].shape}, "
          f"voxel {fields[0].voxel_size_um} um, "
          f"basis {fields[0].segmentation_basis}, transform {transform[0]}")

    rows = []
    for i, f in enumerate(fields):
        # The last frame of one file is the SAME field of view and the same
        # biofilm, correlated with every earlier frame. It can show whether a
        # coarsening rule stays stable in time; it cannot show generalisation
        # to another sample, so it is named for what it is.
        kind = ("temporal_last_frame"
                if len(fields) > 1 and i == len(fields) - 1 else "none")
        rows += spatial_pass(f, factors, occupancy_mapping=args.occupancy_mapping,
                             tau=args.tau, seed=args.seed, tolerances=tolerances,
                             transform=transform, holdout_kind=kind)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        fh.write(f"# Public-data pilot: {REFERENCE_SYSTEM}\n"
                 "# target_calibration = false. Vibrio cholerae is NOT one of the\n"
                 "# seven modelled species, so this validates the MACHINERY against a\n"
                 "# named public reference system and calibrates nothing.\n"
                 "# holdout_kind is never an independent sample: no frame of one\n"
                 "# file is a held-out stack, so select_pitch refuses these rows.\n"
                 f"# source: Dryad doi:10.5061/dryad.zcrjdfnph (CC0), file {path.name}\n")
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    print(f"wrote {len(rows)} rows to {out}")

    temporal = None
    if len(fields) > 1:
        temporal = temporal_pass(fields, fields[0].voxel_size_um,
                                 transform=transform)
        print("\ntemporal observable:")
        print(json.dumps(temporal, indent=2))

    meta = dict(fields[0].meta)
    sidecar = {
        "reference_system_id": REFERENCE_SYSTEM,
        "target_calibration": False,
        "source_id": SOURCE_ID,
        "identifier": "doi:10.5061/dryad.zcrjdfnph",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": {"path": str(path), "name": path.name,
                  "size_bytes": path.stat().st_size, "sha256": sha,
                  "sha256_verified": sha is not None,
                  "manifest": str(MANIFEST) if sha else None},
        "reader": meta,
        "declaration": {"segmentation_basis": args.segmentation_basis,
                        "segmentation_basis_source": "operator_declared_cli",
                        "segmentation_method": transform[0],
                        "transform": list(transform)},
        "parameters": {"occupancy_mapping": args.occupancy_mapping,
                       "tau": args.tau, "seed": args.seed, "factors": factors,
                       "tolerance": args.tolerance, "tolerances": tolerances,
                       "max_timepoints": args.max_timepoints},
        "temporal": temporal,
        "provenance": {**git_provenance(), "python": sys.version.split()[0],
                       "numpy": np.__version__},
        "outputs": {"spatial_csv": str(out), "n_rows": len(rows)},
    }
    side = out.with_suffix(".metadata.json")
    side.write_text(json.dumps(_jsonable(sidecar), indent=2, default=str,
                               allow_nan=False) + "\n")
    print(f"wrote provenance to {side}")

    print(f"\nreference_system_id = {REFERENCE_SYSTEM}; target_calibration = false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
