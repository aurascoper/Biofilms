#!/usr/bin/env python
"""Detectability pilot: can this apparatus resolve the measurement at all?

    python calibration/scripts/detectability_pilot.py --demo-refusals
    python calibration/scripts/detectability_pilot.py \
        --data data/calibration/materials \
        --acceptance config/reference_d_material_acceptance.toml \
        --out data/calibration/materials/pilot_detectability.csv

THE ORDER IS THE DESIGN, and every step can refuse:

    support commensurability  is the volume the same material as the mass?
    segmentation basis        does the volume enclose what the balance weighed?
    blank correction          the biofilm mass is a DIFFERENCE
    detectability             is that difference above the blank's noise?
    invalid-draw accounting   does the distribution cross zero?
    -------------------------------------------------------- only then
    density / water fraction
    replicate sizing

A mass buried under the blank's noise is a SUBSTRATE or BIOMASS problem, not a
replicate-count problem. Replicates shrink the standard error of a mean; they do
not recover a signal lost to a subtraction. So this script refuses to size a
campaign rather than return a number that implies more coupons would help. That
refusal is the whole point of running it first.

Reference D runs read their parameters from the FROZEN acceptance config rather
than from flags, so a pilot cannot drift from the declaration it is testing. The
effective values are recorded in the sidecar either way.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from physical_contract import git_provenance as _git_provenance

from biofilm_calibration.materials.basis_conversion import (BasisError,
                                                            HydratedVolume,
                                                            blank_corrected_mass,
                                                            dry_biomass_concentration,
                                                            water_mass_fraction,
                                                            wet_bulk_density)
from biofilm_calibration.materials.schema import BULK_MEASUREMENTS
from biofilm_calibration.materials.uncertainty import propagate
from biofilm_calibration.acquisition import BLANKS, WHOLE_BIOFILM_BASIS
from biofilm_calibration.precision import (PrecisionError, detectability,
                                           required_replicates)
from biofilm_calibration.schema import read_table

REPO = Path(__file__).resolve().parents[2]
REFERENCE_SYSTEM = "D_engineered_composite"

EXIT_NO_DATA = 2
EXIT_REFUSED = 3

# How much of the uncertainty distribution may sit at or below zero before the
# blank-corrected mass stops being a measurement. A DECLARED engineering target,
# in the same spirit as precision.py's minimum_ratio = 10.0, not a law.
MAX_INVALID_DRAW_FRACTION = 0.01

# What the number in `blank_uncertainty_g` actually is. Wet, dry and blank
# masses weighed on one balance share its calibration, so a single sigma cannot
# be interpreted without saying which of these it represents.
BLANK_UNCERTAINTY_KINDS = ("matched_blank_sd", "blank_mean_se",
                           "balance_calibration", "batch_variance",
                           "combined_model")


def git_provenance() -> dict:
    return _git_provenance(str(REPO))


def _n(x):
    return None if x is None or not np.isfinite(x) else float(x)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float):
        return _n(obj)
    return obj


def support_problems(row: dict) -> list[str]:
    """Is the hydrated volume commensurate with the weighed mass?

    rho_wet = m/V is only a density when the numerator and the denominator
    describe the same material. `volume_basis` says what the volume ENCLOSES;
    `volume_support` says what it COVERS. A perfectly segmented whole-biofilm
    envelope of one confocal field is still not the volume of the coupon the
    balance weighed.
    """
    problems = []
    support = (row.get("volume_support") or "").strip()
    if support in ("", "unresolved"):
        problems.append(
            "volume_support is unresolved — a local field of view cannot be "
            "divided into a whole-coupon mass, and nothing else in the row says "
            "whether it was")
    elif support in ("stereological_scaling", "sibling_batch_model"):
        if not (row.get("scaling_method") or "").strip():
            problems.append(f"volume_support={support} with no scaling_method — "
                            "an undeclared extrapolation is not a measurement")
        if row.get("scaling_uncertainty") is None:
            problems.append(f"volume_support={support} with no "
                            "scaling_uncertainty — a scaled volume carries the "
                            "scaling's error into rho_wet")
    if support == "stereological_scaling":
        imaged, weighed = row.get("imaged_area_cm2"), row.get("weighed_area_cm2")
        if imaged is None or weighed is None:
            problems.append("stereological scaling needs both imaged_area_cm2 "
                            "and weighed_area_cm2 to state what was scaled")
    return problems


def basis_problems(row: dict) -> list[str]:
    """Does the volume enclose what the balance weighed?

    Delegates the judgement to HydratedVolume so the refusal text stays in one
    place, and reports rather than raises so a pilot can list every problem at
    once instead of stopping at the first coupon.
    """
    basis = (row.get("volume_basis") or "").strip()
    volume = row.get("hydrated_volume_cm3")
    if volume is None:
        return ["hydrated_volume_cm3 is unset"]
    try:
        hv = HydratedVolume(float(volume), basis or "unresolved")
        if basis != WHOLE_BIOFILM_BASIS and row.get("pore_volume_fraction") is None:
            hv.require_whole_biofilm()
    except BasisError as exc:
        return [str(exc).splitlines()[0]]
    except ValueError as exc:
        return [str(exc)]
    return []


def evaluate_row(row: dict, blanks: dict, *, minimum_ratio: float,
                 seed: int, draws: int) -> dict:
    """One coupon, through the ordered refusals."""
    out = {
        "sample_id": row["sample_id"],
        "replicate_id": row["replicate_id"],
        "reference_system_id": REFERENCE_SYSTEM,
        "target_calibration": "true",
        "support_status": "PASS",
        "basis_status": "PASS",
        "detectability_status": "NOT_EVALUATED",
        "distribution_status": "NOT_EVALUATED",
        "refusals": "",
    }
    refusals = support_problems(row)
    if refusals:
        out["support_status"] = "REFUSED"
    else:
        basis = basis_problems(row)
        if basis:
            out["basis_status"] = "REFUSED"
            refusals += basis

    blank_id = (row.get("blank_sample_id") or "").strip()
    blank = blanks.get(blank_id)
    if blank is None:
        refusals.append(f"blank_sample_id {blank_id!r} does not resolve in "
                        "blanks.csv — blank subtraction defines the biofilm mass")

    if refusals:
        out["refusals"] = " | ".join(refusals)
        return out

    wet = float(row["wet_mass_sample_plus_substrate_g"])
    blank_wet = blank.get("wet_mass_g")
    sigma_blank = blank.get("wet_mass_uncertainty_g")
    try:
        biofilm_wet = blank_corrected_mass(wet, blank_wet, label="wet")
        det = detectability(biofilm_wet, float(sigma_blank or 0.0),
                            minimum_ratio=minimum_ratio)
    except (BasisError, PrecisionError) as exc:
        out["detectability_status"] = "REFUSED"
        out["refusals"] = str(exc).splitlines()[0]
        return out

    out["biofilm_wet_mass_g"] = round(biofilm_wet, 9)
    out["signal_to_blank_noise"] = round(det.signal_to_blank_noise, 4)
    out["detectability_status"] = "PASS" if det.detectable else "FAIL"
    if not det.detectable:
        out["refusals"] = det.reason.splitlines()[0]
        return out

    # A draw where the blank exceeds its sample is finite and NEGATIVE. It would
    # leave here as a negative density if nobody counted it, so it is counted.
    summary = propagate(
        lambda sample, blank_mass: sample - blank_mass,
        {"sample": (wet, row.get("wet_mass_uncertainty_g") or 0.0),
         "blank_mass": (float(blank_wet), float(sigma_blank or 0.0))},
        seed=seed, draws=draws, invalid_if=lambda v: v <= 0.0)
    invalid = summary.get("invalid_fraction", float("nan"))
    out["invalid_draw_fraction"] = round(float(invalid), 6)
    out["negative_blank_corrected_mass_probability"] = round(float(invalid), 6)
    if not np.isfinite(invalid) or invalid > MAX_INVALID_DRAW_FRACTION:
        out["distribution_status"] = "NOT_IDENTIFIABLE"
        out["refusals"] = (
            f"{invalid:.1%} of the uncertainty distribution puts the "
            "blank-corrected mass at or below zero — the measurement does not "
            "distinguish the biofilm from its blank, and replicates cannot "
            "rescue a subtraction. This is a substrate or biomass problem")
        return out

    out["distribution_status"] = "VALID"
    out["wet_mass_rel_uncertainty"] = round(
        summary["sd"] / abs(summary["mean"]), 6) if summary["mean"] else None
    return out


def size_replicates(rows: list[dict], *, target_rel_uncertainty: float,
                    n_within: int) -> dict | None:
    """Replicate sizing, attempted ONLY on coupons that cleared every refusal.

    Returns None when nothing cleared, which is a real answer: a campaign
    sized from undetectable pilots would be a campaign designed to fail.
    """
    usable = [r for r in rows if r.get("distribution_status") == "VALID"]
    if len(usable) < 2:
        return None
    masses = np.array([r["biofilm_wet_mass_g"] for r in usable], dtype=float)
    batches = {}
    for r, m in zip(usable, masses):
        batches.setdefault(r.get("culture_batch_id") or r["sample_id"], []).append(m)

    within = [np.var(v, ddof=1) for v in batches.values() if len(v) > 1]
    between = (np.var([np.mean(v) for v in batches.values()], ddof=1)
               if len(batches) > 1 else 0.0)
    plan = required_replicates(
        target_rel_uncertainty=target_rel_uncertainty,
        mean_value=float(masses.mean()),
        var_between_batch=float(between),
        var_within_batch=float(np.mean(within)) if within else 0.0,
        n_within=n_within)
    return {
        "required_batches": plan.required_batches,
        "achieved_rel_uncertainty": plan.achieved_rel_uncertainty,
        "dominant_component": plan.dominant_component,
        "components": plan.components,
        "note": plan.note,
        "n_usable_coupons": len(usable),
        "n_batches_observed": len(batches),
    }


def demo_refusals() -> int:
    """The four refusals, on invented rows. Exercises the logic without data."""
    print("The four ways this pilot refuses, demonstrated on invented rows.\n")

    print("1. INCOMMENSURATE SUPPORT — a field of view over a coupon mass")
    for p in support_problems({"volume_support": "unresolved"}):
        print(f"   - {p}")

    print("\n2. INCOMPATIBLE BASIS — a cells-only volume under a whole-biofilm mass")
    for p in basis_problems({"volume_basis": "cells_only",
                             "hydrated_volume_cm3": 0.1,
                             "pore_volume_fraction": None}):
        print(f"   - {p}")

    print("\n3. SIGNAL UNDER THE BLANK'S NOISE")
    det = detectability(0.002, 0.001, minimum_ratio=10.0)
    print(f"   - detectable={det.detectable}  ratio={det.signal_to_blank_noise:.2f}")
    print(f"   - {det.reason.splitlines()[0]}")

    print("\n4. THE DISTRIBUTION CROSSES ZERO")
    summary = propagate(lambda sample, blank_mass: sample - blank_mass,
                        {"sample": (1.00, 0.05), "blank_mass": (0.97, 0.05)},
                        seed=1, draws=20000, invalid_if=lambda v: v <= 0.0)
    print(f"   - invalid_draw_fraction={summary['invalid_fraction']:.1%} "
          f"(ceiling {MAX_INVALID_DRAW_FRACTION:.0%})")
    print("   - replicate sizing is NOT attempted while detectability fails")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path,
                    default=REPO / "data" / "calibration" / "materials")
    ap.add_argument("--acceptance", type=Path,
                    help="frozen acceptance config; Reference D runs take "
                         "minimum_replicates and the uncertainty target from it "
                         "rather than from flags")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--minimum-ratio", type=float,
                    help="signal-to-blank-noise target; declared, not a law")
    ap.add_argument("--target-rel-uncertainty", type=float, default=0.10)
    ap.add_argument("--blank-uncertainty-kind", choices=BLANK_UNCERTAINTY_KINDS,
                    help="what the number in blank_uncertainty_g IS. No default: "
                         "wet, dry and blank masses share a balance calibration, "
                         "so an unlabelled sigma cannot be interpreted")
    ap.add_argument("--n-within", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--demo-refusals", action="store_true")
    args = ap.parse_args(argv)

    if args.demo_refusals:
        return demo_refusals()

    acceptance = {}
    if args.acceptance and args.acceptance.exists():
        with open(args.acceptance, "rb") as fh:
            acceptance = tomllib.load(fh).get("acceptance", {})
    minimum_ratio = args.minimum_ratio if args.minimum_ratio is not None else 10.0
    target_rel = acceptance.get("maximum_density_relative_uncertainty",
                                args.target_rel_uncertainty)

    if not args.blank_uncertainty_kind:
        ap.error("--blank-uncertainty-kind is required: a sigma whose meaning is "
                 "undeclared cannot be propagated honestly")

    rows = read_table(args.data / "bulk_measurements.csv", BULK_MEASUREMENTS)
    blank_rows = read_table(args.data / "blanks.csv", BLANKS)
    if not rows:
        print(f"no bulk measurements in {args.data} — nothing to evaluate. The "
              "pilot exists to be run BEFORE the campaign, on real coupons.",
              file=sys.stderr)
        return EXIT_NO_DATA
    blanks = {b["blank_sample_id"]: b for b in blank_rows}

    results = [evaluate_row(r, blanks, minimum_ratio=minimum_ratio,
                            seed=args.seed, draws=args.draws) for r in rows]
    for r in results:
        status = r["detectability_status"]
        print(f"{r['sample_id']}/{r['replicate_id']}: support={r['support_status']} "
              f"basis={r['basis_status']} detect={status} "
              f"dist={r['distribution_status']}")
        if r["refusals"]:
            print(f"    {r['refusals']}")

    plan = size_replicates(results, target_rel_uncertainty=target_rel,
                           n_within=args.n_within)
    if plan is None:
        print("\nREPLICATE SIZING NOT ATTEMPTED — fewer than two coupons cleared "
              "every refusal. Sizing a campaign from undetectable pilots would "
              "design a campaign to fail.")
    else:
        print(f"\nreplicates: {plan['required_batches']} batches, achieving "
              f"{plan['achieved_rel_uncertainty']:.2%} (target {target_rel:.2%})")
        print(f"dominant variance: {plan['dominant_component']} — {plan['note']}")

    if args.out:
        sidecar = {
            "reference_system_id": REFERENCE_SYSTEM,
            "target_calibration": True,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "input": {"data_dir": str(args.data), "n_rows": len(rows)},
            "declaration": {
                "blank_uncertainty_kind": args.blank_uncertainty_kind,
                "blank_uncertainty_kind_source": "operator_declared_cli",
                "acceptance_config": str(args.acceptance) if args.acceptance else None,
                "independence_assumption": (
                    "propagate() samples inputs independently. Wet, dry and "
                    "blank masses weighed on one balance SHARE its calibration; "
                    "any such correlation must be combined before it reaches "
                    "the sampler, not declared independent for convenience."),
            },
            "parameters": {
                "minimum_ratio": minimum_ratio,
                "target_rel_uncertainty": target_rel,
                "max_invalid_draw_fraction": MAX_INVALID_DRAW_FRACTION,
                "n_within": args.n_within, "seed": args.seed, "draws": args.draws,
            },
            "replicate_plan": plan,
            "provenance": {**git_provenance(), "python": sys.version.split()[0],
                           "numpy": np.__version__},
            "outputs": {"csv": str(args.out), "n_rows": len(results)},
        }
        with open(args.out, "w", encoding="utf-8", newline="") as fh:
            fh.write("# Detectability pilot for Reference D.\n"
                     "# It answers ONE question: can this apparatus resolve the\n"
                     "# measurement at all? A PASS licenses replicate sizing and\n"
                     "# nothing else. It is not a density, not a composition, and\n"
                     "# not a calibration of anything.\n")
            w = csv.DictWriter(fh, fieldnames=list(results[0]),
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            w.writerows(results)
        side = args.out.with_suffix(".metadata.json")
        side.write_text(json.dumps(_jsonable(sidecar), indent=2, default=str,
                                   allow_nan=False) + "\n", encoding="utf-8")
        print(f"\nwrote {args.out}\nwrote {side}")

    refused = [r for r in results if r["distribution_status"] != "VALID"]
    return EXIT_REFUSED if len(refused) == len(results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
