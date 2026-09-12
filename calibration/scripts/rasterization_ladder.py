#!/usr/bin/env python
"""Phase E: what a lattice pitch costs, measured against a known answer.

    python calibration/scripts/rasterization_ladder.py --outdir artifacts/ladder

THE QUESTION, AND WHY IT IS NOT THE ONE PHASE B ANSWERED. The subvoxel study
asked how finely the TRANSPORT must be resolved on a given snapshot, and found
the effect metric is not resolution-invariant. This asks the separate question
underneath it: how finely must the MORPHOLOGY be sampled before the structural
observables the spatial gate declares are trustworthy? A tally can always be
refined after the fact. A snapshot cannot — its pitch is a measurement decision.

WHY AN ANALYTIC OBJECT AND NOT A CPM RUN. Two reasons, and the second is fatal
to the alternative.

Every builder in `spatial/synthetic.py` takes feature sizes in VOXELS, so
`slab(shape=(16,16,32), thickness_voxels=8)` at a different shape is a different
slab. A sweep built from those changes the object and the sampling together and
cannot separate them. `PhysicalSlab` and `PhysicalSpheres` carry the object in
micrometres and rasterise at whatever pitch is asked, so only the sampling moves.

And a CPM sweep would be worse than uninformative. `V_target` in
biofilms_potts.jl is in SITES, so raising N at fixed V_target shrinks every
parcel: the object changes with the grid, in a way that looks like convergence
and is not. That study needs a discretisation-renormalisation contract first --
contact energy, volume penalty, temperature, sweep rate and field diffusion all
carry a scale -- and none of that is needed here, because an analytic shape has
no dynamics to renormalise.

WHAT IT MEASURES. At each pitch, every declared observable against a truth
computed in closed form, and against the tolerance
`config/reference_d_spatial_acceptance.toml` declares for it. The output is the
coarsest pitch at which each observable is still inside its own declared bound
-- which is what `select_pitch` needs and has never had, because on real data
there is no truth to compare against.

NOTHING HERE TOUCHES THE TARGET. The objects are invented, so this calibrates
the METHOD and never the biofilm: it says what a pitch costs a measurement, not
what pitch the specimen needs.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import tomllib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "calibration"))

from biofilm_calibration.spatial.structure import summarise
from biofilm_calibration.spatial.synthetic import (NonTilingPitchError,
                                                   PhysicalSlab, PhysicalSpheres)

ACCEPTANCE = REPO / "config" / "reference_d_spatial_acceptance.toml"

# Which declared tolerance governs which observable. Named explicitly rather
# than matched by string similarity: `porosity` is 1 - biovolume, so an equal
# relative tolerance on it is ~20x tighter, and the config carries a separately
# derived number for exactly that reason.
TOLERANCE_OF = {
    "biovolume_fraction": "maximum_biovolume_fraction_error",
    "porosity": "maximum_porosity_error",
    "specific_interface_area_per_um": "maximum_interface_area_error",
    "correlation_length_um": "maximum_spatial_correlation_length_error",
    "thickness_mean_um": "maximum_thickness_error",
    "component_size_q50_um3": "maximum_component_size_error",
}


def load_tolerances() -> dict:
    with open(ACCEPTANCE, "rb") as fh:
        doc = tomllib.load(fh)
    out = {}
    for section in doc.values():
        if isinstance(section, dict):
            for key, value in section.items():
                if key.startswith(("maximum_", "minimum_")) and isinstance(
                        value, (int, float)):
                    out.setdefault(key, float(value))
    return out


def truth_for(spec) -> dict:
    """Closed-form values. Only what is genuinely exact appears here.

    `correlation_length_um` and the thickness SPREAD are deliberately absent:
    the first has no elementary closed form for a sphere lattice, and the second
    is zero for a slab only in the continuum, so quoting either as truth would
    manufacture a target the geometry does not actually have.
    """
    out = {"biovolume_fraction": spec.true_biovolume_fraction,
           "porosity": 1.0 - spec.true_biovolume_fraction,
           "specific_interface_area_per_um":
               spec.true_specific_interface_area_per_um}
    if isinstance(spec, PhysicalSlab):
        out["thickness_mean_um"] = spec.true_thickness_mean_um
    else:
        out["component_size_q50_um3"] = spec.true_component_size_um3
        out["n_components"] = float(spec.n_spheres)
    return out


def relative_error(got: float, want: float) -> float:
    if want == 0:
        return float("inf") if got else 0.0
    return abs(float(got) - float(want)) / abs(float(want))


def run_ladder(spec, pitches, tolerances) -> list[dict]:
    truth = truth_for(spec)
    rows = []
    for pitch in pitches:
        try:
            B = spec.rasterize(pitch)
        except NonTilingPitchError as err:
            # NAMED, not silently dropped. A pitch that cannot tile this
            # object is not a coarser view of it, and a ladder row that
            # quietly disappeared would look like a pitch nobody tried.
            rows.append({"pitch_um": float(pitch), "skipped": str(err)})
            continue
        # Threshold occupancy at tau = 0.5, the frozen declared mapping. The
        # fields here are already binary, so this is an identity for them and
        # is written out anyway: an occupancy rule that appears only when it
        # changes something is a rule nobody can see.
        mask = B >= 0.5
        got = summarise(B, mask, pitch).as_dict()
        row = {"pitch_um": float(pitch), "shape": list(B.shape),
               "voxels": int(B.size), **{f"got_{k}": v for k, v in got.items()}}
        for name, want in truth.items():
            err = relative_error(got[name], want)
            row[f"err_{name}"] = err
            key = TOLERANCE_OF.get(name)
            if key and key in tolerances:
                row[f"within_{name}"] = bool(err <= tolerances[key])
        # JSON has no NaN, and a report only some parsers can read is not
        # machine-readable — the same rule S0Verdict follows. A NaN here is
        # meaningful (correlation_length returns one when the field never
        # decorrelates, which a slab uniform in-plane never does), so it becomes
        # an explicit null rather than being dropped.
        rows.append({k: (None if isinstance(v, float) and not math.isfinite(v)
                         else v) for k, v in row.items()})
    return rows


def coarsest_passing(rows, observable) -> float | None:
    """The coarsest pitch at which an observable is inside its tolerance AND
    stays inside at every finer pitch.

    The monotonicity requirement is not pedantry: rasterisation error is not
    monotone in pitch — a sphere can be captured well at one pitch and badly at
    the next by where the centres happen to fall — so a single passing coarse
    pitch is luck, not headroom.

    A SKIPPED PITCH IS NOT A FAILED MEASUREMENT. A refused rung carries no
    `within_*` key, so reading it with a default of False made an explicit
    refusal indistinguishable from an observation that missed tolerance -- and
    it broke the tail beneath it. Inserting a non-tiling 1.4 between 1.6 and
    0.8 could therefore drag a passing component-size result from 1.6 down to
    0.8 with neither valid measurement having changed.

    The refusal says this pitch is not a view of the object at all. Absence of
    an observation is not an observation of failure, so the row stays in the
    report and leaves the convergence tail alone.
    """
    key = f"within_{observable}"
    # Walk coarse -> fine and keep the coarsest rung whose tail is unbroken.
    #
    # A previous version computed a running `passing` in a first loop and then
    # DISCARDED it, returning a value from this one -- so the filter I first
    # added to that loop changed nothing at all, and the test still failed with
    # the fix apparently in place. Dead code that looks like the answer is
    # worse than no code: it absorbs a correction and reports success.
    ordered = sorted((r for r in rows if "skipped" not in r),
                     key=lambda r: -r["pitch_um"])
    for i, row in enumerate(ordered):
        if all(r.get(key, False) for r in ordered[i:]):
            return row["pitch_um"]
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--outdir", type=Path, required=True)
    # 3.2 TILES THE SLAB AND NOT THE SPHERES, and the default keeps it for
    # exactly that reason. It does not divide the spheres' 24 um axis (7.5
    # voxels), so that row was measuring a 25.6 um box against a 24 um analytic
    # truth -- a 6.7% denominator shift reported as rasterisation error.
    #
    # The first correction DROPPED 3.2 from this list, which was wrong twice
    # over: the slab's 3.2 control is valid and was lost, and `run_ladder`'s
    # skip-row machinery only runs on a pitch it is actually given, so the
    # promise that the ladder "names the pitch it skipped" quietly became
    # unreachable. A pitch that is absent from the list is indistinguishable
    # from a pitch nobody tried. Keep it; let the spheres refuse it by name.
    ap.add_argument("--pitches", default="3.2,1.6,0.8,0.4,0.2")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    pitches = [float(p) for p in args.pitches.split(",")]
    tolerances = load_tolerances()

    specs = {"slab": PhysicalSlab(), "spheres": PhysicalSpheres()}
    doc = {"schema_version": 1, "tier": "S0", "target_calibration": False,
           "study": "morphology_rasterization_ladder",
           "acceptance_policy": str(ACCEPTANCE.relative_to(REPO)),
           "occupancy_mapping": "threshold_tau_0p5",
           "pitches_um": pitches, "systems": {}}

    for name, spec in specs.items():
        rows = run_ladder(spec, pitches, tolerances)
        truth = truth_for(spec)
        doc["systems"][name] = {
            "spec": {k: (list(v) if isinstance(v, tuple) else v)
                     for k, v in spec.__dict__.items()},
            "truth": truth,
            "coarsest_passing_um": {
                obs: coarsest_passing(rows, obs)
                for obs in truth
                if any(f"within_{obs}" in r for r in rows)},
            "rows": rows,
        }

        print(f"\n=== {name} ===")
        print("  truth: " + "  ".join(f"{k}={v:.5g}" for k, v in truth.items()))
        header = [o for o in truth if any(f"within_{o}" in r for r in rows)]
        print(f"  {'pitch':>7} {'shape':>14} " +
              " ".join(f"{o[:14]:>15}" for o in header))
        for row in rows:
            cells = []
            if row.get("skipped"):
                print(f"  {row['pitch_um']:7.2f} {'SKIPPED':>14}  "
                      f"does not tile this extent")
                continue
            for obs in header:
                mark = "ok " if row.get(f"within_{obs}") else "OVER"
                cells.append(f"{row[f'err_{obs}']:10.2e} {mark}")
            print(f"  {row['pitch_um']:7.2f} {str(tuple(row['shape'])):>14} " +
                  " ".join(cells))
        print("  coarsest pitch still inside tolerance (and at every finer):")
        for obs, pitch in doc["systems"][name]["coarsest_passing_um"].items():
            print(f"    {obs:34} {pitch if pitch else 'never in this range'}")

    (args.outdir / "rasterization_ladder.json").write_text(
        json.dumps(doc, indent=2, allow_nan=False, default=float) + "\n")
    print(f"\nwrote {args.outdir / 'rasterization_ladder.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
