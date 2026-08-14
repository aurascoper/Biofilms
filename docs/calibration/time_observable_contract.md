# The dynamic-observable contract

**Date:** 2026-08-14 · **Status:** no observable selected ·
**Candidates:** `data/calibration/spatial/time_observable.csv` ·
**Code:** `calibration/biofilm_calibration/spatial/time_observable.py`

## The dependency was stated wrong

The integration branch recorded that seconds/MCS is downstream of the lattice pitch, which is
true and incomplete. The mapping is

```
Δt_MCS = a² · S_sim / S_exp
```

so a pitch is necessary. But `S_exp` has to be measured on something **the experiment and the
model both represent**, and under parcel semantics that is not automatic. The real chain is:

```
pitch selected  →  dynamic observable declared  →  seconds/MCS
```

Both `docs/online_coupling_design.md` (gate 5) and the ledger row for
`cpm.seconds_per_mcs_calibration` now say so.

## Why single-cell tracking is not automatically valid

The obvious `S_exp` is a mean-squared displacement from time-lapse tracking of individual cells.
Under the declared entity semantics that is a **category mismatch**: a `cell_id` is a
computational biomass parcel, not an observable bacterium, so a per-cell MSD and a parcel MSD are
different quantities. It remains usable only with a declared mapping from tracked cells to
parcels — which does not exist and would itself need justifying.

This is the same error the biomass-field harness exists to avoid on the spatial side, appearing
again on the temporal side.

## Three conditions, all required

An observable is selectable only if all three hold. The table keeps them as separate columns
because they fail for different reasons and the distinction matters.

**`represented_by_model`** — judged against the four terms actually in `compute_delta_H`
(adhesion, volume, radiation, melanin). The CPM has motility and shape rearrangement driven by
surface minimisation, and nothing else. It has **no nutrient-driven biomass creation**:
`state.nutrient` is written but never read by the acceptance test, and `divide_cell!` has no
trigger.

**`measurable`** — an experimental counterpart exists and can be acquired at the required time
resolution.

**`semantically_compatible_with_parcel`** — the model quantity and the experimental quantity
describe the same kind of object.

## The candidates

| Observable | Model | Measurable | Parcel-compatible | Verdict |
|---|---|---|---|---|
| biomass autocorrelation decay | ✓ | ✓ | ✓ | strongest — field-level on both sides, so it needs no cell-to-parcel mapping, and `structure.correlation_length` already computes the model side |
| aggregate centroid displacement | ✓ | ✓ | ✓ | aggregates are parcel-scale; needs aggregate identity to persist between frames |
| interface rearrangement rate | ✓ | ✓ | ✓ | directly driven by the adhesion and volume terms, which is what the CPM's dynamics *are*; sensitive to the occupancy threshold |
| radial redistribution | ✓ | ✓ | ✓ | **blocked**: needs a gradient the model responds to, and the nutrient field is never read by `compute_delta_H` |
| single-cell MSD | ✓ | ✓ | ✗ | the conventional choice, and the one the earlier plan assumed |
| biofilm thickening rate | ✗ | ✓ | ✓ | **unsupported**: growth-dominated, and the model cannot thicken |
| specific growth rate | ✗ | ✓ | ✗ | **unsupported**, for two independent reasons |

Four satisfy all three conditions. **None is selected**, because choosing among them is an
experimental-design decision, and because two of the four need something the acquisition contract
does not yet require: time-resolved stacks rather than a single stack.

## Growth-dominated observables are excluded on principle

Biofilm thickening and specific growth rate are marked `unsupported_by_current_model`, not
`blocked`. The distinction is the point: a blocked quantity awaits a measurement, an unsupported
one awaits a **model change**. No quantity of thickening data can calibrate a simulation that
cannot thicken.

The same status now applies to `response.growth_survival` in the transport ledger, for the same
reason.

## What this changes

- The spatial gate blocks on a missing observable as well as a missing pitch, so
  `READY_FOR_TIME_CALIBRATION` cannot be reached by acquiring morphology alone.
- `unsupported_by_current_model` joins the shared status vocabulary, and — deliberately — the set
  of statuses that assert a quantity carries no value, so it cannot fall through that check the
  way a new status silently would.
- The acquisition contract gains a consequence: if the selected observable is autocorrelation
  decay or interface rearrangement, the campaign needs **time-resolved** stacks, which is a
  materially harder acquisition than a single stack per condition.
