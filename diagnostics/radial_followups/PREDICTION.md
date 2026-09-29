# Pre-registration: coefficient scales and a longer radial observation

Recorded 2026-09-29 before these production runs. Production requires this file,
the producer, analyser, tests, and ledger entry in a clean, named git revision.
Tiny implementation tests use different seeds, lattice size, and observation time;
their outputs are not evidence for these predictions.

The existing 256-seed experiment at N=40, six parcels per species, and 100 MCS
found small mean-radius contrasts at the shipped coefficients. Its full/null CSV
is `diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv`. The separate
synthetic direction test uses a different configuration and observable. It does
not determine a minimum effective coefficient here. No outcome from the new
coefficient levels or 400-MCS horizon has been inspected for this registration.

## Fixed design

- Seeds: every integer 42 through 297, paired within seed; 256 seeds per cell.
- Geometry: N=40; six initial parcels per species; R=N/2=20 lattice sites.
- Observable: equal-weight mean of parcel centre-of-mass radii for each species,
  divided by R. This is the existing `take_snapshot` observable, not the mean
  radius over occupied lattice sites.
- Random generator and stepping: the original `run_simulation` entry point with
  its `MersenneTwister(seed)` stream, with snapshots every 100 MCS.
- Factors: exactly 1, 10, 100, 1000, 10000, 100000, 1000000. Run every level;
  do not stop when a crossing first appears.
- Production writes one atomic CSV per seed and configuration outside git. Each
  file records the clean revision, hashes of model, producer, prediction and
  baseline, Julia version, machine, parameter vectors, and elapsed time.

P-beta multiplies the entire signed `beta_ion` vector by each factor, preserving
all signs and leaving melanin production at its default. P-alpha instead
multiplies the entire `alpha_M_species` vector, preserving its zeros, and leaves
`beta_ion` at its default. The hard-coded melanin coupling 0.5 is not edited.
Both sweeps end at 100 MCS and compare each scaled full arm with the existing
committed null (both vectors zero) at that seed and time. Thus a full-minus-null
contrast contains both pathways; it does not isolate the scaled pathway alone.
The factor-one runs in both sweeps are independent reproductions of the same
default full trajectory. All their reported rows must match the committed full
rows exactly at six printed decimals before a sweep can be analysed.
Before any scaled production run, reproduce that seed's default full and null
rows and require exact agreement with the committed baseline. A Julia-version
difference alone does not fail this check; a numerical difference does. Preserve
the runtime version in every shard. The null reproduction is a separate
`baseline/null` job; the beta factor-one job provides the full reproduction.

H-400 is a separate run of the default full and null arms to 400 MCS on the same
256 seeds. Both arms' 100-MCS prefixes must match their existing committed rows.
Report full-minus-null at 400 MCS for all species. Compare that paired contrast
with the existing 100-MCS contrast using the within-seed difference of contrasts.
Do not infer a trend at unobserved times from these two endpoints.

## Primary readouts and predictions

For P-beta and P-alpha separately, the primary readout is the smallest **tested**
factor whose C. neoformans (species 1) mean paired full-minus-null difference is
at most -0.05 in r/R: an inward shift of at least one lattice site. The direction
and inclusive boundary are fixed here. An equally large outward shift does not
satisfy this criterion. If no level qualifies, report "no tested factor met the
inward criterion". Do not call it a continuous threshold or a calibrated value.

Prediction P-beta: an inward mean shift of at least 0.05 r/R will occur at one or
more of the seven tested factors. Prediction P-alpha is deliberately weaker:
the sign of its radial response is undetermined, because attraction to a field
made near existing parcels can retain parcels rather than move them inward.
Its primary inward readout is still fixed above, including a no-crossing result.
Prediction H-400: the C. neoformans full-minus-null mean will remain between
-0.05 and +0.05 r/R at 400 MCS. These are model predictions, not organism claims.

The crossing is defined by the paired mean alone, not by a p-value or confidence
interval. Report every mean, paired SD, standard error, normal 95% interval and
seed count for every species and level. Intervals describe uncertainty over
initial seeds; they are not simultaneous intervals over the sweep, and they do
not correct selection of a first crossing. Flag when an interval crosses the
-0.05 boundary. Report nonmonotone responses and later loss of a crossing.

Do not drop seeds whose species disappears. The legacy snapshot writes radius
zero when no parcel remains; that is an absent observable, not a measurement of
position on the axis. Report disappearance counts. A species-level contrast is
not evaluable if either member of any seed pair lacks that species. A primary
first crossing is not evaluable if any tested CN level is incomplete or has an
absent observable. Missing, duplicate, mixed-provenance, nonfinite or mismatched
rows are errors, not evidence for a null result. A failed or interrupted shard
must be rerun; incomplete sweeps have no primary verdict.

## Separation from the manuscript build

These runs do not gate printing the corrected manuscript. No new scientific
result belongs in it before the registered runs finish and their checks pass.
Do not change thresholds, grid, exclusions or hypotheses after inspecting these
outputs. Record any later experiment as a new registration.

Run implementation tests with `julia --project=. diagnostics/radial_followups/test_sweep.jl`
and `python3 -m pytest diagnostics/radial_followups/tests/test_workflow.py`.
The production command and resume scheduler are documented beside the producer.
