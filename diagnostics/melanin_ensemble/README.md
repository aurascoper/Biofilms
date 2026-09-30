# Melanin seed ensemble

A committed seed sweep for the melanin observable, and a paired analysis of it.

This exists because the result it produces was previously run out of band. A write-up
reporting "seeds 42 through 57" cited a commit whose tree contains no seed loop and no
`--seed` flag: `biofilms_potts.jl` has exactly one `ARGS` reference and it is
`--no-radiolysis`. Seed enters only as a keyword defaulting to 42. A number nobody else can
regenerate is a claim rather than a measurement, which is this repository's own standard
applied to its own output.

## Running it

```sh
julia diagnostics/melanin_ensemble/sweep.jl /tmp/sweep.csv      # 16 seeds, about 11 s
julia diagnostics/melanin_ensemble/analyse.jl /tmp/sweep.csv
julia diagnostics/melanin_ensemble/test_statistics.jl           # 150 assertions, no data file
julia --project=diagnostics/melanin_ensemble diagnostics/melanin_ensemble/figure.jl   # the N=40 figure

julia diagnostics/melanin_ensemble/parcel_sweep.jl /tmp/parcel.csv --seeds 42:297   # 256 seeds, about 2 min
julia diagnostics/melanin_ensemble/parcel_sweep.jl --compare /tmp/parcel.csv theirs.csv
```

`figure.jl` needs CairoMakie, which the root `Project.toml` does not declare; this
directory's `Project.toml` and `Manifest.toml` pin the version that rendered the committed
PNG (`julia --project=diagnostics/melanin_ensemble -e 'import Pkg; Pkg.instantiate()'` once).

`sweep.jl` refuses a destination that already exists, an empty seed list, and an `--at`
outside the run, all before the model loads. `test_statistics.jl` takes no argument, no
path and no environment variable, so its assertions cannot be skipped by a missing sweep;
it runs `sweep.jl` itself for four MCS to bind the writer's header to the reader.

No project environment is needed. The model loads through the same `#  13. Figure export`
split-marker sandbox `validate_serial.jl` uses, whose imports are four stdlibs, so this
sidesteps the root `Project.toml`'s AMDGPU entry that cannot instantiate on darwin. It also
means no CairoMakie: this writes a CSV and leaves plotting to whoever wants it.

## What it reproduces

Configuration: `N = 20`, 2 parcels per species, 400 MCS, seeds 42 to 57, observable read at
MCS 100. Committed output is `sweep_seeds42-57.csv` and `analysis.txt`.

| | alpha_M | mean | sd | range |
|---|---|---|---|---|
| *C. sphaerospermum* | 0.140 | 2.192 | 0.274 | 1.603 – 2.587 |
| *C. neoformans* | 0.100 | 1.647 | 0.296 | 1.098 – 2.139 |
| *A. niger* | 0.065 | 1.274 | 0.162 | 0.943 – 1.540 |

11 of 16 seeds display the full `alpha_M` ordering. 4 put *A. niger* above *C. neoformans*;
1 puts *C. neoformans* above *C. sphaerospermum*; the two never coincide, so the counts sum.

## The comparison is paired, and that changes the answer

The three producers share one lattice in every run. A pooled standard deviation across
species therefore treats 16 runs as 32 independent draws. `analyse.jl` reports the paired
statistic, the independence-assuming one, and the assumption-free count side by side:

| pair | mean | paired sd | if independent | within-run r | seeds in the alpha_M direction | sign test |
|---|---|---|---|---|---|---|
| CS − CN | +0.544 | 0.380 | 0.403 | +0.111 | 15 of 16 | p = 0.00052 |
| CN − AN | +0.373 | 0.412 | 0.337 | **−0.585** | 12 of 16 | p = 0.077 |

**The two pairs are not alike, and a single "11 of 16" hides that.** *C. sphaerospermum*
above *C. neoformans* is reliable. *C. neoformans* above *A. niger* is not: at 12 of 16 the
exact two-sided sign test cannot reject a coin flip at the conventional level.

For `CN − AN` the paired spread (0.412) **exceeds** what independence predicts (0.337),
which means the two are negatively correlated within a run and a pooled statistic understates
the uncertainty. For `CS − CN` it falls below, and the pooled statistic overstates. Reporting
one pooled number for both would be wrong in opposite directions.

## What the observable is

Mean melanin over **occupied sites**, read straight from `take_snapshot`
(biofilms_potts.jl@`function take_snapshot(state::CPMState, mcs::Int)`). Every occupied voxel contributes one term, so this is already
volume-weighted: a parcel of 200 sites counts 200 times a parcel of 1.

The mean of **per-parcel means**, weighting each parcel equally, is a different quantity.
`parcel_melanin.jl` computes it, under the name `mel_parcel`. Worth knowing that `mean_r`
in the *same function* uses that other convention, so one snapshot reports two observables
weighted two ways.

*Corrected 2026-09-17:* this section stated that the per-parcel mean "is computed nowhere
in this repository". That was true when written and is false now. The line reference moved
as well: `take_snapshot` is at line 1049, and this file said 857.

`alpha_M` is a declared input. Agreement with it is a statement about how reliably one
trajectory displays an input, and says nothing about whether the input is right.

## Two implementations agree on both estimators, at two time points

F. Fink's Odin port computes `mel_parcel` and supplied 35,840 rows over seeds 42 to 297.
`parcel_sweep.jl` was run against that file twice, reading at MCS 100 and at MCS 400.

| read at | rows compared | volume | ncells | mel_site | mel_parcel |
|---|---|---|---|---|---|
| MCS 100 | 1,792 | 1792 of 1792 | 1792 of 1792 | 1792 of 1792 | 1792 of 1792 |
| MCS 400 | 1,792 | 1792 of 1792 | 1792 of 1792 | 1792 of 1792 | 1792 of 1792 |

Agreement is exact six-decimal string equality. The worst absolute difference on either
melanin column is 0.000e+00.

**What was not compared.** Their file covers 20 MCS values and these runs read two of them.
So 3,584 of 35,840 rows were checked and 32,256 were not. This is a result about MCS 100
and MCS 400, not about the whole file.

**The comparison refuses a one-digit difference.** Changing the last decimal place of a
single `mel_parcel` value out of 1,792 makes `--compare` exit non-zero and name that row.
A comparison that cannot do that reads the same as one that always agrees.

**Every seed checks itself before any of this.** `parcel_sweep.jl` recomputes `mel_site`
from the lattice and requires exact equality with the model's own `mean_melanin`, along
with the site and parcel counts. A misread lattice fails there, not in the comparison.

## A claim this sweep contradicts

A draft write-up states that every species holds between 232 and 234 sites for the whole run.
Measured across the sweep at MCS 100, 200, 300 and 400 (the MCS 0 rows are the seeding
frame, 2 to 246 sites, and are excluded) the range is **231 to 237**, with 152 of 448
samples outside. At
that draft's own seed and frame cadence (42, every 4 MCS) it is **226 to 237**, with 150 of
700 outside, roughly three times the stated width.

The occupancy totals in that draft are exactly right: 1,612 sites at MCS 4 and 1,633 at MCS
400, both reproduced here to the site. The conclusion it supports also survives, since 226 to
237 is still fixed biomass to about ±2.4%. The stated range is simply narrower than what the
model does.

## The radial ensemble: does either radiation term move a parcel?

The manuscript's only radiotropism result is seed 42 alone, and it says of itself that the
exact null must come from replicate simulations. `radial_sweep.jl` runs those replicates
under four arms of the acceptance rule on the same seeds, with the same placement per seed:

| arm | zeroed | dH_rad | dH_mel |
|---|---|---|---|
| `full` | nothing | on | on |
| `no_rad` | `beta_ion`, every species | 0 | on |
| `no_mel` | `alpha_M_species`, every species, so M stays 0 | on | 0 |
| `null` | both radiation-derived terms | 0 | 0 |

```sh
julia diagnostics/melanin_ensemble/radial_sweep.jl diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv \
      --seeds 42:297 --reference diagnostics/melanin_ensemble/sweep_n40_p6_seeds42-57.csv
julia diagnostics/melanin_ensemble/radial_analyse.jl diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv \
      > diagnostics/melanin_ensemble/radial_analysis_n40_p6.txt
```

Four arms of 256 seeds is about four times one `parcel_sweep.jl` pass. The observable is
`take_snapshot`'s `mean_r`, the mean over a species' parcels of the centre-of-mass radius,
each parcel weighted equally, and `mean_r_over_R` divides it by R = N/2. That is a
per-parcel mean; `mean_melanin` in the same row is a per-site mean, as above.

**The tests were registered before the run**, in `PREDICTION_radial_n40_p6.md`, with a
minimum effect of 0.05 r/R (one lattice unit at N = 40) for "moves" and the same value as
the equivalence bound for "does not move". `radial_analyse.jl` carries them as constants
and `test_statistics.jl` asserts the prose and the constants agree. Ledger row RAD-ENS-01
holds the question while it is open.

**The answer (2026-09-24, `radial_analysis_n40_p6.txt`).** Neither term moves a parcel. P1
failed: `full − no_mel` shifts *C. neoformans* by +0.0001 r/R (p = 0.95) and *C. sphaerospermum*
by −0.0004 (p = 0.75). P2 held: `full − no_rad` is inside ±0.05 r/R for all seven species, the
widest interval being [−0.0029, −0.0006]. `full − null` is within ±0.002 r/R everywhere with
within-seed correlation above 0.98, and the inward drift from the seeding frame (0.53 r/R) to about
0.51 r/R appears with both terms off. The *C. neoformans* − *B. subtilis* separation is −0.006
r/R with sd 0.11 over seeds; seed 42's +0.16 was one draw. Run on x86_64-conda-linux-gnu under
Julia 1.12.7, the full arm matching the Mac-generated melanin sweep on all 112 shared rows.
`radial_figure.py` renders the manuscript's Figure 1 from the CSV with matplotlib and writes the
PDF, `pdftotext -layout` sidecar, sha256 and PNG.

**Controls the producer runs on every seed and arm:** the mean radius recomputed from the
cell registry equals the snapshot's exactly; an arm with production zeroed prints a
melanin field of exactly 0; the MCS 0 seeding frame is identical across arms; and, with
`--reference`, the `full` arm's melanin at the read MCS prints identically to the
committed sweep's (six decimals) for every shared seed, and a reference sharing no seed
is refused as a check that never ran. **The header records `julia_version` and
`machine`**: Julia's streams match across platforms, floating-point results on Apple ARM
and x86 need not, and the manuscript will quote this file.

## Files

| | |
|---|---|
| `sweep.jl` | runs the seeds, writes one row per (seed, mcs, species) |
| `analyse.jl` | paired statistics, ordering counts, exact binomial sign test (ties leave the test) |
| `parcel_melanin.jl` | both melanin estimators from one lattice visit; no model, no data |
| `parcel_sweep.jl` | runs the seeds for `mel_parcel`, and compares two files by column |
| `stats.jl` | the sign test and the paired difference, shared by both analyses |
| `radial_sweep.jl` | the four-arm radial ensemble, with its per-seed controls |
| `radial_analyse.jl` | the pre-registered tests and the exploratory cells |
| `PREDICTION_radial_n40_p6.md` | the tests and thresholds, recorded before the run, and the result after it |
| `radial_n40_p6_seeds42-297.csv`, `radial_analysis_n40_p6.txt` | committed receipts of the run |
| `radial_figure.py` | the manuscript's Figure 1, rendered from the CSV |
| `test_statistics.jl` | assertions needing no data file; two testsets run the producers for four MCS |
| `sweep_seeds42-57.csv`, `analysis.txt` | committed receipts of the run above |

The assertion count rose from 114 to 150. Sixteen of the new ones cover the estimator, and
twenty cover the producer and the comparison. The `test_statistics.jl` row read 90 before
2026-09-17. That figure was already 24 short of what the file ran, so the rise is not 60
assertions of new work.
