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
julia diagnostics/melanin_ensemble/test_statistics.jl           # 47 assertions, no data
```

`sweep.jl` refuses a destination that already exists. `test_statistics.jl` takes no
argument, no path and no environment variable, so its assertions cannot be skipped by a
missing sweep.

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

| pair | mean | paired sd | if independent | seeds in the alpha_M direction | sign test |
|---|---|---|---|---|---|
| CS − CN | +0.544 | 0.380 | 0.403 | 15 of 16 | p = 0.00052 |
| CN − AN | +0.373 | 0.412 | 0.337 | 12 of 16 | p = 0.077 |

**The two pairs are not alike, and a single "11 of 16" hides that.** *C. sphaerospermum*
above *C. neoformans* is reliable. *C. neoformans* above *A. niger* is not: at 12 of 16 the
exact two-sided sign test cannot reject a coin flip at the conventional level.

For `CN − AN` the paired spread (0.412) **exceeds** what independence predicts (0.337),
which means the two are negatively correlated within a run and a pooled statistic understates
the uncertainty. For `CS − CN` it falls below, and the pooled statistic overstates. Reporting
one pooled number for both would be wrong in opposite directions.

## What the observable is

Mean melanin over **occupied sites**, read straight from `take_snapshot`
(`biofilms_potts.jl:857`). Every occupied voxel contributes one term, so this is already
volume-weighted: a parcel of 200 sites carries 200 times the weight of a parcel of 1.

The mean of **per-parcel means**, weighting each parcel equally, is a different quantity and
is computed nowhere in this repository. Worth knowing that `mean_r` in the *same function*
uses that other convention, so one snapshot carries two observables weighted two ways.

`alpha_M` is a declared input. Agreement with it is a statement about how reliably one
trajectory displays an input, and says nothing about whether the input is right.

## A claim this sweep contradicts

A draft write-up states that every species holds between 232 and 234 sites for the whole run.
Measured across the sweep the range is **231 to 237**, with 152 of 448 samples outside. At
that draft's own seed and frame cadence (42, every 4 MCS) it is **226 to 237**, with 150 of
700 outside, roughly three times the stated width.

The occupancy totals in that draft are exactly right: 1,612 sites at MCS 4 and 1,633 at MCS
400, both reproduced here to the site. The conclusion it supports also survives, since 226 to
237 is still fixed biomass to about ±2.4%. The stated range is simply narrower than what the
model does.

## Files

| | |
|---|---|
| `sweep.jl` | runs the seeds, writes one row per (seed, mcs, species) |
| `analyse.jl` | paired statistics, ordering counts, exact binomial sign test |
| `test_statistics.jl` | 47 data-free assertions |
| `sweep_seeds42-57.csv`, `analysis.txt` | committed receipts of the run above |
