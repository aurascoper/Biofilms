# Pre-registered tests: the radial ensemble at N=40, six parcels, four arms

**Recorded 2026-09-24, before the sweep was run.** Committed and pushed before
`radial_sweep.jl` executed over seeds 42 to 297, so that the thresholds cannot be moved
to fit what comes back. The same thresholds are constants in `radial_analyse.jl`, and
`test_statistics.jl` asserts that the two agree. If a primary test fails, the prediction
was wrong and the ledger row RAD-ENS-01 says so.

**What had been seen when this was written.** Seed 42 at this configuration (the
manuscript's Figure 1 and its §6 radial paragraph): *C. neoformans* mean r about 13.1,
*B. subtilis* about 9.9, a gap of about 0.15 R, opposite in sign to the direct term's
direction, about 1.4 standard errors from a centre-only null. The 256-seed melanin
ensemble (ledger rows SEED42-SCOPE-01 and MEL-N40-01), which reports the melanin field
and its acceptance bias, and no radial quantity. The term-magnitude audit: the melanin
term reaches 13.4% to 29.2% acceptance bias across seeds; the direct term's reach is
7.505e-2 in ΔH, a 1.5% bias. No radial quantity from any seed other than 42 had been
computed when this file was written.

## The question

Does either radiation-derived term of the CPM acceptance rule move a species' mean
radius? Four arms on the same 256 seeds, same placement per seed:

| arm | what is zeroed | dH_rad | dH_mel |
|---|---|---|---|
| `full` | nothing | on | on |
| `no_rad` | `beta_ion` for every species | 0 | on |
| `no_mel` | `alpha_M_species` for every species, so M stays 0 | on | 0 |
| `null` | both | 0 | 0 |

The observable is `mean_r_over_R`: the mean over a species' parcels of the centre-of-mass
radius, divided by R = N/2 = 20, read at MCS 100. Each comparison is paired within a seed.

## Primary tests, two of them

**P1. The melanin path moves the two species it acts on.** Under `full − no_mel`, for
*C. neoformans* (species 1) and *C. sphaerospermum* (species 3), each separately:
exact two-sided sign test over seeds at α = 0.05, AND |mean Δ r/R| ≥ 0.05. Both
conditions are required, because at n = 256 a sign test detects a consistent shift far
smaller than one the lattice can resolve. 0.05 r/R is one lattice unit at N = 40 and is
smaller than a seeded parcel's own radius (2 to 3 sites), so a mean displacement below it
is below the model's spatial grain.

Predicted: **both move inward** (mean Δ < 0: with the melanin term, the mean radius is
smaller). The melanin term is the only term through which radiation reaches the acceptance
rule at a magnitude above one part in a hundred, and its sign prefers melanin-rich sites,
which the drive field concentrates near the axis.

**P2. The direct term does not move any species by a resolvable amount.** Under
`full − no_rad`, for each of the seven species: two one-sided tests at α = 0.05, which is
the 90% confidence interval of the paired mean lying inside ±0.05 r/R. A sign test is
the wrong instrument here: a 1.5% acceptance bias applied consistently over 256 seeds is
detectable by sign and irrelevant by size, so the test asks about size.

Predicted: **equivalent within the bound, for every species.**

Either answer to either test is publishable. P1 failing means the term-magnitude audit
mis-predicted the dynamics and the melanin coupling does not stratify at this
configuration. P2 failing means the direct term does more than its ΔH reach suggests,
which is a finding about the acceptance rule.

> **Result, added after the run (2026-09-24; everything above is preserved as written).**
> P1 failed: under `full − no_mel` the paired mean shift is +0.0001 r/R for *C. neoformans*
> (129 of 256 seeds positive, p = 0.95) and −0.0004 r/R for *C. sphaerospermum* (125 of 256,
> p = 0.75). The prediction that both move inward was wrong. P2 held: under `full − no_rad`
> every species' 90% interval lies inside ±0.05 r/R, the widest being [−0.0029, −0.0006].
> Exploratory: `full − null` is within ±0.002 r/R for every species with within-seed
> correlation above 0.98, so the arms track one another and the inward drift from 0.53 to
> about 0.51 r/R is present with both terms off. Run on x86_64-conda-linux-gnu under Julia
> 1.12.7; the reference check matched 112 of 112 rows of the Mac-generated melanin sweep.
> Receipts: `radial_n40_p6_seeds42-297.csv`, `radial_analysis_n40_p6.txt`; ledger RAD-ENS-01.

## Exploratory, reported but not tested

Every other species × contrast cell (`full − null`, `no_rad − null`, the five species not
named in P1 under `full − no_mel`), and the *C. neoformans* − *B. subtilis* separation
within each arm and its `full − null` difference. These are reported with the same
paired statistics and no verdict. The seeding frame (MCS 0) is reported as the measured
mean r/R per species over seeds; it is what the manuscript's centre-only null
approximated.

## What a mismatch would and would not mean

The four arms change only which `CPMParams` vectors are zero. A difference between arms
is the model's response to that; it is evidence about the code, and about nothing else.
`alpha_M` and `beta_ion` are declared inputs. The `full` arm's melanin at MCS 100 is
checked against the committed melanin sweep for seeds 42 to 57 at six printed decimals,
which binds this producer to the stream the manuscript's other ensemble numbers came
from. If the check passes on one machine and fails on another, that is the machine, and
the CSV header names it.
