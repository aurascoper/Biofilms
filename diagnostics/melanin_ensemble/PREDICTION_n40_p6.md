# Pre-registered prediction: the melanin ensemble at N=40, six parcels

**Recorded 2026-09-13, before the sweep was run.** Committed and pushed before
`sweep.jl --n 40 --parcels 6` executed, so that the prediction cannot be quietly
adjusted to whatever came back. If the result disagrees, the prediction was wrong
and the ledger row says so.

## The question this answers

`post-biofilms-standalone.html` states it twice, and it is registered nowhere else
in this repository — no ledger row, no doc, no PR body:

> Sixteen seeds at the preprint's own configuration, forty voxels with six parcels.
> Larger parcels average over more sites and the variance should fall. **If eleven of
> sixteen becomes sixteen of sixteen there, then this is a small-lattice artifact and
> the published figure is fine as drawn.**

## What is predicted, and it is the opposite of that

The reasoning is one line: **averaging buys less than the gap loses.**

Occupied sites per species go 234 → 698, a factor of 3.0, so a standard deviation
falling as `1/sqrt(sites)` improves by about 1.7x. But the CN−AN gap at the published
configuration is 0.0347 (from `tests/fixtures/serial_seed42.csv`, seed 42) against
0.373 at N=20 — smaller by 10.7x. The CS−CN gap barely moves, 0.544 → 0.4638.

Per pair, using the PAIRED sd from `analysis.txt` rather than a pooled one, because the
three producers share one lattice per run and the within-run correlation is `-0.585`
for CN−AN against `+0.111` for CS−CN:

| pair | gap N=20 | gap N=40 | separation N=20 | separation N=40 | sign test N=20 | **predicted N=40** |
|---|---|---|---|---|---|---|
| CS − CN | 0.544 | 0.4638 | 1.43 | 2.11 | 15 of 16 | **~16 of 16, hardens** |
| CN − AN | 0.373 | 0.0347 | 0.91 | 0.15 | 12 of 16 | **~9 of 16, a coin flip** |
| full `alpha_M` ordering | | | | | 11 of 16 | **~9 of 16** |

**Predict the pairs, not the total.** The aggregate count moves 11 → about 9, which
reads as "no change", while underneath CS−CN hardens and CN−AN collapses to chance.
A single k-of-16 can stay flat while its composition inverts, and "eleven of sixteen"
is exactly the framing that invites reading the total.

## The free oracle, and what a mismatch would and would not mean

`sweep.jl` calls `run_simulation`; `tests/fixtures/serial_seed42.csv` was produced by
`validate_serial.jl` calling `run_simulation_coupled`, at this same configuration and
seed. So seed 42's three producers are already pinned and must come back
**CS 1.43716, CN 0.97333, AN 0.93862**.

A mismatch would otherwise be ambiguous between "coupling moves melanin" and "the two
entry points consume the RNG differently". **The second is ruled out in advance, by
reading rather than by assuming:** `biofilms_potts.jl:1051` and `:1455` both seed
`MersenneTwister(seed)`, both call `init_state(params; seed = seed)`, and both loops
call `mcs_step!(state, rng)` once per sweep at the same point. The coupled loop adds
`compute_radial_biomass`, `step_radiolysis!`, `radial_to_3d!` and
`update_nutrient_coupled!`, none of which takes `rng`. The streams are identical.

So if seed 42 disagrees, the coupling moves the melanin field, and that is a finding
about the model rather than about this sweep.

## What this cannot establish

A plot of simulation output is evidence about the code, not about organisms. The
observable is the volume-weighted mean over occupied sites, which `sweep.jl`'s own
header notes is a different quantity from the mean of per-parcel means. And `alpha_M`
is a declared input, so an inverted ordering does not contradict the parameter file —
that correction came from external review on PR #35 and stands.
