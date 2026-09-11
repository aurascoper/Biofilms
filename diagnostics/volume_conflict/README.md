# The shared-volume conflict: real, and not the explanation

Deterministic host diagnostics for the acceptance-rate gap measured on Metal
(0.34870 vs 0.17933 on threads, **1.944x**). No GPU is involved here; everything below runs
single-threaded on the host so that read timing is the only variable.

## What was asked and what came back

The conflict is real. `delta_H` reads `vols[sigma_s]` and `vols[sigma_t]`; `cpm_color!`
atomically mutates the same entries on acceptance. The eight-colour decomposition separates
target sites under the Moore-26 stencil, which makes the `lat` reads race-free -- but it does
**not** make donor and recipient parcel ids disjoint.

Measured on a 20^3 lattice, 2 cells per species, one colour pass each:

| colour | proposals | conflicting pairs | donor / recipient / mixed |
|---|---|---|---|
| 3 | 89 | 490 | 142 / 164 / 184 |
| 4 | 90 | 457 | 107 / 183 / 167 |
| 5 | 89 | 450 | 118 / 137 / 195 |
| 6 | 87 | 430 | 89 / 157 / 184 |
| 7 | 99 | 501 | 238 / 76 / 187 |

**Mixed conflicts** -- one proposal's donor is another's recipient -- are consistently among
the largest categories. Any conflict detection that looks at one side only misses roughly a
third of them.

## The isolation test, and what it refutes

Two arms, identical proposal streams, differing only in when `vols` is read:

| | accepted over 8 colours |
|---|---|
| live `vols` (reference semantics) | **228** |
| sweep-start snapshot | **239** |
| ratio | **1.048** |

414 of 715 proposals get a different ΔH; **17 decisions flip**.

**A sweep-start snapshot is the maximum staleness achievable within a colour pass** -- every
read is as stale as it can be. Concurrent GPU execution is bounded by that. So 4.8% is an
**upper bound** on the within-pass volume-read contribution, against an observed gap of 94.4%.

**That is 5.1% of the effect. The volume conflict does not explain the Metal result.**

The direction is right -- staler reads accept more, matching Metal's sign -- which is exactly
how a wrong mechanism looks if you stop at the sign. An earlier revision of
`docs/metal_feasibility_macos_arm64.md` attributed the gap to this race. **Withdrawn.**

## What this does not show

It does not identify what does explain the gap. The leading untested candidate is the
`JACC.@atomic` updates themselves: if increments are lost on device, volumes stay near
`V_target`, the penalty stays small and acceptance stays high, which would be large rather
than marginal. `id_histogram` exists for exactly that test -- recompute a positive-id
histogram from the lattice after each synchronized colour pass and compare every entry with
`vols`. On this host oracle it matches after all eight colours. **It has not been run on
Metal, and that is the decisive next measurement.**

Other untested contributors: the device `exp` and `u01` implementations, and whether the
counter-based RNG produces identical words on both backends. Compare generated words and
draws before attributing anything to dynamics.

## Snapshot freezing is not a candidate fix

For `H(V) = lambda_V (V - V_target)^2` the joint change after two same-parcel growths exceeds
the sum of two changes evaluated at the starting volume by exactly `2*lambda_V`, with the
opposite sign for a growth paired with a shrinkage. Worked case at `lambda_V = 10`,
`V = V_target = 120`, `T = 5`, uniform draw `0.01`:

| evaluation rule | first | second | accepted |
|---|---|---|---|
| sequential, current volumes | ΔH 10, accept | ΔH 30, reject | 1 |
| sweep-start snapshot | ΔH 10, accept | ΔH 10, accept | 2 |

Joint volume-energy change 40 against a sum of evaluated changes of 20. **Freezing removes the
concurrent read by changing which moves are accepted**, so it is a different algorithm, not a
repair. It is included in `oracle.jl` as a diagnostic arm and is labelled as one.

Also withdrawn: "under-counted volume lowers ΔH". The volume contribution is
`lambda_V[2(V_s - V_target) + 1] + lambda_V[-2(V_t - V_target) + 1]`, so with read errors
`e_s`, `e_t` the energy error is `2*lambda_V*(e_s - e_t)`. Undercounting the **donor** lowers
ΔH; undercounting the **recipient** raises it. The net sign is a measurement.

## Files

| | |
|---|---|
| `oracle.jl` | deterministic replay of one colour pass; `snapshot=` selects the arm. Not the random-sequential serial CPM -- that is a different schedule and comparing against it would conflate two differences. |
| `contention_test.jl` | the three measurements above |

Run: `julia --project=. diagnostics/volume_conflict/contention_test.jl`
