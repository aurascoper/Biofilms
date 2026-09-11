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

## The isolation test

Two arms, identical proposal streams, differing only in when `vols` is read. Both run the
**full 50-MCS trajectory** the parity tier runs.

| arm | proposed | accepted | rate |
|---|---|---|---|
| live `vols` (reference semantics) | 47902 | 4540 | **0.09478** |
| sweep-start snapshot | 50807 | 10583 | **0.20830** |
| ratio | | | **2.198x** |

**Measured on Metal vs threads at the same `n_mcs`: 0.34870 / 0.17933 = 1.944x.**

A sweep-start snapshot is the **maximum** staleness achievable within a colour pass -- every
read as stale as it can be. Live is 1.0x by construction. Concurrent GPU execution carries
**partial** staleness and lands at 1.944x, between the two. The volume conflict brackets the
observed effect.

### Why it compounds, and why a single pass does not show it

| MCS | live | snapshot | ratio |
|---|---|---|---|
| 10 | 0.22110 | 0.24874 | 1.125 |
| 20 | 0.15174 | 0.20604 | 1.358 |
| 30 | 0.12395 | 0.20355 | 1.642 |
| 40 | 0.10829 | 0.20324 | 1.877 |
| 50 | 0.09478 | 0.20830 | **2.198** |

The live arm's rate **decays** as the system equilibrates and volumes approach `V_target`,
so the penalty begins to bite. The snapshot arm stays roughly **flat**, because stale reads
never see the accumulated volume. The volume penalty's function is to accumulate; staleness
defeats accumulation, and the gap widens monotonically.

**An earlier revision of this file measured one pass, found ~5%, and concluded the volume
conflict could not explain a 94% gap. That was wrong, and the reason is instructive:
a single pass from a common state cannot exhibit a divergence whose entire character is
compounding.** Measuring an equilibration mechanism at t = 1 understates it by construction.

## What the first pass shows, mechanistically

One colour pass, identical initial state, both backends:

| | |
|---|---|
| proposal sets | **identical**, 99 = 99 sites -- the counter-based RNG and `nb26` agree |
| ΔH | **differs on 62 of 99 sites**, max 180, mean 70.3 |
| decisions differing | 6 |
| of those, with identical ΔH | **0** -- `exp()` and `u01()` are exonerated |
| of those, with differing ΔH | **6** -- the `vols` read is implicated |

The ΔH gaps quantise. Since the volume error is `2*lambda_V*(e_s - e_t)` with `lambda_V = 10`,
observed gaps of 120, 80 and 60 mean volume reads off by **6, 4 and 3 voxels**.

## Atomics are not dropping updates

`vols` is maintained incrementally by `JACC.@atomic` while `lat` carries the labels; the two
are redundant, so a positive-id histogram of `lat` must reproduce `vols` after any
synchronized pass. Run on device:

| backend | passes checked | mismatched |
|---|---|---|
| threads | 40 | **0** |
| **metal** | 40 | **0** |

So lost increments are **refuted**. Atomic addition is doing its job; the defect is that the
read-evaluate-decide-update sequence is not one transaction.

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
