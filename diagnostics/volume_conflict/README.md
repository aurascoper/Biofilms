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

## A conflict-serializing scheduler, verified against the oracle

`scheduler.jl` layers each colour class's proposals by the conflict DAG, in reference order:

```
batch(i) = 1 + max{ batch(j) : j < i, ids(j) ∩ ids(i) ≠ ∅ }      (max over ∅ = 0)
```

`ids(p)` is **both** positive endpoints. A growth of parcel A and a shrinkage of parcel A
conflict exactly as two growths do, and the contention survey shows mixed-role pairs are
roughly a third of all conflicts, so a one-sided rule misses them.

**The enabling property** is that a colour class's proposal *set* is fully determined by the
pass-start lattice, independently of any acceptance in that pass. The RNG is counter-based on
`(seed, step, linear index)`; two distinct class members differ by ≥ 2 in some axis while a
site's source is within ±1 of itself, so one member's source can never *be* another member,
and the only lattice write a proposal makes is to its own target. So `(donor, recipient)` is
knowable up front and the schedule can be built before evaluating anything. Only `vols` is
contended.

### It reproduces the oracle exactly

50 MCS x 8 colours, N=20, 2 cells/species:

| | |
|---|---|
| proposals compared | **47442** |
| proposals differing | **0** |
| final lattice identical | **yes** |
| final volumes identical | **yes** |

Within a batch no two proposals share a parcel id, so evaluating them all against one common
`vols` is identical to evaluating them one at a time. That is the correctness argument, and
it is what a concurrent launch of a batch would actually do.

### Order preservation is load-bearing, and the control shows it

Reversing the within-class enumeration and re-layering produces batches that are **still
conflict-free** — and the trajectory does **not** match the oracle, in either the lattice or
the volumes. **Conflict-freedom alone is insufficient.** Being concurrency-safe is not the
same as being semantics-preserving, and an arbitrary safe partition can reorder conflicting
pairs and change the dynamics.

### The cost: parallelism collapses

The batch count is the critical path of the conflict DAG, and that is what bounds available
concurrency. Measured at the initial state:

| N | cells/species | parcels | proposals/pass | batches | parallelism |
|---|---|---|---|---|---|
| 20 | 2 | 14 | 92.4 | 20.9 | 4.43 |
| 20 | 6 | 42 | 160.9 | 36.6 | 4.39 |
| 40 | 2 | 14 | 117.8 | 14.6 | 8.05 |
| 40 | 6 | 42 | 284.0 | 22.1 | **12.84** |
| 60 | 6 | 42 | 305.8 | 17.5 | 17.47 |
| 40 | 12 | 84 | 527.0 | 26.0 | 20.27 |

It improves with parcel count, and it **degrades as the film densifies**: the same N=20/2
configuration averages **2.61** over 50 MCS against 4.43 at MCS 0, because sites per parcel
rise as the consortium grows.

**The current kernel launches every proposal in a class concurrently** — 92 to 527 of them.
Correct serialization admits 2.6 to 20. That is roughly a **10x to 45x** reduction in
available parallelism, on top of a backend already measured 2.3x slower than saturated CPU
threads on the two field kernels.

### What this settles

A correct concurrent CPM on this kernel is not a scheduling optimisation away. The conflict
structure — many proposals contending over few parcels, worsening as the biofilm grows — does
not admit the parallelism the current launch assumes. **Do not port**, and the reason is now
structural rather than a benchmark result.

The scheduler is not proposed as production code. It is the reference implementation any
candidate concurrent scheme must be checked against, and the answer it gives is that the
speedup such a scheme could earn is small.

### Files

| | |
|---|---|
| `scheduler.jl` | `propose` (geometry only), `schedule` (DAG layering), `scheduled_pass!` |
| `scheduler_test.jl` | exactness against the oracle, plus the reversed-order control |

## The half-window parity effect, attributed

`jacc_parity_tests.jl:209` asserts `s1.V < V_MAX && s2.V < V_MAX`. Julia prints no operand
values for a `&&` compound, so the three Metal failures named neither the offending half nor
its value. `parity_halves.jl` reports both, reusing the shipped `run_tables` / `pooled` /
`parity_stats` / `PERMS` so the statistics are the file's own.

| seed / ordering | rate | pooled V | maxdev | s1.V | s2.V | over |
|---|---|---|---|---|---|---|
| 42 / identity | 0.37039 | 0.01036 | 0.0269 | 0.01564 | 0.02101 | – |
| 42 / perm1 | 0.36775 | 0.02313 | 0.0637 | 0.01351 | **0.03228** | second |
| 42 / perm2 | 0.36665 | 0.01844 | 0.0317 | 0.01597 | 0.02458 | – |
| 43 / identity | 0.34360 | 0.01243 | 0.0423 | 0.01141 | 0.01915 | – |
| 43 / perm1 | 0.34617 | 0.01893 | 0.0405 | 0.01226 | **0.02883** | second |
| 43 / perm2 | 0.34094 | 0.02078 | 0.0551 | 0.01803 | **0.02517** | second |
| 44 / identity | 0.33830 | 0.01453 | 0.0266 | 0.01232 | 0.01921 | – |
| 44 / perm1 | 0.34302 | 0.01761 | 0.0448 | 0.01945 | 0.02300 | – |
| 44 / perm2 | 0.32617 | 0.01442 | 0.0367 | 0.01717 | 0.01954 | – |

**Rate outside band 9 of 9. Pooled V or maxdev over: 0 of 9. First half over: 0 of 9. Second
half over: 3 of 9.** That reconciles the log exactly: 9 failures at line 200 and 3 at line 209.

### Repeatability first -- and it fails

**Metal is not reproducible on identical input.** Same seed, same ordering, three replicates:

| replicate | rate | pooled V | s1.V | s2.V |
|---|---|---|---|---|
| 1 | 0.37143 | 0.00729 | 0.01125 | 0.00986 |
| 2 | 0.36655 | 0.01049 | 0.01127 | 0.01477 |
| 3 | 0.37040 | 0.01421 | 0.01805 | 0.01941 |

Accepted-count L1 difference between replicates 1 and 2: **6069**. Threads, the same probe:
all three replicates **byte-identical**, rate 0.19137 every time.

**This had to come first and did not.** Two claims made before it are now withdrawn:

1. *"All three half-window failures are the SECOND half."* That was one sample. A later
   `runtests.jl` run on the same code and seeds gave **1 failure at `s1.V` and 2 at `s2.V`**,
   and the three replicates above put **none** over threshold. Which half drifts is not a
   property of the configuration.
2. *"s2.V > s1.V in 9 of 9 runs, sign test p ≈ 0.002."* Those nine were nine different
   (seed, ordering) configurations observed **once each**, not nine replicates, and the
   p-value assumed a stability the backend does not have. Three replicates of one
   configuration give 2 of 3. The directional claim needs replicates per configuration before
   it can be made at all.

What **is** reproducible is the rate: 0.366-0.371 across replicates, always far outside the
(0.14, 0.23) band, 9 of 9 configurations in every run observed. The acceptance-rate failure is
a property of the backend. The half-window excursions are sporadic, and attributing them to a
half, a seed or an ordering requires replicates that have not been run.

**Pooling still hides something real.** Pooled V stayed under 0.025 in every run observed
while individual half-windows reached 0.03228. That much survives -- but it is a statement
about what the pooled statistic can see, not a measured time trend.

### A note on what the class index means here

`jacc_parity_tests.jl:26-35` records that `c` indexes both spatial class and sequence
position, and that `vols` accumulating across passes is a documented staleness -- which is why
the file uses random permutations rather than reversal. That effect is *across* passes and is
present on both backends. What concurrent execution adds is staleness *within* a pass, and it
is that addition which the second-half excursion tracks.

### The shipped assertion is now split

`@test s1.V < V_MAX && s2.V < V_MAX` is two assertions in `jacc_parity_tests.jl`. Identical
pass/fail, but a failure names the half and prints its value:

```
  Expression: s1.V < V_MAX
   Evaluated: 0.031991940528837394 < 0.025
```

Threads stays green at **74/74** (was 65 -- exactly +9, one per run). It is not a relaxation.

It earned its keep immediately: the first Metal run after the split reported a **first**-half
failure, which is what showed that the earlier second-half attribution was one sample.

`parity_halves.jl` asserts nothing by design -- it reports, so it can be run against a
failing backend without masking the shipped result.
