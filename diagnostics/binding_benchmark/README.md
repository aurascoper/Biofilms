# Closed diffusion-and-binding benchmark on frozen MCS-0 geometry

Workstream B Stage 1. This is a numerical benchmark and nothing else. It integrates a
dissolved scalar and a bound scalar on the labels and mask of one frozen snapshot, and it
asks whether the implemented scheme conserves what it should, loses only what its ledger
books, and converges at the order it claims. It asserts no chemistry, names no molecule,
no species and no isotope, and it does not establish that the generating CPM ever had a
binding loop of this shape.

The thorium/oxalate model is not here and is not coded anywhere in this branch. It is
deferred to a separate specification packet, which owes species and pool definitions,
balanced stoichiometry, a resolution of `S` as exposed surface area with a flux-to-volume
conversion and a finite source inventory, ligand consumption, units, initial and boundary
conditions, rate laws, and a corrected source audit, before any line of it is written.

## Running it

```sh
P=~/Developer/biofilms-paraview-mac-2026-09-11/evidence/biofilm-signal-evidence/original_evidence
J="julia --project=diagnostics/binding_benchmark"

$J diagnostics/binding_benchmark/test_numerics.jl              # 81 assertions, no data
$J diagnostics/binding_benchmark/run.jl       $P /tmp/bench    # benchmark_receipt.json
$J diagnostics/binding_benchmark/controls.jl  $P /tmp/ctl      # control_verification.json
$J diagnostics/binding_benchmark/mutation_controls.jl /tmp/mut # mutation_verification.json
```

`test_numerics.jl` and `mutation_controls.jl` take no data path, no argument and no
environment variable, so their assertions cannot be skipped by a missing bundle. That is
deliberate: a suite whose assertions are reachable only when a data directory happens to
be mounted reports "passed" on a machine where it never ran, and this repository has
already shipped one file with exactly that defect.

`run.jl` and `controls.jl` refuse a destination that already exists, and check it twice:
once before the run so a typo costs nothing, once at write time because a long run can
finish into a directory that appeared while it was running. The directory is created only
at the second check, so a run that fails leaves nothing behind to block the retry.

This diagnostic carries its own `Project.toml`. The repository root project lists AMDGPU,
which cannot instantiate on darwin, so a benchmark sharing it could not be run on the
machine it was written on. The root `Project.toml` is untouched by this branch.

## The scheme

```
B(x)  = O(x) · [B0 + ΔB · h(A(x))],   h(A) = Aⁿ / (Kⁿ + Aⁿ)
r     = k_on · c · (B − b) − k_off · b
∂c/∂t = D_c ∇²c − r − λ c + q_ext
∂b/∂t = r − λ b
```

`c` is dissolved, `b` is bound, `B` is a finite volumetric capacity carried by occupied
sites, and `O` is occupancy. Closed means `q_ext = 0` and no flux through the boundary, so
the only thing that leaves is decay. Naming follows the specification deliberately, to
avoid two collisions that were live in an earlier draft: `λ` is the loss rate applied to
both pools and `B` is never a loss rate, and `q_ext` is the external source and never the
binding relaxation rate `q_bind = k_on·c + k_off + λ`.

Units are the diagnostic's own. Time is in `dtu`, a declared diagnostic time unit; lengths
are in lattice sites; `c` and `b` are in a declared arbitrary concentration unit. Nothing
here converts to seconds or metres, because that conversion is blocked on D-PITCH and
D-TIMESERIES. In particular the days used by the decay-reference illustration elsewhere in
this repository are that diagnostic's clock, not this one's, and importing them would
attach a physical meaning that nothing establishes.

There is exactly one adjacency rule in this diagnostic: 6-connected faces. The diffusion
stencil, the no-flux boundary and the capacity modulator all use it, so there is no second
rule to confuse with the first.

## The frozen input

| | |
|---|---|
| parent run | `manuscript42_seed42_8cbb8ee` |
| snapshot | `snapshots/snap_mcs000000.h5`, MCS 0 |
| chain of custody | one pin: `run_manifest.json` → `f00e470a…` |
| grid | 40 × 40 × 40, spacing 1 lattice site |
| interior sites | 50,200 |
| occupied sites | 3,038, all inside the interior (checked, not assumed) |

Only one hash has to be trusted. The manifest is hashed against the pin, the snapshot is
hashed against the entry the manifest itself carries for it, and the snapshot's own `mcs`,
`run_id`, `label_state_hash` and `mask_sha256` are then checked against the configuration.
Carrying a `git_sha` would say which code produced the bytes; it would not say that these
are the bytes.

Nothing in the parent bundle is written to. The 101-frame authoritative trajectory is not
regenerated, and no figure, protected CSV or manuscript artifact is touched.

`A` is read off the frozen labels, not off any field carried in the snapshot: it is the
fraction of a site's six face neighbours that are occupied. On this frame it takes six of
its seven possible values, with 193 sites at 1/6, 17 at 2/6, 545 at 3/6, 458 at 4/6, 930
at 5/6 and 895 at 6/6, so capacity spans 1.2 to 2.6 and integrates to 7017.9486. A
constant modulator is refused rather than run, because with a constant `A` the
matched-capacity control compares a field against itself and cannot fail. That guard is
not hypothetical: it fired on the first fixture written for it, a solid 2×2×2 block, in
which every site has exactly three occupied neighbours.

The snapshot's `melanin` field is identically zero at MCS 0, and `radiation_cpm` is the
only varying continuous field it carries. Neither is used. Deriving `A` from the labels
keeps the benchmark's one geometric input geometric.

## Declared constants

Every value below is chosen to exercise the scheme. None of it is measured, fitted, or
read off a source, and the receipt says so.

| | | |
|---|---|---|
| `D_c` | 0.1 | lattice² / dtu |
| `λ` | 0.01 | 1 / dtu, both pools |
| `k_on` | 0.05 | 1 / (concentration · dtu) |
| `k_off` | 0.02 | 1 / dtu |
| `B0`, `ΔB` | 1.0, 2.0 | concentration |
| `K`, `n` | 0.5, 2.0 | units of `A`, dimensionless |
| `q_ext` | 0.0 | closed |
| `c(0)`, `b(0)` | 1.0, 0.0 | uniform on the interior |
| `Δt`, `T` | 0.5, 50.0 | dtu, 100 steps |

## The timestep bound is derived here, not imported

The restriction is taken from this stencil and this right-hand side:

```
Δt · max( 2 D Σⱼ hⱼ⁻² + k_on·max(B − b) + λ ,  k_on·max(c) + k_off + λ ) ≤ 1
```

The first term is the magnitude of the six-point Laplacian's diagonal at a site with all
six faces interior. In 3-D with uniform `h` that is `h²/(6D)`, not the one-dimensional
`h²/(2D)`. The difference is not cosmetic: at `h = D = 1` and `Δt = 0.25`, which the 1-D
bound admits, a unit centre with six zero neighbours lands at −0.5. The suite asserts
that value and asserts that `step!` refuses the step that produces it.

Reaction and decay are not split off, so they enter the same restriction, and each of the
three contributing terms is pinned by its own assertion rather than by one combined
number. The rate is re-evaluated inside every step, because `max(c)` and `max(B − b)` both
move and a run can leave the admissible region after entering it. On the shipped
configuration the initial rate is 0.74, the bound is 1.3514, the declared step is 0.5, and
the margin is 0.37 at every recorded time (both maxima start at their largest and decrease
from there).

`Δt ≤ 1/λ` is a positivity condition and not an accuracy condition. At `Δt = 1/λ` the
scheme returns exactly zero against a continuum value of `1/e`, and the suite asserts
both halves of that.

## What the shipped run measured

| | |
|---|---|
| dissolved | 50200.000000 → 27816.012646 |
| bound | 0.000000 → 2593.663266 |
| decayed | 19790.324088 (18720.268988 dissolved, 1070.055100 bound) |
| bound fraction, of inventory | 0.085291 |
| bound fraction, of capacity | 0.369576 |
| closure residual | −2.91e−11, which is −5.8e−16 relative |
| transport residual | −1.49e−12 |
| min dissolved | 0.199895 |

Both bound-fraction denominators are reported, separately, throughout. They differ by a
factor of 4.33 on this run, so a single number called "the bound fraction" would be
wrong under one reading of it whichever one was meant.

Binding is internal and both pools lose at the same rate, so the total inventory obeys
`dI/dt = −λI` exactly. The final total, 30409.6759, is `50200 · (1 − λΔt)¹⁰⁰` to the
digits printed, which is an independent check on the whole integration that costs nothing.

## Controls

Runnable, not narrated: `controls.jl` writes `control_verification.json` and exits
non-zero if any verdict is unexpected. Two rows are expected to be red.

| | verdict | measured |
|---|---|---|
| K1 matched total capacity | MEASURED | integrals agree to 2.6e−15 at every comparison time; bound inventory differs by up to 49.19, which is 1.89% |
| K2 `ΔB = 0` reduction | PASS | local and uniform runs bit-identical, not merely close |
| K3 analytical binding limit | PASS | observed order 0.990, 0.995 against `b_eq + (b(0) − b_eq)e^{−q_bind t}` |
| K4 pure diffusion | PASS | inventory conserved to 1.3e−14; discrete maximum principle held; spread 1 → 0.4031 toward the interior mean 0.060518 |
| K5 pure decay | PASS | deviation from the exact discrete sequence 1.1e−16; from the continuum 7.60e−4 |
| K6 positivity, `0 ≤ b ≤ B` | PASS | min dissolved 0.199895, min bound 0, max overshoot 0 |
| K7 conservative release | FIRES | capacity cut to 1/4 at step 50, 639.5826 released, closure residual 1.0e−15 relative |
| K8 the same transfer, omitted | FIRES | 639.5826 dropped, closure residual −639.5826 (−1.27%), residual plus dropped −6e−11 |
| K9 timestep refinement | PASS | observed order 1.0013 on the final bound inventory |
| K10 both denominators | MEASURED | 0.085291 of inventory, 0.369576 of capacity, ratio 4.333 |
| K11 ledger closure | PASS | grand residual −5.8e−16 relative; the bound pool closes separately at −5.9e−12 |
| K12 the restriction is enforced | FIRES | a step at 1.0001× the bound is refused on the frozen geometry |

K1 is the control the specification required without further discussion, and it measures
rather than judges. At equal integrated capacity at every comparison time, moving capacity
from a uniform distribution to one that follows `h(A)` changes the bound inventory by up
to 1.89%. The runs are also checked against each other structurally: since the total
inventory must be identical between them, any difference in the bound pool has to be
mirrored exactly in the dissolved pool, and it is.

K8 is the deliberately omitted transfer, and its verdict is FIRES because red is the
correct outcome. `closure_residual` does not include the `dropped` accumulator, on purpose:
a ledger that books the material a bug destroys is a ledger that cannot report the bug.
The residual comes out at exactly minus the dropped amount, to 6e−11 on an inventory of
5e4.

Every ledger term is accumulated from the term the step computes, never by differencing
the fields, so the closure identity has two independent sides. The Laplacian's
floating-point residue is reported separately and is deliberately not a term in that
identity: folding it in would let a real boundary leak close the books.

## Can the suite fail?

`mutation_controls.jl` plants seven single-line defects in a scratch copy of the module and
runs the data-free suite against each. All seven are caught, and the unmutated baseline is
green. Each patch is checked for having applied exactly once before its suite is run, and
a patch that matched zero times or many times is reported as DID-NOT-APPLY rather than as
a result, because a textual patch that quietly matched nothing is indistinguishable in the
report from a defect the suite failed to catch. That is not a hypothetical: a colour-table
control elsewhere in this repository spent an entire run rewriting a digit inside an XML
attribute name instead of the value it was aimed at, and reported a pass.

The defects are: dropping the neighbour's half of each face flux; halving the stencil
diagonal to the 1-D bound; dropping the binding term from the stability rate; dropping the
decay term from it; destroying the capacity overflow instead of transferring it; booking
the dropped material into the closure identity; and removing the degenerate-modulator
guard. They are caught by 6, 11, 2, 1, 2, 1 and 2 assertions respectively.

## Proposed, not ratified

These are choices recorded as choices. No one ratified them on the operator's behalf, and
an earlier draft of this work was wrong to describe an equivalent list as ratified.

Bound inventory is fixed at lattice sites with no advection. Capacity is volumetric and
finite, one capacity per occupied site, independent of how much of the site any label
claims. Overflow above capacity, `max(b − B, 0)`, is released locally to the dissolved
pool, recorded, and tested by a forced capacity decrease. Material removed by decay is
accumulated in a loss ledger and does not occupy capacity. The capacity modulator is the
6-connected occupied face-neighbour fraction.

`clamp_dissolved = true`, used only by K3, is a different scheme and not a view of this
one: holding `c` fixed makes the system open. The material the clamp returns is booked as
`chemostat_input` so that the ledger still closes, and the production run never sets it.

## Not done here

Recomputing the inert signal itself on these frozen labels is a separate scenario owned by
`diagnostics/inert_signal`, which lives on another branch and is not imported. What this
branch does compute on the frozen labels is this benchmark, and it carries its own parent
hash, in its own receipt, for exactly that reason.

## Files

| | |
|---|---|
| `BindingBenchmark.jl` | geometry, capacity, ledger, flux-form Laplacian, the step |
| `benchmark.toml` | the declared configuration; unknown keys are refused |
| `setup.jl` | configuration reading and the parent chain of custody |
| `run.jl` | the production run; writes `benchmark_receipt.json` and `timeseries.csv` |
| `controls.jl` | K1 to K12; writes `control_verification.json` |
| `test_numerics.jl` | 81 data-free assertions |
| `mutation_controls.jl` | seven planted defects; writes `mutation_verification.json` |
