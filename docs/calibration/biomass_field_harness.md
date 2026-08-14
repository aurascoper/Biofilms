# The biomass-field harness: how a lattice pitch gets selected

**Date:** 2026-08-14 · **Code:** `calibration/biofilm_calibration/spatial/{field,structure,pitch_selection,synthetic}.py` ·
**Acceptance:** `config/cpm_spatial_acceptance_template.toml`

This is the pitch-selection path. It replaces the intuition that you calibrate a lattice by
measuring cells, which is a category error under the declared entity semantics.

## 1. Why the target is a field

One CPM ID is a **computational biomass parcel** (`docs/calibration/cpm_entity_semantics.md`).
Segmenting individual bacteria, conidia or hyphae and treating each object as one CPM entity
would calibrate a correspondence the model does not claim. What must be preserved under
coarsening is the **structure of the biomass distribution**, not the shape of any organism.

So the input is a calibrated 3-D biomass field `B(x,y,z) ∈ [0,1]` — a binary mask or a
segmentation probability / volume fraction — and the candidate lattice is obtained by
coarse-graining it:

```
φ_j(a) = (1/V_j) ∫_{V_j} B dV
```

Object-level morphology (`morphology.py`) stays, reframed: it is the **representational-contract
diagnostic** that shows rods and filaments lie outside what the Hamiltonian can hold. It is not
the pitch selector.

## 2. Block mean, not block sum

A biomass volume fraction is **intensive**, so coarse-graining averages.

`coupling/biofilm_openmc/mesh.py:coarsen_field` does the opposite *deliberately* — it block-**sums**,
because deposited energy and mass are extensive. Same concept, opposite operation, and using the
wrong one is off by `factor³`. The two are separate installable packages and the calibration one
must not import the transport one (an undeclared dependency that would also drag in h5py), so the
distinction is restated in `field.py` where it can be seen, and pinned by
`test_a_block_sum_would_be_wrong_by_factor_cubed`.

The property that makes cross-resolution comparison meaningful is that `coarse_grain` **conserves
the total biovolume fraction exactly** — to floating point, at every factor.

## 3. The occupancy map is part of the calibration

The CPM lattice is binary, so `φ` has to become occupied/unoccupied. That mapping is a modelling
decision with measurable consequences and cannot stay implicit:

| Mapping | Conserves | Costs |
|---|---|---|
| `threshold` — occupied iff `φ ≥ τ` | connectivity; deterministic | does **not** conserve total biovolume; the bias depends on τ |
| `mass_preserving` — occupied with probability `φ`, seeded | total biovolume, in expectation | the realisation varies by seed, so observables need replicate handling |

Neither is correct in general. `apply_occupancy()` refuses while `[occupancy] mapping` is unset,
and `mass_preserving` refuses without a declared seed — an unseeded lattice cannot be reproduced
or reviewed. `occupancy_biovolume_error()` reports the bias rather than correcting it, because
the size of that bias is a property of the declared mapping that the calibration record should
carry.

## 4. The declared observables

Field-level, all computed on the field and its occupancy mask:

| Observable | What a coarser pitch does to it |
|---|---|
| biovolume fraction | invariant under coarse-graining, by construction — a control |
| porosity | complement of the above |
| specific biomass–void interface area | falls as the interface is smoothed away |
| two-point correlation length | the characteristic scale; a pitch near it erases structure rather than resolving it |
| thickness distribution | counted as occupied voxels per column, so an internal void reduces thickness rather than being counted as biomass |
| component-size distribution | sizes, not a count: a field that keeps its biovolume while fragmenting has not been preserved |
| chord-length distribution | runs through biomass and void |
| depth profile | mean biomass per depth slice |

Two exclusions are deliberate and are the same rule twice. **The array boundary is not an
interface**, and **an edge-truncated chord is not a chord** — biomass reaching the edge of the
field of view is cut off by the *image*, not by the structure. Counting either would make the
metric a property of how large a crop was taken. The first version of
`specific_interface_area` did count the boundary, and a test caught it: the metric changed when
the same slab was cropped.

## 5. The selection rule

> Choose the **coarsest** pitch that preserves every declared observable within tolerance, on
> training stacks **and** on held-out stacks.

Coarsest, because a finer lattice is always more faithful and always more expensive; the question
a calibration answers is how coarse one may go. Held-out, because a pitch fitted and validated on
the same stacks has not been validated.

`select_pitch()` refuses without declared tolerances, refuses when any tolerance is unset, and
refuses when no evaluation is marked held-out. An **undeclared tolerance can never be passed** —
it is not treated as infinitely permissive. When nothing qualifies it returns `None`, which is a
real answer: that field admits no acceptable coarsening at those tolerances.

## 6. Validation is synthetic, and analytic

Fields with known truth, so a failure says the metric is wrong rather than that the data moved:

| Field | Known property |
|---|---|
| `slab` | exact biovolume fraction; exactly one internal face per column |
| `spheres_on_substrate` | exact component count, computed independently of the builder |
| `correlated_field` | correlation length imposed by construction |
| `graded_field` | exact monotone depth profile |

A real stack can be run through the identical path later as a data-only commit; nothing here
changes for that. What a public dataset **cannot** do is clear the gate, unless its organisms,
conditions and observational mapping match the declared target system.

## 7. What this does not do

It selects no pitch — there are no stacks. It does not decide `τ` or the occupancy mapping, does
not set tolerances, and does not touch the spatial gate's verdict, which remains `PROVISIONAL`
for the reasons already recorded. Its purpose is that when a stack arrives, the pitch is chosen
by a rule that was written down first.
