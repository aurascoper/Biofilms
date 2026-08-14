# What one CPM cell ID represents

**Date:** 2026-08-14 · **Registry:** `data/calibration/spatial/entity_semantics.csv` ·
**Status:** declared

Before any lattice pitch can be chosen, the model has to say what it is a model *of*. A pitch
converts sites to micrometres; it cannot tell you whether the thing occupying those sites is a
bacterium, a hyphal compartment, or a parcel of biomass. This document answers that, and the
answer is constrained by what the code actually does.

## Declaration

> **One CPM cell ID is a computational biomass parcel.** `lineage_id` is parcel ancestry.
> `generation` counts parcel splits. No output may describe either as a biological generation.

Machine-readable in `entity_semantics.csv`, one row per species, with
`literal_generation_claim_allowed = false` for all seven. The schema rejects any row that claims
biological lineage semantics while the model cannot support it.

## Why the alternative is not available

This is not a modelling preference. Four properties of `biofilms_potts.jl` rule out a literal-cell
reading, and the file already says so about the fourth.

**One target volume for seven species.** `V_target::Int = 120` (L58) is a global scalar used in
exactly three executable lines (L509, L515, and a printf at L1067). There is no per-species
volume, radius, or morphology parameter anywhere — the species vectors that do exist are `β_ion`,
`α_M_species`, `D_s` and `uptake`. A single scalar target cannot simultaneously mean a
~1 µm³ bacterium, a ~65 µm³ yeast cell, and a fungal compartment. Under a literal-cell reading
the same 120 sites would have to be three different physical volumes at once.

**No shape term in the Hamiltonian.** `compute_delta_H` (L482–543) returns exactly
`ΔH_adh + ΔH_vol + ΔH_rad + ΔH_mel`. There is no surface-area, perimeter, elongation,
aspect-ratio, second-moment or connectivity term. Adhesion and the volume constraint together
minimise surface at fixed volume, so **every cell of every species relaxes toward the same
isotropic blob**. *B. subtilis* is described in the registry as a motile rod and
*C. sphaerospermum* and *A. niger* as filamentous fungi, but nothing in the energy function can
hold a rod elongated or a hypha extended. Rods and hyphae are not representable at any pitch.

**Founders are spheres, and not even at target.** `place_cell!` (L323–343) fills a sphere of
radius drawn from `2:3` (L380) — roughly 33 sites at r=2 and 123 at r=3, against a target of 120.
The initial condition is a geometric primitive, not a cell.

**There is no growth, and therefore no generation.** `divide_cell!` (L1687) is called from
nowhere; it is a manually-invoked demonstration primitive with a `min_daughter_volume` rejection
guard, not a trigger. There is no cell-cycle model, no age-to-division rule, and no path from
nutrient to volume — `state.nutrient` is written but never read by `compute_delta_H`. The file
states the position directly at L1578–1583:

> `divide_cell!` is a MANUALLY INVOKED demonstration primitive. There is deliberately no
> automatic division scheduler: the CPM has no growth dynamics from which biological generations
> could be claimed (audit §2).

So the repository already holds this position in code. This document promotes it from a comment
to a checked contract.

## What the declaration permits and forbids

**Permitted.** Calibrating a lattice pitch so that one parcel corresponds to a stated physical
volume of hydrated biomass; comparing simulated parcel-size and spatial-occupancy distributions
against measured biofilm structure; tracking parcel ancestry; reporting radiation dose per parcel
and per parcel lineage.

**Forbidden.** Describing `generation` as a biological generation or a doubling; reporting a
division rate, growth rate or lag phase; interpreting a parcel's volume distribution as a
cell-size distribution; comparing parcel counts to colony-forming units; attributing a lineage to
a founding organism.

## What would change the declaration

A literal-cell reading needs all of:

1. per-species target volumes — a one-field change in principle (`Int` → `Vector`) plus the two
   call sites at L509/L515, but only meaningful once (2) exists;
2. a shape constraint in the Hamiltonian, so that a rod stays a rod (surface-area or
   elongation term), and a connectivity term for hyphae;
3. an automatic division rule tied to growth, so `generation` means something;
4. per-species morphology measurements to calibrate all of the above against.

Items 2 and 3 are model changes, not calibration. If literal biological genealogy becomes a
required claim, the spatial gate returns `MODEL_REVISION_REQUIRED` and this branch stops rather
than calibrating a representation that cannot carry the claim.

## Consequence for the material branch

A parcel is not a cell, so "occupied voxel" cannot mean "cell interior". The material assigned to
an occupied voxel must be a **bulk hydrated biofilm** — biomass, its water, and its extracellular
material together — not a cytoplasm composition. The two branches meet on exactly this sentence,
and the integration branch is where it gets completed:

> One occupied CPM voxel of physical volume `a³` represents part of a hydrated computational
> biomass parcel, and is assigned OpenMC material ______ at density ______.
