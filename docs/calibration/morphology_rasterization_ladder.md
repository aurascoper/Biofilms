# Rasterization ladder — what a pitch costs, and one observable that never converges

**Date:** 2026-08-16 · **Tier:** S0 · **`target_calibration = false`** ·
**Reference D verdict:** `NOT_EVALUATED` ·
**Driver:** `calibration/scripts/rasterization_ladder.py` ·
**Artifact:** `artifacts/ladder/rasterization_ladder.json` ·
**Acceptance:** `config/reference_d_spatial_acceptance.toml`

The subvoxel study asked how finely the **transport** must be resolved on a given
snapshot. This asks the question underneath it: how finely must the
**morphology** be sampled before the structural observables the spatial gate
declares are trustworthy? A tally can be refined after the fact. A snapshot
cannot — its pitch is a measurement decision, made once, at the microscope.

## Why an analytic object, and why not a CPM sweep

Every builder in `spatial/synthetic.py` takes feature sizes in **voxels**, so
`slab(shape=(16,16,32), thickness_voxels=8)` at a different shape is a different
slab. A sweep built from those moves the object and the sampling together and
cannot separate them. `PhysicalSlab` and `PhysicalSpheres` carry the object in
micrometres and rasterize at whatever pitch is asked.

A CPM sweep would be worse than uninformative. `V_target` in `biofilms_potts.jl`
is in **sites**, so raising N at fixed `V_target` shrinks every parcel: the
object changes with the grid in a way that *looks* like convergence. That study
needs a discretisation-renormalisation contract first — contact energy, volume
penalty, temperature, sweep rate and field diffusion all carry a scale. An
analytic shape has no dynamics to renormalise, which is why it goes first.

Voxels are occupied by **centre** sampling, not by any-overlap: the latter
systematically over-fills by half a voxel at every boundary, which would appear
as a pitch-dependent biovolume bias that is an artifact of the sampling rule
rather than of the resolution.

## The axis-aligned control is exact everywhere

| pitch (µm) | shape | biovolume | porosity | interface | thickness |
|---|---|---|---|---|---|
| 3.2 → 0.2 | 10³ → 160×160×320 | 0 | 0 | 2e-16 | 0 |

A slab has an axis-aligned boundary, so voxelisation loses nothing. This is the
control: any error below is about **curvature**, not about the ladder machinery.

## Spheres: two observables converge, one does not

Nine spheres of radius 5 µm on a 16 µm lattice in a 48 × 48 × 24 µm box.
Relative error against the closed-form truth, with the declared Reference D
tolerance applied:

| pitch (µm) | biovolume | porosity | **interface area** | component size |
|---|---|---|---|---|
| 3.2 | 1.20e-01 ✗ | 1.12e-02 ✗ | **4.06e-01 ✗** | 6.13e-02 ✓ |
| 1.6 | 6.39e-02 ✗ | 5.95e-03 ✗ | **4.99e-01 ✗** | 6.39e-02 ✓ |
| 0.8 | 4.04e-02 ✓ | 3.77e-03 ✗ | **4.67e-01 ✗** | 4.04e-02 ✓ |
| 0.4 | 4.55e-03 ✓ | 4.24e-04 ✓ | **4.79e-01 ✗** | 4.55e-03 ✓ |
| 0.2 | 4.62e-03 ✓ | 4.30e-04 ✓ | **5.10e-01 ✗** | 4.62e-03 ✓ |

Coarsest pitch inside tolerance, and inside it at every finer pitch:

| observable | coarsest passing pitch |
|---|---|
| component size (q50) | 3.2 µm |
| biovolume fraction | 0.8 µm |
| porosity | 0.4 µm |
| **specific interface area** | **never, at any pitch** |

Porosity needing a 2× finer pitch than biovolume is not a surprise — it is the
`derived` tolerance doing its job. Porosity is `1 − biovolume`, so an equal
*relative* tolerance on it is ~20× tighter, and the config carries a separately
derived 0.0025 for exactly that reason.

## The finding: `specific_interface_area` is not an area

Its error does not fall with pitch. It sits at 0.41–0.51 across a 16× range and
is slightly *worse* at the finest pitch. That is not slow convergence; it is a
systematic factor, and the mechanism is exact.

**Counting axis-aligned faces measures ∫|n|₁ dA, not ∫|n|₂ dA.** The ratio to
the true area is therefore `E[|n|₁]` over the surface's orientation
distribution. Measured against three geometries with different orientation
statistics:

| surface | predicted `E[|n|₁]` | measured (pitch 0.4 / 0.2 / 0.1) |
|---|---|---|
| plane, axis-aligned | 1 | 1.0000 / 1.0000 / 1.0000 |
| plane, 45° | √2 = 1.4142 | 1.4001 / 1.4071 / 1.4107 |
| sphere, isotropic | 3·E\|nₓ\| = 3/2 | 1.5096 / 1.5011 / 1.5006 |

Three orientations, three different constants, each matched to within 2% and
converging *toward the constant* rather than toward 1. The bias is a property of
the estimator and the surface's orientation, and **no refinement removes it**.

### What it means for Reference D

`maximum_interface_area_error = 0.10` is declared **`hard`** in
`config/reference_d_spatial_acceptance.toml`. Read as an error against a physical
area, **no pitch can satisfy it for any curved interface** — and a biofilm
surface is nothing but curved interface.

The tolerance is not wrong so much as **denominated in the wrong thing**. What
`select_pitch` actually does is compare the *same estimator* on a reference field
and on a coarsened copy of it, and the ~1.5 factor largely cancels between them.
That is why the V. cholerae pilot's "specific interface area falls 46% at factor
2" stands: it is one estimator at two pitches, not an estimator against an area.

So the estimator is a valid **comparator between fields**, and is not an area.
Three consequences:

1. **The name invites the error.** `specific_interface_area_per_um` reads as a
   physical area density and is a voxel face-count density. Anything comparing it
   to a literature specific surface area is comparing a number to 1.5× itself.
2. **The cancellation is only partial.** Coarsening changes the orientation
   statistics of the discretised surface, so the factor is not identical between
   the two fields being compared — which is part of what the 46% is measuring.
3. **If a physical area is ever needed**, the estimator has to change — a
   marching-cubes or Cauchy–Crofton area — not the pitch and not the tolerance.

## What this does not establish

The objects are invented, so this calibrates the **method** and never the
specimen: it says what a pitch costs a measurement, not what pitch a biofilm
needs. Reference D stays `NOT_EVALUATED`, no pitch is selected, and no tolerance
is changed here — changing a frozen acceptance value requires a new
`acceptance_policy_id`, and this finding is an argument for that decision rather
than the decision itself.

Correlation length carries no closed form for a sphere lattice, and a slab's
thickness *spread* is zero only in the continuum, so neither appears as a truth
value: quoting either would manufacture a target the geometry does not have.

## Next

The open item is whether `maximum_interface_area_error` should be redenominated
against a same-pitch reference, or the estimator replaced. That is a declaration
change and belongs to whoever owns the acceptance policy, with the number above
as the evidence.

The dynamic CPM grid-convergence study remains deferred behind its
renormalisation contract, and nothing here shortens that.
