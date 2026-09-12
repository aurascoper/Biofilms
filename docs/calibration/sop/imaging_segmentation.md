# SOP — confocal acquisition, voxel calibration and the declared segmentation basis

Produces `D-PITCH`'s artifact: time-resolved 3-D stacks with voxel calibration and a declared
segmentation basis, at least one training and one **independent held-out** stack. Lands in
`data/calibration/spatial/{sample_metadata|object_morphology|biofilm_structure}.csv`.

**Constrained by `D-RHOWET`, which is a different requirement.** `D-PITCH` requires that *a*
segmentation basis be declared. `D-RHOWET` decides *which* basis is acceptable: `volume_basis =
whole_biofilm_envelope`, or a declared `pore_volume_fraction`. Those are two requirements and the
index records them separately, because an SOP whose stated requirement does not constrain it is the
covered-but-uncovered case the index exists to prevent.

## The constraint most likely to be violated by a default setting

**A cell-stain segmentation is `cells_only`, and `D-RHOWET` refuses it.** A `cells_only` mask
encloses neither the extracellular matrix nor the interstitial water, so the volume it reports is
smaller than the volume the balance weighed, and ρ_wet comes out high. The mask basis is therefore
declared **before acquisition**, not chosen afterwards from what segments cleanly.

## Preconditions

1. Confocal booking and account. **Not held** — this is a core-facility dependency, recorded in the
   index as `blocked_by` rather than assumed.
2. Stage and voxel calibration current, with the calibration recorded per session.
3. The paired coupon exists and its `paired_sample_group` is assigned, so the stack and the mass
   measurement can be joined afterwards.

## Procedure

**1. Declare the basis.** `whole_biofilm_envelope`, or `pore_volume_fraction` with its value. Record
it in `sample_metadata` before the first stack.

**2. Acquire.** Voxel dimensions recorded, not inferred from the objective's nominal figure. One
training stack and one **independent** held-out stack — `select_pitch()`'s acceptance criterion is
that the coarsest pitch preserves all six declared observables within tolerance on **both**, so a
single stack cannot satisfy it however good it is.

**3. Segment to the declared basis.** If the stain cannot resolve the matrix, the honest outcome is
that this basis is unavailable with this stain, which is a finding about the method rather than a
reason to fall back to `cells_only`.

**4. Where the film is too thick for confocal depth**, OCT is the alternative. It has its own voxel
calibration and its own basis declaration; it does not inherit these.

## What this does not do

It does not select the pitch. `select_pitch()` does that from these stacks, and `D-LATTICEN` is
derived from the pitch and the domain afterwards. This SOP produces the input.
