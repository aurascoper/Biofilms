# Subvoxel tally refinement — the effect metric is not resolution-invariant

**Date:** 2026-08-16 · **Tier:** S0 · **`target_calibration = false`** ·
**Reference D verdict:** `NOT_EVALUATED` ·
**Driver:** `coupling/scripts/subvoxel_refinement.py` ·
**Artifact:** `artifacts/refinement/subvoxel_refinement.json` ·
**Statistic:** `debiased_relative_l2_squared`

The pilot found the effect rising monotonically as the tally mesh **coarsens** —
0.1224 at the lattice pitch to 0.1821 at factor 5 — and could not say whether the
lattice resolution itself was converged, because factor 1 was the finest tally
the code could build. This study goes the other way.

## "Factor 1 is the finest mesh available" was arithmetic, not physics

An earlier report stated it as a property of the model. It was a property of one
line: `resolve_mesh_dimension` returned `base // factor` with no multiplicative
path. An `openmc.RegularMesh` sets its physical extent and its bin count
independently, `biofilm_mesh_extent_cm()` already ignored the dimension, and
`mesh_material_volumes` raytraces whatever mesh it is handed — the rays and the
histories see the material lattice through the CSG, they do not index it.

A finer tally therefore adds **transport resolution and never biological
detail**. The material lattice stays piecewise constant at the CPM pitch, so
subdividing a voxel resolves how dose varies *within one parcel of uniform
material*. That is exactly the discretisation question, and it is emphatically
not the morphology question — a finer snapshot is a separate study with a
separate confound, since `V_target` is in sites and refining N shrinks every
parcel.

## The measurement

Ratios 1×, 2×, 4× subvoxel bins per CPM voxel, fixed snapshot, fixed geometry,
5 replicates, 3.2 × 10⁵ particles × 20 batches per run. Ω_b is defined **once**
at the lattice resolution and upsampled, so every ratio is evaluated over an
identical physical region by construction — 1 631 lattice voxels, 13 048 at 2×,
104 384 at 4×, exactly ×8 and ×64.

| ratio | mesh | Ω_b bins | debiased E² | z | E | sub-voxel share | median rel err |
|---|---|---|---|---|---|---|---|
| 1 | 20³ | 1 631 | 1.5082e-02 | 249 | 0.12281 | 0 by definition | 0.028 |
| 2 | 40³ | 13 048 | 1.5144e-02 | 264 | 0.12306 | 0.161 | 0.079 |
| 4 | 80³ | 104 384 | 1.1141e-02 | 148 | **0.10555** | **0.723** | 0.222 |

Ω_b is purely geometric — biomass volume fraction ≥ 0.5 and positive mass. It
previously also carried a dose floor taken from replicate 0, which an external
review correctly identified as making the region a random set drawn from one of
the replicates the estimator averages over. Removing it changed these numbers by
**≤ 0.02%**, so nothing here rested on it; it was removed because a region that
depends on the estimator needs an argument about magnitude, and one that does not
needs none.

The null control stays null at every resolution — `noise_floor` returns
z = 2.0, 0.2 and −1.3, so nothing here is an artifact of the estimator losing
its footing on smaller bins.

### It is not a statistics artifact

The obvious objection is that 4× refinement leaves each bin 1/64 of the
histories, so the drop is noise. It is not. Across a **16× history range**, with
the per-bin relative error falling from 0.736 to 0.222, the 4× result does not
move:

| particles/run | rel err at 4× | E² at 4× | sub-voxel share |
|---|---|---|---|
| 2 × 10⁴ | 0.736 | 1.1132e-02 | 0.7235 |
| 8 × 10⁴ | 0.431 | 1.1170e-02 | 0.7224 |
| 3.2 × 10⁵ | 0.222 | 1.1140e-02 | 0.7229 |

Stable to 0.3%. The debiased estimator is doing exactly what it was built for:
the noise grows sharply and the estimate does not follow it.

### Conservation, which is also the strongest check that refinement is sound

A tally mesh is a histogram of the same events. With the seed independent of the
ratio, every ratio transports **bit-identical histories**, so block-summing a
finer histogram back to the lattice grid — heating and mass separately, divided
afterwards — must return the coarse value exactly.

> **Correction.** The first version of this section checked only the E² statistic
> and read its agreement (nine significant figures) as proof that the dose field
> was reproduced. It was not. **E² is a ratio in which the mass denominator
> appears in both sums and cancels to high order**, so it agreed to nine figures
> while the underlying dose field differed by **1.5%** — because each ratio was
> raytracing its own mass denominator, and independent raytraces of the same
> geometry differ by up to 3.0% per voxel. The check passed for a different
> reason than the one given for it.

**Every ratio now shares one raytrace**, taken at the finest ratio and coarsened
down. Mass is extensive, so coarsening the finest estimate is exact, every ratio
gets a mutually consistent denominator by construction, and it costs one raytrace
instead of one per ratio. With that, the field-level identity holds:

| quantity | max relative difference, 4× reduced vs native 1× |
|---|---|
| mass | 4.8e-16 |
| heating | 6.0e-16 |
| **dose field** | **5.6e-16** |

The driver now asserts on the **field** and aborts otherwise — the claim worth
making, and one that only became true after the shared-raytrace fix. It covers a
non-conservative remap, an intensive quantity averaged where it should have been
summed, refinement perturbing the transport, and the mass inconsistency that hid
here. **It is checked on every run rather than argued for once.**

### And it immediately found a third thing

At the production history count the identity does not hold to machine precision.
Measured drift is **0 at 1×, 1.5e-13 at 2×, and 3.0e-09 at 4×**. Float
accumulation would put the 8-term and 64-term sums about 8× apart, not four
orders of magnitude, so that is not the explanation — and it is worth saying that
"accumulation" was the first explanation reached for, and was wrong.

It is the **orphaned-energy discard, and it is resolution-dependent.** A sliver
bin whose raytraced volume is zero has its heating discarded. At 4× that sliver
is its own bin and the energy is dropped; at 1× the same energy sits inside a
coarse voxel with plenty of mass and is kept. The localisation is exact: of the
voxels carrying drift above 1% of the maximum there is **exactly one**, and it
holds **28 of its 64 sub-bins massless**, against 19% of the remaining voxels
holding any at all. Mass itself conserves to 4.8e-16, so the denominator is not
implicated.

The tolerance is therefore set at 1e-7 — above a mechanism that is understood and
bounded, five orders of magnitude below the 1.5e-02 that the mass inconsistency
produced. The measured drift is recorded per row as
`lattice_field_drift_vs_ratio1` rather than merely tested, so the cost of the
sliver tolerance is visible rather than assumed.

The headline numbers below are **unchanged** by that fix — 1.52847e-02 /
1.53695e-02 / 1.11320e-02 before and after, to six digits — which is itself the
evidence that E² is insensitive to the mass denominator, and so strengthens the
finding rather than qualifying it.

## What it means

**The relative-L2 effect metric is not resolution-invariant, and refining does
not converge it.** Between the lattice pitch and half-pitch it is flat (+0.2%),
then it falls 14% by quarter-pitch. Decomposing the ratio says why:

| ratio | numerator ‖E[d]‖²_W | denominator ‖E[D₀]‖²_W | E² |
|---|---|---|---|
| 1 → 4 | ×1.55 | **×2.10** | ×0.73 |

Both grow when finer bins expose sub-voxel variance, because a mass-weighted sum
of squares over sub-bins exceeds the same sum over their mass-weighted means by
exactly that variance. But **the baseline field has more sub-voxel structure than
the difference field**, so the denominator grows faster and the ratio falls.
That is physically sensible: D₀ carries the steep gradients of an axial line
source and of every material boundary, while d — the response to a uniform
composition change — is smoother.

**72% of the baseline field's weighted L2 norm lives below the CPM lattice
pitch**, and most of it below half-pitch: the sub-voxel share is 0.16 at 2× and
0.72 at 4×. This is measured cross-replicate, so it is structure and not noise —
the naive single-replicate version reads 1.05 at 4× and is almost pure noise.

The practical consequence is that **a declared effect threshold is only
meaningful alongside a declared tally resolution**, in the same way it is only
meaningful alongside a `metric_id`. Quoting "the effect is 0.12" without saying
at what resolution is quoting a number that moves by 14% within the range this
study covers, and by 48% across the full 4×-to-coarsening-factor-5 span
(0.1055 to 0.1821).

This also unifies the coarsening result: the effect rises monotonically as bins
grow and falls as they shrink, and both are the same mechanism seen from
opposite ends.

## Three limits of the exact-CSG denominator, found here

**Massless bins are mostly correct, and refinement converges them.** They run
11.25% at 1×, 18.8–23.0% at 2× and 25.3% at 4×, approaching the geometric truth
that a cube circumscribing a cylinder is 1 − π/4 = 21.5% corner, plus the axial
ends. Buying rays does not reduce them (5.12 × 10⁶ rays leaves 4× at 25.25%
against 25.38% at 1.6 × 10⁶) because they are genuinely void.

**A few are slivers, and that limit is real.** At 4×, 1–2 bins in 512 000 receive
heating while the raytrace assigns them zero volume — a material sliver too thin
for any ray to strike. Their share of total heating is at most **6.5 × 10⁻⁷**.
`dose.specific_energy_per_source` used to refuse this outright, which refuses a
correct run; it now refuses **by magnitude**, tolerating up to 1e-4 with a
warning and discarding the orphaned energy rather than redistributing it. A
genuine tally/geometry disagreement — a misaligned mesh, a transposed axis —
puts percent-level or more there and still raises, three orders of magnitude
above the tolerance.

**Independent raytraces per resolution are a trap.** The mass denominator is an
O(1/√n) estimate, so raytracing each ratio separately gives each one a different
denominator — up to 3.0% per voxel here — and dose is heating over mass, so that
lands directly in every dose field. One raytrace at the finest ratio, coarsened
down, is both exact and cheaper.

**Ray counts must be chosen, not scaled.** Each ray crosses `dim` bins along its
axis, so `n_samples = P·dim²` gives about P samples per bin and the per-bin
sampling does not fall as dim³. Raising it carelessly is worse than useless:
`material_volumes` lays rays on a regular grid, and at some counts a ray passes
within floating-point epsilon of a lattice plane and OpenMC **aborts the
process** — not a catchable exception. The failure is neither monotone in ray
count nor related to refinement: 20³ bins abort at 2.56 × 10⁶ rays while 80³ bins
succeed at 5.12 × 10⁶. Setting `max_lost_particles` does not help, because
`material_volumes` runs through `openmc.lib.init` and never reads it.

## Cost

60 runs, 3.84 × 10⁸ histories, 385 s wall, 0.19 GB of statepoints. The 4× tally
is 512 000 bins; refinement costs tally memory and statepoint size, not
transport time.

## What this does not establish

Nothing about morphology. The material lattice is unchanged at every ratio, so
this says how finely the *transport* must be resolved on a given snapshot and
says nothing about whether that snapshot resolves the biofilm. Reference D stays
`NOT_EVALUATED`, online feedback stays disabled, `target_calibration = false`,
and δ stays unset — and note that this study strengthens rather than weakens the
argument against picking δ from observed effects, since those effects now
demonstrably depend on a numerical choice.

## Next

The open question this hands forward is **which resolution the metric should be
declared at.** Finer is not automatically better: the metric weights sub-voxel
variance that a parcel-scale model may have no claim to represent, since the CPM
asserts uniform material within a voxel and dose structure below that scale is a
transport detail the biology cannot see. A defensible declaration would fix the
tally at the lattice pitch and carry the 4× value as a bound — which is what the
uncertainty ledger already asks for, `sampling_role = convergence_axis`.

The snapshot-resolution study remains separate and still needs its
discretisation-renormalisation contract before a CPM grid sweep means anything.
