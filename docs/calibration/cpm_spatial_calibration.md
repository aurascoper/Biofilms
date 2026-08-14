# CPM spatial calibration: identifiability, domain, and representability

**Date:** 2026-08-14 · **Gate:** `PROVISIONAL` ·
**Entity semantics:** `docs/calibration/cpm_entity_semantics.md` ·
**Data:** `data/calibration/spatial/` · **Harness:** `calibration/biofilm_calibration/spatial/`

This branch establishes what would have to be true to choose a CPM lattice pitch, and what the
current model can represent once one is chosen. It selects no pitch, because the evidence to
select one does not exist and because — as set out below — a pitch is not identifiable from the
kind of measurement usually offered for it.

## 1. The pitch is not identifiable from a cell volume

A measured physical volume constrains only a product:

```
V_physical = V_sites · a³
```

so the lattice pitch `a` and the site count `V_sites` are **not separately identifiable** from a
volume measurement. Dividing one reported cell volume by the current `V_target = 120` and taking
the cube root looks like a calibration but is not one: it fixes `V_sites` by accident — at a
value that was never chosen for physical reasons — and then derives `a` from that accident.

Breaking the degeneracy needs a second, independent constraint. Exactly one of:

1. **a declared physical domain size**, which fixes `a = domain / N`;
2. **a resolution criterion** — "the smallest dimension of interest spans at least *k* voxels" —
   which puts a ceiling on `a`;
3. **a declared target site count**, which is a modelling choice and must be labelled as one.

`scale_candidates.enumerate_candidates()` produces the family and reports what each member
implies; `select()` raises unless acceptance thresholds are declared. There is a test that it
raises.

## 2. The domain constraint: L = 2R, always

The Julia geometry is an `N × N × N` cube containing a cylinder of radius `N/2`. `R = N/2.0` is
recomputed as a **local variable at seven separate sites** in `biofilms_potts.jl`; there is no
length, height or z-extent field in `CPMParams` or anywhere else. Under one isotropic pitch:

```
R_physical = (N/2)·a        L_physical = N·a        L = 2R
```

An apparatus with an independently chosen radius and length cannot be represented. This is a
**topology fact, not a resolution one** — no pitch fixes it, and no amount of data changes it.

The consequence for the reference-system programme is direct: **Reference C (the coaxial MABR)
cannot be modelled as a full-length apparatus in the CPM.** Its real aspect ratio is far from 2.
Three futures are permitted, and `config/cpm_spatial_acceptance_template.toml` requires one to be
declared:

| `domain.semantics` | Meaning |
|---|---|
| `representative_segment` | a local axial slice of a longer apparatus, with the axial boundary interpretation stated |
| `microvolume` | an abstracted volume that does not claim to be any apparatus |
| `full_apparatus` | only honest where the real apparatus happens to have L = 2R |

Anything else is `MODEL_TOPOLOGY_CHANGE_REQUIRED`: an anisotropic lattice or an independent axial
extent, which is a change to the Julia model and out of scope here.

Note the asymmetry worth remembering: the **Python** side already separates these —
`coupling/biofilm_openmc/config.py` carries independent `cylinder_radius_cm` and
`cylinder_length_cm`, and the builder only enforces containment. The constraint is Julia-side.

## 3. Resolvable is not maintainable

Two questions get confused in mesh discussions and must not be here.

**Resolvable** — can a shape be written down at this pitch? Arithmetic;
`scale_candidates` answers it via `minor_axis_voxels`.

**Maintainable** — will the model *keep* the shape? `compute_delta_H` (L482–543) is
`ΔH_adh + ΔH_vol + ΔH_rad + ΔH_mel`. No surface-area, elongation, aspect-ratio or connectivity
term exists, so adhesion and the volume constraint together minimise surface at fixed volume and
every object relaxes toward an isotropic blob.

| Shape class | Maintainable | Why |
|---|---|---|
| sphere | yes | the CPM's attractor |
| biomass parcel | yes | the same attractor, and what the entity semantics declare |
| ellipsoid | no | mild anisotropy decays; nothing rewards an elongated state |
| rod | no | needs a surface-area or elongation term — *B. subtilis* is a motile rod |
| filament | no | needs elongation **and** connectivity — *C. sphaerospermum*, *A. niger* |
| disconnected | no | no connectivity constraint; a split object is neither prevented nor penalised |

So a finer pitch buys resolution and buys nothing at all for shape. **Three of the seven species
have morphologies the model cannot hold at any pitch.** That is a recommendation for a future
model revision, not a licence to add Hamiltonian terms here.

`representability.py` builds each class synthetically and measures it. Those masks do not test
the CPM — they validate that the comparison harness *can detect* the blob collapse, by confirming
`morphology.py` separates a rod from a sphere of equal volume. A harness that reported them as
similar would conceal the model's central limitation.

## 4. The comparison harness

`morphology.py` runs the same metrics over a segmented microscopy stack and over a CPM snapshot's
`cell_id` array: volume, equivalent diameter, principal-axis lengths from the inertia tensor,
aspect ratio, 26-connected components, boundary contact. Axis lengths come from the inertia
tensor rather than a bounding box, because a bounding box is orientation-dependent and would
report the wrong aspect ratio for a diagonally-oriented rod — the exact case that matters.

Comparison is by **distribution**, not by mean: two populations with equal mean volume can be
entirely different biology. `size_distribution_summary` reports quantiles and, separately, how
many objects were excluded for touching the image boundary, since those are truncated and would
bias a size distribution downward.

`occupancy()` is deliberately named for what it is — a dimensionless site fraction, the same
quantity `compute_radial_biomass` produces, and **not** a biomass density. See the material
branch's basis contract.

## 5. Gate

**`PROVISIONAL`.** Entity semantics are declared for all seven species; everything else is
missing.

| Requirement | State |
|---|---|
| entity semantics declared | ✓ all seven species, `biomass_parcel` |
| morphology dataset | ✗ `sample_metadata`, `object_morphology`, `biofilm_structure` all empty |
| representative domain selected | ✗ `domain.semantics` unset |
| accepted pitch candidate | ✗ acceptance thresholds unset; `select()` refuses |
| voxel representation acceptable | ✗ not assessable without a dataset |
| uncertainty and held-out validation | ✗ |
| no reliance on normalized R/L | ✓ nothing here uses the legacy `R=1, L=2` |

`READY_FOR_TIME_CALIBRATION` requires all of them. **Seconds per MCS remains out of scope** and
depends on this branch: `Δt_MCS = a²·S_sim / S_exp` needs the calibrated pitch `a`, so the time
calibration cannot start until a pitch is selected.

## 6. Stop conditions

Selection must not proceed while any of these holds:

- the modelled entity is ambiguous for any species;
- a literal-cell reading is required — that returns `MODEL_REVISION_REQUIRED`, and the gate
  enforces it;
- the smallest dimension of interest is unresolved at the candidate pitch;
- the apparatus geometry is incompatible with L = 2R and no domain semantics is declared;
- only normalized graphics or secondary summaries are available;
- the candidate is chosen because it fits the current N = 60;
- fungal lineage is described biologically while the model holds compact parcels.
