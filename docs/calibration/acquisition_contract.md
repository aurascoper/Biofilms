# Calibration acquisition contract

**Date:** 2026-08-14 · **Status:** baseline condition UNSET, protocols declared ·
**Data:** `data/calibration/` · **Binding:** `docs/calibration/integration_contract.md`

The spatial and material gates each need a measurement campaign, and they need the *same* one.
This document freezes what a calibration specimen is and how it is handled, so that a lattice
pitch and a biofilm density end up describing one biofilm rather than two.

## 1. One baseline condition, declared before anything is collected

Every calibration value inherits the conditions it was measured under. The baseline specimen is
declared once, in `data/calibration/baseline_condition.csv`, and it ships **UNSET** — this is a
protocol contract, not a claim about a specimen that exists.

```
consortium or surrogate      strain identities        growth medium
temperature                  pH                       oxygen condition
substrate / membrane         biofilm age              flow or static
irradiation: none (baseline) domain semantics: microvolume
biosafety level per strain   institutional approval   containment facility
```

**Biosafety follows strains, not species.** Several components are commonly BSL-1 while
*C. neoformans* strains are commonly BSL-2, and collection material under one species name can
carry different assigned levels — so a species list cannot determine a facility class. Current
guidance is protocol-driven risk assessment rather than a level label, making the chain:

```
exact strain IDs -> institutional review -> containment decision -> approved procedures
```

The gate blocks on a missing `institutional_approval_id`. A sensible staging is to validate the
whole measurement pipeline — substrates, imaging, blank subtraction, constant-mass drying, sample
linkage, ingestion, uncertainty — on an approved BSL-1 system first, so that scarce
higher-containment time is not spent discovering that the hydrated membrane blank swamps the
biofilm signal.

**`domain.semantics = "microvolume"` for the first calibration.** Not `representative_segment`
and certainly not `full_apparatus`: the CPM forces `L = 2R` (`R = N/2` recomputed as a local at
seven sites in `biofilms_potts.jl`, with no length field anywhere), and no pitch repairs that
topology. A microvolume claims to be an abstracted volume and nothing more, which is the only
honest claim available today.

**A single-species or surrogate dataset exercises the harness; it does not clear the gate.** The
calibrated pitch and density must ultimately refer to the seven-species target system under its
declared conditions. Anything else is a machinery test, and the gate reports it as `PROVISIONAL`
however good the numbers are.

## 2. The morphology target is a biomass field, not individual cells

This follows directly from the entity semantics. One CPM ID is a **computational biomass parcel**
(`docs/calibration/cpm_entity_semantics.md`), so segmenting individual bacteria, conidia or hyphae
and treating each object as one CPM entity would calibrate the wrong thing.

What is needed is a calibrated 3-D biomass field `B(x,y,z) ∈ [0,1]` — a binary mask or a
segmentation probability / volume fraction. Pitch selection then coarse-grains it,

```
φ_j(a) = (1/V_j) ∫_{V_j} B(x,y,z) dV
```

and compares structure. Because the CPM lattice is binary, a **declared occupancy map** turns
`φ_j` into occupied/unoccupied — threshold or mass-preserving. That mapping is part of the
calibration and cannot stay implicit; the acceptance config refuses while it is unset. See
`docs/calibration/biomass_field_harness.md`.

Object-level morphology (`morphology.py`) remains useful as a **diagnostic** showing that rods and
filaments lie outside the model's representational contract. It is not the pitch selector.

## 3. Acquisition package

**The replicate count is derived, not declared.** `n = 3` is a convention, not a calculation.
Run a pilot of a few independent batches, estimate the between-batch, within-batch and blank
variance components, then compute the required count with
`precision.required_replicates(target_rel_uncertainty, mean_value, var_between_batch,
var_within_batch, n_within, var_blank)`. It also names the **dominant** component, which is what
says where effort belongs — the three answers are different actions:

| Dominant | What buys precision |
|---|---|
| between-batch | more independent cultures; extra fields per coupon buy nothing |
| within-batch | more imaging fields per coupon, which is cheaper than more cultures |
| blank | a lighter or more reproducible substrate, worth more than either replicate |

**Detectability comes first, and is a different question.** If the biofilm mass is small compared
with the uncertainty of the blank subtracted from it, no replicate count rescues the measurement:
replicates shrink the standard error of a mean, they do not recover a signal buried under a
subtraction. `precision.detectability()` refuses, and the fix is more biomass per coupon or a
substrate with a lighter, more reproducible blank.

**Time-resolved stacks are required, not optional.** The dynamic observable is selected —
`biomass_autocorrelation_decay`, with `interface_rearrangement_rate` in reserve — so a
longitudinal subset must be part of the same campaign. Discovering this after a static campaign
would mean repeating the culturing.

The deliverables that matter are not renderings:

```
raw or losslessly processed 3-D image     segmentation / probability volume
physical voxel calibration                sample metadata
segmentation provenance                   held-out designation
```

The **held-out designation** is a column, not an afterthought: a pitch fitted and validated on the
same stacks has not been validated.

## 4. Pairing: one culture batch, two coupons

The density must come from the same biological condition as the imaging. The contract is:

- one coupon is imaged;
- a sibling coupon from the same culture batch is weighed;
- both carry the same `paired_sample_group`;
- matched blanks undergo the identical hydration, rinsing, blotting, drying and weighing sequence.

Six cross-domain columns now exist on both branches' sample tables:

| Column | Why |
|---|---|
| `culture_batch_id` | the unit of biological replication |
| `paired_sample_group` | links an imaged coupon to its weighed sibling |
| `blank_sample_id` | which blank corrects this sample |
| `measurement_order` | detects drift within a weighing session |
| `growth_condition_id` | joins to the frozen baseline condition |
| `medium_batch_id` | medium salts end up in the ash and in the dry mass |

Both gates now block on unpaired samples. Using unrelated microscopy and gravimetry datasets
would give a pitch for one biofilm and a density for another, and nothing downstream would notice.

## 5. Gravimetry protocol

### Hydrated volume, and what it must enclose

Preferably from the calibrated 3-D imaging over a known substrate area,
`V_hydrated = ∫B(x,y,z) dV`. A mean-thickness × area approximation is permitted only with its
assumptions and uncertainty recorded — `volume_method` is a required column for that reason.

**The volume and the mass must describe the same material.** A balance weighs the whole biofilm:
cells, extracellular matrix, interstitial water and retained solutes. A cell-only fluorescence
segmentation encloses none of the last three, so `ρ_wet = m_wet / V_hydrated` computed from that
pair is **systematically high** by whatever fraction of the biofilm is matrix and water. Both
numbers are positive floats of plausible size, so nothing about them reveals the mismatch — the
closest published analogue built its total biovolume from live cells, dead cells *and* EPS before
combining it with dry mass, precisely for this reason.

`segmentation_basis` (imaging) and `volume_basis` (gravimetry) therefore record what the mask
encloses, and the converters refuse a mismatch:

| Basis | Encloses | Usable as V_hydrated |
|---|---|---|
| `cells_only` | stained cells | **no** — omits matrix and pore water |
| `cells_and_matrix` | cells + stained EPS | only with a declared `pore_volume_fraction`, via `V_envelope = V_stained / (1 − p)` |
| `whole_biofilm_envelope` | the envelope, internal voids included | **yes** |

A wet mass divided by a cell-only volume **fails** material calibration; the gate says so and
`wet_bulk_density()` raises.

### Wet mass, after a reproducible surface-water procedure

"Wet mass" is not a measurement until the surface water is defined. These are now separate
required columns rather than one free-text note:

```
drain_orientation      drain_time_s        blot_material
blot_contact_time_s    ambient_temperature_C   time_to_weighing_s
```

### Blank correction is the definition, not a refinement

```
m_wet,biofilm = m_wet,sample+substrate − m_wet,blank_substrate
m_dry,biofilm = m_dry,sample+substrate − m_dry,blank_substrate
```

Three blank kinds, in `data/calibration/materials/blanks.csv`:

| Kind | Catches |
|---|---|
| `dry_substrate` | substrate tare and manufacturing variation |
| `hydrated_substrate` | water retained by the membrane or support |
| `medium_exposed_abiotic` | salts or precipitates deposited without biology |

The hydrated substrate blank is the one most easily skipped and the most costly to skip: a
hydrated membrane retains substantial water, which would otherwise be counted as biofilm water,
inflating `ρ_wet` and deflating `w_water`. The medium-exposed blank separates medium salts from
biomass ash.

### Derived quantities, per replicate before pooling

```
ρ_wet   = m_wet,biofilm / V_hydrated
X_total = m_dry,biofilm / V_hydrated
w_water = 1 − X_total / ρ_wet
```

Computed **per replicate, then pooled**. A ratio formed from a pooled numerator and a pooled
denominator is a different estimator; if that aggregation is intended it must be declared, because
the two disagree whenever volume and mass are correlated across replicates — which they are.

## 6. Composition, for the OpenMC material gate

Density alone does not clear it. A closed wet-bulk elemental composition is required, blended from
water, dry organic biomass, ash/minerals, retained medium salts, and explicit metal loading where
present. The strong path is a measured water fraction, a measured ash fraction, CHNS on dry
biomass, oxygen measured or explicitly by difference, the exact medium recipe, and ICP for any
intentionally modelled metal.

Without elemental analysis, a literature dry-biomass composition may back **bounded provisional
scenarios** (`low_solids` / `central` / `high_solids`), and the export gate will still refuse to
ship it as calibrated — `evidence_basis = proxy` is not exportable, by design.

## 7. What this contract does not do

It declares no baseline condition (the fields ship unset), imports no data, and selects nothing.
Its purpose is that when data arrives, the pairing, the blanks and the protocol detail are already
required — rather than reconstructed afterwards from lab notes, which is when a surface-water
procedure quietly becomes "as usual".
