# Radiotrophic calibration map

Sorts every quantity the coupling contract needs into one of four tiers, on the
evidence audit at revision f4e8dbf (2026-08-15). This document POPULATES NOTHING.
No pitch, density, seconds_per_mcs, hamiltonian_scale or melanin_scale appears
here as a value, and none is implied. Where a tier is assigned it is assigned on
the repository's own gate definitions, not on a judgement about data quality.

## Tier definitions, mapped to the repository status vocabulary

| Tier | Meaning | schema.STATUS token |
|---|---|---|
| `calibratable_now` | Every input exists, or is a declaration the repository is the authority for. No new measurement and no model change. | `ready` / `provisional` |
| `calibratable_after_target_measurement` | The model holds the quantity and the term exists. The gate is a MISSING MEASUREMENT on the target system. | `blocked` — "awaits a measurement" |
| `requires_model_revision` | The gate is a defect or an absence in the MODEL. Excellent data do not move it. | `unsupported_by_current_model` — "awaits a MODEL CHANGE" |
| `unsupported` | The entity is not licensed to exist in its current form. The first task is to justify the distinction it asserts, not to measure a value. | `not_applicable` until justified |

The distinction between the middle two is the repository's own, from
`schema.py`: `unsupported_by_current_model` is "deliberately distinct from
`blocked`. A blocked quantity awaits a measurement; an unsupported one awaits a
MODEL CHANGE."

## Summary

| Item | Tier |
|---|---|
| domain semantics | `calibratable_now` |
| parcel temporal observable | `calibratable_now` |
| pitch (`lattice_pitch_um`) | `calibratable_after_target_measurement` |
| seconds/MCS (`seconds_per_mcs`) | `calibratable_after_target_measurement` |
| Hamiltonian radiation response (`normalization.hamiltonian_scale`) | `calibratable_after_target_measurement` |
| melanin response (`normalization.melanin_scale`) | `calibratable_after_target_measurement` |
| wet density (`density_g_cm3`) | `calibratable_after_target_measurement` |
| hydrated composition (elemental mass fractions) | `calibratable_after_target_measurement` |
| metal-loaded material class | `calibratable_after_target_measurement` |
| growth/survival | `requires_model_revision` |
| `X_total` | `requires_model_revision` |
| `X_red` | `requires_model_revision` |
| melanized material class | `unsupported` |

---

## Tier 1 — `calibratable_now`

### domain semantics

**Classification.** `calibratable_now`. `DOMAIN_SEMANTICS = {microvolume,
representative_segment, full_apparatus, unresolved}` is a closed vocabulary over
a MODELLING CHOICE, and the repository is the authority for it.
`evidence_basis=declared` is the correct and honest basis.

**What moves it.** Nothing external. `spatial/report.evaluate()` currently
returns a `PROVISIONAL` blocker reading "domain semantics undeclared — the CPM
forces L = 2R, so a full apparatus cannot be claimed". Writing the declaration
clears that specific blocker.

**Required action.** Declare `microvolume` or `representative_segment`.
`full_apparatus` is refused by construction:
`scale_candidates.enumerate_candidates()` hard-codes
`full_apparatus_compatible=False` for every candidate, "so a real apparatus is
compatible only if its own aspect ratio happens to be 2". Record the declaration
with a `source_id` of the `<TOPIC>_DECL_<YEAR>` form; a `ready` row with an empty
`source_id` is refused by `_provenance_problems()`.

### parcel temporal observable

**Classification.** `calibratable_now`, and already done.
`data/calibration/spatial/time_observable.csv` carries a `# SELECTED:` block and
`time_observable.selectable()` requires exactly one selected row with all three of
`represented_by_model`, `measurable`, `semantically_compatible_with_parcel` true.
The selected observable is `biomass_autocorrelation_decay` — the correlation
length of the segmented biomass field over time.

**What moves it.** Nothing. It is a completed declaration. It moves DOWN only if
the model changes such that the observable is no longer `represented_by_model`.

**Required action.** None for the selection itself. Note the consequence, which is
the whole point of having made the choice: every public time series found in the
audit fails the selection on the SEMANTIC condition, not on timestamps. Hyphal
elongation rate (JRR 2024, ~10 min OD600 and hourly microscopy), normalised mean
grey level over one operator-chosen ROI (ISS 2022, 30 min cadence), colony
diameter (Bland 2022, six-hourly), and mRNA recovery kinetics (GSE80230, 30/60/120
min) all have real physical timestamps and none is parcel-compatible. Timestamps
were never the missing half.

---

## Tier 2 — `calibratable_after_target_measurement`

### pitch (`lattice_pitch_um`)

**Classification.** `calibratable_after_target_measurement`. The model holds a
lattice; the gate is that no target-system stack exists.

**Evidence that would move it up.** Calibrated 3D stacks of the seven-species
consortium under its declared conditions, on a `whole_biofilm_envelope`
segmentation basis, from at least two INDEPENDENT acquisition sessions so that one
can be designated `held_out`.

**Required measurement.** Five things, all of which the public record fails:
1. `has_3d_stacks=true` with a stated physical voxel size in x, y AND z.
   `datasets.problems()` refuses `use_first`/`inspect_first` on
   `has_3d_stacks=false`, and the strongest exact-species candidates found are
   single-plane: SECRaM is 1600 spectral × 63 000 spatial pixels in one focal
   plane (filenames carry the literal token `z0`); the PEDOT:PSS `.oir` files carry
   `Z=1` and `z extent 0.000`.
2. `segmentation_basis=whole_biofilm_envelope`, or `cells_and_matrix` with a
   MEASURED `pore_volume_fraction`. S-BIAD474 is constitutive cytoplasmic GFP with
   no matrix stain anywhere in the study, so it is `cells_only`, which
   `HydratedVolume.require_whole_biofilm()` refuses outright.
3. A designated `held_out` sample. S-BIAD474's only 100 % GFP volumetric condition
   is one acquisition on 2020-02-05 (48 h) and one on 2020-02-06 (72 h);
   `pitch_selection.select_pitch()` refuses with "a pitch fitted and validated on
   the same stacks has not been validated".
4. All five `REQUIRED_THRESHOLDS` and all six `OBSERVABLE_TOLERANCES` declared
   before selection. Both selectors refuse on undeclared thresholds, and
   `evaluate_pitch()` treats an undeclared tolerance as failed.
5. `target_relationship=target_consortium`. B. subtilis and D. radiodurans are
   `exact_species`; the enforcement is in `datasets.problems()` and it is not
   conventional.

Anisotropy is a live constraint: S-BIAD474's 3.87 µm z-step against a 0.869 µm
inferred xy pixel is 4.45:1, which floors any pitch near the z-step. Acquire
isotropically or declare the floor.

### seconds/MCS (`seconds_per_mcs`)

**Classification.** `calibratable_after_target_measurement`. The observable half is
selected; the timestamp half must be measured on the target system.

**Evidence that would move it up.** A live time series of the SEGMENTED BIOMASS
FIELD of the target consortium, at a cadence resolving the autocorrelation decay,
with per-frame physical timestamps.

**Required measurement.** Repeated calibrated 3D acquisition of the same field
under non-perturbing conditions, on the same stacks that fix the pitch — the two
are one acquisition, not two campaigns. Note the standing hazard the audit
surfaced: every candidate time series in the public record is growth-dominated
colony expansion, and the repository has already excluded growth-dominated
observables because the CPM has no nutrient-driven biomass creation. The target
time series must be acquired in a regime where the autocorrelation decay is not
contaminated by net growth, or the contamination must be argued for explicitly
rather than assumed away.

### Hamiltonian radiation response (`normalization.hamiltonian_scale`)

**Classification.** `calibratable_after_target_measurement`. The config key exists
(`config/coupling_template.toml:84-90`, commented out, "no default exists"), the
field is declared required at the `cpm_feedback` stage
(`coupling/biofilm_openmc/config.py:108`), and `response.hamiltonian` is `blocked`
in `data/parameter_provenance.csv:74` — the repository's own token for "awaits a
measurement".

**Evidence that would move it up.** Spatial redistribution of biomass in a
CHARACTERISED radiation gradient. That is the entire mapping constraint, and it is
the one the whole audit failed to satisfy from any source.

**Required measurement.** A dose gradient MAPPED at the specimen plane by physical
dosimeters, varying appreciably across a distance spanning many lattice pitches,
with the observable being BIOMASS REDISTRIBUTION, not areal growth. The audit
found exactly one published gradient of the right shape and it fails on the
observable: Bland 2022 has a real, TLD-anchored 14-fold intra-colony gradient
(85.2 rad at the central plug to 6.1 rad at the growing edge) and measures colony
DIAMETER, which is areal expansion. An areal growth or shielding experiment cannot
identify a directional term. Zhdanova 2004 is the only architecture that could —
directional source, quantified return angle, and *C. sphaerospermum* 3176
confirmed among the responders — and its single remaining blocker is precisely
whether the gradient at the growth front was measured or merely calculated.

**Standing hazard, recorded not resolved.** The declared functional form is
`signal = scale × dose_rate_Gy_s`: linear, modality-independent,
history-independent. Three independent findings contradict all three properties.
Malo 2020: gamma increased colony growth, beta INHIBITED it, alpha responded only
after pre-exposure — same organism, same melanin, three signs. Szarka 2026: MSB
before irradiation and irradiation before MSB give opposite outcomes, so the
response is order-dependent. Yuzon 2023 vs Dadachova 2007: verified dose rates
spanning 250 Gy/min to 0.05 mGy/h, roughly 3×10⁸, in the same literature. A
measurement programme that fits a scale without testing the form will fit a number
to a form the evidence has already falsified. Fit dose rate and cumulative dose
separately, per the standing `response.melanin` note, and test modality
separately.

### melanin response (`normalization.melanin_scale`)

**Classification.** `calibratable_after_target_measurement`. Config key declared
and unpopulated; `response.melanin` is `blocked`.

**Evidence that would move it up.** Quantitative melanin or pigmentation against
dose, dose rate, modality, nutrient state and time, WITH a melanized/non-melanized
contrast. Six axes. The audit's exhaustive result: no single study supplies more
than three, and the two halves — quantitative melanin, and a dose axis — sit in
different papers on different strains under different melanization regimes.

**Required measurement.** One strain, one melanization regime, one nutrient state
held fixed and a second varied, melanin quantified by a mass- or ESR-calibrated
assay at four or more doses, at two or more dose rates, one modality, against an
ISOGENIC non-melanized control.

**Four disqualifiers the design must avoid, each drawn from a real record.**
- Pigmentation as mean grey value is an optical proxy on a rendered image, not a
  melanin quantity (Bland 2022, Sci Rep 2022).
- A between-SPECIES melanized/non-melanized comparison confounds pigment with
  everything else about the two organisms (Bland 2022, *C. cladosporioides* vs
  *P. variotii*).
- Melanin structure is not melanin amount (Khajo 2011: EPR, TBARS and ²¹³Bi
  binding at two doses, no content curve).
- "Melanin" denotes at least three chemically distinct polymers — induced
  DOPA-melanin, constitutive DHN-melanin, allomelanin — with different elemental
  compositions. Name the polymer in the row or the row is meaningless.

**Sign hazard.** Two independent records report pigmentation or the melanin
pathway DECREASING under gamma: Bland 2022 (P < 0.001, both species, while
pigmentation rose under UV) and Jung/Kwon 2016 (LAC1/LAC2 laccase transcript
"gradually decreased" at 1 and 3 kGy, authors concluding "gamma radiation itself
did not trigger melanin formation"). A `melanin_scale` imported from a UV study
would carry the wrong sign. Also note that four of the seven modelled species are
not melanized at all, so this parameter is scoped to CN, CS and AN and must not be
declared consortium-wide.

### wet density (`density_g_cm3`)

**Classification.** `calibratable_after_target_measurement`. The binding sentence
demands `density_basis = wet_bulk` and `integration.evaluate()` refuses on a
`None` value: "a coherent sentence with unmeasured blanks in it is still a
sentence with blanks."

**Evidence that would move it up.** ρ_wet = m_wet / V_hydrated, where V_hydrated
is a `whole_biofilm_envelope` volume and m_wet is blank-corrected.

**Required measurement.** Wet mass and a hydrated envelope volume of the SAME
paired sample group, under a DOCUMENTED SURFACE-WATER REMOVAL PROTOCOL, with a
matched abiotic substrate blank. Every candidate in the audit fails on one of the
three. The single closest public approach is Cortesão 2022 (A. niger, paired wet
and dry filter masses, n=3) and it fails twice: the weighed object is aluminium
foil + filter + colony, so w_water needs two tares and a filter-water term, not
merely a blank; and no volume of any basis was measured. `blank_corrected_mass()`
refuses a missing blank because "blank subtraction defines the biofilm mass", and
`_require_volume()` refuses a bare float because "a bare number cannot say whether
it encloses the same material the balance weighed".

**Refused shortcut, recorded.** The widely cited 1.09 g/cm³ buoyant density with
30 % dry matter is a CELL buoyant density from density-gradient centrifugation. It
excludes the matrix and the interstitial water entirely and overstates the bulk by
roughly the void fraction. It is dimensionally plausible, easy to justify in a
sentence, and wrong. `density_g_cm3` stays blank; blank means not measured and is
never read as zero.

### hydrated composition (elemental mass fractions)

**Classification.** `calibratable_after_target_measurement`.

**Evidence that would move it up.** A CLOSED elemental mass-fraction set on a
wet-bulk basis, summing to 1 within `mixture.TOLERANCE = 1e-6`.

**Required measurement.** Water mass fraction, ash mass fraction and CHNS/O on the
same paired sample group, blended by `mixture.blend()` so that
`sum_k w_k · w_(i|k) = 1`. The tolerance is not pedantry: the error message is
"a missing component would otherwise be absorbed as a silent renormalization".

**Three basis traps the audit caught, all of which would break closure.**
- The Barthel deposit reports chitin per mg of EXTRACTED CELL WALL and glucan per
  mg of LYOPHILIZED WHOLE MYCELIUM, under one shared `[/mg BM]` label in the plate
  files and NO denominator at all in the summary files. Three files, three bases,
  one deposit.
- Wang 1996's 15.4 % is lyophilised-ghost mass over lyophilised-whole-cell mass —
  a dry-cell basis, from a single pair of weighings, at day 10 of a still-rising
  14-day curve, on ATCC 24067 with 1.0 mM L-dopa, against an up-to-eightfold strain
  spread and a 20-fold precursor spread measured in the same paper.
- Genome-scale-model biomass objective functions are DECLARED stoichiometry,
  frequently borrowed from *E. coli*. `export.BLOCKED_EVIDENCE` refuses them
  precisely so shipping one would not "launder it into a measurement".

Converting any dry-basis composition to wet-bulk needs a w_water from the SAME
system under the SAME surface-water removal protocol. That is why the composition
and the gravimetry must come from sibling samples of one culture, not from two
literatures.

### metal-loaded material class (`hydrated_metal_loaded_biofilm`)

**Classification.** `calibratable_after_target_measurement`, conditional on a
class-existence test passing first.

**Evidence that would move it up.** A TRANSPORT-RELEVANT DIFFERENCE — in wet bulk
density, water fraction, dry solids fraction, elemental composition, ash, or metal
loading — measured on the SAME biological system with and without metal, on a
consistent basis. `biofilm_material_calibration.md` §1 requires this before a
second material may be created at all.

**Why it is plausible where the melanized class is not.** Sorbed or precipitated
metal is high-Z mass added to a low-Z hydrated medium; a transport-relevant
difference is physically expected a priori. The audit found the metal in
exact-species biofilms (SERS-mapped biogenic Ag; XRF- and immuno-EM-localised UO₂
in the EPS; XRD-identified uranyl phosphate) and found it NEVER QUANTIFIED as a
mass fraction on a hydrated basis, with no metal-free control biofilm
characterised the same way.

**Required measurement.** ICP-OES/ICP-MS on the same digest that feeds the
elemental analysis, on a metal-amended arm and a metal-free control arm grown
alongside and carried through identically, both on a wet-bulk basis with matched
blanks. `_check_metal_double_count()` will refuse an element determined by
combustion into ash AND assayed separately, so the ash and the metal assay must be
partitioned before blending, not after.

**Representation caveat.** SECRaM shows the metal is spatially heterogeneous at
0.5 µm lateral resolution, and `material_model_kind=explicit_components` is
refused at both the export gate and the binding gate because
"`coupling/biofilm_openmc/materials.py` maps one material class per occupied
voxel" and sub-voxel homogenisation "does not exist". Whether that heterogeneity is
SUB-voxel cannot be stated until a pitch exists. Record the question; do not
pre-empt it.

---

## Tier 3 — `requires_model_revision`

### growth/survival

**Classification.** `requires_model_revision`. Repository status token
`unsupported_by_current_model`, `data/parameter_provenance.csv:76`.

**THIS IS A MODEL FACT, NOT AN EVIDENCE GAP.** The CPM has no baseline growth
dynamics and no automatic division or growth process. There is nothing in the
model for a survival curve or a growth rate to be a parameter OF. The audit found
excellent data for it — D10 ≈ 84 Gy for S. oneidensis under Fricke-dosimetered
X-rays with 3 biological × 9 technical replicates; LD90 366/506/112 Gy for A. niger
N402 across X-rays, He and Fe ions with 6 replicates and traceable ionisation-
chamber dosimetry; a gamma reproductive-death phase separating 69 % germination
from 5.5 % colony formation; D10 ≈ 3000 Gy for E. dermatitidis — and NONE OF IT
MOVES THE TIER. That is the test of whether the rule is honoured, and it is
honoured here.

**What implementing it would entail.** Not a parameter fit. A change to the model's
dynamics, of at least these parts:
1. A biomass creation process. Currently no site is created by growth; the CPM
   moves and reshapes existing parcels. A division or growth rule has to be added
   and its stopping condition declared.
2. A nutrient or resource field, or an explicit declaration that growth is
   resource-unlimited — because a growth rule with no limiting field grows without
   bound, and the repository has already excluded growth-dominated observables on
   exactly this ground.
3. A death or removal process, so survival is representable as something other
   than the absence of growth.
4. A re-derivation of the ENTITY SEMANTICS. `entity_semantics.csv` sets
   `entity_kind=biomass_parcel`, `lineage_semantics=computational`, and
   `integration.coherence_problems()` rejection (4) refuses a `biological` lineage
   on a `biomass_parcel` because "parcel ancestry is not organism ancestry
   (biofilms_potts.jl:1578-1583)". A division process that creates parcels creates
   COMPUTATIONAL descendants, not daughter cells. Any survival curve fitted to it
   is a per-parcel survival, not a per-organism one, and the conversion between the
   two does not exist.
5. A re-selection of the time observable, because `seconds_per_mcs` would then be
   constrained by growth kinetics as well as by relaxation, and
   `time_observable.selectable()` admits exactly one selected row.

Until 1–5 exist, `growth_survival` is `unsupported_by_current_model` and every row
that would carry it reads `candidate_parameter = none`.

### `X_total`

**Classification.** `requires_model_revision`. The gate is a units error in the
model, not missing data.

**The defect.** `biofilms_potts.jl` `compute_radial_biomass` L1341-1342 computes
`X_total = total_cells[i] / counts[i]` — occupied sites per radial bin over interior
sites in that bin, a DIMENSIONLESS OCCUPANCY FRACTION. Its unweighted mean is
assigned at L1445-1446 and L1797 to `RadiolysisParams.X_total`, declared at L1161
as "total dry-mass density (g cm⁻³)", and consumed at L1241 as
`uptake = rp.k_ads * rp.X_total + rp.k_red * rp.X_red` with k in cm³ g⁻¹ s⁻¹. The
product is s⁻¹ only if X is a mass concentration. `Quantity.require()` now refuses
the substitution by name: "site occupancy is a DIMENSIONLESS lattice fraction, not
a biomass concentration. This is the conversion that is wrong in biofilms_potts.jl
L1445/L1797."

**Required model change.** Declare and implement a SPATIAL MAPPING RULE that turns
occupancy into a mass concentration — occupancy × ρ_wet × (1 − w_water), per the
contract chain X_total = ρ_wet · (1 − w_water) pinned by test — and replace the
unweighted radial mean, which under-weights the outer bins of a cylinder, with an
area- or volume-weighted one.

**Ordering.** The model change comes FIRST. `materials/report.evaluate()` returns
`RADIODIALYSIS: BLOCKED` unconditionally and says why: "blocked by a units error,
not by missing data". A measured ρ_wet and w_water are also required, but they
cannot be plumbed into a field that is still a site fraction.

### `X_red`

**Classification.** `requires_model_revision`, and it needs TWO independent fixes
plus a measurement.

**Defect one, units.** Inherits the whole of the `X_total` defect above.

**Defect two, definition.** `X_red = red_cells[i] / counts[i]` counts sites whose
`c.species == SO` and divides by ALL interior sites. It is neither a fraction of
biomass nor a fraction of active reducers. The contract requires
X_red = f_red,dry · X_total where f_red,dry is the ACTIVE-REDUCER dry-mass
fraction, and the model's default `X_red = 0.3` is labelled "(Shewanella proxy)" —
a taxonomic proxy standing in for an activity fraction, which
`active_from_taxonomic()` now refuses outright: "taxonomic abundance is not
functional activity, and assuming they are equal would silently set every present
cell to fully active".

**The evidence that the refusal is right, not merely cautious.** Brown 2015:
Fe(III) reduction MORE THAN DOUBLED in irradiated *S. oneidensis* biomass at 50 Gy
while viability was a small fraction of the non-irradiated control, and the effect
required exogenous riboflavin as a shuttle. Reducing activity is decoupled from
viable cell count, and a fortiori from taxonomic abundance and from CPM site
occupancy. The active fraction is a function of dose and physiological state; a
single static scalar cannot represent it at all.

**Required action, in order.** (1) Fix the units, as for `X_total`. (2) Replace the
species-site-count definition with an explicit f_red,dry field carrying its own
evidence basis. (3) Only then measure f_red,dry — an activity assay (rate of Fe(III)
or metal reduction per unit dry biomass) paired to the same dry mass that feeds
X_total. The audit found f_red,dry unmeasured for every organism and every biofilm
system inspected, and note that even a perfect measurement would not unblock
RADIODIALYSIS on its own.

---

## Tier 4 — `unsupported`

### melanized material class (`hydrated_melanized_biofilm`)

**Classification.** `unsupported`. The vocabulary member exists. What does not
exist is any evidence that a melanized hydrated biofilm DIFFERS from a baseline
one in a transport-relevant way — and the one direct test of that premise found no
difference.

**Why this is not merely "unmeasured".** `biofilm_material_calibration.md` §1:
"Do not create seven materials because there are seven species. Per-species
compositions need evidence that their transport-relevant density or elemental
composition actually differs; otherwise they are seven copies of one number wearing
different labels." Four findings say the difference has not been shown and one says
it was looked for and not found:

1. **Measured null on the premise.** Vasileiou 2020 (Zenodo 3667494, CC-BY,
   13,235,258 bytes): Sr-90 beta transmission through melanin solutions and
   suspensions against a cellulose-nanocrystal control matched on elemental
   composition — "melanin does not provide improved shielding in comparison to
   cellulose from beta-radiation". That is the physically expected result if
   attenuation is governed by areal density and Z, and it strips melanin of special
   status absent a spatial-arrangement argument. (Caveat kept: the measured arm is
   electrons; the 40 kVp photon case in that archive is Geant4 simulation only.)
2. **No elemental composition exists for any fungal melanin polymer** in any record
   in the audit. Wang 1996's CHN on the ghost sums to 32.4 % (30.0 % C + 2.4 % N),
   leaving roughly two thirds unmeasured, plus 6–7 % carbohydrate — and even those
   numbers are unverified at source.
3. **"Melanin" is at least three different polymers** with different elemental
   compositions. A single material class asserts they are one substance.
4. **Mass scale.** The largest melanin fraction found anywhere is 15.4 % of DRY
   CELL mass under maximal artificial induction in one serotype-D strain. On a
   hydrated bulk basis, against a medium that is mostly water, the transport
   perturbation from that fraction is small and has never been computed, let alone
   measured.

**What would move it to `calibratable_after_target_measurement`.** Not a melanin
number. A demonstrated DIFFERENCE: wet bulk density, water fraction, ash, or a
closed elemental mass-fraction set, measured on melanized and non-melanized
biofilm of the SAME system on a consistent basis, with the difference exceeding the
propagated uncertainty. Until that exists, the melanized biofilm is assigned the
same `hydrated_effective_medium` material as the baseline, and the burden sits on
the evidence to move it — "not the other way round".

---

## Standing hazards that cross every tier

- **Modality never collapses.** UV is not ionizing radiation here. Mixed LEO fields,
  ²²⁵Ac (alpha + beta + gamma, distributed in the medium), ¹⁸⁸Re beta at
  0.05 mGy/h, 200 kV unfiltered bremsstrahlung, protons, deuterons, alphas and
  heavy ions are each distinct and none is a monoenergetic gamma benchmark without
  an explicit written reconciliation. No such reconciliation exists anywhere in the
  repository.
- **Dose rate is a separate axis from cumulative dose.** Verified rates in the
  audit span 250 Gy/min to 0.05 mGy/h. One acute point constrains neither a
  dose-rate response nor a chronic one.
- **Phenotype follows the measurement, not the title.** Three isogenic melanin-null
  contrasts in *E. dermatitidis* measured melanin's effect on survival and found it
  null under gamma and under alpha, proton and deuteron exposure. Peer review
  removed "radiotrophic" from the *C. sphaerospermum* ISS title between preprint
  and publication. `biofilms_potts.jl:22` and `:24` still carry "radiotrophic" as
  unsourced inline comments, and `entity_semantics.csv` repeats "Radiotrophic
  melanized yeast" as `evidence_basis=declared`. That is a documentation finding,
  reported here and not fixed.
- **A single-species dataset never clears the target gate,** however good it is.
  `clears_target_gate=true` requires `target_relationship=target_consortium` and
  `datasets.problems()` enforces it.
