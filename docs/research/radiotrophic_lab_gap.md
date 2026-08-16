# The laboratory gap: the smallest experiment public data cannot replace

## The gap, in one sentence

The public record has every piece and pairs none of them: calibrated 3D biofilm
volumes exist, wet and dry masses exist, ash and CHNS exist, melanin assays exist,
dosimetered irradiations exist — and no public record anywhere pairs a calibrated
hydrated volume, a wet mass, a dry mass and matched blanks for ONE system under
ONE growth condition.

That is not a coverage problem. It is a PAIRING problem, and pairing cannot be
retrieved. It has to be produced, from one culture, by splitting one harvest.

Two independent audits reached the same conclusion by different routes: the
spatial branch recorded `PUBLIC_COMPONENTS_FOUND_BUT_NOT_PAIRED`, and the material
branch found no baseline-versus-loaded pair for any organism. The repository gate
that enforces it is in `materials/report.evaluate()`: any sample with no
`paired_sample_group` is a `PROVISIONAL` blocker, because "an unpaired density
describes a different biofilm from the one the pitch was fitted to".

## Why paired siblings, and what a sibling is

One inoculum. One reactor run. N nominally identical coupons harvested at one
time point. Every coupon carries the same `paired_sample_group` identifier, and
each is destined for a different assay. The assays destroy the sample in different
ways — drying, ashing, digesting, staining — so they cannot all be run on one
coupon. The `paired_sample_group` column exists to make the claim that they
describe the same biofilm auditable rather than assumed.

Two constraints on the split are load-bearing and easy to get wrong.

**ρ_wet needs mass and volume of the same object.** `_require_volume()` refuses a
bare number because "a bare number cannot say whether it encloses the same
material the balance weighed". A sibling volume and a sibling mass do not divide
into a density. The design below therefore images ONE coupon label-free and
non-destructively, then weighs that same coupon; the stained volumetric stacks for
pitch selection come from separate siblings.

**Blanks are siblings too.** `blank_corrected_mass()` refuses a missing blank
because "blank subtraction defines the biofilm mass". An abiotic substrate coupon
and a medium-exposed abiotic coupon go through the identical growth vessel,
identical handling and identical surface-water removal protocol, and are weighed
wet and dry alongside the biofilm coupons. Without them every mass in the campaign
is refused at the gate.

---

## Organism, containment and approvals — read this before the design

**BIOSAFETY DETERMINATION IS NOT MADE IN THIS DOCUMENT.** The table below is a set
of POINTERS to the registries that are authoritative, plus the specific hazards
that make each lookup non-routine. The determination belongs to the institutional
biosafety committee, and the repository already enforces that it precede culturing:
`spatial/report.evaluate()` raises the blocker "has no `institutional_approval_id`
— biosafety follows strains, not species, and the approval must precede culturing".

**Where the repository already records this (commit 2dbc154, "Campaign design:
biosafety by strain, derived replicates, observable selected"):**
- `data/calibration/baseline_condition.csv` — the baseline condition row and its
  approval field.
- `docs/calibration/acquisition_contract.md` §1 — the biosafety-by-strain rule in
  prose.
- `calibration/biofilm_calibration/acquisition.py:69-85` — the code that carries it.
- `calibration/biofilm_calibration/spatial/report.py` — the enforcing blocker.

Follow that precedent exactly: one row per STRAIN with its own
`institutional_approval_id`, never one row per species.

**Registries to consult, in this order.** ABSA International Risk Group Database
(risk-group assignments across national frameworks, including where they disagree);
the depositor's own catalogue entry — ATCC, DSMZ, CBS/Westerdijk, NCTC, NBRC — which
states the biosafety level the supplier will ship under; the applicable national
guidance (US: NIH Guidelines Appendix B and the CDC/NIH BMBL; Canada: PHAC
ePATHogen; EU/Germany: TRBA 460 for fungi and TRBA 466 for bacteria; UK: HSE ACDP
Approved List); and, for any *Brucella*-adjacent name, the CDC laboratory
notification on the *Ochrobactrum* reclassification and the Federal Select Agent
Program list.

| Species (strain) | Expected containment — CONFIRM, do not adopt | Facility | Specific hazard driving the lookup |
|---|---|---|---|
| *Bacillus subtilis* NCIB 3610 | RG1 / BSL-1 | Standard microbiology lab | None material. The biofilm workhorse; isogenic matrix mutants and constitutive fluorophore strains are widely available. |
| *Deinococcus radiodurans* R1 | RG1 / BSL-1 | Standard microbiology lab | Wild-type R1 is reported in the literature as a NON-biofilm former. See the strain note below — this is a science blocker, not a safety one. |
| *Shewanella oneidensis* MR-1 (ATCC 700550) | RG1 / BSL-1 | Standard microbiology lab; anaerobic capability | None material. |
| *Cladosporium sphaerospermum* ATCC 11289 (CBS 2) | RG1 / BSL-1 expected | BSL-1 with a Class II BSC | Filamentous fungus: conidia aerosolise. Handle in a BSC regardless of risk group. Allergen exposure is the realistic risk, not infection. |
| *Aspergillus niger* — STRAIN-DEPENDENT | RG1 in some frameworks, RG2 in others; individual ATCC/DSMZ entries differ | BSL-1 or BSL-2 per the determination; Class II BSC MANDATORY either way | Three separate issues: (a) risk group is not consistent across national frameworks for this species; (b) some strains within the *A. niger* aggregate produce ochratoxin A and/or fumonisins, which is a chemical-hazard question independent of risk group; (c) heavy conidiation makes aerosol control the operative control regardless of RG. Pick a strain with a published catalogue biosafety level and a published mycotoxin profile, and record BOTH. |
| *Cryptococcus neoformans* — STRAIN-DEPENDENT | RG2 / BSL-2 expected for H99 and clinical-lineage strains | BSL-2, Class II BSC, BSL-2 waste stream and trained personnel | A human pathogen; the encapsulated yeast is the infectious form and desiccated cells are inhalable. THIS SPECIES DRIVES THE CONTAINMENT DECISION FOR THE WHOLE CONSORTIUM — the repository already records exactly that in `baseline_condition.csv`. Acapsular avirulent derivatives such as cap67 (B3501 background) reduce but do not eliminate the determination; they are still *C. neoformans* and still require the IBC to rule. |
| *Ochrobactrum intermedium* AM7 | **UNRESOLVED, AND THIS IS THE MOST CONSEQUENTIAL ITEM IN THE TABLE** | Do not procure before the determination | NCBI taxid 94625 is now served as `Brucella intermedia`, lineage "…Brucellaceae; Brucella/Ochrobactrum group; Brucella", with *Ochrobactrum intermedium* retained only as a synonym. Some *Brucella* species are select agents and category A infectious substances; CDC issued a laboratory update (2022-12-19) directing labs to handle organisms identified as *Brucella* in a Class II BSC and refer presumptive isolates to public health laboratories. The reclassification is contested in the literature ("If You're Not Confused, You're Not Paying Attention: *Ochrobactrum* Is Not *Brucella*", J Clin Microbiol 2023, doi:10.1128/jcm.00438-23). An IBC may apply *Brucella* handling rules to a shipment whose paperwork carries the new name. **Carry both names on every record, as the Wangiella/Exophiala mapping is carried.** Separately: no culture-collection accession for AM7 was located in any repository. The strain may simply be unobtainable. |

### Nonpathogenic surrogates that would genuinely serve

Surrogate status is enforced, not conventional — nothing below can clear the target
gate, and `datasets.problems()` refuses `clears_target_gate=true` for anything but
the consortium. These are protocol-development and mechanism substitutes.

- **For the melanized-fungus arm: *Knufia petricola* A95.** A non-pathogenic
  rock-inhabiting black fungus with an existing ISOGENIC pigment-mutant panel —
  Δ*pks1* (DHN-melanin null), Δ*phs1* (carotenoid null) and the double mutant. It is
  the only organism in the whole audit offering a clean isogenic melanin contrast at
  RG1, and it separates melanin from a second pigment class, which no other candidate
  does. It is the correct vehicle for developing the melanin quantification and the
  dose-response protocol before any of that is attempted at BSL-2. Confirm its risk
  group with the depositing collection.
- **For *A. niger*: choose the strain, not a different organism.** An *A. niger*
  strain with a published BSL-1 catalogue listing and a published mycotoxin profile
  is the surrogate. "Biosafety follows strains" is the repository's own rule and it
  cuts in the lab's favour here.
- **For *C. neoformans*: there is no adequate surrogate for the target gate,** and
  the honest consequence is that the consortium campaign is a BSL-2 campaign.
  *E. dermatitidis* is the best-deposited melanized-fungus radiation model in the
  literature but is itself an opportunistic pathogen and not RG1. *K. petricola* is
  the mechanism surrogate. Neither substitutes for CN in the consortium.
- **For AM7: no surrogate, and possibly no strain.** Other *O. intermedium* isolates
  (SDCr-5, BCR400, BB12) carry quantified metal phenotypes and are different isolates
  from different environments with different EPS; transferring a chromate-reductase
  phenotype from a tannery isolate to a power-station isolate is exactly the
  strain-collapse the project forbids. The campaign must decide, and RECORD, whether
  it drops the species from the pilot.

### One strain decision that is a science blocker, not a safety one

Wild-type *D. radiodurans* R1 is reported as a non-biofilm former; the only
published *D. radiodurans* biofilm work uses a recombinant biofilm-forming
derivative (gfp/kanR plasmid, SlpA overexpression, 25 mM Ca²⁺). If the campaign
accepts an engineered derivative, the resulting morphology describes an engineered
biofilm-forming construct under a calcium supplement, not the modelled organism,
and that caveat must travel with every row derived from it. Decide before culturing
and record the decision.

---

## The experiment

**One culture. One harvest. Every assay on siblings of it.** Everything below is
the minimum; nothing is included that does not unlock a named parameter.

### Arms

| Arm | Purpose |
|---|---|
| **Sham** | Handled identically, mounted identically, never irradiated. The control for handling, not for radiation. |
| **Irradiated, dose series** | ≥4 doses, one modality, one dose rate. |
| **Irradiated, dose-rate series** | ≥2 dose rates at one matched cumulative dose. Non-negotiable: the standing rule is that dose rate and cumulative dose are fitted SEPARATELY, and one acute point constrains neither. |
| **Gradient arm** | A dose gradient MAPPED across the specimen plane. This is the only arm that can bear on `hamiltonian_scale`. |
| **Metal-amended / metal-free pair** | Only for the metal-loading question; both grown alongside and carried through identically. |
| **Melanized / non-melanized isogenic pair** | Only for the melanin question; a between-species contrast does not count. |

### Sibling allocation, per group

| Sibling | Assay | Order and constraints |
|---|---|---|
| **A** | Label-free hydrated envelope volume, then wet mass, then dry mass | ONE coupon, in this order. OCT or confocal reflectance gives the `whole_biofilm_envelope` directly, adds no stain mass, and does not fix the specimen. Then the surface-water removal protocol, then the balance, then dry to constant mass at a stated temperature. This coupon alone yields ρ_wet. |
| **B** | Calibrated 3D stained stack | A general biomass channel PLUS a matrix channel — not cells-only. If the basis is `cells_and_matrix`, a `pore_volume_fraction` must be measured, not assumed; a non-penetrating fluorescent tracer in the bulk fluid marks the envelope while a penetrating one marks interstitial water, giving the void directly. `to_whole_biofilm()` refuses `pore_volume_fraction=None` because "it is precisely the quantity that separates a stained-voxel volume from the biofilm envelope". |
| **C** | Time-lapse of the same field as B | For `seconds_per_mcs`. Non-perturbing, per-frame timestamps recorded. |
| **D** | Ash, from sibling A's dried mass | Sequential on the same material, so ash is on a known dry basis. |
| **E** | CHNS/O elemental, from sibling A's dried mass | Sequential likewise. |
| **F** | Melanin quantification | Melanized species only. Mass-calibrated extraction with a stated protocol, plus an ESR double integral as a second axis. NAME THE POLYMER (DHN, DOPA, allomelanin) and the precursor and its concentration. Pair to sibling A's w_water so the fraction can be reported on both a dry-cell and a hydrated-bulk basis. |
| **G** | ICP-OES/ICP-MS metal loading | Metal arm only. On the same digest chemistry as E, with the ash/metal partition declared, or `_check_metal_double_count()` refuses the blend. |
| **H** | Reducing-activity assay | Metal arm only. Rate of Fe(III) or target-metal reduction per unit dry biomass, paired to sibling A's dry mass. This is the only route to f_red,dry. |
| **BLK-1** | Abiotic substrate blank | Same coupon material, same vessel, no inoculum, same medium, same handling, same surface-water removal protocol. Weighed wet and dry. `BLANK_KIND = dry_substrate` / `hydrated_substrate`. |
| **BLK-2** | Medium-exposed abiotic blank | `BLANK_KIND = medium_exposed_abiotic`. Distinguishes retained medium solids from biofilm. |

### The surface-water removal protocol

Write it down before the first coupon is harvested, and apply it identically to
blanks. It is the definition of the wet mass, and the repository names it by
requirement. A workable minimum: drain at a fixed angle for a fixed time; contact
a defined absorbent to the coupon EDGE only, a fixed number of times, for a fixed
duration; weigh within a fixed interval of the last contact. Record all five
numbers. `water_mass_fraction()` refuses dry > wet with "check the drying and
surface-water removal protocols" — that refusal is the protocol's own alarm.

### Dosimetry

Physical, traceable, and MAPPED. Absorbed dose to water with a stated uncertainty,
measured at the specimen plane by TLD or OSL at ≥3 positions (≥5 across the
gradient arm), with the field itself calibrated against a traceable ionisation
chamber. Report dose and dose rate as separate quantities with separate
uncertainties. A nominal facility setting is not dosimetry: the audit's records
include an experiment stating outright that "dosimetric data was not obtained", a
transcriptome series naming an isotope and nothing else, a MicroShield model whose
centre TLD read double its prediction, and a beta field characterised by analogy to
Chernobyl. Every one of those is unusable as a dose axis.

### Precision, before replicates

`precision.detectability()` returns `detectable=False` when
`biofilm_mass_g / blank_uncertainty_g < 10`, and states why: "This is a SUBSTRATE
GEOMETRY or BIOMASS problem, not a replicate-count problem: replicates shrink the
standard error of a mean, they do not recover a signal buried under a subtraction."
Choose the coupon area, material and roughness so the biofilm mass clears the blank
uncertainty by an order of magnitude BEFORE choosing n. Then use
`required_replicates()` to set n, and `uncertainty.propagate()` with a declared
seed, because "an uncertainty interval that changes between runs cannot be
reviewed, cited, or regression-tested".

---

## What each measurement unlocks, by parameter name

| Measurement | Unlocks | Gate it clears |
|---|---|---|
| Label-free hydrated envelope volume (sibling A) | The denominator of `density_g_cm3` and of `X_total` | `HydratedVolume.require_whole_biofilm()` — only `whole_biofilm_envelope` passes |
| Wet mass, blank-corrected (A + BLK-1/2) | `density_g_cm3` (`density_basis = wet_bulk`) | `blank_corrected_mass()`; `wet_bulk_density()`; the binding gate's density-basis check |
| Dry mass, blank-corrected (A) | `w_water`; `X_total` = m_dry / V_hydrated | `water_mass_fraction()`; `dry_biomass_concentration()` |
| Ash (D) | Ash mass fraction; the inorganic component of the blend | `ash_mass_fraction()`; `mixture.blend()` closure |
| CHNS/O (E) | The elemental mass-fraction set, closed to 1 on a wet-bulk basis | `mixture.blend()` (`TOLERANCE = 1e-6`); `export.problems_for()`; `to_config_fragment()` |
| Calibrated 3D stacks, ≥2 independent sessions, one designated held-out (B) | `lattice_pitch_um`; the occupancy mapping (`threshold_tau`, or `mass_preserving` with a declared seed); `object_morphology` and `biofilm_structure` | `pitch_selection.select_pitch()` (tolerances + held-out); `scale_candidates.select()` (thresholds); `field.apply_occupancy()` |
| Pore volume fraction, measured (B) | Converts a stained-voxel volume to an envelope volume | `to_whole_biofilm()` |
| Time-lapse of the segmented biomass field (C) | `seconds_per_mcs`, via the already-selected `biomass_autocorrelation_decay` | `time_observable.selectable()` — the observable half is already satisfied |
| Melanin vs dose × dose rate × nutrient state × time, isogenic contrast (F) | `normalization.melanin_scale`; `response.melanin` | The `melanin_production` mapping constraint, all six axes |
| Dose gradient MAPPED at the specimen plane + biomass redistribution measured (gradient arm + B) | `normalization.hamiltonian_scale`; `response.hamiltonian` | The `hamiltonian_radiation` mapping constraint — "spatial redistribution in a CHARACTERISED radiation gradient" |
| Traceable dosimetry, dose and rate reported separately | Makes every dose axis above usable at all | Nothing is fitted against a nominal setting |
| ICP-MS on the same digest, metal-amended vs metal-free pair (G) | Tests whether `hydrated_metal_loaded_biofilm` differs transport-relevantly from baseline; supplies the metal mass fraction if it does | `biofilm_material_calibration.md` §1 (the class-existence test); `_check_metal_double_count()` |
| Reducing activity per unit dry biomass (H) | `f_red,dry`, hence `X_red = f_red,dry · X_total` — AFTER the model units fix | `active_from_taxonomic()` refuses the taxonomic substitution this replaces |
| Blanks, wet and dry, through the identical protocol (BLK-1, BLK-2) | Every mass above. Without them every mass is refused. | `blank_corrected_mass()` |

## What this experiment does NOT unlock, and will not

- **`growth_survival`.** It stays `unsupported_by_current_model` no matter how good
  the survival data are. The gate is the absence of baseline growth dynamics and an
  automatic division process in the CPM. See the calibration map, tier 3.
- **`X_total` and `X_red` as model fields.** The measurements above supply the
  numerators and denominators. They cannot be plumbed into fields that still hold a
  dimensionless site fraction. The model change comes first.
- **The target gate, if the pilot runs on a surrogate or on a BSL-1 subset.**
  `clears_target_gate` can only ever be true for the seven-species consortium under
  its declared conditions. A pilot validates the PROTOCOL — that the pairing closes,
  that the blanks are adequate, that the substrate geometry clears detectability —
  and it clears no parameter gate. Say so in the row.
- **An OpenMC photon-transport validation.** Nothing here transports anything. A
  material composition is an input to transport, not a validation of it.

## Order of operations

1. **Declare, before any culture exists:** the strain list with per-strain
   `institutional_approval_id`; domain semantics; the surface-water removal protocol;
   the acceptance thresholds and observable tolerances; the occupancy mapping and, if
   `mass_preserving`, its seed; the uncertainty-propagation seed. Every one of these
   is a declaration the repository is the authority for, and every one is a refusal
   condition if left blank at evaluation time.
2. **Run the detectability check on the intended substrate geometry.** If the ratio
   is under 10, change the coupon, not the replicate count.
3. **Grow, harvest, split, assay** — in the sibling order above, blanks in parallel
   at every step.
4. **Write the rows with blanks where nothing was measured.** `_parse_numeric`:
   "Blank -> None. Never zero: an unmeasured quantity is not a measured zero, and
   silently conflating them is how a calibration ledger lies."
