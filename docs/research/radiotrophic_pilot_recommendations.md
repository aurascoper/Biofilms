# Radiotrophic pilot recommendations

Revision f4e8dbf · 2026-08-15 · evidence audit, not a calibration.

Three pilots are recommended. None of them calibrates anything. Each is a bounded,
cheap, verifiable retrieval whose purpose is to make a specific existing repository
claim *attributable* — or to demonstrate, against real bytes, that a gate refuses for
the reason the gate says it refuses.

No pilot in this document populates `lattice_pitch_um`, `density_g_cm3`,
`seconds_per_mcs`, `normalization.hamiltonian_scale`, `normalization.melanin_scale`,
any production material composition, or any target-system biofilm transport
configuration. Every pilot is single-species. A single-species public dataset does not
calibrate the seven-species target consortium, and no pilot below is written as if it
could.

---

## 0. The existing V. cholerae pilot is untouched

`data/calibration/spatial/dataset_candidates.csv` row `vcholerae_dryad_zcrjdfnph`
KEEPS `reference_system_id = public_vcholerae_surrogate` and
`target_calibration = false`. It remains the **generic machinery check** — the pilot
that proves the ND2 read path, the field harness and the pitch-evaluation code run
end to end on a named public 3D reference system that is not one of the seven.

None of the three pilots below weakens it, overwrites it, absorbs it, or supersedes
it. In particular: Pilot 3 is an *exact-species* morphology pilot and Pilot 3's
success would still not retire the V. cholerae row, because the two answer different
questions — V. cholerae asks "does the machinery run", B. subtilis asks "does it run
on a species the model names". A row that passes the first is not evidence for the
second and vice versa. Its `access_status` also stays `browser_manual`: the Dryad
Anubis proof-of-work challenge was deliberately not defeated, that decision stands,
and the recorded SHA-256 digests remain the verification route for a manual copy.

---

## 1. Pilot A — *Cladosporium sphaerospermum* radiotrophic growth-and-shielding pilot

This pilot is named **the radiotrophic growth-and-shielding pilot** and is never to be
called a 3D morphology pilot, an OpenMC photon-dose validation, or a wet-density
calibration. Why, in each case:

- **Not a 3D morphology pilot.** There is no z axis. The observable is a normalised
  mean grey value over one operator-selected ROI in a Raspberry Pi camera image, in
  picture coordinates, with no mm-per-pixel scale anywhere in the article or its
  supplement. `raw_3d_images = false`, `physical_voxel_calibration = false`. There is
  no voxel to calibrate and nothing for `pitch_selection.evaluate_pitch()` to evaluate.
- **Not an OpenMC photon-dose validation.** See §1.6. The field is mixed low-Earth-orbit
  radiation which the authors themselves attribute at this detector to trapped-belt
  electrons of hundreds of keV and protons above 100 MeV, and *"dosimetric data was not
  obtained"* (verbatim). A photon-only transport model cannot represent that field.
- **Not a wet-density calibration.** The lawn was never weighed or harvested. No wet
  mass, dry mass, ash, hydrated volume, water fraction, blank or surface-water removal
  protocol exists. The `~1.67 mm` figure that circulates as a lawn thickness is
  `15 mm − 20 mL/75 cm²`, a geometric residual explicitly described in the supplement
  as *"for the fungal growth layer AND/OR GASEOUS HEADSPACE"* — it may be entirely air.
  `HydratedVolume.require_whole_biofilm()` has nothing to refuse because nothing was
  measured.

The pilot's actual purpose is **corrective**: `biofilms_potts.jl:21-31` labels species 3
a *"radiotrophic filamentous fungus"* in an unsourced inline comment, and
`data/calibration/spatial/entity_semantics.csv` repeats declared radiotrophy language
for CN. This pilot makes the flagship study behind that word citable, and makes visible
that peer review removed "radiotrophic" from its own title between the 2020 preprint and
the 2022 article.

### 1.1 The claim success licenses

> The published *C. sphaerospermum* ISS record measured no dose, no dose rate, no
> melanin and no dose-response, its shielding result is 147 vs 151 CPM at n = 1 with
> p = 0.069, and it therefore supports `none_demonstrated` — so the repository's
> unsourced "radiotrophic" annotation for species 3 is a declared modelling choice
> with no measurement behind it.

### 1.2 The claims NOT licensed

1. **Not `radiotrophic`.** There is no dose axis. The comparison is flight-vs-ground
   with microgravity, transit history, temperature regime and dish-to-dish variance
   fully confounded; the authors state microgravity's role *"cannot be assessed and/or
   gauged with the employed experimental setup"*.
2. **Not `radiation_shielding`.** Withdrawn by audit. 147 vs 151 CPM (2.6%), n = 1 dish,
   one sensor pair, phase-3 p = 0.069 in SI-1 model 3.1 against a proper phase-1
   negative control at p = 0.970. Counts dropped; the effect was not established.
3. **Not `radiation_responsive`.** Granting it would launder an unattributed
   flight-vs-ground growth difference into a radiation response — the same error one
   step smaller.
4. **Not `radioresistant`.** No survival assay, no kill curve, no irradiated-vs-sham
   contrast in this study. (Radioresistance *is* earned for this species elsewhere —
   Vember & Zhdanova 2015 — by a different experiment on different strains.)
5. **Not `radiotropic`.** Growth was areal and uniform; no directional source, no
   gradient, no redistribution.
6. **Not `melanized_radioprotective`.** Melanin was never measured by any method. The
   paper says *"presumably constant"* and *"putatively"*.
7. **Not `melanin_production` / `normalization.melanin_scale`.** Zero of the six axes
   the mapping constraint requires exist: no melanin quantity, no dose, no dose rate,
   no modality, no nutrient-state axis, no melanized/non-melanized contrast.
8. **Not `hamiltonian_radiation` / `normalization.hamiltonian_scale`.** No gradient, no
   directional source, no redistribution, no dose. An areal growth measurement plus a
   below-dish count difference cannot identify a directional term in any form.
   `response.hamiltonian` is `blocked` in the ledger independently.
9. **Not `lattice_pitch_um`.** No stacks, no segmentation, no length calibration.
   `datasets.problems()` refuses any row claiming `use_first`/`inspect_first` with
   `has_3d_stacks=false`.
10. **Not `seconds_per_mcs`.** Physical timestamps exist (30 min imaging, ~95 s
    radiation, ~81 s T/RH), but `time_observable.selectable()` requires
    `represented_by_model` AND `measurable` AND `semantically_compatible_with_parcel`,
    and the first and third fail on a normalised grey value over an operator-chosen ROI.
11. **Not `growth_survival`.** `unsupported_by_current_model` until baseline growth
    dynamics and an automatic division/growth process EXIST in the CPM. The k values
    are real; the model change is what is missing.
12. **Not any material class.** No mass, no volume, no blank, no elemental fractions.
    It is not even a biofilm — a fungal lawn on 13.33 mm of potato dextrose agar, where
    the agar dominates the attenuation path.
13. **Not the target gate.** One axenic culture-collection strain (ATCC 11289 = CBS 2,
    explicitly *not* a Chernobyl isolate) on agar is not the seven-species consortium.
14. **Not the "2.17 ± 0.25 %" figure.** It exists in no version of this work. bioRxiv v1
    says 2.17 ± 0.35 %; v7 says approximately 0.84 %; the published article contains
    neither. Anyone quoting 2.17 ± 0.25 % is quoting a mis-transcription of a superseded
    preprint.
15. **Not a single headline growth ratio.** The abstract says 1.21 ± 0.37-fold, SI-1 §D
    says 1.64-fold ± 0.279 SE on a logistic-slope basis, and the article's own
    exponential k values give 1.24 / 1.29 / 1.53. The basis of 1.21 ± 0.37 is
    undocumented. Record the disagreement; do not average, do not choose.

### 1.3 Data size and access

| Artefact | Size | Access |
|---|---|---|
| Supplement `.s002` XLSX (figshare 20225757) | 13,826,234 B, md5 `04263ef054cb605e6e0e0d83168e0c2d` | `direct_download`, CC BY 4.0 |
| Supplement `.s001` PDF (figshare 20225754) | 995,759 B, md5 `d1ad3fe93fab6fda6b0bd6049c0e84a0` | `direct_download`, CC BY 4.0 |
| **Total** | **14,821,993 B** | |

Route note: the `.s002` DOI redirect chain terminates on the *article landing page*, not
on a file record. Reach the file records through the Figshare API. Supporting records:
Vember & Zhdanova 2015 (`doi:10.15407/microbiolj77.05.070`) is a free 469,611-byte
publisher PDF, `direct_download`, in Ukrainian — this is the record that earns
`radioresistant` for the species. Zhdanova 2004 (`doi:10.1017/s0953756204000966`) is
`inaccessible` behind an Elsevier paywall and was not circumvented; it is the only
located experiment whose architecture could ever license `hamiltonian_radiation`, and it
is an author-request or library task, never a scrape.

### 1.4 Preprocessing

Open the XLSX and read four sheets — `rad`, `growth`, `temp`, `dose`. In `rad`, the
columns are `g_cont, g_cont_n, g_cont_t, g_cont_nt, g_exp, g_exp_n, g_exp_t, g_exp_nt`
against `runtime [h]` / `runtime [d]`. **The `_n` / `_nt` channels are the PocketGeiger
NOISE channels and must never be silently merged with the signal channels** — 23,607
counts per sensor are "radiation and noise" pooled. In `growth`, ground is columns
`AW:CV` and flight is `CV:CZ`. That layout is recoverable from
`github.com/chkern/space-radiation` `src/01_preprocess.R` without downloading anything;
the repo has no LICENSE file and is not redistributable, so cite it by commit SHA if at
all. No unit conversion, no dose derivation: the two count-to-dose conversions in SI-1
§B (≈2.2 nGy per count; ≈0.5 nGy per CPM) disagree by 4.4× and both divide someone
else's whole-body dose rate by this detector's count rate, which assumes the answer.

### 1.5 Repository modules exercised

This pilot exercises **no calibration module**, and that is the finding. It exercises:

- `calibration/biofilm_calibration/schema.py` — a new ledger row must survive
  `read_table()` column-order, vocabulary and `_provenance_problems()` checks. Note the
  collision: `status` values in `NO_VALUE_STATUSES` (`blocked`, `unresolved`,
  `unsupported_by_current_model`) refuse any non-key numeric, so the CPM figures and the
  k values must live in prose cells, not numeric columns.
- `data/calibration/spatial/sources.csv` — one new `*_SCREEN_2026`-pattern source row,
  branch-local, pinning url + document_version + accessed_date + sha256.
- Documentation targets, flagged not fixed: `biofilms_potts.jl:21-31` and
  `data/calibration/spatial/entity_semantics.csv` row CS.

### 1.6 The OpenMC modality reconciliation, and what it would require

A defensible reconciliation between this ISS mixed field and an OpenMC photon-transport
benchmark requires four inputs, **all four of which are absent**:

1. an incident particle spectrum by species and energy at the sensor location — absent;
   only a qualitative South Atlantic Anomaly correlation and a literature annual dose
   equivalent exist;
2. a detector response function for the X100-7 SMD PIN photodiode in that field —
   absent; the only published band is 3–20 keV peak sensitivity, 5.5–10 keV maximum,
   pointed at a field whose relevant primaries are hundreds-of-keV electrons and
   >100 MeV protons, with radiation and noise counts pooled;
3. a mass model of everything between space and the diode — absent; only "2U CubeLab,
   4 × 4 × 8 in, 103.4 in³, air-tight", with no wall material, no Destiny module
   shielding, no dish wall thickness;
4. the attenuator's own areal density and elemental composition — absent; the thickness
   is a subtraction and the mass was never taken.

SI-1 §B names OLTARIS, SPENVIS, GEANT4 and CREME96 as the appropriate tools and states
such analyses are *"out of the scope of this study"*. Supplying any of the four would
mean inventing it. **This is a STOP CONDITION, not a difficulty to be worked around.**

### 1.7 Expected outputs

Ledger rows only. Every row carries, verbatim:

```
reference_system_id  = public_csphaerospermum_iss_prior
target_calibration   = false
```

so that surrogate/prior status survives being copied out of context. Rows:

| output row | `candidate_class` | `model_contribution` | `status` | phenotype |
|---|---|---|---|---|
| `csphaerospermum_iss_frontiers2022` | `target_exact_radiation` | `prior_only` | `unsupported_by_current_model` | `none_demonstrated` |
| `csphaerospermum_beta_shelter_vember2015` | `target_exact_radiation` | `prior_only` | `not_applicable` | `radioresistant` |
| `zhdanova_radiotropism_2004` | `target_exact_radiation` | `prior_only` | `blocked` | `radiotropic` |
| `csphaerospermum_iss_analysis_code_github` | `incompatible` | `protocol_only` | `not_applicable` | `none_demonstrated` |

Note that `target_exact_radiation` here means *exact modelled species, radiation axis* —
it does **not** mean the target consortium, and `target_calibration = false` on every row
is what carries that distinction. `clears_target_gate` can never be `true` for any of them.

### 1.8 STOP CONDITIONS

- **Abandon rather than reconcile the modality.** If anyone proposes supplying the
  incident spectrum, the detector response function, the intervening mass model or the
  attenuator areal density from a handbook, a simulation or an analogous mission —
  stop. That is inventing the measurement.
- **Abandon rather than promote the phenotype.** If a later pass wants
  `radiation_shielding` on p = 0.069 at n = 1, or `radiotrophic` on a flight-vs-ground
  growth ratio, stop. Relabelling radioresistance as radiotrophy is a STOP CONDITION,
  not a judgement call.
- **Abandon rather than transcribe 2.17 ± 0.25 %.** It is not in the record.
- **Abandon rather than turn 1.67 mm into a thickness.** It is 15 mm minus a volume over
  an area, and the paper says it may be headspace.
- **Abandon rather than defeat bioRxiv's 403.** The preprint values stay unverified at
  source. Same discipline as the Dryad Anubis decision.
- **If the strain matters, stop.** ATCC 11289 is a culture-collection strain, not a
  radiation-adapted Chernobyl isolate; the supplement says the Chernobyl isolates
  *"were not easily accessible"*. Any argument that leans on adaptation history fails here.

---

## 2. Pilot B — *Shewanella oneidensis* chemical / material-map pilot

The exact-target chemical/material-map pilot, on the surface-enhanced confocal Raman
microscopy (SECRaM) deposit for a silver-nanoparticle-precipitating MR-1 biofilm.

Its purpose is to make the repository's **refusals** concrete against real bytes: this
is the highest-quality exact-species chemical map available and it still yields zero
terms of `sum_k w_k w_(i|k)` and zero density.

### 2.1 The claim success licenses

> The chemistry of a metal-loaded *S. oneidensis* MR-1 biofilm is measurably
> heterogeneous at 0.5 µm lateral sampling across five separately identified component
> classes and changes between 6, 9 and 35 days — which is an argument that the
> one-material-per-voxel hydrated-effective-medium binding is an approximation whose
> adequacy is testable, and is not a value for any term in it.

### 2.2 The claims NOT licensed

1. **No `w_k` and no `w_(i|k)`.** Raman measures inelastic photon scattering. No
   balance, no blank, no drying endpoint, no ashing, no digest exists anywhere in the
   study. A mass fraction from a spectral intensity needs a per-component scattering
   cross-section and a calibration curve; neither exists.
2. **Intensity ratios are forbidden by the source itself.** Verbatim: SERS *"peak
   intensity greatly depends on the distance from the surface and on the vibrational
   mode orientation … rendering data analysis algorithms that depend on peak intensity
   ratio … not applicable."* Deriving `w_k` from band intensities is precisely such an
   algorithm.
3. **The sampled volume is a ~2 nm shell, not the bulk.** *"The enhancement is effective
   up to ca. 2 nm distance from the metal particle surface."* Surface-selected, non-
   uniformly distributed. Hydrated bulk biofilm composition is a distinct measurement
   basis and this is not it.
4. **The component set is not closed.** Five spectral classes: cytochromes, reduced
   flavins, oxidised flavins, polysaccharides, phosphate. Water — the dominant mass of a
   hydrated biofilm and the whole point of the effective-medium basis — is absent, as are
   bulk protein, lipid, nucleic acid, ash and the silver itself. `mixture.blend()` refuses
   any set not summing to 1 within `TOLERANCE = 1e-6`, *"a missing component would
   otherwise be absorbed as a silent renormalization"*.
5. **No volume.** One focal plane, no z extent. `HydratedVolume.require_whole_biofilm()`
   admits only `whole_biofilm_envelope`; there is no volume of any segmentation basis
   here, so `ρ_wet = m_wet / V_hydrated` has neither numerator nor denominator.
6. **No silver quantification.** No ICP-MS, no ICP-OES, no XRF, no gravimetry, no Ag mass
   fraction, no ash content. SEM-EDX elemental mapping appears in the *primary article*
   at days 7 and 14, but it is qualitative distribution mapping, it is not in this
   deposit, and the deposit is Raman only.
7. **No baseline, and none is possible by construction.** The biofilm is grown *on* the
   Ag/AgCl interface, so no silver-free control biofilm exists in the deposit. Without a
   demonstrated transport-relevant *difference*, `hydrated_metal_loaded_biofilm` is a
   name in a controlled vocabulary with no measurement behind it.
8. **No `lattice_pitch_um`.** No z axis, so no volumetric object morphology, so no
   minor-axis-in-voxels and no volume-quantization-error test is even computable.
9. **No `seconds_per_mcs`.** Three day-spaced ages over 35 days, with the changing
   quantity being SERS band intensity — not a parcel morphology observable, and not
   selectable under `time_observable.selectable()`.
10. **Nothing radiation.** 532 nm at 3 mW is a Raman *probe*, not ionizing radiation and
    not even UV. It must never enter a dose or dose-rate field. No
    `hamiltonian_radiation`, no `melanin_production` — *S. oneidensis* is not melanized.
11. **No `X_red`.** Sorption and distribution are not reduction, and no active-reducer
    dry-mass fraction is defined anywhere. `active_from_taxonomic()` already refuses
    substituting abundance for activity, and the empirical support for that refusal is
    elsewhere (irradiated MR-1 more than doubled Fe(III) reduction at a small fraction of
    viability). `RADIODIALYSIS` stays `BLOCKED` by a units error in
    `biofilms_potts.jl` L1341-1342 / L1445-1446 / L1797, not by missing data — even a
    perfect `f_red,dry` would not unblock it.
12. **Not sub-voxel heterogeneity.** The map resolves 0.5 µm laterally in one plane and
    no lattice pitch has been selected, so the comparison scale does not exist yet. The
    heterogeneity is real; calling it sub-voxel is not yet a claim anyone can make.
13. **Not the target gate.** One species of seven.

### 2.3 Data size and access

Two rows, and the access difference is the load-bearing fact:

- **`shewanella_secram_zenodo_4944335`** — Zenodo record 4944335, carrying the Dryad DOI
  `10.5061/dryad.8sc52` (conceptdoi identical, parent record 4944334). CC0-1.0.
  `access_status = direct_download`: the file endpoint returns HTTP 200 with
  `content-disposition: attachment`, no challenge, no token, no registration.
- **`shewanella_secram_dryad_8sc52`** — same DOI, `access_status = browser_manual`. The
  API download endpoint returns HTTP 401 and the web route is Anubis-fronted. **The
  challenge was not probed and not defeated.** Using the Zenodo mirror is not
  circumvention: it is a second, independent, openly serving host of the same CC0
  content, and all seven MD5 digests match the Dryad manifest byte-for-byte, which is
  the integrity check.

Manifest: 7 files, **5,141,559,882 bytes**. Dryad's `storageSize` reports
5,141,603,486 — 43,604 bytes higher, unreconciled; the manifest sum is the download
budget. Recommended fetch is **1.21 GB, not 5.14 GB**: the three `.WID` scan projects
(403,206,108 / 403,206,475 / 403,206,456) plus the 673,643-byte
`Proxy component spectra.WIP`. The three 1,310,422,400-byte ASCII exports are duplicates
of the same scans and add 3.93 GB for nothing.

### 2.4 Preprocessing

The `.WID` files are WITec project containers; the proxy spectra library pins what each
map channel means (horse-heart cytochrome c, riboflavin phosphate, sodium alginate, PBS,
each in dithionite-reduced and hexacyanoferrate-oxidised variants, with colloidal Ag).

If the ASCII route is taken instead, the array geometry is pinned by arithmetic and not
by assumption: every line is exactly 819,014 bytes (CR at 819,012, LF at 819,013), fields
are fixed 12 characters plus tab so 819,013/13 = 63,001 fields exactly, and
1,310,422,400 / 819,014 = 1600 exactly. **Each export is 1600 spectral channels ×
(1 axis column + 63,000 spatial pixels in ONE plane).** The third axis is spectral. A
reader seeing 5 GB of "biofilm image scans" with a 1600-deep array will mistake this for
a z-stack; the `z0` filename token and the 20×/0.4 objective chosen "to keep as many
bacteria as possible in the focal plane" corroborate that it is not. Do not coarse-grain,
reslice or reinterpret a wavenumber axis into a biovolume.

Deposit ages are 6 d, 9 d, 35 d — a 3-of-5 subset of the article's 1/3/6/9/35 d series,
in which the same location was repeatedly revisited.

### 2.5 Repository modules exercised

The point of this pilot is that these modules **refuse**, in writing, against real data:

- `calibration/biofilm_calibration/materials/mixture.py` — `blend()` on the five spectral
  classes: refused for non-closure.
- `calibration/biofilm_calibration/materials/export.py` — `problems_for()` /
  `to_config_fragment()`: refused for no measured density and no closed element set.
- `calibration/biofilm_calibration/materials/basis_conversion.py` —
  `HydratedVolume.require_whole_biofilm()` and `_require_volume()`: refused, no volume of
  any declared segmentation basis exists.
- `calibration/biofilm_calibration/materials/report.py` — `evaluate()`: `PROVISIONAL` for
  OPENMC on every blocker; `RADIODIALYSIS` unconditionally `BLOCKED`.
- `calibration/biofilm_calibration/spatial/datasets.py` — `problems()` on a row with
  `has_3d_stacks=false`.

### 2.6 Expected outputs

Every row carries, verbatim:

```
reference_system_id  = public_soneidensis_agnp_secram
target_calibration   = false
```

| output row | `candidate_class` | `candidate_material_class` | `model_contribution` | `status` |
|---|---|---|---|---|
| `shewanella_secram_zenodo_4944335` | `target_exact_material` | `hydrated_metal_loaded_biofilm` | `model_revision_evidence` | `ready` |
| `shewanella_secram_dryad_8sc52` | `target_exact_material` | `hydrated_metal_loaded_biofilm` | `model_revision_evidence` | `ready` |

`candidate_parameter = none` on both. Cross-reference the two rows so nobody downloads
twice. Naming a `candidate_material_class` here records *which* class the argument bears
on; it is not a claim that the class is populated, and `candidate_parameter = none` plus
`target_calibration = false` is what carries that.

Two ledger corrections fall out of this pilot and should be applied to the existing
`shewanella_dryad_8sc52` row: its `rationale` cell terminates mid-sentence at
*"Do not use it to select a."* (text lost, almost certainly "a lattice pitch"), and its
`size_bytes=5142000000` with `size_basis=current_version` is a rounded `storageSize`,
not the current-version file total of 5,141,559,882.

### 2.7 STOP CONDITIONS

- **Abandon rather than derive a mass fraction from a band intensity.** The source
  forbids it in writing. This is the single most likely failure mode of this pilot.
- **Abandon rather than treat the 1600 axis as z.** It is spectral.
- **Abandon rather than manufacture a baseline.** No silver-free control exists and none
  can, because the biofilm grows on the silver. A "before" state cannot be reconstructed.
- **Abandon rather than convert a dry or per-EPS sorption capacity onto a wet-bulk
  basis.** No water mass fraction exists for this system.
- **Abandon rather than defeat the Dryad challenge.** If the Zenodo mirror ever
  disappears, the row goes to `browser_manual` and stays there.
- **If someone proposes a second material class on this evidence, stop.** The
  demonstration required is a transport-relevant *difference* in wet bulk density, water
  fraction, dry solids fraction, elemental composition, ash or metal loading, measured on
  the same biological system with and without metal, on a consistent basis. Nothing
  public supplies it for this organism. A second class must be earned by a designed
  measurement campaign, not retrieved.

---

## 3. Pilot C — *Bacillus subtilis* exact-target morphology pilot

The exact-species 3D morphology pilot, on `S-BIAD474` (BioImage Archive), the colony-
biofilm matrix-synergy deposit.

Its purpose is to move the spatial pilot from *"a different organism"* (V. cholerae) to
*"one of the seven"* — a strictly stronger position — while making the held-out failure
visible rather than papering over it.

### 3.1 The claim success licenses

> Calibrated exact-species 3D colony-biofilm stacks exist for *Bacillus subtilis*
> NCIB 3610 at a declared 3.87 µm z-step over an 890 µm field at 48 h and 72 h, so the
> spatial harness can be exercised on a modelled species — and the deposit's single
> 100 %-GFP acquisition per timepoint and cells-only segmentation basis are exactly why
> a pitch still cannot be selected from it.

### 3.2 The claims NOT licensed

1. **Not `lattice_pitch_um`.** `pitch_selection.select_pitch()` raises the held-out gate:
   *"no held-out stack among the evaluations — refusing to select. A pitch fitted and
   validated on the same stacks has not been validated."* The 100 %-GFP volumetric
   condition is **one** acquisition at 48 h (dated 200205) and **one** at 72 h (200206).
   The many independent acquisition dates belong to the mixed-label time-lapse
   components, which are not biomass fields.
2. **Not a pitch even with thresholds declared.** `scale_candidates.select()` refuses
   without `minimum_minor_axis_voxels`, `maximum_volume_quantization_error` and
   `maximum_memory_bytes`; `pitch_selection.select_pitch()` refuses without all six
   `OBSERVABLE_TOLERANCES`. Both are unset and this pilot does not set them —
   *"a pitch chosen because it fits the current N=60 is not a calibration."*
3. **Not a hydrated volume.** Contrast is constitutive cytoplasmic GFP and mKate2 only.
   No matrix stain, no membrane dye, no general biomass stain anywhere in the study;
   matrix is inferred from mutant contrast, never imaged. `segmentation_basis =
   cells_only`, which `HydratedVolume.require_whole_biofilm()` refuses outright —
   *"the volume excludes the extracellular matrix and the interstitial water that the
   balance weighed."*
4. **Not a biomass field from the mixed-label components.** Those are
   *"80 µl non-fluorescent strain with 20 µl fluorescent strain"* after OD
   normalisation — a cell subsample, and "20 % of cells" is itself an inference from OD
   equivalence, not a count.
5. **Not `seconds_per_mcs`, and not `biomass_autocorrelation_decay`.** Every time-lapse
   component is 80:20 mixed-label; the two 100 %-GFP files are static timepoints. There
   is no time-resolved biomass field in this deposit at all, so the growth-contamination
   question does not even arise.
6. **Not a verified voxel calibration.** The 3.87 µm z-step / 1024² / 890 µm figures are
   publication-level, in one methods paragraph covering all 10× GFP colony work; no
   `.lif` was opened. `has_voxel_calibration = unverified`. The 0.869 µm lateral pixel is
   my arithmetic, and the 4.45:1 anisotropy floors any pitch near the z-step and
   quantises vertical morphology statistics at ~4 µm.
7. **Not `density_g_cm3`, `X_total` or any OpenMC material.** No wet mass, no dry mass,
   no blanks, no drying protocol, no elemental analysis.
8. **Not `hamiltonian_radiation` or `melanin_production`.** No radiation of any kind was
   applied; no melanin.
9. **Not `growth_survival`.** `unsupported_by_current_model`.
10. **Not the target gate.** `datasets.problems()`: *"only the target consortium under
    its declared conditions can clear the target gate."* `target_relationship =
    exact_species` is not `target_consortium`, and the phrase "exact species clears the
    named-species spatial gate" means a strictly weaker thing.
11. **Not a licence claim.** The BioStudies record carries `Template`, `Title`,
    `ReleaseDate`, `AttachTo` and **no licence attribute**. The article is CC BY; the
    deposit's terms are unverified and must not be transcribed.
12. **Not a full-apparatus geometry.** `scale_candidates.enumerate_candidates()`
    hard-codes `full_apparatus_compatible=False` for every candidate because L = 2R is
    forced. Declaring `domain_semantics` anything but `microvolume` or
    `representative_segment` would be refused.
13. **Not a held-out set from `S-BSST261`.** That deposit is biofilm *edges* at 24 h and
    5 days on an **unstated medium**; `S-BIAD474` Fig 7 is edge+middle at 48 h and 72 h on
    MSgg. A split across the two tests **condition transfer**, not pitch stability. It can
    serve as an independent-study check and the distinction must be declared in
    `sample_metadata`, never assumed.

### 3.3 Data size and access

Whole study is ~1 TB (one Movie 6 file alone is 63 G). **Fetch one file, 3.0 G**, at the
corrected path — the previously circulated path omits a directory level:

```
https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/474/S-BIAD474/Files/
  Systematic-microscopical-analysis-of-Bacillus-subtilis-colony-biofilm-development/
  Fig%207/200205_100pc-GFP_Edge_Middle_48h/100pc-GFP_Edge_Middle_48h.lif
```

`access_status = direct_download`: unauthenticated EBI FTP/HTTPS, no challenge. The
`Files/` root holds only the 24 filelist manifests plus that one directory. Fig 7 is six
files at 1.9 / 1.8 / 3.7 / 4.3 / **3.0** / 3.2 G.

Second deposit, for the independent-study check: `S-BSST261`, 5 Leica `.lif` files,
**18,529,647,313 B** total; the cheapest inspection point is
`WT_Pulcherrimin_edges_24h_2.lif` at 730,131,371 B. Licence unverified there too. Note
pulcherrimin is an iron-binding pigment, **not melanin** — different molecule, different
pathway, no bearing on the melanin branch.

### 3.4 Preprocessing

Read the `.lif` and answer two questions before anything else: does the container carry
physical voxel size, and does it carry per-frame timestamps. Both are expected and both
are unverified. Then segment the 100 %-GFP channel to a `[0,1]` field —
`field.check_field()` requires 3-D, non-empty, NaN-free, in `[0,1]`, *"a segmentation
probability or volume fraction, not a raw intensity"*. Occupancy mapping must be declared
(`field.apply_occupancy()` refuses an undeclared mapping; `mass_preserving` additionally
requires a seed). Coarse-graining uses block-**MEAN**, not block-sum, and the factor must
divide the shape evenly.

Record `segmentation_basis = cells_only` honestly. Do not convert it with
`to_whole_biofilm()` — that needs a measured `pore_volume_fraction`, *"precisely the
quantity that separates a stained-voxel volume from the biofilm envelope"*, and it does
not exist here.

### 3.5 Repository modules exercised

- `calibration/biofilm_calibration/spatial/field.py` — `check_field()`, `coarse_grain()`,
  `apply_occupancy()`.
- `calibration/biofilm_calibration/spatial/scale_candidates.py` —
  `enumerate_candidates()` runs; `select()` refuses (thresholds unset).
- `calibration/biofilm_calibration/spatial/pitch_selection.py` — `evaluate_pitch()` runs;
  `select_pitch()` refuses (tolerances unset **and** no held-out stack).
- `calibration/biofilm_calibration/spatial/datasets.py` — `problems()`,
  `usable_for_pitch()`.
- `calibration/biofilm_calibration/spatial/report.py` — `evaluate()` returns
  `PROVISIONAL` with the held-out, surrogate, paired-sample and threshold blockers named.
- `calibration/biofilm_calibration/acquisition.py` — the `institutional_approval_id` and
  biosafety-follows-strains rules, for the campaign this pilot is designing.
- `calibration/biofilm_calibration/spatial/representability.py:49-53` — B. subtilis is
  already named there as a morphology the model cannot hold at any pitch.

### 3.6 Expected outputs

Every row carries, verbatim:

```
reference_system_id  = public_bsubtilis_colony_biofilm
target_calibration   = false
```

| output row | `target_relationship` | `verdict` | `clears_target_gate` | `has_3d_stacks` | `has_voxel_calibration` | `status` |
|---|---|---|---|---|---|---|
| `bsubtilis_biad474_matrix_synergy` | `exact_species` | `use_first` | `false` | `true` | `unverified` | `provisional` |
| `bsubtilis_bsst261_pulcherrimin` | `exact_species` | `inspect_first` | `false` | `unverified` | `unverified` | `provisional` |

Plus a `spatial/report.evaluate()` verdict of `PROVISIONAL` with its blocker list quoted
verbatim — that list *is* the deliverable, because it enumerates what the target campaign
must acquire.

### 3.7 STOP CONDITIONS

- **Abandon rather than fit a pitch on one acquisition.** If the pilot's output is a
  number for `lattice_pitch_um`, it has failed. The correct output is
  `select_pitch()` refusing with the held-out message.
- **Abandon rather than declare acceptance thresholds to make the gate pass.** Thresholds
  are a calibration decision with an owner; inventing them to unblock a pilot is the
  N=60 error wearing a biological hat.
- **Abandon rather than treat `cells_only` as `whole_biofilm_envelope`.** No
  `pore_volume_fraction` exists.
- **Abandon rather than call `S-BSST261` a held-out set.** Different condition, unstated
  medium. If a within-condition held-out stack is needed, it must be acquired.
- **Abandon rather than transcribe a licence.** If the BioStudies terms cannot be
  verified from an authoritative page, `license` stays `unverified`.
- **If the heiData soft-X-ray tomography alternative (`doi:10.11588/data/KH6NDD`) is
  reached for, stop and read this first.** It is CC BY 4.0, 31 files,
  13,799,884,261 B, label-free and hydrated — and: 8 of the 31 files are **cells in
  suspension at t = 0, not biofilm**; the deposited TIFF voxel calibration is internally
  inconsistent (`unit=micron` with `spacing=2.105e-07`, i.e. 0.21 pm, which is not a
  credible voxel); the article's field of view is a **per-acquisition 15–17 µm range**, so
  the FOV cannot pin a voxel even in principle; and the API is **user-agent-dependently
  Anubis-fronted** (v1.26.2, the same proof-of-work system as Dryad). It also images at
  530 eV — an imaging probe that deposits dose but measures no radiation endpoint, so it
  is `none_demonstrated`, and it must never be reconciled against the OpenMC photon
  benchmark as an irradiation. Until the voxel question is answered by the authors, every
  length derived from that deposit is unusable.

---

## 4. What none of the three pilots does

- None populates a calibration value. `candidate_parameter` is `none` on every row in
  Pilots A and B, and Pilot C's contribution is a refusal message, not a pitch.
- None clears the target gate. All three carry `target_calibration = false`.
- None reconciles a mixed or non-photon field against OpenMC. The ISS field, the 532 nm
  Raman probe and the 530 eV soft X-ray beam are three different things and none of them
  is a monoenergetic gamma benchmark.
- None defeats an access control. Dryad Anubis, bioRxiv 403, ScienceDirect 403, Wiley
  402, Elsevier paywalls and the heiData proof-of-work were all encountered and none was
  circumvented. Where a browser is the intended public route, the row reads
  `browser_manual` and the block is recorded as a value, not omitted.
- None supplies a wet bulk density, a closed elemental mass-fraction set, or a
  water fraction for any organism. No public record anywhere in this audit pairs a
  calibrated hydrated volume, a wet mass, a dry mass and matched blanks under one growth
  condition for any of the seven species. That measurement must be **designed**, not
  retrieved.
