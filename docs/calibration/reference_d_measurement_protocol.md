# Reference D — the paired measurement protocol

**Frozen:** 2026-08-15 · **Requirements:** `data/calibration/reference_d_requirements.csv` ·
**Status:** `calibration/scripts/reference_d_status.py` · **Target:** D — engineered
radiotropic composite (`docs/physical_reference_system.md` §2)

This is the specification for the campaign that fills Reference D. It exists because "run the
lab campaign" is not one, and because the cost of discovering a missing column after the
coupons are dry is a repeated campaign.

**What this document is not.** It does not restate the contracts. Each requirement points at
the document that already specifies it — `acquisition_contract.md`, `material_basis_contract.md`,
`cpm_spatial_calibration.md`, `time_observable_contract.md`, `integration_contract.md`. What is
here is the part none of those hold: the **ordering**, the **physical procedure**, and the
**declarations that must be made before any data is looked at**.

**Where the authority lives.** Not here. The requirements table is checked against the live
gates and against `emit_transport_config`'s own refusal, in both directions — a field the code
requires with no row is a failure, and so is a row justified by nothing the code refuses on.
Run it before planning anything:

```sh
python calibration/scripts/reference_d_status.py
```

Today: **30 requirements, 10 satisfied.** Spatial gate `PROVISIONAL`, material OpenMC gate
`PROVISIONAL`, radiodialysis gate `BLOCKED`. `CAMPAIGN_READY: no` — blocked by exactly one
thing, `D-APPROVAL`, which only a biosafety committee can issue.

---

## 0. The declarations — MADE, not pending

All twelve modelling decisions are frozen as `reference_d_acceptance_v1`, before any target
data exists, in `config/reference_d_spatial_acceptance.toml` and
`config/reference_d_material_acceptance.toml`. Every value carries its reason on the line
above it, and every threshold carries an `enforcement` class, because a gate must not close
on a criterion nothing compares against.

| Declaration | Value | Requirement |
|---|---|---|
| Domain semantics | `representative_segment` | `D-DOMAIN` |
| Occupancy map | `threshold`, `tau = 0.5` | `D-OCCUPANCY` |
| Source idealization | `line_z_axis`, `isotropic` | `D-SRCSHAPE` |
| Boundaries | radial `vacuum`, axial `reflective` | `D-BOUNDS` |
| Material model kind | `hydrated_effective_medium` | `D-MODELKIND` |
| Five representation thresholds | see config | `D-THRESH` |
| Six observable tolerances | see config | `D-OBSTOL` |
| Four material acceptance thresholds | see config | `D-MATACCEPT` |
| **Institutional biosafety approval** | **outstanding** | `D-APPROVAL` |

**0.1 The one that is still open, and gates culturing itself.** Biosafety follows *strains*,
not species — the consortium mixes BSL-1 organisms with *C. neoformans* at BSL-2 — so
`biosafety_level_by_strain` is per-strain (`'SO:BSL1;CN:BSL2'`), never one level for the
consortium. `institutional_approval_id` must precede culturing, and `is_target_system` must
be `true`: a surrogate exercises the harness and can never clear this gate. Nothing in this
repository can supply it.

**0.2 The occupancy map is a calibration decision, and its cost is measured.** `threshold`
is deterministic and needs no seed, so structural observables need no replicate handling.
It does **not** conserve biovolume — caught by its own tolerance — and it does **not**
preserve connectivity: a majority threshold severs thin bridges and deletes narrow
filaments. That loss is reported by `field.occupancy_topology_report()` rather than assumed
away, and on a two-blob fixture joined by a sub-tau bridge it registers one component
becoming two.

**0.3 The axial boundary is conditional.** `reflective` imposes mirror symmetry, which
represents a longer uniform reactor only if the material near both axial faces supports that
continuation. An arbitrary CPM snapshot may vary strongly along z. The declaration is
therefore stamped `axial_boundary_validation_required = "axial_face_compatibility"` and
`D-AXIALFACE` records the check, which needs a real snapshot. Fallback if it fails: axial
vacuum, scored over a central region distant from the ends.

**0.4 A sealed capsule still cannot be represented.** `line_z_axis` is an idealized
zero-radius axial line with no encapsulation and no self-attenuation. Reference A1 is
declared and unsupported, so if the built source is a capsule, this is a recorded modelling
limitation and not a description of the apparatus.

**0.5 Three declared criteria are decoration, and say so.** `maximum_axis_ratio_error`,
`minimum_objects_per_species` and `require_independent_validation_sample` have no consumer.
They are `enforcement = not_implemented` and may not close a gate. `minimum_objects_per_species`
is additionally questionable under parcel semantics: pitch selection operates on a
coarse-grained field, not on segmented organism objects, and mixed-consortium imaging may
not separate the seven species at all.

---

## 0.6 Institutional authorization — before any of this

Everything below this line is laboratory work on live organisms, and none of it
may begin until an institution says so. That authorization is **evidence
produced by a biosafety process**, not a declaration anyone in this repository
can make. It was once asserted here as the words "d approved"; that string has
no evidentiary force, and the gate now refuses it by value.

The milestone is `REFERENCE_D_CAMPAIGN_INSTITUTIONALLY_AUTHORIZED`, printed by
`calibration/scripts/reference_d_status.py` with the date it judged against:

1. exact strain identities frozen
2. exact growth conditions frozen
3. containment facility identified
4. risk assessment referenced
5. institutional biosafety approval issued, with an issuing authority
6. approval scope matches the strains, procedures, facility and conditions
7. approval effective before culturing began, and not expired
8. approval provenance recorded and retrievable
9. `baseline_condition.csv` carries the approved target row

`CAMPAIGN_READY` requires it, so the two cannot disagree. A tenth criterion of
the form "the status script agrees" would be circular and is deliberately absent.

**BIOSAFETY FOLLOWS STRAINS, NOT SPECIES**, which is why criterion 1 is first:
*C. neoformans* strains are commonly BSL-2 where the other six are commonly
BSL-1, so a species list cannot determine a containment class.

**The approval expires**, so this verdict depends on when it is asked. That is
correct — a lapsed approval authorises nothing — and the evaluation date is
printed rather than assumed.

**What the checks cannot do.** They verify the fields a reviewer would ask for:
that the identifier is not filler, that an authority is named, that the document
resolves and is the right *kind* of document, that the dates bracket the
culturing, and that the conditions on the row still match what was approved.
They cannot verify the document is genuine. No check local to this repository
can. The scope digest detects **drift, not forgery** — it makes an edited growth
condition visible, which is the failure that happens by accident.

Preparatory work belongs in `data/calibration/reference_d_condition_proposal.csv`
(`authoritative_for_campaign = false`) and
`docs/calibration/reference_d_ibc_submission_draft.md`.
`baseline_condition.csv` stays at zero rows until an approval exists: it is an
authoritative input, not a planning worksheet.

---

## 1. The pairing discipline, which is the whole design

Everything else in this protocol is ordinary laboratory work. This is the part that is easy to
get wrong and impossible to repair afterwards.

> **The imaged stack and the weighed coupon must be the same biofilm, and the record must
> prove it.**

A pitch fitted to one biofilm and a density measured on another are two numbers about two
systems, and the voxel binding joins them into one sentence regardless. `paired_sample_group`
is what makes the join true, and it is enforced on **both** sides — the spatial gate and the
material OpenMC gate each refuse an unpaired sample independently.

The six linkage columns travel with every sample row and are deliberately *optional at the
schema and enforced at the gates*, so a table can honestly record an unpaired pilot sample
while the gate still refuses to calibrate from it:

`culture_batch_id` · `paired_sample_group` · `blank_sample_id` · `measurement_order` ·
`growth_condition_id` · `medium_batch_id`

At least one imaged sample must carry `held_out = true`. **Held out means an independent
sample** — not another frame, channel or crop of the same file. This is exactly what stopped
the *V. cholerae* pilot: one 3-D stack cannot validate a pitch fitted to itself, and
`select_pitch()` refuses without a held-out evaluation.

---

## 2. Imaging — what the pitch is fitted to

**Requirements:** `D-PITCH`, `D-TIMESERIES` → `data/calibration/spatial/{sample_metadata,
object_morphology, biofilm_structure}.csv` (all three currently header-only).

Per sample: voxel calibration in all three axes, the segmentation method **and** its
`segmentation_basis`, and the linkage columns.

**The stacks must be time-resolved.** The selected dynamic observable is
`biomass_autocorrelation_decay`, whose experimental quantity is *the correlation length of the
segmented biomass field over time*. Selecting it before the spatial calibration was deliberate:
it makes time-resolution a requirement of this campaign rather than of a second one.

> **No minimum timepoint count exists anywhere in this repository.** Two is the arithmetic
> floor for a decay rate and is not a design. The campaign must **declare** the count and
> interval and justify them. Do not infer a number from this repo — there is none to infer.

**A pitch is not identifiable from a cell volume.** `V_physical = V_sites · a³` fixes only the
product, so the morphology dataset is necessary and not sufficient; a second constraint must be
declared. And note what the gate already accepts: all seven species are declared
`biomass_parcel`. A CPM cell ID is a *computational biomass parcel*, not an organism — if any
species is re-declared with `literal_generation_claim_allowed = true` the gate returns
`MODEL_REVISION_REQUIRED` and outranks everything, because `compute_delta_H` has no shape or
connectivity term, `V_target` is one global scalar, and `divide_cell!` has no trigger.

---

## 3. Mass and volume — what the density is divided by

**Requirements:** `D-RHOWET`, `D-BLANK`, `D-DRY` → `data/calibration/materials/
{bulk_measurements, blanks}.csv`.

### 3.1 The volume must enclose what the balance weighed

A balance weighs cells **plus** matrix **plus** interstitial water. A `cells_only` mask
encloses none of the last two, so `ρ_wet = m/V` comes out systematically **high**.

| `volume_basis` | Accepted? |
|---|---|
| `whole_biofilm_envelope` | yes — segment the envelope |
| `cells_and_matrix` | only with a declared `pore_volume_fraction`; `V_envelope = V_stained/(1−p)` |
| `cells_only` | **refused** |
| `unresolved` | refused — what the volume encloses must be declared before it is divided into a mass |

### 3.2 Surface water is the measurement, not a nuisance

A wet mass without its handling record is not reproducible, so six columns are **schema-required
on every row**: `drain_orientation`, `drain_time_s`, `blot_material`, `blot_contact_time_s`,
`ambient_temperature_C`, `time_to_weighing_s`. Fix the procedure once and apply it identically
to samples and blanks.

Likewise `drying_protocol` and `drying_endpoint`: *"dried overnight"* and *"dried to constant
mass"* are different measurements and only one of them is a specification.

### 3.3 Blanks define the biofilm mass

    m_wet,biofilm = m_wet,sample+substrate − m_wet,blank_substrate

Blank subtraction is not a correction applied to a measurement; it **is** the definition. An
unblanked sample attributes the substrate and its retained water to the biofilm. Every bulk
measurement must name a `blank_sample_id`, and the blank is the same substrate carried through
identical handling and drying with no biofilm (`blank_kind` ∈ `dry_substrate`,
`hydrated_substrate`, `medium_exposed_abiotic`).

### 3.4 Detectability before replicates — the order is the point

Run `precision.detectability(biofilm_mass_g, blank_uncertainty_g)` **first**.

> A mass buried under the blank's noise is a **substrate or biomass problem, not a
> replicate-count problem.** Replicates shrink the standard error of a mean; they do not
> recover a signal lost to a subtraction. Grow more biomass per coupon, or use a substrate
> whose blank is lighter or more reproducible.

Only once it clears does `required_replicates()` apply, and its useful output is not the number
but `dominant_component` — it says where effort is worth spending:

| Dominant | What actually helps |
|---|---|
| `between_batch` | more independent cultures; more fields per coupon buy nothing |
| `within_batch_effective` | more imaging fields per coupon, which are cheaper than more cultures |
| `blank` | a lighter or more reproducible substrate, worth more than either replicate |

The floor is 2 batches — a variance estimate needs at least two. A blank uncertainty entered as
zero is reported as such: check it was measured and not assumed.

---

## 4. Composition — closed, not complete

**Requirement:** `D-COMP` → `data/calibration/materials/{elemental_analysis,
component_definitions}.csv`.

The transport loader needs elemental **mass fractions summing to 1**. The blend is
`w_i = Σ_k w_k · w_i|k` with `Σ_k w_k = 1` over components — water, dry organic (CHNS), ash,
salts, metals — each declared as a row with its own `mass_fraction`.

**Closed is the requirement; complete is not.** Anything unmeasured must be carried explicitly
as a residual component with a declared composition. It cannot be left out: `blend()` refuses
component fractions that do not sum to 1, precisely so a missing component is not absorbed as a
silent renormalization.

**Do not count a metal twice.** An element determined by combustion into the ash *and* assayed
separately is counted twice; `contributes_metals` and `counted_in_ash` exist to disambiguate,
and `blend()` refuses when two components claim the same element.

**Oxygen by difference is allowed and must be labelled.** `allow_oxygen_by_difference = true`,
but the fraction must be recorded as `analysis_method = by_difference` — otherwise a closure
that was *assumed* reads as a closure that was *measured*.

**The density must be measured, not blended.** `ideal_mixture_density()` returns evidence basis
`derived_proxy`, which is in `BLOCKED_EVIDENCE`, so a config built on it is refused. Use it as a
cross-check on §3, never as the value.

**Evidence basis decides admissibility.** Exportable: `direct_measurement`,
`manufacturer_datasheet`, `primary_literature`, `derived`. Refused: `proxy`, `derived_proxy`,
`synthetic`, `declared`, `unresolved` — *shipping a placeholder would launder it into a
measurement.* Every exported material also needs a non-empty `source_id` and
`system_conditions`; a composition without its growth medium, temperature and hydration state
cannot be judged applicable.

**Uncertainty is Monte Carlo and seeded** (`uncertainty.propagate`), because an interval that
changes between runs cannot be reviewed or regression-tested. It assumes independence, and the
wet and dry mass of one sample **share a balance calibration** — combine correlated inputs
before they reach it rather than declaring them independent for convenience.

---

## 5. Apparatus and source — metrology, not biology

**Requirements:** `D-GEOM`, `D-MEMBGEOM`, `D-MEMBMAT`, `D-MEDIUM`, `D-SRCPOS`, `D-ACTIVITY`,
`D-SPECTRUM`.

As-built dimensions of the biological domain and the membrane annulus; membrane and medium
density and composition (a manufacturer datasheet is exportable evidence); the source position
relative to the domain origin.

Two constraints the config enforces and the build must satisfy:

- `cylinder_radius_cm + membrane_thickness_cm` must stay inside the tallied lattice cube, or the
  membrane is partly untallied and its mass accounting is silently incomplete.
- The source position must sit **off every lattice plane** and inside the domain. The derived
  centre `n·pitch/2` lands on a plane for every even lattice size, so an explicit
  `position_cm` is required in practice.

**Activity is a dosimetry input only.** Transport reports eV per source particle and reads no
activity at all — `transport_state_hash` excludes it. For Gy/s:

    S = A(t) · Σ_j Y_j ,   A(t) = A₀ · 2^(−(t−t₀)/T_half)

Keep that separate from the PMF `p_j = Y_j / Σ_k Y_k`. The loader does **not** normalize, and
conflating the two silently mis-normalizes any nuclide emitting more than one photon per decay.

### What the model cannot represent, whatever is built

> **A real sealed capsule cannot be modelled.** Reference A1 is declared and *unsupported*:
> there is no active-core or capsule CSG, so no encapsulation and no self-attenuation.
> `source_spatial` accepts only `point_origin` or `line_z_axis`.

If the apparatus uses a sealed source, representing it as a point or line is a **declared
modelling limitation**, not a measurement, and `D-SRCSHAPE` requires it to be recorded as one.
Do not let a metrology report imply the capsule was modelled.

---

## 6. What is deliberately out of scope

**Metal reduction (`D-XRED`) is blocked by a units error, not by missing data.**
`compute_radial_biomass` returns dimensionless site occupancy and its unweighted mean is
assigned to `X_total` / `X_red`, which is neither a biomass fraction nor an active-reducer
fraction. Until that is replaced by a mass concentration the radiodialysis mapping stays
`BLOCKED`, and **measuring metal loading will not unblock it.** Taxonomic abundance is also not
functional activity: an `active_fraction` is required per row and may not be substituted.

Put it in the campaign only if it is wanted for its own sake. It is not on Reference D's
critical path.

**Growth and survival responses are unsupported by the current model**, not merely unmeasured:
the CPM has no birth and no death, so there is no simulated growth rate to fit against. No
quantity of data changes that without a model revision.

**The numerics are not measured.** `batches`, `particles`, `seed` and `mesh_coarsening_factor`
come from the **biofilm-specific resolution and history study**, which is a separate threshold
*after* the pitch exists. A0's sweep decided the phantom's numbers and its ledger rows carry
`sensitivity_domain = transport_numerical` so it cannot appear to have ranked biology. The
biofilm study must additionally report cell-, lineage- and generation-dose rank stability,
which A0 has no labels to measure.

---

## 7. Done

The campaign is complete when this prints a config instead of a refusal:

```sh
python calibration/scripts/reference_d_status.py --check
```

with the spatial gate at `READY_FOR_TIME_CALIBRATION`, the material OpenMC gate at
`READY_FOR_OPENMC_BIOFILM_TRANSPORT`, and `lattice_pitch_um` and `density_g_cm3` filled in
`data/calibration/voxel_binding.csv`.

Then `emit_transport_config(binding, spec)` under `evidence_policy = "measured_only"` and
`execution_class = "target_calibration"` yields Reference D's first transport config — and the
path that consumes it has already been demonstrated end to end on the synthetic system
(`docs/physical_reference_system.md` §2, reference S).

The radiodialysis gate will still read `BLOCKED`. That is the expected and correct result, and
it is not a failure of this campaign.

What happens after that is specified in `docs/feedback_gate_spec.md`: the target one-way dose
model, then the **offline** material-feedback gate, which asks whether a plausible change in
biological material state moves the radiation field by more than uncertainty and a declared
threshold. That gate needs no calibrated biological response, so it can be answered before the
response fitting is paid for — and if it fails, expensive two-way coupling is not justified.
