# SOP — coupon harvest, surface-water removal, wet mass and dry mass

Covers `D-RHOWET` (wet mass and hydrated volume) and `D-DRY` (dry mass to constant weight).
Lands in `data/calibration/materials/bulk_measurements.csv`.

**Written from the schema, not beside it.** Every step below fills a named column of
`bulk_measurements.csv`. If a step here has no column, the step is wrong or the schema is
incomplete — resolve that before running, never afterwards.

## Preconditions — none of this runs before all four hold

1. `D-APPROVAL` — institutional biosafety approval, strain- and protocol-specific.
2. Strain identity verified against an accession record. **For *O. intermedium* AM7 no
   culture-collection accession was located in any repository** — `radiotrophic_lab_gap.md:87`,
   which also records that the strain may simply be unobtainable and that CDC's 2022-12-19 update
   directs *Brucella*-identified isolates to a Class II BSC with public-health referral. For that
   strain this step may **terminate** the SOP rather than delay it.
3. Balance calibration log current, and pipette calibration where volumes are dispensed.
4. **`D-BLANK` has been run and the detectability floor is known.** See
   `matched_blank.md`. A biofilm mass under the blank's noise is a substrate problem no
   replicate count fixes, and §3.4 of `reference_d_measurement_protocol.md` puts detectability
   before replicates for that reason.

## Procedure

**1. Harvest.** Record `coupon_id`, `culture_batch_id`, `growth_condition_id`, `medium_batch_id`,
`sample_id`, `replicate_id`, and the `paired_sample_group` that links this coupon to its imaging
counterpart. Record `measurement_order` — the position in the session, which is what lets a drift
in the balance be separated from a difference between samples.

**2. Surface-water removal — the measurement, not a nuisance.** §3.2 of the protocol is explicit
that this is part of the measurement. Fix and record all six:

| Column | What to record |
|---|---|
| `drain_orientation` | the fixed angle the coupon is held at |
| `drain_time_s` | held constant across the batch |
| `blot_material` | one material for the whole study |
| `blot_contact_time_s` | held constant |
| `ambient_temperature_C` | measured, not assumed |
| `time_to_weighing_s` | from end of blot to balance reading |

**3. Wet mass.** `wet_mass_sample_plus_substrate_g` with `wet_mass_uncertainty_g`. The substrate is
weighed with the sample; the blank supplies the subtraction.

**4. Hydrated volume.** `hydrated_volume_cm3` with `volume_uncertainty_cm3`, `volume_method`, and
`volume_support`. **`volume_basis` must be `whole_biofilm_envelope`, or `pore_volume_fraction` must
be declared** — this is `D-RHOWET`'s acceptance criterion, not a preference. See
`imaging_segmentation.md`, where the same constraint decides the mask.

**5. Area scaling, if the imaged and weighed areas differ.** `imaged_area_cm2`, `weighed_area_cm2`,
`scaling_method`, `scaling_uncertainty`. The uncertainty on the bridge is declared, never assumed
to be zero.

**6. Dry mass to constant weight** — `D-DRY`. Record `drying_protocol` and `drying_endpoint`:
successive weighings agreeing within a stated tolerance, at a stated temperature.
`dry_mass_sample_plus_substrate_g` with `dry_mass_uncertainty_g`. Ash, if taken, uses
`ash_protocol` and `ash_mass_g`.

**7. Close the row.** `blank_sample_id` naming this sample's matched blank, `quality_flag`, and
`source_id`.

## A protocol-variance measurement this SOP is also for

Run the wet/dry cycle on **n replicates of one biomass** and report the CV. **A large CV retires wet
mass as a normaliser** — an output that invalidates the procedure it characterises, which is why it
is stated as an expected outcome rather than a risk. If more than one operator will run this,
repeat across operators and report between-operator CV separately.

## What this SOP does not establish

It produces a density. It does not produce a growth rate, a survival response, or a dose response,
because the CPM has neither birth nor death — no quantity of data changes that without a model
revision (protocol §6).
