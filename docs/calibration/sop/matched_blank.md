# SOP — the matched blank

Covers `D-BLANK`. Lands in `data/calibration/materials/blanks.csv`, and every row of
`bulk_measurements.csv` names one through `blank_sample_id`.

**This SOP runs before biomass production, not after it.** §3.4 of
`reference_d_measurement_protocol.md` — "Detectability before replicates — the order is the point" —
is the argument: a biofilm mass under the blank's noise is a substrate problem, and no replicate
count fixes a substrate problem. A discipline-ordered index would file blanks under analytical
chemistry and run them fourth.

## Preconditions

Balance calibration log current. Water quality recorded as resistivity **and TOC** — demineralisation
removes ions and passes neutral organics, so resistivity alone does not characterise the water a
low-biomass measurement sits in.

## Procedure

**1. One blank per bulk measurement.** Not one per batch. `D-RHOWET`'s acceptance criterion requires
every row to name a `blank_sample_id`, so a shared blank fails the criterion whatever its quality.

**2. The blank is the identical cycle with no biomass.** Same substrate lot, same wetting, same
`drain_orientation`, `drain_time_s`, `blot_material`, `blot_contact_time_s`, `ambient_temperature_C`
and `time_to_weighing_s` as `harvest_dewatering.md` specifies. A blank that differs in any of those
measures a different thing from the sample it is subtracted from.

**3. Abiotic control, separately.** A no-biomass vessel through the same handling, to separate
plasticware wall adsorption from biofilm uptake. This is not the same as the substrate blank and
does not replace it.

**4. Report the floor, not just the value.** The output is a mean **and a standard deviation across
blanks**. The standard deviation is the detectability floor, and it is what the preconditions of
`harvest_dewatering.md` require before biomass is grown.

## The disposition, stated before the run

If the expected biofilm mass is **within the blank's standard deviation**, the measurement does not
work on this substrate at this scale. The response is a different substrate or a larger coupon, not
more replicates. Writing that here means it is decided before a number exists to argue with.
