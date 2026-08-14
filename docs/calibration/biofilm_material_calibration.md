# Biofilm material calibration

**Date:** 2026-08-14 · **Gates:** `OPENMC: PROVISIONAL` · `RADIODIALYSIS: BLOCKED` ·
**Unit contract:** `docs/calibration/material_basis_contract.md` ·
**Data:** `data/calibration/materials/` · **Harness:** `calibration/biofilm_calibration/materials/`

This branch defines what material occupies an occupied voxel and in what units, and builds the
conversions and closure checks that would let a measured composition be exported. It exports
nothing, because nothing has been measured — and it is designed so that when something is, the
easy mistakes are refused rather than shipped.

## 1. Material model kind

`material_model_kind = hydrated_effective_medium`: **one** OpenMC material representing the wet
biofilm bulk — water, dry organic biomass, ash, and any sorbed or precipitated metal blended by
mass.

This is not a preference. `coupling/biofilm_openmc/materials.py` maps one material class per
occupied voxel; representing dry solid and interstitial water as separate phases would need
sub-voxel homogenization or additional classes, and neither exists. `explicit_components` is a
declared alternative that the export gate refuses today, with that reason.

It also follows from the spatial branch: a cell ID is a computational biomass *parcel*, not a
cell, so "occupied voxel" cannot mean "cytoplasm". A parcel bundles cells and their extracellular
material, and a bulk hydrated composition is what matches it.

**Do not create seven materials because there are seven species.** Per-species compositions need
evidence that their transport-relevant density or elemental composition actually differs;
otherwise they are seven copies of one number wearing different labels.

## 2. Four quantities that are not each other

Set out in full in the basis contract. In brief: site occupancy (dimensionless), ρ_wet, X_total
and X_red (all g cm⁻³). Three of them can plausibly read `0.05`.

The defect this branch exists to stop is live. `compute_radial_biomass`
(`biofilms_potts.jl:1341–1342`) returns occupancy fractions; their **unweighted** mean is assigned
to `RadiolysisParams.X_total` and `.X_red` — declared `g cm⁻³` — at L1445–1446 and L1797, then
multiplied by rates in `cm³ g⁻¹ s⁻¹` at L1241. The product is `s⁻¹` only if X is a mass
concentration.

Two further problems sit in the same lines: `X_red` counts sites of species `SO` divided by *all*
interior sites, so it is neither a biomass fraction nor an active-reducer fraction; and an
unweighted mean over radial bins under-weights the outside of a cylinder, where most of the sites
are.

**The fix is basis tagging, not bounds.** `Quantity` carries `(value, unit, basis)` and every
converter checks what it was handed. `site_occupancy` is accepted nowhere a concentration is
required, and the error names the offending lines.

## 3. Derived quantities

```
ρ_wet     = m_wet / V_hydrated
X_total   = m_dry / V_hydrated          dry mass, HYDRATED volume
w_water   = (m_wet − m_dry) / m_wet
f_red,dry = m_dry,active_reducer / m_dry,total
X_red     = f_red,dry · X_total
```

with the round-trip `X_total = ρ_wet · (1 − w_water)` pinned by test. `X_total` mixes bases on
purpose — dry mass over hydrated volume is what the kinetics need — and that is exactly the thing
that gets "simplified" to dry mass over dry volume, which is larger by roughly 1/(1 − w_water).

`active_from_taxonomic` refuses without a declared activity fraction, rather than assuming every
cell of a metal-reducing species is reducing metal.

## 4. Measurement schemas, and the columns that carry the weight

Six tables under `data/calibration/materials/`, all shipping header-only. The columns that are
easy to omit and expensive to lack:

- **`surface_water_removal_protocol`** — a wet mass without a stated drainage or blotting
  procedure is not reproducible, and the surface water it includes goes straight into ρ_wet and
  w_water.
- **`drying_endpoint`** — "dried overnight" and "dried to constant mass" are different
  measurements.
- **`basis`** on every elemental row — so dry-basis CHNS values and wet-basis water fractions
  cannot be summed without an explicit conversion.
- **`analysis_method`** distinguishing `measured` from `by_difference` — oxygen is normally
  obtained by difference, which is legitimate and absorbs the error of every other element and of
  the ash determination.
- **`taxonomic_fraction` and `active_fraction` as separate columns** — see above.
- **`counted_in_ash`** — metal determined by combustion into ash *and* assayed separately is
  counted twice. `blend()` rejects it.

Bounded scenarios (`low_solids`, `central`, `high_solids`, and metal-loading variants) are
expressed as separate `scenario_id`s in `component_definitions`, each retaining its own sources,
conditions, uncertainty and evidence basis — not collapsed into one false central value.

## 5. Mixing and closure

`w_i_wet = Σ_k w_k · w_i|k`, with rejections for negative fractions, component fractions not
summing to 1 (a forgotten component would otherwise be absorbed as a silent renormalization),
elemental fractions not closing within a component, and metal double-counted between ash and an
explicit loading.

Density is measured where possible. `ideal_mixture_density` returns its value tagged
`derived_proxy` and the export gate refuses to ship a proxy as calibrated: real biofilm has
excluded volume and pore structure that volume additivity does not model.

Uncertainty propagation is Monte Carlo and **seeded** — an interval that changes between runs
cannot be reviewed or cited. Independence between inputs is an assumption and often wrong (wet and
dry mass of one sample share a balance calibration), so correlated inputs must be combined before
they reach the propagator.

## 6. Medium and membrane

**Medium.** The exact growth-medium recipe, once selected. A0's pure water cannot silently stand
in for a salt-bearing biological medium.

**Membrane.** Fixed for this phase: geometry, density and composition constant, no dose response.
A stoichiometric PDMS composition remains a **proxy** for an unknown commercial silicone — fillers
and additives are unresolved without the exact grade — and the Nafion irradiation envelope stays
scoped to Nafion, as recorded in `docs/physical_reference_system.md` §3.

## 7. Gates

Two, because they need different things and one verdict would hide which half was missing.

**`OPENMC: PROVISIONAL`** — needs a hydrated bulk density or declared distribution, a closed
elemental composition, the water/solids basis resolved, the exact medium or a declared proxy,
membrane status, uncertainty, and a mapping compatible with the spatial branch. All missing;
none blocked in principle.

**`RADIODIALYSIS: BLOCKED`** — needs everything above plus physical `X_total` and `X_red`, a
definition of the active-reducer fraction, and a spatial mapping rule. This one is **blocked by a
units error in the model, not by missing data**: until `compute_radial_biomass` supplies a mass
concentration rather than a site fraction, no measurement makes the mapping correct. Declaring a
`spatial_mapping_rule` in the acceptance template does not fix it either.

## 8. Stop conditions

Do not export while any holds: wet and dry bases unspecified; hydrated volume unknown; wet density
inferred from dry density without a declared model; elemental composition incomplete;
oxygen-by-difference undocumented; medium recipe unknown; a PDMS proxy labelled as exact silicone;
metal loading inferred from a simulation parameter rather than an assay; species occupancy used as
a dry-mass fraction; or literature values from incompatible growth conditions averaged into a
false central value.
