# Material basis and unit contract

**Date:** 2026-08-14 · **Enforced by:** `calibration/biofilm_calibration/materials/`

Four quantities in this project are routinely called "biomass" and are not interchangeable.
Conflating any two of them is a units error that no array shape will catch, because all four are
plain floats and three of them can plausibly read `0.05`.

| Quantity | Symbol | Unit | What it is |
|---|---|---|---|
| CPM site occupancy | — | dimensionless, [0,1] | fraction of lattice sites occupied |
| Bulk wet density | ρ_wet | g cm⁻³ | wet mass per hydrated biofilm volume |
| Dry biomass concentration | X_total | g cm⁻³ | dry mass per **hydrated** volume |
| Metal-reducing biomass concentration | X_red | g cm⁻³ | active-reducer dry mass per hydrated volume |

## The live defect this contract exists to stop

`compute_radial_biomass` (`biofilms_potts.jl:1309–1344`) counts occupied lattice sites per radial
bin and divides by the number of interior sites in that bin:

```julia
X_total = [counts[i] > 0 ? total_cells[i] / counts[i] : 0.0 for i in 1:Nr]
X_red   = [counts[i] > 0 ? red_cells[i]   / counts[i] : 0.0 for i in 1:Nr]
```

Those are **dimensionless occupancy fractions**. The unweighted mean of them is then assigned to
`RadiolysisParams.X_total` and `.X_red` at L1445–1446 and again at L1797 — fields declared
`g cm⁻³` at L1162–1163 — and consumed at L1241:

```julia
uptake = rp.k_ads * rp.X_total + rp.k_red * rp.X_red
```

with `k_ads`, `k_red` in `cm³ g⁻¹ s⁻¹` (L1156–1157). The product is `s⁻¹` **only if X is a mass
concentration**. Fed an occupancy fraction, the uptake rate is wrong by whatever ρ happens to be,
and nothing anywhere reports a problem.

Two further defects sit in the same three lines:

- **`X_red` is a species site fraction, not a biomass fraction.** It counts sites whose species is
  `SO` and divides by *all* interior sites (L1334–1336, L1342). So it is neither a fraction of
  biomass nor a fraction of active reducers.
- **The radial profile is collapsed with an unweighted `mean`.** In a cylinder an outer radial bin
  contains far more sites than an inner one, so an unweighted mean under-weights the outside. If a
  scalar is wanted it should be volume-weighted; the radial structure is computed and then thrown
  away either way.

## The rule

**Every quantity carries a declared basis, and conversions refuse untagged or wrong-basis input.**

```
site_occupancy      dimensionless    NEVER a concentration
wet_bulk            per wet mass or per hydrated volume
dry_biomass         per dry mass
ash                 per ash mass
medium              per medium mass
membrane_dry / membrane_hydrated
```

Tagging is the actual fix, not a range check. A biofilm dry-biomass concentration around
0.05 g cm⁻³ and a site occupancy of 0.05 are both entirely plausible numbers; no bound
distinguishes them. Only the label does.

`basis_conversion.Quantity` carries `(value, unit, basis)` and every converter checks the basis it
was handed. `site_occupancy` is accepted **nowhere** that a concentration is required, and the
error says so explicitly.

## Definitions

```
ρ_wet     = m_wet / V_hydrated                        [g cm⁻³]
X_total   = m_dry / V_hydrated                        [g cm⁻³]   dry mass, HYDRATED volume
w_water   = (m_wet − m_dry) / m_wet                   [—]
f_red,dry = m_dry,active_reducer / m_dry,total        [—]
X_red     = f_red,dry · X_total                       [g cm⁻³]
```

`X_total` is a dry mass over a **hydrated** volume. That mixed basis is deliberate — it is what
the kinetics need — and it is exactly the kind of thing that gets silently "simplified" to dry
mass over dry volume, which would be larger by roughly 1/(1−w_water).

## Two things that are not the same as each other

**Taxonomic abundance is not functional activity.** `f_red` is the *active reducer* dry-mass
fraction. The fraction of biomass that is *S. oneidensis* is an upper bound on it at best: cells
present but not respiring the metal do not reduce it. The schema keeps `taxonomic_fraction` and
`active_fraction` as separate columns, and the conversion from one to the other requires a
declared activity fraction or is marked blocked.

**Elemental composition is not isotopic composition.** The schema records elements. An OpenMC
material built by element uses natural isotopic abundance, which is right for biomass and wrong
the moment enriched or depleted material is involved.

## Oxygen by difference

CHNS analysers do not measure oxygen. It is routinely obtained as
`w_O = 1 − (w_C + w_H + w_N + w_S + w_ash)`, which is legitimate and must be **recorded as such**:
the `analysis_method` column distinguishes `measured` from `by_difference`, because by-difference
oxygen absorbs the error of every other element and of the ash determination.

## Mass closure

For components *k* with mass fractions *w_k* and elemental fractions *w_i|k*:

```
w_i_wet = Σ_k w_k · w_i|k          with Σ_k w_k = 1
```

The mixer rejects negative component fractions, sums outside tolerance, and **metal counted in
both the ash component and an explicit metal-loading component** — a double count that is easy to
create when ash is determined by combustion and the metal is also assayed directly.

## Density is measured or it is a proxy

An ideal-mixture density computed from component densities is labelled `derived_proxy`, never
`measured`. `evidence_basis = proxy` cannot be exported as a calibrated material.
