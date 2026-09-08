# The voxel binding: where the two calibration branches meet

**Date:** 2026-08-14 · **Status:** coherent, incomplete · **Emission:** REFUSED ·
**Contract:** `data/calibration/voxel_binding.csv` ·
**Checker:** `calibration/biofilm_calibration/integration.py`

The spatial branch says what an occupied CPM voxel *represents*. The material branch says what a
voxel is *made of*. Each passed its own gate without ever needing the other. This branch joins
them, and the join is where a specific class of error becomes visible.

## The contract

> **One occupied CPM voxel of physical volume `a³` represents part of a hydrated computational
> biomass parcel, and is assigned the hydrated-effective-biofilm material at a wet bulk density.**

Machine-readable as one row in `data/calibration/voxel_binding.csv`, with `lattice_pitch_um` and
`density_g_cm3` **blank** — blank means not measured, and is never read as zero.

## Why the check has to live here

Consider the binding the review named as incoherent:

> a literal cell fragment, whose density comes from bulk wet biofilm, whose lineage is claimed as
> a biological cell lineage.

Every clause passes a branch gate on its own. The spatial branch would accept `biological_cell`
as an entity kind — it is in the vocabulary. The material branch would accept a bulk wet
composition — that is what it is built to produce. **Neither branch can see the contradiction**,
because neither holds both halves. Only the join does:

- a cell interior is not bulk biofilm — bulk biofilm is cells *plus* their water *plus* their
  extracellular material, blended;
- a lineage of parcels is not a lineage of organisms.

`coherence_problems()` rejects that combination and five others:

| Rejected combination | Why |
|---|---|
| organism-scale entity + `hydrated_bulk_biofilm` material | a cell interior is not bulk biofilm; pick one scale |
| `biomass_parcel` + `cytoplasm` material | a parcel bundles cells with water and EPS, so its contents are not cytoplasm |
| `biological` lineage + `biomass_parcel` | parcel ancestry is not organism ancestry (`biofilms_potts.jl:1578-1583`) |
| `computational` lineage + organism-scale entity | if the entity is an organism, say whether its lineage is too |
| `dry_solid` + `hydrated_effective_medium` | the effective medium includes the water |
| `explicit_components` bound to a voxel | the mapper assigns one class per voxel; separate phases need sub-voxel homogenization that does not exist |
| `density_basis = site_occupancy` | not a density at all — the live model defect, caught here too |

## Coherent is not the same as ready

The declared binding **is** coherent. It still cannot produce a configuration, because a coherent
sentence with blanks in it is still a sentence with blanks. `evaluate()` refuses, and running it
against the repository's actual state prints:

```
REFUSING to emit a biofilm transport config
  binding coherent: yes
  + voxel binding is internally coherent
  - spatial gate is PROVISIONAL, not READY_FOR_TIME_CALIBRATION — without a
    selected lattice pitch there is no voxel volume, so no mass per voxel
  - material OpenMC gate is PROVISIONAL, not READY_FOR_OPENMC_BIOFILM_TRANSPORT
    — no measured density or closed composition exists to assign
  - lattice_pitch_um is unset
  - density_g_cm3 is unset
```

Three independent conditions, all required: the binding must be coherent, **and** both gates must
be ready, **and** the two numbers must exist. Coherence alone was never sufficient — that is
precisely the failure mode where a well-formed but unmeasured model gets shipped.

### The binding is not the apparatus

Emitting a config needs a second input, and it lives in
`biofilm_calibration.reference_system`. The binding's subject is one voxel; the cylinder, the
membrane thickness, the boundaries, the source and the medium are properties of the *system*, and a
loadable `biofilm_cylinder` config requires all of them plus three material classes. Treating the
binding as if it held the whole apparatus is why the emitter could previously produce only a pitch
and one density — a fragment, not a configuration.

`emit_transport_config(binding, spec)` now renders a complete config or refuses naming every
missing field, and its checks come in three groups that are deliberately not interchangeable:

| group | question | bypassable? |
|---|---|---|
| `coherence_problems` | does the binding contradict itself? | never |
| `structural_problems` | is this physically constructible — closure, congruence, source placement, units? | never |
| `evidence_problems` | is anyone entitled to believe these numbers? | only by `evidence_policy = "synthetic"` |

A synthetic system may invent its values. It may not invent a geometry that cannot be built, a
composition that does not close, or a source born on a lattice plane — otherwise "synthetic" becomes
a universal escape hatch around the contracts this document exists to describe.

## What each blank waits on

| Blank | Waits on | Gate |
|---|---|---|
| `lattice_pitch_um` | morphology dataset, declared acceptance thresholds, declared domain semantics | spatial, `PROVISIONAL` |
| `density_g_cm3` | bulk wet/dry measurements with a stated surface-water protocol, closed elemental composition | material OpenMC, `PROVISIONAL` |

The pitch is additionally **not identifiable** from a cell volume alone (`V_physical = V_sites·a³`
constrains only the product), so a morphology dataset is necessary but not sufficient — a second
constraint has to be declared.

## Downstream consequences now recorded in the ledger

- **`cpm.seconds_per_mcs_calibration` is downstream of the pitch.** `Δt_MCS = a²·S_sim / S_exp`
  needs `a`, so the time calibration cannot begin until the spatial gate clears. The ledger row
  now says so.
- **`response.growth_survival` is currently unfittable in principle**, not merely unmeasured: the
  CPM has no growth dynamics and `divide_cell!` has no trigger, so there is no simulated growth
  rate to fit a survival response against.
- **`response.metal_reduction` is blocked by a units error upstream**, not only by missing data:
  `X_red` is a site fraction of one species divided by *all* interior sites, which is neither a
  biomass fraction nor an active-reducer fraction.
- Every response is a **per-parcel** response, not a per-organism one, because that is what a cell
  ID is.

## Which radiation phenotype this model actually represents

Five phenomena travel under one loose word, and they are not interchangeable. A survival curve
cannot calibrate a directional term, and a transcriptome cannot set a material density.

| phenotype | what it means | represented here? |
|-----------|---------------|-------------------|
| **radiotropic** | directional growth or redistribution toward or away from a source | **yes** — this is `response.hamiltonian` |
| radiotrophic | ionizing radiation measurably **enhances growth or metabolism** | no — requires growth |
| radioresistant | elevated **survival** after exposure | no — requires death |
| melanized_radioprotective | melanin measurably alters damage, survival or shielding | no — requires damage |
| radiation_responsive | molecular, morphological or material change after exposure | partly — `response.melanin` |

`compute_delta_H` has four terms: adhesion, volume, radiation, melanin. The radiation term adds
`β_ion[s]·I_γ` when a parcel gains a site, so the sign of `β_ion` decides whether that parcel
drifts up or down the dose gradient. **That is a tropism.** The Hamiltonian is a Metropolis
acceptance functional in arbitrary units, not an energy budget, so a negative `β_ion` is a spatial
preference and not metabolic gain — the model creates no biomass for such gain to act on.

Radiotrophy and radioresistance are therefore not weakly supported here; they are **structurally
absent**, for the same reason `response.growth_survival` is `unsupported_by_current_model`. A
model with no birth and no death cannot express either, and no quantity of data changes that
without a model revision.

This matters beyond bookkeeping. The published D10 evidence — which is strong, and is what anchors
the per-species `β_ion` magnitudes — is *radioresistance* evidence being spent on a *tropism*
coefficient. That substitution is legitimate only as a declared modelling choice, and it is
declared here rather than left implicit in a sign convention.

Radiotrophy additionally fails on evidence, not only on representability: it is not established
for any of the seven species. See `docs/research/radiotrophic_compatibility_audit.md`, verdict
`TARGETED_LAB_EXPERIMENT_REQUIRED`.

### …and the tropism is not driven by the term that is tabulated

Naming the radiation term correctly is not the end of it. At the shipped `I0 = 1.0` and
`T_cpm = 5.0`:

| term | ΔH | Metropolis acceptance bias |
|------|-----|---------------------------|
| `ΔH_rad`, radiotropic species (β = −5e-5) | −5.0e-5 | 1.000010 |
| `ΔH_rad`, most radiosensitive species (β = 7.5e-2) | +7.5e-2 | 0.985 |
| `ΔH_mel` at the reported M = 1.44 | −0.720 | **1.155** |

`β_ion` — the one parameter Table 2 tabulates per species, and the one the sign convention is
written around — biases acceptance by **one part in 10⁵** *for one role of the two negatively
signed species occupying a site*, which is not the term's reach: signed by role it reaches
7.505e-2. The melanin term biases acceptance by 15.5%, larger by about an order of
magnitude more, through a coefficient of `0.5` hard-coded at its call site and appearing in no
table and no configuration file.

So the model's tropism is **melanin-mediated**. Radiation still drives it, but indirectly:
`melanin_drive` is copied from the radiation field, giving

```
radiation → production (α_M, tabulated) → M → ΔH_mel (0.5, untabulated) → tropism
```

Two consequences. A sensitivity analysis over `β_ion` would report the radiation response as
negligible and be right about the coefficient while wrong about the mechanism. And a parameter
that decides the headline result must not live only at its call site: it is now ledgered as
`cpm.melanin_coupling` in `data/parameter_provenance.csv`, `status=blocked`,
`sensitivity_rank=high`. It is one of five rows in that ledger not ranked `unknown`, and the only
one outside the A0 transport sweep: `transport.mesh.base_dimension` and
`transport.mesh.coarsening_factor` are also `high`, `transport.batches` and `transport.particles`
`medium`, and all four carry `sensitivity_domain = transport_numerical`, which ranks numerical
convergence rather than a biological response. Every biological-response row is still `unknown`.
This one was measured here rather than awaiting a sweep.

It also cannot be fitted alone. Only the product `α_M · k` reaches the dynamics, so
`cpm.melanin_coupling` and `response.melanin` are jointly identifiable at best.

### And the coefficient vector is identifiable only up to a common scale

That is the second of three degeneracies, not the last. The Metropolis rule
(`biofilms_potts.jl:806`) reads only `ΔH/T`, and each active term is linear in exactly one
coefficient — `ΔH_adh` in `J`, `ΔH_vol` in `λ_V`, `ΔH_rad` in `β_ion`. So scaling every
coefficient *and* `T_cpm` by a common `c > 0` leaves every acceptance decision unchanged,
including which branch of the `ΔH <= 0` ternary is taken, so the generator is consulted at the
same proposals and the random stream is untouched. **The shipped `T_cpm = 5.0` is one redundant
degree of freedom**, and no fit against a trajectory can recover the scale of the coefficient
vector: a fit must fix one coefficient by independent measurement, or report the ratio and say
so. Ledgered as `SCALE-01`.

With `c` a power of two this is bit-exact rather than statistical — scaling by a power of two is
exact in IEEE-754 and commutes with rounding through the sums and the division, so `(2ΔH)/(2T)`
is exactly `ΔH/T`. `tests/hamiltonian_scale_invariance.jl` asserts lattice, registry and volume
equality after 30 MCS, having first asserted `f(2θ) == 2·f(θ)` per term so that the premise is
tested before the conclusion.

Two consequences worth stating. `biofilms_potts.jl:81` calls `T_cpm` a *"CPM Boltzmann
temperature"*, which is wrong rather than incomplete: a quantity fixed only up to a shared scale
has no absolute units to be measured in, and the repository denies exactly this reading for `H`
in four other places (`SCALE-02`). And the melanin term is linear in a hard-coded `0.5`
(`:593`, `:596`) rather than in a parameter, so the scaling cannot reach it — with a non-zero
melanin field the symmetry provably breaks. That is asserted as a control, which makes the test
a standing detector for **any** future coefficient in `H` placed outside the parameter
interface (`SCALE-03`, rediscovering `PP-T2-29` from the opposite direction).

The third degeneracy is spatial: `V_physical = V_sites · a³` fixes only a product, so the
lattice pitch is not identifiable from a cell volume — `docs/calibration/cpm_spatial_calibration.md`.

#### The reporting convention, and why the shipped constants are not changed

Because only ratios reach the dynamics, **every coefficient is quoted in units of `T`**.
`T_cpm = 5.0` is a choice of scale, not a measurement; only `J/T`, `λ_V/T` and `β_ion/T` have
content, and the melanin statement above is cleanest read the same way — the product `α_M · k`
in units of `T`. Any future fit should work in that gauge, pinning `T = 1` as a *fitting
coordinate* and optimising `ln(J/T)`, `ln(λ_V/T)`, `ln(β_ion/T)`, so that estimator variance is
never injected into a direction the data cannot constrain.

**That gauge is a coordinate choice, not an edit to `biofilms_potts.jl`, and the distinction is
numerical rather than pedantic.** Setting the shipped `T_cpm` to 1.0 and dividing the
coefficients by 5 is free algebraically and is *not* free in floating point:

| rescaling | `ΔH/T` unchanged | trajectory after 30 MCS, N=12 |
|---|---|---|
| × 2 (a power of two) | **7659 / 7659** | bit-identical |
| × 1/5 (to reach `T = 1`) | **1957 / 7659 — 74.45% differ** | **365 sites differ** |

Every shipped constant survives the round trip `(x/5)·5 == x`, and no proposal changes the sign
of `ΔH`, so the `ΔH <= 0` branch is untouched. The divergence is entirely in the last bits of
the comparison: `ΔH` is a **sum**, and `Σ fl(x_i/5) ≠ fl(Σ x_i)/5`. `1/5` is not a power of two,
so there is no exact route from `T_cpm = 5` to `T_cpm = 1`.

Changing the constants would therefore break the source→snapshot chain for
`manuscript42_seed42_8cbb8ee`, whose manifest pins `cpm_parameters.T_cpm = 5`, `λ_V = 10` and
`β_ion` alongside 122 hashed artifacts, and with it Figure 5, every `label_state_hash` and the
derived inert-signal run. The convention is documentary; the constants stay. Ledgered as
`SCALE-04`.

## What this branch merged

Both branch-local source registries are folded into `data/sources.csv` (only `SPATIAL_DECL_2026`
had content; the material registry ships empty because nothing has been cited yet), plus
`INTEGRATION_DECL_2026` for this declaration. The provenance ledger gains the binding as a row
and has its calibration rows repointed at the harnesses that now own them.

Neither calibration branch touched these shared files, which is why both merged without conflict.
