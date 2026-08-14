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
sentence with blanks in it is still a sentence with blanks. `emit_transport_config()` refuses, and
running it against the repository's actual state prints:

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

## What this branch merged

Both branch-local source registries are folded into `data/sources.csv` (only `SPATIAL_DECL_2026`
had content; the material registry ships empty because nothing has been cited yet), plus
`INTEGRATION_DECL_2026` for this declaration. The provenance ledger gains the binding as a row
and has its calibration rows repointed at the harnesses that now own them.

Neither calibration branch touched these shared files, which is why both merged without conflict.
