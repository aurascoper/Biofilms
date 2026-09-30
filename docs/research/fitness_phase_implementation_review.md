# Revised implementation contract

The brief is implementable after the changes below. These are model and acceptance decisions, not results. Base revision: `8ad9439afc32e5a0db606ddf99ea2827a054d1c4`. The supplied prompt digest and local tools are in `diagnostics/fitness_phase/preflight.json`.

## Red-team findings

1. **There is no frozen V1 production protocol yet.** The current specifications explicitly say so. Retain their seeds 42–297, N=40, six parcels per species, and MCS 100 endpoint; freeze synthetic parameters and analysis before production. Do not call this replication of a completed study or external preregistration.
2. **Resorption requires a per-parcel target.** The serial and JACC kernels currently use a scalar. Add an explicit serial energy callback and per-parcel policy outside the legacy state. Keep the default path bit-identical. Refuse accelerated execution of this extension; existing JACC parity covers only the legacy model.
3. **The existing nutrient field already consumes nutrient.** Growth-OFF must retain that exact update. Growth-ON uses a separate conservative reaction/diffusion closure with a no-flux mask, named synthetic concentration, biomass and length units, and a transfer budget. The legacy wall reservoir and new closed domain must not be conflated.
4. **Mortality changes geometry as well as membership.** The existing phase specification uses a moving Boolean 26-neighbour contact graph. Keep it for V1. V3 may explicitly choose a face-area graph; the geometry path and row normalisation must be reported. The present implementation keeps the Boolean graph across both stages to isolate growth, with face-area graph supported as a separately declared option.
5. **Physical time and model time are different contracts.** Use MCS for growth, phase and display. The declared 60 s/MCS and 0.01 Gy/s are synthetic exposure settings only; the lattice spacing is one synthetic length unit, not micrometres. No biological calibration is claimed.
6. **A target radius is not biomass.** Initialise synthetic biomass to the global target in synthetic mass units, density one per synthetic volume. Target voxel count is B/(rho*voxel_volume); actual volume may lag. Death retains biomass in resorbing material until physical removal, and its loss is a named budget sink. Division conserves biomass and archives ancestral dose.
7. **Existing division can refuse fragmented daughters.** Use its deterministic principal-axis split and 26-connectivity refusal. A deferred division is recorded without mutation. CPM itself does not guarantee connected parcels; do not claim it does. The enabled public division API rejects dying parents.
8. **Membership can manufacture coherence.** Preserve frozen-phase removal accounting, fixed-window survivor IDs, unwrapped pair paths and zero-coupling replay. For V3, separately record phase correlation introduced by daughter inheritance. A replay retains every parcel until its recorded death; measurement cohorts do not alter the dynamics.
9. **The finite all-to-all reference is not a universal graph threshold.** A normal-frequency continuum threshold cannot certify a finite, moving, degree-normalised graph. Verify two-node analytic locking, synchronous updates, graph normalisation and time refinement; report moving-graph locking only under frozen drift and frequency criteria. A null result is permitted.
10. **Raw checkpoint hashes compare irrelevant metadata.** Compare canonical nonphase state at every MCS, including both process RNGs, histories and targets. Compare complete state for same-arm restart. Use a versioned extension checkpoint with explicit runtime/source/config identity, separate from the legacy exchange schema.
11. **Growth-OFF has no defined biomass dynamics.** Export biomass as unavailable with a validity mask when that extension is off. Export the declared zero growth contribution and mode flag. Conditional survival is a model probability conditional on each recorded exposure history, not the realised alive state.
12. **Visualization acceptance is independent of a file write.** Version the new CellData schema; independently read arrays and masks with VTK, reopen ParaView from a relocated directory, and compare times and source identity. Rendering and analysis cannot consume simulation RNGs. Missing views or receipts must remain incomplete.

## Implementation boundary

Model IDs are `cpm-sp-v1` (survival and passive phase), export schema `cpm-sp-vti-1` (the same trajectories with export) and `cpm-sp-growth-v1` (synthetic growth/division). No phase feedback, periodic topology, phase noise, forced entrainment, ART, Dsup factor, repair or dormancy is enabled. Unsupported modes refuse input.

V1 keeps the specified mechanical → fields → occupied-site dose → conditional Bernoulli hazard → cleanup → phase order. New hazard draws use a separate MersenneTwister, in increasing parcel-ID order, skipping zero hazard. A dying parcel immediately leaves the phase graph and ceases production and uptake on later steps. It retains contact energies and radiation/melanin acceptance until resorbed below two voxels. The legacy zero-volume cleanup remains unchanged.

V3 inserts conservative nutrient/biomass updates at the field stage and deterministic eligible divisions after hazard and cleanup, before phase. It retains the predecessor mechanical-first split. Daughters inherit the parent's phase and frequency exactly; independent future conditional Bernoulli draws express the declared memoryless hazard. Targets follow biomass; biomass is independent of accepted voxel copies. Minimum age, division mass, daughter size, substeps and rates are configuration fields frozen before outcomes.

The serial reference loads the existing model and reuses its proposal schedule, mechanics, field updates and division geometry. New state lives in an extension owner; legacy checkpoint and CellInfo layouts are preserved. The extension cannot be passed to JACC.

## Acceptance and reporting

The required paired identities compare baseline/phase-only and survival-only/combined after each MCS and restart, for growth OFF and ON. They include CPM arrays, dose, parcel histories, biomass, target volumes, event order and nonphase RNG states. Export ON/OFF must preserve complete simulation state. The legacy disabled path has its own golden comparison.

Run controls before a committed production PREDICTION. Record tests separately from outcomes. Freeze seed ranges, horizons, synthetic parameter vectors, phase substeps, analysis windows, thresholds, output schedules and resource limits. Commit source and configuration before production; retain all trial cases. Original manuscript results remain tied to their original configuration. New results describe only the separately named synthetic models.

VTI uses explicit CellData, Int64 IDs, x-fastest ordering, one extra point per voxel axis, masks for undefined quantities, MCS PVD times and a hashed sidecar manifest. A versioned schema and independent reader define the scientific export. The legacy exporter remains available.

Official API references: [VTK XML](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html), [WriteVTK associations](https://juliavtk.github.io/WriteVTK.jl/stable/grids/datasets/), and [ParaView state](https://docs.paraview.org/en/latest/UsersGuide/savingResults.html). Installed versions, not the latest online labels, govern acceptance.
