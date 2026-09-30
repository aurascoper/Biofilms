# Synthetic survival, passive phase and growth

This standalone serial runner exercises defined synthetic mechanisms. It does not solve the manuscript's original PSDE, fit biological survival, implement transport feedback, or change the legacy movement experiment.

`docs/research/fitness_phase_implementation_review.md` records the red-team findings and revised contract. `PREDICTION.md`, `protocol.toml` and `config.toml` fix the local experiment. Read the execution report for completed checks and results; source availability alone is not execution evidence.

## Equations and state

| Quantity | Equation or rule | Code and state |
| --- | --- | --- |
| Incremental exposure | occupied-voxel mean dose rate times explicit exposure clock | `model.jl`, `Parcel.increment`; cumulative `dose` remains separate |
| Hazard | ΔH = alpha_s ΔD; commitment probability = -expm1(-ΔH) | `hazard_draw`, `Parcel.hazard`; independent mortality RNG |
| Viability | irreversible alive → dying at commitment; dying exits phase, production and uptake immediately | `step!`; target zero; no further hazard trials |
| Physical removal | mechanical shrinkage followed by cleanup below two voxels | core energy callback and `step!`; dose continues until removal |
| Conditional survival | exp(-cumulative hazard), frozen at commitment | `export.jl`; not the realized alive flag |
| Phase | dθ_i/dt = ω_i + K Σ_j A_ij sin(θ_j−θ_i)/Σ_j A_ij | `phase.jl`; synchronous explicit Euler, isolated coupling zero |
| Graph | Boolean contact through the existing 26-neighbour adhesion stencil | `contact_graph`; rebuilt after geometry and lifecycle changes |
| Uptake | u_i = B_i q_s mean[C/(K_s+C)] over occupied voxels | `growth.jl`; one implicitly computed nutrient transfer supplies biomass |
| Biomass | ΔB_i = Y_s ΔU_i − m_s B_i Δt | growth substeps; source/sink ledger includes removed corpses |
| Target | V*_i = B_i/(rho_s voxel_volume) | side table `Parcel.target`; enabled serial volume energy only |
| Division | minimum age, threshold biomass and two connected daughters | existing deterministic divider; conservative volume-proportional biomass split |
| Net proxy | conditional growth rate minus instantaneous mortality hazard | `export.jl`; not realized biomass change or population fitness |

Within one MCS: validate → mechanics → fields/growth → incremental dose → commitment → physical cleanup → division when enabled → current graph → phase. Enabled mechanics records removals in ascending parcel-ID order so a restored dictionary cannot reorder events. No phase value enters mechanics, nutrient, melanin, survival or random streams. Daughter IDs are new; their phases and frequencies inherit the parent's values. Dose since birth resets and ancestral dose is retained separately. Equal daughter phases are an imposed inheritance rule.

Growth uses closed no-flux diffusion, positive implicit Monod uptake and explicit maintenance under its positivity bound. Dying parcels keep biomass until removed, when it enters the sink ledger. The local growth and net rates are per MCS; physical absorbed-dose bookkeeping uses the explicitly synthetic clock. Spacing is synthetic, with no inferred micrometre calibration. Growth OFF preserves the existing nutrient update and has no defined biomass state.

## Parameter registry

The exact species order and all inputs are in `config.toml`. The seven alpha values `[.01,.002,.015,.04,.02,.06,.03]` Gy⁻¹, all growth parameters, phase frequencies, coupling, density, exposure rate, clock and spacing are **assumed synthetic inputs**. None is a measured seven-species calibration. Existing adhesion and direct radiation constants come from unchanged `CPMParams`; source SHA records their version. No D10 conversion is required by this runner. The existing `ResponseConfig` validator and unset-clock refusal remain in force. JACC, forcing, phase noise and signed coupling modes refuse input.

## Commands

Run from the repository root. `OUTPUT` is a new external directory. Production refuses dirty source, a changed configuration, seeds/horizons outside the frozen plan, or a changed export schedule.

```sh
julia --project=. diagnostics/fitness_phase/tests.jl
julia --project=. diagnostics/fitness_phase/run.jl diagnostics/fitness_phase/config.toml OUTPUT 42 297 100 0
julia --project=. diagnostics/fitness_phase/run.jl diagnostics/fitness_phase/config.toml OUTPUT 42 57 200 1
julia --project=. diagnostics/fitness_phase/controls.jl diagnostics/fitness_phase/config.toml OUTPUT/controls
julia --project=. diagnostics/fitness_phase/resume.jl OUTPUT/growth_0_seed_42/combined 50 100 OUTPUT/resume.json
python diagnostics/fitness_phase/validate_exports.py OUTPUT/growth_0_seed_42/combined
python diagnostics/fitness_phase/analyse.py OUTPUT 42 297 100 0 OUTPUT/v1_analysis.json
python diagnostics/fitness_phase/analyse.py OUTPUT 42 57 200 1 OUTPUT/v3_analysis.json
pvpython --force-offscreen-rendering diagnostics/fitness_phase/pv_workflow.py build OUTPUT/growth_0_seed_42/combined VIEW
pvpython --force-offscreen-rendering diagnostics/fitness_phase/pv_workflow.py reopen RELOCATED_RUN RELOCATED_VIEW
pvpython --force-offscreen-rendering diagnostics/fitness_phase/pv_workflow.py build OUTPUT/growth_1_seed_42/combined VIEW_CLIP --clip
```

Here `python` means the existing coupling environment's interpreter, tested locally at `/home/aurascoper/Developer/Biofilms/coupling/.venv/bin/python`. Julia uses the repository project and existing local Manifest. Checkpoints use Julia Serialization and require the recorded Julia version, source SHA and config hash; they are not an interchange format. `resume.jl` verifies continuation against every stored reference state; it refuses mismatched identity. Failed output is retained for diagnosis and cannot be overwritten by the runner.

Use `--development --dense` only for a separate smoke seed; it marks results nonproduction. Export tests require `CPM_SMOKE_RUN` to point at a complete five-frame smoke output. VTK tests include real file corruption and asymmetric coordinate/Int64 fixtures. Relocation uses a fresh ParaView process and explicit reader remapping. No GUI is needed for simulation or readback.

## VTI schema and interpretation

Schema `cpm-sp-vti-1` uses Float64 scientific values, Int64 IDs and UInt8 masks, all **CellData**. N³ voxels require (N+1)³ points, x-fastest tuple order, origin at the lower corner and centers at origin + (index−0.5) spacing. Every frame sidecar binds dimensions, units, source, config, run identity, source-array hashes and file hash. The PVD records actual MCS with relative paths.

`cell_id`, `species_id`, `parent_id`, `lineage_id` use zero off-parcel; founders have parent zero. `interior_mask`, `occupied`, `alive`, `phase_valid`, `biomass_valid`, `local_coherence_valid`, `growth_enabled`, `survival_enabled` declare membership and modes. Environmental `dose_rate`, `ambient_cum_dose`, `nutrient`, `melanin`, `radiation_cpm` are fixed-voxel fields. Melanin and CPM radiation are legacy scores, not Gy. Dose arrays use Gy and Gy/s.

`parcel_cum_dose`, `cumulative_hazard`, `conditional_survival`, `biomass`, `target_volume`, rates and phase quantities repeat a parcel value on each occupied voxel. **Never sum repeated biomass.** Integrate `biomass_density` times voxel volume or sum the parcel table. Undefined floating values are NaN and have a validity mask; growth OFF biomass is unavailable, while the disabled growth contribution is explicitly zero. `phase_rad` is wrapped only for display; trajectory metrics retain unwrapped phases. Local coherence is the magnitude of the contact-weighted mean of neighbour phasors, excluding self. Isolates have undefined local coherence. Graph degree counts neighbours and contact strength sums weights.

The 2×2 ParaView layout shows species, net proxy, valid living phase, and dose. Named hidden toggles expose conditional survival and hazard. Fixed limits are species 1–7, net −0.06–0.05/MCS, cyclic phase 0–2π and dose 0–0.01 Gy/s. Cameras and time are common. `--clip` removes the camera-facing half with a whole-voxel plane cut through the axis, so the interior and the peak-dose axis are visible; the receipt records `axis_clip`. Voxels remain unsmoothed; IDs never pass through point interpolation. Companion plots show count, biomass, coherence, components, selection and unwrapped relative-phase change.

The analytic zero-K replay has its own ID and receipt. It uses identical exogenous membership and inheritance; selecting a measurement cohort does not remove other parcels from the coupled dynamics. Finite-window locking criteria and missing-data rules are in the preserved prediction. Optional radiation entrainment and feedback into survival are not implemented.
