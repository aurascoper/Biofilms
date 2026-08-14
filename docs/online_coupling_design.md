# Online coupling: design and gate status (NOT IMPLEMENTED)

The online (two-way) driver is deliberately absent from this branch. This
document records the agreed design so implementation is mechanical once the
gates below clear — and records why implementing it now would be premature.

## Design (windowed API, Julia as state authority)

The restart-file round-trip architecture was rejected in review: a portable
snapshot cannot resume the serial stochastic process, and shuttling full
restart checkpoints per window is wasteful. Instead, one long-lived Julia
process owns the exact simulation state (`CoupledSimulation`: CPM state,
radiodialysis state, contaminant field, the MersenneTwister, clocks), all of
which exists and is tested as of commit 3:

```
loop (Python orchestrator, Julia state authority):
    Julia:  export_transport_snapshot(sim, t0)          # labels at t0
    Python: needs_rerun? → OpenMC(state at t0)          # exact-hash policy
    Python: transport_result + dose_attribution
    Julia:  import_dose_field!(state, dose_mean, dose_sd;
                hamiltonian_transform, melanin_transform)   # config-declared
    Julia:  advance_window!(sim, n_window)              # accrues dose per-MCS
    repeat with t1 = t0 + n_window · seconds_per_mcs
```

Operator order is first-order explicit and starts with transport — the first
biological window is never advanced on a stale or legacy field (review
amendment 5). Dose is attributed to *current* voxel labels every MCS inside
`advance_window!` (`accrue_dose!`), so exposure during a held field is never
assigned to endpoint geometry.

Membrane update: `membrane_dose_rate_Gy_s` is computed by the config-declared
statistic (mass-weighted volume mean or area-weighted surface measure) and
would replace the legacy constant `Ddot_R` — but see gate 2.

## Gates (status by stage, 2026-08-14)

Gates are reported per stage, because requirements are staged
(`config.STAGES`: transport → dosimetry → cpm_feedback → membrane_feedback).
A transport-only experiment is not "missing" a CPM clock it never uses, and
reporting it that way hid the fact that the transport blocker was artificial.

| Gate | Fixed-membrane / transport stage | Dose-responsive stage |
|---|---|---|
| 0 · Restart/window semantics proven | **CLEAR** — windowed runs bit-identical to unbroken runs; restart resume bit-identical (`tests/checkpoint_io_tests.jl`) | CLEAR (same evidence) |
| 1 · Membrane constitutive choice | **CLEAR** — membrane declared fixed: `m_mech = 1`, `P_j = P_j0`, no dose→membrane update (`docs/physical_reference_system.md` §3) | **BLOCKED** — `m` vs `P_eff` unresolved in code (audit §6); **STOP CONDITION RETAINED: no OpenMC dose enters the radiodialysis model until the constitutive choice is explicitly selected** |
| 2 · Membrane dose statistic | **NOT APPLICABLE** — diagnostic only while the membrane is fixed; the config key stays required at `membrane_feedback` and unset | **BLOCKED** — `membrane_statistic` unset |
| 3 · Transport geometry / materials / spectrum | **IN PROGRESS** — per reference system; `data/parameter_provenance.csv` tracks each key. A0 is buildable as of 2026-08-14 (`config/reference_a0_water_phantom.toml`): spectrum and geometry ready, biomass/membrane now `not_applicable` rather than blocked (water-phantom gap closed). Mesh resolution and history count are `provisional` pending the A0 sweep | Required |
| 4 · Source activity | **NOT NEEDED** — heating is tallied per source particle and every sensitivity metric is a ratio, so `biofilm-dose-scan --stage transport` runs without it | **BLOCKED** — needs an assay certificate with a reference date; a catalogue activity range is not an activity |
| 5 · Seconds per MCS | **NOT NEEDED** — no CPM is advanced | **BLOCKED TWICE OVER** — `Δt_MCS = a²·S_sim/S_exp` needs the calibrated pitch **and** a declared dynamic observable the model represents (gate 9). Single-cell tracking is not automatically valid: a cell ID is a biomass parcel, so a per-cell MSD is not a parcel MSD |
| 6 · Biological response transforms | **NOT NEEDED** | **BLOCKED** — four distinct response functions, none calibrated. Each is a PER-PARCEL response, not per-organism (gate 8). `response.growth_survival` is `unsupported_by_current_model` rather than blocked: it awaits a MODEL CHANGE, not a measurement, because the CPM has no growth dynamics and no quantity of survival data would clear it |
| 7 · Offline feedback gate exceeds uncertainty + threshold | **NOT APPLICABLE** to transport validation | **NOT EVALUATED** — reported while gates above are open |
| 8 · Voxel binding: what an occupied voxel is and is made of | **NOT APPLICABLE** — a water phantom has no occupied voxels | **PROVISIONAL** — the binding is declared and *coherent* (one occupied voxel is part of a hydrated computational biomass parcel, assigned the hydrated-effective-biofilm material at a wet bulk density) but both its numbers are unset. Spatial gate `PROVISIONAL`, material OpenMC gate `PROVISIONAL`, material radiodialysis gate `BLOCKED` by the site-occupancy units defect. `emit_transport_config()` refuses. See `docs/calibration/integration_contract.md` |

The reason the stop condition survives the fixed-membrane declaration: fixing
the membrane removes the feedback path, it does not choose between the two
contradictory laws. Both still stand in `biofilms_potts.jl` and would
reactivate the moment feedback is wired.

| 9 · Dynamic observable for the time calibration | **NOT APPLICABLE** — no CPM is advanced | **BLOCKED** — seven candidates declared, four satisfy all three conditions (represented by the model, measurable, semantically compatible with a parcel), none selected. Growth-dominated observables are excluded on principle: the CPM has motility and shape rearrangement but no nutrient-driven biomass creation. See `docs/calibration/time_observable_contract.md` |

Gate 8 is where the two calibration branches meet, and it exists because a
binding can be incoherent in ways neither branch can detect alone: an
organism-scale entity assigned a bulk-biofilm material passes both gates
separately and contradicts itself jointly.

Until the dose-responsive column clears, running "coupled" dynamics would mean
cells responding to an OpenMC field while membrane integrity responds to the
hard-coded `Ddot_R = 1.0 Gy/s` — two incompatible radiation models in one
simulation.
