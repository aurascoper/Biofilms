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

## Gates (ALL must clear before implementation)

| # | Gate | Status (2026-08-13) |
|---|---|---|
| 1 | Restart/window semantics proven | **CLEAR** — windowed runs bit-identical to unbroken runs; restart checkpoint resume bit-identical (tests/checkpoint_io_tests.jl) |
| 2 | Membrane m-vs-P_eff constitutive choice declared | **BLOCKED** — open contradiction, audit §6: `P_eff` grows with dose while `m` decays; the intended coupling is undeclared. STOP CONDITION: no OpenMC dose into radiodialysis until the user selects it |
| 3 | Membrane dose statistic declared & tested | **BLOCKED** — config key exists (`membrane_statistic`), value unset |
| 4 | Physical inputs supplied | **BLOCKED** — voxel pitch, s/MCS, densities/compositions, spectrum, source rate all unset (intentional; audit §5) |
| 5 | Offline feedback gate exceeded uncertainty + declared threshold | **NOT EVALUATED** — `biofilm-dose-scan` reports NOT EVALUATED while gate 4 is blocked |

Until gates 2–5 clear, running "coupled" dynamics would mean cells responding
to an OpenMC field while membrane integrity responds to the hard-coded
`Ddot_R = 1.0 Gy/s` — two incompatible radiation models in one simulation.
