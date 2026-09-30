# Prospective survival and phase factorial

This specification defines how to exercise the synthetic survival analogue and the phase passenger separately and together. It does not implement either layer or report a result.

Status: design specification, 2026-09-30; not a production preregistration. Production requires a clean revision containing both implementations, their tests, a fixed parameter manifest, an analysis program and a `PREDICTION` file written before any production output is inspected.

## 1. Four matched arms

| Arm | Synthetic survival | Phase passenger |
|---|---|---|
| Baseline | Off | Off |
| Survival only | On | Off |
| Phase only | Off | On |
| Combined | On | On |

All arms use the same initial lattice, parcel IDs, mechanical seed, physical dose map, movement parameters and observation times. The survival-only and combined arms also share the same hazard seed. The phase update is deterministic in the first model. If phase noise is later added, it receives a third generator whose state is independent of the mechanical and hazard generators.

The baseline and disabled-module controls must reproduce the original CPM state, event output and mechanical random-number-generator bytes. Attaching the phase passenger cannot consume the mechanical or hazard streams.

## 2. The first combined model has one directed connection

A survival commitment makes a parcel nonviable and removes its ID from the active phase graph at the declared lifecycle point. Its phase may be archived for reporting but cannot influence another update after removal. This is the only connection between the modules in the first combined model.

Phase cannot change absorbed dose, hazard probability, survival status, target volume, lattice moves, field production or Metropolis acceptance. Therefore, under matched initial state, physical field, mechanical seed and hazard seed:

- the combined and survival-only arms must have identical hazard commitments, removals, dose histories and mechanical trajectories; and
- the combined phase result may differ from phase-only because selective loss changes the active node set and contact graph.

Failure of the first identity is an implementation error. A phase-dependent hazard or a phase term in mechanics is a new coupled model. It requires an explicit equation, sign, units, update order, evidence boundary, positive control and separate preregistration.

## 3. Outcomes remain separate

Report the survival specification's viable fraction, survivor-weighted radius, initial radius of the eventual survivors and identified-parcel displacement. Report the phase specification's order parameter, effective-frequency spread, active node count and graph components. Never substitute survivor radius for movement or phase coherence for survival.

For the combined arm, report phase summaries on the active set and state how node loss changes the denominator. If no active phase nodes remain, phase summaries are missing and `zero_active_phase` is recorded. Do not assign zero coherence.

The first factorial tests three questions:

1. Does the survival implementation change survival outcomes while preserving the declared movement decomposition?
2. Does the phase passenger reproduce its synthetic locking controls without changing the CPM?
3. Does selective loss alter the phase diagnostic through the declared active-graph update while leaving survival and mechanics identical to survival only?

It does not test radiation entrainment of phase, a biological phase-survival interaction, or the manuscript's fitness-field equation.

## 4. Required controls before production

- Zero dose, zero sensitivity and disabled survival produce no hazard commitments.
- The fixed-position survival control can change survivor radius while identified-parcel displacement stays zero.
- The two-oscillator, star-graph, isolated-node and synchronous-update phase controls pass.
- Phase-only and baseline have identical CPM state and mechanical generator bytes.
- Combined and survival-only have identical survival events, mechanical trajectories and hazard-generator state.
- Restart at a declared checkpoint reproduces every later survival and phase record.
- Missing, duplicate or inconsistent parcel IDs stop the run before either extension mutates state.
- The production analysis rejects incomplete arms, mixed revisions, changed thresholds and nonfinite values.

The production `PREDICTION` must fix the seed range, horizon, survival vector, clock, dose-rate regime, phase frequency sample, coupling grid, graph rule, missing-data rules and decision thresholds. Synthetic values must remain labelled assumed.
