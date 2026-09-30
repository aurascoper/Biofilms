# Prospective survival and phase factorial

Implementation status, 2026-09-30: the original design below is retained as a historical specification. The closed synthetic implementation, changes to its checkpoint/graph contract, frozen prediction and execution receipts are under `diagnostics/fitness_phase/` and `fitness_phase_implementation_review.md`. The legacy movement study and the original undefined PSDE terms retain their stated boundaries. Production completed from source 841d052 on 2026-09-30; results are appended in `diagnostics/fitness_phase/PREDICTION.md` and summarised in the manuscript appendix on execution.

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

| Paired arms | Required identity |
|---|---|
| Baseline versus phase only | Identical per-ID survival status and event histories, dose histories and CPM trajectories, including lattice, fields, parcel records and mechanical RNG bytes |
| Survival only versus combined | Identical per-ID deaths (hazard commitments), physical removals, incremental and cumulative dose histories, CPM trajectories, and mechanical and hazard RNG bytes |

Compare the paired records exactly at each completed MCS and after restart, not just final summaries. Phase state is excluded from the CPM comparison. Different seed values alone do not establish stream isolation: attachment and updates must leave both process generators untouched. Each identity needs an enabled phase test, including a nonzero-death case for survival only versus combined.

Failure of either identity is an implementation error. Combined phase output may differ from phase only because mortality changes the active graph. A phase-dependent hazard or a phase term in mechanics is a new coupled model requiring an equation, sign, units, update order, evidence boundary, positive control and separate preregistration.

## 3. Outcomes remain separate

Report viable fraction, survivor-weighted radius, initial radius of eventual survivors and identified-parcel displacement as defined in the survival specification. At every phase observation, record survivor IDs and count, coherence, relative-phase drift and effective-frequency spread. Record graph edges, component membership and sizes (including isolates), isolate count and largest-component fraction. Report coherence within components as well as globally; disconnected components need not share a phase. Never substitute survivor radius for movement or coherence for survival.

### Separate selection from phase evolution

After survival commitment and physical cleanup, refresh geometry and the viable graph, then advance phase synchronously. A dying parcel leaves the phase graph at commitment even while its biomass remains. Record every active-ID loss and its cause, including ordinary mechanical removal.

For each interval let `S_old` be the previous active set, `S_new` the remaining set and `theta_old` the phases before the next update. With `r(S, theta) = abs(sum(exp(im*theta_i) for i in S))/length(S)`, record:

- selection: `r(S_new, theta_old) - r(S_old, theta_old)`;
- evolution: `r(S_new, theta_new) - r(S_new, theta_old)`.

Their sum is the total coherence change. The first holds all phases fixed. The second measures evolution on the retained set, but is not by itself evidence of locking. Additions are outside this no-birth design and must be refused.

For each preregistered analysis window, also compute coherence and pairwise relative-phase drift on the fixed IDs present throughout that window. Save unwrapped phase paths: for each fixed pair report `(delta_ij(t1)-delta_ij(t0))/(t1-t0)` and the range of `delta_ij(t)`, where `delta_ij = theta_i_unwrapped - theta_j_unwrapped`. Keep the pair set fixed across edge changes. Endpoint drift alone can hide intervening slips. This is a survivor-conditioned diagnostic, not an unbiased estimate for the starting population.

Replay phase with zero coupling on the same recorded geometry and removal schedule, with identical initial phases and frequencies. Compare the same survivor cohort and pairs. A locking interpretation requires the declared drift and frequency-spread criteria over the declared window, and comparison with this control. Increased active-set coherence alone is insufficient; equal intrinsic frequencies can also give zero drift without coupling.

If no active nodes remain, phase summaries and the decomposition are missing with `zero_active_phase`. For one node, report coherence as 1 but pairwise drift and locking as not assessable. An empty or singleton fixed cohort has the same limits. Never replace missing values with zero.

The first factorial tests three questions:

1. Does the survival implementation change survival outcomes while preserving the declared movement decomposition?
2. Does the phase passenger reproduce its synthetic locking controls without changing the CPM?
3. Does selective loss alter the phase diagnostic through the declared active-graph update while leaving survival and mechanics identical to survival only?

It does not test radiation entrainment of phase, a biological phase-survival interaction, or the manuscript's fitness-field equation.

## 4. Required controls before production

- Zero dose, zero sensitivity and disabled survival produce no hazard commitments.
- The fixed-position survival control can change survivor radius while identified-parcel displacement stays zero.
- The two-oscillator, star-graph, isolated-node and synchronous-update phase controls pass.
- Both paired identities in Section 2 pass; mutants consuming either process RNG or changing a survival, dose or CPM record must fail through the integrated runner.
- A frozen-phase fixture `(0, 0, pi)` losing its third parcel raises coherence from `1/3` to `1`, with selection contribution `2/3` and zero evolution contribution. It cannot establish coupling-induced locking.
- A known locking oscillator control with a fixed cohort is detected; a drifting-pair control with wrapped endpoint agreement is rejected.
- The analysis rejects changing pair denominators, omitted isolates, singleton locking and substitutions of zero for missing phase values.
- Restart at a declared checkpoint reproduces every later survival and phase record.
- Missing, duplicate or inconsistent parcel IDs stop the run before either extension mutates state.
- The production analysis rejects incomplete arms, mixed revisions, changed thresholds and nonfinite values.

The production `PREDICTION` must fix the seed range, horizon, survival vector, clock, dose-rate regime, phase frequency sample, coupling grid, graph rule, analysis windows, cohort and pair rules, missing-data rules and decision thresholds. Synthetic values must remain labelled assumed. These checks are required contracts; none is claimed implemented or passed here.
