# Phase locking as a passenger diagnostic

This specification defines a phase calculation that can be attached to a CPM run without changing that run. It does not implement the phase layer, report a phase result, or discretise the manuscript's fitness-field equation.

Status: design specification, 2026-09-29; paired invariants and selection analysis refined 2026-09-30. The reported CPM has no phase state. The first implementation must remain diagnostic: it reads parcel geometry after each completed CPM step and cannot write the lattice, fields, acceptance terms, parcel records, checkpoints or mechanical random-number generator.

## 1. State and configuration

Store one phase `theta_i` and one intrinsic angular frequency `omega_i` in a side table keyed by persistent parcel ID. Do not add either value to `CellInfo`. A parcel ID missing from the table is an error while the layer is active. A removed parcel leaves the active table and may be copied to a diagnostic archive; its value never affects the CPM.

`PhaseConfig` is immutable and has no biological defaults. It must require an explicit time step, initial phases, intrinsic frequencies, graph rule, coupling matrix and radiation-coupling mode. Time, frequency and coupling must use one declared unit system. Every scalar and matrix entry must be finite. The time step must be positive. Couplings may be zero but cannot be inferred from the CPM's adhesion matrix. The only accepted radiation-coupling mode in this version is `none`.

Attach and update operations receive the CPM state as read-only input. With the layer attached or absent, a run with the same initial state and mechanical seed must produce byte-identical CPM state, event output and mechanical RNG state at every checkpoint.

The factorial requires two explicit comparisons at each completed MCS and after restart: baseline versus phase only must have identical survival and CPM trajectories; survival only versus combined must have identical deaths, dose histories and CPM trajectories. Phase must not consume either the mechanical or hazard RNG. The full per-ID comparison and selection controls are specified in `survival_phase_factorial_spec.md`.

## 2. The update is synchronous and degree normalised

For old phases `theta(t)`, use

`theta_i(t+dt) = wrap(theta_i(t) + dt * [omega_i + sum_j A_ij K_s(i),s(j) sin(theta_j(t)-theta_i(t)) / d_i])`,

where `d_i = sum_j A_ij`. Every right-hand-side phase is read from the old vector, and all new phases are committed together. An isolated parcel has zero coupling and advances only at its intrinsic frequency. The graph normalization is `D^-1 A`; it is not the raw adjacency matrix. Noise is absent. A later noise term needs its own distribution, seed, update order and checkpoint state.

The reference graph is the undirected parcel contact graph induced by the CPM adhesion neighbourhood after a completed mechanical step. Two distinct parcels are adjacent when at least one lattice-neighbour pair crosses their boundary. Multiple touching faces do not create multiple edges in this first model.

A centroid k-nearest-neighbour graph is a separate, declared fallback for data without contact geometry. Symmetrise it before updating: an edge exists when either endpoint selects the other. The graph rule is fixed before a run. It cannot change because a result is weak or a graph is disconnected.

The order parameter is `r(t) = abs(sum_i exp(im*theta_i(t))) / N` over the declared active population. Report effective frequencies over a declared late-time window as phase increments before wrapping divided by elapsed time. Report their spread together with `r`; either one alone can mislabel a transient alignment as locking.

Record survivor IDs and count, graph edges, component membership and sizes, isolates and largest-component fraction alongside coherence and pairwise relative-phase drift. Preserve unwrapped phase paths for drift and slip analysis. Report both global and component coherence. When mortality is enabled, compute coherence before removal, after removal with phases held fixed, and after phase evolution. Also analyse fixed survivor cohorts against a zero-coupling replay using the same geometry and removals. The factorial specification defines these estimators and their missing-data rules. A coherence increase caused by removal of discordant parcels cannot establish dynamical locking.

## 3. Reference calculation and positive controls

For an all-to-all graph and a declared frequency distribution, compute the finite-sample self-consistency reference from the exact frequencies used in the run. An infinite-population threshold is context, not the acceptance test for a small moving graph.

The two-oscillator control uses

`d(delta)/dt = delta_omega - 2 K sin(delta)`.

A fixed point exists when `abs(delta_omega) <= 2 K`. At equality it is marginal, so the numerical control must distinguish convergence below the boundary from slow or neutral behaviour at it. A second control uses a star graph to reject raw-degree coupling: leaves and hub must receive the declared degree-normalised sums. A third uses phases for which an in-place Gauss--Seidel update differs from the synchronous result.

The prospective exercise varies a scalar multiplier `K` over a fixed grid while holding frequencies, initial phases, graph rule and CPM trajectories fixed. Its pre-registered checks are:

- the rank correlation between seed-averaged late-time `r` and `K` has the declared positive sign;
- on the all-to-all reference, the observed onset lies within a factor of two of the finite-sample self-consistency prediction;
- on the contact graph, the late-time spread of effective frequencies falls by at least 50% from the zero-coupling control.

These are synthetic implementation and behaviour checks. They do not establish a microbial oscillator, entrainment by radiation, or a biological phase-locking coefficient.

The evidence gates behind that boundary are recorded in `phase_evidence_brief.md`. Liu et al. (2017) confirm a two-biofilm coupled-phase-oscillator precedent in *Bacillus subtilis*, but the accessible main article does not expose the coupling equation. Until its supplementary equation is retrieved, the specific label `Kuramoto-type` is unverified. The same brief derives the infinite all-to-all normal-frequency result `K_c = sigma sqrt(8/pi)` and keeps it separate from the exact finite-sample reference required here.

## 4. Tests that must fail for defective implementations

The implementation test set must include controls that reject:

1. a mutable or default-filled configuration;
2. nonfinite phases, frequencies, couplings or time steps;
3. zero or negative time step;
4. any radiation-coupling mode other than `none`;
5. a missing or duplicate parcel ID;
6. self-edges or a directed contact graph;
7. raw adjacency in place of `D^-1 A`, using the star graph;
8. nonzero coupling for an isolated parcel;
9. an in-place Gauss--Seidel update;
10. a phase wrap that changes the unwrapped increment used for effective frequency;
11. an incorrect two-oscillator locking boundary;
12. a k-nearest-neighbour fallback that is not symmetrised;
13. a graph-rule switch during a run;
14. reuse of the mechanical RNG by the diagnostic;
15. mutation of CPM state, fields, events or checkpoint bytes;
16. a status or manuscript guard that calls the layer exercised before a registered run exists;
17. reuse of the hazard RNG or a phase-induced change to deaths or dose histories in the combined arm;
18. a locking verdict based only on rising coherence after frozen-phase removal;
19. omitted isolates or changing survivor/pair sets within a fixed-cohort comparison; and
20. a singleton locking verdict or a drift calculation that misses an unwrapped phase slip.

Each test needs a mutant or bad input that reaches the production entry point and is rejected. A control that calls only a helper does not protect attachment to the CPM.

## 5. Checkpoint and reporting boundary

The first passenger implementation need not alter CPM checkpoints. Its own restart record must contain the source revision, phase schema, full `PhaseConfig`, graph rule, active and archived ID tables, old phase vector, intrinsic frequencies, diagnostic step index and any diagnostic RNG state added later. Include unwrapped phases and a verified reference to the phase, graph and removal records needed to reproduce windowed analysis across the restart. Restore must validate IDs against the CPM checkpoint before either state becomes runnable.

The manuscript may describe this specification as a prospective diagnostic. It must continue to say that the reported simulation did not execute phase dynamics. The words `exercised`, `validated` and `radiation entrainment` remain unavailable until their separate evidence requirements are met.
