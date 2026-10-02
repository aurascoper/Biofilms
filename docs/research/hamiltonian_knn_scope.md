# H_kNN semantic audit and observer-only feasibility

Implementation status, 2026-10-02: this document scopes work on the manuscript's `H_kNN` expression and implements none of it. The executed CPM acceptance path, the manuscript, the claims ledger and every diagnostic are unchanged. The substantive change is this document; the one other diff beside it is the citation-context metadata the repository requires for the code and text cited here.

This specification freezes the expression as written, audits the three names it carries, separates the manuscript's phase symbols from the one executable phase model in the tree, and states what each of three terminal outcomes would require. It does not define the missing symbols, commit to a replacement equation, train a model or report a result.

Status: design specification, 2026-10-02; not a production preregistration.

## 1. The written expression, frozen without repair

The manuscript introduces the expression with one sentence, preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`To model species transition probabilities modulated by neighbourhood phase states`, and writes it as

```
H_kNN = (1 / sigma_n) * sum_{j=1}^{n} P_sj(t) * F_j(t, x) * Gamma_s(t, x)
```

This is the historical specification under audit. The table records what the manuscript says about each symbol and what it leaves unsaid. It infers nothing: a blank cell is a finding, not an omission of this document.

| Symbol | Manuscript wording | Type implied | Dimensions | Clock | Spatial support | Audit status |
|---|---|---|---|---|---|---|
| `P_sj(t)` | a permutation matrix of species-s to species-j transitions in the symbol table; "the transition matrix" in the PSDE section | a 0/1 selection matrix, or a row of transition weights | none stated | depends on `t`, law unstated | none stated | type contradiction |
| `F_j(t, x)` | the fitness field of the PSDE section | scalar field per species | none stated | none stated | a field over `x` | written in the manuscript, absent from executed code |
| `Gamma_s(t, x)` | "phase-locked kernel for species s" | scalar per species and position | none stated | no evolution law | none stated | dynamics undefined |
| `n` | summation limit only | integer | none | none | none | undefined |
| `sigma_n` | appears only in this equation | normaliser | none stated | none | none | undefined |
| `j` | second index of `P_sj`, a species label in the symbol table | species label | none | none | none | not a geometric-neighbour index |
| `s` | the species the functional is written for | species label | none | none | none | defined |
| `t` | time argument | none | none | relation to the MCS unstated | none | clock undefined |
| `x` | position argument | none | none | none | lattice or continuum unstated | support undefined |

The symbol table calls `P_sj` a permutation matrix, preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`Permutation matrix describing species-$s$ to species-$j$ transitions`, while the PSDE section describes directed motility as preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`controlled by $\mu_s$ and the transition matrix $P_{sj}$`. Read literally, a permutation matrix has exactly one selected `j` in each row, so the sum collapses to one transformed fitness term; the manuscript's own audit paragraph says that it preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`would select one entry for each species`. Read as transition probabilities or mixing weights, it is not a permutation matrix. This document records the contradiction and chooses neither reading.

`Gamma_s` carries no `j` index. It is therefore a common multiplier of the whole sum, `H_kNN = (Gamma_s / sigma_n) * sum_j P_sj F_j`, and not a weight on individual terms. Unless `Gamma_s` itself already summarises the neighbourhood, which nothing in the manuscript specifies, the written equation does not weight individual `j` by phase state. And because `j` indexes species, the phrase "neighbour sum" names an intended interpretation that the equation does not define.

`F_j` is the fitness field of the PSDE section, which the manuscript records as preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`reported there as absent from this repository's simulation sources`. The aggregated quantity is therefore itself specified rather than implemented. `sigma_n` and `n` are defined nowhere in the manuscript, and no symbol in the expression carries a dimension, an evaluation clock or a lattice support.

## 2. Three names, audited separately

### kNN

Not descriptive of the written expression. The manuscript concedes this itself, preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`names an operation defined by a distance metric, a ranking and a cutoff`, and the claims ledger records the coinage as declared, data/claims_ledger.csv@`COINAGE DECLARED 2026-09-01`. The expression supplies no feature space, no distance metric, no ranking, no `k`, no nearest-neighbour selection, no training or reference observations and no classification or regression rule. The only k-nearest-neighbour code in the tree is `reactor_decision_tree.R`, which runs a three-neighbour classifier, reactor_decision_tree.R@`k <- 3`, over the cost and energy density of four fuel isotopes; it has no species, phase or lattice input and bears on none of this. This document records the conclusion rather than rediscovering it.

### Decision tree

Not descriptive. The expression has no predictors, no split predicates, no thresholds, no nodes, no leaves and no training rule. It is a single product under a sum, and there is nothing in it to split on.

### Hamiltonian

Not established by the written equation. The executed rule has four terms, preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`The active Metropolis rule contains four terms`, which the code sums in a fixed order, biofilms_potts.jl@`return t.adh + t.vol + t.rad + t.mel`. Those terms carry no physical unit: the adhesion matrix, the volume stiffness, the tropism coefficients and the melanin coefficient are model-scale numbers whose only common scale is the acceptance temperature `T_cpm`. A fifth term is therefore not asked for energy units. It is asked for five things. First, a scalar state functional `H_5(z, theta)` of the configuration `z`, with the phase `theta` held fixed during a proposal. Second, values in the same arbitrary acceptance-functional units as the existing terms, commensurate with `T_cpm`. Third, a site-copy change derived as an actual difference, `dH_5 = H_5(z_after, theta) - H_5(z_before, theta)`, and not written down as a score. Fourth, a declared ordering of the phase update and the mechanical proposal. Fifth, path independence: with phases frozen, a closed reversible sequence of site changes returns the accumulated `dH_5` to zero, up to a declared numerical tolerance. A proposal score that depends on the direction of the move or on proposal history fails the fifth requirement and must not be called Hamiltonian.

Until those exist, the neutral description is a candidate phase-weighted species-transition functional. This document may recommend that a later manuscript revision rename the expression; it performs no manuscript edit.

## 3. Two phase systems, not one

The manuscript layer is species-indexed. Its symbol table gives `theta_s` as preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`Adaptive phase-locking coefficient`, `Gamma_s` as preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`Phase-locked kernel for species $s$` and `Delta_s` as preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`Phase-locking synchronization deviation`, and its own audit records preprint/modeling_radioresistance_and_radiotropic_fitness.tex@`conflicting stated roles for $P_{sj}$ and $\theta_s$`. No evolution law for any of them is written.

The synthetic diagnostic layer is parcel-indexed and lives in `diagnostics/fitness_phase/`. Each persistent parcel carries one unwrapped phase and one intrinsic frequency, advanced synchronously on a degree-normalised contact graph by diagnostics/fitness_phase/phase.jl@`v+=sum(A[i,j]*K[i,j]*sin(theta[j]-theta[i]) for j in 1:n)/d`, in explicit Euler substeps within each MCS, beside an analytic zero-coupling replay, diagnostics/fitness_phase/phase.jl@`exact zero-coupling integral for one MCS`. The law is Kuramoto-type, with the repository's standing qualification that the attribution of such a law to the cited literature is unverified. The layer is read-only with respect to the CPM: the run records diagnostics/fitness_phase/run.jl@`phase_feedback=false`, and its paired arms must agree on every non-phase state hash. Its specification admits a centroid graph only as docs/research/phase_locking_spec.md@`A centroid k-nearest-neighbour graph is a separate, declared fallback` for data without contact geometry; that is a graph construction for the passenger and has no relation to the manuscript's `H_kNN`.

These two systems must not be identified merely because each contains an angle. Before any coupling into mechanics, five choices are declared explicitly: the state entity (parcel or species), the graph entity (parcel contact, species interaction or other), the phase value (`theta_i` or `theta_s`), the coupling clock relative to the mechanical MCS, and the rules for a new parcel and a removed parcel. The parcel passenger is the only executable candidate in the tree. Its recorded behaviour with growth off and with growth on, in `diagnostics/fitness_phase/PREDICTION.md`, is evidence about that diagnostic alone and supplies none of the semantics missing from the written expression.

## 4. Outcome A: semantic correction

Outcome A requires no new experiment. Its conditions are already met: the expression is not a k-nearest-neighbour operation, it is not a decision tree, its neighbourhood semantics are unresolved and its Hamiltonian status is not established. The scope verdict is that the current name must not be used as evidence about the mathematical operation it labels, and that a later manuscript-only revision may replace it with a neutral descriptive name once the intended object has been chosen. No rename occurs on this branch.

## 5. Outcome B: observer-only learning, null by construction

### 5.1 The formal null

For every evaluated serial proposal, let `D` be the sum of the four executed terms and `Y` the accept decision. The stepper draws biofilms_potts.jl@`u = ΔH <= 0 ? NaN : rand(rng)` and decides biofilms_potts.jl@`accept = ΔH <= 0 ? true : u < exp(-ΔH / p.T_cpm)`. The source site and its target are drawn uniformly from the same generator before any state feature is consulted, biofilms_potts.jl@`dx, dy, dz = NEIGHBORS_26[rand(rng, 1:26)]`. So for any feature `X` computed from the state before the decision,

```
P(Y = 1 | D, X, T_cpm) = p_0(D, T_cpm) = 1            if D <= 0
                                        = exp(-D / T_cpm) if D > 0
```

and therefore `Y` is independent of `X` given `(D, T_cpm)`; the conditional mutual information `I(Y; X | D, T_cpm)` is zero under the simulator's own rule. This is stronger than "a learned model probably will not help". `p_0` is the Bayes predictor of acceptance under any proper scoring rule, so a model that improves on it systematically on held-out simulator proposals indicates finite-sample fluctuation or one of the defects listed in 5.4, never phase physics. The same statement holds for the JACC port, where the recorded summed change and a fresh counter-based draw decide acceptance; a within-pass stale volume is inside that recorded value, not a feature beside it.

Excluded from `X` by construction: the uniform draw itself, the generator state, the seed, the JACC counter keys, any state after the decision, and any acceptance or driver label derived from the draw. The proposal index is retained as a join key only and is never a predictor.

### 5.2 Admissible predictors

Candidate diagnostic features are the parcel phase, the relative phase across the contact, the contact degree, the graph component, the species pair, the local contact composition and any other quantity computed from the state before the decision.

### 5.3 Baseline and scoring

The baseline is the analytic `p_0`, not a logistic model fitted to the four terms. The primary comparison is the Brier score difference on held-out seeds, because it stays finite where `p_0 = 1`. Secondary readouts are the log score where it is defined, the calibration of the residual `r = Y - p_0`, and held-out prediction of that residual. Under the null, `E[r | D, X] = 0`, so a residual predictor learned from `X` has zero held-out signal in expectation, which tests the additional-information question directly rather than asking one classifier to compete with a known generative law. Classification accuracy is unsuitable: thresholding discards the probability information the null specifies.

### 5.4 Interpretation

On ordinary simulator data the predicted result is an expected improvement of zero. Small positive or negative held-out differences are sampling variation. A reproducible held-out improvement is a diagnostic failure condition, investigated in this order: outcome leakage; generator or counter leakage; seed leakage between training and test; incorrect reconstruction of `D`; an incorrect `T_cpm`; an observer that changed the trajectory; drift in the acceptance rule; and the scoring implementation.

### 5.5 The observer that would be needed

The current hooks cannot produce the dataset. The proposal hook receives the terms and their sum only, biofilms_potts.jl@`on_proposal === nothing || on_proposal(terms, ΔH)`, fires before the draw, and is forbidden to consult the generator because biofilms_potts.jl@`The hook is called before the draw, so a closure that` touches it moves the trajectory. The accepted-copy hook carries the proposal index, biofilms_potts.jl@`proposal_index = proposal_index,`, but fires only for accepted copies. The replaceable energy function sees the state and both sites but not the decision. Repeated site pairs within a sweep make retrospective matching of these channels unsafe, so a wrapper around the energy function joined to the accepted channel does not yield an unambiguous accepted-or-rejected row for every evaluated proposal.

A future implementation therefore needs one immutable post-decision record for every evaluated proposal: the MCS, the proposal index, donor and recipient sites and species, the four terms, their sum, the draw or its recorded absence, the accept flag, and the pre-decision features captured separately. The observer receives no mutable state and no reference to the generator. Its controls inherit the pattern of `tests/rad_proposals_tests.jl`: bare and passive-observer runs agree exactly in state and generator bytes; the observer is shown to fire; a different seed changes the trajectory; and a deliberately generator-consuming observer changes the trajectory, as that suite already shows for the existing hook, tests/rad_proposals_tests.jl@`@test drew.lattice != bare_lat`, so the inertness test is one that can go red. The recorder's row count must equal the evaluated-proposal count, and its accepted rows must reconcile exactly with the accepted-event channel. None of this is implemented by this document.

## 6. Controls for outcome B

The planted positive control is built only in a scratch or mutant path, never in the shipped stepper. A synthetic feature `X*` is generated from its own declared deterministic source, independent of the four terms, and acceptance in the mutant uses `D* = D + g(X*)` for a known `g`, while the baseline under test is given `D` only. A diagnostic learner must recover held-out predictive value from `X*`, and a global permutation of `X*` must destroy that advantage by construction. This establishes that the analysis pipeline can see a fifth term when one exists.

On real simulator data there is no planted term and the expected result is no held-out residual signal. A stratified permutation of real phase or neighbourhood features is a secondary diagnostic only, because those features can correlate with `D` and a continuous `D` leaves residual information inside any bin. If it is used, the strata and the bins of `D` are fixed before any result is inspected, and the within-stratum imbalance of `D` is reported beside the result.

## 7. A scientific outcome B is gated on measurements

A learned relationship can be interpreted scientifically only on observations whose outcome was not generated by the four-term Metropolis rule. The minimum is measured transition or movement observations with a declared observational unit, a mapping from parcels to those observations, a spatial scale, a time mapping, a target variable that is not the simulator's own acceptance decision, and held-out biological observations. The repository supplies none of these in combination: the lattice pitch and the MCS-to-time mapping are both awaiting measurement, and `seconds_per_mcs` is undefined by contract. Simulator-based work under outcome B can therefore verify the recorder and the inference pipeline. It cannot discover an omitted biological acceptance mechanism from labels that the existing acceptance mechanism generated.

## 8. Outcome C: a mechanistic coupling is a new model

The factorial specification already states the bar: docs/research/survival_phase_factorial_spec.md@`a phase term in mechanics is a new coupled model requiring an equation, sign, units, update order, evidence boundary, positive control and separate preregistration`. Before implementation, outcome C requires the state variable and its entity level; a global scalar functional; the acceptance-unit scale; a coupling coefficient `kappa` with its sign; the spatial or contact support; the ordering of mechanical and phase updates; division and removal semantics; a restart schema; a serial kernel and a JACC kernel; an exact negative control at `kappa = 0`; a planted positive control at nonzero `kappa`; the closed-loop consequence stated; and a prospective effect criterion.

One complete example, offered only to make the missing choices visible and not as the selected model, is a symmetric phase-modulated contact term over unordered contacting parcel pairs:

```
H_phi = -kappa * sum_{<i,j>} w_{s_i s_j} * cos(theta_i - theta_j)
```

A site-copy proposal then changes only the boundary contacts of that site, with phase frozen during the decision. The example exposes every choice the written expression leaves open: what constitutes a contact; whether contacts count per face or once per parcel pair; the species weight `w`; the units of `kappa`; the update order; serial and JACC parity; and the global functional from which the site-copy change actually derives. Its integrability control is the one named in section 2: with phases frozen, a closed reversible sequence of configuration changes returns the accumulated change to zero. And the moment contacts influence phase while phase influences the acceptance of contacts, the former passenger is a feedback system, and it needs a separate preregistration.

## 9. Terminal decisions

| Outcome | Current status | What this scope permits |
|---|---|---|
| A, semantic rename | conditions already established | recommend a later manuscript-only correction |
| B, diagnostic learning | the null is a theorem on simulator labels | a future recorder and pipeline check only |
| B, scientific learning | blocked | requires independently measured outcomes |
| C, fifth acceptance term | unspecified | requires the complete mechanistic specification of section 8 and a new preregistration |

Nothing in this scope enters `compute_delta_H_terms`, `compute_delta_H`, `mcs_step!` or the JACC acceptance kernel.
