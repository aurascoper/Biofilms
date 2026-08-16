# The feedback gates — offline counterfactual, online authorization

**Date:** 2026-08-15 · **Status:** specified and synthetically testable; both verdicts
`NOT_EVALUATED` · **Engine:** `coupling/biofilm_openmc/feedback_gate.py` ·
**Effects:** `coupling/biofilm_openmc/feedback_uq.py` ·
**Ledger:** `data/uncertainty/feedback_parameter_distributions.csv` ·
**Verdicts:** `data/uncertainty/feedback_gate_verdicts.csv` (empty) ·
**Policy:** `config/feedback_effect_policy_template.toml` (threshold unset)

One document rather than two, because the two gates share one verdict vocabulary, one
decision rule and one uncertainty taxonomy, and splitting them would duplicate all three.
They are separate *gates*, not separate *specifications*.

## Why there are two

Conflating them would make **"feedback matters"** indistinguishable from **"the software
can perform feedback"**, which is the last distinction this repository has left to protect.

| | Offline | Online |
|---|---|---|
| Question | Does changing a material state move the field more than noise? | May the loop change state? |
| Kind | A counterfactual transport experiment | An authorization condition |
| Operates on | Immutable snapshots | The live simulation state |
| Touches biology | Never | Only through a licensed transform |
| Needs source activity | **No** — see below | Yes, for absolute Gy/s |
| Needs a response posterior | No | Yes |

The offline gate is worth running **before** any dose→biology response function exists,
because it answers a question that does not need one: *if the biofilm moved from its
measured baseline to a physically plausible melanized or metal-loaded state, would OpenMC
predict a change large enough to matter at all?* If the answer is no across the plausible
material envelope, expensive two-way coupling is not justified, and that is worth knowing
before the response fitting is paid for rather than after.

## The decision rule

$$q_{0.01}(s\,\Delta Q) > \delta \qquad\text{equivalently}\qquad P(s\,\Delta Q > \delta) \ge 0.99$$

- $\Delta Q$ is a **paired** effect, $Q(M_1) - Q(M_0)$, on the same snapshot. Not two
  independent results with error bars to compare — subtracting separately-reported means
  and adding their uncertainties in quadrature discards exactly the correlation a paired
  design exists to exploit.
- $s$ is the **declared direction**. A two-sided test passes on an effect of the wrong sign.
- $\delta$ is the **declared scientific effect threshold**, and it **has no default**.
- The rule is on the **whole distribution**, not a Gaussian half-width. The model is
  nonlinear, several inputs are bounded or compositional, and the output need not be
  symmetric. Quadrature is not offered because several of this project's uncertainties are
  correlated, and `uncertainty.propagate` says so of its own inputs.

**Nobody may default $\delta$.** No uncertainty calculation can decide what magnitude of
feedback would matter biologically or operationally. Until it is declared, the verdict is
`NOT_EVALUATED` — and a threshold declared *after* results are inspected is refused
outright, because it is not a threshold.

Cohen-style standardized effect size is reported and is **not** the acceptance rule: a
physically irrelevant effect can have a huge $d$ when uncertainty is small, and an
operationally important one can have a modest $d$ under a wide posterior. The gate needs
magnitude *and* certainty, which is why $\delta$ exists separately.

## Not every uncertainty is a probability

The load-bearing column in the ledger is `sampling_role`, not `distribution_family`.
Sampling every uncertain thing from one undifferentiated design produces a probability with
no physical interpretation.

| Uncertainty type | Role | Why |
|---|---|---|
| `physical_measurement`, `biological_posterior`, `source_metrology` | `outer_random` | Genuine epistemic uncertainty |
| `monte_carlo_estimator` | `inner_transport_replicate` | Estimates $\operatorname{Var}_\xi(\Delta Q \mid \theta)$; a seed is not a prior |
| `numerical_convergence` | `convergence_axis` | A mesh factor is an engineering decision. Bound it; never draw it |
| `discrete_model_choice`, `nuclear_data_model` | `scenario_branch` | A categorical choice gets a named branch, not an invented categorical probability |

`physical_contract.distribution_row_problems()` **refuses** `outer_random` for the last
three, so this is enforced rather than written down and hoped for.

Three consequences worth stating plainly:

- **Composition is compositional data.** Mass fractions must stay non-negative and sum to
  one, so independent Gaussian perturbation is invalid. Oxygen-by-difference creates strong
  negative correlations that must survive into the transport ensemble.
- **Segmentation is a field uncertainty.** Alternative plausible masks move biovolume,
  surface geometry, the selected pitch, CSG placement and per-label dose simultaneously, so
  the propagation object is a *set of fields* sharing specimen identifiers, not a scalar.
- **No Gaussian cross-section uncertainty is invented.** No usable photo-atomic covariance
  package was identified, so the pinned library stays for reproducibility and a library
  change is a named sensitivity scenario.

## Where variance comes from, and what fixes it

$$\operatorname{Var}(\Delta Q) = \underbrace{E_\theta[\operatorname{Var}_\xi(\Delta Q \mid \theta)]}_{\text{more seeds and histories}} + \underbrace{\operatorname{Var}_\theta[E_\xi(\Delta Q \mid \theta)]}_{\text{only better measurements}}$$

Reported rather than collapsed, because the two terms are reduced by completely different
and differently priced actions. A budget that reports one number cannot say which to buy.

Mesh bins are **not** independent — one history deposits in several — so per-label
uncertainty comes from independent replicate runs, exactly as `lineage.py` already insists.
Nothing sums per-bin variances as if they were independent.

**Seeds.** The relative error of a standard deviation estimated from $r$ replicates is about
$1/\sqrt{2(r-1)}$: ~35% at 5, ~16% at 20. Five spot a large effect; they are weak evidence
for a tight budget. `seed_sufficiency()` reports this so a number never stands in for its
own precision.

**Common random numbers are adopted only if measured.** Matched seeds cancel common noise
when the two material states transport similarly, but a large material change makes
histories diverge. Assuming the benefit would understate exactly the effects that matter
most, so `common_random_number_benefit()` measures it.

## The numerical budget

A pass driven by barely resolved numerics is a pass about the mesh. Fractions of $\delta$,
and conservative **project defaults** rather than metrological constants:

| Contribution | Ceiling |
|---|---|
| transport + numerical, 99% half-width | $0.25\,\delta$ |
| material-volume raytrace | $0.10\,\delta$ — the CSG mass is denominator-critical |
| surrogate approximation, if used | $0.10\,\delta$ |

**The mass denominator is not an ensemble member.** Production must use
`exact_csg_raytrace`. The full-bin method overstates A0's cylinder by exactly $4/\pi$
because it weighs the circumscribing cube; it is a **known negative control**, and averaging
it with the exact method averages in a known error.

An uncharacterised budget **blocks**. It is not assumed to be zero.

## Source activity cancels — sometimes

For fixed spectrum and geometry, $D_{\text{rate}} = D_{\text{per source}} \cdot \dot N_\gamma$,
so a relative effect $(Q_1 - Q_0)/Q_0$ cancels a common source-rate multiplier. **The offline
material screen can therefore run per source particle, before any assay certificate exists**
— which `drivers.compare_fields` already records of its own metrics.

Absolute Gy/s, time-integrated dose and biological response fitting get no such
cancellation and need the activity distribution with its reference date and decay.

## The state machine

One vocabulary covers both gates, so a blocked verdict says *which* thing is wrong:

```
NOT_EVALUATED                     BLOCKED_ON_PHYSICAL_CALIBRATION
BLOCKED_ON_NUMERICAL_RESOLUTION   BLOCKED_ON_BIOLOGICAL_CALIBRATION
OFFLINE_EFFECT_BELOW_THRESHOLD    OFFLINE_UNCERTAINTY_TOO_LARGE
OFFLINE_PASS                      ONLINE_OUT_OF_DOMAIN
ONLINE_UNCERTAINTY_TOO_LARGE      ONLINE_ENABLED
UNSUPPORTED_BY_CURRENT_MODEL
```

`OFFLINE_EFFECT_BELOW_THRESHOLD` and `OFFLINE_UNCERTAINTY_TOO_LARGE` are deliberately
distinct: the first is small however well resolved, the second may not be. Conflating them
hides which one more histories can fix.

**Everything fails closed.** Every `OnlinePrerequisites` field defaults to its refusing
value, so a caller who forgets one cannot authorize feedback by omission. A missing
threshold, an empty distribution, an unresolved model branch, or a state outside the
validated envelope all block. `ONLINE_ENABLED` requires, independently: the offline gate
passed, physical calibration and the production resolution decision ready, the biological
posterior and `seconds_per_mcs` calibrated, the current state inside the validated envelope,
the numerical budget met, and the effect *still* clearing $\delta$ at the current state.

**Julia cannot act on a decision nobody made.** The authorization lives in Python; the CPM's
obligation is to be structurally unable to apply an unauthorized transform.
`import_dose_field!` with no transform changes no radiation signal, no melanin drive, no
lattice and no accumulated dose, and `accrue_dose!` refuses without a clock. Importing is
not authorizing, and `tests/genealogy_tests.jl` asserts it.

## What is deliberately excluded

- **Growth and survival.** The CPM has no birth and no death, so radiotrophy as enhanced
  growth and radioresistance as survival are structurally absent. A model-revision problem,
  not a missing-data one; `evaluate_online` refuses the pathway outright.
- **Membrane response.** The m-versus-$P_{\text{eff}}$ constitutive contradiction is
  unresolved. Fixed-membrane one-way dosimetry does not need it; any dose-responsive
  feedback does.
- **Melanin as a fixed coefficient.** Radioprotection is not a universal monotonic
  multiplier — it varies with chemical composition, spatial arrangement, metabolic state and
  radiation modality. It enters as a bounded distribution over a scenario. Note also that
  only the *product* of melanin production and coupling reaches the dynamics, so those two
  are jointly unidentifiable and are grouped as such in the ledger.

## Order of work

```
specify (done)
  → exercise offline UQ on synthetic_biofilm_e2e
  → integrate detectability outputs as the upstream producer of distributions
  → acquire Reference D measurements
  → target numerical sweep
  → scientific offline gate
  → fit supported response posteriors
  → only then may the online gate become eligible
```

A physical row may carry a distribution only once its pilot cleared detectability, basis and
pairing. Converting a non-detectable measurement into a very wide Gaussian is not
conservatism; it is inventing information, and such rows take
`distribution_status = NOT_IDENTIFIABLE` instead.

## Terminal state of this pass

```
OFFLINE_GATE_IMPLEMENTATION = SYNTHETICALLY_VALIDATED
OFFLINE_TARGET_VERDICT      = NOT_EVALUATED
ONLINE_GATE_LOGIC           = FAIL_CLOSED_AND_TESTED
ONLINE_FEEDBACK             = DISABLED
```

That is a successful pass, not an incomplete one. The machinery exists and is tested; no
target verdict has been rendered, because no measured configuration and no declared
threshold exist to render one from.
