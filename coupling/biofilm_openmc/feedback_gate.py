"""The feedback gates, as a fail-closed decision engine.

TWO GATES, TWO QUESTIONS.

    OFFLINE  Does changing a biologically controlled material state produce a
             transport effect distinguishable from Monte Carlo noise, numerical
             discretization and calibration uncertainty? A counterfactual on
             immutable snapshots. Nothing is fed back into biology.
    ONLINE   May the real two-way loop change state? An authorization
             condition, not a measurement.

THE DECISION RULE is a statement about the whole distribution, not a
comparison of error bars:

    q_0.01( s * dQ ) > delta        equivalently   P(s*dQ > delta) >= 0.99

`s` is the predeclared direction, because a two-sided test passes on an effect
of the wrong sign. `delta` is the predeclared scientific effect threshold and
HAS NO DEFAULT — until someone declares what magnitude of feedback would
matter, the honest verdict is NOT_EVALUATED. Quadrature over independent error
bars is not offered: several of the project's uncertainties are correlated, and
`uncertainty.propagate` says so of its own inputs.

EVERYTHING FAILS CLOSED. A missing threshold, an empty distribution, a NaN, an
unresolved model branch, or a state outside the validated envelope all return a
blocking verdict rather than a permissive default. This module consumes
versioned distributions and thresholds; it does not decide what a plausible
melanin range is, and it must never learn to.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from physical_contract import (DEFAULT_EXCEEDANCE_PROBABILITY,
                               DEFAULT_NUMERICAL_BUDGET_FRACTION,
                               DEFAULT_SURROGATE_BUDGET_FRACTION,
                               DEFAULT_VOLUME_BUDGET_FRACTION,
                               EFFECT_DIRECTIONS, FEEDBACK_GATE_VERDICTS)

NOT_EVALUATED = "NOT_EVALUATED"
BLOCKED_PHYSICAL = "BLOCKED_ON_PHYSICAL_CALIBRATION"
BLOCKED_NUMERICAL = "BLOCKED_ON_NUMERICAL_RESOLUTION"
BLOCKED_BIOLOGICAL = "BLOCKED_ON_BIOLOGICAL_CALIBRATION"
BELOW_THRESHOLD = "OFFLINE_EFFECT_BELOW_THRESHOLD"
OFFLINE_UNCERTAIN = "OFFLINE_UNCERTAINTY_TOO_LARGE"
OFFLINE_PASS = "OFFLINE_PASS"
ONLINE_OUT_OF_DOMAIN = "ONLINE_OUT_OF_DOMAIN"
ONLINE_UNCERTAIN = "ONLINE_UNCERTAINTY_TOO_LARGE"
ONLINE_ENABLED = "ONLINE_ENABLED"
UNSUPPORTED = "UNSUPPORTED_BY_CURRENT_MODEL"

# The production mass denominator. The full-bin method overstates A0's cylinder
# by exactly 4/pi because it weighs the circumscribing cube; it is a known
# negative control, never one member of an equally weighted model ensemble.
PRODUCTION_MASS_DENOMINATOR = "exact_csg_raytrace"


@dataclass(frozen=True)
class EffectPolicy:
    """What counts as an effect worth acting on, declared BEFORE results.

    `effect_threshold` is deliberately not defaulted. A gate that supplies its
    own relevance threshold has decided the science it was meant to test.
    """
    qoi_id: str
    effect_direction: str
    effect_threshold: float | None = None
    exceedance_probability: float = DEFAULT_EXCEEDANCE_PROBABILITY
    numerical_budget_fraction: float = DEFAULT_NUMERICAL_BUDGET_FRACTION
    volume_budget_fraction: float = DEFAULT_VOLUME_BUDGET_FRACTION
    surrogate_budget_fraction: float = DEFAULT_SURROGATE_BUDGET_FRACTION
    policy_id: str = ""
    declared_before_results: bool = False

    def problems(self) -> list[str]:
        out = []
        if self.effect_threshold is None:
            out.append(
                "effect_threshold is undeclared — no uncertainty calculation "
                "can choose what magnitude of feedback matters, and a gate that "
                "supplies its own relevance threshold has decided the science "
                "it was meant to test")
        elif self.effect_threshold <= 0:
            out.append(f"effect_threshold {self.effect_threshold} is not positive")
        if self.effect_direction not in EFFECT_DIRECTIONS:
            out.append(f"effect_direction {self.effect_direction!r} not in "
                       f"{sorted(EFFECT_DIRECTIONS)} — a two-sided test passes "
                       "on an effect of the wrong sign")
        if not 0.5 < self.exceedance_probability < 1.0:
            out.append("exceedance_probability must be in (0.5, 1)")
        if not self.policy_id:
            out.append("policy_id is empty — a threshold with no version cannot "
                       "be revised prospectively, only silently")
        if not self.declared_before_results:
            out.append("declared_before_results is false — a threshold chosen "
                       "after the effect was seen is not a threshold")
        return out

    @property
    def sign(self) -> float:
        return 1.0 if self.effect_direction == "increase" else -1.0


@dataclass(frozen=True)
class UncertaintyBudget:
    """The controllable half of the uncertainty, which must be small relative
    to the decision margin or a pass is driven by barely resolved numerics."""
    transport_numerical_u99: float = float("nan")
    material_volume_u99: float = float("nan")
    surrogate_error: float = 0.0
    mass_denominator_method: str = ""
    unresolved_model_branches: tuple = ()

    def problems(self, policy: EffectPolicy) -> list[str]:
        out = []
        delta = policy.effect_threshold
        if delta is None:
            return ["no threshold to size the numerical budget against"]
        checks = (
            ("transport+numerical", self.transport_numerical_u99,
             policy.numerical_budget_fraction),
            ("material-volume", self.material_volume_u99,
             policy.volume_budget_fraction),
            ("surrogate", self.surrogate_error, policy.surrogate_budget_fraction),
        )
        for name, value, fraction in checks:
            if value is None or (isinstance(value, float) and math.isnan(value)):
                out.append(f"{name} uncertainty is not characterised")
            elif value > fraction * delta:
                out.append(
                    f"{name} uncertainty {value:.4g} exceeds {fraction:.0%} of "
                    f"the effect threshold {delta:.4g} — a pass here would be "
                    "driven by how well the numerics were resolved")
        if self.mass_denominator_method != PRODUCTION_MASS_DENOMINATOR:
            out.append(
                f"mass_denominator_method is {self.mass_denominator_method!r}, "
                f"not {PRODUCTION_MASS_DENOMINATOR!r} — the full-bin denominator "
                "weighs the circumscribing cube and is a known negative "
                "control, not an alternative")
        if self.unresolved_model_branches:
            out.append(
                "unresolved model branches: "
                f"{', '.join(self.unresolved_model_branches)} — averaging over "
                "an undeclared choice is not a result")
        return out


@dataclass(frozen=True)
class OfflineVerdict:
    verdict: str
    reasons: tuple = ()
    effect_mean: float = float("nan")
    effect_q001: float = float("nan")
    probability_exceeds: float = float("nan")
    n_draws: int = 0

    @property
    def passed(self) -> bool:
        return self.verdict == OFFLINE_PASS


def directed_exceedance(effect_draws, policy: EffectPolicy) -> tuple:
    """(q_0.01 of the directed effect, P(directed effect > delta)).

    The full joint distribution, not a Gaussian half-width: the model is
    nonlinear, several inputs are bounded or compositional, and the output
    need not be symmetric.
    """
    draws = np.asarray(effect_draws, dtype=float)
    finite = draws[np.isfinite(draws)]
    if finite.size == 0:
        return float("nan"), float("nan"), 0
    directed = policy.sign * finite
    q = float(np.quantile(directed, 1.0 - policy.exceedance_probability))
    p = float(np.mean(directed > policy.effect_threshold)) \
        if policy.effect_threshold is not None else float("nan")
    return q, p, int(finite.size)


def evaluate_offline(effect_draws, policy: EffectPolicy,
                     budget: UncertaintyBudget, *,
                     physical_calibration_ready: bool = False,
                     supported_by_model: bool = True) -> OfflineVerdict:
    """The counterfactual verdict. Fails closed at every step."""
    if not supported_by_model:
        return OfflineVerdict(UNSUPPORTED, (
            "the response this gate would license is not representable by the "
            "current model — no quantity of data changes that",))

    policy_problems = policy.problems()
    if policy_problems:
        return OfflineVerdict(NOT_EVALUATED, tuple(policy_problems))

    if not physical_calibration_ready:
        return OfflineVerdict(BLOCKED_PHYSICAL, (
            "Reference D physical calibration is not READY — an effect measured "
            "on uncalibrated material is an effect on a material nobody has",))

    budget_problems = budget.problems(policy)
    if budget_problems:
        return OfflineVerdict(BLOCKED_NUMERICAL, tuple(budget_problems))

    q, p, n = directed_exceedance(effect_draws, policy)
    if n == 0:
        return OfflineVerdict(NOT_EVALUATED,
                              ("no finite effect draws to decide on",))
    mean = float(np.mean(np.asarray(effect_draws, dtype=float)[
        np.isfinite(np.asarray(effect_draws, dtype=float))]))

    if q > policy.effect_threshold:
        return OfflineVerdict(OFFLINE_PASS, (
            f"q{1 - policy.exceedance_probability:.2f} of the directed effect "
            f"({q:.4g}) exceeds the declared threshold "
            f"({policy.effect_threshold:.4g})",), mean, q, p, n)

    # Distinguish "too small to matter" from "too uncertain to tell". Both
    # block, and conflating them would hide which one more histories can fix.
    spread = float(np.std(np.asarray(effect_draws, dtype=float)[
        np.isfinite(np.asarray(effect_draws, dtype=float))], ddof=1)) if n > 1 else 0.0
    directed_mean = policy.sign * mean
    if directed_mean > policy.effect_threshold:
        return OfflineVerdict(OFFLINE_UNCERTAIN, (
            f"the directed effect mean ({directed_mean:.4g}) exceeds the "
            f"threshold but its lower bound ({q:.4g}) does not — the "
            f"distribution is too wide (sd {spread:.4g}) to decide",),
            mean, q, p, n)
    return OfflineVerdict(BELOW_THRESHOLD, (
        f"the directed effect ({directed_mean:.4g}) does not reach the declared "
        f"threshold ({policy.effect_threshold:.4g}); more histories would "
        "sharpen it, not enlarge it",), mean, q, p, n)


@dataclass(frozen=True)
class OnlinePrerequisites:
    """Everything that must independently be true before feedback may change
    state. Every field defaults to the blocking value."""
    offline: OfflineVerdict | None = None
    physical_calibration_ready: bool = False
    numerical_resolution_ready: bool = False
    biological_posterior_ready: bool = False
    seconds_per_mcs_ready: bool = False
    state_within_validity_domain: bool = False
    response_pathway: str = ""
    unsupported_pathways: tuple = ("growth_survival", "membrane_response")
    extras: dict = field(default_factory=dict)


def evaluate_online(prereq: OnlinePrerequisites, policy: EffectPolicy,
                    budget: UncertaintyBudget,
                    current_effect_draws=None) -> OfflineVerdict:
    """Authorization for state-changing feedback. Returns a verdict from the
    same vocabulary so one state machine covers both gates.

    Ordering is deliberate: the cheapest and most fundamental refusals come
    first, so a blocked gate names the earliest thing that is wrong rather than
    the last thing checked.
    """
    if prereq.response_pathway in prereq.unsupported_pathways:
        return OfflineVerdict(UNSUPPORTED, (
            f"response pathway {prereq.response_pathway!r} is unsupported by the "
            "current model; authorizing it would licence a response the model "
            "cannot represent",))

    if prereq.offline is None or not prereq.offline.passed:
        got = prereq.offline.verdict if prereq.offline else NOT_EVALUATED
        return OfflineVerdict(got if got != OFFLINE_PASS else NOT_EVALUATED, (
            "the offline gate has not passed — online feedback may not be "
            f"enabled on the strength of {got}",))

    if not (prereq.physical_calibration_ready and prereq.numerical_resolution_ready):
        return OfflineVerdict(BLOCKED_PHYSICAL, (
            "physical calibration or the production resolution decision is not "
            "READY",))
    if not (prereq.biological_posterior_ready and prereq.seconds_per_mcs_ready):
        return OfflineVerdict(BLOCKED_BIOLOGICAL, (
            "the biological response posterior or seconds_per_mcs is not "
            "calibrated — a dose cannot be turned into a biological signal, "
            "nor accrued over time, without them",))
    if not prereq.state_within_validity_domain:
        return OfflineVerdict(ONLINE_OUT_OF_DOMAIN, (
            "the current simulation state is outside the envelope the offline "
            "gate validated — offline evidence does not license extrapolation",))

    budget_problems = budget.problems(policy)
    if budget_problems:
        return OfflineVerdict(BLOCKED_NUMERICAL, tuple(budget_problems))

    if current_effect_draws is not None:
        live = evaluate_offline(current_effect_draws, policy, budget,
                                physical_calibration_ready=True)
        if not live.passed:
            return OfflineVerdict(ONLINE_UNCERTAIN, (
                "the effect at the CURRENT state no longer clears the "
                f"threshold: {live.verdict}",) + live.reasons,
                live.effect_mean, live.effect_q001, live.probability_exceeds,
                live.n_draws)

    return OfflineVerdict(ONLINE_ENABLED, (
        "offline pass, calibration ready, state inside the validated envelope, "
        "and the effect still exceeds its declared threshold",))


assert set(FEEDBACK_GATE_VERDICTS) >= {
    NOT_EVALUATED, BLOCKED_PHYSICAL, BLOCKED_NUMERICAL, BLOCKED_BIOLOGICAL,
    BELOW_THRESHOLD, OFFLINE_UNCERTAIN, OFFLINE_PASS, ONLINE_OUT_OF_DOMAIN,
    ONLINE_UNCERTAIN, ONLINE_ENABLED, UNSUPPORTED}, \
    "a verdict exists here that the shared vocabulary does not know"
