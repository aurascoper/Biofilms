"""Tier S0: the synthetic gate, which decides verdicts without deciding science.

This is a SEPARATE TIER from the scientific offline gate in `feedback_gate.py`,
and the separation is structural rather than conventional: `decide()` cannot
return `OFFLINE_PASS` because that constant is not in its vocabulary. Synthetic
evidence may pass a synthetic gate. It may never pass the offline feedback gate,
whatever its effect size, and no configuration of this module makes it possible.

WHY THE INDETERMINATE VERDICTS ARE SPLIT. "Too uncertain to tell" is not one
outcome, because the four ways of being uncertain call for four different and
differently priced actions:

    INDETERMINATE_TRANSPORT    more histories, more replicates.       Reducible.
    INDETERMINATE_NUMERICS     finer mesh, more raytrace samples.
                               A FLOOR: histories cannot remove it — A0 showed
                               statistical error falling while resolution loss
                               did not.
    INDETERMINATE_CALIBRATION  better measurements. No amount of compute helps.
    INDETERMINATE_MODEL_FORM   a decision. Two admissible scenarios disagree
                               across the threshold, and averaging them under
                               arbitrary weights would manufacture a verdict
                               nobody is entitled to.

Collapsing them into one number is how a compute budget gets spent on a
calibration problem.

The verdict is a function of a DISTRIBUTION and a declared threshold, never of a
point estimate and an error bar.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# Tier S0 vocabulary. `OFFLINE_PASS` is deliberately absent.
EFFECT_BELOW_THRESHOLD = "EFFECT_BELOW_THRESHOLD"
PASS_SYNTHETIC_GATE = "PASS_SYNTHETIC_GATE"
INDETERMINATE_TRANSPORT = "INDETERMINATE_TRANSPORT"
INDETERMINATE_NUMERICS = "INDETERMINATE_NUMERICS"
INDETERMINATE_CALIBRATION = "INDETERMINATE_CALIBRATION"
INDETERMINATE_MODEL_FORM = "INDETERMINATE_MODEL_FORM"
EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT = (
    "EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT")
NOT_EVALUATED = "NOT_EVALUATED"

S0_VERDICTS = frozenset({
    EFFECT_BELOW_THRESHOLD, PASS_SYNTHETIC_GATE, INDETERMINATE_TRANSPORT,
    INDETERMINATE_NUMERICS, INDETERMINATE_CALIBRATION, INDETERMINATE_MODEL_FORM,
    EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT, NOT_EVALUATED,
})

# Verdicts this tier must never be able to produce, asserted rather than trusted.
FORBIDDEN_IN_S0 = frozenset({"OFFLINE_PASS", "PASS_OFFLINE_FEEDBACK",
                             "ONLINE_ENABLED"})

_SOURCE_TO_VERDICT = {
    "transport": INDETERMINATE_TRANSPORT,
    "numerics": INDETERMINATE_NUMERICS,
    "calibration": INDETERMINATE_CALIBRATION,
    "model_form": INDETERMINATE_MODEL_FORM,
}


@dataclass(frozen=True)
class VarianceBudget:
    """Variance by SOURCE, never summed into one number before the verdict.

    `numerics` is a floor and is labelled as one: mesh-resolution loss does not
    fall with histories, so a run that is numerics-dominated cannot be rescued
    by buying more of the thing that is cheap.
    """
    transport: float = 0.0       # OpenMC replicate spread
    numerics: float = 0.0        # mesh resolution loss + raytrace volume error
    calibration: float = 0.0     # outer correlated parameter draws
    model_form: float = 0.0      # spread across admissible discrete scenarios

    def as_dict(self) -> dict:
        return {"transport": self.transport, "numerics": self.numerics,
                "calibration": self.calibration, "model_form": self.model_form}

    @property
    def total(self) -> float:
        return sum(self.as_dict().values())

    @property
    def dominant(self) -> str:
        d = self.as_dict()
        return max(d, key=lambda k: d[k]) if self.total > 0 else "transport"


@dataclass(frozen=True)
class ThresholdPolicy:
    """Declared BEFORE the effect is seen. Tier S0 values are contract-test
    values and are explicitly not the eventual Reference D numbers, which stay
    unset until a deep-research pass declares them prospectively."""
    threshold_policy_id: str = "synthetic_gate_contract_v1"
    tier: str = "S0"
    target_calibration: bool = False
    # WHICH METRIC THE THRESHOLDS ARE DENOMINATED IN. A threshold is meaningless
    # without one: 0.10 declared for a relative L2 effect and 0.10 compared
    # against a SQUARED effect differ by a factor of ten, and nothing about the
    # number itself reveals the mistake. `decide()` refuses a mismatch rather
    # than trusting call sites to keep them aligned.
    metric_id: str = "relative_l2"
    # The effect magnitude that would matter, in `metric_id`'s units. At tier S0
    # this is a fixture constant used to exercise the logic, NOT a scientific
    # relevance claim.
    effect_threshold: float = 0.10
    # Below this the effect is real but not worth acting on.
    practical_importance_floor: float = 0.02
    # One-sided credibility for the lower bound.
    credibility: float = 0.99
    declared_before_results: bool = True

    def problems(self) -> list[str]:
        out = []
        if self.target_calibration:
            out.append("tier S0 may not carry target_calibration = true")
        if self.tier != "S0":
            out.append(f"tier {self.tier!r} is not S0")
        if not self.metric_id:
            out.append("a threshold with no metric_id is not interpretable")
        if not (0.0 < self.practical_importance_floor <= self.effect_threshold):
            out.append("practical_importance_floor must be in (0, effect_threshold]")
        if not 0.5 < self.credibility < 1.0:
            out.append("credibility must be in (0.5, 1)")
        if not self.declared_before_results:
            out.append("a threshold chosen after the effect was seen is not a threshold")
        return out


@dataclass(frozen=True)
class S0Verdict:
    verdict: str
    reason: str
    effect_median: float = float("nan")
    effect_low: float = float("nan")
    effect_high: float = float("nan")
    dominant_variance: str = ""
    variance: dict = field(default_factory=dict)
    n_draws: int = 0
    threshold_policy_id: str = ""
    tier: str = "S0"
    target_calibration: bool = False

    def as_dict(self) -> dict:
        """JSON-ready and NaN-free: a non-finite number is not JSON, and a
        verdict that only some parsers can read is not machine-readable."""
        def clean(x):
            return None if x is None or not math.isfinite(x) else float(x)
        return {
            "verdict": self.verdict, "reason": self.reason,
            "effect_median": clean(self.effect_median),
            "effect_low": clean(self.effect_low),
            "effect_high": clean(self.effect_high),
            "dominant_variance": self.dominant_variance,
            "variance": {k: clean(v) for k, v in self.variance.items()},
            "n_draws": int(self.n_draws),
            "threshold_policy_id": self.threshold_policy_id,
            "tier": self.tier,
            "target_calibration": bool(self.target_calibration),
        }


def scenario_envelope(scenario_effects: dict[str, np.ndarray]) -> dict:
    """Spread ACROSS admissible model-form scenarios, as an envelope.

    Deliberately not a weighted average. Two admissible discrete scenarios are
    not two samples from a distribution over scenarios, and assigning them equal
    probability invents information — a number that looks like evidence and is
    an assumption. The envelope says what the scenarios permit, and if they
    permit both sides of the threshold the answer is that nobody knows yet.
    """
    means = {k: float(np.mean(v)) for k, v in scenario_effects.items()}
    lo, hi = (min(means.values()), max(means.values())) if means else (
        float("nan"), float("nan"))
    return {"per_scenario_mean": means, "envelope_low": lo, "envelope_high": hi,
            "spread": hi - lo if means else float("nan")}


def decide(effect_draws, budget: VarianceBudget, policy: ThresholdPolicy, *,
           scenario_effects: dict[str, np.ndarray] | None = None,
           metric_id: str | None = None) -> S0Verdict:
    """The tier-S0 verdict. Exactly one outcome per input, by construction.

    `metric_id` names the metric the draws were computed with. Supplying it is
    optional but strongly advised: a threshold declared for one metric compared
    against draws from another is a silent factor-of-ten error, and it is the
    exact mistake that moving from a relative L2 effect to a debiased SQUARED
    effect invites.
    """
    if metric_id is not None and metric_id != policy.metric_id:
        return S0Verdict(
            NOT_EVALUATED,
            f"draws are {metric_id!r} but the policy declares thresholds in "
            f"{policy.metric_id!r}; a threshold is not transferable between "
            "metrics",
            threshold_policy_id=policy.threshold_policy_id, tier=policy.tier,
            target_calibration=policy.target_calibration)
    problems = policy.problems()
    if problems:
        return S0Verdict(NOT_EVALUATED, "; ".join(problems),
                         threshold_policy_id=policy.threshold_policy_id,
                         tier=policy.tier,
                         target_calibration=policy.target_calibration)

    draws = np.asarray(effect_draws, dtype=float)
    draws = draws[np.isfinite(draws)]
    common = dict(threshold_policy_id=policy.threshold_policy_id,
                  tier=policy.tier, variance=budget.as_dict(),
                  dominant_variance=budget.dominant,
                  target_calibration=policy.target_calibration)
    if draws.size == 0:
        return S0Verdict(NOT_EVALUATED, "no finite effect draws", **common)

    tail = 1.0 - policy.credibility
    low = float(np.quantile(draws, tail))
    high = float(np.quantile(draws, policy.credibility))
    median = float(np.median(draws))
    common.update(effect_median=median, effect_low=low, effect_high=high,
                  n_draws=int(draws.size))

    # MODEL FORM FIRST. If two admissible scenarios sit on opposite sides of the
    # threshold, no amount of sampling within either one settles it, and the
    # pooled distribution would hide the disagreement it is made of.
    if scenario_effects:
        env = scenario_envelope(scenario_effects)
        if env["envelope_low"] < policy.effect_threshold <= env["envelope_high"]:
            return S0Verdict(
                INDETERMINATE_MODEL_FORM,
                f"admissible scenarios straddle the threshold "
                f"({env['envelope_low']:.4g} to {env['envelope_high']:.4g} vs "
                f"{policy.effect_threshold:.4g}); they are not samples from a "
                "distribution and averaging them would invent a weight",
                **{**common, "dominant_variance": "model_form"})

    # Resolved above the threshold: the only passing outcome this tier has.
    if low > policy.effect_threshold:
        return S0Verdict(
            PASS_SYNTHETIC_GATE,
            f"the {policy.credibility:.0%} lower bound ({low:.4g}) exceeds the "
            f"declared threshold ({policy.effect_threshold:.4g})", **common)

    # Entirely below the threshold. Two different meanings, kept apart.
    if high < policy.effect_threshold:
        if high < policy.practical_importance_floor:
            return S0Verdict(
                EFFECT_BELOW_THRESHOLD,
                f"the whole distribution ({high:.4g} at the upper bound) sits "
                f"below the practical floor ({policy.practical_importance_floor:.4g})",
                **common)
        return S0Verdict(
            EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT,
            f"resolved above the practical floor "
            f"({policy.practical_importance_floor:.4g}) but entirely below the "
            f"threshold ({policy.effect_threshold:.4g}): real, and not worth "
            "acting on", **common)

    # Straddles the threshold. WHICH source of spread dominates decides what to
    # buy, so the verdict names it instead of saying "uncertain".
    return S0Verdict(
        _SOURCE_TO_VERDICT[budget.dominant],
        f"the distribution spans the threshold ({low:.4g} to {high:.4g} vs "
        f"{policy.effect_threshold:.4g}) and {budget.dominant} variance "
        f"dominates ({budget.as_dict()[budget.dominant]:.4g} of "
        f"{budget.total:.4g})"
        + (" — a resolution floor that more histories cannot remove"
           if budget.dominant == "numerics" else ""),
        **common)


assert not (S0_VERDICTS & FORBIDDEN_IN_S0), \
    "tier S0 must not be able to express an offline-feedback pass"
