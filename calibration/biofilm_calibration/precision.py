"""How many biological replicates, and whether the measurement is possible.

Two different questions, and conflating them wastes a campaign.

DETECTABILITY comes first. If the biofilm mass is small compared with the
uncertainty of the blank that must be subtracted from it, no number of
replicates rescues the measurement — the substrate geometry or the amount of
biomass has to change. Replicates shrink the standard error of a mean; they do
nothing about a signal buried under a subtraction.

REPLICATE COUNT comes second, and is derived rather than declared. `n = 3` is a
convention, not a calculation. The required count follows from the target
uncertainty and the observed variance components, and the dominant component is
reported alongside it because that is what says where to spend effort: more
batches, more fields per batch, or a better balance.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


class PrecisionError(ValueError):
    """Raised when a precision question is asked of unusable inputs."""


@dataclass(frozen=True)
class Detectability:
    detectable: bool
    signal_to_blank_noise: float
    reason: str


def detectability(biofilm_mass_g: float, blank_uncertainty_g: float, *,
                  minimum_ratio: float = 10.0) -> Detectability:
    """Is the biofilm mass measurable above the blank's own uncertainty?

    `minimum_ratio` is a DECLARED engineering target, not a law: at 10 the
    blank contributes about 10% of the biofilm mass in the worst case. Raise it
    for a tighter density, but note that the fix for a failure is never more
    replicates.
    """
    if biofilm_mass_g is None or biofilm_mass_g <= 0:
        raise PrecisionError(f"biofilm mass must be positive, got {biofilm_mass_g}")
    if blank_uncertainty_g is None or blank_uncertainty_g < 0:
        raise PrecisionError(
            f"blank uncertainty must be non-negative, got {blank_uncertainty_g}")
    if blank_uncertainty_g == 0:
        return Detectability(True, float("inf"),
                             "blank uncertainty declared as zero — check that "
                             "it was measured rather than assumed")
    ratio = float(biofilm_mass_g) / float(blank_uncertainty_g)
    if ratio >= minimum_ratio:
        return Detectability(True, ratio,
                             f"biofilm mass is {ratio:.1f}x the blank "
                             "uncertainty")
    return Detectability(
        False, ratio,
        f"biofilm mass is only {ratio:.1f}x the blank uncertainty (target "
        f"{minimum_ratio:.0f}x). This is a SUBSTRATE GEOMETRY or BIOMASS "
        "problem, not a replicate-count problem: replicates shrink the standard "
        "error of a mean, they do not recover a signal buried under a "
        "subtraction. Grow more biomass per coupon, or use a substrate whose "
        "blank is lighter or more reproducible.")


@dataclass(frozen=True)
class ReplicatePlan:
    required_batches: int
    achieved_rel_uncertainty: float
    dominant_component: str
    components: dict
    note: str


def required_replicates(*, target_rel_uncertainty: float, mean_value: float,
                        var_between_batch: float, var_within_batch: float,
                        n_within: int = 1, var_blank: float = 0.0,
                        max_batches: int = 200) -> ReplicatePlan:
    """Biological replicates needed to reach a target relative uncertainty.

    The standard error of the grand mean over `n` batches, each with `n_within`
    technical measurements, is

        SE^2 = var_between/n + var_within/(n * n_within) + var_blank/n

    so every component falls as 1/n except that `var_within` also falls with
    the number of measurements per batch. That asymmetry is the reason the
    dominant component is reported: if between-batch variance dominates, more
    fields per coupon buy nothing.
    """
    if not 0 < target_rel_uncertainty < 1:
        raise PrecisionError(
            f"target relative uncertainty must be in (0, 1), got "
            f"{target_rel_uncertainty}")
    if mean_value is None or mean_value <= 0:
        raise PrecisionError(f"mean value must be positive, got {mean_value}")
    for name, v in (("var_between_batch", var_between_batch),
                    ("var_within_batch", var_within_batch),
                    ("var_blank", var_blank)):
        if v is None or v < 0:
            raise PrecisionError(f"{name} must be non-negative, got {v}")
    if n_within < 1:
        raise PrecisionError(f"n_within must be >= 1, got {n_within}")

    target_se = float(target_rel_uncertainty) * float(mean_value)
    per_batch = (float(var_between_batch)
                 + float(var_within_batch) / n_within
                 + float(var_blank))
    if per_batch == 0:
        return ReplicatePlan(1, 0.0, "none", {}, "no variance declared")

    n = math.ceil(per_batch / target_se ** 2)
    n = max(2, n)          # a variance estimate needs at least two batches
    capped = n > max_batches
    n_final = min(n, max_batches)
    achieved = math.sqrt(per_batch / n_final) / float(mean_value)

    components = {
        "between_batch": float(var_between_batch),
        "within_batch_effective": float(var_within_batch) / n_within,
        "blank": float(var_blank),
    }
    dominant = max(components, key=components.get)

    note = ""
    if dominant == "between_batch":
        note = ("between-batch variance dominates: more imaging fields per "
                "coupon buy nothing, only more independent cultures do")
    elif dominant == "within_batch_effective":
        note = (f"within-batch variance dominates at {n_within} measurement(s) "
                "per batch: more fields per coupon are cheaper than more "
                "cultures")
    elif dominant == "blank":
        note = ("the blank dominates: a lighter or more reproducible substrate "
                "is worth more than either kind of replicate")
    if capped:
        note += (f" | CAPPED at {max_batches} batches; the target is not "
                 "reachable by replication at these variances")
    return ReplicatePlan(n_final, achieved, dominant, components, note)


def plan_problems(plan: ReplicatePlan, target_rel_uncertainty: float) -> list[str]:
    """Gate-facing: does this plan actually meet its target?"""
    out = []
    if plan.achieved_rel_uncertainty > target_rel_uncertainty:
        out.append(
            f"achieved relative uncertainty {plan.achieved_rel_uncertainty:.3%} "
            f"exceeds the target {target_rel_uncertainty:.3%} at "
            f"{plan.required_batches} batches — {plan.note}")
    return out
