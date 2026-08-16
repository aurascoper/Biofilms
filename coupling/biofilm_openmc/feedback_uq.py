"""Paired counterfactual effects, and where their variance comes from.

The offline gate asks what changing a material state does to the transport
field. That is a PAIRED question, and it must be computed as one:

    dQ = Q(M1) - Q(M0)

evaluated on the SAME immutable snapshot, not as two independent results with
error bars to be compared. Subtracting two separately-reported means and adding
their uncertainties in quadrature throws away exactly the correlation that a
paired design exists to exploit.

VARIANCE HAS TWO SOURCES and the law of total variance separates them:

    Var(dQ) = E_theta[ Var_xi(dQ | theta) ]  +  Var_theta[ E_xi(dQ | theta) ]
              \\_______ transport _______/     \\____ calibration/biology ____/

The first term is what more histories and more seeds reduce. The second is what
only better measurements reduce. A budget that reports one number cannot say
which of those to spend on, which is the practical reason the decomposition is
worth the bookkeeping.

Mesh bins are NOT independent — one history deposits in several — so per-label
uncertainty comes from independent replicate runs, exactly as `lineage.py`
already insists. Nothing here sums per-bin variances as if they were.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Every metric here is a RATIO or a paired difference of quantities carrying the
# same units, so a common source-rate multiplier cancels. That is why the
# offline material screen can run per source particle, BEFORE the source
# activity is known — the same reasoning `drivers.compare_fields` already
# records. Absolute Gy/s and time-integrated dose get no such cancellation.
SOURCE_RATE_INVARIANT = ("relative_effect", "mass_weighted_mean_dose",
                         "absorbed_energy_rate")


def absorbed_energy_rate(dose_field, mass_kg) -> float:
    """Q_E = sum_v D_v m_v — the whole-biomass transport quantity.

    Dose rate times mass is energy per time, so this is the most fundamental
    thing a material change can move, and it is the same product the lineage
    attribution conserves.
    """
    d = np.asarray(dose_field, dtype=float)
    m = np.asarray(mass_kg, dtype=float)
    return float(np.nansum(d * m))


def mass_weighted_mean_dose(dose_field, mass_kg, mask=None) -> float:
    d = np.asarray(dose_field, dtype=float)
    m = np.asarray(mass_kg, dtype=float)
    if mask is not None:
        sel = np.asarray(mask, dtype=bool)
        d, m = d[sel], m[sel]
    total = float(np.nansum(m))
    return float("nan") if total == 0 else float(np.nansum(d * m) / total)


def upper_tail_dose(dose_field, mask=None, percentile: float = 95.0) -> float:
    """A predeclared upper percentile. A localized effect can be real and
    invisible in a whole-system average, so the hierarchy carries one."""
    d = np.asarray(dose_field, dtype=float)
    if mask is not None:
        d = d[np.asarray(mask, dtype=bool)]
    d = d[np.isfinite(d)]
    return float("nan") if d.size == 0 else float(np.percentile(d, percentile))


@dataclass(frozen=True)
class PairedEffect:
    """One outer parameter point, with its transport replicates."""
    qoi_id: str
    baseline: np.ndarray        # per-seed Q(M0)
    feedback: np.ndarray        # per-seed Q(M1)

    @property
    def per_seed_effect(self) -> np.ndarray:
        """Paired per seed where the seeds match, so a common random stream can
        cancel. Whether that reduces variance is an EMPIRICAL question — a
        material change makes histories diverge — so it is measured by
        `common_random_number_benefit`, never assumed."""
        n = min(len(self.baseline), len(self.feedback))
        return np.asarray(self.feedback[:n]) - np.asarray(self.baseline[:n])

    @property
    def mean_effect(self) -> float:
        e = self.per_seed_effect
        return float(np.mean(e)) if e.size else float("nan")

    @property
    def transport_variance(self) -> float:
        """Var_xi(dQ | theta): what more seeds and more histories reduce."""
        e = self.per_seed_effect
        return float(np.var(e, ddof=1)) if e.size > 1 else 0.0

    @property
    def relative_effect(self) -> float:
        b = float(np.mean(self.baseline)) if len(self.baseline) else float("nan")
        return float("nan") if not b else self.mean_effect / b


@dataclass(frozen=True)
class VarianceDecomposition:
    transport: float            # E_theta[Var_xi]
    parameter: float            # Var_theta[E_xi]
    total: float
    n_outer: int
    n_seeds_min: int

    @property
    def transport_share(self) -> float:
        return float("nan") if self.total == 0 else self.transport / self.total

    def render(self) -> str:
        return (f"Var(dQ) = {self.total:.6g}\n"
                f"  transport (more seeds/histories help):  {self.transport:.6g}"
                f"  [{self.transport_share:.1%}]\n"
                f"  parameter (only measurements help):     {self.parameter:.6g}"
                f"  [{1 - self.transport_share:.1%}]\n"
                f"  {self.n_outer} outer points, >= {self.n_seeds_min} seeds each")


def decompose_variance(effects: list[PairedEffect]) -> VarianceDecomposition:
    """Law of total variance over outer parameter points.

    Reported rather than collapsed, because the two terms are reduced by
    completely different and differently priced actions.
    """
    if not effects:
        return VarianceDecomposition(float("nan"), float("nan"), float("nan"), 0, 0)
    means = np.array([e.mean_effect for e in effects], dtype=float)
    within = np.array([e.transport_variance for e in effects], dtype=float)
    finite = np.isfinite(means)
    transport = float(np.mean(within[finite])) if finite.any() else float("nan")
    parameter = (float(np.var(means[finite], ddof=1))
                 if finite.sum() > 1 else 0.0)
    return VarianceDecomposition(
        transport, parameter, transport + parameter, int(finite.sum()),
        min((e.per_seed_effect.size for e in effects), default=0))


def effect_draws(effects: list[PairedEffect]) -> np.ndarray:
    """The joint distribution the gate decides on: every outer point's mean
    effect, carrying parameter uncertainty, with transport noise already
    averaged within each point."""
    return np.array([e.mean_effect for e in effects], dtype=float)


def seed_sufficiency(n_seeds: int) -> dict:
    """How well a variance estimated from `n_seeds` replicates is itself known.

    For an approximately normal replicate distribution the relative error of an
    estimated standard deviation is about 1/sqrt(2(r-1)) — roughly 35% at 5
    seeds and 16% at 20. Five seeds spot a large effect; they are weak evidence
    for a tight budget, and the gate should say so rather than let a number
    stand in for its own precision.
    """
    if n_seeds < 2:
        return {"n_seeds": n_seeds, "sd_relative_error": float("inf"),
                "adequate_for_budget": False,
                "note": "a variance estimate needs at least two replicates"}
    rel = 1.0 / np.sqrt(2.0 * (n_seeds - 1))
    return {
        "n_seeds": n_seeds,
        "sd_relative_error": float(rel),
        # 20 seeds -> ~16%, the recommended floor for characterising a
        # production condition rather than merely screening one.
        "adequate_for_budget": n_seeds >= 20,
        "note": ("adequate for screening, not for a tight uncertainty budget"
                 if n_seeds < 20 else "adequate for a production budget"),
    }


def common_random_number_benefit(paired: PairedEffect,
                                 independent: PairedEffect) -> dict:
    """Did matched seeds actually reduce Var(dQ)?

    Adopted only if measured. Matched streams cancel common noise when the two
    material states transport similarly, but a large material change makes
    histories diverge and the pairing can stop helping. Assuming the benefit
    would understate the uncertainty of exactly the effects that matter most.
    """
    v_paired = paired.transport_variance
    v_indep = independent.transport_variance
    ratio = float("nan") if v_indep == 0 else v_paired / v_indep
    return {
        "variance_paired": v_paired,
        "variance_independent": v_indep,
        "variance_ratio": ratio,
        "adopt_common_random_numbers": bool(np.isfinite(ratio) and ratio < 1.0),
        "note": ("matched seeds reduced the paired variance"
                 if np.isfinite(ratio) and ratio < 1.0 else
                 "matched seeds did not help; use independent streams"),
    }
