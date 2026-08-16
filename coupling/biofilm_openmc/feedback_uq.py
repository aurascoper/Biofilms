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

A NORM OF A NOISY DIFFERENCE IS BIASED, AND ITS SCATTER WILL NOT SAY SO. This is
the subtlest thing in this module. `relative_l2_effect` is unsigned, so with zero
true effect it returns |noise| rather than 0, and it does so *stably* — the bias
is systematic, so replicate variance stays small and reassuring while the
estimate sits well above zero. The pilot measured a transport standard deviation
of 0.0023 against a bias of 0.0545, a factor of 24, and reported only the former.

`debiased_squared_effect` removes it by construction rather than by subtraction:
across independent replicates the cross terms have no noise contribution, so the
estimator is centred on the true squared effect and is negative about half the
time when that effect is zero. It is what a gate should consume;
`relative_l2_effect` is what a diagnostic should.
"""

from __future__ import annotations

import math
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


def relative_l2_effect(d0, d1, mass, mask=None) -> float:
    """Mass-weighted relative L2 field effect. THE RAW, BIASED DIAGNOSTIC.

        E = sqrt( sum_b m_b (D1 - D0)^2 ) / sqrt( sum_b m_b D0^2 )

    Not a raw maximum relative difference: near-zero bins would make that a
    measure of the dose floor rather than of the material change.

    IT IS AN UNSIGNED NORM, so it CANNOT return zero when the true effect is
    zero — residual Monte Carlo noise enters as positive bias that never
    cancels, and it does so stably, so replicate scatter never reveals it. The
    pilot measured that floor at 0.0545 (finest mesh, 1e5 histories per run),
    which is larger than four of six material levers of interest.

    Kept because it remains the right diagnostic for three jobs: checking the
    expected N^-1/2 history scaling, measuring what common random numbers buy,
    and catching gross estimator regressions. For a GATE, use
    `debiased_squared_effect` instead.
    """
    a = np.asarray(d0, dtype=float)
    b = np.asarray(d1, dtype=float)
    m = np.asarray(mass, dtype=float)
    if mask is not None:
        sel = np.asarray(mask, dtype=bool)
        a, b, m = a[sel], b[sel], m[sel]
    denom = math.sqrt(float(np.sum(m * a * a)))
    if denom == 0.0:
        return float("nan")
    return math.sqrt(float(np.sum(m * (b - a) ** 2))) / denom


def cross_replicate_inner(X, w) -> float:
    """The U-statistic kernel: mean of x_r^T W x_s over ORDERED PAIRS r != s.

        (1 / (R(R-1))) * sum_{r != s} x_r^T W x_s

    Because replicates are independent, E[x_r^T W x_s] = E[x_r]^T W E[x_s] for
    r != s, so this estimates ||E[x]||^2_W with NO noise term. The diagonal
    r == s is exactly where the noise would enter — E[x_r^T W x_r] =
    ||E[x]||^2_W + E[||noise||^2_W] — which is why it is excluded rather than
    subtracted afterwards.

    THE DIAGONAL IS ZEROED, NOT SUBTRACTED. `G.sum() - trace(G)` is
    algebraically identical but forms and then cancels the large diagonal sum,
    and near the null the off-diagonal total is tiny by comparison — precisely
    the regime this estimator exists for. Zeroing keeps the cancellation out of
    the arithmetic.

    `X` is (R, B): one row per replicate, restricted to the region of interest.
    `w` is the (B,) weight vector, i.e. the diagonal of W.
    """
    x = np.asarray(X, dtype=float)
    weights = np.asarray(w, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"X must be (R, B), got shape {x.shape}")
    r = x.shape[0]
    if r < 2:
        return float("nan")
    gram = (x * weights) @ x.T
    np.fill_diagonal(gram, 0.0)
    return float(gram.sum() / (r * (r - 1)))


@dataclass(frozen=True)
class DebiasedEffect:
    """A squared effect that is zero in expectation when nothing changed."""
    s_hat: float                # unbiased ||E[d]||^2_W, MAY BE NEGATIVE
    denominator: float          # unbiased ||E[D0]||^2_W
    e_squared: float            # s_hat / denominator, the gate statistic
    n_replicates: int
    jackknife_var: float        # delete-one-replicate variance of e_squared
    metric_id: str = "debiased_relative_l2_squared"

    @property
    def signed_root(self) -> float:
        """FOR HUMAN DISPLAY ONLY. Never feed this to a gate: taking a root
        near zero reintroduces exactly the positive bias the estimator removes,
        and the sign carries no directional meaning — it is an artifact of an
        unbiased estimate of a non-negative quantity landing below zero."""
        e = self.e_squared
        if not math.isfinite(e):
            return float("nan")
        return math.copysign(math.sqrt(abs(e)), e)

    def as_dict(self) -> dict:
        def clean(x):
            return None if x is None or not math.isfinite(x) else float(x)
        return {"s_hat": clean(self.s_hat), "denominator": clean(self.denominator),
                "e_squared": clean(self.e_squared),
                "n_replicates": int(self.n_replicates),
                "jackknife_var": clean(self.jackknife_var),
                "metric_id": self.metric_id}


def _ratio(baseline_rows, difference_rows, weights) -> float:
    num = cross_replicate_inner(difference_rows, weights)
    den = cross_replicate_inner(baseline_rows, weights)
    if not math.isfinite(den) or den == 0.0:
        return float("nan")
    return num / den


def debiased_squared_effect(baseline, feedback, mass, mask=None) -> DebiasedEffect:
    """The gate statistic: an unbiased estimate of the SQUARED relative effect.

    `baseline` and `feedback` are sequences of R replicate fields. Within a
    replicate the two states may share a random stream (common random numbers,
    which reduce Var(d_r) and are therefore welcome); ACROSS replicates they
    must be independent, which is what makes the r != s terms unbiased.

    Returned squared, and compared against delta^2, because taking a square
    root near zero reintroduces the positive bias this exists to remove.

    THE DENOMINATOR IS DEBIASED TOO. E[||D0_r||^2_W] carries the baseline's own
    noise power, so a naive ||D0||^2 denominator is biased HIGH and would drag
    every ratio low. The same cross-replicate construction removes it.
    """
    b_rows = [np.asarray(x, dtype=float).ravel() for x in baseline]
    f_rows = [np.asarray(x, dtype=float).ravel() for x in feedback]
    r = min(len(b_rows), len(f_rows))
    if r < 3:
        # 2 replicates admit the estimator but not its jackknife, and a variance
        # nobody can estimate is not a usable gate input.
        return DebiasedEffect(float("nan"), float("nan"), float("nan"), r,
                              float("nan"))

    m = np.asarray(mass, dtype=float).ravel()
    if mask is None:
        sel = slice(None)
    else:
        sel = np.asarray(mask, dtype=bool).ravel()
    weights = m[sel]

    d0 = np.array([row[sel] for row in b_rows[:r]])
    d1 = np.array([row[sel] for row in f_rows[:r]])
    diff = d1 - d0

    s_hat = cross_replicate_inner(diff, weights)
    denom = cross_replicate_inner(d0, weights)
    e2 = _ratio(d0, diff, weights)

    # Delete-one-replicate jackknife on the RATIO, since the ratio is what the
    # gate consumes and its variance is not the numerator's.
    keep = [i for i in range(r)]
    loo = np.array([_ratio(np.delete(d0, i, axis=0), np.delete(diff, i, axis=0),
                           weights) for i in keep], dtype=float)
    finite = loo[np.isfinite(loo)]
    jack = (float((r - 1) / r * np.sum((finite - finite.mean()) ** 2))
            if finite.size > 1 else float("nan"))

    return DebiasedEffect(s_hat, denom, e2, r, jack)


def sobolev_smooth(field, alpha: float, spacing=None) -> np.ndarray:
    """Solve (I - alpha * Laplacian) u = f exactly, homogeneous Neumann.

    An H^1 regularisation: `alpha` has units of length squared and sets the
    smoothing length sqrt(alpha). Zero-flux boundaries are the physically right
    choice for a bounded specimen — a periodic solve would wrap dose across the
    cylinder — and they match the 6-connected Neumann stencil the CPM already
    uses (`biofilms_potts.jl:668`, `biofilms_potts_jacc.jl:215`).

    DIAGNOSTIC AND DISPLAY ONLY. This must never enter a gate statistic.
    Smoothing a difference field changes any norm computed on it: it can hide
    real localized structure and it can manufacture apparent structure out of
    noise. Nothing in `feedback_gate` or `synthetic_gate` may import it, and a
    test asserts exactly that.

    Method: even-symmetric extension on every axis realises a DCT, which
    diagonalises the Neumann Laplacian, so the solve is one FFT, one divide and
    one inverse FFT — exact, deterministic, no solver tolerance, and numpy-only.
    scipy is deliberately absent from this package's dependency tier.

    NOTE THE TRAP: applying 1-D solves axis by axis computes
    prod_i (I - alpha d_i^2)^-1, NOT (I - alpha sum_i d_i^2)^-1. Both preserve
    constants and both reduce to the identity at alpha = 0, so only the PDE
    residual distinguishes them. The transform is taken over all axes at once.
    """
    f = np.asarray(field, dtype=float)
    if alpha < 0:
        raise ValueError(f"alpha must be non-negative, got {alpha}")
    if alpha == 0:
        return f.copy()
    h = ((1.0,) * f.ndim if spacing is None
         else tuple(float(s) for s in np.broadcast_to(np.asarray(spacing, float),
                                                      (f.ndim,))))
    if any(s <= 0 for s in h):
        raise ValueError(f"spacing must be positive, got {h}")

    ext = f
    for axis in range(f.ndim):
        ext = np.concatenate([ext, np.flip(ext, axis)], axis=axis)

    lam = np.zeros(ext.shape, dtype=float)
    for axis, step in enumerate(h):
        k = np.fft.fftfreq(ext.shape[axis])
        eig = (4.0 / step ** 2) * np.sin(np.pi * k) ** 2
        shape = [-1 if i == axis else 1 for i in range(f.ndim)]
        lam = lam + eig.reshape(shape)

    u = np.fft.ifftn(np.fft.fftn(ext) / (1.0 + alpha * lam)).real
    return u[tuple(slice(0, n) for n in f.shape)]


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
