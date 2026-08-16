"""The debiased effect estimator, and the smoother that must stay out of it.

WHY A DEBIASED ESTIMATOR AT ALL. The primary metric was a mass-weighted
relative L2 norm, and a norm is UNSIGNED: with zero true effect and residual
Monte Carlo noise it returns |noise|, not 0. Worse, it returns it *stably*, so
replicate scatter never reveals the bias — the pilot's own budget reported a
transport standard deviation of 0.0023 while the metric sat 0.0545 above zero.

The cross-replicate U-statistic removes it by construction rather than by
subtraction: for independent replicates r != s,

    E[d_r^T W d_s] = E[d_r]^T W E[d_s] = ||E[d]||^2_W

with no noise term to estimate and subtract. These tests pin that it is
actually unbiased, that the raw norm on the SAME draws is not, and that the
estimator refuses inputs it cannot support.

Bare tier: pure numpy, no OpenMC, no scipy, no marker.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_openmc.feedback_uq import (DebiasedEffect, cross_replicate_inner,
                                        debiased_squared_effect,
                                        relative_l2_effect, sobolev_smooth)
from biofilm_openmc.synthetic_gate import (NOT_EVALUATED, ThresholdPolicy,
                                           VarianceBudget, decide)

SHAPE = (8, 8, 8)


def _world(seed=0):
    rng = np.random.default_rng(seed)
    mass = rng.uniform(0.5, 1.5, SHAPE)
    mask = rng.random(SHAPE) < 0.6
    mu0 = rng.uniform(0.5, 1.5, SHAPE)
    return rng, mass, mask, mu0


def _replicates(rng, mu0, mass, mask, rel_effect, n_rep, sigma=0.05, rho=0.9):
    """Paired replicates with a known true relative effect.

    `rho` is the common-random-number correlation: within a replicate the two
    states share most of their histories, across replicates they share none.
    """
    scale = np.sqrt((mass[mask] * mu0[mask] ** 2).sum() / mass[mask].sum())
    eff = rel_effect * scale
    e0 = rng.normal(0, sigma, (n_rep,) + SHAPE)
    e1 = rho * e0 + np.sqrt(1 - rho ** 2) * rng.normal(0, sigma, (n_rep,) + SHAPE)
    return ([mu0 + e0[r] for r in range(n_rep)],
            [mu0 + eff + e1[r] for r in range(n_rep)])


# ---------------------------------------------------------------- the estimator

def test_under_the_null_the_estimator_is_centred_on_zero():
    """The property the whole thing exists for. An unbiased estimate of a
    non-negative quantity MUST fall below zero about half the time when the
    truth is zero; a statistic that never goes negative is still biased."""
    rng, mass, mask, mu0 = _world(1)
    vals = []
    for _ in range(200):
        b, f = _replicates(rng, mu0, mass, mask, 0.0, 3)
        vals.append(debiased_squared_effect(b, f, mass, mask).e_squared)
    vals = np.array(vals)
    assert abs(np.mean(vals)) < 3e-4, "the null estimate is not centred on zero"
    assert 0.35 < np.mean(vals > 0) < 0.65, (
        "an unbiased null statistic must straddle zero, not sit above it")


def test_the_raw_norm_is_biased_high_on_the_very_same_draws():
    """Not a claim about norms in general — a measurement on identical data, so
    the comparison cannot be explained by a difference in inputs."""
    rng, mass, mask, mu0 = _world(2)
    raw, debiased = [], []
    for _ in range(120):
        b, f = _replicates(rng, mu0, mass, mask, 0.0, 3)
        raw.append(np.mean([relative_l2_effect(b[r], f[r], mass, mask)
                            for r in range(3)]))
        debiased.append(debiased_squared_effect(b, f, mass, mask).signed_root)
    assert np.median(raw) > 0.01, "the raw norm should show its noise bias here"
    assert abs(np.median(debiased)) < 0.1 * np.median(raw), (
        "the debiased statistic should be far closer to zero than the raw norm")


@pytest.mark.parametrize("n_rep", [3, 5, 20])
@pytest.mark.parametrize("rel_effect", [0.02, 0.05])
def test_a_known_effect_is_recovered_at_every_replicate_count(n_rep, rel_effect):
    """UNBIASEDNESS DOES NOT NEED MANY REPLICATES. R = 3 recovers the truth as
    well as R = 20; what more replicates buy is PRECISION. This is why the
    estimator supersedes rather than confirms the report's ~40x-histories
    recommendation, which was sized to push the raw norm's bias down."""
    rng, mass, mask, mu0 = _world(3)
    got = []
    for _ in range(60):
        b, f = _replicates(rng, mu0, mass, mask, rel_effect, n_rep)
        got.append(debiased_squared_effect(b, f, mass, mask).signed_root)
    assert np.median(got) == pytest.approx(rel_effect, rel=0.08)


def test_more_replicates_buy_precision_not_accuracy():
    rng, mass, mask, mu0 = _world(4)

    def spread(n_rep):
        vals = [debiased_squared_effect(
            *_replicates(rng, mu0, mass, mask, 0.05, n_rep), mass, mask).e_squared
            for _ in range(60)]
        return float(np.std(vals))

    assert spread(20) < spread(3), "the jackknife-era claim: R reduces variance"


def test_the_denominator_is_debiased_too():
    """E[||D0_r||^2_W] carries the baseline's own noise power, so a naive
    squared-norm denominator is biased HIGH and drags every ratio low. With a
    large baseline noise the two constructions must visibly disagree."""
    rng, mass, mask, mu0 = _world(5)
    b, f = _replicates(rng, mu0, mass, mask, 0.05, 6, sigma=0.35)
    got = debiased_squared_effect(b, f, mass, mask)
    naive_den = float(np.sum(mass[mask] * np.asarray(b[0])[mask] ** 2))
    assert got.denominator < naive_den, (
        "the debiased denominator must exclude the baseline noise power that "
        "a single-replicate squared norm includes")
    assert got.e_squared == pytest.approx(0.05 ** 2, rel=0.30)


def test_the_jackknife_tracks_the_actual_sampling_variance():
    """The gate's spread comes from this, so it must not be decorative."""
    rng, mass, mask, mu0 = _world(6)
    estimates, jack = [], []
    for _ in range(120):
        b, f = _replicates(rng, mu0, mass, mask, 0.03, 6)
        d = debiased_squared_effect(b, f, mass, mask)
        estimates.append(d.e_squared)
        jack.append(d.jackknife_var)
    empirical = float(np.var(estimates, ddof=1))
    predicted = float(np.mean(jack))
    assert 0.25 < predicted / empirical < 4.0, (
        f"jackknife {predicted:.3g} should be within a factor of a few of the "
        f"empirical variance {empirical:.3g}")


def test_two_replicates_are_refused_because_their_spread_is_unknowable():
    """The U-statistic itself needs only 2, but the delete-one jackknife then
    has a single leave-one-out value, and a gate input with no estimable
    variance is not a gate input."""
    rng, mass, mask, mu0 = _world(7)
    b, f = _replicates(rng, mu0, mass, mask, 0.05, 2)
    got = debiased_squared_effect(b, f, mass, mask)
    assert np.isnan(got.e_squared) and got.n_replicates == 2


def test_common_random_numbers_reduce_variance_without_creating_bias():
    """Pairing is a variance-reduction technique, not a correctness one — it
    must not move the expectation, only tighten it."""
    rng, mass, mask, mu0 = _world(8)

    def sample(rho):
        return np.array([debiased_squared_effect(
            *_replicates(rng, mu0, mass, mask, 0.04, 4, rho=rho), mass, mask
        ).e_squared for _ in range(80)])

    paired, independent = sample(0.95), sample(0.0)
    assert np.var(paired) < np.var(independent)
    for got in (paired, independent):
        assert np.median(got) == pytest.approx(0.04 ** 2, rel=0.35)


def test_cross_replicate_inner_excludes_the_diagonal():
    """The diagonal is exactly where the noise power lives. With identical
    rows the off-diagonal mean equals the row's own weighted square, so an
    included diagonal would be invisible here — hence the second case, where
    a zero-mean set must give ~0 rather than the positive diagonal total."""
    w = np.ones(4)
    same = np.tile(np.array([1.0, 2.0, 3.0, 4.0]), (3, 1))
    assert cross_replicate_inner(same, w) == pytest.approx(30.0)

    rng = np.random.default_rng(9)
    noise = rng.normal(0, 1, (400, 20))
    assert abs(cross_replicate_inner(noise, np.ones(20))) < 1.0, (
        "pure noise has a large diagonal and a near-zero off-diagonal mean")
    assert np.isnan(cross_replicate_inner(np.ones((1, 4)), w))


# ------------------------------------------------------------------- the gate

def test_a_threshold_declared_for_one_metric_is_refused_for_another():
    """A delta meant for E compared against E^2 is a silent factor-of-ten
    error, and nothing about the number itself reveals it."""
    policy = ThresholdPolicy(metric_id="debiased_relative_l2_squared",
                             effect_threshold=0.01,
                             practical_importance_floor=0.0004)
    draws = np.full(50, 0.02)
    ok = decide(draws, VarianceBudget(), policy,
                metric_id="debiased_relative_l2_squared")
    assert ok.verdict != NOT_EVALUATED

    bad = decide(draws, VarianceBudget(), policy, metric_id="relative_l2")
    assert bad.verdict == NOT_EVALUATED
    assert "not transferable" in bad.reason


def test_a_policy_without_a_metric_identity_is_refused():
    problems = ThresholdPolicy(metric_id="").problems()
    assert any("metric_id" in p for p in problems)


def test_negative_effect_draws_survive_the_verdict_machinery():
    """Null-region draws are legitimately negative, so the gate must handle
    them as data rather than as an error."""
    policy = ThresholdPolicy(metric_id="debiased_relative_l2_squared",
                             effect_threshold=0.01,
                             practical_importance_floor=0.0004)
    draws = np.random.default_rng(10).normal(0.0, 1e-4, 500)
    got = decide(draws, VarianceBudget(transport=1e-8), policy,
                 metric_id=policy.metric_id)
    assert got.verdict != NOT_EVALUATED
    assert got.as_dict()["effect_median"] is not None


# --------------------------------------------------------------- the smoother

def _laplacian(u, spacing):
    out = np.zeros_like(u)
    for axis, step in enumerate(spacing):
        padded = np.concatenate([u.take([0], axis), u, u.take([-1], axis)],
                                axis=axis)
        out += np.diff(padded, 2, axis=axis) / step ** 2
    return out


@pytest.mark.parametrize("spacing", [(1.0, 1.0, 1.0), (0.4, 1.3, 0.8)])
def test_sobolev_actually_solves_the_pde(spacing):
    """THE ONLY TEST THAT DISTINGUISHES THE RIGHT OPERATOR FROM A PLAUSIBLE
    WRONG ONE. Applying 1-D solves axis by axis computes
    prod_i (I - a d_i^2)^-1 instead of (I - a sum_i d_i^2)^-1. That mistake
    preserves constants and reduces to the identity at alpha = 0, so it passes
    both of the tests below; only the residual catches it."""
    f = np.random.default_rng(11).normal(size=(9, 7, 5))
    alpha = 2.0
    u = sobolev_smooth(f, alpha, spacing)
    residual = np.abs(u - alpha * _laplacian(u, spacing) - f).max()
    assert residual < 1e-10, f"residual {residual:.3g} — wrong operator"


def test_sobolev_preserves_constants_and_the_total():
    """Zero-flux boundaries mean nothing leaves the domain, so the mean is
    conserved exactly and a constant field is a fixed point."""
    const = np.full((6, 5, 4), 3.7)
    assert np.allclose(sobolev_smooth(const, 5.0), 3.7)
    f = np.random.default_rng(12).normal(size=(6, 5, 4))
    assert sobolev_smooth(f, 3.0).mean() == pytest.approx(f.mean())


def test_sobolev_edge_cases():
    f = np.random.default_rng(13).normal(size=(5, 5, 5))
    assert np.allclose(sobolev_smooth(f, 0.0), f)
    assert sobolev_smooth(f, 1.0).var() < f.var()
    with pytest.raises(ValueError):
        sobolev_smooth(f, -1.0)
    with pytest.raises(ValueError):
        sobolev_smooth(f, 1.0, spacing=(1.0, 0.0, 1.0))


def test_the_smoother_is_structurally_barred_from_every_gate():
    """Smoothing a difference field changes any norm computed on it: it can
    hide real localized structure and manufacture apparent structure out of
    noise. It is a display and diagnostic tool, and the gates must not be able
    to reach it even by accident."""
    from pathlib import Path
    pkg = Path(__file__).resolve().parents[1] / "biofilm_openmc"
    for module in ("synthetic_gate.py", "feedback_gate.py"):
        src = (pkg / module).read_text()
        assert "sobolev" not in src, (
            f"{module} must not reference the smoother — a gate that can "
            "smooth its own input can be tuned toward a verdict")


def test_debiased_effect_serialises_without_nan():
    """JSON has no NaN, and a verdict only some parsers can read is not
    machine-readable — the same rule S0Verdict already follows."""
    got = DebiasedEffect(float("nan"), 1.0, float("inf"), 3, 2.0)
    as_dict = got.as_dict()
    assert as_dict["s_hat"] is None and as_dict["e_squared"] is None
    assert as_dict["metric_id"] == "debiased_relative_l2_squared"
    assert np.isnan(got.signed_root)
