"""Monod fitter: code verification on synthetic curves only.

There is no measured depletion curve in this repository. These tests generate
curves from known mu_max and monod_half_saturation, add Gaussian noise, and
check that the fitter recovers them, reports non-identifiability when the
design cannot bound K, and that its asymptotic intervals cover at a rate
consistent with nominal 95%. None of this is biological calibration.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest

from biofilm_calibration.monod import (CHI2_1_95, MonodFitError, fit_monod,
                                       simulate_monod_batch)

# Truth and a design whose substrate range spans K: S0 = 25 K, sampled through
# depletion (exhausted near t ~ 8).
TRUE_MU, TRUE_K = 0.5, 0.2
Y, X0, S0 = 0.5, 0.05, 5.0
SIGMA = 0.05
BOUNDS = dict(mu_max_bounds=(0.05, 5.0), half_saturation_bounds=(2e-3, 20.0))


def _curve(t, rng):
    S, _ = simulate_monod_batch(t, TRUE_MU, TRUE_K, Y=Y, X0=X0, S0=S0)
    return S + rng.normal(0.0, SIGMA, size=S.shape)


def _contains(ci, value):
    lo, hi = ci
    return (lo is None or lo <= value) and (hi is None or value <= hi)


def test_linear_invariant_is_conserved_to_roundoff():
    # d(S + X/Y)/dt = 0 identically, and every Runge-Kutta method conserves
    # linear invariants EXACTLY, whatever the step. So this checks that the
    # right-hand side is wired consistently, at roundoff tolerance; it cannot
    # check integration accuracy (see the convergence test for that).
    t = np.linspace(0.0, 12.0, 25)
    S, X = simulate_monod_batch(t, TRUE_MU, TRUE_K, Y=Y, X0=X0, S0=S0, h_max=0.2)
    invariant0 = S0 + X0 / Y
    np.testing.assert_allclose(S + X / Y, invariant0, rtol=1e-12)


def test_rk4_converges_at_fourth_order():
    # Integration accuracy, separately: halving h should cut the error ~16x.
    # Window ends before depletion so the S >= 0 clamp is never active.
    t = np.array([0.0, 6.0])
    S = [simulate_monod_batch(t, TRUE_MU, TRUE_K, Y=Y, X0=X0, S0=S0,
                              h_max=h)[0][-1] for h in (0.4, 0.2, 0.1)]
    ratio = (S[0] - S[1]) / (S[1] - S[2])
    assert 12.0 < ratio < 20.0, ratio


def test_recovers_parameters_when_design_spans_K():
    rng = np.random.default_rng(20261002)
    t = np.linspace(0.0, 12.0, 25)
    fit = fit_monod(t, _curve(t, rng), Y=Y, X0=X0, S0=S0, **BOUNDS)

    assert fit.identifiable == {"mu_max": True, "monod_half_saturation": True}
    # No "interval contains the truth" here: one draw misses ~5% of the time by
    # construction, and this seed's K interval (0.205, 0.476) does. Containment
    # is a rate, and the coverage test below measures it.
    assert abs(fit.mu_max / TRUE_MU - 1) < 0.10
    # K is the weakly identified one; a factor-of-2 band is the honest claim.
    assert 0.5 < fit.monod_half_saturation / TRUE_K < 2.0


def test_half_saturation_unbounded_when_substrate_never_approaches_K():
    # Stop sampling at t = 5, where S ~ 3.9 = 20 K. The curve is effectively
    # zero-order in S, so it cannot rule out arbitrarily small K.
    rng = np.random.default_rng(7)
    t = np.linspace(0.0, 5.0, 21)
    fit = fit_monod(t, _curve(t, rng), Y=Y, X0=X0, S0=S0, **BOUNDS)

    assert fit.half_saturation_ci[0] is None
    assert fit.identifiable["monod_half_saturation"] is False


def _binomial_acceptance(n: int, p: float, alpha: float) -> tuple[int, int]:
    """Smallest lo and largest hi with P(X < lo) <= alpha/2, P(X > hi) <= alpha/2."""
    pmf = [math.comb(n, k) * p**k * (1 - p) ** (n - k) for k in range(n + 1)]
    cdf = np.cumsum(pmf)
    lo = int(np.searchsorted(cdf, alpha / 2, side="right"))
    sf = 1.0 - np.concatenate([[0.0], cdf[:-1]])  # sf[k] = P(X >= k)
    hi = int(max(k for k in range(n + 1) if sf[k] > alpha / 2))
    return lo, hi


# Predeclared BEFORE looking at any coverage result: exact Binomial(100, 0.95),
# two-sided, alpha = 0.01. Failing means the asymptotic interval misbehaves for
# this design, not that the test needs a wider band.
N_REPLICATES, NOMINAL, ALPHA = 100, 0.95, 0.01


def test_asymptotic_interval_coverage_is_consistent_with_nominal():
    lo, hi = _binomial_acceptance(N_REPLICATES, NOMINAL, ALPHA)
    rng = np.random.default_rng(1)
    t = np.linspace(0.0, 12.0, 25)
    hits = {"mu_max": 0, "monod_half_saturation": 0}
    start = time.perf_counter()
    for _ in range(N_REPLICATES):
        fit = fit_monod(t, _curve(t, rng), Y=Y, X0=X0, S0=S0, n_grid=120,
                        **BOUNDS)
        hits["mu_max"] += _contains(fit.mu_max_ci, TRUE_MU)
        hits["monod_half_saturation"] += _contains(fit.half_saturation_ci, TRUE_K)
    elapsed = time.perf_counter() - start
    print(f"\ncoverage over {N_REPLICATES} replicates: {hits}; "
          f"acceptance [{lo}, {hi}]; {elapsed:.1f} s")
    for name, count in hits.items():
        assert lo <= count <= hi, (name, count, (lo, hi))


def test_cutoff_is_chi2_1_quantile():
    assert CHI2_1_95 == pytest.approx(3.8415, abs=1e-4)


def test_rejects_unusable_inputs():
    t = np.linspace(0.0, 1.0, 5)
    with pytest.raises(MonodFitError):
        fit_monod(t, np.ones(4), Y=Y, X0=X0, S0=S0, **BOUNDS)
    with pytest.raises(MonodFitError):
        fit_monod(t, np.ones(5), Y=Y, X0=X0, S0=S0,
                  mu_max_bounds=(1.0, 0.5), half_saturation_bounds=(0.1, 1.0))
