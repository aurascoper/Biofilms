"""Monod parameters from a batch substrate-depletion curve, with intervals.

Substrate S and biomass X in a closed batch, yield Y and initial biomass X0
known (measured separately), S0 known:

    dS/dt = -(mu_max / Y) * S / (K + S) * X
    dX/dt =   mu_max      * S / (K + S) * X

K here is the MONOD HALF-SATURATION, named `monod_half_saturation` everywhere
it leaves this module. It is not the carrying capacity `K_s` of the preprint's
Table 2 (claims_ledger PP-T2-21..24); the two share a letter and nothing else.

mu_max and K are strongly correlated and poorly identified from one curve
unless the substrate range actually spans K: while S >> K the curve is
zero-order and says nothing about K. A fit that returns a point estimate
without saying so is the same failure as a D10 quoted without its endpoint, so
this returns asymptotic one-parameter PROFILE-LIKELIHOOD intervals and marks a
parameter non-identifiable when its interval runs off the grid.

NumPy only, by the same rule as the rest of this package: no optimizer. Two
parameters are profiled exactly on a log grid; the ODE is integrated by
fixed-step RK4 over the whole grid at once.

NOT CONNECTED TO THE CPM. `seconds_per_mcs` is NaN by contract and C_wall has
no physical concentration, so a fitted mu_max [1/time] and K [concentration]
have no lattice counterpart. The CPM's `K_C_monod` is a declared, dimensionless
value; nothing here feeds it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# chi-square(1) 0.95 quantile: the asymptotic likelihood-ratio cutoff for one
# parameter. The resulting intervals are asymptotic, not exact finite-sample
# 95% intervals; tests/test_monod_fit.py measures their actual coverage.
CHI2_1_95 = 3.841458820694124


class MonodFitError(ValueError):
    """Raised when a Monod fit is asked of unusable inputs."""


@dataclass(frozen=True)
class MonodFit:
    mu_max: float
    monod_half_saturation: float
    # (lo, hi); None on a side where the interval reached the grid edge.
    mu_max_ci: tuple[float | None, float | None]
    half_saturation_ci: tuple[float | None, float | None]
    identifiable: dict[str, bool]
    n: int
    rss: float


def simulate_monod_batch(t, mu_max, K, *, Y: float, X0: float, S0: float,
                         h_max: float = 0.05):
    """Integrate the batch Monod system and return (S, X) at times `t`.

    `mu_max` and `K` may be arrays (broadcast together); the outputs then have
    shape broadcast(mu_max, K).shape + (len(t),). `t` must be non-decreasing
    and start at the time S0 and X0 refer to. Each interval between samples is
    split into equal RK4 steps no longer than `h_max`.
    """
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.size < 2 or np.any(np.diff(t) < 0):
        raise MonodFitError("t must be a non-decreasing 1-D array of >= 2 times")
    mu, k = np.broadcast_arrays(np.asarray(mu_max, float), np.asarray(K, float))
    S = np.full(mu.shape, float(S0))
    X = np.full(mu.shape, float(X0))

    def rates(S, X):
        # Clamp only inside the rate: dS and dX share the factor, so the
        # linear invariant S + X/Y is still conserved exactly.
        Sp = np.maximum(S, 0.0)
        growth = mu * Sp / (k + Sp) * X
        return -growth / Y, growth

    S_out = np.empty(mu.shape + (t.size,))
    X_out = np.empty_like(S_out)
    S_out[..., 0], X_out[..., 0] = S, X
    for i in range(1, t.size):
        dt = t[i] - t[i - 1]
        steps = max(1, int(np.ceil(dt / h_max)))
        h = dt / steps
        for _ in range(steps):
            a1, b1 = rates(S, X)
            a2, b2 = rates(S + h / 2 * a1, X + h / 2 * b1)
            a3, b3 = rates(S + h / 2 * a2, X + h / 2 * b2)
            a4, b4 = rates(S + h * a3, X + h * b3)
            S = S + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            X = X + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        S_out[..., i], X_out[..., i] = S, X
    return S_out, X_out


def _interval_from_profile(grid: np.ndarray, delta: np.ndarray,
                           threshold: float) -> tuple[float | None, float | None]:
    """Turn a 1-D profile into an interval on `grid`.

    `grid` is ascending; `delta[i]` is the profile -2logL at grid[i] minus its
    minimum, so min(delta) == 0. The interval is the set of grid values with
    delta <= threshold. A side that reaches the first or last grid point is
    returned as None: the data did not bound the parameter there, and the
    grid edge is a choice of the caller, not a result.
    """
    admitted = np.flatnonzero(delta <= threshold)
    if admitted.size == 0:
        raise MonodFitError("profile has no points at or below threshold")
    lo_i, hi_i = admitted[0], admitted[-1]
    # Outer span: a (lo, hi) pair cannot represent a disjoint admitted set, and
    # keeping only the MLE's component would silently discard a second mode.
    # The hull loses information in the conservative direction. If
    # multimodality matters, change the API to return several intervals.
    return (None if lo_i == 0 else float(grid[lo_i]),
            None if hi_i == len(grid) - 1 else float(grid[hi_i]))


def fit_monod(t, S_obs, *, Y: float, X0: float, S0: float,
              mu_max_bounds: tuple[float, float],
              half_saturation_bounds: tuple[float, float],
              n_grid: int = 200, h_max: float = 0.05,
              max_widen: int = 8) -> MonodFit:
    """Fit mu_max and the Monod half-saturation to a depletion curve.

    Gaussian errors of unknown, constant variance; sigma is profiled out
    analytically, so -2logL = n log(RSS/n) + const. Both parameters are
    profiled on log-spaced n_grid x n_grid grids: a coarse pass over the given
    bounds, then re-grids of a box around the region it admits, widened until
    the region no longer touches a box edge short of the caller's bounds. The
    MLE is the final grid minimum. An interval side reported None reached a
    caller bound: the data did not bound it within what was searched.
    """
    t = np.asarray(t, dtype=float)
    S_obs = np.asarray(S_obs, dtype=float)
    if S_obs.shape != t.shape:
        raise MonodFitError("t and S_obs must have the same shape")
    if not (Y > 0 and X0 > 0 and S0 > 0):
        raise MonodFitError("Y, X0 and S0 must be positive")
    for lo, hi in (mu_max_bounds, half_saturation_bounds):
        if not 0 < lo < hi:
            raise MonodFitError("bounds must satisfy 0 < lo < hi")

    n = t.size

    def surface(mu_grid, k_grid):
        S_pred, _ = simulate_monod_batch(t, mu_grid[:, None], k_grid[None, :],
                                         Y=Y, X0=X0, S0=S0, h_max=h_max)
        rss = np.sum((S_pred - S_obs) ** 2, axis=-1)
        return rss, n * np.log(rss / n)

    # A coarse grid over the caller's bounds only LOCATES the admitted region.
    # Its intervals are not usable: the (mu_max, K) ridge is thin and diagonal
    # and falls between grid points, so a profile read off a grid spaced like
    # the interval width is too narrow (measured: 61/100 coverage at 4%
    # spacing) and it under-reports the region itself. So re-grid a box around
    # the region, and widen any side the region still touches until it either
    # stops touching or reaches the caller's bound. Only a caller bound can
    # then produce an unbounded (None) interval side.
    bounds = np.log([mu_max_bounds, half_saturation_bounds])  # (axis, lo/hi)
    box = bounds.copy()
    for attempt in range(max_widen):
        grids = [np.geomspace(*np.exp(b), n_grid) for b in box]
        rss, m2ll = surface(*grids)
        admitted = np.nonzero(m2ll - m2ll.min() <= CHI2_1_95)
        lo_hit = [idx.min() == 0 and box[ax, 0] > bounds[ax, 0]
                  for ax, idx in enumerate(admitted)]
        hi_hit = [idx.max() == n_grid - 1 and box[ax, 1] < bounds[ax, 1]
                  for ax, idx in enumerate(admitted)]
        if attempt > 0 and not any(lo_hit + hi_hit):
            break
        new = np.empty_like(box)
        for ax, idx in enumerate(admitted):
            lo, hi = np.log(grids[ax][[idx.min(), idx.max()]])
            pad = max(hi - lo, (box[ax, 1] - box[ax, 0]) / n_grid)
            box_width = box[ax, 1] - box[ax, 0]
            new[ax] = (lo - (box_width if lo_hit[ax] else pad),
                       hi + (box_width if hi_hit[ax] else pad))
        box = np.clip(new, bounds[:, :1], bounds[:, 1:])
    else:
        raise MonodFitError("admitted region did not settle within "
                            f"{max_widen} re-grids")
    mu_grid, k_grid = grids

    i, j = np.unravel_index(np.argmin(m2ll), m2ll.shape)
    prof_mu = m2ll.min(axis=1)
    prof_k = m2ll.min(axis=0)
    mu_ci = _interval_from_profile(mu_grid, prof_mu - prof_mu.min(), CHI2_1_95)
    k_ci = _interval_from_profile(k_grid, prof_k - prof_k.min(), CHI2_1_95)
    return MonodFit(
        mu_max=float(mu_grid[i]),
        monod_half_saturation=float(k_grid[j]),
        mu_max_ci=mu_ci,
        half_saturation_ci=k_ci,
        identifiable={"mu_max": None not in mu_ci,
                      "monod_half_saturation": None not in k_ci},
        n=n,
        rss=float(rss[i, j]),
    )
