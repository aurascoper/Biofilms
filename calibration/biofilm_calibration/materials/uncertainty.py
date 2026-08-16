"""Seeded, reproducible uncertainty propagation.

Monte Carlo rather than a derivative expansion, because the quantities of
interest are ratios of measured masses and volumes and are not well
approximated by a first-order expansion when the relative uncertainties are
large — which, for a wet biofilm mass, they are.

SEEDED ON PURPOSE. An uncertainty interval that changes between runs cannot be
reviewed, cited, or regression-tested. Every function here takes a seed and the
same seed gives the same interval.
"""

from __future__ import annotations

import numpy as np


def propagate(fn, inputs: dict, *, seed: int, draws: int = 20000,
              quantiles=(0.025, 0.5, 0.975), invalid_if=None) -> dict:
    """Propagate independent Gaussian uncertainties through `fn`.

    `inputs` maps a name to (value, sigma). A sigma of 0 or None is treated as
    exact. Independence is an ASSUMPTION and often wrong — wet and dry mass of
    the same sample share a balance calibration — so correlated inputs must be
    combined before they get here, not declared independent for convenience.
    """
    rng = np.random.default_rng(seed)
    names = list(inputs)
    samples = {}
    for name in names:
        value, sigma = inputs[name]
        if sigma in (None, 0):
            samples[name] = np.full(draws, float(value))
        else:
            samples[name] = rng.normal(float(value), float(sigma), draws)

    out = np.array([fn(**{n: samples[n][i] for n in names}) for i in range(draws)])
    finite = out[np.isfinite(out)]

    # A draw can be perfectly finite and still physically impossible — a blank
    # exceeding its own sample gives a NEGATIVE biofilm mass, which would sail
    # through the finiteness filter and out the other side as a negative
    # density. That is information about the measurement, not noise to discard:
    # it says the distribution crosses zero, which no number of replicates can
    # fix. `invalid_if` names the impossibility; the fraction is reported.
    invalid_fraction = 0.0
    if invalid_if is not None and finite.size:
        invalid_fraction = float(np.mean([bool(invalid_if(v)) for v in finite]))

    if finite.size == 0:
        return {"n_valid": 0, "n_rejected": int(draws),
                "invalid_fraction": float("nan")}
    qs = np.quantile(finite, quantiles)
    return {
        "n_valid": int(finite.size),
        "n_rejected": int(draws - finite.size),
        "invalid_fraction": invalid_fraction,
        "mean": float(finite.mean()),
        "sd": float(finite.std(ddof=1)) if finite.size > 1 else 0.0,
        **{f"q{int(q * 1000):03d}": float(v) for q, v in zip(quantiles, qs)},
    }


def relative_uncertainty(summary: dict) -> float:
    """sd/mean, or nan when the mean is zero or nothing was valid."""
    if not summary.get("n_valid") or summary.get("mean", 0.0) == 0.0:
        return float("nan")
    return summary["sd"] / abs(summary["mean"])
