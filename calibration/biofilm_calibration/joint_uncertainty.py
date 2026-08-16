"""Joint sampling for correlated and constrained calibration inputs.

`materials.uncertainty.propagate` samples every input independently and says so
in its own docstring. That is fine for a single ratio and wrong for the feedback
experiment, where two constraints are structural rather than incidental:

    CORRELATION IS NOT OPTIONAL. Wet and dry mass of one coupon are weighed on
    one balance and share its calibration error. Sampling them independently
    inflates or deflates the water fraction depending on which way the error
    happens to fall, and the resulting interval describes a measurement nobody
    made.

    COMPOSITION IS COMPOSITIONAL DATA. Elemental mass fractions are
    non-negative and sum to one. Perturbing them with independent Gaussians
    leaves the simplex immediately — negative fractions, sums drifting off 1 —
    and the transport loader then refuses the config, or worse, a renormalising
    caller hides the drift. Oxygen-by-difference makes it worse: it induces
    strong NEGATIVE correlations by construction, since every other element's
    error lands in it.

So this module samples groups jointly and keeps constrained quantities on their
constraint set, rather than sampling freely and repairing afterwards. Repair is
what loses the correlation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class CorrelationGroup:
    """Quantities that share a measurement error, sampled together.

    `shared_relative_sigma` is the error common to the whole group — the balance
    calibration, the segmentation threshold, the instrument drift. Each member
    additionally carries its own independent component. Sampling the shared part
    ONCE per draw is the entire point: it is what makes the members correlated.
    """
    name: str
    values: dict                      # member -> central value
    shared_relative_sigma: float = 0.0
    independent_relative_sigma: dict = field(default_factory=dict)

    def sample(self, rng: np.random.Generator, draws: int) -> dict:
        shared = rng.normal(0.0, self.shared_relative_sigma, draws) \
            if self.shared_relative_sigma else np.zeros(draws)
        out = {}
        for name, value in self.values.items():
            own_sigma = self.independent_relative_sigma.get(name, 0.0)
            own = rng.normal(0.0, own_sigma, draws) if own_sigma else np.zeros(draws)
            out[name] = float(value) * (1.0 + shared + own)
        return out


def sample_composition(central: dict, rng: np.random.Generator, draws: int, *,
                       log_sigma: float = 0.05) -> dict:
    """Draws from the simplex, closed by construction rather than by repair.

    Logistic-normal: perturb in log space, then normalise. Every draw is
    non-negative and sums to exactly 1 up to floating point, so no draw can
    leave the constraint set and none needs renormalising after the fact —
    which is what would silently destroy the induced negative correlations.

    Those correlations are real: raising one element's fraction must lower the
    others, because they sum to one. A sampler that produced independent
    fractions would be describing a different material at every draw.
    """
    keys = sorted(central)
    base = np.array([float(central[k]) for k in keys], dtype=float)
    if np.any(base < 0):
        raise ValueError("central composition has a negative mass fraction")
    total = base.sum()
    if total <= 0:
        raise ValueError("central composition sums to zero")
    base = base / total

    # Perturb only the strictly positive components; a declared zero stays zero.
    positive = base > 0
    logs = np.zeros((draws, base.size))
    logs[:, positive] = (np.log(base[positive])
                         + rng.normal(0.0, log_sigma, (draws, positive.sum())))
    weights = np.zeros((draws, base.size))
    weights[:, positive] = np.exp(logs[:, positive])
    weights /= weights.sum(axis=1, keepdims=True)
    return {k: weights[:, i] for i, k in enumerate(keys)}


def closure_error(sampled: dict) -> np.ndarray:
    """|sum of fractions - 1| per draw. Zero to floating point, or the sampler
    is broken and the transport loader will refuse every config it feeds."""
    stacked = np.vstack([np.asarray(v, dtype=float) for v in sampled.values()])
    return np.abs(stacked.sum(axis=0) - 1.0)


def observed_correlation(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


@dataclass(frozen=True)
class JointSample:
    """One outer design, with its provenance. Seeded, because an interval that
    changes between runs cannot be reviewed, cited or regression-tested."""
    draws: int
    seed: int
    groups: dict = field(default_factory=dict)
    compositions: dict = field(default_factory=dict)

    def flat(self) -> dict:
        out = dict(self.groups)
        for name, comp in self.compositions.items():
            for element, values in comp.items():
                out[f"{name}.{element}"] = values
        return out


def sample_joint(groups: list[CorrelationGroup],
                 compositions: dict, *, seed: int, draws: int = 4,
                 log_sigma: float = 0.05) -> JointSample:
    """Sample every group and composition from ONE seeded stream.

    One stream rather than one per quantity, so the whole design is reproducible
    from a single recorded integer.
    """
    rng = np.random.default_rng(seed)
    sampled_groups: dict = {}
    for g in groups:
        sampled_groups.update(g.sample(rng, draws))
    sampled_comps = {name: sample_composition(central, rng, draws,
                                              log_sigma=log_sigma)
                     for name, central in compositions.items()}
    return JointSample(draws=draws, seed=seed, groups=sampled_groups,
                       compositions=sampled_comps)
