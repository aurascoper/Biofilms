"""Biomass-field coarse-graining and the occupancy map.

This is the pitch-selection path. It exists because one CPM ID is a
*computational biomass parcel*: segmenting individual cells and treating each
object as one CPM entity would calibrate the wrong thing entirely. What gets
compared is a calibrated 3-D biomass field B(x,y,z) in [0,1] — a binary mask
or a segmentation probability / volume fraction — coarse-grained onto a
candidate lattice:

    phi_j(a) = (1/V_j) * integral_{V_j} B dV

BLOCK MEAN, NOT BLOCK SUM. A biomass volume fraction is INTENSIVE, so
coarse-graining averages. Note that `coupling/biofilm_openmc/mesh.py` does the
opposite deliberately — it block-SUMS, because deposited energy and mass are
extensive. Same word, opposite operation, and using the wrong one is off by
factor**3. The two are separate installable packages and this one must not
import that one (it would be an undeclared dependency and would drag in h5py),
so the distinction is restated here where it can be seen.

THE OCCUPANCY MAP IS PART OF THE CALIBRATION. The CPM lattice is binary, so
turning phi into occupied/unoccupied is a modelling decision with consequences,
and it cannot stay implicit:

    threshold        occupied iff phi >= tau. Deterministic, preserves
                     connectivity, does NOT conserve total biovolume.
    mass_preserving  occupied with probability phi, seeded. Conserves total
                     biovolume in expectation, varies by seed.

Neither is correct in general; one must be declared, and the acceptance config
refuses while it is unset.
"""

from __future__ import annotations

import numpy as np

OCCUPANCY_MAPPINGS = ("threshold", "mass_preserving")


class FieldError(ValueError):
    """Raised on an invalid field, factor, or occupancy declaration."""


def check_field(B: np.ndarray) -> np.ndarray:
    """A biomass field is 3-D and lies in [0, 1]."""
    a = np.asarray(B, dtype=float)
    if a.ndim != 3:
        raise FieldError(f"biomass field must be 3-D, got {a.ndim}-D")
    if a.size == 0:
        raise FieldError("biomass field is empty")
    if np.isnan(a).any():
        raise FieldError("biomass field contains NaN")
    lo, hi = float(a.min()), float(a.max())
    if lo < 0.0 or hi > 1.0:
        raise FieldError(
            f"biomass field must lie in [0, 1], got [{lo}, {hi}] — a "
            "segmentation probability or volume fraction, not a raw intensity")
    return a


def coarse_grain(B: np.ndarray, factor: int) -> np.ndarray:
    """Block-MEAN a biomass volume fraction by an integer factor.

    Conserves the total biovolume fraction exactly, which is the property that
    makes cross-resolution comparison meaningful. See the module docstring for
    why this is a mean where the transport package uses a sum.
    """
    a = check_field(B)
    f = int(factor)
    if f < 1:
        raise FieldError(f"coarsening factor must be >= 1, got {f}")
    if f == 1:
        return a
    bad = [d for d in a.shape if d % f]
    if bad:
        raise FieldError(
            f"coarsening factor {f} does not divide field shape {a.shape} "
            "evenly — partial blocks would be averaged over fewer voxels and "
            "silently reweight the edges")
    nx, ny, nz = (d // f for d in a.shape)
    return a.reshape(nx, f, ny, f, nz, f).mean(axis=(1, 3, 5))


def biovolume_fraction(B: np.ndarray) -> float:
    """Mean biomass fraction over the whole field. Invariant under
    coarse_grain, by construction."""
    return float(check_field(B).mean())


def occupancy_threshold(phi: np.ndarray, tau: float) -> np.ndarray:
    """Occupied iff phi >= tau. Deterministic; does not conserve biovolume."""
    a = check_field(phi)
    if not 0.0 < float(tau) <= 1.0:
        raise FieldError(f"threshold tau must be in (0, 1], got {tau}")
    return a >= float(tau)


def occupancy_mass_preserving(phi: np.ndarray, seed: int) -> np.ndarray:
    """Occupied with probability phi. Conserves total biovolume in
    expectation; the realisation varies by seed, so structural observables
    need replicate handling."""
    a = check_field(phi)
    if seed is None:
        raise FieldError(
            "mass_preserving occupancy is stochastic and requires a declared "
            "seed — an unseeded lattice cannot be reproduced or reviewed")
    rng = np.random.default_rng(int(seed))
    return rng.random(a.shape) < a


def apply_occupancy(phi: np.ndarray, mapping: str, *, tau: float | None = None,
                    seed: int | None = None) -> np.ndarray:
    """Dispatch on a DECLARED mapping. Refuses when it is unset, exactly as
    the pitch thresholds do."""
    if not mapping:
        raise FieldError(
            "no occupancy mapping declared — turning a biomass fraction into a "
            "binary lattice is part of the calibration, not an implementation "
            f"detail. Declare one of {list(OCCUPANCY_MAPPINGS)} in "
            "config/cpm_spatial_acceptance_template.toml")
    if mapping == "threshold":
        if tau is None:
            raise FieldError("occupancy mapping 'threshold' requires threshold_tau")
        return occupancy_threshold(phi, tau)
    if mapping == "mass_preserving":
        return occupancy_mass_preserving(phi, seed)
    raise FieldError(
        f"unknown occupancy mapping {mapping!r}; expected one of "
        f"{list(OCCUPANCY_MAPPINGS)}")


def occupancy_biovolume_error(phi: np.ndarray, mapping: str, *,
                              tau: float | None = None,
                              seed: int | None = None) -> float:
    """Relative change in total biovolume caused by the occupancy map.

    Reported rather than corrected: a threshold map systematically gains or
    loses biomass depending on tau, and the size of that bias is a property of
    the declared mapping that the calibration record should carry.
    """
    a = check_field(phi)
    target = float(a.mean())
    got = float(apply_occupancy(a, mapping, tau=tau, seed=seed).mean())
    if target == 0.0:
        return float("nan")
    return (got - target) / target
