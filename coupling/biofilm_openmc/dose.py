"""Dose normalization: heating tally (eV/source-particle) -> Gy/s.

OpenMC documents `heating` in eV per source particle. The dimensioned dose
rate for voxel v is

    dose_Gy_s[v] = H_v [eV/src] * 1.602176634e-19 [J/eV] * S [src/s] / m_v [kg]

Results are voxel-averaged absorbed-dose estimates under OpenMC's
charged-particle local-deposition approximation — NOT single-cell
microdosimetry.

Statepoint access is duck-typed on `.get_tally(name=...)` returning an
object with `.mean` / `.std_dev`, so the same code path serves the real
`openmc.StatePoint` and the test fixtures' FakeStatepoint.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

EV_TO_J = 1.602176634e-19


@dataclass
class DoseResult:
    dose_rate_mean_Gy_s: np.ndarray   # logical (x,y,z)
    dose_rate_sd_Gy_s: np.ndarray
    rel_err: np.ndarray               # sd/mean where mean > 0, else inf
    source_rate: float


def normalize_heating(heating_eV_per_src: np.ndarray, source_rate_per_s: float,
                      mass_kg: np.ndarray) -> np.ndarray:
    """The single-line formula, kept pure for exact unit tests."""
    if np.any(mass_kg <= 0):
        raise ValueError("non-positive voxel mass — dose undefined")
    return heating_eV_per_src * EV_TO_J * source_rate_per_s / mass_kg


def extract_heating(statepoint, mesh_shape) -> tuple[np.ndarray, np.ndarray]:
    """Heating tally mean/sd reshaped to logical (x,y,z).

    OpenMC mesh-filter bins run x fastest, z slowest (Fortran-like), so a
    C-order reshape to (z,y,x) then transpose recovers logical order.
    """
    tally = statepoint.get_tally(name="heating")
    n = int(np.prod(mesh_shape))
    mean = np.asarray(tally.mean).reshape(mesh_shape[::-1]).transpose(2, 1, 0)
    sd = np.asarray(tally.std_dev).reshape(mesh_shape[::-1]).transpose(2, 1, 0)
    if mean.size != n:  # pragma: no cover
        raise ValueError("tally size does not match mesh")
    return mean, sd


def dose_from_statepoint(statepoint, mass_kg: np.ndarray,
                         source_rate_per_s: float) -> DoseResult:
    mean_eV, sd_eV = extract_heating(statepoint, mass_kg.shape)
    mean = normalize_heating(mean_eV, source_rate_per_s, mass_kg)
    sd = normalize_heating(sd_eV, source_rate_per_s, mass_kg)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(mean > 0, sd / mean, np.inf)
    return DoseResult(mean, sd, rel, source_rate_per_s)


def sparsity_report(result: DoseResult, occupied_mask: np.ndarray) -> dict:
    """Feasibility numbers the integration benchmark must report (review
    amendment 11): scoring sparsity and uncertainty on occupied voxels."""
    occ = result.dose_rate_mean_Gy_s[occupied_mask]
    rel = result.rel_err[occupied_mask]
    finite = rel[np.isfinite(rel)]
    return {
        "occupied_voxels": int(occupied_mask.sum()),
        "occupied_unscored_fraction": float((occ == 0).mean()) if occ.size else float("nan"),
        "median_rel_err": float(np.median(finite)) if finite.size else float("nan"),
        "max_rel_err": float(finite.max()) if finite.size else float("nan"),
    }
