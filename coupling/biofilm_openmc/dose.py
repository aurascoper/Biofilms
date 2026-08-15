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
    """Dosimetry-stage field: absorbed dose RATE, Gy/s."""
    dose_rate_mean_Gy_s: np.ndarray   # logical (x,y,z)
    dose_rate_sd_Gy_s: np.ndarray
    rel_err: np.ndarray               # sd/mean where mean > 0, else inf
    source_rate: float

    unit = "Gy/s"

    # Unit-neutral accessors so scale-invariant comparisons (effect/noise
    # ratios) run identically on this and on PerSourceResult.
    @property
    def field(self) -> np.ndarray:
        return self.dose_rate_mean_Gy_s

    @property
    def field_sd(self) -> np.ndarray:
        return self.dose_rate_sd_Gy_s


@dataclass
class PerSourceResult:
    """Transport-stage field: specific energy per source particle, Gy per
    source particle (J/kg per source particle).

    This is what transport actually yields. Multiplying by the emission rate
    (photons/s) gives a DoseResult — so a source activity is NEVER needed to
    run or compare transport. Do not fabricate `photons_per_second = 1.0` to
    reach this: that would silently relabel Gy/source as Gy/s.
    """
    specific_energy_mean_Gy_per_src: np.ndarray
    specific_energy_sd_Gy_per_src: np.ndarray
    rel_err: np.ndarray

    unit = "Gy/source-particle"

    @property
    def field(self) -> np.ndarray:
        return self.specific_energy_mean_Gy_per_src

    @property
    def field_sd(self) -> np.ndarray:
        return self.specific_energy_sd_Gy_per_src


def specific_energy_per_source(heating_eV_per_src: np.ndarray,
                               mass_kg: np.ndarray) -> np.ndarray:
    """Gy per source particle. The activity-free half of the dose formula.

    An EMPTY bin — no mass and no heating — is dose zero, not an error. Exact
    CSG masses (`materials.mesh_material_masses_kg`) produce them wherever the
    mesh cube reaches outside the geometry, which is every corner bin. The
    full-bin approximations never did, which is why this guard could once be
    unconditional.

    Mass with the WRONG SIGN is still refused, and so is energy deposited in a
    bin holding nothing — that combination is a real inconsistency between the
    tally and the geometry, and silently returning zero would hide it.
    """
    if np.any(mass_kg < 0):
        raise ValueError("negative voxel mass — dose undefined")
    empty = mass_kg == 0
    if np.any(empty & (heating_eV_per_src != 0)):
        raise ValueError(
            "heating scored in a bin with zero mass — the tally and the "
            "geometry disagree about where material is")
    with np.errstate(divide="ignore", invalid="ignore"):
        out = heating_eV_per_src * EV_TO_J / mass_kg
    return np.where(empty, 0.0, out)


def normalize_heating(heating_eV_per_src: np.ndarray, source_rate_per_s: float,
                      mass_kg: np.ndarray) -> np.ndarray:
    """The single-line formula, kept pure for exact unit tests."""
    return specific_energy_per_source(heating_eV_per_src, mass_kg) * source_rate_per_s


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


def per_source_from_statepoint(statepoint, mass_kg: np.ndarray) -> PerSourceResult:
    """Transport-stage result: no source activity required."""
    mean_eV, sd_eV = extract_heating(statepoint, mass_kg.shape)
    mean = specific_energy_per_source(mean_eV, mass_kg)
    sd = specific_energy_per_source(sd_eV, mass_kg)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(mean > 0, sd / mean, np.inf)
    return PerSourceResult(mean, sd, rel)


def dose_from_per_source(result: PerSourceResult, source_rate_per_s: float) -> DoseResult:
    """Apply the emission rate to a transport-stage field. The relative error
    is carried over unchanged: it is scale-invariant."""
    return DoseResult(result.field * source_rate_per_s,
                      result.field_sd * source_rate_per_s,
                      result.rel_err, source_rate_per_s)


def sparsity_report(result, occupied_mask: np.ndarray) -> dict:
    """Feasibility numbers the integration benchmark must report (review
    amendment 11): scoring sparsity and uncertainty on occupied voxels.
    Accepts either result type — these statistics are unit-independent."""
    occ = result.field[occupied_mask]
    rel = result.rel_err[occupied_mask]
    finite = rel[np.isfinite(rel)]
    return {
        "occupied_voxels": int(occupied_mask.sum()),
        "occupied_unscored_fraction": float((occ == 0).mean()) if occ.size else float("nan"),
        "median_rel_err": float(np.median(finite)) if finite.size else float("nan"),
        "max_rel_err": float(finite.max()) if finite.size else float("nan"),
    }
