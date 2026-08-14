"""Lattice-scale candidate enumeration, and why it cannot select one.

THE IDENTIFIABILITY PROBLEM. A measured physical volume constrains only the
product:

    V_physical = V_sites * a**3

so the lattice pitch `a` and the target site count `V_sites` are NOT separately
identifiable from a volume measurement alone. Dividing one reported cell volume
by the current `V_target = 120` and calling the cube root a calibrated pitch is
choosing `V_sites` by accident and then deriving `a` from that accident.

Breaking the degeneracy needs a second, independent constraint — one of:

  * a declared physical domain size, which fixes `a = domain / N`;
  * a resolution criterion, e.g. "the smallest dimension of interest must span
    at least k voxels", which puts a ceiling on `a`;
  * a declared target site count, which is a modelling choice, not a measurement.

This module enumerates the family and reports what each member implies. It
never picks one, and there is a test that it never picks one.

DOMAIN CONSTRAINT. The Julia model's geometry is an N x N x N cube containing a
cylinder of radius N/2 (`R = N / 2.0`, recomputed as a local at seven sites in
biofilms_potts.jl; there is no length, height or z-extent field). So under one
isotropic pitch:

    R_physical = (N / 2) * a        L_physical = N * a        L = 2R, always.

An apparatus with an independent radius and length cannot be represented. That
is a topology fact, not a resolution one, and no pitch fixes it.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

# What the modelled domain is claimed to be. Chosen, not derived.
DOMAIN_SEMANTICS = frozenset({
    "representative_segment",   # a local axial slice of a longer apparatus
    "microvolume",              # an abstracted volume, not a whole apparatus
    "full_apparatus",           # only valid where the real apparatus has L = 2R
    "unresolved",
})

REPRESENTABILITY = frozenset({"representable", "marginal", "unresolved",
                              "not_representable"})


@dataclass(frozen=True)
class ScaleCandidate:
    candidate_id: str
    pitch_um: float
    N: int
    physical_radius_um: float
    physical_length_um: float
    species: str
    target_volume_sites: float
    entity_volume_um3: float
    volume_quantization_error: float
    minor_axis_voxels: float
    major_axis_voxels: float
    estimated_grid_memory_bytes: int
    full_apparatus_compatible: bool
    representative_segment_compatible: bool
    representability_status: str

    def as_row(self) -> dict:
        return asdict(self)


def domain_from_pitch(pitch_um: float, N: int) -> tuple[float, float]:
    """(radius, length) in micrometres. L = 2R is structural, not a choice."""
    return (N / 2.0) * pitch_um, N * pitch_um


def sites_for_volume(volume_um3: float, pitch_um: float) -> float:
    """How many lattice sites a physical volume occupies at this pitch."""
    return volume_um3 / pitch_um ** 3


def grid_memory_bytes(N: int, bytes_per_site: int = 4) -> int:
    """One Int32 lattice. The CPM carries several float64 fields alongside it,
    so treat this as a floor, not a budget."""
    return int(N) ** 3 * bytes_per_site


def enumerate_candidates(*, pitches_um, grid_sizes, entities,
                         target_volume_sites=None,
                         min_axis_voxels: float | None = None) -> list[ScaleCandidate]:
    """Every (pitch, N, entity) combination, with what each implies.

    `entities` maps a species/entity name to a dict with at least
    `volume_um3`, and optionally `minor_axis_um` / `major_axis_um` so the
    resolution of its smallest dimension can be reported.

    `target_volume_sites` is optional and, when given, is a DECLARED modelling
    choice: supplying it is what breaks the a/V_sites degeneracy, and the
    quantization error column then says what that choice costs.
    """
    out: list[ScaleCandidate] = []
    for a in pitches_um:
        for N in grid_sizes:
            radius, length = domain_from_pitch(a, N)
            for name, spec in entities.items():
                vol = float(spec["volume_um3"])
                exact_sites = sites_for_volume(vol, a)
                target = float(target_volume_sites) if target_volume_sites \
                    else exact_sites
                # error from forcing the entity into an integer-ish site target
                q_err = abs(target - exact_sites) / exact_sites if exact_sites else float("nan")

                minor_um = spec.get("minor_axis_um")
                major_um = spec.get("major_axis_um")
                minor_vox = (minor_um / a) if minor_um else float("nan")
                major_vox = (major_um / a) if major_um else float("nan")

                out.append(ScaleCandidate(
                    candidate_id=f"a{a:g}um_N{N}_{name}",
                    pitch_um=float(a), N=int(N),
                    physical_radius_um=radius, physical_length_um=length,
                    species=name,
                    target_volume_sites=target,
                    entity_volume_um3=vol,
                    volume_quantization_error=q_err,
                    minor_axis_voxels=minor_vox,
                    major_axis_voxels=major_vox,
                    estimated_grid_memory_bytes=grid_memory_bytes(N),
                    # L = 2R is forced, so a real apparatus is compatible only
                    # if its own aspect ratio happens to be 2.
                    full_apparatus_compatible=False,
                    representative_segment_compatible=True,
                    representability_status=_status(minor_vox, min_axis_voxels),
                ))
    return out


def _status(minor_axis_voxels: float, threshold: float | None) -> str:
    if threshold is None:
        # No declared acceptance threshold means no verdict is available.
        return "unresolved"
    if not np.isfinite(minor_axis_voxels):
        return "unresolved"
    if minor_axis_voxels >= threshold:
        return "representable"
    if minor_axis_voxels >= threshold / 2.0:
        return "marginal"
    return "not_representable"


def select(candidates, acceptance: dict | None):
    """Deliberately refuses without declared acceptance thresholds.

    Returns the candidates that satisfy every declared threshold — never a
    single 'the' answer, and never anything at all when the thresholds are
    unset. Choosing a pitch because it happens to fit the current N = 60 is
    exactly the failure this guards against.
    """
    if not acceptance:
        raise ValueError(
            "no acceptance thresholds declared — refusing to select a lattice "
            "pitch. Fill config/cpm_spatial_acceptance_template.toml; a pitch "
            "chosen because it fits the current N=60 is not a calibration.")
    required = ("minimum_minor_axis_voxels", "maximum_volume_quantization_error",
                "maximum_memory_bytes")
    missing = [k for k in required if acceptance.get(k) is None]
    if missing:
        raise ValueError(
            "acceptance thresholds are incomplete — refusing to select. "
            f"Unset: {', '.join(missing)}")

    return [c for c in candidates
            if np.isfinite(c.minor_axis_voxels)
            and c.minor_axis_voxels >= acceptance["minimum_minor_axis_voxels"]
            and c.volume_quantization_error <= acceptance["maximum_volume_quantization_error"]
            and c.estimated_grid_memory_bytes <= acceptance["maximum_memory_bytes"]]
