"""Pitch selection from coarse-grained biomass fields, and its refusals.

The rule:

    Choose the COARSEST pitch that preserves every declared calibration
    observable within its tolerance, on training stacks AND on held-out
    stacks.

Coarsest, because a finer lattice is always more faithful and always more
expensive; the calibration question is how coarse you may go. Held-out,
because a pitch fitted and validated on the same stacks has not been validated.

This module selects nothing without declared tolerances and refuses without a
held-out set, mirroring `scale_candidates.select()`.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .field import FieldError, apply_occupancy, coarse_grain
from .structure import summarise

# The observables a candidate pitch must preserve, and the tolerance key each
# is judged by. Kept in step with config/cpm_spatial_acceptance_template.toml.
OBSERVABLE_TOLERANCES = {
    "biovolume_fraction": "maximum_biovolume_fraction_error",
    "porosity": "maximum_porosity_error",
    "specific_interface_area_per_um": "maximum_interface_area_error",
    "correlation_length_um": "maximum_correlation_length_error",
    "thickness_mean_um": "maximum_thickness_error",
    "component_size_q50_um3": "maximum_component_size_error",
}


@dataclass(frozen=True)
class PitchEvaluation:
    factor: int
    pitch_um: float
    sample_id: str
    held_out: bool
    passed: bool
    errors: dict
    reference: dict
    coarse: dict

    def as_row(self) -> dict:
        d = asdict(self)
        d["errors"] = dict(sorted(d["errors"].items()))
        return d


def _relative_error(reference: float, got: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(got):
        return float("nan")
    if reference == 0.0:
        return 0.0 if got == 0.0 else float("inf")
    return abs(got - reference) / abs(reference)


def evaluate_pitch(B: np.ndarray, base_voxel_um, factor: int, *,
                   sample_id: str, held_out: bool, occupancy_mapping: str,
                   tolerances: dict, tau: float | None = None,
                   seed: int | None = None) -> PitchEvaluation:
    """Compare one candidate coarsening against the full-resolution field."""
    base = np.asarray(base_voxel_um, dtype=float)
    if base.shape == ():
        base = np.repeat(base, 3)

    ref_mask = apply_occupancy(B, occupancy_mapping, tau=tau, seed=seed)
    reference = summarise(B, ref_mask, base).as_dict()

    phi = coarse_grain(B, factor)
    coarse_voxel = base * factor
    mask = apply_occupancy(phi, occupancy_mapping, tau=tau, seed=seed)
    coarse = summarise(phi, mask, coarse_voxel).as_dict()

    errors, passed = {}, True
    for observable, tol_key in OBSERVABLE_TOLERANCES.items():
        err = _relative_error(reference[observable], coarse[observable])
        errors[observable] = err
        tol = tolerances.get(tol_key)
        if tol is None:
            passed = False           # an undeclared tolerance cannot be passed
        elif not np.isfinite(err) or err > float(tol):
            passed = False

    return PitchEvaluation(factor=int(factor), pitch_um=float(coarse_voxel[0]),
                           sample_id=sample_id, held_out=bool(held_out),
                           passed=passed, errors=errors,
                           reference=reference, coarse=coarse)


def select_pitch(evaluations, tolerances: dict):
    """The coarsest factor passing on every stack, training and held-out.

    Refuses without declared tolerances and without at least one held-out
    stack. Returns the winning factor, or None when nothing qualifies — which
    is a real answer, not an error.
    """
    if not tolerances:
        raise FieldError(
            "no acceptance tolerances declared — refusing to select a pitch. "
            "Fill config/cpm_spatial_acceptance_template.toml; a pitch chosen "
            "because it looks reasonable is not a calibration.")
    missing = [k for k in OBSERVABLE_TOLERANCES.values()
               if tolerances.get(k) is None]
    if missing:
        raise FieldError(
            "acceptance tolerances are incomplete — refusing to select. "
            f"Unset: {', '.join(sorted(missing))}")

    evaluations = list(evaluations)
    if not evaluations:
        raise FieldError("no pitch evaluations supplied")
    if not any(e.held_out for e in evaluations):
        raise FieldError(
            "no held-out stack among the evaluations — refusing to select. A "
            "pitch fitted and validated on the same stacks has not been "
            "validated.")

    by_factor: dict[int, list] = {}
    for e in evaluations:
        by_factor.setdefault(e.factor, []).append(e)

    qualifying = [f for f, group in by_factor.items()
                  if all(e.passed for e in group)
                  and any(e.held_out for e in group)]
    return max(qualifying) if qualifying else None
