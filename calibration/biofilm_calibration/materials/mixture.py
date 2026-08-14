"""Component mixing with mass closure, for a hydrated effective medium.

The hydrated biofilm is represented as ONE OpenMC material whose elemental
composition is the mass-weighted blend of its components:

    w_i_wet = sum_k w_k * w_i|k        with sum_k w_k = 1

This is `material_model_kind = hydrated_effective_medium`, the only kind
compatible with the present mapper, which assigns one material class per
occupied voxel (`coupling/biofilm_openmc/materials.py`). Representing dry solid
and interstitial water as separate phases would need sub-voxel homogenization
or additional classes, and neither exists.

The rejections here are the point. Each corresponds to a way a plausible-looking
composition is wrong:

  * negative or out-of-range component fractions;
  * component fractions not summing to 1 (a forgotten component reads as a
    silent renormalization);
  * elemental fractions within a component not summing to 1;
  * the same element counted in both the ash component and an explicit metal
    loading — easy to create when ash comes from combustion and the metal is
    also assayed directly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .basis_conversion import BasisError

TOLERANCE = 1e-6


@dataclass(frozen=True)
class Component:
    """One phase of the hydrated biofilm."""
    name: str
    mass_fraction: float                 # of the wet bulk
    elements: dict                       # element -> mass fraction within this component
    density_g_cm3: float | None = None
    contributes_metals: frozenset = field(default_factory=frozenset)

    def __post_init__(self):
        if not 0.0 <= self.mass_fraction <= 1.0:
            raise BasisError(
                f"component {self.name!r} mass fraction {self.mass_fraction} "
                "is outside [0, 1]")
        if not self.elements:
            raise BasisError(f"component {self.name!r} has no elemental composition")
        for el, frac in self.elements.items():
            if frac < 0:
                raise BasisError(
                    f"component {self.name!r}: negative mass fraction for {el}")
        total = sum(self.elements.values())
        if not math.isclose(total, 1.0, rel_tol=1e-6):
            raise BasisError(
                f"component {self.name!r}: elemental mass fractions sum to "
                f"{total}, not 1")


def _check_metal_double_count(components) -> list[str]:
    """An element claimed by more than one component's `contributes_metals` is
    being counted twice — the classic ash-plus-assay double count."""
    problems = []
    seen: dict[str, str] = {}
    for c in components:
        for el in c.contributes_metals:
            if el in seen:
                problems.append(
                    f"{el} is claimed by both {seen[el]!r} and {c.name!r} — "
                    "metal determined by combustion into ash AND assayed "
                    "separately is counted twice")
            else:
                seen[el] = c.name
    return problems


def blend(components) -> dict:
    """Elemental mass fractions of the mixture. Raises on any closure failure."""
    components = list(components)
    problems: list[str] = []

    total_mass = sum(c.mass_fraction for c in components)
    if not math.isclose(total_mass, 1.0, rel_tol=TOLERANCE):
        problems.append(
            f"component mass fractions sum to {total_mass}, not 1 — a missing "
            "component would otherwise be absorbed as a silent renormalization")

    problems.extend(_check_metal_double_count(components))
    if problems:
        raise BasisError("mixture does not close:\n  - " + "\n  - ".join(problems))

    out: dict[str, float] = {}
    for c in components:
        for el, frac in c.elements.items():
            out[el] = out.get(el, 0.0) + c.mass_fraction * frac

    blended = sum(out.values())
    if not math.isclose(blended, 1.0, rel_tol=1e-6):     # arithmetic guard
        raise BasisError(f"blended elemental fractions sum to {blended}, not 1")
    return dict(sorted(out.items()))


def ideal_mixture_density(components) -> tuple[float, str]:
    """Volume-additive density estimate.

    Returns (value, evidence_basis) with evidence_basis fixed at
    'derived_proxy'. An ideal-mixture density is NOT a measured density: real
    biofilm has excluded volume and pore structure that volume additivity does
    not model. `export` refuses to ship a proxy as calibrated.
    """
    components = list(components)
    missing = [c.name for c in components if c.density_g_cm3 is None]
    if missing:
        raise BasisError(
            "cannot estimate a mixture density without every component "
            f"density; missing: {', '.join(missing)}")
    inv = sum(c.mass_fraction / c.density_g_cm3 for c in components)
    if inv <= 0:
        raise BasisError("non-positive specific volume")
    return 1.0 / inv, "derived_proxy"
