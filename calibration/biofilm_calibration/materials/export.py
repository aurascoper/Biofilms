"""Export to an OpenMC material class — and the refusals that guard it.

The output shape is the one `coupling/biofilm_openmc/config.py` already
requires: a density in g/cm3 plus elemental MASS fractions summing to 1. That
loader will reject a bad sum on its own; what it cannot possibly know is
whether the number it was handed was measured, estimated, or a proxy for a
material nobody has characterised. That judgement lives here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .basis_conversion import BasisError

# Evidence that may back an exported, calibrated material.
EXPORTABLE_EVIDENCE = frozenset({"direct_measurement", "manufacturer_datasheet",
                                 "primary_literature", "derived"})
# Evidence that may NOT, however plausible the numbers look.
BLOCKED_EVIDENCE = frozenset({"proxy", "derived_proxy", "synthetic", "declared",
                              "unresolved"})


@dataclass(frozen=True)
class MaterialClass:
    """A candidate OpenMC material class."""
    name: str
    density_g_cm3: float
    elements: dict
    evidence_basis: str
    source_id: str
    system_conditions: str
    material_model_kind: str = "hydrated_effective_medium"


def problems_for(material: MaterialClass) -> list[str]:
    """Everything that stops this material being exported as calibrated."""
    problems: list[str] = []

    if material.evidence_basis in BLOCKED_EVIDENCE:
        problems.append(
            f"evidence_basis={material.evidence_basis!r} cannot be exported as "
            "a calibrated material — a proxy or declared composition is a "
            "placeholder, and shipping it would launder it into a measurement")
    elif material.evidence_basis not in EXPORTABLE_EVIDENCE:
        problems.append(f"unknown evidence_basis {material.evidence_basis!r}")

    if not material.source_id.strip():
        problems.append("no source_id — an exported material must be attributable")
    if not material.system_conditions.strip():
        problems.append(
            "no system_conditions — a composition without its growth medium, "
            "temperature and hydration state cannot be judged applicable")

    if material.density_g_cm3 is None or material.density_g_cm3 <= 0:
        problems.append(f"density {material.density_g_cm3} is not positive")

    if not material.elements:
        problems.append("no elemental composition")
    else:
        for el, frac in material.elements.items():
            if frac < 0:
                problems.append(f"negative mass fraction for {el}")
        total = sum(material.elements.values())
        if not math.isclose(total, 1.0, rel_tol=1e-6):
            problems.append(
                f"elemental mass fractions sum to {total}, not 1 — the "
                "coupling loader requires a closed composition")

    if material.material_model_kind != "hydrated_effective_medium":
        problems.append(
            f"material_model_kind={material.material_model_kind!r} is not "
            "compatible with the current mapper, which assigns one material "
            "class per occupied voxel; explicit components would need sub-voxel "
            "homogenization that does not exist")
    return problems


def to_config_fragment(material: MaterialClass, class_name: str) -> str:
    """Render a TOML fragment for `[materials.<class_name>]`, or refuse."""
    problems = problems_for(material)
    if problems:
        raise BasisError(
            f"refusing to export {material.name!r} as a calibrated material:\n  - "
            + "\n  - ".join(problems))

    lines = [f"# {material.name} — {material.evidence_basis}, source "
             f"{material.source_id}",
             f"# conditions: {material.system_conditions}",
             f"[materials.{class_name}]",
             f"density_g_cm3 = {material.density_g_cm3!r}",
             f"  [materials.{class_name}.elements]"]
    for el, frac in sorted(material.elements.items()):
        lines.append(f"  {el} = {frac!r}")
    return "\n".join(lines) + "\n"
