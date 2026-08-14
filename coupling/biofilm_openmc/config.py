"""Physical coupling config: load, validate, fail loudly.

Every REQUIRED key missing from the TOML is collected and reported in ONE
error — never silently defaulted (audit §5: no physical value may be
invented). tomllib is stdlib (Python >= 3.11).
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field


class ConfigError(ValueError):
    """Raised with the complete list of missing/invalid keys."""


_REQUIRED_SCALARS = [
    ("geometry", "voxel_pitch_cm"),
    ("geometry", "cylinder_radius_cm"),
    ("geometry", "cylinder_length_cm"),
    ("geometry", "membrane_thickness_cm"),
    ("boundaries", "axial"),
    ("boundaries", "radial_outer"),
    ("time", "seconds_per_mcs"),
    ("source", "photons_per_second"),
    ("source", "spatial"),
    ("source", "angular"),
    ("normalization", "hamiltonian_scale"),
    ("normalization", "melanin_scale"),
    ("normalization", "membrane_statistic"),
]
_REQUIRED_LISTS = [
    ("geometry", "origin_cm"),
    ("source", "spectrum_energies_eV"),
    ("source", "spectrum_probabilities"),
]
_ALLOWED_BC = {"vacuum", "reflective"}
_ALLOWED_SPATIAL = {"line_z_axis", "point_origin"}
_ALLOWED_ANGULAR = {"isotropic"}
_ALLOWED_MEMBRANE_STAT = {"mass_weighted", "area_weighted"}


@dataclass(frozen=True)
class MaterialClass:
    name: str
    density_g_cm3: float
    elements: tuple  # sorted ((symbol, mass_fraction), ...) — hashable for dedup


@dataclass(frozen=True)
class Config:
    voxel_pitch_cm: float
    origin_cm: tuple
    cylinder_radius_cm: float
    cylinder_length_cm: float
    membrane_thickness_cm: float
    axial_bc: str
    radial_outer_bc: str
    seconds_per_mcs: float
    photons_per_second: float
    spectrum_energies_eV: tuple
    spectrum_probabilities: tuple
    source_spatial: str
    source_angular: str
    hamiltonian_scale: float
    melanin_scale: float
    membrane_statistic: str
    materials: dict = field(default_factory=dict)  # name -> MaterialClass
    batches: int = 10
    particles: int = 10000
    seed: int = 1
    raw_toml: str = ""


def load_config(path_or_str) -> Config:
    """Load from a file path or a raw TOML string; raise ConfigError listing
    every missing/invalid REQUIRED key at once."""
    if isinstance(path_or_str, str) and "\n" in path_or_str:
        raw = path_or_str
    else:
        with open(path_or_str, "rb") as fh:
            raw = fh.read().decode()
    return config_from_dict(tomllib.loads(raw), raw)


def config_from_dict(data: dict, raw: str = "") -> Config:
    """Validate an already-parsed config dict (used by scenario overrides)."""
    problems: list[str] = []

    def get(section, key):
        val = data.get(section, {}).get(key)
        if val is None:
            problems.append(f"[{section}] {key}  (REQUIRED, unset)")
        return val

    scalars = {(s, k): get(s, k) for s, k in _REQUIRED_SCALARS}
    lists = {(s, k): get(s, k) for s, k in _REQUIRED_LISTS}

    e = lists[("source", "spectrum_energies_eV")]
    p = lists[("source", "spectrum_probabilities")]
    if e is not None and p is not None:
        if len(e) == 0 or len(p) == 0:
            problems.append("[source] spectrum arrays are empty (REQUIRED, unset)")
        elif len(e) != len(p):
            problems.append("[source] spectrum energies/probabilities length mismatch")

    for (s, k), allowed in [
        (("boundaries", "axial"), _ALLOWED_BC),
        (("boundaries", "radial_outer"), _ALLOWED_BC),
        (("source", "spatial"), _ALLOWED_SPATIAL),
        (("source", "angular"), _ALLOWED_ANGULAR),
        (("normalization", "membrane_statistic"), _ALLOWED_MEMBRANE_STAT),
    ]:
        v = scalars[(s, k)]
        if v is not None and v not in allowed:
            problems.append(f"[{s}] {k} = {v!r} not in {sorted(allowed)}")

    mats: dict[str, MaterialClass] = {}
    mat_table = data.get("materials", {})
    if not mat_table:
        problems.append("[materials] class table  (REQUIRED, unset)")
    for name, spec in mat_table.items():
        density = spec.get("density_g_cm3")
        elements = spec.get("elements", {})
        if density is None:
            problems.append(f"[materials.{name}] density_g_cm3  (REQUIRED, unset)")
        if not elements:
            problems.append(f"[materials.{name}] elements  (REQUIRED, unset)")
        else:
            total = sum(elements.values())
            if not math.isclose(total, 1.0, rel_tol=1e-6):
                problems.append(
                    f"[materials.{name}] element mass fractions sum to {total}, not 1")
        if density is not None and elements:
            mats[name] = MaterialClass(
                name, float(density), tuple(sorted(elements.items())))

    if problems:
        raise ConfigError(
            "physical coupling config is incomplete — refusing to invent values "
            "(audit §5). Problems:\n  - " + "\n  - ".join(problems))

    return Config(
        voxel_pitch_cm=float(scalars[("geometry", "voxel_pitch_cm")]),
        origin_cm=tuple(lists[("geometry", "origin_cm")]),
        cylinder_radius_cm=float(scalars[("geometry", "cylinder_radius_cm")]),
        cylinder_length_cm=float(scalars[("geometry", "cylinder_length_cm")]),
        membrane_thickness_cm=float(scalars[("geometry", "membrane_thickness_cm")]),
        axial_bc=scalars[("boundaries", "axial")],
        radial_outer_bc=scalars[("boundaries", "radial_outer")],
        seconds_per_mcs=float(scalars[("time", "seconds_per_mcs")]),
        photons_per_second=float(scalars[("source", "photons_per_second")]),
        spectrum_energies_eV=tuple(float(x) for x in e),
        spectrum_probabilities=tuple(float(x) for x in p),
        source_spatial=scalars[("source", "spatial")],
        source_angular=scalars[("source", "angular")],
        hamiltonian_scale=float(scalars[("normalization", "hamiltonian_scale")]),
        melanin_scale=float(scalars[("normalization", "melanin_scale")]),
        membrane_statistic=scalars[("normalization", "membrane_statistic")],
        materials=mats,
        batches=int(data.get("transport", {}).get("batches", 10)),
        particles=int(data.get("transport", {}).get("particles", 10000)),
        seed=int(data.get("transport", {}).get("seed", 1)),
        raw_toml=raw,
    )
