"""Physical coupling config: load, validate, fail loudly — STAGED.

Every REQUIRED key missing from the TOML is collected and reported in ONE
error — never silently defaulted (audit §5: no physical value may be
invented). tomllib is stdlib (Python >= 3.11).

Requirements are staged, because the pipeline is staged. Transport does not
depend on the source activity, the CPM clock, or the biological response
scales — `build_model()` reads none of them — so demanding them before an
OpenMC run is an artificial blocker, not a physical contract:

    transport          geometry, boundaries, source shape + normalized
                       spectrum, materials, MC controls
                       -> heating in eV/source-particle
    dosimetry          + source.photons_per_second        -> Gy/s
    cpm_feedback       + time.seconds_per_mcs, normalization.{hamiltonian,
                         melanin}_scale                   -> CPM response
    membrane_feedback  + normalization.membrane_statistic -> membrane state

Stages are cumulative and ordered, and the config types mirror that by
inheritance: a stronger stage IS-A weaker one, so anything accepting a
`TransportConfig` also accepts a full config. The stage is chosen by the
CALLER (which loader you call) — it is never declared in the TOML.

`FIELD_SPECS` and `required_keys()` are public on purpose: the provenance
ledger's schema test consumes them, so the ledger cannot drift from the
loader.
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field


class ConfigError(ValueError):
    """Raised with the complete list of missing/invalid keys."""


# Ordered and cumulative: a field required at stage i is required at every j > i.
STAGES = ("transport", "dosimetry", "cpm_feedback", "membrane_feedback")

_ALLOWED_BC = frozenset({"vacuum", "reflective"})
_ALLOWED_SPATIAL = frozenset({"line_z_axis", "point_origin"})
_ALLOWED_ANGULAR = frozenset({"isotropic"})
_ALLOWED_MEMBRANE_STAT = frozenset({"mass_weighted", "area_weighted"})


@dataclass(frozen=True)
class FieldSpec:
    """One required TOML key: where it lives, what it becomes, when it is needed."""
    section: str
    key: str
    attr: str                       # constructor kwarg on the config dataclass
    cast: str                       # "float" | "str" | "tuple" | "float_tuple"
    min_stage: str                  # first stage that requires it
    allowed: frozenset | None = None

    @property
    def dotted(self) -> str:
        return f"{self.section}.{self.key}"


FIELD_SPECS: tuple[FieldSpec, ...] = (
    # --- transport ---------------------------------------------------------
    FieldSpec("geometry", "voxel_pitch_cm", "voxel_pitch_cm", "float", "transport"),
    FieldSpec("geometry", "origin_cm", "origin_cm", "tuple", "transport"),
    FieldSpec("geometry", "cylinder_radius_cm", "cylinder_radius_cm", "float", "transport"),
    FieldSpec("geometry", "cylinder_length_cm", "cylinder_length_cm", "float", "transport"),
    FieldSpec("geometry", "membrane_thickness_cm", "membrane_thickness_cm", "float", "transport"),
    FieldSpec("boundaries", "axial", "axial_bc", "str", "transport", _ALLOWED_BC),
    FieldSpec("boundaries", "radial_outer", "radial_outer_bc", "str", "transport", _ALLOWED_BC),
    FieldSpec("source", "spectrum_energies_eV", "spectrum_energies_eV", "float_tuple", "transport"),
    FieldSpec("source", "spectrum_probabilities", "spectrum_probabilities", "float_tuple", "transport"),
    FieldSpec("source", "spatial", "source_spatial", "str", "transport", _ALLOWED_SPATIAL),
    FieldSpec("source", "angular", "source_angular", "str", "transport", _ALLOWED_ANGULAR),
    # --- dosimetry ---------------------------------------------------------
    FieldSpec("source", "photons_per_second", "photons_per_second", "float", "dosimetry"),
    # --- cpm_feedback ------------------------------------------------------
    FieldSpec("time", "seconds_per_mcs", "seconds_per_mcs", "float", "cpm_feedback"),
    FieldSpec("normalization", "hamiltonian_scale", "hamiltonian_scale", "float", "cpm_feedback"),
    FieldSpec("normalization", "melanin_scale", "melanin_scale", "float", "cpm_feedback"),
    # --- membrane_feedback -------------------------------------------------
    FieldSpec("normalization", "membrane_statistic", "membrane_statistic", "str",
              "membrane_feedback", _ALLOWED_MEMBRANE_STAT),
)


def _stage_index(stage: str) -> int:
    try:
        return STAGES.index(stage)
    except ValueError:
        raise ValueError(f"unknown stage {stage!r}; expected one of {list(STAGES)}") from None


def specs_for(stage: str) -> tuple[FieldSpec, ...]:
    """Every FieldSpec required at `stage` (cumulative)."""
    i = _stage_index(stage)
    return tuple(f for f in FIELD_SPECS if _stage_index(f.min_stage) <= i)


def required_keys(stage: str) -> frozenset[str]:
    """Dotted "section.key" paths required at `stage`. The `[materials]` class
    table is required at every stage and is not a FieldSpec — see
    `MATERIALS_REQUIRED_FROM_STAGE`."""
    return frozenset(f.dotted for f in specs_for(stage))


MATERIALS_REQUIRED_FROM_STAGE = "transport"


@dataclass(frozen=True)
class MaterialClass:
    name: str
    density_g_cm3: float
    elements: tuple  # sorted ((symbol, mass_fraction), ...) — hashable for dedup


@dataclass(frozen=True, kw_only=True)
class TransportConfig:
    """Everything `build_model()` actually reads. Sufficient for a transport
    run reporting heating in eV per source particle."""
    voxel_pitch_cm: float
    origin_cm: tuple
    cylinder_radius_cm: float
    cylinder_length_cm: float
    membrane_thickness_cm: float
    axial_bc: str
    radial_outer_bc: str
    spectrum_energies_eV: tuple
    spectrum_probabilities: tuple
    source_spatial: str
    source_angular: str
    materials: dict = field(default_factory=dict)  # name -> MaterialClass
    batches: int = 10
    particles: int = 10000
    seed: int = 1
    raw_toml: str = ""


@dataclass(frozen=True, kw_only=True)
class DoseRateConfig(TransportConfig):
    """Adds the physical emission rate, turning eV/source into Gy/s."""
    photons_per_second: float


@dataclass(frozen=True, kw_only=True)
class CPMFeedbackConfig(DoseRateConfig):
    """Adds the single time authority and the per-consumer response scales."""
    seconds_per_mcs: float
    hamiltonian_scale: float
    melanin_scale: float


@dataclass(frozen=True, kw_only=True)
class MembraneFeedbackConfig(CPMFeedbackConfig):
    """Adds the voxel-field -> membrane scalar statistic.

    STOP CONDITION (audit §6) is unchanged by this type existing: even a fully
    populated membrane_feedback config must not drive the radiodialysis model
    until the m-vs-P_eff constitutive choice is declared.
    """
    membrane_statistic: str


Config = MembraneFeedbackConfig  # back-compat alias: the full contract

_STAGE_TYPE = {
    "transport": TransportConfig,
    "dosimetry": DoseRateConfig,
    "cpm_feedback": CPMFeedbackConfig,
    "membrane_feedback": MembraneFeedbackConfig,
}


def _read(path_or_str) -> str:
    if isinstance(path_or_str, str) and "\n" in path_or_str:
        return path_or_str
    with open(path_or_str, "rb") as fh:
        return fh.read().decode()


def _cast(spec: FieldSpec, value):
    if spec.cast == "float":
        return float(value)
    if spec.cast == "str":
        return value
    if spec.cast == "tuple":
        return tuple(value)
    return tuple(float(x) for x in value)  # "float_tuple"


def _check_spectrum(energies, probabilities, problems: list[str]) -> None:
    """The loader does NOT normalize. Probabilities are a discrete PMF and must
    arrive as one; raw per-decay yields must be normalized by the caller.

    For photon line j with yield Y_j (photons per decay) and activity A:
        source.photons_per_second = A * sum_j Y_j      (total emission rate)
        spectrum_probabilities[j] = Y_j / sum_k Y_k    (the PMF)
    Keeping these separate is what stops a nuclide emitting more than one
    photon per decay from being silently mis-normalized.
    """
    if energies is None or probabilities is None:
        return
    if len(energies) == 0 or len(probabilities) == 0:
        problems.append("[source] spectrum arrays are empty (REQUIRED, unset)")
        return
    if len(energies) != len(probabilities):
        problems.append("[source] spectrum energies/probabilities length mismatch")
        return
    if any(e <= 0 for e in energies):
        problems.append("[source] spectrum_energies_eV must all be positive")
    if any(p < 0 for p in probabilities):
        problems.append("[source] spectrum_probabilities must all be non-negative")
    total = sum(probabilities)
    if total <= 0:
        problems.append("[source] spectrum_probabilities are all zero")
    elif not math.isclose(total, 1.0, rel_tol=1e-6):
        problems.append(
            f"[source] spectrum_probabilities sum to {total}, not 1 — supply a "
            "normalized PMF (p_j = Y_j / sum_k Y_k); the loader does not "
            "normalize, and the total emission rate belongs in "
            "source.photons_per_second")


def _check_materials(data: dict, problems: list[str]) -> dict:
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
    return mats


def config_from_dict(data: dict, raw: str = "", stage: str = "membrane_feedback"):
    """Validate an already-parsed config dict against `stage`'s requirements.

    Returns the dataclass for that stage. Every problem is collected and
    reported in one ConfigError — the loader never partially succeeds.
    """
    cls = _STAGE_TYPE[stage] if stage in _STAGE_TYPE else None
    if cls is None:
        raise ValueError(f"unknown stage {stage!r}; expected one of {list(STAGES)}")

    problems: list[str] = []
    values: dict[str, object] = {}
    raw_values: dict[str, object] = {}

    for spec in specs_for(stage):
        val = data.get(spec.section, {}).get(spec.key)
        raw_values[spec.dotted] = val
        if val is None:
            problems.append(f"[{spec.section}] {spec.key}  (REQUIRED, unset)")
            continue
        if spec.allowed is not None and val not in spec.allowed:
            problems.append(
                f"[{spec.section}] {spec.key} = {val!r} not in {sorted(spec.allowed)}")
            continue
        values[spec.attr] = val

    _check_spectrum(raw_values.get("source.spectrum_energies_eV"),
                    raw_values.get("source.spectrum_probabilities"), problems)
    mats = _check_materials(data, problems)

    if problems:
        raise ConfigError(
            "physical coupling config is incomplete — refusing to invent values "
            f"(audit §5). Stage: {stage}. Problems:\n  - " + "\n  - ".join(problems))

    kwargs = {spec.attr: _cast(spec, values[spec.attr]) for spec in specs_for(stage)}
    transport = data.get("transport", {})
    return cls(
        materials=mats,
        batches=int(transport.get("batches", 10)),
        particles=int(transport.get("particles", 10000)),
        seed=int(transport.get("seed", 1)),
        raw_toml=raw,
        **kwargs,
    )


def load_staged(path_or_str, stage: str):
    """Load from a file path or a raw TOML string at an explicit stage."""
    raw = _read(path_or_str)
    return config_from_dict(tomllib.loads(raw), raw, stage)


def load_transport_config(path_or_str) -> TransportConfig:
    return load_staged(path_or_str, "transport")


def load_dose_rate_config(path_or_str) -> DoseRateConfig:
    return load_staged(path_or_str, "dosimetry")


def load_cpm_feedback_config(path_or_str) -> CPMFeedbackConfig:
    return load_staged(path_or_str, "cpm_feedback")


def load_membrane_feedback_config(path_or_str) -> MembraneFeedbackConfig:
    return load_staged(path_or_str, "membrane_feedback")


def load_config(path_or_str) -> Config:
    """Back-compat alias for the strongest (full) contract."""
    return load_membrane_feedback_config(path_or_str)
