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

from physical_contract import (ALLOWED_BOUNDARY_CONDITIONS,
                               ALLOWED_SOURCE_ANGULAR, ALLOWED_SOURCE_SPATIAL,
                               EVIDENCE_POLICIES, EXECUTION_CLASSES,
                               SYSTEM_PROVENANCE)


class ConfigError(ValueError):
    """Raised with the complete list of missing/invalid keys."""


# Ordered and cumulative: a field required at stage i is required at every j > i.
STAGES = ("transport", "dosimetry", "cpm_feedback", "membrane_feedback")

# Shared with calibration/ via `physical_contract`, which owns them: the
# emitting side must validate against the same vocabulary the loader enforces,
# and it cannot reach into another package's private names to do it. Aliased
# here so no call site in this package moves.
_ALLOWED_BC = ALLOWED_BOUNDARY_CONDITIONS
_ALLOWED_SPATIAL = ALLOWED_SOURCE_SPATIAL
_ALLOWED_ANGULAR = ALLOWED_SOURCE_ANGULAR
# Not shared: nothing on the calibration side emits a membrane statistic, and
# the m-vs-P_eff stop condition means nothing should until it is resolved.
_ALLOWED_MEMBRANE_STAT = frozenset({"mass_weighted", "area_weighted"})

# Which model a config is loaded for. Like the stage, this is chosen by the
# CALLER (which builder/loader you call) and never declared in the TOML, so a
# declared kind can never contradict the builder that was actually invoked.
WATER_PHANTOM = "water_phantom"
BIOFILM_CYLINDER = "biofilm_cylinder"
MODEL_KINDS = (WATER_PHANTOM, BIOFILM_CYLINDER)
_ALL_KINDS = frozenset(MODEL_KINDS)


@dataclass(frozen=True)
class FieldSpec:
    """One required TOML key: where it lives, what it becomes, when it is needed.

    `section` may be dotted (`"transport.mesh"`) to address a TOML subtable.
    `kinds` restricts the key to particular model kinds — a water phantom has
    no CPM lattice and no membrane, so demanding their parameters would be the
    same artificial blocking the staged contract removed.
    """
    section: str
    key: str
    attr: str                       # constructor kwarg on the config dataclass
    cast: str                       # "float" | "str" | "tuple" | "float_tuple" | "int_tuple"
    min_stage: str                  # first stage that requires it
    allowed: frozenset | None = None
    kinds: frozenset = _ALL_KINDS

    @property
    def dotted(self) -> str:
        return f"{self.section}.{self.key}"


_BIOFILM_ONLY = frozenset({BIOFILM_CYLINDER})
_PHANTOM_ONLY = frozenset({WATER_PHANTOM})

_ALLOWED_TARGET_CALIBRATION = frozenset({True, False})

FIELD_SPECS: tuple[FieldSpec, ...] = (
    # --- transport ---------------------------------------------------------
    # WHAT SYSTEM THIS IS. Required, not optional: a config that cannot say what
    # it represents has no business building a model, and `target_calibration =
    # false` alone cannot distinguish an A0 benchmark from a public surrogate,
    # a synthetic fixture, an uncalibrated Reference D, or a sensitivity case.
    # These four axes are independent — an engineered composite can hold
    # measured, certified, derived and declared components at once, so
    # "may this run use unmeasured values" cannot be a provenance value.
    FieldSpec("provenance", "reference_system_id", "reference_system_id", "str",
              "transport"),
    FieldSpec("provenance", "system_provenance", "system_provenance", "str",
              "transport", SYSTEM_PROVENANCE),
    FieldSpec("provenance", "evidence_policy", "evidence_policy", "str",
              "transport", EVIDENCE_POLICIES),
    FieldSpec("provenance", "execution_class", "execution_class", "str",
              "transport", EXECUTION_CLASSES),
    FieldSpec("provenance", "target_calibration", "target_calibration", "bool",
              "transport", _ALLOWED_TARGET_CALIBRATION),
    # The CPM lattice pitch: a biofilm-model parameter. A water phantom has no
    # lattice, and its tally resolution comes from transport.mesh instead.
    FieldSpec("geometry", "voxel_pitch_cm", "voxel_pitch_cm", "float", "transport",
              kinds=_BIOFILM_ONLY),
    FieldSpec("geometry", "origin_cm", "origin_cm", "tuple", "transport"),
    FieldSpec("geometry", "cylinder_radius_cm", "cylinder_radius_cm", "float", "transport"),
    FieldSpec("geometry", "cylinder_length_cm", "cylinder_length_cm", "float", "transport"),
    FieldSpec("geometry", "membrane_thickness_cm", "membrane_thickness_cm", "float", "transport",
              kinds=_BIOFILM_ONLY),
    FieldSpec("boundaries", "axial", "axial_bc", "str", "transport", _ALLOWED_BC),
    FieldSpec("boundaries", "radial_outer", "radial_outer_bc", "str", "transport", _ALLOWED_BC),
    FieldSpec("source", "spectrum_energies_eV", "spectrum_energies_eV", "float_tuple", "transport"),
    FieldSpec("source", "spectrum_probabilities", "spectrum_probabilities", "float_tuple", "transport"),
    FieldSpec("source", "spatial", "source_spatial", "str", "transport", _ALLOWED_SPATIAL),
    FieldSpec("source", "angular", "source_angular", "str", "transport", _ALLOWED_ANGULAR),
    # Absolute tally resolution. Required only where there is no lattice to
    # inherit one from; the biofilm model derives its base from the snapshot.
    FieldSpec("transport.mesh", "base_dimension", "mesh_base_dimension", "int_tuple",
              "transport", kinds=_PHANTOM_ONLY),
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


@dataclass(frozen=True)
class OptionalSpec:
    """A key with a safe default, so it never blocks a run — but still a real
    config key, and the numerical ones are exactly what a resolution/history
    sweep varies, so the provenance ledger must be able to name them."""
    section: str
    key: str
    attr: str
    default: object

    @property
    def dotted(self) -> str:
        return f"{self.section}.{self.key}"


OPTIONAL_SPECS: tuple[OptionalSpec, ...] = (
    OptionalSpec("transport", "batches", "batches", 10),
    OptionalSpec("transport", "particles", "particles", 10000),
    OptionalSpec("transport", "seed", "seed", 1),
    OptionalSpec("transport.mesh", "coarsening_factor", "mesh_coarsening_factor", 1),
    # Absent means "let the builder centre it", which the biofilm builder
    # refuses for an even lattice — so this is optional at the loader and
    # required in practice for that kind.
    OptionalSpec("source", "position_cm", "source_position_cm", None),
)


def known_keys() -> frozenset[str]:
    """Every config key the loader recognises, required or defaulted."""
    return frozenset(f.dotted for f in FIELD_SPECS) | \
        frozenset(o.dotted for o in OPTIONAL_SPECS)


def _stage_index(stage: str) -> int:
    try:
        return STAGES.index(stage)
    except ValueError:
        raise ValueError(f"unknown stage {stage!r}; expected one of {list(STAGES)}") from None


def _check_kind(kind: str) -> str:
    if kind not in _ALL_KINDS:
        raise ValueError(f"unknown model kind {kind!r}; expected one of {list(MODEL_KINDS)}")
    return kind


def specs_for(stage: str, kind: str = BIOFILM_CYLINDER) -> tuple[FieldSpec, ...]:
    """Every FieldSpec required at `stage` for `kind` (stages are cumulative)."""
    i = _stage_index(stage)
    _check_kind(kind)
    return tuple(f for f in FIELD_SPECS
                 if _stage_index(f.min_stage) <= i and kind in f.kinds)


def required_keys(stage: str, kind: str = BIOFILM_CYLINDER) -> frozenset[str]:
    """Dotted "section.key" paths required at `stage` for `kind`. The
    `[materials]` class table is required at every stage and is not a FieldSpec
    — see `MATERIALS_REQUIRED_FROM_STAGE` and `required_material_classes`."""
    return frozenset(f.dotted for f in specs_for(stage, kind))


MATERIALS_REQUIRED_FROM_STAGE = "transport"


@dataclass(frozen=True)
class MaterialClass:
    name: str
    density_g_cm3: float
    elements: tuple  # sorted ((symbol, mass_fraction), ...) — hashable for dedup


@dataclass(frozen=True, kw_only=True)
class TransportConfig:
    """Everything a transport builder actually reads. Sufficient for a run
    reporting heating in eV per source particle.

    `model_kind` records which builder this config was loaded FOR, so a builder
    can refuse a config prepared for the other model instead of silently
    modeling the wrong thing. The kind-specific fields default to None and are
    guaranteed present by the loader whenever the kind requires them.
    """
    reference_system_id: str
    system_provenance: str
    evidence_policy: str
    execution_class: str
    target_calibration: bool
    origin_cm: tuple
    cylinder_radius_cm: float
    cylinder_length_cm: float
    axial_bc: str
    radial_outer_bc: str
    spectrum_energies_eV: tuple
    spectrum_probabilities: tuple
    source_spatial: str
    source_angular: str
    model_kind: str = BIOFILM_CYLINDER
    # Explicit source placement. Optional because the builders derive a centre
    # when it is absent — but that derived centre lands exactly on internal
    # lattice planes for every EVEN lattice size, so the biofilm builder
    # requires it. See `check_source_placement`.
    source_position_cm: tuple | None = None
    # biofilm_cylinder only
    voxel_pitch_cm: float | None = None
    membrane_thickness_cm: float | None = None
    # water_phantom only
    mesh_base_dimension: tuple | None = None
    # tally coarsening: applies to both kinds; 1 means "as fine as the base"
    mesh_coarsening_factor: int = 1
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
    if spec.cast in ("str", "bool"):
        return value
    if spec.cast == "tuple":
        return tuple(value)
    if spec.cast == "int_tuple":
        return tuple(int(x) for x in value)
    return tuple(float(x) for x in value)  # "float_tuple"


def _lookup(data: dict, section: str, key: str):
    """Read `key` from a possibly-dotted TOML section ("transport.mesh")."""
    node = data
    for part in section.split("."):
        node = node.get(part, {})
        if not isinstance(node, dict):
            return None
    return node.get(key)


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


def required_material_classes(kind: str = BIOFILM_CYLINDER) -> tuple[str, ...]:
    """Material classes a model kind cannot be built without.

    A water phantom needs only the medium. Demanding biomass and a membrane for
    it was the water-phantom gap: it would have forced fabricated entries into a
    config just to satisfy the biofilm builder
    (`docs/physical_reference_system.md` §5).
    """
    _check_kind(kind)
    return ("medium",) if kind == WATER_PHANTOM else (
        "medium", "baseline_biomass", "membrane")


def config_from_dict(data: dict, raw: str = "", stage: str = "membrane_feedback",
                     kind: str = BIOFILM_CYLINDER):
    """Validate an already-parsed config dict against `stage` and `kind`.

    Returns the dataclass for that stage. Every problem is collected and
    reported in one ConfigError — the loader never partially succeeds.
    """
    cls = _STAGE_TYPE[stage] if stage in _STAGE_TYPE else None
    if cls is None:
        raise ValueError(f"unknown stage {stage!r}; expected one of {list(STAGES)}")
    _check_kind(kind)

    problems: list[str] = []
    values: dict[str, object] = {}
    raw_values: dict[str, object] = {}

    for spec in specs_for(stage, kind):
        val = _lookup(data, spec.section, spec.key)
        raw_values[spec.dotted] = val
        if val is None:
            problems.append(f"[{spec.section}] {spec.key}  (REQUIRED, unset)")
            continue
        if spec.cast == "bool" and not isinstance(val, bool):
            # `target_calibration = 1` is truthy and would pass an `in {True,
            # False}` check, because 1 == True. A provenance flag is not a
            # number.
            problems.append(
                f"[{spec.section}] {spec.key} = {val!r} must be a TOML boolean")
            continue
        if spec.allowed is not None and val not in spec.allowed:
            problems.append(
                f"[{spec.section}] {spec.key} = {val!r} not in "
                f"{sorted(spec.allowed, key=repr)}")
            continue
        values[spec.attr] = val

    _check_spectrum(raw_values.get("source.spectrum_energies_eV"),
                    raw_values.get("source.spectrum_probabilities"), problems)
    mats = _check_materials(data, problems)
    missing_classes = [n for n in required_material_classes(kind) if n not in mats]
    if mats and missing_classes:
        problems.append(
            f"[materials] {kind} requires the class(es): {', '.join(missing_classes)}")

    coarsening = data.get("transport", {}).get("mesh", {}).get("coarsening_factor", 1)
    if not isinstance(coarsening, int) or isinstance(coarsening, bool) or coarsening < 1:
        problems.append(
            f"[transport.mesh] coarsening_factor = {coarsening!r} must be an "
            "integer >= 1 (1 means the tally mesh is as fine as its base)")

    position = _lookup(data, "source", "position_cm")
    if position is not None and (not isinstance(position, (list, tuple))
                                 or len(position) != 3):
        problems.append(
            f"[source] position_cm = {position!r} must be three coordinates")

    base_dim = values.get("mesh_base_dimension")
    if base_dim is not None and (len(base_dim) != 3 or any(int(d) < 1 for d in base_dim)):
        problems.append(
            f"[transport.mesh] base_dimension = {list(base_dim)} must be three "
            "positive integers")

    if problems:
        raise ConfigError(
            "physical coupling config is incomplete — refusing to invent values "
            f"(audit §5). Stage: {stage}, model kind: {kind}. Problems:\n  - "
            + "\n  - ".join(problems))

    kwargs = {spec.attr: _cast(spec, values[spec.attr])
              for spec in specs_for(stage, kind)}
    transport = data.get("transport", {})
    position = _lookup(data, "source", "position_cm")
    return cls(
        model_kind=kind,
        source_position_cm=None if position is None else tuple(float(v) for v in position),
        mesh_coarsening_factor=coarsening,
        materials=mats,
        batches=int(transport.get("batches", 10)),
        particles=int(transport.get("particles", 10000)),
        seed=int(transport.get("seed", 1)),
        raw_toml=raw,
        **kwargs,
    )


def load_staged(path_or_str, stage: str, kind: str = BIOFILM_CYLINDER):
    """Load from a file path or a raw TOML string at an explicit stage/kind."""
    raw = _read(path_or_str)
    return config_from_dict(tomllib.loads(raw), raw, stage, kind)


def load_transport_config(path_or_str, kind: str = BIOFILM_CYLINDER) -> TransportConfig:
    return load_staged(path_or_str, "transport", kind)


def load_dose_rate_config(path_or_str, kind: str = BIOFILM_CYLINDER) -> DoseRateConfig:
    return load_staged(path_or_str, "dosimetry", kind)


def load_cpm_feedback_config(path_or_str, kind: str = BIOFILM_CYLINDER) -> CPMFeedbackConfig:
    return load_staged(path_or_str, "cpm_feedback", kind)


def load_membrane_feedback_config(path_or_str,
                                  kind: str = BIOFILM_CYLINDER) -> MembraneFeedbackConfig:
    return load_staged(path_or_str, "membrane_feedback", kind)


def load_config(path_or_str) -> Config:
    """Back-compat alias for the strongest (full) contract, biofilm model."""
    return load_membrane_feedback_config(path_or_str)
