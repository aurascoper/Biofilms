"""The apparatus half of a transport config, and the checks that gate it.

The voxel binding says what one occupied CPM voxel IS and is MADE OF. That is
the biological voxel and nothing else — it does not know the cylinder it sits
in, the membrane around it, the source irradiating it, or how many histories to
run. Treating the binding as if it held the whole apparatus is why
`emit_transport_config()` could only ever produce a two-key fragment.

So the apparatus is a second, separate input, and the emitter needs both.

THREE GROUPS OF PROBLEMS, and the split is load-bearing:

    coherence_problems()   does the binding contradict itself?
    structural_problems()  is this a physically constructible configuration?
    evidence_problems()    is anyone entitled to believe these numbers?

Only the third may ever be bypassed, and only by `evidence_policy = "synthetic"`.
A synthetic run is allowed to invent values; it is NOT allowed to invent a
geometry that cannot be built, a composition that does not close, or a source
born on a lattice plane. Without that split, "synthetic" becomes a universal
escape hatch around the very contracts that exist to stop incoherent configs.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

from physical_contract import (ALLOWED_BOUNDARY_CONDITIONS,
                               ALLOWED_SOURCE_ANGULAR, ALLOWED_SOURCE_SPATIAL,
                               EVIDENCE_POLICIES, EXECUTION_CLASSES,
                               MaterialSpec, SYSTEM_PROVENANCE,
                               closed_composition_problems,
                               render_material_toml, source_placement_problems)

from .integration import VoxelBinding
from .materials.export import MaterialClass, problems_for

SYNTHETIC = "synthetic"

# Fields that only the dosimetry stage needs. Demanding a source activity for a
# transport run was the same artificial blocking the staged contract removed:
# transport reports eV per source particle and reads no activity at all.
_DOSIMETRY_ONLY = frozenset({"photons_per_second"})
_OPTIONAL = frozenset({"notes"}) | _DOSIMETRY_ONLY


@dataclass(frozen=True)
class ReferenceSystemSpec:
    """Everything a biofilm transport config needs that is NOT the voxel.

    `lattice_n` is here because the geometry checks are meaningless without it:
    whether a cylinder fits, and whether a source sits on a lattice plane, are
    both questions about `n * pitch`. It is a property of the snapshot the
    config will be run against, and a spec that does not know it cannot be
    validated before the run.
    """
    reference_system_id: str
    system_provenance: str
    evidence_policy: str
    execution_class: str
    target_calibration: bool

    lattice_n: int
    origin_cm: tuple
    cylinder_radius_cm: float
    cylinder_length_cm: float
    membrane_thickness_cm: float

    axial_bc: str
    radial_outer_bc: str

    source_spatial: str
    source_angular: str
    source_position_cm: tuple
    spectrum_energies_eV: tuple
    spectrum_probabilities: tuple

    medium: MaterialClass
    membrane: MaterialClass
    biomass_elements: dict

    batches: int
    particles: int
    seed: int
    mesh_coarsening_factor: int

    photons_per_second: float | None = None
    notes: str = ""

    @property
    def is_synthetic(self) -> bool:
        return self.evidence_policy == SYNTHETIC


def _missing(spec: ReferenceSystemSpec, stage: str) -> list[str]:
    """Every unset field, named. A blank is never read as a default."""
    out = []
    for f in fields(spec):
        if f.name in _OPTIONAL:
            continue
        if getattr(spec, f.name) is None:
            out.append(f"{f.name} is unset")
    if stage == "dosimetry" and spec.photons_per_second is None:
        out.append("photons_per_second is unset (required at the dosimetry "
                   "stage: Gy/s needs a source rate; transport does not)")
    return out


def _vocabulary(spec: ReferenceSystemSpec) -> list[str]:
    checks = [("system_provenance", spec.system_provenance, SYSTEM_PROVENANCE),
              ("evidence_policy", spec.evidence_policy, EVIDENCE_POLICIES),
              ("execution_class", spec.execution_class, EXECUTION_CLASSES),
              ("axial_bc", spec.axial_bc, ALLOWED_BOUNDARY_CONDITIONS),
              ("radial_outer_bc", spec.radial_outer_bc, ALLOWED_BOUNDARY_CONDITIONS),
              ("source_spatial", spec.source_spatial, ALLOWED_SOURCE_SPATIAL),
              ("source_angular", spec.source_angular, ALLOWED_SOURCE_ANGULAR)]
    return [f"{name} = {value!r} not in {sorted(allowed)}"
            for name, value, allowed in checks
            if value is not None and value not in allowed]


def _spectrum(spec: ReferenceSystemSpec) -> list[str]:
    """The loader does NOT normalize, so neither may the emitter: probabilities
    are a PMF, and the total emission rate belongs in photons_per_second."""
    e, p = spec.spectrum_energies_eV, spec.spectrum_probabilities
    if e is None or p is None:
        return []
    if not e or not p:
        return ["spectrum arrays are empty"]
    if len(e) != len(p):
        return [f"spectrum has {len(e)} energies but {len(p)} probabilities"]
    problems = []
    if any(v <= 0 for v in e):
        problems.append("spectrum_energies_eV must all be positive")
    if any(v < 0 for v in p):
        problems.append("spectrum_probabilities must all be non-negative")
    total = sum(p)
    if abs(total - 1.0) > 1e-6:
        problems.append(
            f"spectrum_probabilities sum to {total}, not 1 — supply a "
            "normalized PMF (p_j = Y_j / sum_k Y_k)")
    return problems


def _geometry(binding: VoxelBinding, spec: ReferenceSystemSpec) -> list[str]:
    """Congruence, membrane containment and source placement — the questions a
    synthetic system must answer exactly as a measured one does."""
    # Every input here is separately reported as a blank by `_missing`, so a
    # geometry check on a missing value would only add noise — and crash.
    needed = (binding.lattice_pitch_um, spec.lattice_n, spec.origin_cm,
              spec.cylinder_radius_cm, spec.cylinder_length_cm,
              spec.membrane_thickness_cm, spec.source_position_cm)
    if any(v is None for v in needed):
        return []
    pitch_cm = binding.lattice_pitch_um * 1e-4
    side = spec.lattice_n * pitch_cm
    problems = []

    if spec.cylinder_radius_cm > side / 2.0 * (1 + 1e-9):
        problems.append(
            f"cylinder_radius_cm = {spec.cylinder_radius_cm} exceeds the lattice "
            f"half-width {side / 2.0} (n={spec.lattice_n}, pitch={pitch_cm} cm)")
    if spec.cylinder_length_cm > side * (1 + 1e-9):
        problems.append(
            f"cylinder_length_cm = {spec.cylinder_length_cm} exceeds the lattice "
            f"extent {side}")
    # The membrane annulus lies OUTSIDE the biological cylinder but the tally
    # mesh covers only the lattice cube, so a membrane reaching past the cube is
    # partly untallied and its mass accounting is silently incomplete.
    r_out = spec.cylinder_radius_cm + spec.membrane_thickness_cm
    if r_out > side / 2.0 * (1 + 1e-9):
        problems.append(
            f"membrane outer radius {r_out} reaches outside the tallied lattice "
            f"cube (half-width {side / 2.0}) — its dose would be partly untallied")

    problems += source_placement_problems(
        spec.source_position_cm, spec.origin_cm, pitch_cm, spec.lattice_n,
        spec.cylinder_radius_cm, spec.cylinder_length_cm)
    return problems


def biomass_material(binding: VoxelBinding, spec: ReferenceSystemSpec,
                     evidence_basis: str | None = None) -> MaterialClass:
    """The occupied-voxel material: composition from the spec, density from the
    BINDING. The density belongs to the binding because it is a property of what
    a voxel is made of, not of the apparatus around it.

    Its evidence basis follows the run's policy by default, so a synthetic
    system cannot end up shipping one material labelled `declared` beside two
    labelled `synthetic` and inviting the reader to believe the odd one out.
    """
    if evidence_basis is None:
        evidence_basis = SYNTHETIC if spec.is_synthetic else "declared"
    return MaterialClass(
        name=binding.material_class,
        density_g_cm3=binding.density_g_cm3,
        elements=dict(spec.biomass_elements or {}),
        evidence_basis=evidence_basis,
        source_id=spec.reference_system_id,
        system_conditions=spec.notes or spec.reference_system_id,
        material_model_kind=binding.material_model_kind)


def structural_problems(binding: VoxelBinding, spec: ReferenceSystemSpec, *,
                        stage: str = "transport") -> list[str]:
    """Is this a physically constructible configuration? NEVER bypassed.

    A synthetic system invents its numbers; it does not get to invent a
    geometry that cannot be built or a composition that does not close.
    """
    problems = _missing(spec, stage)
    problems += _vocabulary(spec)
    problems += _spectrum(spec)

    if binding.lattice_pitch_um is None:
        problems.append("lattice_pitch_um is unset")
    elif binding.lattice_pitch_um <= 0:
        problems.append(f"lattice_pitch_um = {binding.lattice_pitch_um} is not positive")
    if binding.density_g_cm3 is None:
        problems.append("density_g_cm3 is unset")
    elif binding.density_g_cm3 <= 0:
        problems.append(f"density_g_cm3 = {binding.density_g_cm3} is not positive")

    for name, mat in (("medium", spec.medium), ("membrane", spec.membrane)):
        if mat is None:
            continue
        if mat.density_g_cm3 is None or mat.density_g_cm3 <= 0:
            problems.append(f"{name} density {mat.density_g_cm3} is not positive")
        problems += [f"{name}: {p}" for p in closed_composition_problems(mat.elements)]
    problems += [f"baseline_biomass: {p}"
                 for p in closed_composition_problems(spec.biomass_elements or {})]

    for name, value in (("batches", spec.batches), ("particles", spec.particles),
                        ("mesh_coarsening_factor", spec.mesh_coarsening_factor)):
        if value is not None and (not isinstance(value, int) or value < 1):
            problems.append(f"{name} = {value!r} must be an integer >= 1")
    if (spec.lattice_n and spec.mesh_coarsening_factor
            and spec.lattice_n % spec.mesh_coarsening_factor):
        problems.append(
            f"mesh_coarsening_factor {spec.mesh_coarsening_factor} does not "
            f"divide the lattice size {spec.lattice_n} — the meshes would not nest")

    if stage == "dosimetry" and spec.photons_per_second is not None \
            and spec.photons_per_second <= 0:
        problems.append(
            f"photons_per_second = {spec.photons_per_second} is not positive")

    if spec.is_synthetic and spec.target_calibration:
        problems.append(
            "evidence_policy = 'synthetic' with target_calibration = true — a run "
            "on invented values cannot calibrate the target system")

    problems += _geometry(binding, spec)
    return problems


def evidence_problems(binding: VoxelBinding, spec: ReferenceSystemSpec, *,
                      spatial_verdict: str, material_openmc_verdict: str
                      ) -> list[str]:
    """Is anyone entitled to believe these numbers? The ONLY bypassable group.

    `evidence_policy = "synthetic"` skips exactly this and nothing else.
    """
    problems = []
    if spatial_verdict != "READY_FOR_TIME_CALIBRATION":
        problems.append(
            f"spatial gate is {spatial_verdict}, not READY_FOR_TIME_CALIBRATION "
            "— without a selected lattice pitch there is no voxel volume, so no "
            "mass per voxel")
    if material_openmc_verdict != "READY_FOR_OPENMC_BIOFILM_TRANSPORT":
        problems.append(
            f"material OpenMC gate is {material_openmc_verdict}, not "
            "READY_FOR_OPENMC_BIOFILM_TRANSPORT — no measured density or closed "
            "composition exists to assign")

    for name, mat in (("medium", spec.medium), ("membrane", spec.membrane),
                      ("baseline_biomass", biomass_material(binding, spec))):
        if mat is not None:
            problems += [f"{name}: {p}" for p in problems_for(mat)]
    return problems


def emit_problems(binding: VoxelBinding, spec: ReferenceSystemSpec, *,
                  spatial_verdict: str, material_openmc_verdict: str,
                  stage: str = "transport") -> list[str]:
    """Every reason this configuration may not be emitted, in one list."""
    from .integration import coherence_problems

    problems = list(coherence_problems(binding))
    problems += structural_problems(binding, spec, stage=stage)
    if not spec.is_synthetic:
        problems += evidence_problems(
            binding, spec, spatial_verdict=spatial_verdict,
            material_openmc_verdict=material_openmc_verdict)
    return problems


_STAGE_NOTE = {
    "transport": ("Transport stage: reports heating in eV per source particle. "
                  "No source activity is read, so none is written."),
    "dosimetry": ("Dosimetry stage: Gy/s, so a source rate is written. "
                  "S = A(t) * sum_j Y_j — keep it separate from the PMF."),
}


def emit_transport_config(binding: VoxelBinding, spec: ReferenceSystemSpec, *,
                          spatial_verdict: str = "",
                          material_openmc_verdict: str = "",
                          stage: str = "transport") -> str:
    """Render a COMPLETE biofilm transport config, or refuse naming everything.

    Complete means loadable: `load_transport_config(..., kind="biofilm_cylinder")`
    accepts the result. The previous version emitted a pitch and one density,
    which is not a config — the loader also requires the membrane thickness, the
    cylinder, the boundaries, the source, and all three material classes.
    """
    problems = emit_problems(
        binding, spec, spatial_verdict=spatial_verdict,
        material_openmc_verdict=material_openmc_verdict, stage=stage)
    if problems:
        raise ValueError("refusing to emit a biofilm transport config:\n  - "
                         + "\n  - ".join(problems))

    pitch_cm = binding.lattice_pitch_um * 1e-4
    synthetic_banner = ("\n# SYNTHETIC: every physical value below is invented to "
                        "exercise the\n# software path. It calibrates nothing and "
                        "may not be cited as a\n# measurement of anything.\n"
                        if spec.is_synthetic else "\n")

    out = [f"# {spec.reference_system_id} — emitted by "
           f"biofilm_calibration.reference_system.",
           f"# {binding.sentence()}",
           f"# {_STAGE_NOTE.get(stage, stage)}",
           synthetic_banner.rstrip("\n") if spec.is_synthetic else "",
           f"# {spec.notes}" if spec.notes else "",
           "",
           "schema_version = 1",
           "",
           "[provenance]",
           f'reference_system_id = "{spec.reference_system_id}"',
           f'system_provenance = "{spec.system_provenance}"',
           f'evidence_policy = "{spec.evidence_policy}"',
           f'execution_class = "{spec.execution_class}"',
           f"target_calibration = {str(bool(spec.target_calibration)).lower()}",
           "",
           "[geometry]",
           f"voxel_pitch_cm = {pitch_cm!r}",
           f"origin_cm = {list(spec.origin_cm)!r}",
           f"cylinder_radius_cm = {spec.cylinder_radius_cm!r}",
           f"cylinder_length_cm = {spec.cylinder_length_cm!r}",
           f"membrane_thickness_cm = {spec.membrane_thickness_cm!r}",
           "",
           "[boundaries]",
           f'axial = "{spec.axial_bc}"',
           f'radial_outer = "{spec.radial_outer_bc}"',
           "",
           "[source]",
           f"spectrum_energies_eV = {list(spec.spectrum_energies_eV)!r}",
           f"spectrum_probabilities = {list(spec.spectrum_probabilities)!r}",
           f'spatial = "{spec.source_spatial}"',
           f'angular = "{spec.source_angular}"',
           "# Off every lattice plane: a source born on a plane is born on a "
           "surface.",
           f"position_cm = {list(spec.source_position_cm)!r}"]
    if stage == "dosimetry":
        out.append(f"photons_per_second = {spec.photons_per_second!r}")
    out += ["",
            "[transport]",
            f"batches = {spec.batches!r}",
            f"particles = {spec.particles!r}",
            f"seed = {spec.seed!r}",
            "",
            "  [transport.mesh]",
            f"  coarsening_factor = {spec.mesh_coarsening_factor!r}",
            ""]

    biomass = biomass_material(binding, spec)
    for class_name, mat in (("medium", spec.medium),
                            ("baseline_biomass", biomass),
                            ("membrane", spec.membrane)):
        out.append(render_material_toml(
            MaterialSpec(name=mat.name, density_g_cm3=mat.density_g_cm3,
                         elements=dict(mat.elements),
                         material_model_kind=mat.material_model_kind),
            class_name,
            header=f"{mat.name} — {mat.evidence_basis}, source {mat.source_id}\n"
                   f"conditions: {mat.system_conditions}"))
    return "\n".join(l for l in out if l is not None) + "\n"
