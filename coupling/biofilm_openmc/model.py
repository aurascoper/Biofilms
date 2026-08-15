"""OpenMC model builders — one per model kind.

    build_water_phantom_model(config)             # snapshot-free
    build_biofilm_cylinder_model(snapshot, config)

The kind is chosen by which builder you call, exactly as the stage is chosen by
which loader you call; it is never declared in the TOML, so a declared kind can
never contradict the builder actually invoked. Each builder asserts the config
was loaded for it, so a mismatch fails loudly rather than modeling the wrong
thing.

`import openmc` happens INSIDE functions so every schema/postprocessing path
runs without the transport stack. Geometry semantics (audit §6): the CPM's
-1 voxels mean "outside the biological domain"; the physical membrane is a
finite CSG annulus with configured thickness, never lattice voxels.

The tally mesh is no longer welded to the lattice — see `mesh.py`.
"""

from __future__ import annotations

import numpy as np

from physical_contract import source_placement_problems

from .config import (BIOFILM_CYLINDER, WATER_PHANTOM, ConfigError,
                     TransportConfig)
from .materials import BIOMASS, MEDIUM, MEMBRANE, required_classes, voxel_class_array
from .mesh import (biofilm_mesh_extent_cm, phantom_mesh_extent_cm,
                   resolve_mesh_dimension)
from .snapshot import Snapshot, to_openmc_lattice_order


def _check_kind(config: TransportConfig, expected: str) -> None:
    if config.model_kind != expected:
        raise ConfigError(
            f"config was loaded for model kind {config.model_kind!r} but the "
            f"{expected!r} builder was called. Load it with "
            f"kind={expected!r} — the two kinds require different keys and "
            "different material classes.")


def _make_material(openmc, mc):
    m = openmc.Material(name=mc.name)
    m.set_density("g/cm3", mc.density_g_cm3)
    for symbol, frac in mc.elements:
        m.add_element(symbol, frac, percent_type="wo")
    return m


def source_position_cm(config, cx, cy, z_lo, z_hi) -> tuple:
    """Where the source actually sits: the explicit `source.position_cm` when
    given, otherwise the geometric centre the builder derives.

    Public because the placement guard and the emitting side both need the same
    answer, and recomputing a centre in two places is how they come to disagree.
    """
    if config.source_position_cm is not None:
        return tuple(float(v) for v in config.source_position_cm)
    return (cx, cy, (z_lo + z_hi) / 2.0)


def _photon_source(openmc, config, cx, cy, z_lo, z_hi):
    sx, sy, sz = source_position_cm(config, cx, cy, z_lo, z_hi)
    if config.source_spatial == "line_z_axis":
        space = openmc.stats.CartesianIndependent(
            openmc.stats.Discrete([sx], [1.0]),
            openmc.stats.Discrete([sy], [1.0]),
            openmc.stats.Uniform(z_lo, z_hi))
    elif config.source_spatial == "point_origin":
        space = openmc.stats.Point((sx, sy, sz))
    else:  # pragma: no cover — config validation forbids this
        raise ValueError(f"unknown source spatial {config.source_spatial}")
    return openmc.IndependentSource(
        space=space,
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete(list(config.spectrum_energies_eV),
                                     list(config.spectrum_probabilities)),
        particle="photon")


def _settings(openmc, config, source):
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.photon_transport = True
    settings.source = source
    settings.batches = config.batches
    settings.particles = config.particles
    settings.seed = config.seed
    return settings


def _heating_tally(openmc, dimension, lower_left, upper_right):
    mesh = openmc.RegularMesh(name="voxel_mesh")
    mesh.dimension = tuple(dimension)
    mesh.lower_left = tuple(lower_left)
    mesh.upper_right = tuple(upper_right)
    heating = openmc.Tally(name="heating")
    heating.filters = [openmc.MeshFilter(mesh)]
    heating.scores = ["heating"]
    return openmc.Tallies([heating])


def check_lattice_congruence(n: int, config: TransportConfig) -> None:
    """The lattice is cubic (n,n,n) with a single isotropic pitch, so the
    modeled domain is a cube of side n*pitch. The CSG cylinder is configured
    independently — if it does not fit inside that cube, the excess is filled
    by `lattice.outer` (medium) and lies outside the tally mesh, so the run
    silently models something other than the configured apparatus.

    This is the guard that makes an apparatus of the wrong shape fail loudly.
    An extreme aspect ratio (e.g. a 3 mm bore, 1100 mm long) satisfies
    containment only at an infeasible voxel count: representing it needs
    anisotropic pitch or a declared axial segment, not a bigger cube.
    """
    side = n * config.voxel_pitch_cm
    problems = []
    if config.cylinder_radius_cm > side / 2.0 * (1 + 1e-9):
        problems.append(
            f"cylinder_radius_cm = {config.cylinder_radius_cm} exceeds the "
            f"lattice half-width {side / 2.0} (n={n}, pitch={config.voxel_pitch_cm})")
    if config.cylinder_length_cm > side * (1 + 1e-9):
        problems.append(
            f"cylinder_length_cm = {config.cylinder_length_cm} exceeds the "
            f"lattice extent {side} (n={n}, pitch={config.voxel_pitch_cm})")
    if problems:
        raise ConfigError(
            "CSG cylinder does not fit inside the voxel lattice it contains — "
            "the excess would be filled with medium and never tallied. "
            "Problems:\n  - " + "\n  - ".join(problems))


def check_source_placement(n: int, config: TransportConfig) -> None:
    """A source born exactly on a lattice plane is born on a surface, and which
    cell it starts in is then decided by floating-point tie-breaking rather than
    by the geometry. That is a real fragility, and the DEFAULT placement walks
    straight into it: the builder centres the source at `n*pitch/2`, and lattice
    planes fall at `k*pitch`, so for every EVEN `n` the centre coincides with a
    plane in all three dimensions at once.

    So this is not a property of one fixture. Rather than special-casing the
    fixture, the biofilm builder demands an explicit `source.position_cm`
    whenever the derived centre would land on a plane, and checks that whatever
    position it gets is inside the biological domain.

    The water phantom has no lattice and is not subject to any of this.
    """
    x0, y0, z0 = config.origin_cm
    pitch = config.voxel_pitch_cm
    side = n * pitch
    cx, cy = x0 + side / 2.0, y0 + side / 2.0
    z_lo, z_hi = z0, z0 + config.cylinder_length_cm
    pos = source_position_cm(config, cx, cy, z_lo, z_hi)

    # Shared with the emitting side: the builder must refuse a bad placement
    # and the emitter must not write one, from one implementation.
    problems = source_placement_problems(
        pos, config.origin_cm, pitch, n,
        config.cylinder_radius_cm, config.cylinder_length_cm)

    if problems:
        raise ConfigError(
            "source placement is degenerate — a particle born on a surface has "
            "no well-defined starting cell. Set [source] position_cm "
            "explicitly, off every lattice plane and inside the biological "
            "domain. Problems:\n  - " + "\n  - ".join(problems))


def build_water_phantom_model(config: TransportConfig):
    """Reference A0: a homogeneous water cylinder with an idealized photon
    source. No snapshot, no lattice, no membrane, no biomass.

    This is a TRANSPORT VALIDATION geometry: it exists so the chain from CSG
    through the tally to the eV/source normalization can be checked against an
    independent expectation. It is deliberately dosimetry-scale — at cell scale
    a 661.7 keV photon would traverse the whole domain without interacting, and
    the run would validate nothing.

    IDEALIZED source only: a point or axial line, with no active core, no
    encapsulation and no self-attenuation. This must not be described as a
    sealed-capsule model (reference A1, unsupported).
    """
    import openmc

    _check_kind(config, WATER_PHANTOM)
    medium = required_classes(config, (MEDIUM,))[MEDIUM]
    water = _make_material(openmc, medium)

    x0, y0, z0 = config.origin_cm
    extent = phantom_mesh_extent_cm(config)
    cx = x0 + extent[0] / 2.0
    cy = y0 + extent[1] / 2.0
    z_lo, z_hi = z0, z0 + config.cylinder_length_cm

    cyl = openmc.ZCylinder(x0=cx, y0=cy, r=config.cylinder_radius_cm,
                           boundary_type=config.radial_outer_bc)
    bottom = openmc.ZPlane(z0=z_lo, boundary_type=config.axial_bc)
    top = openmc.ZPlane(z0=z_hi, boundary_type=config.axial_bc)
    phantom = openmc.Cell(name="water_phantom", fill=water,
                          region=-cyl & +bottom & -top)
    geometry = openmc.Geometry([phantom])

    dimension = resolve_mesh_dimension(config.mesh_base_dimension,
                                       config.mesh_coarsening_factor)
    tallies = _heating_tally(openmc, dimension, (x0, y0, z0),
                             (x0 + extent[0], y0 + extent[1], z0 + extent[2]))
    source = _photon_source(openmc, config, cx, cy, z_lo, z_hi)
    return openmc.Model(geometry=geometry,
                        settings=_settings(openmc, config, source),
                        tallies=tallies)


def build_biofilm_cylinder_model(snapshot: Snapshot, config: TransportConfig):
    """Return an `openmc.Model`: fixed photon source, photon transport on,
    RegularMesh heating tally over the voxel lattice, optionally coarsened.

    Needs only a TransportConfig: the source activity, CPM clock, and
    biological response scales are not transport inputs.
    """
    import openmc

    _check_kind(config, BIOFILM_CYLINDER)
    n = snapshot.cell_id.shape[0]
    check_lattice_congruence(n, config)
    check_source_placement(n, config)
    classes = required_classes(config)

    def make_material(mc):
        return _make_material(openmc, mc)

    # Composition-level dedup: identical (density, elements) share one material
    by_spec: dict[tuple, object] = {}
    mat_for: dict[str, object] = {}
    for name, mc in classes.items():
        key = (mc.density_g_cm3, mc.elements)
        if key not in by_spec:
            by_spec[key] = make_material(mc)
        mat_for[name] = by_spec[key]

    # One universe per material class present in the voxel array
    uni_for = {}
    for name in (MEDIUM, BIOMASS):
        cell = openmc.Cell(fill=mat_for[name], name=f"class_{name}")
        uni_for[name] = openmc.Universe(cells=[cell], name=name)

    class_arr = voxel_class_array(snapshot, config)          # logical xyz
    lattice_arr = to_openmc_lattice_order(class_arr)         # (z, y-inv, x)
    pitch = config.voxel_pitch_cm
    x0, y0, z0 = config.origin_cm

    lattice = openmc.RectLattice(name="biofilm_voxels")
    lattice.lower_left = (x0, y0, z0)
    lattice.pitch = (pitch, pitch, pitch)
    lattice.universes = [
        [[uni_for[class_arr_name] for class_arr_name in row] for row in plane]
        for plane in lattice_arr
    ]
    lattice.outer = uni_for[MEDIUM]

    # CSG: interior cylinder holds the lattice; membrane is a finite annulus.
    cx = x0 + n * pitch / 2.0
    cy = y0 + n * pitch / 2.0
    r_in = config.cylinder_radius_cm
    r_out = r_in + config.membrane_thickness_cm
    z_lo, z_hi = z0, z0 + config.cylinder_length_cm

    inner_cyl = openmc.ZCylinder(x0=cx, y0=cy, r=r_in)
    outer_cyl = openmc.ZCylinder(x0=cx, y0=cy, r=r_out,
                                 boundary_type=config.radial_outer_bc)
    bottom = openmc.ZPlane(z0=z_lo, boundary_type=config.axial_bc)
    top = openmc.ZPlane(z0=z_hi, boundary_type=config.axial_bc)

    interior_cell = openmc.Cell(name="biological_domain", fill=lattice,
                                region=-inner_cyl & +bottom & -top)
    membrane_cell = openmc.Cell(name="membrane", fill=mat_for[MEMBRANE],
                                region=+inner_cyl & -outer_cyl & +bottom & -top)
    geometry = openmc.Geometry([interior_cell, membrane_cell])

    source = _photon_source(openmc, config, cx, cy, z_lo, z_hi)

    # Heating tally over the lattice cube, at the configured coarsening. The
    # extent is fixed by the geometry; only the bin count moves.
    extent = biofilm_mesh_extent_cm(config, n)
    dimension = resolve_mesh_dimension((n, n, n), config.mesh_coarsening_factor)
    tallies = _heating_tally(openmc, dimension, (x0, y0, z0),
                             (x0 + extent[0], y0 + extent[1], z0 + extent[2]))
    return openmc.Model(geometry=geometry,
                        settings=_settings(openmc, config, source),
                        tallies=tallies)


# The old name built the biofilm model; keep it pointing there so existing
# callers keep the behaviour they had.
build_model = build_biofilm_cylinder_model
