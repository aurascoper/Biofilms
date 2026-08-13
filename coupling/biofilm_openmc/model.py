"""OpenMC model builder: CSG cylinder + membrane annulus + voxel RectLattice.

`import openmc` happens INSIDE functions so every schema/postprocessing path
runs without the transport stack. Geometry semantics (audit §6): the CPM's
-1 voxels mean "outside the biological domain"; the physical membrane is a
finite CSG annulus with configured thickness, never lattice voxels.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .materials import BIOMASS, MEDIUM, MEMBRANE, required_classes, voxel_class_array
from .snapshot import Snapshot, to_openmc_lattice_order


def build_model(snapshot: Snapshot, config: Config):
    """Return an `openmc.Model`: fixed photon source, photon transport on,
    RegularMesh heating tally over the voxel lattice."""
    import openmc

    classes = required_classes(config)

    def make_material(mc):
        m = openmc.Material(name=mc.name)
        m.set_density("g/cm3", mc.density_g_cm3)
        for symbol, frac in mc.elements:
            m.add_element(symbol, frac, percent_type="wo")
        return m

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
    n = snapshot.cell_id.shape[0]
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

    # Fixed photon source
    if config.source_spatial == "line_z_axis":
        space = openmc.stats.CartesianIndependent(
            openmc.stats.Discrete([cx], [1.0]),
            openmc.stats.Discrete([cy], [1.0]),
            openmc.stats.Uniform(z_lo, z_hi))
    elif config.source_spatial == "point_origin":
        space = openmc.stats.Point((cx, cy, (z_lo + z_hi) / 2.0))
    else:  # pragma: no cover — config validation forbids this
        raise ValueError(f"unknown source spatial {config.source_spatial}")

    source = openmc.IndependentSource(
        space=space,
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete(list(config.spectrum_energies_eV),
                                     list(config.spectrum_probabilities)),
        particle="photon")

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.photon_transport = True
    settings.source = source
    settings.batches = config.batches
    settings.particles = config.particles
    settings.seed = config.seed

    # Heating tally on a mesh congruent with the voxel lattice
    mesh = openmc.RegularMesh(name="voxel_mesh")
    mesh.dimension = (n, n, n)
    mesh.lower_left = (x0, y0, z0)
    mesh.upper_right = (x0 + n * pitch, y0 + n * pitch, z0 + n * pitch)

    heating = openmc.Tally(name="heating")
    heating.filters = [openmc.MeshFilter(mesh)]
    heating.scores = ["heating"]

    tallies = openmc.Tallies([heating])
    return openmc.Model(geometry=geometry, settings=settings, tallies=tallies)
