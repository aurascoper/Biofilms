"""Material-class mapping: config-defined table, composition-deduplicated.

Genealogy/cell labels NEVER become materials: the voxel class array is a
function of occupancy + config only, and unique materials are bounded by the
config's class table regardless of how many cells/lineages exist.
"""

from __future__ import annotations

import numpy as np

from .config import ConfigError, MaterialClass, TransportConfig
from .snapshot import Snapshot

# Voxel class codes (indices into the per-snapshot class list)
MEDIUM = "medium"
BIOMASS = "baseline_biomass"
MEMBRANE = "membrane"  # CSG region, never a lattice voxel class


def required_classes(config: TransportConfig,
                     names: tuple = (MEDIUM, BIOMASS, MEMBRANE)) -> dict[str, MaterialClass]:
    """The classes a caller needs; loud error if the config lacks one.

    `names` is a parameter because a water phantom needs only the medium.
    Hard-coding all three here was the water-phantom gap: it forced fabricated
    biomass and membrane entries into any config that just wanted water
    (`docs/physical_reference_system.md` §5).
    """
    missing = [n for n in names if n not in config.materials]
    if missing:
        raise ConfigError(
            "material class table lacks required classes: " + ", ".join(missing))
    return {n: config.materials[n] for n in names}


def voxel_class_array(snapshot: Snapshot, config: TransportConfig) -> np.ndarray:
    """Logical (x,y,z) array of class NAMES (object dtype indexing into the
    config table). Occupied voxels -> baseline_biomass; everything else ->
    medium. Voxels outside the biological domain (cell_id == -1) get medium
    as filler: they are clipped by the CSG cylinder and their lattice
    material is never seen by transport (the membrane is CSG, audit §6).
    """
    required_classes(config)
    classes = np.full(snapshot.cell_id.shape, MEDIUM, dtype=object)
    classes[snapshot.cell_id > 0] = BIOMASS
    return classes


def unique_material_specs(snapshot: Snapshot, config: TransportConfig) -> list[MaterialClass]:
    """Deduplicated (by frozen composition spec) materials actually present in
    the voxel array plus the CSG membrane. Bounded by the config table size —
    NEVER by cell or lineage count."""
    present = set(np.unique(voxel_class_array(snapshot, config)))
    present.add(MEMBRANE)
    # dedup by value: identical (density, elements) collapse to one spec
    specs: dict[tuple, MaterialClass] = {}
    for name in sorted(present):
        mc = config.materials[name]
        specs.setdefault((mc.density_g_cm3, mc.elements), mc)
    return list(specs.values())


def phantom_mass_kg(config: TransportConfig, dimension, extent_cm) -> np.ndarray:
    """Uniform medium masses (kg) over the tally mesh, logical (x,y,z).

    Snapshot-free by construction, and the bin volume comes from the PHYSICAL
    extent divided by the bin count — never from `voxel_pitch_cm`, which is the
    lattice pitch and would rescale every dose by the coarsening factor cubed.

    ponytail: full rectangular-bin volume, so bins in the corners of the
    circumscribing cube are counted as water although the CSG cylinder leaves
    them void. Their heating is zero, so dose there is zero and unaffected;
    only mass-weighted aggregates carry the bias. Exact volumes need openmc
    RegularMesh.material_volumes() — the same approximation, and the same
    upgrade path, as voxel_mass_kg below.
    """
    from .mesh import mesh_bin_volume_cm3

    medium = required_classes(config, (MEDIUM,))[MEDIUM]
    bin_volume_cm3 = mesh_bin_volume_cm3(extent_cm, dimension)
    return np.full(tuple(int(d) for d in dimension),
                   medium.density_g_cm3 * bin_volume_cm3 * 1e-3)  # g -> kg


def voxel_mass_kg(snapshot: Snapshot, config: TransportConfig) -> np.ndarray:
    """Full-voxel masses (kg), logical (x,y,z), at LATTICE resolution.

    Coarsen with `mesh.coarsen_field` (a block SUM) when the tally mesh is
    coarser than the lattice — mass is extensive, so the coarse bin's mass is
    the sum of the voxel masses it contains, and dividing summed energy by
    summed mass reproduces the mass-weighted mean of the fine dose.

    ponytail: full rectangular-voxel volume; bins intersected by the CSG
    cylinder need openmc RegularMesh.material_volumes() — applied in the
    integration path where openmc is available (see dose.py notes).
    """
    density = np.empty(snapshot.cell_id.shape, dtype=float)
    for name in (MEDIUM, BIOMASS):
        cls = config.materials[name]
        density[voxel_class_array(snapshot, config) == name] = cls.density_g_cm3
    voxel_volume_cm3 = config.voxel_pitch_cm ** 3
    return density * voxel_volume_cm3 * 1e-3  # g -> kg


def _volume_raytrace_model(model, mesh):
    """The same CSG, wrapped in a void box that spans the mesh.

    `material_volumes` raytraces the geometry and refuses with "Mesh not fully
    contained in geometry" when a ray leaves the defined region — which the
    cube mesh always does, because it circumscribes a cylinder and the corners
    are outside it.

    The wrapper is a MEASURING INSTRUMENT, never a transport model: it is built
    here, used for raytracing, and discarded. The real model's boundary
    conditions are untouched, so this cannot change any physics — and the void
    fill contributes no mass, so it cannot change any volume either. Validated
    against the A0 cylinder: pi*r^2*L to within 0.06% at 1e6 rays.
    """
    import openmc

    # ROOT-universe cells only. `get_all_cells()` also returns the cells that
    # fill lattice universes, and those carry no region — they are bounded by
    # the lattice, not by surfaces — so unioning them raises. The nested
    # universes still travel with the lattice fill.
    cells = list(model.geometry.root_universe.cells.values())
    occupied = cells[0].region
    for cell in cells[1:]:
        occupied = occupied | cell.region
    bounds = [v for pair in zip(mesh.lower_left, mesh.upper_right) for v in pair]
    box = openmc.model.RectangularParallelepiped(*bounds, boundary_type="vacuum")
    void = openmc.Cell(name="raytrace_void", fill=None, region=-box & ~occupied)
    return openmc.Model(geometry=openmc.Geometry(cells + [void]),
                        materials=model.materials, settings=model.settings)


def mesh_material_masses_kg(mesh, model, dimension, *, n_samples: int = 10_000_000,
                            max_materials: int = 4) -> np.ndarray:
    """EXACT per-bin masses (kg), logical (x,y,z), from the CSG itself.

        m_v = sum_k rho_k * V_vk

    `voxel_mass_kg` and `phantom_mass_kg` both give every bin the full
    rectangular volume of one material. That is wrong wherever a bin is cut by
    a curved surface: the cylinder wall, and the membrane annulus, which
    `voxel_mass_kg` cannot see at all because the membrane is CSG and never a
    lattice class. It is the mechanism behind A0's ~11% coarse-mesh aggregate
    drift, and both functions name this as their upgrade path.

    TWO PROPERTIES OF THE ANSWER, neither of which may be hidden:

    - It is a RAYTRACING ESTIMATE, so the denominator carries its own O(1/sqrt
      n) sampling error. Record `n_samples` next to any result that uses it,
      and size a mass-conservation tolerance from it rather than asserting
      exact agreement.
    - It needs `openmc.lib.init`, hence cross sections. There is no offline
      path to an exact CSG volume here.

    `max_materials` must cover the most materials in any one bin: biomass,
    medium, membrane and void is already four.
    """
    volumes = mesh.material_volumes(_volume_raytrace_model(model, mesh),
                                    n_samples=n_samples,
                                    max_materials=max_materials)
    # From the GEOMETRY, not `model.materials`: the builders pass materials
    # implicitly through the cells, so `model.materials` is empty and a lookup
    # against it raises KeyError on the first real bin.
    density_by_id = {m.id: m.get_mass_density()
                     for m in model.geometry.get_all_materials().values()}

    dim = tuple(int(d) for d in dimension)
    grams = np.zeros(int(np.prod(dim)), dtype=float)
    for material_id, per_element in volumes.items():
        if material_id is None:            # void contributes no mass
            continue
        # openmc returns numpy integer ids; `Material.id` is a Python int.
        grams += density_by_id[int(material_id)] * np.asarray(per_element, dtype=float)
    # Mesh elements run x fastest, z slowest — the same convention
    # `dose.extract_heating` undoes for the tally.
    return (grams * 1e-3).reshape(dim[::-1]).transpose(2, 1, 0)
