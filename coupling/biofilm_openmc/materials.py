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
