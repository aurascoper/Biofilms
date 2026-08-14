import numpy as np
import pytest

from biofilm_openmc.config import (WATER_PHANTOM, ConfigError, load_config,
                                   load_transport_config)
from biofilm_openmc.materials import (BIOMASS, MEDIUM, unique_material_specs,
                                      voxel_class_array, voxel_mass_kg)
from biofilm_openmc.model import check_lattice_congruence

from conftest import VALID_CONFIG, WATER_PHANTOM_CONFIG


def test_cylinder_must_fit_inside_the_voxel_lattice():
    """The lattice is a cube of side n*pitch; a CSG cylinder larger than it
    would be silently padded with medium and fall outside the tally mesh."""
    n = 8                                     # VALID_CONFIG: pitch 0.001 -> 0.008 cm
    check_lattice_congruence(n, load_transport_config(VALID_CONFIG))

    too_wide = load_transport_config(
        VALID_CONFIG.replace("cylinder_radius_cm = 0.004",
                             "cylinder_radius_cm = 0.010"))
    with pytest.raises(ConfigError, match="cylinder_radius_cm"):
        check_lattice_congruence(n, too_wide)

    # a real apparatus aspect ratio (3 mm bore, 1100 mm long) cannot fit
    too_long = load_transport_config(
        VALID_CONFIG.replace("cylinder_length_cm = 0.008",
                             "cylinder_length_cm = 110.0"))
    with pytest.raises(ConfigError, match="cylinder_length_cm"):
        check_lattice_congruence(n, too_long)


def test_voxel_classes_from_occupancy_only(snapshot, config):
    classes = voxel_class_array(snapshot, config)
    assert set(np.unique(classes)) == {MEDIUM, BIOMASS}
    assert (classes[snapshot.cell_id > 0] == BIOMASS).all()
    assert (classes[snapshot.cell_id == 0] == MEDIUM).all()


def test_no_material_per_lineage(snapshot, config):
    # blow the lineage label array up to many distinct values: material
    # count must not move (labels never become materials)
    snapshot.lineage_id[snapshot.cell_id > 0] = np.arange(
        1, 1 + (snapshot.cell_id > 0).sum())
    specs = unique_material_specs(snapshot, config)
    assert len(specs) <= len(config.materials)


def test_composition_dedup(snapshot):
    # membrane duplicated as an exact copy of medium -> collapses to one spec
    dup = VALID_CONFIG.replace(
        """[materials.membrane]
density_g_cm3 = 0.94
  [materials.membrane.elements]
  C = 0.856
  H = 0.144""",
        """[materials.membrane]
density_g_cm3 = 1.0
  [materials.membrane.elements]
  H = 0.111894
  O = 0.888106""")
    cfg = load_config(dup)
    specs = unique_material_specs(snapshot, cfg)
    assert len(specs) == 2  # medium==membrane composition + biomass


_NO_MEMBRANE = VALID_CONFIG.replace(
    """[materials.membrane]
density_g_cm3 = 0.94
  [materials.membrane.elements]
  C = 0.856
  H = 0.144""", "")


def test_missing_required_class_is_loud():
    """The biofilm model needs all three classes, and the loader now says so
    before a model is ever built."""
    with pytest.raises(ConfigError, match="membrane"):
        load_config(_NO_MEMBRANE)


def test_required_classes_still_guards_a_directly_built_config(snapshot, config):
    """Defence in depth: the loader is not the only way to make a config, so
    the mapper keeps its own check."""
    from dataclasses import replace
    stripped = replace(config, materials={k: v for k, v in config.materials.items()
                                          if k != "membrane"})
    with pytest.raises(ConfigError, match="membrane"):
        voxel_class_array(snapshot, stripped)


def test_water_phantom_needs_only_the_medium():
    """The water-phantom gap, closed: a config with no biomass and no membrane
    loads for the phantom and is refused for the biofilm."""
    cfg = load_transport_config(WATER_PHANTOM_CONFIG, kind=WATER_PHANTOM)
    assert set(cfg.materials) == {"medium"}
    with pytest.raises(ConfigError) as exc:
        load_transport_config(WATER_PHANTOM_CONFIG, kind="biofilm_cylinder")
    msg = str(exc.value)
    assert "baseline_biomass" in msg and "membrane" in msg


def test_voxel_mass_positive_and_scaled(snapshot, config):
    m = voxel_mass_kg(snapshot, config)
    assert (m > 0).all()
    # medium voxel: 1.0 g/cm3 * (1e-3 cm)^3 = 1e-9 g = 1e-12 kg
    assert np.isclose(m[snapshot.cell_id == 0][0], 1e-12)
    # biomass voxel is 1.1x denser
    assert np.isclose(m[snapshot.cell_id > 0][0], 1.1e-12)
