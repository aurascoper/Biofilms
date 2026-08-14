import numpy as np
import pytest

from biofilm_openmc.config import ConfigError, load_transport_config
from biofilm_openmc.model import check_lattice_congruence

from conftest import VALID_CONFIG


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
import pytest

from biofilm_openmc.config import ConfigError, load_config
from biofilm_openmc.materials import (BIOMASS, MEDIUM, unique_material_specs,
                                      voxel_class_array, voxel_mass_kg)

from conftest import VALID_CONFIG


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


def test_missing_required_class_is_loud(snapshot, config):
    no_membrane = VALID_CONFIG.replace(
        """[materials.membrane]
density_g_cm3 = 0.94
  [materials.membrane.elements]
  C = 0.856
  H = 0.144""", "")
    cfg_missing = load_config(no_membrane)
    with pytest.raises(ConfigError, match="membrane"):
        voxel_class_array(snapshot, cfg_missing)


def test_voxel_mass_positive_and_scaled(snapshot, config):
    m = voxel_mass_kg(snapshot, config)
    assert (m > 0).all()
    # medium voxel: 1.0 g/cm3 * (1e-3 cm)^3 = 1e-9 g = 1e-12 kg
    assert np.isclose(m[snapshot.cell_id == 0][0], 1e-12)
    # biomass voxel is 1.1x denser
    assert np.isclose(m[snapshot.cell_id > 0][0], 1.1e-12)
