from pathlib import Path

import pytest

from biofilm_openmc.config import ConfigError, load_config

from conftest import VALID_CONFIG

TEMPLATE = Path(__file__).parents[2] / "config" / "coupling_template.toml"


def test_shipped_template_fails_loudly_naming_every_required_key():
    with pytest.raises(ConfigError) as exc:
        load_config(str(TEMPLATE))
    msg = str(exc.value)
    for key in ["voxel_pitch_cm", "origin_cm", "cylinder_radius_cm",
                "cylinder_length_cm", "membrane_thickness_cm", "axial",
                "radial_outer", "seconds_per_mcs", "photons_per_second",
                "spectrum_energies_eV", "spectrum_probabilities", "spatial",
                "angular", "hamiltonian_scale", "melanin_scale",
                "membrane_statistic", "[materials]"]:
        assert key in msg, f"missing key {key!r} not named in the error"


def test_valid_config_loads(config):
    assert config.voxel_pitch_cm == 0.001
    assert config.membrane_statistic == "mass_weighted"
    assert set(config.materials) == {"medium", "baseline_biomass", "membrane"}
    assert config.batches == 5


def test_single_missing_key_is_reported():
    broken = VALID_CONFIG.replace("seconds_per_mcs = 1.0", "")
    with pytest.raises(ConfigError, match="seconds_per_mcs"):
        load_config(broken)


def test_spectrum_length_mismatch():
    broken = VALID_CONFIG.replace("spectrum_probabilities = [1.0]",
                                  "spectrum_probabilities = [0.5, 0.5]")
    with pytest.raises(ConfigError, match="length mismatch"):
        load_config(broken)


def test_element_fractions_must_sum_to_one():
    broken = VALID_CONFIG.replace("H = 0.111894", "H = 0.5")
    with pytest.raises(ConfigError, match="sum to"):
        load_config(broken)


def test_invalid_enum_value():
    broken = VALID_CONFIG.replace('axial = "vacuum"', 'axial = "mirror"')
    with pytest.raises(ConfigError, match="axial"):
        load_config(broken)
