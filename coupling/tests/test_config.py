from pathlib import Path

import pytest

from biofilm_openmc.config import (FIELD_SPECS, STAGES, ConfigError,
                                   DoseRateConfig, TransportConfig, load_config,
                                   load_cpm_feedback_config,
                                   load_dose_rate_config,
                                   load_transport_config, required_keys)

from conftest import VALID_CONFIG

TEMPLATE = Path(__file__).parents[2] / "config" / "coupling_template.toml"

# Keys the CPM/membrane stages add — transport must never demand these.
_BIOLOGICAL_KEYS = ["seconds_per_mcs", "hamiltonian_scale", "melanin_scale",
                    "membrane_statistic"]


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


# --- staged contract -------------------------------------------------------

def test_required_keys_are_cumulative_across_stages():
    prev = frozenset()
    for stage in STAGES:
        keys = required_keys(stage)
        assert prev <= keys, f"{stage} dropped a key required earlier"
        prev = keys
    assert required_keys(STAGES[-1]) == {f.dotted for f in FIELD_SPECS}


def test_transport_stage_does_not_require_time_or_biological_scales():
    keys = required_keys("transport")
    for key in _BIOLOGICAL_KEYS + ["photons_per_second"]:
        assert not any(k.endswith(key) for k in keys), \
            f"transport must not require {key}"


def test_transport_config_loads_without_time_or_biological_scales():
    """The whole point: a config with only transport inputs must build a model,
    so the resolution/history sweep can run before any assay certificate or
    biological calibration exists."""
    toml = VALID_CONFIG
    for line in ["seconds_per_mcs = 1.0", "hamiltonian_scale = 1.0",
                 "melanin_scale = 1.0", 'membrane_statistic = "mass_weighted"',
                 "photons_per_second = 1.0e9"]:
        toml = toml.replace(line, "")

    cfg = load_transport_config(toml)
    assert isinstance(cfg, TransportConfig)
    assert cfg.voxel_pitch_cm == 0.001
    assert not hasattr(cfg, "photons_per_second")

    with pytest.raises(ConfigError) as exc:
        load_dose_rate_config(toml)
    assert "photons_per_second" in str(exc.value)

    with pytest.raises(ConfigError) as exc:
        load_cpm_feedback_config(toml)
    msg = str(exc.value)
    for key in ["seconds_per_mcs", "hamiltonian_scale", "melanin_scale"]:
        assert key in msg
    assert "membrane_statistic" not in msg   # that one is membrane_feedback


def test_stronger_stage_is_a_weaker_stage(config):
    """Inheritance direction: anything accepting a TransportConfig accepts a
    full config, so build_model needs no unwrapping."""
    assert isinstance(config, TransportConfig)
    assert isinstance(config, DoseRateConfig)
    assert isinstance(load_transport_config(VALID_CONFIG), TransportConfig)


def test_template_still_fails_loudly_at_every_stage():
    for stage in STAGES:
        with pytest.raises(ConfigError, match="voxel_pitch_cm"):
            load_config(str(TEMPLATE)) if stage == STAGES[-1] else \
                load_transport_config(str(TEMPLATE))


# --- spectrum contract -----------------------------------------------------

def test_spectrum_probabilities_must_sum_to_one():
    """The loader does not normalize; raw per-decay yields must be converted
    to a PMF by the caller, keeping S = A*sum(Y) separate from p_j."""
    broken = VALID_CONFIG.replace("spectrum_probabilities = [1.0]",
                                  "spectrum_probabilities = [0.85]")
    with pytest.raises(ConfigError, match="sum to 0.85, not 1"):
        load_config(broken)


def test_spectrum_rejects_negative_and_zero_energy():
    neg = VALID_CONFIG.replace("spectrum_probabilities = [1.0]",
                               "spectrum_probabilities = [-1.0]")
    with pytest.raises(ConfigError, match="non-negative"):
        load_config(neg)
    zero = VALID_CONFIG.replace("spectrum_energies_eV = [1.0e6]",
                                "spectrum_energies_eV = [0.0]")
    with pytest.raises(ConfigError, match="positive"):
        load_config(zero)
