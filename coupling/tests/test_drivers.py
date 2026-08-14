"""Driver tests without OpenMC: gate reporting, scenario merge, comparison
math, and the scan loop with an injected fake runner."""

from pathlib import Path

import numpy as np
import pytest

from biofilm_openmc.dose import DoseResult
from biofilm_openmc.drivers import (GATE_FAILED, GATE_NOT_EVALUATED, GATE_PASSED,
                                    compare_fields, deep_merge, evaluate_gate,
                                    main_scan, scan, scenario_config)

from conftest import VALID_CONFIG, heating_statepoint_from_logical

TEMPLATE = Path(__file__).parents[2] / "config" / "coupling_template.toml"


def test_unset_config_reports_not_evaluated(tmp_path, snapshot_path, capsys):
    rc = main_scan(["--snapshot", str(snapshot_path),
                    "--config", str(TEMPLATE)])
    out = capsys.readouterr().out
    assert rc == 2
    assert GATE_NOT_EVALUATED in out
    assert "voxel_pitch_cm" in out          # the ConfigError text is surfaced


def test_scenario_override_merges(tmp_path):
    override = tmp_path / "denser.toml"
    override.write_text("[materials.baseline_biomass]\ndensity_g_cm3 = 1.3\n"
                        "  [materials.baseline_biomass.elements]\n"
                        "  H = 0.10\n  C = 0.50\n  N = 0.10\n  O = 0.30\n")
    cfg = scenario_config(VALID_CONFIG, override)
    assert cfg.materials["baseline_biomass"].density_g_cm3 == 1.3
    assert cfg.materials["medium"].density_g_cm3 == 1.0   # untouched


def test_deep_merge_nested():
    assert deep_merge({"a": {"b": 1, "c": 2}}, {"a": {"b": 9}}) == \
        {"a": {"b": 9, "c": 2}}


def _res(mean, sd_frac=0.01):
    mean = np.asarray(mean, dtype=float)
    sd = mean * sd_frac
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.where(mean > 0, sd / mean, np.inf)
    return DoseResult(mean, sd, rel, 1.0)


def test_compare_fields_floor_and_mask():
    base = _res([[[1.0, 1e-12]]])           # second voxel under the floor
    var = _res([[[1.2, 5.0]]])
    mass = np.ones((1, 1, 2))
    m = compare_fields(base, var, mass, dose_floor_Gy_s=1e-6)
    assert m["masked_voxels"] == 1          # floored voxel excluded
    assert np.isclose(m["effect_rel"], 0.2)


def test_gate_decision_paths():
    passed = {"hot": {"effect_rel": 0.30, "noise_rel": 0.01,
                      "total_deposited_ratio": 1.3, "masked_voxels": 10}}
    below_threshold = {"hot": {"effect_rel": 0.01, "noise_rel": 0.001,
                               "total_deposited_ratio": 1.0, "masked_voxels": 10}}
    below_noise = {"hot": {"effect_rel": 0.10, "noise_rel": 0.50,
                           "total_deposited_ratio": 1.1, "masked_voxels": 10}}
    assert evaluate_gate(passed, 0.05)[0] == GATE_PASSED
    assert evaluate_gate(below_threshold, 0.05)[0] == GATE_FAILED
    assert evaluate_gate(below_noise, 0.05)[0] == GATE_FAILED
    assert evaluate_gate({}, 0.05)[0] == GATE_NOT_EVALUATED


def test_scan_loop_with_fake_runner(tmp_path, snapshot, config, monkeypatch):
    """Full scan loop: fake transport returns a hotter field for the denser
    scenario; identical physical state is run once (exact-hash reuse)."""
    from biofilm_openmc.config import load_config

    denser = load_config(VALID_CONFIG.replace(
        'density_g_cm3 = 1.1', 'density_g_cm3 = 1.2'))
    same_as_base = load_config(VALID_CONFIG)

    n = snapshot.cell_id.shape[0]
    calls = []

    def fake_runner(model, name):
        calls.append(name)
        hot = 2.0 if name == "denser" else 1.0
        field = np.full((n, n, n), hot)
        return heating_statepoint_from_logical(field, field * 0.01)

    # build_model imports openmc — stub it out for the fake path
    import biofilm_openmc.drivers as drv
    monkeypatch.setattr(drv, "build_model", lambda s, c: None, raising=False)
    import biofilm_openmc.model as model_mod
    monkeypatch.setattr(model_mod, "build_model", lambda s, c: None)

    report = scan(snapshot, config,
                  {"denser": denser, "duplicate": same_as_base},
                  fake_runner, tmp_path, effect_threshold=0.05,
                  dose_floor_Gy_s=0.0, nuclear_data_id="test",
                  openmc_version="fake")

    # duplicate scenario matched baseline's transport hash: never re-run
    assert sorted(calls) == ["baseline", "denser"]
    assert report["gate"] == GATE_PASSED
    assert (tmp_path / "transport_result_baseline.h5").exists()
    assert (tmp_path / "transport_result_denser.h5").exists()
    assert not (tmp_path / "transport_result_duplicate.h5").exists()
    assert report["metrics"]["duplicate"]["effect_rel"] == 0.0
