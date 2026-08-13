import numpy as np
import h5py

from biofilm_openmc.dose import DoseResult
from biofilm_openmc.results import (load_transport_result, write_dose_attribution,
                                    write_transport_result)


def _result(n=4):
    rng = np.random.default_rng(2)
    mean = rng.random((n, n, n))
    sd = mean * 0.1
    return DoseResult(mean, sd, sd / mean, 1.0e9)


def test_transport_result_roundtrip(tmp_path):
    res = _result()
    path = tmp_path / "transport_result.h5"
    write_transport_result(path, res, transport_hash="t" * 64,
                           openmc_version="0.15.3",
                           nuclear_data_id="endfb-viii.0",
                           batches=3, particles=1000, seed=7)
    back, attrs = load_transport_result(path)
    assert np.array_equal(back.dose_rate_mean_Gy_s, res.dose_rate_mean_Gy_s)
    assert np.array_equal(back.rel_err, res.rel_err)
    assert attrs["transport_state_hash"] == "t" * 64
    assert attrs["openmc_version"] == "0.15.3"
    assert attrs["seed"] == 7


def test_attribution_file_keys_both_hashes(tmp_path):
    path = tmp_path / "dose_attribution.h5"
    aggregates = {"lineage": {
        1: {"dose_rate_Gy_s": 2.5, "mass_kg": 4.0, "n_voxels": 2},
        3: {"dose_rate_Gy_s": 1.0, "mass_kg": 1.0, "n_voxels": 1},
    }}
    write_dose_attribution(path, aggregates, transport_hash="t" * 64,
                           label_hash="l" * 64,
                           uncertainty_method="replicates")
    with h5py.File(path, "r") as f:
        assert f.attrs["transport_state_hash_used"] == "t" * 64
        assert f.attrs["label_state_hash"] == "l" * 64
        assert f.attrs["uncertainty_method"] == "replicates"
        assert f["lineage/label"][()].tolist() == [1, 3]
        assert np.allclose(f["lineage/dose_rate_Gy_s"][()], [2.5, 1.0])
