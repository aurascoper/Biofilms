import numpy as np
import pytest

from biofilm_openmc.dose import (EV_TO_J, dose_from_statepoint, normalize_heating,
                                 sparsity_report)

from conftest import heating_statepoint_from_logical


def test_normalize_heating_hand_computed():
    # 1 eV/src * 1.602176634e-19 J/eV * 1e9 src/s / 1e-6 kg = 1.602176634e-4 Gy/s
    h = np.array([1.0])
    out = normalize_heating(h, 1.0e9, np.array([1.0e-6]))
    assert np.isclose(out[0], 1.602176634e-4, rtol=1e-12)
    assert EV_TO_J == 1.602176634e-19


def test_nonpositive_mass_is_refused():
    with pytest.raises(ValueError, match="mass"):
        normalize_heating(np.ones(2), 1.0, np.array([1.0, 0.0]))


def test_statepoint_bin_order_recovers_logical_field():
    # distinct value per voxel: L[x,y,z] = x + 10y + 100z
    n = 3
    x, y, z = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    L = (x + 10 * y + 100 * z).astype(float)
    sp = heating_statepoint_from_logical(L)
    mass = np.full((n, n, n), 1.0)
    res = dose_from_statepoint(sp, mass, 1.0 / EV_TO_J)  # scale cancels units
    assert np.allclose(res.dose_rate_mean_Gy_s, L)


def test_rel_err_and_sparsity_report():
    n = 4
    mean = np.zeros((n, n, n))
    mean[0, 0, 0] = 2.0
    sd = np.zeros_like(mean)
    sd[0, 0, 0] = 0.5
    sp = heating_statepoint_from_logical(mean, sd)
    res = dose_from_statepoint(sp, np.ones_like(mean), 1.0 / EV_TO_J)
    assert np.isclose(res.rel_err[0, 0, 0], 0.25)
    assert np.isinf(res.rel_err[1, 1, 1])

    occupied = np.zeros_like(mean, dtype=bool)
    occupied[0, 0, 0] = True
    occupied[1, 1, 1] = True
    rep = sparsity_report(res, occupied)
    assert rep["occupied_voxels"] == 2
    assert rep["occupied_unscored_fraction"] == 0.5
    assert np.isclose(rep["median_rel_err"], 0.25)
