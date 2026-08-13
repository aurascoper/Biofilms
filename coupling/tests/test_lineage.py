import numpy as np
import pytest

from biofilm_openmc.lineage import (aggregate_by_label, quadrature_uncertainty_approx,
                                    replicate_uncertainty)


def test_mass_weighted_aggregation_exact():
    dose = np.array([[[1.0, 3.0]]])          # two voxels
    labels = np.array([[[1, 1]]])
    mass = np.array([[[1.0, 3.0]]])          # unequal masses
    out = aggregate_by_label(dose, labels, mass)
    # mass-weighted: (1*1 + 3*3) / 4 = 2.5 — NOT the unweighted 2.0
    assert np.isclose(out[1]["dose_rate_Gy_s"], 2.5)
    assert out[1]["mass_kg"] == 4.0
    assert out[1]["n_voxels"] == 2


def test_label_zero_is_skipped_and_labels_separate():
    dose = np.array([[[1.0, 2.0, 8.0]]])
    labels = np.array([[[0, 2, 5]]])
    mass = np.ones_like(dose)
    out = aggregate_by_label(dose, labels, mass)
    assert set(out) == {2, 5}
    assert out[5]["dose_rate_Gy_s"] == 8.0


def test_replicate_uncertainty():
    r1 = {1: {"dose_rate_Gy_s": 1.0}}
    r2 = {1: {"dose_rate_Gy_s": 3.0}}
    out = replicate_uncertainty([r1, r2])
    assert out[1]["dose_rate_Gy_s"] == 2.0
    assert np.isclose(out[1]["sd_Gy_s"], np.sqrt(2.0))
    with pytest.raises(ValueError, match="replicates"):
        replicate_uncertainty([r1])


def test_quadrature_is_available_but_separate():
    sd = np.array([[[0.3, 0.4]]])
    labels = np.array([[[7, 7]]])
    mass = np.ones_like(sd)
    out = quadrature_uncertainty_approx(sd, labels, mass)
    assert np.isclose(out[7], 0.25)  # sqrt(.09+.16)/2
