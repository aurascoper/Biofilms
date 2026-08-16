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


def test_generation_zero_is_a_founder_not_background():
    # Four founder parcels (generation 0) and four background voxels, which
    # also read 0. The `labels > 0` fallback drops the founders entirely; the
    # occupancy mask is what separates the two meanings of zero.
    cell_id = np.array([[[0, 0], [0, 0]], [[1, 2], [3, 4]]])
    generation = np.array([[[0, 0], [0, 0]], [[0, 0], [1, 1]]])
    dose = np.full((2, 2, 2), 5.0)
    mass = np.ones((2, 2, 2))

    assert set(aggregate_by_label(dose, generation, mass)) == {1}

    occupied = cell_id > 0
    out = aggregate_by_label(dose, generation, mass, occupied=occupied)
    assert set(out) == {0, 1}
    assert out[0]["n_voxels"] == 2      # the two generation-0 parcels, not the background
    assert out[0]["dose_rate_Gy_s"] == 5.0


def test_wall_voxels_are_refused_rather_than_miscounted():
    dose = np.ones((1, 1, 3))
    labels = np.array([[[-1, 0, 1]]])          # -1 is the wall sentinel
    with pytest.raises(ValueError, match="negative label"):
        aggregate_by_label(dose, labels, np.ones_like(dose),
                           occupied=np.ones((1, 1, 3), dtype=bool))


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
