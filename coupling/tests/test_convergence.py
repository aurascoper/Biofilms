"""Convergence-metric unit tests (no OpenMC).

The two accounting rules that change the answer are pinned here: the region of
interest excludes structurally-void bins, and cross-resolution comparison
aggregates extensive energy rather than averaging intensive dose.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_openmc.config import WATER_PHANTOM, load_transport_config
from biofilm_openmc.convergence import (field_difference, in_cylinder_mask,
                                        point_metrics, seed_spread)
from biofilm_openmc.dose import EV_TO_J
from biofilm_openmc.mesh import resolve_mesh_dimension

from conftest import WATER_PHANTOM_CONFIG


@pytest.fixture
def cfg():
    return load_transport_config(WATER_PHANTOM_CONFIG, kind=WATER_PHANTOM)


def test_roi_is_the_inscribed_cylinder_not_the_cube(cfg):
    """The mesh spans the circumscribing cube, so ~1 - pi/4 of bins are void by
    construction. Counting them would add a fixed offset to every sparsity
    number and make each mesh look worse than it is."""
    for factor in (1, 2, 4):
        dim = resolve_mesh_dimension((48, 48, 48), factor)
        mask = in_cylinder_mask(cfg, dim)
        frac = mask.mean()
        assert np.isclose(frac, np.pi / 4, atol=0.02), (factor, frac)


def test_zero_score_fraction_ignores_void_bins(cfg):
    dim = resolve_mesh_dimension((48, 48, 48), 4)
    mask = in_cylinder_mask(cfg, dim)
    heating = np.where(mask, 1.0, 0.0)          # everything in the ROI scores
    mass = np.full(dim, 1e-3)
    m = point_metrics(heating, heating * 0.01, mass, mask,
                      e_src_eV=661655.0)
    assert m["zero_score_fraction_roi"] == 0.0
    # ... even though most-of-a-quarter of the full cube scored nothing
    assert (heating == 0).mean() > 0.2


def test_absorbed_fraction_and_total_energy(cfg):
    dim = resolve_mesh_dimension((48, 48, 48), 4)
    mask = in_cylinder_mask(cfg, dim)
    heating = np.zeros(dim)
    heating[0, 0, 0] = 100.0
    mass = np.full(dim, 2e-3)
    m = point_metrics(heating, heating * 0.1, mass, mask, e_src_eV=1000.0)
    assert m["total_heating_eV_per_src"] == 100.0
    assert np.isclose(m["absorbed_fraction"], 0.1)


def test_field_difference_is_zero_for_a_consistent_refinement(cfg):
    """A fine field that aggregates exactly onto the coarse one is not a
    disagreement — the metric must see zero."""
    factor = 2
    fine_dim = resolve_mesh_dimension((48, 48, 48), 1)
    coarse_dim = resolve_mesh_dimension((48, 48, 48), factor)
    rng = np.random.default_rng(0)
    fine = rng.random(fine_dim) * 10.0

    from biofilm_openmc.mesh import coarsen_field
    coarse = coarsen_field(fine, factor)          # exactly consistent
    mass_coarse = np.full(coarse_dim, 1e-3)
    mask = in_cylinder_mask(cfg, coarse_dim)

    assert np.isclose(field_difference(coarse, fine, factor, mass_coarse, mask), 0.0)


def test_field_difference_scales_with_a_uniform_bias(cfg):
    factor = 2
    fine_dim = resolve_mesh_dimension((48, 48, 48), 1)
    coarse_dim = resolve_mesh_dimension((48, 48, 48), factor)
    fine = np.full(fine_dim, 1.0)
    from biofilm_openmc.mesh import coarsen_field
    coarse = coarsen_field(fine, factor) * 1.10   # 10% hot
    mass_coarse = np.full(coarse_dim, 1e-3)
    mask = in_cylinder_mask(cfg, coarse_dim)
    assert np.isclose(
        field_difference(coarse, fine, factor, mass_coarse, mask), 0.10)


def test_field_difference_uses_energy_not_a_dose_average(cfg):
    """With non-uniform mass, comparing dose averages and comparing deposited
    energy give different answers; the metric must use energy."""
    factor = 2
    fine_dim = resolve_mesh_dimension((48, 48, 48), 1)
    coarse_dim = resolve_mesh_dimension((48, 48, 48), factor)
    from biofilm_openmc.mesh import coarsen_field

    rng = np.random.default_rng(3)
    fine = rng.random(fine_dim) * 10.0
    coarse = coarsen_field(fine, factor)
    mass_coarse = rng.random(coarse_dim) * 1e-3 + 1e-4
    mask = in_cylinder_mask(cfg, coarse_dim)
    # energy-consistent refinement -> exactly zero regardless of the mass field
    assert np.isclose(
        field_difference(coarse, fine, factor, mass_coarse, mask), 0.0)


def test_seed_spread_needs_replicates():
    assert np.isnan(seed_spread([1.0]))
    assert np.isclose(seed_spread([1.0, 1.0, 1.0]), 0.0)
    assert seed_spread([1.0, 1.1, 0.9]) > 0


def test_mass_weighted_mean_dose_matches_the_hand_formula(cfg):
    dim = resolve_mesh_dimension((48, 48, 48), 8)
    mask = in_cylinder_mask(cfg, dim)
    rng = np.random.default_rng(5)
    heating = rng.random(dim) * 100
    mass = rng.random(dim) * 1e-3 + 1e-4
    m = point_metrics(heating, heating * 0.05, mass, mask, e_src_eV=661655.0)
    dose = heating * EV_TO_J / mass
    expected = (dose[mask] * mass[mask]).sum() / mass[mask].sum()
    assert np.isclose(m["mass_weighted_mean_dose_Gy_per_src"], expected)


def test_resolution_loss_sees_what_field_difference_cannot(cfg):
    """Nested meshes make field_difference blind to discretization: summing a
    fine field onto a coarse grid IS that grid's exact answer. resolution_loss
    measures the smearing instead, on the fine grid where the structure is."""
    from biofilm_openmc.convergence import resolution_loss
    from biofilm_openmc.mesh import coarsen_field

    factor = 4
    fine_dim = resolve_mesh_dimension((48, 48, 48), 1)
    coarse_dim = resolve_mesh_dimension((48, 48, 48), factor)

    # a strong radial gradient: exactly the structure a coarse mesh destroys
    xs = np.arange(fine_dim[0]) - fine_dim[0] / 2
    r = np.sqrt(xs[:, None, None] ** 2 + xs[None, :, None] ** 2
                + np.zeros((1, 1, fine_dim[2])))
    fine = np.exp(-r / 4.0) + 1e-3

    coarse = coarsen_field(fine, factor)
    fine_mass = np.full(fine_dim, 1e-3)
    coarse_mass = coarsen_field(fine_mass, factor)
    fine_mask = in_cylinder_mask(cfg, fine_dim)
    coarse_mask = in_cylinder_mask(cfg, coarse_dim)

    # energy-consistent, so the aggregate comparison sees nothing at all
    assert np.isclose(
        field_difference(coarse, fine, factor, coarse_mass, coarse_mask), 0.0)
    # ... while the gradient really has been smeared
    loss = resolution_loss(coarse, coarse_mass, factor, fine, fine_mass, fine_mask)
    assert loss > 0.05, loss


def test_resolution_loss_is_zero_at_factor_one(cfg):
    """No coarsening, no loss — the metric must not manufacture one."""
    from biofilm_openmc.convergence import resolution_loss

    dim = resolve_mesh_dimension((48, 48, 48), 1)
    rng = np.random.default_rng(7)
    field = rng.random(dim) + 0.1
    mass = np.full(dim, 1e-3)
    mask = in_cylinder_mask(cfg, dim)
    assert np.isclose(
        resolution_loss(field, mass, 1, field, mass, mask), 0.0)
