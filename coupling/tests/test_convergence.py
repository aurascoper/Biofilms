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


def test_the_three_column_classes_cover_the_csv_exactly():
    """A column that belongs to no class is a column nobody decided how to
    compare — it would be silently skipped by compare_sweeps."""
    from biofilm_openmc.convergence import (CSV_COLUMNS, DETERMINISTIC_COLUMNS,
                                            ENVIRONMENTAL_COLUMNS,
                                            REFERENCE_RELATIVE_COLUMNS,
                                            STOCHASTIC_COLUMNS)
    classes = (DETERMINISTIC_COLUMNS, STOCHASTIC_COLUMNS,
               REFERENCE_RELATIVE_COLUMNS, ENVIRONMENTAL_COLUMNS)
    assert set().union(*classes) == set(CSV_COLUMNS)
    assert sum(len(c) for c in classes) == len(CSV_COLUMNS)


def test_compare_sweeps_ignores_the_machine_and_not_the_physics():
    from biofilm_openmc.convergence import compare_sweeps

    base = {"coarsening_factor": "2", "histories": "4000000", "seed": "1",
            "mesh_nx": "24", "mesh_ny": "24", "mesh_nz": "24",
            "bin_volume_cm3": "7.4", "batches": "80", "particles": "50000",
            "roi_bins": "6000", "heating_floor_eV": "1e-9",
            "total_heating_eV_per_src": "100.0", "absorbed_fraction": "0.5",
            "zero_score_fraction_roi": "0.0", "bins_above_floor": "6000",
            "median_rel_err": "0.067", "p90_rel_err": "0.101",
            "max_rel_err": "0.2", "mass_weighted_mean_dose_Gy_per_src": "1e-12",
            "total_mass_kg": "21.2", "field_diff_vs_reference": "0.0",
            "resolution_loss_vs_reference": "0.177",
            "runtime_s": "2.9", "histories_per_s": "1379310",
            "statepoint_bytes": "500000", "peak_child_rss_kb_cumulative": "80000"}

    # A three-times-slower host on a bigger machine: not a difference.
    slower = {**base, "runtime_s": "9.1", "histories_per_s": "439560",
              "peak_child_rss_kb_cumulative": "240000",
              "statepoint_bytes": "500123"}
    assert compare_sweeps([base], [slower]) == []

    # Monte Carlo jitter inside the tolerance: not a difference either.
    jitter = {**base, "median_rel_err": "0.0672"}
    assert compare_sweeps([base], [jitter]) == []

    # A changed mesh IS a difference, and an exact one.
    remeshed = {**base, "mesh_nx": "12"}
    assert any("mesh_nx changed" in p for p in compare_sweeps([base], [remeshed]))

    # So is a dose that moved more than the declared tolerance.
    moved = {**base, "mass_weighted_mean_dose_Gy_per_src": "1.5e-12"}
    problems = compare_sweeps([base], [moved])
    assert any("mass_weighted_mean_dose_Gy_per_src moved" in p for p in problems)

    # And a missing row is never silently tolerated.
    assert any("missing" in p for p in compare_sweeps([base], []))


def test_reference_relative_columns_are_only_compared_within_one_sweep():
    """field_diff_vs_reference is computed against the sweep's OWN reference
    row, so a sweep over factors 1,2 and one over 1,2,4,8,16 give the same row
    different values. Measured: re-running A0's factor-1 row alone reproduces
    every other column exactly and moves this one by 8%."""
    from biofilm_openmc.convergence import compare_sweeps

    a = {"coarsening_factor": "1", "histories": "4000000", "seed": "1",
         "field_diff_vs_reference": "0.08960055712264413"}
    b = {**a, "field_diff_vs_reference": "0.09700926009429321"}
    assert compare_sweeps([a], [b]) == []
    assert any("field_diff_vs_reference moved" in p
               for p in compare_sweeps([a], [b], same_reference=True))
