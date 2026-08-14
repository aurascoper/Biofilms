"""Tally-mesh resolution: nesting, extensive aggregation, and the f-cubed trap.

The dangerous bug this file exists to catch has no shape mismatch to give it
away: if the tally mesh is coarsened by f while the dose denominator still
comes from `voxel_pitch_cm`, every dose is wrong by exactly f**3 and every
array still lines up.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_openmc.config import ConfigError, load_transport_config
from biofilm_openmc.dose import EV_TO_J, specific_energy_per_source
from biofilm_openmc.mesh import (biofilm_mesh_extent_cm, coarsen_field,
                                 mesh_bin_volume_cm3, phantom_mesh_extent_cm,
                                 resolve_mesh_dimension, upsample_field)

from conftest import VALID_CONFIG, WATER_PHANTOM_CONFIG


def test_resolve_requires_nesting():
    assert resolve_mesh_dimension((60, 60, 60), 4) == (15, 15, 15)
    assert resolve_mesh_dimension((12, 12, 12), 1) == (12, 12, 12)
    with pytest.raises(ConfigError, match="does not divide"):
        resolve_mesh_dimension((60, 60, 60), 7)
    with pytest.raises(ConfigError, match=">= 1"):
        resolve_mesh_dimension((60, 60, 60), 0)


def test_coarsen_is_a_sum_and_conserves_the_total():
    rng = np.random.default_rng(0)
    fine = rng.random((8, 8, 8))
    coarse = coarsen_field(fine, 2)
    assert coarse.shape == (4, 4, 4)
    assert np.isclose(coarse.sum(), fine.sum())          # extensive: conserved
    assert np.isclose(coarse[0, 0, 0], fine[0:2, 0:2, 0:2].sum())
    with pytest.raises(ConfigError, match="does not divide"):
        coarsen_field(fine, 3)


def test_coarse_dose_equals_the_mass_weighted_mean_of_the_fine_dose():
    """THE invariant. Coarsening energy and mass separately, then dividing,
    must reproduce the mass-weighted mean of the fine dose over each block —
    which is what makes a coarse run a legitimate estimate of the fine one."""
    rng = np.random.default_rng(1)
    f = 2
    heating_eV = rng.random((8, 8, 8)) * 1e3
    mass_kg = rng.random((8, 8, 8)) * 1e-9 + 1e-10     # non-uniform on purpose

    fine_dose = specific_energy_per_source(heating_eV, mass_kg)
    coarse_dose = specific_energy_per_source(coarsen_field(heating_eV, f),
                                             coarsen_field(mass_kg, f))

    # block-wise mass-weighted mean of the fine dose
    expected = (coarsen_field(fine_dose * mass_kg, f) / coarsen_field(mass_kg, f))
    assert np.allclose(coarse_dose, expected)

    # and total deposited energy is identical at both resolutions
    assert np.isclose((fine_dose * mass_kg).sum(),
                      (coarse_dose * coarsen_field(mass_kg, f)).sum())


def test_using_the_lattice_pitch_after_coarsening_is_wrong_by_f_cubed():
    """The regression guard. A coarsened mesh whose mass still comes from the
    LATTICE pitch understates dose by exactly the coarsening factor cubed."""
    f = 2
    pitch_cm = 0.001
    n = 8
    heating_eV = np.full((n, n, n), 10.0)
    density_g_cm3 = 1.0

    voxel_mass = density_g_cm3 * pitch_cm ** 3 * 1e-3
    fine_mass = np.full((n, n, n), voxel_mass)

    right = specific_energy_per_source(coarsen_field(heating_eV, f),
                                       coarsen_field(fine_mass, f))
    wrong = specific_energy_per_source(coarsen_field(heating_eV, f),
                                       np.full((n // f,) * 3, voxel_mass))
    assert np.allclose(wrong / right, f ** 3)


def test_bin_volume_comes_from_extent_not_pitch():
    cfg = load_transport_config(WATER_PHANTOM_CONFIG, kind="water_phantom")
    extent = phantom_mesh_extent_cm(cfg)
    assert extent == (30.0, 30.0, 30.0)                 # 2r x 2r x L
    fine = mesh_bin_volume_cm3(extent, resolve_mesh_dimension((12, 12, 12), 1))
    coarse = mesh_bin_volume_cm3(extent, resolve_mesh_dimension((12, 12, 12), 2))
    assert np.isclose(fine, 2.5 ** 3)
    assert np.isclose(coarse / fine, 8.0)               # exactly f**3


def test_biofilm_extent_is_the_lattice_cube_and_does_not_move_with_coarsening():
    """Physical size fixed, resolution free — the whole point of the split."""
    cfg = load_transport_config(VALID_CONFIG)
    assert biofilm_mesh_extent_cm(cfg, 8) == (0.008, 0.008, 0.008)
    coarsened = load_transport_config(
        VALID_CONFIG.replace("[transport]",
                             "[transport]\n  [transport.mesh]\n  coarsening_factor = 2\n"))
    assert biofilm_mesh_extent_cm(coarsened, 8) == biofilm_mesh_extent_cm(cfg, 8)
    assert resolve_mesh_dimension((8, 8, 8), coarsened.mesh_coarsening_factor) == (4, 4, 4)


def test_upsample_broadcasts_an_intensive_field():
    """Dose rate is intensive, so a coarse bin's rate applies to every fine
    voxel inside it. This is what keeps lineage attribution meaningful."""
    coarse = np.arange(8, dtype=float).reshape(2, 2, 2)
    fine = upsample_field(coarse, 2)
    assert fine.shape == (4, 4, 4)
    assert np.allclose(fine[0:2, 0:2, 0:2], coarse[0, 0, 0])
    # round trip: coarsening an upsampled intensive field scales by f**3
    assert np.allclose(coarsen_field(fine, 2), coarse * 8)


def test_coarsening_factor_must_be_a_positive_integer():
    for bad in ("2", 0, -1, 1.5):
        broken = VALID_CONFIG.replace(
            "[transport]",
            f"[transport]\n  [transport.mesh]\n  coarsening_factor = {bad!r}\n")
        with pytest.raises(ConfigError, match="coarsening_factor"):
            load_transport_config(broken)
