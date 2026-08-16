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
                                 coarsen_ratio, mesh_bin_volume_cm3,
                                 phantom_mesh_extent_cm,
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


# --- tally refinement: the other end of the resolution axis ------------------

def test_refinement_multiplies_the_base_dimension():
    """The finest tally was never limited by OpenMC. A RegularMesh sets its
    extent and its bin count independently, and `biofilm_mesh_extent_cm` already
    ignores the dimension — what blocked refinement was this function returning
    `d // f` with no multiplicative path."""
    assert resolve_mesh_dimension((20, 20, 20), 1, 1) == (20, 20, 20)
    assert resolve_mesh_dimension((20, 20, 20), 1, 2) == (40, 40, 40)
    assert resolve_mesh_dimension((20, 20, 20), 1, 4) == (80, 80, 80)
    # coarsening is untouched by the new parameter
    assert resolve_mesh_dimension((20, 20, 20), 5) == (4, 4, 4)


def test_refinement_and_coarsening_may_not_both_be_set():
    """They are opposite ends of one axis. Composing them silently would give a
    ratio nobody declared, and `mesh_coarsening_factor` carries
    sampling_role = convergence_axis in the uncertainty ledger."""
    with pytest.raises(ConfigError, match="opposite ends"):
        resolve_mesh_dimension((20, 20, 20), 2, 2)
    with pytest.raises(ConfigError, match=">= 1"):
        resolve_mesh_dimension((20, 20, 20), 1, 0)


def _with_mesh(**keys):
    body = "".join(f"  {k} = {v!r}\n" for k, v in keys.items())
    return VALID_CONFIG.replace("[transport]",
                                "[transport]\n  [transport.mesh]\n" + body)


def test_refinement_keeps_the_physical_extent_fixed():
    """Refinement must add bins, never move the apparatus — the whole reason
    the mesh module exists is that the two used to be welded together."""
    cfg = load_transport_config(VALID_CONFIG)
    refined = load_transport_config(_with_mesh(refinement_factor=4))
    assert refined.mesh_refinement_factor == 4
    assert biofilm_mesh_extent_cm(refined, 8) == biofilm_mesh_extent_cm(cfg, 8)
    assert resolve_mesh_dimension((8, 8, 8), 1, refined.mesh_refinement_factor) \
        == (32, 32, 32)


def test_config_refuses_a_bad_or_contradictory_refinement():
    for bad in ("2", 0, -1, 1.5):
        with pytest.raises(ConfigError, match="refinement_factor"):
            load_transport_config(_with_mesh(refinement_factor=bad))
    with pytest.raises(ConfigError, match="opposite ends"):
        load_transport_config(_with_mesh(coarsening_factor=2, refinement_factor=2))


def test_coarsen_ratio_is_the_mass_weighted_mean_not_the_plain_mean():
    """Dose is intensive, so reducing it means summing the extensive numerator
    and denominator separately. A plain block-average would weight a bin that
    is 2% biomass the same as one that is full."""
    rng = np.random.default_rng(3)
    energy = rng.uniform(1.0, 5.0, (4, 4, 4))
    mass = rng.uniform(0.1, 2.0, (4, 4, 4))
    got = coarsen_ratio(energy, mass, 2)

    blocks = lambda a: (a.reshape(2, 2, 2, 2, 2, 2)
                        .transpose(0, 2, 4, 1, 3, 5).reshape(2, 2, 2, 8))
    w = blocks(mass)
    expected = (blocks(energy / mass) * w).sum(-1) / w.sum(-1)
    assert np.allclose(got, expected)
    # and it is NOT the unweighted mean, or this test would prove nothing
    assert not np.allclose(got, blocks(energy / mass).mean(-1))


def test_coarsen_ratio_returns_zero_where_there_is_no_mass():
    """Empty corner bins are normal once the denominator is the exact CSG
    volume rather than the full bin, so this must not be a division error."""
    got = coarsen_ratio(np.ones((2, 2, 2)), np.zeros((2, 2, 2)), 2)
    assert got.shape == (1, 1, 1) and got[0, 0, 0] == 0.0


def test_refining_then_summing_back_recovers_the_coarse_field_exactly():
    """The conservation identity the refinement study rests on: a tally mesh is
    a histogram of the same events, so block-summing a finer histogram returns
    the coarse bin exactly. If this ever fails, either the remap is not
    conservative or refinement perturbed the transport."""
    rng = np.random.default_rng(4)
    fine_energy = rng.uniform(1.0, 5.0, (8, 8, 8))
    fine_mass = rng.uniform(0.1, 2.0, (8, 8, 8))
    coarse_dose = coarsen_ratio(fine_energy, fine_mass, 2)
    # the same reduction done in two steps must agree with one step
    two_step = coarsen_ratio(coarsen_field(fine_energy, 2),
                             coarsen_field(fine_mass, 2), 1)
    assert np.allclose(coarse_dose, two_step)


def test_the_scan_path_refuses_refinement_before_spending_the_histories():
    """`voxel_mass_kg` is defined at the lattice and its only reducer is a block
    sum, so a refined biofilm tally has no denominator in the scan path. Left
    alone the mismatch surfaces inside `per_source_from_statepoint` as a reshape
    failure — AFTER the transport run has been paid for."""
    from biofilm_openmc.drivers import _biofilm_scan_mass

    cfg = load_transport_config(_with_mesh(refinement_factor=2))
    with pytest.raises(ConfigError, match="not supported by the scan path"):
        _biofilm_scan_mass(object(), cfg)


def test_the_water_phantom_honours_refinement_because_its_hash_claims_it():
    """`water_phantom_state_hash` records the refined dimension. A builder that
    ignored refinement would give a run an identity describing a mesh it never
    used — a cache key for a different model."""
    from biofilm_openmc.fingerprint import water_phantom_state_hash

    base = load_transport_config(WATER_PHANTOM_CONFIG, kind="water_phantom")
    refined = load_transport_config(
        WATER_PHANTOM_CONFIG.replace(
            "[transport.mesh]", "[transport.mesh]\n  refinement_factor = 2\n"),
        kind="water_phantom")
    assert refined.mesh_refinement_factor == 2
    assert resolve_mesh_dimension(refined.mesh_base_dimension, 1, 2) == \
        tuple(2 * d for d in resolve_mesh_dimension(base.mesh_base_dimension, 1))
    # and the identity must actually differ, or a refined run reuses a cache
    assert water_phantom_state_hash(base) != water_phantom_state_hash(refined)
