import numpy as np

from biofilm_openmc.config import load_config
from biofilm_openmc.fingerprint import (dose_state_hash, label_state_hash,
                                        needs_rerun, transport_state_hash)

from conftest import VALID_CONFIG


def test_identical_state_reuses(snapshot, config):
    h = transport_state_hash(snapshot, config, "endfb-viii.0")
    assert h == transport_state_hash(snapshot, config, "endfb-viii.0")
    assert needs_rerun(h, snapshot, config, "endfb-viii.0") is False
    assert needs_rerun(None, snapshot, config, "endfb-viii.0") is True


def test_label_only_change_reuses_field_but_changes_label_hash(snapshot, config):
    h_before = transport_state_hash(snapshot, config)
    l_before = label_state_hash(snapshot)
    # relabel lineages/generations only — physical state untouched
    snapshot.lineage_id[snapshot.cell_id > 0] += 100
    snapshot.generation[snapshot.cell_id > 0] += 1
    assert transport_state_hash(snapshot, config) == h_before
    assert needs_rerun(h_before, snapshot, config) is False
    assert label_state_hash(snapshot) != l_before


def test_composition_change_forces_rerun(snapshot, config):
    h = transport_state_hash(snapshot, config)
    denser = load_config(VALID_CONFIG.replace("density_g_cm3 = 1.1",
                                              "density_g_cm3 = 1.2"))
    assert needs_rerun(h, snapshot, denser) is True


def test_occupancy_change_forces_rerun_no_threshold(snapshot, config):
    # a SINGLE voxel changing class must force a rerun — exact policy,
    # no churn threshold (review amendment 7)
    h = transport_state_hash(snapshot, config)
    x = np.argwhere(snapshot.cell_id == 0)[0]
    snapshot.cell_id[tuple(x)] = 99
    assert needs_rerun(h, snapshot, config) is True


def test_source_and_data_identity_in_hash(snapshot, config):
    h = transport_state_hash(snapshot, config, "endfb-viii.0")
    assert transport_state_hash(snapshot, config, "endfb-viii.1") != h
    softer = load_config(VALID_CONFIG.replace("spectrum_energies_eV = [1.0e6]",
                                              "spectrum_energies_eV = [6.6e5]"))
    assert transport_state_hash(snapshot, softer, "endfb-viii.0") != h


def test_activity_changes_dose_identity_but_not_transport(snapshot, config):
    """Heating is tallied per source particle, so the activity cannot change
    the transport result and must not invalidate it — but it does change the
    Gy/s field, so the dose identity must move."""
    h = transport_state_hash(snapshot, config, "endfb-viii.0")
    hotter = load_config(VALID_CONFIG.replace("photons_per_second = 1.0e9",
                                              "photons_per_second = 2.0e9"))
    assert transport_state_hash(snapshot, hotter, "endfb-viii.0") == h
    assert needs_rerun(h, snapshot, hotter, "endfb-viii.0") is False
    assert dose_state_hash(h, hotter.photons_per_second) != \
        dose_state_hash(h, config.photons_per_second)


def test_coarsening_changes_the_transport_identity(snapshot, config):
    """Changing the tally resolution changes the run, so a coarse result must
    never be served from a fine run's cache."""
    from biofilm_openmc.config import load_transport_config

    coarse = load_transport_config(VALID_CONFIG.replace(
        "[transport]", "[transport]\n  [transport.mesh]\n  coarsening_factor = 2\n"))
    assert transport_state_hash(snapshot, coarse) != \
        transport_state_hash(snapshot, config)


def test_water_phantom_hash_is_snapshot_free_and_mesh_sensitive():
    from biofilm_openmc.config import WATER_PHANTOM, load_transport_config
    from biofilm_openmc.fingerprint import water_phantom_state_hash

    from conftest import WATER_PHANTOM_CONFIG

    base = load_transport_config(WATER_PHANTOM_CONFIG, kind=WATER_PHANTOM)
    h = water_phantom_state_hash(base, "endfb-viii.0")
    assert h == water_phantom_state_hash(base, "endfb-viii.0")
    assert water_phantom_state_hash(base, "endfb-viii.1") != h

    coarse = load_transport_config(
        WATER_PHANTOM_CONFIG.replace("base_dimension = [12, 12, 12]",
                                     "base_dimension = [12, 12, 12]\n  coarsening_factor = 2"),
        kind=WATER_PHANTOM)
    assert water_phantom_state_hash(coarse, "endfb-viii.0") != h
