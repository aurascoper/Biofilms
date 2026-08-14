import numpy as np

from biofilm_openmc.config import load_config
from biofilm_openmc.fingerprint import (label_state_hash, needs_rerun,
                                        transport_state_hash)

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
    hotter = load_config(VALID_CONFIG.replace("photons_per_second = 1.0e9",
                                              "photons_per_second = 2.0e9"))
    assert transport_state_hash(snapshot, hotter, "endfb-viii.0") != h
