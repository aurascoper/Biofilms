"""Scale identifiability and the refusal to select.

The two things this file pins: a pitch is not identifiable from a volume
alone, and the engine will not choose one without declared thresholds.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_calibration.spatial.scale_candidates import (
    domain_from_pitch, enumerate_candidates, grid_memory_bytes, select,
    sites_for_volume)

ENTITIES = {
    "coccus": {"volume_um3": 1.0, "minor_axis_um": 1.2, "major_axis_um": 1.2},
    "yeast": {"volume_um3": 65.4, "minor_axis_um": 5.0, "major_axis_um": 5.0},
}


def test_pitch_and_site_count_are_not_separately_identifiable():
    """V_physical = V_sites * a^3 constrains only the product. Many (a, sites)
    pairs reproduce the same physical volume exactly, so a volume measurement
    alone cannot pick one."""
    volume = 65.4
    pairs = [(a, sites_for_volume(volume, a)) for a in (0.2, 0.5, 0.82, 1.0)]
    for a, sites in pairs:
        assert np.isclose(sites * a ** 3, volume)
    # ... and they disagree wildly about how many sites a cell is
    counts = [s for _, s in pairs]
    assert max(counts) / min(counts) > 100


def test_the_famous_shortcut_is_an_accident_not_a_calibration():
    """Dividing a cell volume by the incumbent V_target = 120 and taking the
    cube root yields a pitch — but the 120 was never chosen physically, so the
    pitch inherits that arbitrariness."""
    a_from_shortcut = (65.4 / 120) ** (1 / 3)
    assert a_from_shortcut == pytest.approx(0.82, abs=0.01)
    # a different, equally unjustified site target moves the pitch a lot
    a_alt = (65.4 / 1000) ** (1 / 3)
    assert a_alt / a_from_shortcut == pytest.approx(0.4914, abs=0.01)


def test_domain_is_forced_to_L_equals_2R():
    """R = N/2 is structural in biofilms_potts.jl; there is no length field."""
    for a in (0.1, 0.5, 2.0):
        for N in (30, 60, 120):
            r, L = domain_from_pitch(a, N)
            assert np.isclose(L, 2 * r)


def test_no_candidate_claims_full_apparatus_compatibility():
    cands = enumerate_candidates(pitches_um=[0.5, 1.0], grid_sizes=[60],
                                 entities=ENTITIES)
    assert cands and not any(c.full_apparatus_compatible for c in cands)
    assert all(c.representative_segment_compatible for c in cands)


def test_status_is_unresolved_without_a_declared_threshold():
    cands = enumerate_candidates(pitches_um=[1.0], grid_sizes=[60],
                                 entities=ENTITIES)
    assert {c.representability_status for c in cands} == {"unresolved"}


def test_status_responds_to_a_declared_threshold():
    cands = enumerate_candidates(pitches_um=[0.2, 2.0], grid_sizes=[60],
                                 entities=ENTITIES, min_axis_voxels=4.0)
    by = {(c.pitch_um, c.species): c.representability_status for c in cands}
    # a 1.2 um coccus is 6 voxels across at 0.2 um, and under 1 at 2.0 um
    assert by[(0.2, "coccus")] == "representable"
    assert by[(2.0, "coccus")] == "not_representable"


def test_select_refuses_without_thresholds():
    cands = enumerate_candidates(pitches_um=[0.5], grid_sizes=[60],
                                 entities=ENTITIES, min_axis_voxels=4.0)
    with pytest.raises(ValueError, match="refusing to select"):
        select(cands, None)
    with pytest.raises(ValueError, match="refusing to select"):
        select(cands, {})
    with pytest.raises(ValueError, match="incomplete"):
        select(cands, {"minimum_minor_axis_voxels": 4.0})


def test_select_applies_every_declared_threshold():
    cands = enumerate_candidates(pitches_um=[0.2, 1.0], grid_sizes=[60, 200],
                                 entities=ENTITIES, min_axis_voxels=4.0)
    chosen = select(cands, {"minimum_minor_axis_voxels": 4.0,
                            "maximum_volume_quantization_error": 1e-9,
                            "maximum_memory_bytes": grid_memory_bytes(64)})
    assert chosen, "expected at least one admissible candidate"
    for c in chosen:
        assert c.minor_axis_voxels >= 4.0
        assert c.estimated_grid_memory_bytes <= grid_memory_bytes(64)
    # N=200 is over the memory ceiling and must be gone
    assert {c.N for c in chosen} == {60}


def test_memory_grows_as_N_cubed():
    assert grid_memory_bytes(120) == 8 * grid_memory_bytes(60)


def test_quantization_error_reported_against_a_declared_target():
    cands = enumerate_candidates(pitches_um=[0.82], grid_sizes=[60],
                                 entities={"yeast": ENTITIES["yeast"]},
                                 target_volume_sites=120)
    c = cands[0]
    assert c.target_volume_sites == 120
    assert c.volume_quantization_error < 0.05      # 0.82 um was reverse-engineered
    # with the same pitch but a different declared target, the cost is visible
    worse = enumerate_candidates(pitches_um=[0.82], grid_sizes=[60],
                                 entities={"yeast": ENTITIES["yeast"]},
                                 target_volume_sites=500)[0]
    assert worse.volume_quantization_error > 3.0
