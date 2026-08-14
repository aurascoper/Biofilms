"""Morphology metrics, and the shape classes the CPM cannot hold.

The load-bearing test here is that the harness SEPARATES a rod from a sphere of
equal volume. compute_delta_H has no elongation term, so every object relaxes
to a blob; if these metrics could not tell a rod from a blob, that collapse
would be invisible in any comparison against microscopy.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_calibration.spatial.morphology import (all_objects, object_morphology,
                                                    occupancy,
                                                    size_distribution_summary)
from biofilm_calibration.spatial import representability as rep


def test_sphere_volume_and_isotropy():
    arr = rep.make_sphere(radius=6.0)
    m = object_morphology(arr, 1, 1.0)
    # discretized sphere volume approaches 4/3 pi r^3
    assert np.isclose(m.volume_um3, 4 / 3 * np.pi * 6.0 ** 3, rtol=0.06)
    assert m.is_isotropic
    assert m.aspect_ratio_major_minor == pytest.approx(1.0, abs=0.1)
    assert m.connected_components == 1


def test_axis_lengths_recover_a_known_ellipsoid():
    """Inertia-tensor axes, so the answer must not depend on orientation."""
    arr = rep.make_ellipsoid(shape=(40, 40, 40), semi=(12.0, 6.0, 6.0))
    m = object_morphology(arr, 1, 1.0)
    assert m.major_axis_um == pytest.approx(24.0, rel=0.08)
    assert m.minor_axis_um == pytest.approx(12.0, rel=0.08)
    assert m.aspect_ratio_major_minor == pytest.approx(2.0, rel=0.10)


def test_rod_is_distinguishable_from_a_sphere_of_equal_volume():
    """THE test. If this failed, the CPM's collapse of every species to an
    isotropic blob would be undetectable in a morphology comparison."""
    rod = object_morphology(rep.make_rod(), 1, 1.0)
    radius = (3.0 * rod.volume_um3 / (4.0 * np.pi)) ** (1 / 3)
    sphere = object_morphology(
        rep.make_sphere(shape=(32, 32, 32), radius=radius), 1, 1.0)

    assert np.isclose(rod.volume_um3, sphere.volume_um3, rtol=0.10)
    assert not rod.is_isotropic
    assert sphere.is_isotropic
    assert rod.aspect_ratio_major_minor > 2.0 * sphere.aspect_ratio_major_minor


def test_disconnected_object_is_detected():
    """No connectivity term exists in the energy, so a split object is neither
    prevented nor penalised — the harness must at least see it."""
    m = object_morphology(rep.make_disconnected(), 1, 1.0)
    assert m.connected_components == 2


def test_voxel_size_scales_volume_cubically_and_axes_linearly():
    arr = rep.make_sphere(radius=5.0)
    a = object_morphology(arr, 1, 1.0)
    b = object_morphology(arr, 1, 0.5)
    assert np.isclose(b.volume_um3, a.volume_um3 / 8.0)
    assert np.isclose(b.major_axis_um, a.major_axis_um / 2.0)


def test_anisotropic_voxels_are_honoured():
    """Confocal z-spacing is routinely coarser than xy; ignoring that would
    silently elongate every object along z."""
    arr = rep.make_sphere(radius=5.0)
    m = object_morphology(arr, 1, (0.1, 0.1, 0.5))
    assert m.major_axis_um > m.minor_axis_um * 2

    # order matters: the long axis must be the one with the coarse spacing
    assert m.aspect_ratio_major_minor == pytest.approx(5.0, rel=0.15)


def test_single_voxel_has_no_resolvable_shape():
    arr = np.zeros((5, 5, 5), dtype=np.int32)
    arr[2, 2, 2] = 1
    m = object_morphology(arr, 1, 1.0)
    assert m.n_voxels == 1
    assert m.major_axis_um == 0.0
    assert not np.isfinite(m.aspect_ratio_major_minor)


def test_wall_sites_are_not_objects():
    """The CPM uses -1 for 'outside the biological domain'. It is a wall
    marker, and counting it as an object would invent a huge phantom cell."""
    arr = np.zeros((8, 8, 8), dtype=np.int32)
    arr[:2] = -1
    arr[4:6, 4:6, 4:6] = 7
    objs = all_objects(arr, 1.0)
    assert [o.label for o in objs] == [7]


def test_occupancy_is_a_dimensionless_fraction_of_non_wall_sites():
    arr = np.zeros((10, 10, 10), dtype=np.int32)
    arr[0] = -1                       # 100 wall sites
    arr[5:7, :, :] = 3                # 200 occupied
    assert occupancy(arr) == pytest.approx(200 / 900)
    assert 0.0 <= occupancy(arr) <= 1.0


def test_boundary_objects_are_excluded_from_the_size_distribution():
    """Objects clipped by the image edge are truncated and would bias a size
    distribution downward."""
    arr = np.zeros((20, 20, 20), dtype=np.int32)
    arr[0:3, 0:3, 0:3] = 1            # touches the edge
    arr[8:12, 8:12, 8:12] = 2         # interior
    objs = all_objects(arr, 1.0)
    summary = size_distribution_summary(objs)
    assert summary["n"] == 1
    assert summary["n_excluded_boundary"] == 1


def test_representability_audit_matches_the_hamiltonian():
    """Every anisotropic or disconnected class must be both detectable by the
    harness and unmaintainable by the CPM."""
    rows = {r["shape_class"]: r for r in rep.audit()}
    for name in ("rod", "filament", "disconnected"):
        assert rows[name]["detectable_as_non_blob"], name
        assert not rows[name]["maintainable_by_cpm"], name
    for name in ("sphere", "biomass_parcel"):
        assert rows[name]["maintainable_by_cpm"], name


def test_bad_input_is_refused():
    with pytest.raises(ValueError, match="3-D"):
        object_morphology(np.zeros((4, 4)), 1, 1.0)
    with pytest.raises(ValueError, match="not present"):
        object_morphology(np.zeros((4, 4, 4), dtype=np.int32), 9, 1.0)
    with pytest.raises(ValueError, match="positive"):
        object_morphology(rep.make_sphere(), 1, (1.0, 0.0, 1.0))
