"""One physical object at many pitches, and what that reveals about the metrics.

The builders in `spatial/synthetic.py` take feature sizes in VOXELS, which is
right for checking a metric against a hand-computed property at one resolution
and useless for asking what a pitch costs: at a different shape they describe a
different object, so a sweep built from them moves the object and the sampling
together. `PhysicalSlab` and `PhysicalSpheres` carry the object in micrometres
so only the sampling moves.

The interesting result is negative, and it is pinned here rather than left in a
report: `specific_interface_area` does not converge to a surface area.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_calibration.spatial.structure import (specific_interface_area,
                                                   summarise)
from biofilm_calibration.spatial.synthetic import PhysicalSlab, PhysicalSpheres

PITCHES = (3.2, 1.6, 0.8, 0.4)


def _rel(got, want):
    return abs(float(got) - float(want)) / abs(float(want))


# ------------------------------------------------- the object stays the object

def test_the_slab_is_the_same_slab_at_every_pitch():
    """The whole premise. If the physical object moved with the pitch, every
    number downstream would be measuring two things at once."""
    spec = PhysicalSlab()
    for pitch in PITCHES:
        B = spec.rasterize(pitch)
        assert B.shape == tuple(round(e / pitch) for e in spec.extent_um)
        assert B.mean() == pytest.approx(spec.true_biovolume_fraction, abs=1e-12)


def test_the_sphere_lattice_keeps_its_count_and_converges_in_volume():
    spec = PhysicalSpheres()
    errors = []
    for pitch in PITCHES:
        B = spec.rasterize(pitch)
        errors.append(_rel(B.mean(), spec.true_biovolume_fraction))
        assert summarise(B, B >= 0.5, pitch).n_components == spec.n_spheres, (
            "the component count is exact at every pitch or the spheres have "
            "merged, which changes the object")
    assert errors[-1] < errors[0] / 5, (
        f"biovolume should converge with pitch; got {errors}")


def test_an_axis_aligned_slab_is_exact_at_every_pitch():
    """Voxelisation of an axis-aligned boundary loses nothing, so this is the
    control: any error seen for the spheres is about curvature, not about the
    ladder machinery."""
    spec = PhysicalSlab()
    for pitch in PITCHES:
        B = spec.rasterize(pitch)
        got = summarise(B, B >= 0.5, pitch)
        assert got.biovolume_fraction == pytest.approx(
            spec.true_biovolume_fraction, abs=1e-12)
        assert got.thickness_mean_um == pytest.approx(spec.thickness_um, abs=1e-12)
        assert got.specific_interface_area_per_um == pytest.approx(
            spec.true_specific_interface_area_per_um, rel=1e-9)


# ------------------------------------- the observable that does not converge

@pytest.mark.parametrize("pitch,expected", [(0.4, 1.5), (0.2, 1.5)])
def test_voxel_interface_area_converges_to_three_halves_of_a_sphere_not_to_one(
        pitch, expected):
    """SPECIFIC INTERFACE AREA IS NOT AN AREA, AND REFINING DOES NOT FIX IT.

    Counting axis-aligned faces measures the integral of |n|_1 over the surface,
    not |n|_2, so the ratio to the true area is E[|n|_1] for the surface's
    orientation distribution — 1 for an axis-aligned plane, sqrt(2) at 45
    degrees, and 3*E|n_x| = 3/2 for an isotropic surface. It is a systematic
    orientation-dependent factor, not a discretisation error, so it survives
    every refinement.

    This matters beyond bookkeeping: `maximum_interface_area_error = 0.10` is
    declared `hard` in the Reference D acceptance config, and NO pitch can meet
    it against a continuum area for a curved interface.
    """
    n = int(round(40 / pitch))
    centre, radius = (n - 1) / 2, n / 4
    idx = np.indices((n, n, n)).astype(float)
    ball = (((idx[0] - centre) ** 2 + (idx[1] - centre) ** 2
             + (idx[2] - centre) ** 2) <= radius ** 2)

    measured = specific_interface_area(ball, pitch)
    true_area_density = 4 * np.pi * (radius * pitch) ** 2 / (n * pitch) ** 3
    assert measured / true_area_density == pytest.approx(expected, rel=0.02)


def test_the_same_estimator_is_exact_for_an_axis_aligned_plane():
    """The other end of the same law: |n|_1 == |n|_2 when the normal is on an
    axis, so the factor is exactly 1. Together with the sphere case this shows
    the bias is orientation, not resolution."""
    n, pitch = 100, 0.2
    B = np.zeros((n, n, n), dtype=bool)
    B[:, :, :n // 2] = True
    assert specific_interface_area(B, pitch) == pytest.approx(
        1.0 / (n * pitch), rel=1e-9)


def test_a_forty_five_degree_plane_lands_on_root_two():
    """The intermediate case, which is what makes the mechanism a law rather
    than two coincidences."""
    n, pitch = 200, 0.2
    idx = np.indices((n, n, n)).astype(float)
    B = (idx[0] + idx[1]) < n
    true_area_density = np.sqrt(2) * (n * pitch) ** 2 / (n * pitch) ** 3
    assert specific_interface_area(B, pitch) / true_area_density == \
        pytest.approx(np.sqrt(2), rel=0.02)


def test_the_bias_does_not_shrink_with_refinement():
    """The claim that makes it a defect rather than a tolerance question."""
    ratios = []
    for pitch in (0.4, 0.2, 0.1):
        n = int(round(20 / pitch))
        centre, radius = (n - 1) / 2, n / 4
        idx = np.indices((n, n, n)).astype(float)
        ball = (((idx[0] - centre) ** 2 + (idx[1] - centre) ** 2
                 + (idx[2] - centre) ** 2) <= radius ** 2)
        true = 4 * np.pi * (radius * pitch) ** 2 / (n * pitch) ** 3
        ratios.append(specific_interface_area(ball, pitch) / true)
    assert min(ratios) > 1.4, f"the bias vanished, which it must not: {ratios}"
    assert abs(ratios[-1] - 1.5) < 0.02
