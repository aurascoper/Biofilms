"""Biomass-field coarse-graining, occupancy maps, and structural observables.

Validated against synthetic fields with analytic truth, so a failure says the
metric is wrong rather than that the data moved. The load-bearing tests are
that coarse-graining CONSERVES the biovolume fraction (it is a mean, not a
sum), and that the two occupancy maps differ in exactly the way that makes the
choice a calibration decision.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_calibration.spatial import structure, synthetic
from biofilm_calibration.spatial.field import (FieldError, apply_occupancy,
                                               biovolume_fraction,
                                               coarse_grain,
                                               occupancy_biovolume_error,
                                               occupancy_mass_preserving,
                                               occupancy_threshold)
from biofilm_calibration.spatial.pitch_selection import (OBSERVABLE_TOLERANCES,
                                                         evaluate_pitch,
                                                         select_pitch)

TOL = {k: 0.05 for k in OBSERVABLE_TOLERANCES.values()}


# --- coarse-graining -------------------------------------------------------

def test_coarse_graining_conserves_the_biovolume_fraction():
    """THE test. A biomass volume fraction is INTENSIVE, so coarsening is a
    block mean. The transport package's coarsen_field block-SUMS deliberately
    because energy and mass are extensive; using that here would be off by
    factor**3."""
    B = synthetic.correlated_field(shape=(32, 32, 32), seed=1)
    fine = biovolume_fraction(B)
    for f in (2, 4, 8):
        assert biovolume_fraction(coarse_grain(B, f)) == pytest.approx(fine, rel=1e-12)


def test_a_block_sum_would_be_wrong_by_factor_cubed():
    """Stated as a test so the distinction cannot quietly rot."""
    B = synthetic.slab(shape=(8, 8, 8), thickness_voxels=4)
    f = 2
    mean = coarse_grain(B, f)
    summed = B.reshape(4, f, 4, f, 4, f).sum(axis=(1, 3, 5))
    assert np.allclose(summed, mean * f ** 3)


def test_slab_biovolume_is_exact():
    B = synthetic.slab(shape=(16, 16, 32), thickness_voxels=8)
    assert biovolume_fraction(B) == pytest.approx(8 / 32)


def test_coarsening_refuses_a_non_dividing_factor():
    B = synthetic.slab(shape=(16, 16, 32))
    with pytest.raises(FieldError, match="does not divide"):
        coarse_grain(B, 5)


def test_field_must_be_a_fraction():
    with pytest.raises(FieldError, match=r"\[0, 1\]"):
        biovolume_fraction(np.full((4, 4, 4), 2.5))
    with pytest.raises(FieldError, match="3-D"):
        biovolume_fraction(np.zeros((4, 4)))
    with pytest.raises(FieldError, match="NaN"):
        biovolume_fraction(np.full((4, 4, 4), np.nan))


# --- occupancy maps --------------------------------------------------------

def test_mass_preserving_conserves_biovolume_and_threshold_does_not():
    """The two maps differ in what they conserve, which is exactly why one has
    to be declared rather than assumed."""
    phi = synthetic.correlated_field(shape=(48, 48, 48), seed=2,
                                     target_fraction=0.35)
    mass_err = occupancy_biovolume_error(phi, "mass_preserving", seed=7)
    assert abs(mass_err) < 0.02, mass_err

    # a threshold at 0.5 discards every partially-filled bin below it
    thr_err = occupancy_biovolume_error(phi, "threshold", tau=0.5)
    assert abs(thr_err) > abs(mass_err)


def test_threshold_is_deterministic_and_mass_preserving_is_seeded():
    phi = synthetic.correlated_field(shape=(24, 24, 24), seed=3)
    assert np.array_equal(occupancy_threshold(phi, 0.5),
                          occupancy_threshold(phi, 0.5))
    assert np.array_equal(occupancy_mass_preserving(phi, 11),
                          occupancy_mass_preserving(phi, 11))
    assert not np.array_equal(occupancy_mass_preserving(phi, 11),
                              occupancy_mass_preserving(phi, 12))


def test_occupancy_refuses_when_undeclared():
    phi = synthetic.slab()
    with pytest.raises(FieldError, match="no occupancy mapping declared"):
        apply_occupancy(phi, "")
    with pytest.raises(FieldError, match="requires threshold_tau"):
        apply_occupancy(phi, "threshold")
    with pytest.raises(FieldError, match="requires a declared seed"):
        apply_occupancy(phi, "mass_preserving", seed=None)
    with pytest.raises(FieldError, match="unknown occupancy mapping"):
        apply_occupancy(phi, "vibes")


# --- structural observables ------------------------------------------------

def test_porosity_is_the_complement_of_biovolume():
    B = synthetic.slab(shape=(16, 16, 32), thickness_voxels=8)
    assert structure.porosity(B) == pytest.approx(1 - 8 / 32)


def test_slab_interface_area_is_analytic():
    """A slab exposes exactly one INTERNAL face per column — its top. The
    array boundary is not an interface."""
    shape, t = (10, 10, 20), 5
    B = synthetic.slab(shape=shape, thickness_voxels=t)
    got = structure.specific_interface_area(B > 0.5, 1.0)
    expected = (shape[0] * shape[1]) / float(np.prod(shape))
    assert got == pytest.approx(expected, rel=1e-9)


def test_interface_area_does_not_depend_on_the_crop():
    """Biomass reaching the edge of the field of view is cut off by the IMAGE,
    not by void. Counting the array boundary would make the metric a property
    of the crop rather than of the biofilm."""
    big = synthetic.slab(shape=(20, 20, 20), thickness_voxels=5)
    small = big[:10, :10, :]          # same structure, half the field of view
    assert structure.specific_interface_area(big > 0.5, 1.0) == pytest.approx(
        structure.specific_interface_area(small > 0.5, 1.0), rel=1e-9)


def test_thickness_recovers_the_slab():
    B = synthetic.slab(shape=(12, 12, 30), thickness_voxels=9)
    thick = structure.thickness_profile(B > 0.5, (0.2, 0.2, 0.5))
    assert thick.shape == (12, 12)
    assert np.allclose(thick, 9 * 0.5)


def test_thickness_counts_biomass_not_extent():
    """An internal void must reduce the thickness, not be counted as biomass."""
    B = synthetic.slab(shape=(4, 4, 10), thickness_voxels=8)
    B[:, :, 3:5] = 0.0                     # a void inside the slab
    thick = structure.thickness_profile(B > 0.5, 1.0)
    assert np.allclose(thick, 6.0)


def test_component_sizes_recover_the_sphere_count():
    B = synthetic.spheres_on_substrate()
    sizes = structure.components_26(B > 0.5)
    assert len(sizes) == synthetic.n_spheres()
    assert min(sizes) > 0 and max(sizes) == sizes[0]     # sorted descending


def test_correlation_length_recovers_the_synthetic_scale():
    """The field is smoothed to a known scale, so the measured decorrelation
    must track it — a pitch comparable to this erases structure rather than
    resolving it."""
    short = synthetic.correlated_field(shape=(64, 64, 64),
                                       correlation_length_voxels=4.0, seed=4)
    long = synthetic.correlated_field(shape=(64, 64, 64),
                                      correlation_length_voxels=16.0, seed=4)
    l_short = structure.correlation_length(short, 1.0)
    l_long = structure.correlation_length(long, 1.0)
    assert np.isfinite(l_short) and np.isfinite(l_long)
    assert l_long > 2.0 * l_short, (l_short, l_long)


def test_depth_profile_recovers_a_gradient():
    B = synthetic.graded_field(shape=(8, 8, 16), top=0.9, bottom=0.1)
    prof = structure.depth_profile(B)
    assert prof[0] == pytest.approx(0.9)
    assert prof[-1] == pytest.approx(0.1)
    assert np.all(np.diff(prof) < 0)


def test_chord_lengths_drop_edge_truncated_runs():
    """Runs clipped by the field of view are truncated by the image, not by
    the structure, and would bias the distribution short."""
    mask = np.zeros((1, 1, 10), dtype=bool)
    mask[0, 0, 0:3] = True        # touches the start: dropped
    mask[0, 0, 5:8] = True        # interior: kept
    chords = structure.chord_lengths(mask, 1.0, axis=2)
    assert chords.tolist() == [3.0]


def test_summarise_returns_every_declared_observable():
    B = synthetic.correlated_field(shape=(32, 32, 32), seed=5)
    s = structure.summarise(B, B > 0.5, 0.5)
    for key in OBSERVABLE_TOLERANCES:
        assert key in s.as_dict(), key


# --- pitch selection -------------------------------------------------------

def _evals(B, factors, tolerances, held_out_factor_ok=True):
    out = []
    for f in factors:
        for sample, held in (("train1", False), ("holdout1", True)):
            out.append(evaluate_pitch(B, 0.25, f, sample_id=sample,
                                      held_out=held,
                                      occupancy_mapping="threshold",
                                      tolerances=tolerances, tau=0.5))
    return out


def test_select_refuses_without_declared_tolerances():
    B = synthetic.correlated_field(shape=(32, 32, 32), seed=6)
    evals = _evals(B, (1, 2), TOL)
    with pytest.raises(FieldError, match="refusing to select"):
        select_pitch(evals, {})
    with pytest.raises(FieldError, match="incomplete"):
        select_pitch(evals, {"maximum_porosity_error": 0.05})


def test_select_refuses_without_a_held_out_stack():
    B = synthetic.correlated_field(shape=(32, 32, 32), seed=6)
    train_only = [e for e in _evals(B, (1, 2), TOL) if not e.held_out]
    with pytest.raises(FieldError, match="has not been validated"):
        select_pitch(train_only, TOL)


def test_factor_one_always_passes_and_coarser_eventually_fails():
    """Sanity for the rule: no coarsening is exact, and enough coarsening
    destroys the structure."""
    B = synthetic.correlated_field(shape=(32, 32, 32),
                                   correlation_length_voxels=4.0, seed=7)
    evals = _evals(B, (1, 16), TOL)
    at1 = [e for e in evals if e.factor == 1]
    at16 = [e for e in evals if e.factor == 16]
    assert all(e.passed for e in at1)
    assert all(max(e.errors.values()) > 0 for e in at16)


def test_select_takes_the_coarsest_qualifying_factor():
    B = synthetic.correlated_field(shape=(32, 32, 32), seed=8)
    evals = _evals(B, (1, 2, 16), TOL)
    chosen = select_pitch(evals, TOL)
    assert chosen is not None
    # never coarser than something that failed
    failed = {e.factor for e in evals if not e.passed}
    assert chosen not in failed


def test_select_returns_none_when_nothing_qualifies():
    """A real answer, not an error: this field admits no acceptable coarsening
    at these tolerances."""
    B = synthetic.correlated_field(shape=(32, 32, 32),
                                   correlation_length_voxels=2.0, seed=9)
    strict = {k: 1e-6 for k in OBSERVABLE_TOLERANCES.values()}
    evals = _evals(B, (2, 4), strict)
    assert select_pitch(evals, strict) is None


def test_an_undeclared_tolerance_cannot_be_passed():
    B = synthetic.slab(shape=(16, 16, 16), thickness_voxels=8)
    partial = dict(TOL)
    partial.pop("maximum_porosity_error")
    e = evaluate_pitch(B, 0.25, 1, sample_id="s", held_out=True,
                       occupancy_mapping="threshold", tolerances=partial, tau=0.5)
    assert not e.passed
