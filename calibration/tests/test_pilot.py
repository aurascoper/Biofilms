"""The public-data pilot path, tested without the public data.

The parts that can go wrong silently — axis selection, intensity rescaling, the
temporal statistic — are exercised on synthetic arrays. The parts that need a
real ND2 skip cleanly, marked `needs_measured_data`, so the report distinguishes
"absent" from "passed".
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from biofilm_calibration.spatial import synthetic
from biofilm_calibration.spatial.field import FieldError
from biofilm_calibration.spatial.ingest import ingest
from biofilm_calibration.spatial.readers import SMOKE, _select_volume, _to_fraction


# --- axis handling: the part that silently ruins a volume ------------------

def test_axes_are_selected_by_name_not_position():
    """ND2 axis order is metadata, not convention. Guessing by position is how
    a time index becomes a z index and the stack quietly becomes nonsense."""
    # (T, C, Z, Y, X) — the common confocal layout
    arr = np.zeros((3, 2, 4, 5, 6))
    arr[1, 1] = 7.0                      # only t=1, c=1 is hot
    vol = _select_volume(arr, ["T", "C", "Z", "Y", "X"], channel=1, timepoint=1)
    assert vol.shape == (6, 5, 4)        # returned as logical (x, y, z)
    assert (vol == 7.0).all()

    # a different declared order must give a different, correct answer
    arr2 = np.zeros((2, 3, 6, 5, 4))     # (C, T, X, Y, Z)
    arr2[1, 1] = 7.0
    vol2 = _select_volume(arr2, ["C", "T", "X", "Y", "Z"], channel=1, timepoint=1)
    assert vol2.shape == (6, 5, 4)
    assert (vol2 == 7.0).all()


def test_a_time_axis_without_a_timepoint_is_refused():
    """Silently taking t=0 from a time course would discard the dynamics the
    selected observable depends on."""
    arr = np.zeros((3, 4, 5, 6))
    with pytest.raises(FieldError, match="pass timepoint="):
        _select_volume(arr, ["T", "Z", "Y", "X"], channel=0, timepoint=None)


def test_a_volume_with_no_z_axis_is_refused():
    arr = np.zeros((5, 6))
    with pytest.raises(FieldError, match="expected a 3-D volume"):
        _select_volume(arr, ["Y", "X"], channel=0, timepoint=None)


# --- intensity to biomass fraction -----------------------------------------

def test_threshold_gives_a_mask_and_smoke_gives_a_probability_field():
    from biofilm_calibration.spatial.readers import SMOKE

    v = np.array([[[0.0, 10.0], [20.0, 30.0]]])
    masked = _to_fraction(v, ("fixed_threshold", 15.0, float("nan")))
    assert set(np.unique(masked)) == {0.0, 1.0}
    assert masked.sum() == 2

    scaled = _to_fraction(v, SMOKE)
    assert scaled.min() == 0.0 and scaled.max() == 1.0
    assert len(np.unique(scaled)) > 2     # not a mask


def test_a_flat_volume_cannot_be_rescaled():
    from biofilm_calibration.spatial.readers import SMOKE

    with pytest.raises(FieldError, match="no dynamic range"):
        _to_fraction(np.ones((2, 2, 2)), SMOKE)


def test_a_declared_scale_is_applied_as_given_not_re_estimated():
    """The point of fixed_background_and_scale: two volumes with different
    dynamic ranges must map onto the SAME scale, which per-volume min-max
    cannot do."""
    from biofilm_calibration.spatial.readers import SMOKE

    dim = np.array([[[0.0, 25.0], [50.0, 75.0]]])
    bright = dim * 2.0
    t = ("fixed_background_and_scale", 0.0, 200.0)
    assert _to_fraction(bright, t).mean() > _to_fraction(dim, t).mean()
    # ... whereas the smoke path erases exactly that difference
    assert _to_fraction(bright, SMOKE).mean() == _to_fraction(dim, SMOKE).mean()


def test_an_unknown_transform_is_refused():
    with pytest.raises(FieldError, match="unknown transform"):
        _to_fraction(np.ones((2, 2, 2)), ("otsu", 0.0, 1.0))


# --- the pilot's own passes ------------------------------------------------

def test_spatial_pass_skips_non_dividing_factors_rather_than_failing():
    """Nesting is required for cross-resolution comparison, so a factor that
    does not divide the stack is skipped, not forced."""
    from vcholerae_pilot import spatial_pass
    from biofilm_calibration.spatial.pitch_selection import OBSERVABLE_TOLERANCES

    B = synthetic.correlated_field(shape=(24, 24, 24), seed=1)
    f = ingest(B, 0.2, "cells_only", sample_id="s", source_id="SRC")
    tol = {k: 0.05 for k in OBSERVABLE_TOLERANCES.values()}
    rows = spatial_pass(f, [1, 2, 4, 5, 16], occupancy_mapping="threshold",
                        tau=0.5, seed=1, tolerances=tol, transform=SMOKE)
    factors = [r["coarsening_factor"] for r in rows]
    assert 5 not in factors and 16 not in factors     # 24 % 5, 24 % 16
    assert factors == [1, 2, 4]


def test_every_pilot_row_is_labelled_a_surrogate():
    """The label must survive being copied out of context."""
    from vcholerae_pilot import REFERENCE_SYSTEM, spatial_pass
    from biofilm_calibration.spatial.pitch_selection import OBSERVABLE_TOLERANCES

    B = synthetic.correlated_field(shape=(16, 16, 16), seed=2)
    f = ingest(B, 0.2, "cells_only", sample_id="s", source_id="SRC")
    rows = spatial_pass(f, [1, 2], occupancy_mapping="threshold", tau=0.5,
                        seed=1, transform=SMOKE,
                        tolerances={k: 0.05 for k in OBSERVABLE_TOLERANCES.values()})
    assert rows
    for r in rows:
        assert r["target_calibration"] == "false"
        assert r["reference_system_id"] == REFERENCE_SYSTEM == "public_vcholerae_surrogate"


def test_temporal_pass_recovers_a_known_decorrelation():
    """The selected observable, on a series that decorrelates by construction:
    a field blended from an initial state toward independent noise."""
    from vcholerae_pilot import temporal_pass

    rng = np.random.default_rng(0)
    base = synthetic.correlated_field(shape=(16, 16, 16), seed=3)
    fields = []
    for t in range(8):
        w = np.exp(-t / 2.0)                     # correlation with t=0 decays
        other = synthetic.correlated_field(shape=(16, 16, 16), seed=100 + t)
        mixed = np.clip(w * base + (1 - w) * other, 0.0, 1.0)
        fields.append(ingest(mixed, 0.2, "cells_only",
                             sample_id=f"s_t{t}", source_id="SRC"))

    out = temporal_pass(fields, (0.2, 0.2, 0.2), transform=SMOKE)
    assert out["temporal_observable_computable"]
    assert out["n_timepoints"] == 8
    assert out["correlation"][0] == pytest.approx(1.0)
    assert out["correlation"][-1] < out["correlation"][0]

    # the SELECTED observable is the primary result and is per-frame
    assert out["selected_observable"] == "biomass_autocorrelation_decay"
    assert len(out["spatial_correlation_length_um_by_t"]) == 8

    # the substituted one is present but cannot be mistaken for the selection
    assert out["secondary_observable"] == "biomass_temporal_autocorrelation_decay"
    assert out["secondary_observable_selected"] == "false"

    # a smoke transform can never claim temporal calibration eligibility
    assert out["temporal_calibration_eligible"] == "false"
    assert out["seconds_per_mcs_eligible"] == "false"
    assert json.dumps(out)                       # NaN would not be valid JSON


def test_temporal_pass_reports_a_static_series_rather_than_dividing_by_zero():
    from vcholerae_pilot import temporal_pass
    B = synthetic.slab(shape=(8, 8, 8), thickness_voxels=4)
    fields = [ingest(B, 0.2, "cells_only", sample_id=f"s{t}", source_id="S")
              for t in range(3)]
    out = temporal_pass(fields, (0.2, 0.2, 0.2), transform=SMOKE)
    assert not out["temporal_observable_computable"]
    assert out["reference_system_id"] == "public_vcholerae_surrogate"


# --- the real file ---------------------------------------------------------

_DATA = Path.home() / "biofilm-pilot-data"
# Dryad doi:10.5061/dryad.zcrjdfnph. Preflight 2026-08-15 established that the
# three files differ in exactly the ways that matter, so each one pins a
# different guard:
#   Dataset_1  T=47, Y, X, S=3   uncalibrated RGB brightfield time-lapse
#   Dataset_2  T=49, C=3, Y, X   calibrated 2-D time-lapse (NOT a stack)
#   Dataset_3  Z=94, C=3, Y, X   calibrated 3-D stack, no time axis
_STACK = _DATA / "Dataset_3_Fig3c_VPS-Bap1.nd2"
_PLANAR = _DATA / "Dataset_2_Fig4a_Timecourse.nd2"
_UNCAL = _DATA / "Dataset_1_Fig1b_Timecourse.nd2"

_why = ("not present; Dryad serves it behind a proof-of-work anti-bot "
        "challenge, so fetch it with a browser")


@pytest.mark.needs_measured_data
@pytest.mark.skipif(not _STACK.exists(), reason=f"{_STACK.name} {_why}")
def test_real_stack_ingests_with_its_own_calibration():
    """When the file exists: the voxel size must come from the FILE. A
    hand-typed calibration is the easiest thing to lose in a conversion."""
    pytest.importorskip("nd2")
    from biofilm_calibration.spatial.readers import read_nd2

    f = read_nd2(_STACK, segmentation_basis="cells_only", sample_id="vc",
                 source_id="DRYAD_ZCRJDFNPH", channel=0)
    assert len(f.shape) == 3
    assert all(v > 0 for v in f.voxel_size_um)
    assert f.reader == "nd2"
    # read from the file, not defaulted to an isotropic guess
    assert f.voxel_size_um[2] != f.voxel_size_um[0]


@pytest.mark.needs_measured_data
@pytest.mark.skipif(not _UNCAL.exists(), reason=f"{_UNCAL.name} {_why}")
def test_an_uncalibrated_file_is_refused_rather_than_read_as_one_micron():
    """The defect this guard exists for, in a real file. nd2 returns
    (1.0, 1.0, 1.0) for a file carrying no calibration, which a `v > 0` check
    passes; every pitch downstream would then be 1.0 x factor with nothing
    behind it."""
    pytest.importorskip("nd2")
    from biofilm_calibration.spatial.readers import read_nd2

    with pytest.raises(FieldError, match="axesCalibrated"):
        read_nd2(_UNCAL, segmentation_basis="cells_only", sample_id="vc",
                 source_id="DRYAD_ZCRJDFNPH", channel=0, timepoint=0)


@pytest.mark.needs_measured_data
@pytest.mark.skipif(not _PLANAR.exists(), reason=f"{_PLANAR.name} {_why}")
def test_a_single_plane_time_lapse_is_not_a_stack():
    """A biomass field is 3-D by contract. Dataset_2 is 1200x1201 over 49
    timepoints with no Z axis, so it cannot enter the harness however many
    frames it has."""
    pytest.importorskip("nd2")
    from biofilm_calibration.spatial.readers import read_nd2

    with pytest.raises(FieldError, match="expected a 3-D volume"):
        read_nd2(_PLANAR, segmentation_basis="cells_only", sample_id="vc",
                 source_id="DRYAD_ZCRJDFNPH", channel=0, timepoint=0)


@pytest.mark.needs_measured_data
@pytest.mark.skipif(not _STACK.exists(), reason=f"{_STACK.name} {_why}")
def test_three_channels_force_an_explicit_choice():
    """Which of 488nm / 561nm / 640nm is the biomass channel is an
    experimental judgement. Index 0 is not an answer to it."""
    pytest.importorskip("nd2")
    from biofilm_calibration.spatial.readers import read_nd2

    with pytest.raises(FieldError, match="no choice was made"):
        read_nd2(_STACK, segmentation_basis="cells_only", sample_id="vc",
                 source_id="DRYAD_ZCRJDFNPH")


@pytest.mark.needs_measured_data
@pytest.mark.skipif(not _STACK.exists(), reason=f"{_STACK.name} {_why}")
def test_a_stack_with_no_time_axis_raises_no_time_axis_error_by_name():
    """So a caller can tell 'this is a single volume' from every other reason
    a read fails, instead of catching Exception and retrying a broken read."""
    pytest.importorskip("nd2")
    from biofilm_calibration.spatial.readers import (NoTimeAxisError,
                                                     read_nd2_timecourse)

    with pytest.raises(NoTimeAxisError):
        read_nd2_timecourse(_STACK, segmentation_basis="cells_only",
                            sample_id="vc", source_id="DRYAD_ZCRJDFNPH",
                            channel=0)
