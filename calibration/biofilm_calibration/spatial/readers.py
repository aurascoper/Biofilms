"""Optional microscopy readers. NOT imported by the core.

Everything here needs the `dataset` extra (`pip install -e "calibration[dataset]"`).
`ingest.py` defines the format-agnostic contract; this module only turns a file
into the three things that contract requires — an array, a voxel calibration,
and a declared segmentation basis.

THE BASIS IS NEVER READ FROM THE FILE. A microscope records what a channel
detected, not what the biologist intends it to represent. Whether a stain marks
cells only, cells plus matrix, or the whole biofilm envelope is a judgement
about the experiment, so it is always a caller-supplied argument here — and it
is the argument that decides whether the resulting volume may be divided into a
measured mass.

THE CALIBRATION IS ALWAYS READ FROM THE FILE, AND ITS ABSENCE IS NOT A ZERO.
`nd2.ND2File.voxel_size()` returns `(1.0, 1.0, 1.0)` for a file that carries no
spatial calibration — it does not raise, and it does not return zero or None.
A `v <= 0` check therefore passes an uncalibrated file straight through and
every downstream length silently becomes "per pixel" wearing micrometre units.
`axesCalibrated` is the signal that distinguishes a real 1 um voxel from an
absent one, and it is checked before any pixel is read. This is not
hypothetical: Dataset_1 of Dryad doi:10.5061/dryad.zcrjdfnph reports
axesCalibrated=(False, False, False) with voxel_size()==(1.0, 1.0, 1.0).
"""

from __future__ import annotations

import numpy as np

from ..acquisition import SEGMENTATION_BASIS
from .field import FieldError
from .ingest import BiomassField, ingest

# How raw intensity becomes a biomass fraction. The mode is DECLARED, never
# inferred from the histogram: where the background ends and the signal begins
# is a segmentation judgement of exactly the same kind as the basis, and an
# argument default must not make it.
TRANSFORMS = ("fixed_threshold", "fixed_background_and_scale",
              "smoke_minmax_per_frame")

# The honest label for "rescale this volume by its own min and max". Adequate
# to prove the array handling works; NOT comparable across volumes, because
# every volume is forced to span [0, 1] whatever its real dynamic range.
SMOKE = ("smoke_minmax_per_frame", float("nan"), float("nan"))

# Axes this reader knows how to reduce. Anything else with a size above one is
# refused rather than silently indexed at zero — 'S' (RGB samples) and 'P'
# (stage positions) both appear in real files and neither is a biomass channel.
_SPATIAL = ("X", "Y", "Z")


class NoTimeAxisError(FieldError):
    """The file has no T axis, so there is no time course to iterate.

    Named so a caller can distinguish "this is a single volume" from every
    other reason a read fails. Catching bare Exception here would turn a bad
    channel, a missing calibration or a corrupt file into "not a time course"
    and retry the same broken read.
    """


def _import_nd2():
    try:
        import nd2
    except ImportError as e:            # pragma: no cover - exercised by skipif
        raise FieldError(
            "reading ND2 needs the optional dataset extra: "
            'pip install -e "calibration[dataset]"') from e
    return nd2


def _voxel_from(f) -> tuple:
    """(x, y, z) in micrometres, or a refusal.

    Checked BEFORE any array is materialised, so an uncalibrated file costs a
    metadata read rather than half a gigabyte.
    """
    vox = f.voxel_size()
    voxel = (float(vox.x), float(vox.y), float(vox.z))

    calibrated = None
    if not f.is_legacy:
        try:
            calibrated = tuple(f.metadata.channels[0].volume.axesCalibrated)
        except (AttributeError, IndexError, TypeError):
            calibrated = None           # absent metadata is not a false claim

    if calibrated is not None and not all(calibrated):
        raise FieldError(
            f"ND2 reports axesCalibrated={calibrated} — this file carries no "
            f"spatial calibration, and nd2 substitutes {voxel} for it. A "
            "voxel size that was never measured must not be used as one.")
    if any(v <= 0 for v in voxel):
        raise FieldError(
            f"ND2 reports a non-positive voxel size {voxel} — the physical "
            "calibration is missing and must not be guessed")
    return voxel, calibrated


def _nd2_meta(f, sizes, voxel, calibrated) -> dict:
    """Structured provenance, for the run's own output to carry."""
    import nd2

    names = []
    try:
        names = [c.channel.name for c in f.metadata.channels]
    except (AttributeError, IndexError, TypeError):
        pass
    return {"nd2_version": nd2.__version__, "sizes": dict(sizes),
            "dtype": str(f.dtype), "is_legacy": bool(f.is_legacy),
            "is_rgb": bool(f.is_rgb), "voxel_size_um": list(voxel),
            "voxel_calibrated": list(calibrated) if calibrated else None,
            "channel_names": names, "n_frames": len(f.loop_indices),
            "nbytes": int(f.nbytes)}


def preflight(path) -> dict:
    """Everything the file declares about itself, reading no pixels.

    Run this BEFORE allocating anything. It is what decides whether a file can
    support the pass you intend, and it answers that for the cost of a header
    read. It deliberately does not raise on an uncalibrated file: reporting
    that a file is uncalibrated is the whole point of a preflight.
    """
    nd2 = _import_nd2()
    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        vox = f.voxel_size()
        voxel = (float(vox.x), float(vox.y), float(vox.z))
        calibrated = None
        interpretation = None
        if not f.is_legacy:
            try:
                vol = f.metadata.channels[0].volume
                calibrated = tuple(vol.axesCalibrated)
                interpretation = tuple(vol.axesInterpretation)
            except (AttributeError, IndexError, TypeError):
                pass
        out = _nd2_meta(f, sizes, voxel, calibrated)
        out.update({
            "path": str(path),
            "axes": list(sizes.keys()),
            "axes_interpretation": list(interpretation) if interpretation else None,
            "n_timepoints": sizes.get("T"),
            "n_z": sizes.get("Z"),
            "n_channels": sizes.get("C"),
            "n_positions": sizes.get("P"),
            "experiment": [(e.type, e.count) for e in f.experiment],
            # The question the whole pass turns on. A biomass field is 3-D by
            # contract, and an axis of size 1 is dropped from f.sizes entirely.
            "is_3d": all(a in sizes for a in _SPATIAL),
            "has_time_axis": "T" in sizes,
        })
        frame_times, time_source = _frame_times(f, sizes.get("T"))
        out["frame_times_s"] = frame_times
        out["time_metadata_source"] = time_source
        return out


def _frame_times(f, n_t):
    """(seconds per timepoint, source) — or a declared absence.

    Never synthesises `t = i * dt`. A file that does not record when its frames
    were taken has not recorded it, and an interval invented here would be
    indistinguishable downstream from a measured one.
    """
    if not n_t:
        return None, "absent"
    try:
        events = f.events()
    except Exception:                   # pragma: no cover - reader-specific
        return None, "absent"
    times = []
    for e in events or ():
        t = e.get("Time [s]") if isinstance(e, dict) else None
        if t is not None:
            times.append(float(t))
    if len(times) >= n_t:
        return times[:n_t], "nd2_events_time_s"
    return None, "absent"


def read_nd2(path, *, segmentation_basis: str, sample_id: str, source_id: str,
             channel: int | None = None, timepoint: int | None = None,
             position: int | None = None,
             transform: tuple = SMOKE, notes: str = "") -> BiomassField:
    """Read one 3-D volume from an ND2 file into a validated biomass field.

    Voxel calibration is taken from the FILE, not from the caller: a hand-typed
    voxel size is the single easiest thing to lose in a conversion, and ND2
    carries it. If the file does not carry it, this raises rather than guessing
    — see the module docstring for why `v > 0` is not enough of a check.

    `transform` is a declared (mode, lo, hi) triple. Which mode is used is a
    segmentation decision, it is recorded in the output, and it is never
    estimated from the data here.
    """
    nd2 = _import_nd2()

    if segmentation_basis not in SEGMENTATION_BASIS:
        raise FieldError(
            f"unknown segmentation basis {segmentation_basis!r}; expected one "
            f"of {sorted(SEGMENTATION_BASIS)}")

    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        voxel, calibrated = _voxel_from(f)      # before any pixel is read
        axes = list(sizes.keys())
        _check_selection(sizes, channel=channel, timepoint=timepoint,
                         position=position)
        arr = f.asarray()
        meta = _nd2_meta(f, sizes, voxel, calibrated)

    vol = _select_volume(arr, axes, channel=channel, timepoint=timepoint,
                         position=position)
    field = _to_fraction(vol, transform)

    meta.update({"channel": channel, "timepoint": timepoint,
                 "position": position, "transform": list(transform)})
    return ingest(field, voxel, segmentation_basis, sample_id=sample_id,
                  source_id=source_id, reader="nd2", meta=meta,
                  notes=(f"sizes={sizes} channel={channel} "
                         f"timepoint={timepoint} transform={transform[0]} "
                         f"{notes}").strip())


def _check_selection(sizes, *, channel, timepoint, position) -> None:
    """Refuse to index an axis the caller did not choose.

    The old reader took index 0 for every axis it did not recognise. That is
    how three fluorescence channels become "the first one", and how an RGB
    brightfield sample axis becomes a biomass channel.
    """
    chosen = {"C": channel, "T": timepoint, "P": position}
    for name, n in sizes.items():
        if name in _SPATIAL:
            continue
        if n <= 1:
            continue                    # only one thing to choose
        if name in chosen:
            if chosen[name] is None:
                flag = {"C": "--channel", "T": "--timepoint",
                        "P": "--position"}[name]
                raise FieldError(
                    f"this ND2 has {n} values on the {name!r} axis and no "
                    f"choice was made; pass {flag} to select one. Which of "
                    "them is the biomass channel is an experimental "
                    "judgement, not a default.")
            if not 0 <= int(chosen[name]) < n:
                raise FieldError(
                    f"{name!r} index {chosen[name]} is out of range for an "
                    f"axis of size {n}")
        else:
            raise FieldError(
                f"this ND2 has {n} values on the {name!r} axis, which this "
                "reader does not know how to choose between. 'S' is an RGB "
                "sample axis and is not a biomass channel; refusing rather "
                "than taking index 0.")


def _select_volume(arr: np.ndarray, axes, *, channel: int | None,
                   timepoint: int | None = None,
                   position: int | None = None) -> np.ndarray:
    """Reduce an ND2 array to one 3-D (x, y, z) volume.

    ND2 axis order is metadata, not convention, so the axes are used by NAME.
    Guessing by position is how a time index becomes a z index.
    """
    a = np.asarray(arr)
    idx = [slice(None)] * a.ndim
    for pos, name in enumerate(axes):
        if name == "C":
            idx[pos] = 0 if channel is None else int(channel)
        elif name == "T":
            if timepoint is None:
                raise FieldError(
                    "this ND2 has a time axis; pass timepoint= to choose one "
                    "volume, or iterate with read_nd2_timecourse()")
            idx[pos] = int(timepoint)
        elif name == "P":
            idx[pos] = 0 if position is None else int(position)
        elif name not in _SPATIAL:
            idx[pos] = 0
    vol = a[tuple(idx)]
    if vol.ndim != 3:
        raise FieldError(
            f"expected a 3-D volume after axis selection, got {vol.ndim}-D "
            f"from axes {axes}. A biomass field is 3-D by contract; a "
            "single-plane acquisition is not a stack.")
    spatial = [n for n in axes if n in _SPATIAL]
    # to logical (x, y, z)
    order = [spatial.index(n) for n in _SPATIAL if n in spatial]
    return np.transpose(vol, order) if len(order) == 3 else vol


def _to_fraction(vol: np.ndarray, transform: tuple) -> np.ndarray:
    """Raw intensity -> biomass fraction in [0, 1], by a DECLARED rule."""
    mode, lo, hi = transform
    if mode not in TRANSFORMS:
        raise FieldError(f"unknown transform {mode!r}; expected one of "
                         f"{list(TRANSFORMS)}")
    v = np.asarray(vol, dtype=float)

    if mode == "fixed_threshold":
        return (v >= float(lo)).astype(float)

    if mode == "fixed_background_and_scale":
        lo, hi = float(lo), float(hi)
        if not hi > lo:
            raise FieldError(
                f"fixed_background_and_scale needs hi > lo, got lo={lo} hi={hi}")
        return np.clip((v - lo) / (hi - lo), 0.0, 1.0)

    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        raise FieldError("volume has no dynamic range; cannot rescale to [0, 1]")
    return (v - lo) / (hi - lo)


def read_nd2_timecourse(path, *, segmentation_basis: str, sample_id: str,
                        source_id: str, channel: int | None = None,
                        position: int | None = None, transform: tuple = SMOKE,
                        max_timepoints: int | None = None) -> list[BiomassField]:
    """Every timepoint as its own validated field, for a temporal observable.

    One transform for every frame, supplied by the caller. Rescaling each frame
    by its own min and max would force them all to span [0, 1] and manufacture
    temporal similarity in a field whose real dynamic range is growing.

    ponytail: opens once and materialises the whole array once, then slices.
    Peak is f.nbytes rather than the old N x f.nbytes, which is enough for the
    0.4-0.6 GB files this has been run against. Switch to loop_indices +
    read_frame if a file arrives that cannot be held whole.
    """
    nd2 = _import_nd2()

    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        n_t = sizes.get("T")
        if not n_t:
            raise NoTimeAxisError(f"{path} has no time axis")
        voxel, calibrated = _voxel_from(f)
        axes = list(sizes.keys())
        _check_selection(sizes, channel=channel, timepoint=0, position=position)
        n = n_t if max_timepoints is None else min(n_t, int(max_timepoints))
        arr = f.asarray()
        meta = _nd2_meta(f, sizes, voxel, calibrated)
        frame_times, time_source = _frame_times(f, n_t)

    out = []
    for t in range(n):
        vol = _select_volume(arr, axes, channel=channel, timepoint=t,
                             position=position)
        m = dict(meta)
        m.update({"channel": channel, "timepoint": t, "position": position,
                  "transform": list(transform),
                  "time_s": frame_times[t] if frame_times else None,
                  "time_metadata_source": time_source})
        out.append(ingest(_to_fraction(vol, transform), voxel,
                          segmentation_basis, sample_id=f"{sample_id}_t{t}",
                          source_id=source_id, reader="nd2", meta=m,
                          notes=(f"sizes={sizes} channel={channel} "
                                 f"timepoint={t} transform={transform[0]}")))
    return out
