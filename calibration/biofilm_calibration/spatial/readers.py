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
"""

from __future__ import annotations

import numpy as np

from ..acquisition import SEGMENTATION_BASIS
from .field import FieldError
from .ingest import BiomassField, ingest


def read_nd2(path, *, segmentation_basis: str, sample_id: str, source_id: str,
             channel: int = 0, timepoint: int | None = None,
             threshold: float | None = None, notes: str = "") -> BiomassField:
    """Read one 3-D volume from an ND2 file into a validated biomass field.

    Voxel calibration is taken from the FILE, not from the caller: a hand-typed
    voxel size is the single easiest thing to lose in a conversion, and ND2
    carries it. If the file does not, this raises rather than guessing.

    `threshold` converts raw intensity to a biomass fraction in [0,1]. Supplying
    it is a segmentation decision and is recorded; omitting it rescales the
    volume linearly to [0,1], which is a probability-like field, not a mask.
    """
    try:
        import nd2
    except ImportError as e:            # pragma: no cover - exercised by skipif
        raise FieldError(
            "reading ND2 needs the optional dataset extra: "
            'pip install -e "calibration[dataset]"') from e

    if segmentation_basis not in SEGMENTATION_BASIS:
        raise FieldError(
            f"unknown segmentation basis {segmentation_basis!r}; expected one "
            f"of {sorted(SEGMENTATION_BASIS)}")

    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        vox = f.voxel_size()            # (x, y, z) in micrometres
        arr = f.asarray()
        axes = list(f.sizes.keys())

    voxel = (float(vox.x), float(vox.y), float(vox.z))
    if any(v <= 0 for v in voxel):
        raise FieldError(
            f"ND2 reports a non-positive voxel size {voxel} — the physical "
            "calibration is missing and must not be guessed")

    vol = _select_volume(arr, axes, channel=channel, timepoint=timepoint)
    field = _to_fraction(vol, threshold)

    return ingest(field, voxel, segmentation_basis, sample_id=sample_id,
                  source_id=source_id, reader="nd2",
                  notes=(f"sizes={sizes} channel={channel} "
                         f"timepoint={timepoint} threshold={threshold} {notes}").strip())


def _select_volume(arr: np.ndarray, axes, *, channel: int,
                   timepoint: int | None) -> np.ndarray:
    """Reduce an ND2 array to one 3-D (x, y, z) volume.

    ND2 axis order is metadata, not convention, so the axes are used by NAME.
    Guessing by position is how a time index becomes a z index.
    """
    a = np.asarray(arr)
    idx = [slice(None)] * a.ndim
    for pos, name in enumerate(axes):
        if name == "C":
            idx[pos] = channel
        elif name == "T":
            if timepoint is None:
                raise FieldError(
                    "this ND2 has a time axis; pass timepoint= to choose one "
                    "volume, or iterate with read_nd2_timecourse()")
            idx[pos] = timepoint
        elif name not in ("X", "Y", "Z"):
            idx[pos] = 0
    vol = a[tuple(idx)]
    if vol.ndim != 3:
        raise FieldError(
            f"expected a 3-D volume after axis selection, got {vol.ndim}-D "
            f"from axes {axes}")
    spatial = [n for n in axes if n in ("X", "Y", "Z")]
    # to logical (x, y, z)
    order = [spatial.index(n) for n in ("X", "Y", "Z") if n in spatial]
    return np.transpose(vol, order) if len(order) == 3 else vol


def _to_fraction(vol: np.ndarray, threshold: float | None) -> np.ndarray:
    """Raw intensity -> biomass fraction in [0, 1]."""
    v = np.asarray(vol, dtype=float)
    if threshold is not None:
        return (v >= float(threshold)).astype(float)
    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        raise FieldError("volume has no dynamic range; cannot rescale to [0, 1]")
    return (v - lo) / (hi - lo)


def read_nd2_timecourse(path, *, segmentation_basis: str, sample_id: str,
                        source_id: str, channel: int = 0,
                        threshold: float | None = None) -> list[BiomassField]:
    """Every timepoint as its own validated field, for the temporal observable."""
    try:
        import nd2
    except ImportError as e:            # pragma: no cover
        raise FieldError("reading ND2 needs the optional dataset extra") from e

    with nd2.ND2File(str(path)) as f:
        n_t = dict(f.sizes).get("T")
    if not n_t:
        raise FieldError(f"{path} has no time axis")
    return [read_nd2(path, segmentation_basis=segmentation_basis,
                     sample_id=f"{sample_id}_t{t}", source_id=source_id,
                     channel=channel, timepoint=t, threshold=threshold)
            for t in range(n_t)]
