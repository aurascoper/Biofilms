"""Structural observables for a biomass field.

These are the declared calibration observables: the quantities a coarse-grained
lattice must preserve for its pitch to be acceptable. They are field-level
rather than object-level because a CPM ID is a biomass parcel, so what has to
survive coarsening is the STRUCTURE of the biomass distribution, not the shape
of any individual organism.

Pure numpy. Two-point correlation uses an FFT so it stays cheap on a full
stack; everything else is direct.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .field import check_field

# Axis convention: z is the depth axis (normal to the substrate), matching the
# CPM's cylinder axis and the usual confocal stack orientation.
DEPTH_AXIS = 2


@dataclass(frozen=True)
class StructureSummary:
    biovolume_fraction: float
    porosity: float
    specific_interface_area_per_um: float
    correlation_length_um: float
    thickness_mean_um: float
    thickness_sd_um: float
    n_components: int
    component_size_q50_um3: float
    component_size_q95_um3: float

    def as_dict(self) -> dict:
        return asdict(self)


def _voxel(voxel_size_um) -> np.ndarray:
    v = np.asarray(voxel_size_um, dtype=float)
    if v.shape == ():
        v = np.repeat(v, 3)
    if v.shape != (3,) or np.any(v <= 0):
        raise ValueError(f"voxel_size_um must be 3 positive values, got {voxel_size_um!r}")
    return v


def porosity(B: np.ndarray) -> float:
    """Void fraction. For a fractional field this is 1 - <B>, i.e. the
    complement of the biovolume fraction, not a pore-network measure."""
    return 1.0 - float(check_field(B).mean())


def specific_interface_area(mask: np.ndarray, voxel_size_um) -> float:
    """Biomass-void interface area per unit total volume, in 1/um.

    Counts INTERNAL transitions only — faces where an occupied voxel meets an
    unoccupied one inside the array. The array's outer boundary is deliberately
    excluded: biomass reaching the edge of the field of view is cut off by the
    IMAGE, not by an interface with void, and counting it would make the metric
    depend on how large a crop was taken rather than on the structure. Same
    rule as the chord-length truncation below.

    It remains a discretisation-dependent estimate, which is precisely why it
    is a calibration observable: a coarser lattice smooths the interface away,
    and the rate at which it does so is the thing being measured.
    """
    m = np.asarray(mask).astype(bool)
    if m.ndim != 3:
        raise ValueError(f"expected a 3-D mask, got {m.ndim}-D")
    v = _voxel(voxel_size_um)
    face_area = (v[1] * v[2], v[0] * v[2], v[0] * v[1])

    area = 0.0
    for axis in range(3):
        a = np.swapaxes(m, 0, axis)
        internal = int(np.logical_xor(a[:-1], a[1:]).sum())
        area += internal * face_area[axis]
    total_volume = float(np.prod(np.array(m.shape) * v))
    return area / total_volume if total_volume else float("nan")


def two_point_correlation(B: np.ndarray, voxel_size_um,
                          axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Autocorrelation of the biomass fluctuation along one axis.

    Returns (lags_um, correlation) with correlation normalised to 1 at zero
    lag. Computed by FFT on the mean-subtracted field, so it is the covariance
    function rather than the raw two-point probability.
    """
    a = check_field(B)
    v = _voxel(voxel_size_um)
    x = a - a.mean()
    n = x.shape[axis]
    # correlate along `axis`, averaging over the other two
    xf = np.fft.rfft(x, n=2 * n, axis=axis)
    acf = np.fft.irfft(xf * np.conjugate(xf), n=2 * n, axis=axis)
    acf = np.moveaxis(acf, axis, 0)[:n]
    acf = acf.reshape(n, -1).mean(axis=1)
    if acf[0] <= 0:
        return np.zeros(n), np.zeros(n)
    # unbiased for the shrinking overlap
    acf = acf / np.arange(n, 0, -1)
    acf = acf / acf[0]
    lags = np.arange(n) * v[axis]
    return lags, acf


def correlation_length(B: np.ndarray, voxel_size_um, axis: int = 0) -> float:
    """Lag at which the autocorrelation first falls to 1/e, interpolated.

    A single number standing in for the structure's characteristic scale. If a
    candidate pitch is comparable to this, the structure is being erased rather
    than resolved.
    """
    lags, acf = two_point_correlation(B, voxel_size_um, axis)
    if acf.size < 2 or not np.isfinite(acf).all():
        return float("nan")
    target = 1.0 / np.e
    below = np.nonzero(acf < target)[0]
    if below.size == 0:
        return float("nan")            # never decorrelates within the field
    i = int(below[0])
    if i == 0:
        return 0.0
    y0, y1 = acf[i - 1], acf[i]
    frac = (y0 - target) / (y0 - y1) if y0 != y1 else 0.0
    return float(lags[i - 1] + frac * (lags[i] - lags[i - 1]))


def thickness_profile(mask: np.ndarray, voxel_size_um) -> np.ndarray:
    """Biofilm thickness at each (x, y) column, in um.

    Occupied-voxel count along the depth axis rather than top-minus-bottom, so
    an internal void reduces the thickness instead of being counted as biomass.
    """
    m = np.asarray(mask).astype(bool)
    v = _voxel(voxel_size_um)
    return m.sum(axis=DEPTH_AXIS).astype(float) * v[DEPTH_AXIS]


def depth_profile(B: np.ndarray) -> np.ndarray:
    """Mean biomass fraction at each depth slice."""
    a = check_field(B)
    return a.mean(axis=(0, 1))


def components_26(mask: np.ndarray) -> list[int]:
    """Sizes, in voxels, of every 26-connected component.

    Sizes rather than a count: a field that keeps its total biovolume while
    fragmenting into many small pieces has not been preserved, and only the
    size distribution shows that.
    """
    m = np.asarray(mask).astype(bool)
    if m.ndim != 3:
        raise ValueError(f"expected a 3-D mask, got {m.ndim}-D")
    remaining = m.copy()
    offsets = [(dx, dy, dz)
               for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
               if (dx, dy, dz) != (0, 0, 0)]
    shape = m.shape
    sizes: list[int] = []
    for start in map(tuple, np.argwhere(remaining)):
        if not remaining[start]:
            continue
        size = 0
        stack = [start]
        remaining[start] = False
        while stack:
            x, y, z = stack.pop()
            size += 1
            for dx, dy, dz in offsets:
                p = (x + dx, y + dy, z + dz)
                if all(0 <= p[i] < shape[i] for i in range(3)) and remaining[p]:
                    remaining[p] = False
                    stack.append(p)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def chord_lengths(mask: np.ndarray, voxel_size_um, axis: int = 0,
                  phase: bool = True) -> np.ndarray:
    """Lengths of uninterrupted runs through one phase, along one axis, in um.

    Runs touching the array edge are dropped: they are truncated by the field
    of view, not by the structure, and would bias the distribution short.
    """
    m = np.asarray(mask).astype(bool)
    if not phase:
        m = ~m
    v = _voxel(voxel_size_um)
    a = np.moveaxis(m, axis, -1)
    out: list[int] = []
    for line in a.reshape(-1, a.shape[-1]):
        run = 0
        touched_start = False
        for i, occupied in enumerate(line):
            if occupied:
                if run == 0 and i == 0:
                    touched_start = True
                run += 1
            elif run:
                if not touched_start:
                    out.append(run)
                run, touched_start = 0, False
        # a run reaching the end is truncated by the field of view
    return np.asarray(out, dtype=float) * v[axis]


def summarise(B: np.ndarray, mask: np.ndarray, voxel_size_um) -> StructureSummary:
    """Every declared observable, from a field and its occupancy mask."""
    v = _voxel(voxel_size_um)
    voxel_volume = float(np.prod(v))
    sizes = components_26(mask)
    thick = thickness_profile(mask, v)
    sizes_um3 = np.asarray(sizes, dtype=float) * voxel_volume
    return StructureSummary(
        biovolume_fraction=float(check_field(B).mean()),
        porosity=porosity(B),
        specific_interface_area_per_um=specific_interface_area(mask, v),
        correlation_length_um=correlation_length(B, v),
        thickness_mean_um=float(thick.mean()),
        thickness_sd_um=float(thick.std()),
        n_components=len(sizes),
        component_size_q50_um3=float(np.quantile(sizes_um3, 0.5)) if sizes else 0.0,
        component_size_q95_um3=float(np.quantile(sizes_um3, 0.95)) if sizes else 0.0,
    )
