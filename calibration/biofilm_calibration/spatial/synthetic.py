"""Synthetic biomass fields with analytically known structure.

These validate the metrics, not a download. Each builder has a property that
can be computed by hand, so a failing test says the metric is wrong rather than
that the data changed:

    slab                 exact biovolume fraction, exact interface area
    spheres_on_substrate exact component count, near-exact biovolume
    correlated_field     known correlation length, controllable biovolume
    graded_field         monotone depth profile with a known mean

Real biofilm data can be run through the same path later as a data-only
commit; nothing here needs to change for that.
"""

from __future__ import annotations

import numpy as np


def slab(shape=(16, 16, 32), thickness_voxels: int = 8,
         value: float = 1.0) -> np.ndarray:
    """A uniform slab on the substrate, occupying the first `thickness_voxels`
    of the depth axis.

    Biovolume fraction is exactly `value * thickness / shape[2]`, and the only
    internal interface is the one face at the top of the slab.
    """
    B = np.zeros(shape, dtype=float)
    B[:, :, :thickness_voxels] = float(value)
    return B


def spheres_on_substrate(shape=(40, 40, 20), radius: float = 4.0,
                         spacing: int = 12, z_centre: int = 6) -> np.ndarray:
    """A regular grid of non-touching spheres. The component count is exactly
    the number of grid positions."""
    B = np.zeros(shape, dtype=float)
    idx = np.indices(shape).astype(float)
    centres = [(x, y) for x in range(spacing // 2, shape[0], spacing)
               for y in range(spacing // 2, shape[1], spacing)]
    for cx, cy in centres:
        r2 = ((idx[0] - cx) ** 2 + (idx[1] - cy) ** 2
              + (idx[2] - z_centre) ** 2)
        B[r2 <= radius ** 2] = 1.0
    return B


def n_spheres(shape=(40, 40, 20), spacing: int = 12) -> int:
    """How many spheres `spheres_on_substrate` places — the expected component
    count, computed independently of the builder's loop."""
    return len(range(spacing // 2, shape[0], spacing)) * \
        len(range(spacing // 2, shape[1], spacing))


def correlated_field(shape=(64, 64, 64), correlation_length_voxels: float = 8.0,
                     seed: int = 0, target_fraction: float = 0.3) -> np.ndarray:
    """A Gaussian random field smoothed to a known correlation length, then
    mapped to [0,1] with a chosen mean biomass fraction.

    Smoothing is a separable Gaussian convolution by FFT, so the imposed
    correlation length is a property of the construction rather than of the
    measurement.
    """
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(shape)

    sigma = float(correlation_length_voxels) / 2.0
    smoothed = noise
    for axis in range(3):
        n = shape[axis]
        freq = np.fft.fftfreq(n)
        kernel = np.exp(-2.0 * (np.pi * sigma * freq) ** 2)
        shape_k = [1, 1, 1]
        shape_k[axis] = n
        smoothed = np.real(np.fft.ifft(
            np.fft.fft(smoothed, axis=axis) * kernel.reshape(shape_k), axis=axis))

    # map to [0,1] with the requested mean via a quantile shift
    flat = smoothed.ravel()
    cut = np.quantile(flat, 1.0 - float(target_fraction))
    scale = flat.std() or 1.0
    B = 1.0 / (1.0 + np.exp(-(smoothed - cut) / (0.25 * scale)))
    return np.clip(B, 0.0, 1.0)


def graded_field(shape=(16, 16, 32), top: float = 0.9,
                 bottom: float = 0.1) -> np.ndarray:
    """Biomass decreasing linearly with depth. The depth profile is exact and
    the mean is (top + bottom) / 2."""
    B = np.zeros(shape, dtype=float)
    ramp = np.linspace(top, bottom, shape[2])
    B[:, :] = ramp
    return B
