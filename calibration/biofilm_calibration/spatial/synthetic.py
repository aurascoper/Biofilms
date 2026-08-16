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


# --------------------------------------------------------------------------
# PHYSICAL SPECIFICATIONS: one object, many pitches.
#
# Everything above takes feature sizes in VOXELS, which is right for validating
# a metric against a hand-computed property at one resolution and useless for
# asking what a pitch costs. `slab(shape=(16,16,32), thickness_voxels=8)` at a
# different shape is a different slab, so a "resolution sweep" built from those
# builders sweeps the object and the sampling together and cannot separate them.
#
# These carry the object in MICROMETRES and rasterise it at whatever pitch is
# asked for. The physical object is fixed; only the sampling moves. That is the
# whole point, and it is the same discipline `mesh.py` applies on the transport
# side, where the tally extent is fixed by geometry and only the bin count
# changes.
#
# THIS IS NOT A CPM CONVERGENCE STUDY, and must not be mistaken for one. There
# are no dynamics here. `V_target` in biofilms_potts.jl is in SITES, so refining
# a CPM lattice at fixed V_target shrinks every parcel — a different physical
# system, not a finer view of the same one. Rasterising a fixed analytic shape
# has no such confound, which is exactly why it is the first study to run.

from dataclasses import dataclass


def _grid_shape(extent_um, pitch_um) -> tuple:
    return tuple(max(1, int(round(float(e) / float(pitch_um)))) for e in extent_um)


@dataclass(frozen=True)
class PhysicalSlab:
    """A slab of fixed physical thickness in a fixed physical box.

    Ground truth is exact and pitch-independent: the biovolume fraction is
    `thickness_um / extent_um[2]`, and there is exactly one internal interface
    plane, so the specific interface area is `1 / extent_um[2]`.
    """
    extent_um: tuple = (32.0, 32.0, 64.0)
    thickness_um: float = 16.0

    def rasterize(self, pitch_um: float) -> np.ndarray:
        shape = _grid_shape(self.extent_um, pitch_um)
        depth = np.arange(shape[2], dtype=float) * float(pitch_um)
        # A voxel is occupied when its CENTRE lies inside the slab. Centre
        # sampling rather than "any overlap" because the latter systematically
        # over-fills by half a voxel at every boundary, which would show up as
        # a pitch-dependent biovolume bias that is an artifact of the sampling
        # rule rather than of the resolution.
        occupied = (depth + 0.5 * float(pitch_um)) < float(self.thickness_um)
        B = np.zeros(shape, dtype=float)
        B[:, :, occupied] = 1.0
        return B

    @property
    def true_biovolume_fraction(self) -> float:
        return float(self.thickness_um) / float(self.extent_um[2])

    @property
    def true_specific_interface_area_per_um(self) -> float:
        """One internal face of area X*Y in a box of volume X*Y*Z."""
        return 1.0 / float(self.extent_um[2])

    @property
    def true_thickness_mean_um(self) -> float:
        return float(self.thickness_um)


@dataclass(frozen=True)
class PhysicalSpheres:
    """A regular lattice of non-touching spheres at fixed physical size."""
    extent_um: tuple = (48.0, 48.0, 24.0)
    radius_um: float = 5.0
    spacing_um: float = 16.0
    z_centre_um: float = 12.0

    def _centres(self) -> list:
        step = float(self.spacing_um)
        xs = np.arange(step / 2, float(self.extent_um[0]), step)
        ys = np.arange(step / 2, float(self.extent_um[1]), step)
        return [(x, y) for x in xs for y in ys]

    def rasterize(self, pitch_um: float) -> np.ndarray:
        shape = _grid_shape(self.extent_um, pitch_um)
        p = float(pitch_um)
        idx = np.indices(shape).astype(float)
        for axis in range(3):
            idx[axis] = (idx[axis] + 0.5) * p          # voxel centres, um
        B = np.zeros(shape, dtype=float)
        for cx, cy in self._centres():
            r2 = ((idx[0] - cx) ** 2 + (idx[1] - cy) ** 2
                  + (idx[2] - float(self.z_centre_um)) ** 2)
            B[r2 <= float(self.radius_um) ** 2] = 1.0
        return B

    @property
    def n_spheres(self) -> int:
        return len(self._centres())

    @property
    def true_biovolume_fraction(self) -> float:
        volume = self.n_spheres * (4.0 / 3.0) * np.pi * float(self.radius_um) ** 3
        return float(volume / np.prod([float(e) for e in self.extent_um]))

    @property
    def true_component_size_um3(self) -> float:
        return float((4.0 / 3.0) * np.pi * float(self.radius_um) ** 3)

    @property
    def true_specific_interface_area_per_um(self) -> float:
        area = self.n_spheres * 4.0 * np.pi * float(self.radius_um) ** 2
        return float(area / np.prod([float(e) for e in self.extent_um]))
