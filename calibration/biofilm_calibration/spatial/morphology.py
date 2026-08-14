"""Morphology metrics for a 3-D labelled volume.

The same function must run over a segmented microscopy stack and over a CPM
snapshot's `cell_id` array, or the comparison between them is not a comparison.
That is the whole design constraint here: metrics take a label array and a
voxel size, and know nothing about where the labels came from.

Axis lengths come from the inertia tensor rather than a bounding box, because a
bounding box is orientation-dependent and a rod lying diagonally would report
the wrong aspect ratio — precisely the shape question this harness exists to
answer.

Pure numpy: no scipy, no scikit-image. Connected components are needed only for
a 26-connectivity check on objects that are already small, so a compact
flood-fill is cheaper than a dependency.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ObjectMorphology:
    label: int
    n_voxels: int
    volume_um3: float
    equivalent_diameter_um: float
    major_axis_um: float
    intermediate_axis_um: float
    minor_axis_um: float
    aspect_ratio_major_minor: float
    connected_components: int
    boundary_contact: bool
    centroid_um: tuple

    @property
    def is_isotropic(self) -> bool:
        """A blob by the usual convention: no axis more than 1.5x another."""
        return self.aspect_ratio_major_minor <= 1.5


def _principal_axis_lengths(coords: np.ndarray, voxel_um: np.ndarray) -> tuple:
    """Equivalent-ellipsoid axis LENGTHS (full, not semi-) in micrometres.

    From the eigenvalues of the covariance of the voxel centres. For a solid
    ellipsoid with semi-axes a>=b>=c the covariance eigenvalues are a^2/5 etc.,
    so length = 2*sqrt(5*lambda). A single voxel gives zero variance and hence
    zero length, which is correct: one voxel has no resolvable shape.
    """
    if coords.shape[0] < 2:
        return (0.0, 0.0, 0.0)
    pts = coords * voxel_um
    centred = pts - pts.mean(axis=0)
    cov = (centred.T @ centred) / coords.shape[0]
    eig = np.linalg.eigvalsh(cov)                 # ascending
    lengths = 2.0 * np.sqrt(5.0 * np.clip(eig, 0.0, None))
    return tuple(float(x) for x in lengths[::-1])  # major, intermediate, minor


def _count_components_26(mask: np.ndarray) -> int:
    """Number of 26-connected components in a small boolean sub-volume."""
    remaining = mask.copy()
    n = 0
    idx = np.argwhere(remaining)
    offsets = [(dx, dy, dz)
               for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
               if (dx, dy, dz) != (0, 0, 0)]
    shape = mask.shape
    for start in map(tuple, idx):
        if not remaining[start]:
            continue
        n += 1
        stack = [start]
        remaining[start] = False
        while stack:
            x, y, z = stack.pop()
            for dx, dy, dz in offsets:
                p = (x + dx, y + dy, z + dz)
                if all(0 <= p[i] < shape[i] for i in range(3)) and remaining[p]:
                    remaining[p] = False
                    stack.append(p)
    return n


def object_morphology(labels: np.ndarray, label: int,
                      voxel_size_um) -> ObjectMorphology:
    """Morphology of one label in a 3-D label array."""
    labels = np.asarray(labels)
    if labels.ndim != 3:
        raise ValueError(f"expected a 3-D label array, got {labels.ndim}-D")
    voxel = np.asarray(voxel_size_um, dtype=float)
    if voxel.shape == ():
        voxel = np.repeat(voxel, 3)
    if voxel.shape != (3,) or np.any(voxel <= 0):
        raise ValueError(f"voxel_size_um must be 3 positive values, got {voxel_size_um!r}")

    coords = np.argwhere(labels == label)
    if coords.size == 0:
        raise ValueError(f"label {label} is not present")

    n = int(coords.shape[0])
    voxel_volume = float(np.prod(voxel))
    volume = n * voxel_volume
    equiv_d = float((6.0 * volume / np.pi) ** (1.0 / 3.0))
    major, inter, minor = _principal_axis_lengths(coords, voxel)
    aspect = float(major / minor) if minor > 0 else float("inf")

    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    touching = bool(np.any(lo == 0) or np.any(hi == np.array(labels.shape) - 1))

    sub = np.zeros(hi - lo + 1, dtype=bool)
    sub[tuple((coords - lo).T)] = True

    return ObjectMorphology(
        label=int(label), n_voxels=n, volume_um3=volume,
        equivalent_diameter_um=equiv_d, major_axis_um=major,
        intermediate_axis_um=inter, minor_axis_um=minor,
        aspect_ratio_major_minor=aspect,
        connected_components=_count_components_26(sub),
        boundary_contact=touching,
        centroid_um=tuple(float(v) for v in (coords.mean(axis=0) * voxel)),
    )


def all_objects(labels: np.ndarray, voxel_size_um,
                background=(0, -1)) -> list[ObjectMorphology]:
    """Morphology of every label except the background ones.

    `-1` is excluded by default because the CPM uses it for "outside the
    biological domain" — it is a wall marker, not an object.
    """
    labels = np.asarray(labels)
    present = [int(v) for v in np.unique(labels) if int(v) not in set(background)]
    return [object_morphology(labels, v, voxel_size_um) for v in present]


def occupancy(labels: np.ndarray, background=(0, -1)) -> float:
    """Fraction of non-wall sites that are occupied.

    DIMENSIONLESS BY CONSTRUCTION, and named so. This is the same quantity
    `compute_radial_biomass` produces, and it is *not* a biomass density — see
    docs/calibration/material_basis_contract.md.
    """
    labels = np.asarray(labels)
    wall = np.isin(labels, [b for b in background if b < 0])
    considered = int((~wall).sum())
    if considered == 0:
        return float("nan")
    occupied = int((~wall & ~np.isin(labels, list(background))).sum())
    return occupied / considered


def size_distribution_summary(objects, quantiles=(0.05, 0.25, 0.5, 0.75, 0.95)) -> dict:
    """Quantiles of the volume distribution, plus the counts a reader needs to
    judge them. Distributions, not means: two populations with the same mean
    volume can be completely different biology."""
    vols = np.array([o.volume_um3 for o in objects if not o.boundary_contact])
    excluded = sum(1 for o in objects if o.boundary_contact)
    if vols.size == 0:
        return {"n": 0, "n_excluded_boundary": excluded}
    out = {"n": int(vols.size), "n_excluded_boundary": excluded,
           "mean_um3": float(vols.mean())}
    for q in quantiles:
        out[f"q{int(q * 100):02d}_um3"] = float(np.quantile(vols, q))
    return out
