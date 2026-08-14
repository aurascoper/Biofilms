"""Synthetic shape classes, and what the CPM can hold.

Two different questions live here and must not be confused.

RESOLVABLE — can a shape be *written down* at this pitch? A rod 1 um across is
not resolvable at a 2 um pitch, whatever the model does. This is arithmetic and
`scale_candidates` answers it.

MAINTAINABLE — will the model *keep* the shape once it has it? The CPM's energy
is adhesion + volume + radiation + melanin (`compute_delta_H`, L482-543). There
is no surface-area, elongation, aspect-ratio or connectivity term, so adhesion
and the volume constraint together minimise surface at fixed volume and every
object relaxes toward an isotropic blob. Nothing can hold a rod elongated or a
hypha extended. No pitch changes this.

So the synthetic masks below are not a test of the CPM. They validate that the
comparison harness can *detect* the blob collapse when it happens — that
`morphology.py` distinguishes a rod from a sphere of equal volume. A harness
that reported them as similar would hide the model's central limitation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .morphology import ObjectMorphology, object_morphology


@dataclass(frozen=True)
class ShapeClass:
    name: str
    maintainable_by_cpm: bool
    reason: str


# What the current Hamiltonian can hold, independent of pitch.
SHAPE_CLASSES = (
    ShapeClass("sphere", True,
               "isotropic; adhesion + volume already minimise surface at fixed "
               "volume, so this is the CPM's attractor"),
    ShapeClass("biomass_parcel", True,
               "a compact parcel is the same attractor, and is what the entity "
               "semantics declare a cell ID to be"),
    ShapeClass("ellipsoid", False,
               "mild anisotropy decays: no term rewards an elongated state"),
    ShapeClass("rod", False,
               "needs a surface-area or elongation term; B. subtilis is "
               "described as a motile rod but cannot be held as one"),
    ShapeClass("filament", False,
               "needs elongation AND connectivity terms; the filamentous fungi "
               "C. sphaerospermum and A. niger cannot be represented"),
    ShapeClass("disconnected", False,
               "no connectivity constraint exists in the energy; a split object "
               "is neither prevented nor penalised"),
)


def shape_class(name: str) -> ShapeClass:
    for s in SHAPE_CLASSES:
        if s.name == name:
            return s
    raise KeyError(f"unknown shape class {name!r}")


def _blank(shape) -> np.ndarray:
    return np.zeros(tuple(shape), dtype=np.int32)


def make_sphere(shape=(24, 24, 24), radius=6.0, label=1) -> np.ndarray:
    arr = _blank(shape)
    c = (np.array(shape) - 1) / 2.0
    idx = np.indices(shape).astype(float)
    r2 = sum((idx[i] - c[i]) ** 2 for i in range(3))
    arr[r2 <= radius ** 2] = label
    return arr


def make_ellipsoid(shape=(24, 24, 24), semi=(10.0, 5.0, 5.0), label=1) -> np.ndarray:
    arr = _blank(shape)
    c = (np.array(shape) - 1) / 2.0
    idx = np.indices(shape).astype(float)
    q = sum(((idx[i] - c[i]) / semi[i]) ** 2 for i in range(3))
    arr[q <= 1.0] = label
    return arr


def make_rod(shape=(32, 16, 16), length=22.0, radius=3.0, label=1) -> np.ndarray:
    """A capsule along x: cylinder plus hemispherical caps."""
    arr = _blank(shape)
    c = (np.array(shape) - 1) / 2.0
    idx = np.indices(shape).astype(float)
    dx = idx[0] - c[0]
    r2 = (idx[1] - c[1]) ** 2 + (idx[2] - c[2]) ** 2
    half = length / 2.0 - radius
    body = (np.abs(dx) <= half) & (r2 <= radius ** 2)
    caps = ((np.abs(dx) - half) ** 2 + r2 <= radius ** 2) & (np.abs(dx) > half)
    arr[body | caps] = label
    return arr


def make_filament(shape=(40, 24, 24), radius=2.0, label=1) -> np.ndarray:
    """A branched filament: a main hypha with one lateral branch."""
    arr = _blank(shape)
    idx = np.indices(shape).astype(float)
    cy, cz = (shape[1] - 1) / 2.0, (shape[2] - 1) / 2.0
    main = ((idx[1] - cy) ** 2 + (idx[2] - cz) ** 2 <= radius ** 2) & \
           (idx[0] >= 3) & (idx[0] <= shape[0] - 4)
    bx = shape[0] * 0.6
    branch = ((idx[0] - bx) ** 2 + (idx[2] - cz) ** 2 <= radius ** 2) & \
             (idx[1] >= cy) & (idx[1] <= shape[1] - 4)
    arr[main | branch] = label
    return arr


def make_disconnected(shape=(24, 24, 24), radius=4.0, gap=6, label=1) -> np.ndarray:
    """Two spheres sharing one label and not touching."""
    arr = _blank(shape)
    idx = np.indices(shape).astype(float)
    c = (np.array(shape) - 1) / 2.0
    for sign in (-1, 1):
        cx = c[0] + sign * (radius + gap / 2.0)
        r2 = (idx[0] - cx) ** 2 + (idx[1] - c[1]) ** 2 + (idx[2] - c[2]) ** 2
        arr[r2 <= radius ** 2] = label
    return arr


def make_biomass_parcel(shape=(24, 24, 24), radius=6.0, seed=0,
                        label=1) -> np.ndarray:
    """A compact but irregular blob — what a relaxed CPM object looks like."""
    arr = make_sphere(shape, radius, label)
    rng = np.random.default_rng(seed)
    surface = np.argwhere(arr == label)
    keep = rng.random(surface.shape[0]) > 0.12
    arr[:] = 0
    arr[tuple(surface[keep].T)] = label
    return arr


BUILDERS = {
    "sphere": make_sphere,
    "ellipsoid": make_ellipsoid,
    "rod": make_rod,
    "filament": make_filament,
    "disconnected": make_disconnected,
    "biomass_parcel": make_biomass_parcel,
}


def audit(voxel_size_um=1.0) -> list[dict]:
    """Measure every synthetic class and pair it with what the CPM can hold.

    The `detectable` column is the one that matters for the harness: it says
    whether the metrics separate this class from an isotropic blob, i.e.
    whether a collapse to blobs would be *visible* in a comparison.
    """
    rows = []
    for name, build in BUILDERS.items():
        m: ObjectMorphology = object_morphology(build(), 1, voxel_size_um)
        cls = shape_class(name)
        rows.append({
            "shape_class": name,
            "n_voxels": m.n_voxels,
            "volume_um3": m.volume_um3,
            "major_axis_um": m.major_axis_um,
            "minor_axis_um": m.minor_axis_um,
            "aspect_ratio": m.aspect_ratio_major_minor,
            "connected_components": m.connected_components,
            "detectable_as_non_blob": (not m.is_isotropic)
                                      or m.connected_components > 1,
            "maintainable_by_cpm": cls.maintainable_by_cpm,
            "reason": cls.reason,
        })
    return rows
