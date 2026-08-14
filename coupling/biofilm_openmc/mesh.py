"""Tally-mesh resolution, decoupled from physical geometry.

The tally mesh used to be welded to the CPM lattice: dimension `(n,n,n)` over
an extent of `n * voxel_pitch_cm`. Changing the pitch to sweep resolution also
changed the size of the modeled apparatus, so a convergence study could not
separate discretization error from a change in geometry. Here the physical
domain is fixed and only `coarsening_factor` moves.

Meshes are nested by construction: the factor must divide every base dimension
exactly, so each coarse bin is a whole block of fine bins and a fine result can
be aggregated onto a coarser grid for comparison.

WHY THE AGGREGATION IS A SUM. The quantities coarsened here are EXTENSIVE —
deposited energy and mass. Coarsening each separately and dividing afterwards
is what makes a coarse dose equal the mass-weighted mean of the fine dose it
replaces:

    D_coarse = sum_v E_v / sum_v m_v  ==  weighted_mean(E_v / m_v, weights=m_v)

Averaging dose directly would silently weight every bin equally, and taking the
bin volume from `voxel_pitch_cm` instead of the mesh extent would understate
dose by exactly the coarsening factor cubed. `tests/test_mesh.py` pins both.

Pure numpy: no `openmc` import, so this is testable in the bare tier.
"""

from __future__ import annotations

import numpy as np

from .config import ConfigError


def resolve_mesh_dimension(base_dimension, coarsening_factor: int) -> tuple[int, int, int]:
    """Tally dimension after coarsening. Refuses a factor that does not divide
    the base evenly — nesting is what makes cross-resolution comparison valid."""
    base = tuple(int(d) for d in base_dimension)
    f = int(coarsening_factor)
    if len(base) != 3:
        raise ConfigError(f"mesh base dimension must be 3-D, got {list(base)}")
    if f < 1:
        raise ConfigError(f"mesh coarsening_factor must be >= 1, got {f}")
    bad = [d for d in base if d % f]
    if bad:
        raise ConfigError(
            f"mesh coarsening_factor {f} does not divide base dimension "
            f"{list(base)} evenly — coarse bins would straddle fine ones and "
            "results at different resolutions could not be compared")
    return tuple(d // f for d in base)


def coarsen_field(arr: np.ndarray, factor: int) -> np.ndarray:
    """Block-SUM an extensive field (energy, mass) by an integer factor.

    Not a mean: see the module docstring. Use on energy and mass separately,
    then divide, to obtain the coarse intensive quantity.
    """
    f = int(factor)
    if f == 1:
        return np.asarray(arr)
    a = np.asarray(arr)
    if a.ndim != 3:
        raise ConfigError(f"coarsen_field expects a 3-D field, got {a.ndim}-D")
    resolve_mesh_dimension(a.shape, f)          # reuse the divisibility refusal
    nx, ny, nz = (d // f for d in a.shape)
    return a.reshape(nx, f, ny, f, nz, f).sum(axis=(1, 3, 5))


def upsample_field(arr: np.ndarray, factor: int) -> np.ndarray:
    """Broadcast a coarse INTENSIVE field back onto the fine grid.

    Valid for dose RATE precisely because it is intensive: every fine voxel
    inside a coarse bin carries that bin's average rate. This is how lineage
    attribution keeps its meaning under coarsening — labels stay at lattice
    resolution while the radiation estimate is coarser. Never use it on an
    extensive quantity; energy and mass would be multiplied by factor cubed.
    """
    f = int(factor)
    if f == 1:
        return np.asarray(arr)
    a = np.asarray(arr)
    for axis in range(3):
        a = np.repeat(a, f, axis=axis)
    return a


def biofilm_mesh_extent_cm(config, n: int) -> tuple[float, float, float]:
    """The biofilm mesh spans the LATTICE cube, side `n * voxel_pitch_cm` — not
    the cylinder, which `check_lattice_congruence` only requires to fit inside
    it. Keeping this as the lattice cube preserves exactly what the foundation
    tallied; only the number of bins across it changes."""
    side = n * config.voxel_pitch_cm
    return (side, side, side)


def phantom_mesh_extent_cm(config) -> tuple[float, float, float]:
    """The phantom mesh spans the cube circumscribing the water cylinder."""
    d = 2.0 * config.cylinder_radius_cm
    return (d, d, float(config.cylinder_length_cm))


def mesh_bin_volume_cm3(extent_cm, dimension) -> float:
    """Volume of one tally bin, from the PHYSICAL extent and the bin count.

    The dose denominator must come from here, never from `voxel_pitch_cm`:
    that is the lattice pitch, and using it once the mesh is coarsened
    understates every dose by the coarsening factor cubed.
    """
    ex, ey, ez = (float(e) for e in extent_cm)
    nx, ny, nz = (int(d) for d in dimension)
    return (ex / nx) * (ey / ny) * (ez / nz)
