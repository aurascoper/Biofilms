"""Do the dose field and the CPM labels describe the same place?

WHAT IS ALREADY COVERED, so this file does not repeat it. Axis permutation is
tested in both directions: `test_lattice_orientation_against_probes` checks
four asymmetric probe points through `lat.universes[z][n-1-y][x]`,
`test_corrupted_axis_order_is_refused` transposes two file axes and requires
the snapshot probes to reject it, and `test_dose_written_by_python_reads_back_in_julia`
round-trips a field encoding `x + 10y + 100z` through real Julia. Those are
index-space checks and they are good ones.

WHAT IS NOT COVERED IS WORLD-SPACE REGISTRATION BETWEEN TWO GRIDS, and the
reason is structural: there are two different transforms wearing one name in
the comments.

    RectLattice.universes   transpose(2,1,0)[:, ::-1, :]   HAS a y-flip
    RegularMesh tally       reshape(dim[::-1]).transpose(2,1,0)   NO y-flip

The first lives in `snapshot.to_openmc_lattice_order` and is probe-verified.
The second is reimplemented independently in `dose.extract_heating` and again
in `materials.mesh_material_masses_kg`, and neither is probed. Nothing asserts
that the voxel which is hot in the tally is the voxel which holds biomass in
the lattice. `synthetic_e2e.py`'s invariant 7 comes closest and cannot close
it: a >3sigma near/far grouping about a near-central source still looks
monotonic under a flip.

WHY NOT A CENTROID. The obvious check -- dose-weighted world centroid against
occupied-label centroid -- has the same blind spot one level up: a y-flip
PRESERVES the first moment whenever the field is symmetric about the flip
axis. It would pass exactly the case it exists to catch. So this file compares
located values at ASYMMETRIC probe points, the way the lattice test already
does, and asserts the asymmetry of its own fixture as a precondition -- because
otherwise the fixture silently becomes the load-bearing part and a later edit
disarms the check with nothing failing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from biofilm_openmc.observer import grid_geometry, read_layer, read_manifest
from biofilm_openmc.viewer import Grid, Layer, write_bundle

# Deliberately unequal on every axis and unequal under a flip of any one of
# them, so a permutation and a reflection are both detectable.
SHAPE = (4, 6, 5)
PROBES = ((0, 0, 0), (1, 4, 2), (3, 1, 4), (2, 5, 1))


def _world_centre(geom: dict, idx) -> np.ndarray:
    """Physical centre of one voxel, in cm, from the grid's own declaration."""
    origin = np.asarray(geom["origin_cm"], dtype=float)
    spacing = np.asarray(geom["spacing_cm"], dtype=float)
    return origin + (np.asarray(idx, dtype=float) + 0.5) * spacing


def _asymmetric_field(shape) -> np.ndarray:
    """A field whose value encodes its own index, so any permutation or
    reflection changes at least one probe. Same idea as the Julia return-leg
    test's `x + 10y + 100z`, which is the precedent for this shape."""
    x, y, z = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    return (x + 10.0 * y + 100.0 * z) + 1.0


def _bundle(tmp_path, dose=None, name="coreg.h5"):
    labels = np.zeros(SHAPE, np.int32)
    labels[1, 4, 2] = 7
    labels[3, 1, 4] = 9
    field = _asymmetric_field(SHAPE) if dose is None else dose
    grid = Grid("cpm_labels", SHAPE, (0.5, -1.0, 2.0), (0.25, 0.5, 0.125))
    dose_grid = Grid("dose_tally", SHAPE, (0.5, -1.0, 2.0), (0.25, 0.5, 0.125),
                     material_resolution_grid="cpm_labels")
    path = tmp_path / name
    write_bundle(path, [grid, dose_grid], [
        Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
              labels, background=0),
        Layer("dose_rate_Gy_s", "dose_tally", "Gy/s", "intensive", field),
    ], provenance={"reference_system_id": "synthetic",
                   "target_calibration": False,
                   "evidence_policy": "synthetic",
                   "openmc_version": "0.15.3"})
    return path


def test_the_fixture_is_asymmetric_on_every_axis():
    """THE PRECONDITION, STATED. Every assertion below is only meaningful
    because a flip changes the field. If the fixture ever becomes symmetric,
    the co-registration checks keep passing while detecting nothing -- the
    failure mode is silent, so it gets its own test rather than a comment.
    """
    field = _asymmetric_field(SHAPE)
    assert len(set(SHAPE)) == len(SHAPE), (
        f"{SHAPE} has a repeated extent, so an axis SWAP is undetectable")
    for axis in range(3):
        assert not np.array_equal(field, np.flip(field, axis=axis)), (
            f"the fixture is symmetric under a flip of axis {axis}, so a "
            "reflection on that axis would pass every check in this file")
        assert not np.array_equal(field, np.moveaxis(field, axis, (axis + 1) % 3)), (
            f"the fixture survives a permutation involving axis {axis}")


def test_dose_and_labels_agree_in_world_coordinates_at_asymmetric_probes(
        tmp_path):
    """The seam: two grids, one coordinate system, checked where a flip shows.

    NOT a centroid -- a y-flip preserves the first moment under symmetry, so a
    moment-based check passes the exact case it is for. These are located
    probes, and their world coordinates are recomputed from each grid's own
    declared `origin_cm`/`spacing_cm` rather than assumed equal, so an origin
    offset between the grids fails here even though both grids would still
    span the same extent.
    """
    path = _bundle(tmp_path)
    doc = read_manifest(path)
    lab_geom = grid_geometry(doc, "cpm_labels")
    dose_geom = grid_geometry(doc, "dose_tally")
    dose, _ = read_layer(path, "dose_rate_Gy_s")
    expected = _asymmetric_field(SHAPE)

    for idx in PROBES:
        assert np.allclose(_world_centre(lab_geom, idx),
                           _world_centre(dose_geom, idx)), (
            f"probe {idx} sits at different physical points on the two grids")
        assert dose[idx] == pytest.approx(expected[idx]), (
            f"dose at probe {idx} is {dose[idx]}, expected {expected[idx]}: "
            "the field is not where the grid says it is")


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_a_flipped_dose_field_is_caught(tmp_path, axis):
    """THE NEGATIVE CONTROL, on every axis rather than the one I happened to
    worry about. The y-flip is the interesting one because `RectLattice`
    genuinely has a y-inversion that the tally transform does not -- so a
    reimplementation copying the wrong convention produces exactly this input.
    """
    flipped = np.flip(_asymmetric_field(SHAPE), axis=axis)
    path = _bundle(tmp_path, dose=flipped, name=f"flip{axis}.h5")
    dose, _ = read_layer(path, "dose_rate_Gy_s")
    expected = _asymmetric_field(SHAPE)

    mismatched = [idx for idx in PROBES if dose[idx] != expected[idx]]
    assert mismatched, (
        f"a flip of axis {axis} changed nothing at any probe in {PROBES}; the "
        "probes are not asymmetric enough to detect it")


def test_a_shifted_grid_origin_is_caught(tmp_path):
    """The failure index arithmetic cannot see at all. Both grids keep their
    shape and spacing, so every index still resolves and every extent still
    matches -- only the physical placement differs. This is why the check is
    written in world coordinates rather than on indices.
    """
    labels = np.zeros(SHAPE, np.int32)
    labels[1, 4, 2] = 7
    good = Grid("cpm_labels", SHAPE, (0.5, -1.0, 2.0), (0.25, 0.5, 0.125))
    shifted = Grid("dose_tally", SHAPE, (0.75, -1.0, 2.0), (0.25, 0.5, 0.125))

    # The bundle refuses it outright -- one coordinate system means one origin.
    with pytest.raises(ValueError, match="origin|coordinate|volume"):
        write_bundle(tmp_path / "shift.h5", [good, shifted], [
            Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                  labels, background=0),
            Layer("dose_rate_Gy_s", "dose_tally", "Gy/s", "intensive",
                  _asymmetric_field(SHAPE)),
        ], provenance={"reference_system_id": "synthetic",
                       "target_calibration": False,
                       "evidence_policy": "synthetic",
                       "openmc_version": "0.15.3"})
