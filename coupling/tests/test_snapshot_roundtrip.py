import h5py
import numpy as np
import pytest

from biofilm_openmc.snapshot import (SnapshotError, from_openmc_lattice_order,
                                     load_snapshot, to_openmc_lattice_order)

from conftest import make_snapshot_file


def test_schema_and_probe_verification(snapshot):
    assert snapshot.cell_id.shape == (8, 8, 8)
    assert snapshot.cell_id.dtype == np.int32
    assert snapshot.material_class_source == "absent"
    assert snapshot.cells["id"].tolist() == [1, 2, 3]
    # registry volumes consistent with the label arrays
    assert snapshot.cells["volume"].sum() == (snapshot.cell_id > 0).sum()


def test_corrupted_axis_order_is_refused(tmp_path):
    path = make_snapshot_file(tmp_path / "bad.h5")
    # swap two file axes of cell_id: probes must catch it
    with h5py.File(path, "r+") as f:
        arr = f["lattice/cell_id"][()]
        del f["lattice/cell_id"]
        f["lattice/cell_id"] = arr.transpose(1, 0, 2)
    with pytest.raises(SnapshotError, match="orientation probe failed"):
        load_snapshot(path)


def test_lattice_order_roundtrip():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 99, size=(4, 5, 6))
    assert np.array_equal(from_openmc_lattice_order(to_openmc_lattice_order(a)), a)


def test_lattice_order_semantics():
    # universes[z][y_index][x] with y_index increasing at DECREASING physical y
    nx, ny, nz = 3, 4, 5
    a = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    lat = to_openmc_lattice_order(a)
    assert lat.shape == (nz, ny, nx)
    for x, y, z in [(0, 0, 0), (2, 3, 4), (1, 2, 3)]:
        assert lat[z, ny - 1 - y, x] == a[x, y, z]
