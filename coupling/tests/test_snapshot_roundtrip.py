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


# The CPM run seed. THREE states, because a file with no cpm_seed_source predates
# the field, and that silence is not a statement by the writer.


def test_a_declared_cpm_seed_is_read(tmp_path):
    snap = load_snapshot(make_snapshot_file(tmp_path / "seeded.h5", cpm_seed=42))
    assert snap.cpm_seed_source == "declared"
    assert snap.cpm_seed == 42


def test_a_writer_with_no_seed_declares_absent(tmp_path):
    snap = load_snapshot(make_snapshot_file(tmp_path / "unseeded.h5"))
    assert snap.cpm_seed_source == "absent"
    assert snap.cpm_seed is None


def test_a_file_predating_the_field_reads_as_unrecorded(tmp_path):
    path = make_snapshot_file(tmp_path / "old.h5", write_seed_attrs=False)
    snap = load_snapshot(path)
    assert snap.cpm_seed_source == "unrecorded"
    assert snap.cpm_seed is None


def test_a_declaration_without_a_seed_is_refused(tmp_path):
    path = make_snapshot_file(tmp_path / "bad_declared.h5", cpm_seed=7)
    with h5py.File(path, "r+") as f:
        del f.attrs["cpm_seed"]
    with pytest.raises(SnapshotError, match="no cpm_seed attribute"):
        load_snapshot(path)


def test_a_seed_without_a_declaration_is_refused(tmp_path):
    path = make_snapshot_file(tmp_path / "bad_absent.h5")
    with h5py.File(path, "r+") as f:
        f.attrs["cpm_seed"] = 7
    with pytest.raises(SnapshotError, match="cpm_seed_source is 'absent'"):
        load_snapshot(path)


def test_an_unknown_seed_source_is_refused(tmp_path):
    path = make_snapshot_file(tmp_path / "bad_source.h5")
    with h5py.File(path, "r+") as f:
        f.attrs["cpm_seed_source"] = "maybe"
    with pytest.raises(SnapshotError, match="unknown cpm_seed_source"):
        load_snapshot(path)


def test_a_stray_seed_on_a_file_predating_the_field_is_refused(tmp_path):
    # The control that proves "unrecorded" is not a free pass. A seed with no
    # declaration beside it still refuses, whichever way the declaration is missing.
    path = make_snapshot_file(tmp_path / "bad_old.h5", write_seed_attrs=False)
    with h5py.File(path, "r+") as f:
        f.attrs["cpm_seed"] = 7
    with pytest.raises(SnapshotError, match="cpm_seed_source is 'unrecorded'"):
        load_snapshot(path)
