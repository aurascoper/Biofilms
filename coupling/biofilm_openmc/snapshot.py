"""Transport-snapshot reader + the ONE place axis math lives.

Layout facts (docs/exchange_schema.md):
- Julia (column-major) writes logical (x,y,z); h5py (C-order) therefore sees
  shape (Nz,Ny,Nx) with element [z,y,x]. `_h5py_to_logical` recovers (x,y,z).
  This is a storage-layout fact, not a coordinate transform.
- OpenMC RectLattice.universes is indexed (z,y,x) with increasing y-index at
  DECREASING physical y. `to_openmc_lattice_order` performs exactly that
  mapping from a logical (x,y,z) array; `from_openmc_lattice_order` inverts it.

Every loaded snapshot is probe-verified: the file's asymmetric
`/orientation_probes` rows must match the recovered logical array.
"""

from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np


class SnapshotError(ValueError):
    pass


def _h5py_to_logical(arr: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(arr.transpose(2, 1, 0))


def _logical_to_h5py(arr: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(arr.transpose(2, 1, 0))


def to_openmc_lattice_order(arr_xyz: np.ndarray) -> np.ndarray:
    """Logical (x,y,z) -> RectLattice.universes order (z, y-inverted, x)."""
    return np.ascontiguousarray(arr_xyz.transpose(2, 1, 0)[:, ::-1, :])


def from_openmc_lattice_order(arr_zyx: np.ndarray) -> np.ndarray:
    """Inverse of `to_openmc_lattice_order`."""
    return np.ascontiguousarray(arr_zyx[:, ::-1, :].transpose(2, 1, 0))


@dataclass
class Snapshot:
    cell_id: np.ndarray          # logical (x,y,z), Int32
    species_id: np.ndarray
    lineage_id: np.ndarray
    generation: np.ndarray
    interior_mask: np.ndarray    # bool
    radiation_cpm: np.ndarray
    melanin: np.ndarray
    accumulated_dose_Gy: np.ndarray
    cells: dict                  # column name -> np.ndarray
    mcs: int
    physical_time_s: float
    label_state_hash: str
    material_class_source: str
    config_toml: str
    attrs: dict


def load_snapshot(path) -> Snapshot:
    with h5py.File(path, "r") as f:
        attrs = {k: f.attrs[k] for k in f.attrs}
        if int(attrs["schema_version"]) != 1:
            raise SnapshotError(f"unsupported schema_version {attrs['schema_version']}")
        if _as_str(attrs["logical_axis_order"]) != "xyz":
            raise SnapshotError("unexpected logical_axis_order")
        if int(attrs["coordinate_index_base"]) != 0:
            raise SnapshotError("orientation probes must be 0-based")

        def arr(name):
            return _h5py_to_logical(f[name][()])

        cell_id = arr("lattice/cell_id")

        # Julia wrote a (P,4) matrix column-major; the h5py (C-order) view is
        # ALWAYS its transpose — no shape heuristic (ambiguous when P == 4).
        probes = f["orientation_probes"][()].T
        for x, y, z, expect in probes:
            got = cell_id[int(x), int(y), int(z)]
            if int(got) != int(expect):
                raise SnapshotError(
                    f"orientation probe failed at (x,y,z)=({x},{y},{z}): "
                    f"expected {expect}, got {got} — axis order is broken, refusing to proceed")

        cells = {k: f[f"cells/{k}"][()] for k in f["cells"]}
        return Snapshot(
            cell_id=cell_id,
            species_id=arr("lattice/species_id"),
            lineage_id=arr("lattice/lineage_id"),
            generation=arr("lattice/generation"),
            interior_mask=arr("lattice/interior_mask").astype(bool),
            radiation_cpm=arr("fields/radiation_cpm"),
            melanin=arr("fields/melanin"),
            accumulated_dose_Gy=arr("dose/accumulated_Gy"),
            cells=cells,
            mcs=int(attrs["mcs"]),
            physical_time_s=float(attrs["physical_time_s"]),
            label_state_hash=_as_str(attrs["label_state_hash"]),
            material_class_source=_as_str(attrs["material_class_source"]),
            config_toml=_as_str(f["config_toml"][()]),
            attrs=attrs,
        )


def write_dose_field(path, dose_mean_xyz: np.ndarray, dose_sd_xyz: np.ndarray,
                     attrs: dict) -> None:
    """Write a dose field back in the Julia-readable layout (logical xyz in
    file dims, i.e. transposed for the C-order writer)."""
    with h5py.File(path, "w") as f:
        f["mesh/dose_rate_mean_Gy_s"] = _logical_to_h5py(dose_mean_xyz)
        f["mesh/dose_rate_sd_Gy_s"] = _logical_to_h5py(dose_sd_xyz)
        f.attrs["logical_axis_order"] = "xyz"
        f.attrs["dataset_axis_order_h5py"] = "zyx"
        for k, v in attrs.items():
            f.attrs[k] = v


def _as_str(v) -> str:
    return v.decode() if isinstance(v, bytes) else str(v)
