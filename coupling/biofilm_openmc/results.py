"""Result files, kept deliberately separate (review amendment 8):

- transport_result.h5 — the physical field only, keyed by
  transport_state_hash. Reusable across label changes.
- dose_attribution.h5 — the join of a (possibly cached) field with the
  CURRENT labels, keyed by BOTH hashes, so a cached field never falsely
  claims the full state matched.
"""

from __future__ import annotations

import h5py
import numpy as np

from .dose import DoseResult
from .snapshot import _logical_to_h5py, _h5py_to_logical


def write_transport_result(path, result: DoseResult, *, transport_hash: str,
                           openmc_version: str, nuclear_data_id: str,
                           batches: int, particles: int, seed: int) -> None:
    with h5py.File(path, "w") as f:
        f["mesh/dose_rate_mean_Gy_s"] = _logical_to_h5py(result.dose_rate_mean_Gy_s)
        f["mesh/dose_rate_sd_Gy_s"] = _logical_to_h5py(result.dose_rate_sd_Gy_s)
        f["mesh/rel_err"] = _logical_to_h5py(result.rel_err)
        f.attrs.update({
            "schema_version": 1,
            "logical_axis_order": "xyz",
            "dataset_axis_order_h5py": "zyx",
            "transport_state_hash": transport_hash,
            "source_rate_photons_per_s": result.source_rate,
            "openmc_version": openmc_version,
            "nuclear_data": nuclear_data_id,
            "batches": batches,
            "particles": particles,
            "seed": seed,
        })


def load_transport_result(path):
    with h5py.File(path, "r") as f:
        attrs = dict(f.attrs)
        result = DoseResult(
            _h5py_to_logical(f["mesh/dose_rate_mean_Gy_s"][()]),
            _h5py_to_logical(f["mesh/dose_rate_sd_Gy_s"][()]),
            _h5py_to_logical(f["mesh/rel_err"][()]),
            float(attrs["source_rate_photons_per_s"]))
        return result, attrs


def write_dose_attribution(path, aggregates: dict[str, dict[int, dict]], *,
                           transport_hash: str, label_hash: str,
                           uncertainty_method: str) -> None:
    """aggregates: {"lineage": {label: {...}}, "cell": {...}, ...}.
    uncertainty_method must be declared: "replicates" or
    "quadrature_independent_bins_approx" (labeled approximation)."""
    with h5py.File(path, "w") as f:
        f.attrs.update({
            "schema_version": 1,
            "transport_state_hash_used": transport_hash,
            "label_state_hash": label_hash,
            "uncertainty_method": uncertainty_method,
        })
        for kind, table in aggregates.items():
            labels = sorted(table)
            f[f"{kind}/label"] = np.array(labels, dtype=np.int64)
            for col in ("dose_rate_Gy_s", "mass_kg"):
                if labels and col in table[labels[0]]:
                    f[f"{kind}/{col}"] = np.array(
                        [table[l][col] for l in labels], dtype=float)
            if labels and "n_voxels" in table[labels[0]]:
                f[f"{kind}/n_voxels"] = np.array(
                    [table[l]["n_voxels"] for l in labels], dtype=np.int64)
            if labels and "sd_Gy_s" in table[labels[0]]:
                f[f"{kind}/sd_Gy_s"] = np.array(
                    [table[l]["sd_Gy_s"] for l in labels], dtype=float)
