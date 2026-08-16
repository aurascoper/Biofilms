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


def provenance_attrs(config, *, config_sha256: str = "",
                     snapshot_sha256: str = "", repo_root: str | None = None) -> dict:
    """What a result must carry to be reproducible from its own metadata.

    `target_calibration = false` alone cannot tell an A0 benchmark from a public
    surrogate, a synthetic fixture, an uncalibrated Reference D or a sensitivity
    case — and a file that cannot say which of those it is will eventually be
    read as whichever one the reader expected.
    """
    from physical_contract import git_provenance

    attrs = {
        "reference_system_id": config.reference_system_id,
        "system_provenance": config.system_provenance,
        "evidence_policy": config.evidence_policy,
        "execution_class": config.execution_class,
        "target_calibration": bool(config.target_calibration),
        "model_kind": config.model_kind,
        "config_sha256": config_sha256,
        "snapshot_sha256": snapshot_sha256,
    }
    for key, value in git_provenance(repo_root).items():
        attrs[key] = "" if value is None else value
    return attrs


def write_transport_result(path, result: DoseResult, *, transport_hash: str,
                           openmc_version: str, nuclear_data_id: str,
                           batches: int, particles: int, seed: int,
                           provenance: dict | None = None,
                           extra: dict | None = None) -> None:
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
        f.attrs.update(provenance or {})
        # Mesh dimension, coarsening factor, spectrum hash, dose hash,
        # material_volume_n_samples — whatever the caller had that the writer
        # cannot derive from a DoseResult.
        f.attrs.update(extra or {})


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
