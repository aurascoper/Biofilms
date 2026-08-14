"""Mass-weighted dose aggregation by lineage/cell/generation labels.

Never unweighted-averages voxel doses (voxel masses differ), and never sums
per-bin variances as if independent: mesh-bin covariances are unavailable
(one history scores several bins), so per-lineage uncertainty comes from
independent replicate transport runs. A quadrature approximation is offered
but explicitly labeled as assuming independent bins.
"""

from __future__ import annotations

import numpy as np


def aggregate_by_label(dose_Gy_s: np.ndarray, labels: np.ndarray,
                       mass_kg: np.ndarray) -> dict[int, dict]:
    """Mass-weighted mean dose rate per label value (label 0 = unlabeled,
    skipped). Exact: sum(m_v * D_v) / sum(m_v) over each label's voxels."""
    lab = labels.ravel()
    d = dose_Gy_s.ravel()
    m = mass_kg.ravel()
    sel = lab > 0
    lab, d, m = lab[sel], d[sel], m[sel]
    out: dict[int, dict] = {}
    if lab.size == 0:
        return out
    energy = np.bincount(lab, weights=d * m)
    mass = np.bincount(lab, weights=m)
    count = np.bincount(lab)
    for label in np.nonzero(mass)[0]:
        out[int(label)] = {
            "dose_rate_Gy_s": float(energy[label] / mass[label]),
            "mass_kg": float(mass[label]),
            "n_voxels": int(count[label]),
        }
    return out


def replicate_uncertainty(per_replicate: list[dict[int, dict]]) -> dict[int, dict]:
    """Primary uncertainty path: aggregate each independent replicate run,
    then report mean and sample-sd across replicates per label."""
    if len(per_replicate) < 2:
        raise ValueError("need >= 2 independent replicates for an uncertainty estimate")
    labels = sorted(set().union(*[set(r) for r in per_replicate]))
    out = {}
    for label in labels:
        vals = np.array([r[label]["dose_rate_Gy_s"]
                         for r in per_replicate if label in r])
        out[label] = {
            "dose_rate_Gy_s": float(vals.mean()),
            "sd_Gy_s": float(vals.std(ddof=1)) if vals.size > 1 else float("nan"),
            "n_replicates": int(vals.size),
        }
    return out


def quadrature_uncertainty_approx(dose_sd_Gy_s: np.ndarray, labels: np.ndarray,
                                  mass_kg: np.ndarray) -> dict[int, float]:
    """APPROXIMATION — assumes independent mesh bins, which OpenMC does not
    guarantee (shared histories). Labeled as such wherever reported; prefer
    `replicate_uncertainty`."""
    lab = labels.ravel()
    sd = dose_sd_Gy_s.ravel()
    m = mass_kg.ravel()
    sel = lab > 0
    lab, sd, m = lab[sel], sd[sel], m[sel]
    var = np.bincount(lab, weights=(sd * m) ** 2)
    mass = np.bincount(lab, weights=m)
    return {int(l): float(np.sqrt(var[l]) / mass[l]) for l in np.nonzero(mass)[0]}
