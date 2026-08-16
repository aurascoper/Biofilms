"""A0 resolution/history convergence sweep.

Answers one question, for a FIXED physical geometry and spectrum: how many
histories and what tally resolution are needed for a stable per-source
transport estimate? It does not identify the biofilm's production mesh — that
depends on the calibrated geometry, materials and gradients, and is the second
pass of the transport program.

Everything here is per source particle, so no source activity is needed
(`heating` is scored in eV per source particle). Fixed-history runs with
explicit zero-score accounting are used rather than tally triggers: a trigger
that ignores unscored bins can fire early precisely when sparsity is the
problem being measured.

TWO ACCOUNTING RULES THAT CHANGE THE ANSWER:

1. Zero-score fraction is computed INSIDE THE CYLINDER only. The mesh spans the
   circumscribing cube, so about 1 - pi/4 = 21.5% of bins are void by
   construction and must score zero. Counting them would inflate the sparsity
   metric by a fixed offset and make every mesh look worse than it is.

2. Cross-resolution comparison aggregates DEPOSITED ENERGY, never a mean of
   per-bin dose. Energy is extensive, so a fine result is summed onto the
   coarse grid (E_j = sum_{v in j} E_v) before comparing. Averaging dose would
   weight every bin equally regardless of mass.
"""

from __future__ import annotations

import numpy as np

from .dose import EV_TO_J
from .mesh import (coarsen_field, mesh_bin_volume_cm3,
                   phantom_mesh_extent_cm, upsample_field)


def in_cylinder_mask(config, dimension) -> np.ndarray:
    """Bins whose centre lies inside the water cylinder, logical (x,y,z).

    The complement is void: it must score zero, so including it in a sparsity
    metric measures the geometry rather than the statistics.
    """
    nx, ny, nz = (int(d) for d in dimension)
    ex, ey, _ = phantom_mesh_extent_cm(config)
    xs = (np.arange(nx) + 0.5) * ex / nx - ex / 2.0
    ys = (np.arange(ny) + 0.5) * ey / ny - ey / 2.0
    r2 = xs[:, None] ** 2 + ys[None, :] ** 2
    inside = r2 <= config.cylinder_radius_cm ** 2
    return np.broadcast_to(inside[:, :, None], (nx, ny, nz))


def deposited_energy_J(heating_eV: np.ndarray) -> np.ndarray:
    """Energy per source particle in each bin. Extensive — this is what may be
    summed across bins and across resolutions."""
    return np.asarray(heating_eV) * EV_TO_J


def point_metrics(heating_mean_eV: np.ndarray, heating_sd_eV: np.ndarray,
                  mass_kg: np.ndarray, mask: np.ndarray, *,
                  e_src_eV: float, heating_floor_frac: float = 0.0) -> dict:
    """Convergence statistics for one run. `mask` selects the region of
    interest (bins inside the cylinder)."""
    mean = np.asarray(heating_mean_eV, dtype=float)
    sd = np.asarray(heating_sd_eV, dtype=float)

    roi = mean[mask]
    scoring = roi > 0
    total_eV = float(mean.sum())

    floor = heating_floor_frac * (roi[scoring].mean() if scoring.any() else 0.0)
    above = roi > floor
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(mean > 0, sd / mean, np.inf)
    rel_roi = rel[mask][above]
    finite = rel_roi[np.isfinite(rel_roi)]

    dose = deposited_energy_J(mean) / mass_kg          # Gy per source particle
    m = mass_kg[mask]
    return {
        "roi_bins": int(mask.sum()),
        "total_heating_eV_per_src": total_eV,
        # vacuum boundaries: the fraction of emitted energy stopped in the phantom
        "absorbed_fraction": total_eV / e_src_eV,
        "zero_score_fraction_roi": float((~scoring).mean()) if roi.size else float("nan"),
        "heating_floor_eV": float(floor),
        "bins_above_floor": int(above.sum()),
        "median_rel_err": float(np.median(finite)) if finite.size else float("nan"),
        "p90_rel_err": float(np.percentile(finite, 90)) if finite.size else float("nan"),
        "max_rel_err": float(finite.max()) if finite.size else float("nan"),
        "mass_weighted_mean_dose_Gy_per_src":
            float((dose[mask] * m).sum() / m.sum()) if m.sum() else float("nan"),
        "total_mass_kg": float(mass_kg.sum()),
    }


def field_difference(coarse_heating_eV: np.ndarray, fine_heating_eV: np.ndarray,
                     factor: int, mass_kg_coarse: np.ndarray,
                     mask_coarse: np.ndarray) -> float:
    """Mass-weighted relative field difference between a coarse run and a finer
    one, after aggregating the fine DEPOSITED ENERGY onto the coarse grid.

    Returns the mass-weighted mean |dose_coarse - dose_fine| divided by the
    mass-weighted mean fine dose, over the region of interest.
    """
    fine_on_coarse = coarsen_field(np.asarray(fine_heating_eV, dtype=float), factor)
    d_coarse = deposited_energy_J(coarse_heating_eV) / mass_kg_coarse
    d_fine = deposited_energy_J(fine_on_coarse) / mass_kg_coarse
    m = mass_kg_coarse[mask_coarse]
    if m.sum() == 0:
        return float("nan")
    ref = float((d_fine[mask_coarse] * m).sum() / m.sum())
    if ref == 0:
        return float("nan")
    diff = float((np.abs(d_coarse - d_fine)[mask_coarse] * m).sum() / m.sum())
    return diff / ref


def resolution_loss(coarse_heating_eV: np.ndarray, coarse_mass_kg: np.ndarray,
                    factor: int, ref_heating_eV: np.ndarray,
                    ref_mass_kg: np.ndarray, ref_mask: np.ndarray) -> float:
    """Information destroyed by tallying coarsely, measured on the FINE grid.

    `field_difference` cannot answer this. Because the meshes are nested,
    summing the reference's energy onto a coarse grid reproduces that grid's
    exact true answer, so the discretization error there is zero by
    construction and only Monte Carlo noise survives. The resolution question
    is what the coarse field fails to represent, which is only visible where
    the structure lives: broadcast the coarse DOSE RATE (intensive, so each
    fine bin inherits its coarse bin's rate) back onto the reference grid and
    compare mass-weighted there.

    This grows with the factor as a real gradient gets smeared flat.
    """
    d_coarse = deposited_energy_J(coarse_heating_eV) / coarse_mass_kg
    d_up = upsample_field(d_coarse, factor)
    d_ref = deposited_energy_J(ref_heating_eV) / ref_mass_kg
    if d_up.shape != d_ref.shape:                       # pragma: no cover
        raise ValueError(f"upsampled {d_up.shape} != reference {d_ref.shape}")
    m = ref_mass_kg[ref_mask]
    if m.sum() == 0:
        return float("nan")
    ref_mean = float((d_ref[ref_mask] * m).sum() / m.sum())
    if ref_mean == 0:
        return float("nan")
    diff = float((np.abs(d_up - d_ref)[ref_mask] * m).sum() / m.sum())
    return diff / ref_mean


def seed_spread(values) -> float:
    """Relative spread of an aggregate across independent seeds.

    Independent replicates, not per-bin standard deviations: a per-bin sigma
    does not capture the covariance that matters for a spatial aggregate.
    """
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if v.size < 2 or v.mean() == 0:
        return float("nan")
    return float(v.std(ddof=1) / v.mean())


CSV_COLUMNS = (
    "coarsening_factor", "mesh_nx", "mesh_ny", "mesh_nz", "bin_volume_cm3",
    "batches", "particles", "histories", "seed",
    "roi_bins", "total_heating_eV_per_src", "absorbed_fraction",
    "zero_score_fraction_roi", "heating_floor_eV", "bins_above_floor",
    "median_rel_err", "p90_rel_err", "max_rel_err",
    "mass_weighted_mean_dose_Gy_per_src", "total_mass_kg",
    "field_diff_vs_reference", "resolution_loss_vs_reference",
    "runtime_s", "histories_per_s",
    "statepoint_bytes", "peak_child_rss_kb_cumulative",
)


# How a sweep CSV may be compared against a committed one. A byte-diff is the
# wrong instrument: four of these columns are machine-specific, and the A0
# report says so of its own runtimes. Partitioning them is what makes "A0 is
# unchanged" a checkable claim rather than a hope.
DETERMINISTIC_COLUMNS = (
    "coarsening_factor", "mesh_nx", "mesh_ny", "mesh_nz", "bin_volume_cm3",
    "batches", "particles", "histories", "seed", "roi_bins", "heating_floor_eV",
)
# Seeded Monte Carlo output. Bit-identical only for the same OpenMC build,
# nuclear data and thread count; compare within a declared tolerance otherwise.
STOCHASTIC_COLUMNS = (
    "total_heating_eV_per_src", "absorbed_fraction", "zero_score_fraction_roi",
    "bins_above_floor", "median_rel_err", "p90_rel_err", "max_rel_err",
    "mass_weighted_mean_dose_Gy_per_src", "total_mass_kg",
)
# Computed against the sweep's OWN reference row, so they are a property of the
# sweep's composition and not of the row. A sweep over factors 1,2 and a sweep
# over 1,2,4,8,16 give the same row different values here, which is correct and
# is why they cannot be compared across sweeps that do not share a reference.
REFERENCE_RELATIVE_COLUMNS = (
    "field_diff_vs_reference", "resolution_loss_vs_reference",
)
# Machine-specific. Never compared: they say what the host was, not what the
# physics is.
ENVIRONMENTAL_COLUMNS = (
    "runtime_s", "histories_per_s", "statepoint_bytes",
    "peak_child_rss_kb_cumulative",
)


def compare_sweeps(reference_rows, candidate_rows, *, rel_tol: float = 0.02,
                   same_reference: bool = False) -> list[str]:
    """Differences that matter between two sweep CSVs. Empty means unchanged.

    Deterministic columns must match exactly; stochastic ones within `rel_tol`;
    environmental ones are ignored entirely. The reference-relative columns are
    compared only when `same_reference` says the two sweeps were run over the
    same factor set — otherwise they legitimately differ and flagging them would
    train the reader to ignore this function.
    """
    problems: list[str] = []
    key = lambda r: (r["coarsening_factor"], r["histories"], r["seed"])  # noqa: E731
    ref = {key(r): r for r in reference_rows}
    cand = {key(r): r for r in candidate_rows}

    for missing in sorted(set(ref) - set(cand)):
        problems.append(f"row {missing} is missing from the candidate sweep")
    for added in sorted(set(cand) - set(ref)):
        problems.append(f"row {added} is not in the reference sweep")

    for k in sorted(set(ref) & set(cand)):
        for col in DETERMINISTIC_COLUMNS:
            a, b = ref[k].get(col, ""), cand[k].get(col, "")
            if str(a) != str(b):
                problems.append(f"row {k}: {col} changed {a!r} -> {b!r}")
        compared = STOCHASTIC_COLUMNS + (
            REFERENCE_RELATIVE_COLUMNS if same_reference else ())
        for col in compared:
            a, b = ref[k].get(col, ""), cand[k].get(col, "")
            if a in ("", None) or b in ("", None):
                continue
            a, b = float(a), float(b)
            if a == 0.0:
                if b != 0.0:
                    problems.append(f"row {k}: {col} moved off zero -> {b!r}")
            elif abs(b - a) / abs(a) > rel_tol:
                problems.append(
                    f"row {k}: {col} moved {a!r} -> {b!r} "
                    f"({abs(b - a) / abs(a):.1%} > {rel_tol:.0%})")
    return problems


def row_for(config, factor: int, dimension, seed: int, metrics: dict, *,
            runtime_s: float, statepoint_bytes: int, peak_rss_kb: int,
            field_diff: float, resolution_loss_: float) -> dict:
    histories = config.batches * config.particles
    row = {
        "coarsening_factor": factor,
        "mesh_nx": dimension[0], "mesh_ny": dimension[1], "mesh_nz": dimension[2],
        "bin_volume_cm3": mesh_bin_volume_cm3(phantom_mesh_extent_cm(config),
                                              dimension),
        "batches": config.batches, "particles": config.particles,
        "histories": histories, "seed": seed,
        "field_diff_vs_reference": field_diff,
        "resolution_loss_vs_reference": resolution_loss_,
        "runtime_s": round(runtime_s, 3),
        "histories_per_s": round(histories / runtime_s) if runtime_s else "",
        "statepoint_bytes": statepoint_bytes,
        "peak_child_rss_kb_cumulative": peak_rss_kb,
    }
    row.update(metrics)
    return {c: row.get(c, "") for c in CSV_COLUMNS}
