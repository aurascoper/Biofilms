"""Offline sensitivity driver (immutable-snapshot replay).

Replays a saved transport snapshot under alternative MATERIAL assumptions
(scenario override files deep-merged onto the base config), runs transport
per scenario (exact fingerprints — identical physical state is never re-run),
and evaluates the two-way-feedback gate:

  effect(scenario) = |mass-weighted mean dose - baseline| / baseline

  gate PASSED       effect > declared threshold AND effect > uncertainty
  gate FAILED       transport ran; effect below threshold or below noise
  gate NOT EVALUATED  config unset / transport stack or data unavailable

Comparison uses a dose floor + uncertainty mask + global metrics — never raw
max relative error in near-zero-dose voxels (review amendment 12).
"""

from __future__ import annotations

import argparse
import sys
import tomllib
from pathlib import Path

import numpy as np

from .config import Config, ConfigError, config_from_dict
from .dose import DoseResult, dose_from_statepoint
from .fingerprint import label_state_hash, transport_state_hash
from .lineage import aggregate_by_label
from .materials import voxel_mass_kg
from .results import write_transport_result
from .snapshot import Snapshot, load_snapshot

GATE_PASSED = "PASSED"
GATE_FAILED = "FAILED"
GATE_NOT_EVALUATED = "NOT EVALUATED"


def deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in override.items():
        out[k] = deep_merge(out[k], v) if (
            isinstance(v, dict) and isinstance(out.get(k), dict)) else v
    return out


def scenario_config(base_toml: str, override_path) -> Config:
    base = tomllib.loads(base_toml)
    with open(override_path, "rb") as fh:
        override = tomllib.load(fh)
    return config_from_dict(deep_merge(base, override),
                            raw=f"{base_toml}\n# merged: {override_path}")


def compare_fields(baseline: DoseResult, variant: DoseResult,
                   mass_kg: np.ndarray, dose_floor_Gy_s: float) -> dict:
    """Global effect metrics with floor + uncertainty mask (no raw max-rel-err
    on near-zero voxels)."""
    mask = (baseline.dose_rate_mean_Gy_s > dose_floor_Gy_s) & \
           np.isfinite(baseline.rel_err)
    m = mass_kg[mask]
    b = baseline.dose_rate_mean_Gy_s[mask]
    v = variant.dose_rate_mean_Gy_s[mask]
    sd = np.hypot(baseline.dose_rate_sd_Gy_s[mask],
                  variant.dose_rate_sd_Gy_s[mask])
    if m.sum() == 0:
        return {"masked_voxels": 0, "effect_rel": float("nan"),
                "noise_rel": float("nan"),
                "total_deposited_ratio": float("nan")}
    wmean = lambda x: float((x * m).sum() / m.sum())
    base_mean = wmean(b)
    return {
        "masked_voxels": int(mask.sum()),
        # mass-weighted mean absolute field difference, relative to baseline
        "effect_rel": wmean(np.abs(v - b)) / base_mean,
        # transport noise at the same aggregation
        "noise_rel": wmean(sd) / base_mean,
        # total deposited power ratio (global energy metric)
        "total_deposited_ratio": float((v * m).sum() / (b * m).sum()),
    }


def evaluate_gate(metrics_by_scenario: dict[str, dict],
                  effect_threshold: float) -> tuple[str, str]:
    effects = {k: m["effect_rel"] for k, m in metrics_by_scenario.items()}
    if not effects or any(np.isnan(e) for e in effects.values()):
        return GATE_NOT_EVALUATED, "no valid masked voxels to compare"
    worst = max(effects, key=lambda k: effects[k])
    eff = effects[worst]
    noise = metrics_by_scenario[worst]["noise_rel"]
    detail = (f"max effect {eff:.3%} (scenario '{worst}'), transport noise "
              f"{noise:.3%}, declared threshold {effect_threshold:.3%}")
    if eff > effect_threshold and eff > noise:
        return GATE_PASSED, detail
    return GATE_FAILED, detail


def _openmc_runner(outdir: Path):
    import openmc  # noqa: deferred import — gate reports NOT EVALUATED if absent

    def run(model, name):
        cwd = outdir / name
        cwd.mkdir(parents=True, exist_ok=True)
        sp_path = model.run(cwd=str(cwd), output=False)
        return openmc.StatePoint(sp_path)

    return run


def scan(snapshot: Snapshot, base_config: Config,
         scenarios: dict[str, Config], runner, outdir: Path, *,
         effect_threshold: float, dose_floor_Gy_s: float,
         nuclear_data_id: str, openmc_version: str) -> dict:
    """Core loop, runner-injectable for tests. Returns the report dict."""
    from .model import build_model

    results: dict[str, DoseResult] = {}
    hashes: dict[str, str] = {}
    seen: dict[str, str] = {}
    for name, cfg in {"baseline": base_config, **scenarios}.items():
        th = transport_state_hash(snapshot, cfg, nuclear_data_id)
        hashes[name] = th
        if th in seen:  # exact-match reuse only
            results[name] = results[seen[th]]
            continue
        sp = runner(build_model(snapshot, cfg), name)
        mass = voxel_mass_kg(snapshot, cfg)
        res = dose_from_statepoint(sp, mass, cfg.photons_per_second)
        results[name] = res
        seen[th] = name
        write_transport_result(
            outdir / f"transport_result_{name}.h5", res,
            transport_hash=th, openmc_version=openmc_version,
            nuclear_data_id=nuclear_data_id, batches=cfg.batches,
            particles=cfg.particles, seed=cfg.seed)

    mass = voxel_mass_kg(snapshot, base_config)
    metrics = {name: compare_fields(results["baseline"], res, mass,
                                    dose_floor_Gy_s)
               for name, res in results.items() if name != "baseline"}
    gate, detail = evaluate_gate(metrics, effect_threshold)
    return {
        "gate": gate,
        "gate_detail": detail,
        "metrics": metrics,
        "transport_hashes": hashes,
        "label_state_hash": label_state_hash(snapshot),
        "baseline_lineage_dose": aggregate_by_label(
            results["baseline"].dose_rate_mean_Gy_s, snapshot.lineage_id, mass),
    }


def main_scan(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="biofilm-dose-scan",
        description="Offline immutable-snapshot dose sensitivity scan")
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--scenario", action="append", default=[],
                    metavar="NAME=OVERRIDE.toml")
    ap.add_argument("--outdir", default="dose_scan_out")
    ap.add_argument("--effect-threshold", type=float, default=0.05)
    ap.add_argument("--dose-floor-Gy-s", type=float, default=0.0)
    ap.add_argument("--nuclear-data-id", default="unset")
    args = ap.parse_args(argv)

    def not_evaluated(reason: str) -> int:
        print(f"feedback gate: {GATE_NOT_EVALUATED} — {reason}")
        return 2

    try:
        with open(args.config, "rb") as fh:
            base_toml = fh.read().decode()
        base_config = config_from_dict(tomllib.loads(base_toml), base_toml)
        scenarios = {}
        for spec in args.scenario:
            name, _, path = spec.partition("=")
            scenarios[name] = scenario_config(base_toml, path)
    except ConfigError as e:
        return not_evaluated(f"physical config unset/invalid:\n{e}")

    try:
        import openmc
    except ImportError:
        return not_evaluated("openmc not importable (pinned conda env required)")
    import os
    xs = os.environ.get("OPENMC_CROSS_SECTIONS", "")
    if not xs or not Path(xs).exists():
        return not_evaluated("OPENMC_CROSS_SECTIONS not available")

    snapshot = load_snapshot(args.snapshot)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    report = scan(snapshot, base_config, scenarios, _openmc_runner(outdir),
                  outdir, effect_threshold=args.effect_threshold,
                  dose_floor_Gy_s=args.dose_floor_Gy_s,
                  nuclear_data_id=args.nuclear_data_id,
                  openmc_version=openmc.__version__)

    print(f"feedback gate: {report['gate']} — {report['gate_detail']}")
    for name, m in report["metrics"].items():
        print(f"  {name}: effect {m['effect_rel']:.3%}, noise {m['noise_rel']:.3%}, "
              f"deposited ratio {m['total_deposited_ratio']:.4f} "
              f"({m['masked_voxels']} masked voxels)")
    return 0 if report["gate"] != GATE_NOT_EVALUATED else 2


if __name__ == "__main__":
    sys.exit(main_scan())
