#!/usr/bin/env python3
"""Validate atomic shards and analyse the registered grid; no simulations here."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import re
import statistics
import subprocess
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASELINE = ROOT / "diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv"
FACTORS = (1, 10, 100, 1000, 10000, 100000, 1000000)
SEEDS = tuple(range(42, 298))
SPECIES = tuple(range(1, 8))
COLUMNS = ("seed", "arm", "mcs", "species", "species_name", "ncells", "volume",
           "mean_r", "mean_r_over_R", "mean_melanin")
FILES = {"model_sha256": "biofilms_potts.jl",
         "producer_sha256": "diagnostics/radial_followups/sweep.jl",
         "prediction_sha256": "diagnostics/radial_followups/PREDICTION.md",
         "baseline_sha256": str(BASELINE.relative_to(ROOT))}


def jobs(root):
    """All 4352 jobs, with both historical baselines first. Paths are absolute."""
    root = Path(root).resolve()
    out = []
    for seed in SEEDS:
        for axis, arm, rel in (("beta", "full", "beta/factor_1"),
                               ("baseline", "null", "baseline/null")):
            out.append(dict(path=str(root / rel / f"seed_{seed}.csv"), axis=axis,
                            arm=arm, factor=1, seed=seed, stage="baseline"))
    for axis, factor, seed in itertools.product(("beta", "alpha"), FACTORS, SEEDS):
        if axis == "beta" and factor == 1:
            continue
        out.append(dict(path=str(root / axis / f"factor_{factor}" / f"seed_{seed}.csv"),
                        axis=axis, arm="full", factor=factor, seed=seed, stage="followup"))
    for arm, seed in itertools.product(("full", "null"), SEEDS):
        out.append(dict(path=str(root / "horizon" / arm / f"seed_{seed}.csv"),
                        axis="horizon", arm=arm, factor=1, seed=seed, stage="followup"))
    return out


@lru_cache(maxsize=4)
def committed_inputs(revision):
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("full 40-character revision required")
    blobs = {key: subprocess.run(["git", "-C", str(ROOT), "show", f"{revision}:{path}"],
                                 check=True, capture_output=True).stdout
             for key, path in FILES.items()}
    return {key: hashlib.sha256(value).hexdigest() for key, value in blobs.items()}, blobs


def read_table(source):
    lines = source.splitlines()
    meta = {}
    for line in lines:
        if line.startswith("# ") and "=" in line:
            key, value = line[2:].split("=", 1)
            if key in meta:
                raise ValueError(f"duplicate metadata {key}")
            meta[key] = value
    reader = csv.DictReader(line for line in lines if line and not line.startswith("#"))
    if tuple(reader.fieldnames or ()) != COLUMNS:
        raise ValueError("wrong CSV columns")
    rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError("wrong CSV row width")
    return meta, rows


@lru_cache(maxsize=4)
def baseline_rows(blob):
    _, rows = read_table(blob.decode())
    keyed = {(int(r["seed"]), r["arm"], int(r["mcs"]), int(r["species"])): r for r in rows}
    expected = set(itertools.product(SEEDS, ("full", "no_rad", "no_mel", "null"), (0, 100), SPECIES))
    if len(keyed) != len(rows) or set(keyed) != expected:
        raise ValueError("committed baseline is incomplete or duplicated")
    return keyed


@lru_cache(maxsize=4)
def vectors(model):
    out = {}
    for name, symbol in (("beta", "β_ion"), ("alpha", "α_M_species")):
        match = re.search(rf"{symbol}::Vector\{{Float64\}}\s*=\s*\[(.*?)\]", model, re.S)
        if not match:
            raise ValueError(f"missing model vector {symbol}")
        body = re.sub(r"#.*", "", match[1])
        out[name] = [float(x.strip()) for x in body.split(",") if x.strip()]
        if len(out[name]) != 7:
            raise ValueError(f"wrong model vector length {symbol}")
    return out


def validate_shard(path, revision, expected_job=None, *, inputs=None):
    """Refuse malformed, unpaired or wrongly attributed files before resume/analysis."""
    hashes, blobs = committed_inputs(revision) if inputs is None else inputs
    meta, rows = read_table(Path(path).read_text())
    if meta.get("format") != "radial_followups_v1" or meta.get("mode") != "production":
        raise ValueError("not a production shard")
    for key, value in {"revision": revision, "N": "40", "parcels": "6", **hashes}.items():
        if meta.get(key) != value:
            raise ValueError(f"wrong provenance/configuration: {key}")
    axis, arm = meta.get("axis"), meta.get("arm")
    factor, seed = int(meta["factor"]), int(meta["seed"])
    if axis not in ("beta", "alpha", "baseline", "horizon") or seed not in SEEDS or factor not in FACTORS:
        raise ValueError("outside registered grid")
    if arm not in ("full", "null") or (axis in ("beta", "alpha") and arm != "full"):
        raise ValueError("wrong arm")
    if (axis == "baseline" and arm != "null") or (axis in ("baseline", "horizon") and factor != 1):
        raise ValueError("wrong baseline/horizon factor")
    if expected_job and any(str(meta[k]) != str(expected_job[k]) for k in ("axis", "arm", "factor", "seed")):
        raise ValueError("shard does not match requested job")
    at = 400 if axis == "horizon" else 100
    times = (0, 100, 400) if axis == "horizon" else (0, 100)
    if meta.get("at") != str(at) or not meta.get("julia_version") or not meta.get("machine"):
        raise ValueError("missing time/runtime provenance")
    expected_vectors = vectors(blobs["model_sha256"].decode())
    for name, vector in expected_vectors.items():
        want = [0.0 if arm == "null" else x * (factor if axis == name else 1) for x in vector]
        if [float(x) for x in meta.get(name, "").split(";")] != want:
            raise ValueError(f"wrong declared {name} vector")
    keyed = {}
    for row in rows:
        key = int(row["mcs"]), int(row["species"])
        if key in keyed:
            raise ValueError(f"duplicate row {key}")
        if int(row["seed"]) != seed or row["arm"] != arm:
            raise ValueError("row does not match seed/arm header")
        if any(int(row[k]) < 0 for k in ("ncells", "volume")):
            raise ValueError("negative counts")
        if not all(math.isfinite(float(row[k])) for k in ("mean_r", "mean_r_over_R", "mean_melanin")):
            raise ValueError("nonfinite observable")
        if abs(float(row["mean_r"]) / 20 - float(row["mean_r_over_R"])) > 0.0000006:
            raise ValueError("radius normalization mismatch")
        if arm == "null" and float(row["mean_melanin"]) != 0:
            raise ValueError("null has melanin")
        keyed[key] = row
    if set(keyed) != set(itertools.product(times, SPECIES)):
        raise ValueError("missing or unexpected rows")
    reference = baseline_rows(blobs["baseline_sha256"])
    checks = (0, 100) if factor == 1 or axis == "horizon" else (0,)
    for mcs, species in itertools.product(checks, SPECIES):
        if keyed[mcs, species] != reference[seed, arm, mcs, species]:
            raise ValueError(f"historical baseline mismatch at {mcs}, {species}")
    if meta.get("baseline_checked") != str(7 * len(checks)):
        raise ValueError("wrong baseline check count")
    return meta, keyed


def paired_summary(pairs):
    absent = sum(int(a["ncells"]) == 0 or int(b["ncells"]) == 0 for a, b in pairs)
    if absent:
        return dict(status="not_evaluable_absent_species", n_pairs=len(pairs), absent_pairs=absent)
    differences = [float(a["mean_r_over_R"]) - float(b["mean_r_over_R"]) for a, b in pairs]
    return difference_summary(differences)


def difference_summary(differences):
    if len(differences) < 2:
        raise ValueError("paired uncertainty requires at least two seeds")
    mean, sd = statistics.mean(differences), statistics.stdev(differences)
    se = sd / math.sqrt(len(differences))
    ci = [mean - 1.959963984540054 * se, mean + 1.959963984540054 * se]
    return dict(status="evaluated", n_pairs=len(differences), absent_pairs=0, mean=mean,
                paired_sd=sd, standard_error=se, ci95=ci, inward_criterion=mean <= -.05,
                interval_crosses_inward_boundary=ci[0] <= -.05 <= ci[1])


def first_crossing(levels):
    if set(levels) != set(FACTORS) or any(levels[f]["status"] != "evaluated" for f in FACTORS):
        return dict(status="not_evaluable_incomplete_or_absent_species", factor=None)
    crossed = [f for f in FACTORS if levels[f]["mean"] <= -.05]
    first = min(crossed) if crossed else None
    return dict(status="crossing" if first else "no_tested_factor_met_inward_criterion", factor=first,
                nonmonotone_mean=any(levels[b]["mean"] > levels[a]["mean"] for a, b in zip(FACTORS, FACTORS[1:])),
                later_loss_of_crossing=bool(first and any(levels[f]["mean"] > -.05 for f in FACTORS if f > first)))


def analyse(root, revision):
    inputs = committed_inputs(revision)
    reference = baseline_rows(inputs[1]["baseline_sha256"])
    records, runtimes = {}, set()
    for job in jobs(root):
        meta, data = validate_shard(job["path"], revision, job, inputs=inputs)
        runtimes.add((meta["julia_version"], meta["machine"]))
        records[job["axis"], job["factor"], job["arm"], job["seed"]] = data
    if len(runtimes) != 1:
        raise ValueError(f"mixed runtime provenance: {sorted(runtimes)}")
    result = dict(revision=revision, input_hashes=inputs[0], runtime=next(iter(runtimes)),
                  criterion="CN mean(full_scaled-null) <= -0.05 r/R; smallest tested factor",
                  uncertainty="normal 95% paired-mean intervals; not simultaneous or selection-adjusted",
                  sweeps={}, horizon={})
    for axis in ("beta", "alpha"):
        levels = {factor: {sp: paired_summary([(records[axis, factor, "full", s][100, sp],
                                              reference[s, "null", 100, sp]) for s in SEEDS])
                           for sp in SPECIES} for factor in FACTORS}
        result["sweeps"][axis] = dict(levels=levels, primary=first_crossing({f: levels[f][1] for f in FACTORS}))
    for sp in SPECIES:
        pairs = [(records["horizon", 1, "full", s][400, sp], records["horizon", 1, "null", s][400, sp]) for s in SEEDS]
        summary = paired_summary(pairs)
        if summary["status"] == "evaluated":
            changes = [float(a["mean_r_over_R"]) - float(b["mean_r_over_R"]) -
                       (float(reference[s, "full", 100, sp]["mean_r_over_R"]) -
                        float(reference[s, "null", 100, sp]["mean_r_over_R"]))
                       for s, (a, b) in zip(SEEDS, pairs)]
            summary["change_from_100"] = difference_summary(changes)
            summary["inside_declared_size_band"] = -.05 <= summary["mean"] <= .05
        result["horizon"][sp] = summary
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    manifest = commands.add_parser("jobs"); manifest.add_argument("root")
    validate = commands.add_parser("validate-shard"); validate.add_argument("path"); validate.add_argument("--revision", required=True)
    report = commands.add_parser("analyse"); report.add_argument("root"); report.add_argument("--revision", required=True); report.add_argument("--out", required=True)
    args = parser.parse_args()
    if args.command == "jobs":
        print(json.dumps(jobs(args.root), indent=2))
    elif args.command == "validate-shard":
        meta, rows = validate_shard(args.path, args.revision)
        print(json.dumps(dict(status="valid", rows=len(rows), seed=meta["seed"])))
    else:
        output = Path(args.out)
        if output.exists():
            raise ValueError("report already exists")
        result = analyse(args.root, args.revision)
        with output.open("x") as fh:
            json.dump(result, fh, indent=2, allow_nan=False); fh.write("\n")


if __name__ == "__main__":
    main()
