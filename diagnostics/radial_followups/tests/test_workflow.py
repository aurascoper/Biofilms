"""Known-bad controls for pairing, provenance, exclusions and primary selection.

Temporary validation fixtures copy historical rows. They are not simulation
outputs at the scaled coefficients, and never leave pytest's temporary tree.
"""
import csv
import hashlib
import importlib.util
import io
import math
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("radial_followup_workflow", Path(__file__).parents[1] / "workflow.py")
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
REV = "a" * 40


@pytest.fixture
def inputs():
    blobs = {key: (w.ROOT / path).read_bytes() for key, path in w.FILES.items()}
    return {k: hashlib.sha256(v).hexdigest() for k, v in blobs.items()}, blobs


def fixture_text(job, inputs):
    hashes, blobs = inputs
    refs = w.baseline_rows(blobs["baseline_sha256"])
    axis, factor, arm, seed = (job[k] for k in ("axis", "factor", "arm", "seed"))
    meta = dict(format="radial_followups_v1", mode="production", revision=REV,
                axis=axis, factor=factor, arm=arm, seed=seed, N=40, parcels=6,
                at=400 if axis == "horizon" else 100, julia_version="test-fixture",
                machine="test-fixture", elapsed_seconds=0,
                baseline_checked=14 if factor == 1 or axis == "horizon" else 7, **hashes)
    for name, vector in w.vectors(blobs["model_sha256"].decode()).items():
        meta[name] = ";".join(str(0 if arm == "null" else x * (factor if axis == name else 1)) for x in vector)
    result = io.StringIO()
    for key, value in meta.items():
        result.write(f"# {key}={value}\n")
    writer = csv.DictWriter(result, fieldnames=w.COLUMNS); writer.writeheader()
    for mcs in ((0, 100, 400) if axis == "horizon" else (0, 100)):
        for sp in w.SPECIES:
            row = dict(refs[seed, arm, min(mcs, 100), sp]); row["mcs"] = str(mcs)
            writer.writerow(row)
    return result.getvalue()


def test_manifest_is_complete_and_baselines_precede_scaled_jobs(tmp_path):
    jobs = w.jobs(tmp_path)
    assert len(jobs) == 4352 == len({j["path"] for j in jobs})
    assert all(j["stage"] == "baseline" for j in jobs[:512])
    assert all(j["stage"] == "followup" for j in jobs[512:])
    assert {j["seed"] for j in jobs[:512]} == set(range(42, 298))


def test_good_archive_fixture_and_bad_row_controls(tmp_path, inputs):
    job = w.jobs(tmp_path)[0]
    path = Path(job["path"]); path.parent.mkdir(parents=True)
    good = fixture_text(job, inputs)
    path.write_text(good)
    assert len(w.validate_shard(path, REV, job, inputs=inputs)[1]) == 14
    lines = good.splitlines()
    first = next(i for i, line in enumerate(lines) if line.startswith("42,full,"))
    variants = {
        "duplicate": lines + [lines[first]],
        "missing": lines[:first] + lines[first + 1:],
        "wrong seed": [line.replace("42,full,", "43,full,") if i == first else line for i, line in enumerate(lines)],
        "wrong arm": [line.replace("42,full,", "42,null,") if i == first else line for i, line in enumerate(lines)],
        "wrong mcs": [line.replace(",0,1,", ",1,1,") if i == first else line for i, line in enumerate(lines)],
    }
    for name, changed in variants.items():
        path.write_text("\n".join(changed) + "\n")
        with pytest.raises(ValueError):
            w.validate_shard(path, REV, job, inputs=inputs)
    for index, value in ((7, "nan"), (7, "99.000000"), (5, "-1")):
        changed = lines[:]; row = changed[first].split(","); row[index] = value; changed[first] = ",".join(row)
        path.write_text("\n".join(changed) + "\n")
        with pytest.raises(ValueError):
            w.validate_shard(path, REV, job, inputs=inputs)


@pytest.mark.parametrize("before,after", [("# mode=production", "# mode=test"),
    ("# revision=" + REV, "# revision=" + "b" * 40), ("# baseline_checked=14", "# baseline_checked=0"),
    ("# N=40", "# N=20"), ("# at=100", "# at=400"), ("# beta=", "# wrong_beta=")])
def test_provenance_and_configuration_controls(tmp_path, inputs, before, after):
    job = w.jobs(tmp_path)[0]; path = tmp_path / "fixture.csv"
    good = fixture_text(job, inputs)
    assert before in good
    path.write_text(good.replace(before, after))
    with pytest.raises(ValueError):
        w.validate_shard(path, REV, job, inputs=inputs)


def test_signed_inclusive_crossing_and_missing_level_controls():
    levels = {f: dict(status="evaluated", mean=.10) for f in w.FACTORS}
    assert w.first_crossing(levels)["factor"] is None  # outward is not inward
    levels[100] = dict(status="evaluated", mean=-.05)
    assert w.first_crossing(levels)["factor"] == 100
    assert w.first_crossing(levels)["later_loss_of_crossing"]
    assert w.first_crossing(levels)["nonmonotone_mean"]
    incomplete = dict(levels); incomplete.pop(10)
    assert w.first_crossing(incomplete)["status"].startswith("not_evaluable")
    levels[10] = dict(status="not_evaluable_absent_species")
    assert w.first_crossing(levels)["status"].startswith("not_evaluable")


def test_absent_species_is_not_zero_radius_and_pairs_determine_uncertainty():
    row = lambda radius, count=6: dict(mean_r_over_R=str(radius), ncells=str(count))
    pairs = [(row(.20), row(.21)), (row(.80), row(.79))]
    summary = w.paired_summary(pairs)
    assert summary["mean"] == pytest.approx(0, abs=1e-15)
    assert summary["paired_sd"] == pytest.approx(math.sqrt(.0002))
    assert summary["standard_error"] == pytest.approx(.01)
    assert w.paired_summary([(row(0, 0), row(.21)), pairs[1]])["absent_pairs"] == 1
    with pytest.raises(ValueError):
        w.paired_summary(pairs[:1])


def test_analysis_wiring_refuses_missing_shard_and_mixed_runtime(tmp_path, inputs, monkeypatch):
    # A reduced TEST grid with archived rows; no new coefficient results are used.
    monkeypatch.setattr(w, "SEEDS", (42, 43)); monkeypatch.setattr(w, "FACTORS", (1, 10))
    hashes, blobs = inputs
    header = ",".join(w.COLUMNS)
    rows = [line for line in blobs["baseline_sha256"].decode().splitlines()
            if line.startswith(("42,", "43,"))]
    blobs = dict(blobs, baseline_sha256=(header + "\n" + "\n".join(rows) + "\n").encode())
    inputs = ({k: hashlib.sha256(v).hexdigest() for k, v in blobs.items()}, blobs)
    monkeypatch.setattr(w, "committed_inputs", lambda revision: inputs)
    jobs = w.jobs(tmp_path)
    for job in jobs:
        path = Path(job["path"]); path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(fixture_text(job, inputs))
    result = w.analyse(tmp_path, REV)
    assert result["sweeps"]["beta"]["primary"]["factor"] is None
    assert all(r["change_from_100"]["mean"] == 0 for r in result["horizon"].values())
    victim = Path(jobs[-1]["path"]); original = victim.read_text()
    victim.unlink()
    with pytest.raises(FileNotFoundError):
        w.analyse(tmp_path, REV)
    victim.write_text(original.replace("# machine=test-fixture", "# machine=other-fixture"))
    with pytest.raises(ValueError, match="mixed runtime"):
        w.analyse(tmp_path, REV)
