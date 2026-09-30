"""Scope and provenance checks for the two figures corrected after 24cd3ec.

The negative controls are the actual PDFs committed at 24cd3ec, copied without
editing to fixtures/presend_24cd3ec. They enter through the same PDF extraction
as the manuscript figures. These checks do not assert absence of parcel motion.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest

REPO = Path(__file__).resolve().parents[2]
NAMES = ("fig1_radial_ensemble", "fig_consortium_seed42")
WITHDRAWN = {
    NAMES[0]: ("neither radiation-derived term moves a parcel", "switching a term off moves nothing"),
    NAMES[1]: ("every reported number comes from the 256-seed ensemble", "commit de8d88e-dirty"),
}
REQUIRED = {
    NAMES[0]: ("small ensemble mean-radius shifts across four paired arms",
               "paired mean-radius contrasts", "90% intervals"),
    NAMES[1]: ("figure 1 reports the 256-seed ensemble; these frames show seed 42",
               "frame archive", "render revision"),
}


def figure_findings(root):
    """Read both PDFs through Poppler; return each missing or withdrawn claim."""
    if not shutil.which("pdftotext"):
        pytest.skip("pdftotext unavailable; rendered claims are not checked")
    findings = set()
    for name in NAMES:
        pdf = root / f"{name}.pdf"
        assert pdf.is_file() and pdf.stat().st_size, f"missing or empty {pdf}"
        text = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
        assert text.strip(), f"empty PDF extraction: {pdf}"
        normal = " ".join(text.lower().split())
        findings.update((name, "withdrawn", phrase) for phrase in WITHDRAWN[name] if phrase in normal)
        findings.update((name, "missing", phrase) for phrase in REQUIRED[name] if phrase not in normal)
    return findings


def test_rendered_figures_have_bounded_claims_and_separate_provenance():
    assert figure_findings(REPO / "preprint" / "figures") == set()


def test_control_rejects_both_actual_pre_correction_pdfs():
    findings = figure_findings(REPO / "calibration" / "tests" / "fixtures" / "presend_24cd3ec")
    expected = {(name, "withdrawn", phrase) for name, phrases in WITHDRAWN.items() for phrase in phrases}
    expected |= {(name, "missing", phrase) for name, phrases in REQUIRED.items() for phrase in phrases}
    assert findings == expected


def render_helpers(repo):
    """Exercise production stdlib helpers without importing optional plotting packages."""
    source = (REPO / "tools" / "render_consortium.py").read_text()
    names = {"manifest", "clean_render_record", "external_output", "provenance"}
    nodes = [node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == names
    namespace = {"REPO": repo, "Path": Path, "subprocess": subprocess, "hashlib": hashlib}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "render_consortium.py", "exec"), namespace)
    return namespace


@pytest.fixture
def clean_render_repo(tmp_path):
    repo = tmp_path / "source"
    repo.mkdir()
    for relative in ("tools/render_consortium.py", "tools/render_label_trajectory.py",
                     "calibration/tests/test_claims_ledger.py",
                     "diagnostics/melanin_ensemble/consortium_seed42/manifest.txt",
                     "diagnostics/melanin_ensemble/consortium_seed42/full_mcs000.csv.gz"):
        dest = repo / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / relative, dest)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                    "-c", "commit.gpgsign=false", "commit", "-qm", "Fixture"], check=True)
    return repo, render_helpers(repo)


def test_provenance_preserves_original_run_and_names_archive_and_renderer(clean_render_repo):
    repo, helpers = clean_render_repo
    frames = repo / "diagnostics/melanin_ensemble/consortium_seed42"
    record = helpers["clean_render_record"](frames)
    assert record["original_run"]["git_commit"] == "de8d88e-dirty"
    assert record["render_git_commit"] == subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    assert record["frame_archive_tree"] == subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD:diagnostics/melanin_ensemble/consortium_seed42"], text=True).strip()
    caption = helpers["provenance"](helpers["manifest"](frames), record)
    assert f"frame archive {record['frame_archive_tree'][:10]}" in caption
    assert f"render revision {record['render_git_commit'][:10]}" in caption
    assert "de8d88e-dirty" not in caption
    for field in ("source_sha256", "input_sha256"):
        assert record[field]
        assert all(hashlib.sha256((repo / p).read_bytes()).hexdigest() == digest
                   for p, digest in record[field].items())


@pytest.mark.parametrize("staged", [False, True])
def test_control_dirty_tracked_source_refuses_clean_render_revision(clean_render_repo, staged):
    repo, helpers = clean_render_repo
    source = repo / "tools/render_consortium.py"
    before = source.read_bytes()
    source.write_bytes(before + b"\n# changed source\n")
    assert source.read_bytes() != before
    if staged:
        subprocess.run(["git", "-C", str(repo), "add", str(source)], check=True)
    with pytest.raises(SystemExit, match="tracked source files differ"):
        helpers["clean_render_record"](repo / "diagnostics/melanin_ensemble/consortium_seed42")


def test_control_untracked_frame_refuses_archive_claim(clean_render_repo):
    repo, helpers = clean_render_repo
    frames = repo / "diagnostics/melanin_ensemble/consortium_seed42"
    (frames / "unarchived.csv").write_text("not committed\n")
    with pytest.raises(SystemExit, match="exactly its committed files"):
        helpers["clean_render_record"](frames)


def test_control_output_inside_checkout_is_refused(tmp_path):
    helpers = render_helpers(REPO)
    with pytest.raises(SystemExit, match="outside the source checkout"):
        helpers["external_output"](REPO / "preprint/figures")
    assert helpers["external_output"](tmp_path) == tmp_path.resolve()


def test_consortium_receipt_matches_shipped_sources_inputs_and_outputs():
    path = REPO / "preprint/figures/fig_consortium_seed42.render.json"
    assert path.is_file(), "run render_consortium.py from a clean commit before shipping the figures"
    record = json.loads(path.read_text())
    original = render_helpers(REPO)["manifest"](REPO / "diagnostics/melanin_ensemble/consortium_seed42")
    assert record["original_run"] == {key: original[key] for key in (
        "git_commit", "julia_version", "machine", "seed", "N", "parcels_per_species", "n_mcs", "arms")}
    assert len(record["render_git_commit"]) == 40 and "dirty" not in record["render_git_commit"]
    assert len(record["frame_archive_tree"]) == 40
    assert record["frame_archive_tree"] == subprocess.check_output(
        ["git", "-C", str(REPO), "rev-parse", "HEAD:diagnostics/melanin_ensemble/consortium_seed42"], text=True).strip()
    assert set(record["source_sha256"]) == {
        "tools/render_consortium.py", "tools/render_label_trajectory.py", "calibration/tests/test_claims_ledger.py"}
    expected_outputs = {f"preprint/figures/{NAMES[1]}.{ext}" for ext in ("pdf", "txt", "png", "sha256")}
    expected_outputs |= {f"assets/consortium_seed42_full_vs_null.{ext}"
                         for ext in ("mp4", "gif", "mp4.sha256", "gif.sha256")}
    assert set(record["output_sha256"]) == expected_outputs
    expected_inputs = {p.relative_to(REPO).as_posix()
                       for p in (REPO / "diagnostics/melanin_ensemble/consortium_seed42").glob("*") if p.is_file()}
    assert set(record["input_sha256"]) == expected_inputs
    for field in ("source_sha256", "input_sha256", "output_sha256"):
        assert record[field]
        assert all(hashlib.sha256((REPO / p).read_bytes()).hexdigest() == digest
                   for p, digest in record[field].items()), field
    text = " ".join((REPO / "preprint/figures/fig_consortium_seed42.txt").read_text().split())
    assert f"frame archive {record['frame_archive_tree'][:10]}" in text
    assert f"render revision {record['render_git_commit'][:10]}" in text
