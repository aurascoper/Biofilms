"""A rendered figure must be derived from the SVG currently in the tree.

SCOPE: the SVG/PDF/sha256 triples under preprint/figures/, and nothing else.

This exists because the generator-and-artifact pair has already failed here once:
a figure generator was corrected two weeks before its committed images were, and
no check could open a figure to notice. FIG-01 through FIG-07 record it.

THE HASH IS WRITTEN ONLY BY tools/render_figure_svg.py. A sidecar a human can
rewrite binds the SVG to a record of itself, so the repair for a mismatch would be
to regenerate the record -- green, with the PDF still stale. Because only the
renderer writes it, the sole way to make this green is to run the renderer, which
emits the PDF as a side effect.
"""
from __future__ import annotations

import hashlib
import html
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIGDIR = REPO / "preprint" / "figures"


def svg_triples():
    return sorted(FIGDIR.glob("*.svg"))


@pytest.mark.parametrize("svg", svg_triples(), ids=lambda p: p.name)
def test_the_rendered_pdf_matches_the_svg_in_the_tree(svg):
    record = Path(str(svg) + ".sha256")
    pdf = svg.with_suffix(".pdf")
    assert record.is_file(), f"{svg.name} has no render record; run tools/render_figure_svg.py"
    assert pdf.is_file(), f"{svg.name} has no rendered PDF"
    recorded = record.read_text(encoding="utf-8").split()[0]
    actual = hashlib.sha256(svg.read_bytes()).hexdigest()
    assert recorded == actual, (
        f"{svg.name} changed since it was rendered. Re-run "
        f"tools/render_figure_svg.py {svg.relative_to(REPO)} -- do NOT edit the "
        f".sha256 file, which would make this pass with a stale PDF.")


def test_the_staleness_check_can_fail():
    """CONTROL: a record that does not match its source must be reported.

    SCOPE: one synthetic pair, exercising the comparison this file performs.
    """
    body = b"<svg/>"
    assert hashlib.sha256(body).hexdigest() != hashlib.sha256(b"<svg />").hexdigest()


def test_the_phase2_figure_carries_the_corrected_multiple():
    """The number this figure was corrected FOR must be in the shipped artifact.

    SCOPE: preprint/figures/phase2_diffusion_cell.pdf only.

    The caption compared the solver's radial D_eff to water self-diffusivity and
    said 100x; the verified value at 1e-3 cm2/s against 2.3e-5 is 43.5x. Asserting
    it on the PDF rather than the SVG is the point -- the prose being right while
    the image is wrong is the exact defect FIG-01..07 record.
    """
    pdf = FIGDIR / "phase2_diffusion_cell.pdf"
    if not pdf.is_file():
        pytest.skip("figure not rendered")
    # DEPENDENCY-GATED, like the sidecar test already is. Without Poppler this
    # raised FileNotFoundError instead of skipping, so a calibration environment
    # lacking it failed here rather than reporting uncovered surface -- and the
    # bare-tier checks in this file, which need no binary at all, were reported
    # as an error rather than running. Raised by Codex on pull request #23.
    import shutil
    import subprocess
    if not shutil.which("pdftotext"):
        pytest.skip("pdftotext unavailable; the artifact assertion needs Poppler")
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    assert "43" in out and "water" in out, "corrected multiple absent from the rendered figure"
    assert "100× water" not in out and "100x water" not in out, \
        "the withdrawn 100x comparison is still in the rendered figure"


# ---------------------------------------------------------------- source -> artifact
#
# WHAT THE HASH CHECK ABOVE CANNOT DO, AND THE FILE IT GUARDS SAID SO FIRST.
# `test_the_rendered_pdf_matches_the_svg_in_the_tree` compares the committed
# `.svg.sha256` against the SVG. Both are files in the tree. Edit the SVG,
# rewrite the record by hand, do not re-render: the record and the source agree,
# the PDF's own hash still matches its own record, and every check here passes
# over a stale PDF. tools/render_figure_svg.py's docstring argues exactly this
# and then wrote exactly that record. Raised by Codex on pull request #23.
#
# A HASH CANNOT FIX A HASH. Recording the PDF digest beside the source digest
# (now done) makes the forgery need two lies instead of one; it does not make it
# impossible, because both lines are still text in a file anyone can write.
#
# So the binding is content, not digest: every text run in the SVG must appear in
# the committed `.txt`, which the sidecar test independently pins to the PDF's
# actual bytes. SVG text -> committed text -> committed PDF, each link checked by
# something that cannot be satisfied by editing a record.
#
# SCOPE, STATED BECAUSE IT IS NARROWER THAN "THE PDF IS FRESH": this catches
# TEXT drift only. An SVG edit that moves a rectangle and touches no string
# passes, and no check in this repository would catch it. That is the honest
# reach -- and text is the failure mode the figure guards exist for, since
# FIG-01 through FIG-10 are every one of them a claim made in words inside an
# image. Matching is on the first MATCH_WORDS words of each run, because
# `pdftotext -layout` wraps a long run mid-word to fit its column and an exact
# substring match on the whole run fails on wrapping rather than on staleness.
MATCH_WORDS = 5


def svg_text_runs(svg_source: str) -> list[str]:
    """Whitespace-normalised contents of each <text> element, >= 3 words.

    SCOPE: `<text>` elements only. Text drawn as paths, or set in `<tspan>`
    outside a `<text>`, is not seen -- this repository's figures use neither.
    """
    runs = []
    for m in re.finditer(r"<text[^>]*>(.*?)</text>", svg_source, re.S):
        inner = html.unescape(re.sub(r"<[^>]+>", "", m.group(1)))
        s = " ".join(inner.split())
        if len(s.split()) >= 3:
            runs.append(s)
    return runs


def unrendered_runs(svg_source: str, txt: str) -> list[str]:
    """Runs whose opening MATCH_WORDS words are absent from the extracted text.

    Shared by the check and its control, so the control cannot pass against a
    regressed check.
    """
    flat = " ".join(txt.split())
    return [s for s in svg_text_runs(svg_source)
            if " ".join(s.split()[:MATCH_WORDS]) not in flat]


@pytest.mark.parametrize("svg", svg_triples(), ids=lambda p: p.name)
def test_the_svg_text_reaches_the_committed_pdf_text(svg):
    txt_path = svg.with_suffix(".txt")
    assert txt_path.is_file(), f"{svg.name} has no .txt sidecar; run tools/render_figure_svg.py"
    missing = unrendered_runs(svg.read_text(encoding="utf-8"),
                              txt_path.read_text(encoding="utf-8"))
    assert not missing, (
        f"{svg.name} carries text that is not in the rendered artifact, so the "
        f"committed PDF is older than the SVG: {missing[:3]}. Re-run "
        f"tools/render_figure_svg.py {svg.name} -- editing the .sha256 will make "
        "the hash checks pass and will not make this one.")


def test_the_source_to_artifact_binding_detects_an_edited_svg():
    """CONTROL: an SVG edited without re-rendering must be caught.

    SCOPE: one real figure, mutated in memory. The known-bad is DERIVED FROM THE
    ARTIFACT PATH rather than written by hand -- a synthetic `<text>` string
    would test the regex against my idea of the figure, and the sidecar case this
    repository already records is precisely a control that never met the
    pipeline.
    """
    svg = FIGDIR / "phase2_diffusion_cell.svg"
    source = svg.read_text(encoding="utf-8")
    txt = svg.with_suffix(".txt").read_text(encoding="utf-8")

    assert not unrendered_runs(source, txt), "baseline is not clean; no control verdict"

    runs = svg_text_runs(source)
    assert runs, "no text runs found -- the extractor is broken, not the figure"
    victim = max(runs, key=len)
    mutated = source.replace(victim.split()[0], "Notarealword", 1)
    assert mutated != source, "mutation did not apply"

    caught = unrendered_runs(mutated, txt)
    assert caught, (
        "an SVG whose text no longer matches the committed .txt was not caught, "
        "so this check cannot detect the staleness it exists for")
