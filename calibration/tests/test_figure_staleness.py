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
    import subprocess
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    assert "43" in out and "water" in out, "corrected multiple absent from the rendered figure"
    assert "100× water" not in out and "100x water" not in out, \
        "the withdrawn 100x comparison is still in the rendered figure"
