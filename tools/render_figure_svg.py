#!/usr/bin/env python3
"""Render an SVG figure to PDF and write the source hash THE CONVERSION SAW.

WHY THE HASH IS WRITTEN HERE AND NOWHERE ELSE. An SVG source-of-truth plus a
committed derived PDF is the generator-and-artifact pair that already failed in
this repository: the generator was corrected two weeks before the committed
images were, and nothing could open a figure to notice.

A plain `svg.sha256` sidecar does NOT fix that. Compared against the SVG it binds
the SVG to a record of itself -- edit the SVG, the hash mismatches, and the
obvious repair is to regenerate the hash. That is green again with the PDF still
stale, so the check meant to catch staleness has been satisfied without catching
it.

So the hash file is written ONLY by this script, as a side effect of actually
converting. The only way to make the guard green is to run the conversion, which
produces the PDF. This is contract_csv.jl's compare-never-regenerate reasoning
inverted: there the fixture must never be rewritten by the thing it checks; here
the record must be writable only by the generator.

    python3 tools/render_figure_svg.py preprint/figures/phase2_diffusion_cell.svg
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

CONVERTERS = [
    ("rsvg-convert", lambda src, dst: ["rsvg-convert", "-f", "pdf", "-o", str(dst), str(src)]),
]
CAIRO_PY = Path("/home/aurascoper/.local/share/mamba/envs/openmc-biofilms/bin/python")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render(svg: Path) -> Path:
    pdf = svg.with_suffix(".pdf")
    import shutil
    for name, argv in CONVERTERS:
        if shutil.which(name):
            subprocess.run(argv(svg, pdf), check=True)
            return pdf
    if CAIRO_PY.is_file():
        subprocess.run(
            [str(CAIRO_PY), "-c",
             "import sys,cairosvg;cairosvg.svg2pdf(url=sys.argv[1],write_to=sys.argv[2])",
             str(svg), str(pdf)], check=True)
        return pdf
    raise SystemExit("no SVG converter available; install librsvg2-bin or cairosvg")


def sidecars(pdf: Path) -> None:
    """The sidecar set preprint/figures/ already requires of every committed PDF.

    A figure entering that directory joins an existing contract: `<name>.sha256`
    of the PDF's own bytes, `<name>.txt` from pdftotext -layout, and `<name>.png`.
    The .txt is what RETRACTED_IN_FIGURES scans, so producing it is not
    bookkeeping -- it is what puts this figure's text under the phrase guard that
    exists because a withdrawn claim once survived inside an image.
    """
    # BARE HASH, matching what preprint/figures/ already uses. The first version
    # wrote "<hash>  <name>", a format invented beside the convention rather than
    # mirroring it, and the existing comparison reads the whole file.
    pdf.with_suffix(".sha256").write_text(f"{sha256(pdf)}\n", encoding="utf-8")
    txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                         capture_output=True, text=True, check=True).stdout
    pdf.with_suffix(".txt").write_text(txt, encoding="utf-8")
    import shutil
    if shutil.which("pdftoppm"):
        subprocess.run(["pdftoppm", "-png", "-r", "150", "-singlefile",
                        str(pdf), str(pdf.with_suffix(""))], check=True)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(__doc__)
    svg = Path(argv[1])
    pdf = render(svg)
    sidecars(pdf)
    # THE PAIR, NOT THE SOURCE ALONE. The docstring above argues that a record
    # comparing the SVG to a record of itself can be satisfied by regenerating
    # the record, leaving the PDF stale -- and then the first version of this
    # line wrote exactly that record. Codex raised it on pull request #23.
    #
    # Recording the digest of the PDF this conversion produced states the pair
    # the guard actually cares about. BE CLEAR ABOUT WHAT THIS DOES AND DOES NOT
    # BUY: both lines are still hand-writable, so it does not make forgery
    # impossible, it makes it require two deliberate falsehoods instead of one.
    # What actually binds source to artifact is the text check in
    # calibration/tests/test_figure_staleness.py, which reads the SVG's own text
    # runs out of the committed .txt and cannot be satisfied by editing a hash.
    #
    # `split()[0]` still yields the source digest, so every existing reader of
    # this file is unaffected.
    digest = sha256(svg)
    Path(str(svg) + ".sha256").write_text(
        f"{digest}\npdf={sha256(pdf)}\n", encoding="utf-8")
    print(f"rendered {pdf} from {svg}\n  source sha256 {digest}\n"
          f"  sidecars: {pdf.with_suffix('.sha256').name}, {pdf.with_suffix('.txt').name}, "
          f"{pdf.with_suffix('.png').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
