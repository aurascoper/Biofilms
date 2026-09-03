#!/usr/bin/env python3
"""Render an SVG figure to PDF and record THE SOURCE/OUTPUT PAIR IT SAW.

WHY THE HASH IS WRITTEN HERE AND NOWHERE ELSE. An SVG source-of-truth plus a
committed derived PDF is the generator-and-artifact pair that already failed in
this repository: the generator was corrected two weeks before the committed
images were, and nothing could open a figure to notice.

A plain `svg.sha256` sidecar does NOT fix that. Compared against the SVG it binds
the SVG to a record of itself -- edit the SVG, the hash mismatches, and the
obvious repair is to regenerate the hash. That is green again with the PDF still
stale, so the check meant to catch staleness has been satisfied without catching
it.

The record therefore contains source, output, and a framed joint digest. The
component hashes diagnose which file moved; the joint digest makes their exact
combination the unit the guard accepts. The file remains writable, so this turns
an accidental "update the source hash" repair into an explicit forged pair rather
than pretending hashes prove provenance. The PDF text is independently checked
through its own extraction paths.

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


def pair_sha256(svg: Path, pdf: Path) -> str:
    """Digest the ordered source/output pair with explicit framing."""
    digest = hashlib.sha256(b"biofilms-svg-pdf-pair-v1\0")
    for path in (svg, pdf):
        body = path.read_bytes()
        digest.update(len(body).to_bytes(8, "big"))
        digest.update(body)
    return digest.hexdigest()


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
    # THE PAIR, NOT TWO INDEPENDENT FACTS. The docstring above argues that a
    # record comparing the SVG to itself can be satisfied by regenerating the
    # record, leaving the PDF stale. The second version put independent source
    # and PDF hashes in the same file: editing the source hash still left a valid
    # old PDF hash, so the stale combination remained admissible. Codex raised
    # both versions on pull request #23.
    #
    # A framed joint digest makes the exact combination the renderer saw a
    # first-class value. Component hashes remain for diagnosis, but neither can
    # substitute for the pair. This does not make a text file unforgeable; it
    # turns the ordinary "update the source hash" repair into an explicit false
    # claim about the source/output pair. Exact text binding is checked through
    # the PDF itself in test_figure_staleness.py.
    source_digest = sha256(svg)
    pdf_digest = sha256(pdf)
    Path(str(svg) + ".sha256").write_text(
        f"source={source_digest}\npdf={pdf_digest}\npair={pair_sha256(svg, pdf)}\n",
        encoding="utf-8")
    print(f"rendered {pdf} from {svg}\n  source sha256 {source_digest}\n"
          f"  sidecars: {pdf.with_suffix('.sha256').name}, {pdf.with_suffix('.txt').name}, "
          f"{pdf.with_suffix('.png').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
