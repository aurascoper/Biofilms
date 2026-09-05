#!/usr/bin/env python3
"""Render a guide's LaTeX to PDF and record the .md hash THE RENDER SAW.

    python3 tools/render_guide.py docs/guides/calculus_in_this_code.md

WHY THIS EXISTS. docs/guides/<name>.md is the declared source of record;
<name>.tex is a hand-maintained rendering of it and <name>.pdf is built from the
.tex. Nothing bound the three. On 2026-08-31 the .md's prose was edited and the
.tex and .pdf kept the old wording while every check stayed green, because the
citation binding covers line numbers only -- which that test states about itself.
The .md was hand-synced to the .tex twice in one session.

WHAT THIS BINDS, AND IT IS NARROWER THAN "THE RENDER IS CURRENT". Measured before
choosing: matching the .md's prose paragraphs against the rendered text layer
gives 38 of 48 at a six-word prefix, and the ten failures are MARKUP differences
rather than divergence -- `r` in backticks renders as math, a literal `##` inside
prose renders as emphasis. THE .tex IS A RE-AUTHORING, NOT A RENDER, so the two
differ by design and a prose binding would be 20% false failures. That is not a
threshold to tune; it is the wrong instrument for a hand-authored .tex.

So this binds what is real:

  1. it renders the .tex with tectonic, so a PDF exists that came from the .tex;
  2. it verifies every `file:line` citation in the .md appears in that PDF's text
     layer, which is the check that has already caught three stale renders;
  3. and only then does it write `<name>.md.sha256`.

The record is written ONLY as a side effect of a successful render, the same
inversion tools/render_figure_svg.py uses: the sole way to make the staleness
check green is to run this, which produces the PDF. Editing the record by hand
makes the hash check pass and leaves the citation check to fail on the next run.

WHAT IT DOES NOT DO, STATED SO NOBODY READS MORE INTO A GREEN RUN: it does not
verify that the .tex says what the .md says in prose. A prose edit to the .md
that touches no citation will fail the hash check -- which is the point, it
forces a human to re-sync and re-render -- but nothing here checks that the
re-sync was faithful. That remains a human obligation and an open gap.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

ROW = re.compile(
    r"^\|(?P<claim>[^|]*)\|\s*`(?P<path>[^`:]+):(?P<line>\d+)`\s*\|\s*`(?P<frag>.+?)`\s*\|\s*$",
    re.M)


def citations(md_text: str) -> list[tuple[str, int]]:
    return [(m.group("path").strip(), int(m.group("line")))
            for m in ROW.finditer(md_text)]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(__doc__)
    md = Path(argv[1]).resolve()
    tex, pdf = md.with_suffix(".tex"), md.with_suffix(".pdf")
    for f in (md, tex):
        if not f.is_file():
            raise SystemExit(f"missing: {f}")
    if not shutil.which("tectonic"):
        raise SystemExit("tectonic not found; it is what renders the guide")

    subprocess.run(["tectonic", tex.name, "--outdir", str(tex.parent)],
                   cwd=tex.parent, check=True)

    if not shutil.which("pdftotext"):
        raise SystemExit("pdftotext not found; the citation check needs Poppler")
    text = " ".join(subprocess.run(
        ["pdftotext", "-layout", str(pdf), "-"],
        capture_output=True, text=True, check=True).stdout.split())

    # The tail after the final underscore is what survives the template's
    # line-breaking of long paths -- measured, see test_guide_citations.py.
    cites = citations(md.read_text(encoding="utf-8"))
    if not cites:
        raise SystemExit(f"{md.name} yields no citations; refusing to record a hash")
    missing = [f"{p}:{n}" for p, n in cites
               if f"{Path(p).name.split('_')[-1]}:{n}" not in text]
    if missing:
        raise SystemExit(
            f"REFUSING TO RECORD. {len(missing)} citation(s) in {md.name} are not "
            f"in the rendered PDF: {missing}. The .tex is not a render of this "
            ".md -- sync it before re-running.")

    digest = hashlib.sha256(md.read_bytes()).hexdigest()
    Path(str(md) + ".sha256").write_text(f"{digest}\n", encoding="utf-8")
    print(f"rendered {pdf.name} from {tex.name}\n"
          f"  {len(cites)} citations verified in the artifact\n"
          f"  source-of-record sha256 {digest}\n"
          f"  recorded in {md.name}.sha256")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
