"""Every `file:line` citation in docs/guides/ resolves to what it claims.

WHY THIS EXISTS. A guide whose premise is that each claim is checkable against a
line of source is, without this, a tour with footnotes. Nothing in this
repository verified a file:line citation before it: `manuscript_claims_tests.jl`
does a regex SEARCH (find every line matching a pattern), and
`test_claims_ledger.py` resolves a git SHA. Neither takes a pre-named citation
and asks whether it still points at what it says.

WHAT MAKES A CITATION RESOLVE RATHER THAN MERELY POINT. File-exists and
line-exists are nearly free to satisfy -- a citation to any line of a file that
happens to be long enough passes both. The fragment is what carries the claim,
so every citation states one and the fragment must be on that line.

AND THE FAILURE MESSAGE HAS TO SPLIT TWO CASES, or the repair becomes reflexive.
When an edit shifts a line the guard fails, and the cheap fix is to bump the
number. That is right when the code is unchanged and moved, and wrong when the
fragment moved because the code changed -- those need different responses. So on
failure the guard searches the file for the fragment and says which case it is.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
GUIDES = REPO / "docs" / "guides"

# A row of a "Checks against the source" table:
#   | claim | `path:line` | `fragment` |
# The fragment is taken verbatim between the last pair of backticks, so it may
# contain pipes escaped as \| -- which the Robin condition's |_{r=R} needs.
ROW = re.compile(
    r"^[ \t]*\|(?P<claim>[^|]*)\|\s*`(?P<path>[^`:]+):(?P<line>\d+)`\s*\|\s*`(?P<frag>.+?)`\s*\|\s*$",
    re.M)


# A LOOSE PATTERN FOR THE SAME THING, USED ONLY TO COUNT. `ROW` is strict on
# purpose, and a strict parser silently drops what it cannot read -- so the count
# it produces is a SUBSET while the assertion below it would be about the SET.
# That is the defect this repository has now had caught three times, and putting
# it in the guard written to check citations would be the joke completing itself.
# So anything shaped like a citation row is counted, and the strict parser must
# account for every one.
# LEADING WHITESPACE IS ALLOWED IN BOTH, and it has to be both. A row indented
# two spaces -- still a valid Markdown table row -- matched NEITHER pattern, so
# `cites` and `loose` shrank together and the count-equality check still held.
# The safety net added so an unparsed row could not hide was itself blind to the
# one way a row hides. A citation smuggled in indented was never resolved and
# never reported.
LOOSE_ROW = re.compile(r"^[ \t]*\|.*`[^`]+:\d+`.*\|", re.M)


def guide_files() -> list[Path]:
    return sorted(GUIDES.glob("*.md")) if GUIDES.is_dir() else []


def citations(text: str) -> list[dict]:
    out = []
    for m in ROW.finditer(text):
        out.append({
            "claim": m.group("claim").strip(),
            "path": m.group("path").strip(),
            "line": int(m.group("line")),
            # `\|` is markdown's escape for a pipe inside a table cell; the
            # source line contains a bare pipe.
            "frag": m.group("frag").replace(r"\|", "|"),
        })
    return out


def resolve(cite: dict) -> tuple[bool, str]:
    """(ok, message). The message names WHICH repair a failure needs."""
    target = REPO / cite["path"]
    if not target.is_file():
        return False, f"no such file: {cite['path']}"
    lines = target.read_text(encoding="utf-8", errors="replace").split("\n")
    if cite["line"] > len(lines):
        return False, (f"{cite['path']} has {len(lines)} lines; cited "
                       f"{cite['line']}")
    if cite["frag"] in lines[cite["line"] - 1]:
        return True, ""

    # THE SPLIT. Where else does the fragment live?
    elsewhere = [i for i, l in enumerate(lines, 1) if cite["frag"] in l]
    if elsewhere:
        return False, (
            f"THE CODE MOVED, NOT THE CLAIM. {cite['path']}:{cite['line']} no "
            f"longer contains this fragment, but line(s) {elsewhere} do. "
            f"Renumber the citation to :{elsewhere[0]}.")
    return False, (
        f"THE CODE CHANGED. {cite['path']} no longer contains this fragment "
        f"anywhere, so the claim it supports may no longer be true. RE-READ the "
        f"source before renumbering -- bumping the line number here would be "
        f"the reflexive repair this message exists to prevent. "
        f"Claim: {cite['claim']!r}")


@pytest.mark.skipif(not guide_files(), reason="no guides committed yet")
@pytest.mark.parametrize("guide", guide_files(), ids=lambda p: p.name)
def test_every_guide_citation_resolves(guide):
    text = guide.read_text(encoding="utf-8")
    cites = citations(text)
    assert cites, (
        f"{guide.name} has no parseable citations. A guide in docs/guides/ that "
        "cites nothing is the shape this check exists to prevent; if it "
        "genuinely makes no claim about source, say so in the file.")

    # EVERY citation-shaped row, not just the ones the strict parser liked. A
    # row that fails ROW would otherwise be checked by nothing while this test
    # reported green over the rows it happened to read.
    loose = len(LOOSE_ROW.findall(text))
    assert len(cites) == loose, (
        f"{guide.name}: {loose} rows look like citations and only {len(cites)} "
        f"parsed, so {loose - len(cites)} are being checked by nothing. Fix the "
        "row's format or widen ROW -- do not leave it silently unread.")

    failures = []
    for c in cites:
        ok, msg = resolve(c)
        if not ok:
            failures.append(msg)
    assert not failures, "\n\n".join(failures)


def test_the_citation_check_detects_a_shifted_line():
    """CONTROL: shift a real citation and require the guard to catch it.

    SCOPE: one citation drawn from a committed guide, shifted by ONE and by TEN.

    BOTH SHIFTS, AND ONE IS NOT ENOUGH. Off-by-one is exactly the case a fragment
    match can survive by accident: `w_plus` and `w_minus` sit on consecutive
    lines of face_weights(), as do `rtol` and `atol` in the ode() call, so a
    control that only shifted by one could pass while the guard was blind to
    small drift. Ten crosses out of any such neighbourhood.

    The known-bad is DERIVED FROM THE ARTIFACT -- a real citation from a real
    guide -- rather than hand-written, because a synthetic citation would test
    the regex against an idea of the format instead of the format in use.
    """
    guides = guide_files()
    assert guides, "no guide to draw a control from"
    cites = citations(guides[0].read_text(encoding="utf-8"))
    assert cites, "no citations to draw a control from"

    real = cites[0]
    ok, msg = resolve(real)
    assert ok, f"baseline citation does not resolve, so no control verdict: {msg}"

    for delta in (1, 10):
        shifted = dict(real, line=real["line"] + delta)
        ok, msg = resolve(shifted)
        assert not ok, (
            f"shifting {real['path']}:{real['line']} by {delta} was NOT caught, "
            "so this guard cannot detect a stale line number")
        assert "THE CODE MOVED" in msg, (
            f"shift by {delta} was caught but misdiagnosed; the fragment still "
            f"exists in the file, so the message must say the code moved: {msg}")


def test_the_two_failure_cases_are_distinguished():
    """CONTROL: a fragment that exists nowhere must not be reported as moved.

    SCOPE: one synthetic fragment against a real file. This half MUST be
    synthetic -- the case is 'no line contains this', which cannot be drawn from
    a citation that resolves.
    """
    guides = guide_files()
    assert guides, "no guide to draw a control from"
    real = citations(guides[0].read_text(encoding="utf-8"))[0]

    gone = dict(real, frag="Notarealfragment_zzz")
    ok, msg = resolve(gone)
    assert not ok
    assert "THE CODE CHANGED" in msg, (
        "a fragment absent from the whole file was not diagnosed as a changed "
        f"claim, so the two repairs are not being distinguished: {msg}")


# ------------------------------------------- the rendered artifact, not just the source
#
# THE SOURCE-OF-RECORD IS THE .md AND A READER RECEIVES THE .pdf. Everything above
# reads `docs/guides/*.md`. On 2026-08-31 the committed PDF cited
# biofilms_potts.jl:1444, biofilms_..._fitness.tex:660 and :912 -- three line
# numbers the .md had already been corrected away from, because edits elsewhere had
# shifted them. Every check above passed. That is the figure-staleness defect
# exactly: prose right, artifact wrong, and the guard reading only the prose.
#
# It is also the rendered-versus-source gap this session recorded one commit
# earlier while writing about it, which is why the binding is asserted here rather
# than left to the discipline of remembering to re-render.
#
# SCOPE, NARROWER THAN "THE PDF IS CURRENT": line numbers only. A prose edit to the
# .md that never touches a citation leaves the PDF stale and this check green. The
# .tex is a hand-maintained rendering, not generated, so nothing binds its prose to
# the .md's -- that is a real and open gap, named here rather than implied away.
def guide_pdfs() -> list[Path]:
    return sorted(GUIDES.glob("*.pdf")) if GUIDES.is_dir() else []


def pdf_text(pdf: Path) -> str:
    import subprocess
    return subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                          capture_output=True, text=True).stdout


@pytest.mark.skipif(not guide_pdfs(), reason="no rendered guide committed")
@pytest.mark.parametrize("pdf", guide_pdfs(), ids=lambda p: p.name)
def test_the_rendered_guide_cites_what_the_source_cites(pdf):
    import shutil
    md = pdf.with_suffix(".md")
    assert md.is_file(), f"{pdf.name} has no .md source of record"
    if not shutil.which("pdftotext"):
        pytest.skip("pdftotext unavailable; the rendered-artifact check needs Poppler")

    cites = citations(md.read_text(encoding="utf-8"))
    assert cites, f"{md.name} yields no citations to compare against"
    text = " ".join(pdf_text(pdf).split())

    # MATCH WHAT THE EXTRACTION PRODUCES, NOT WHAT THE SOURCE WROTE. The template
    # breaks long paths at underscores to fit the table column, and `pdftotext
    # -layout` puts the two halves in DIFFERENT COLUMNS -- `biofilms_` on one side,
    # `radiodialysis.R:223` on the other -- so no amount of whitespace-joining
    # rejoins them. Checked against the real text layer rather than assumed: the
    # segment after the final underscore survives intact for all seven paths, and
    # full-path matching reported all fourteen citations missing from a PDF that
    # in fact contained every one. That is the extraction-shape lesson the figure
    # sidecars taught, arriving in a second artifact.
    def _tail(path: str) -> str:
        return Path(path).name.split("_")[-1]

    missing = [f"{c['path']}:{c['line']}" for c in cites
               if f"{_tail(c['path'])}:{c['line']}" not in text]
    # ...and no line number the source has moved away from may survive in it. This
    # half is what catches a stale render: the first is satisfied by a PDF that
    # merely contains the right numbers somewhere among the wrong ones.
    current = {str(c["line"]) for c in cites}
    import re as _re
    stale = sorted({m for m in _re.findall(r"\.(?:jl|R|py|tex|md):(\d+)", text)
                    if m not in current})

    assert not missing and not stale, (
        f"{pdf.name} is not the render of {md.name}. Missing from the artifact: "
        f"{missing}. Line numbers in the artifact the source no longer cites: "
        f"{stale}. Re-render the guide -- the .md is the source of record, and "
        "editing it without re-rendering leaves a reader holding wrong citations.")


# ------------------------------------------------- the .md moved without a render
#
# THE GAP THIS CLOSES, AND THE ONE IT LEAVES OPEN. Everything above binds line
# NUMBERS. A prose edit to the .md that touches no citation left the .tex and the
# PDF holding the old wording with every check green -- which happened here on
# 2026-08-31 and was repaired by hand-syncing three files twice in one session.
#
# tools/render_guide.py writes `<name>.md.sha256` ONLY as a side effect of a
# successful render that verified every citation in the artifact. So the record
# cannot be satisfied by editing it: a hand-written hash passes this check and
# leaves the citation check to fail on the next real run. Same inversion as
# tools/render_figure_svg.py.
#
# IT DOES NOT VERIFY THAT THE .tex SAYS WHAT THE .md SAYS. Measured before
# choosing: the .md's prose paragraphs match the rendered text layer 38 of 48 at a
# six-word prefix, and the ten failures are MARKUP differences rather than
# divergence -- backticked `r` renders as math, a literal `##` in prose renders as
# emphasis. The .tex is a re-authoring, not a render, so the two differ by design
# and a prose binding would be a fifth false failures. That is the wrong
# instrument for a hand-authored .tex, not a threshold to tune. A faithful
# re-sync remains a human obligation and an open gap.
def guide_records() -> list[Path]:
    return sorted(GUIDES.glob("*.md.sha256")) if GUIDES.is_dir() else []


@pytest.mark.skipif(not guide_records(), reason="no rendered guide recorded")
@pytest.mark.parametrize("record", guide_records(), ids=lambda p: p.name)
def test_the_guide_source_matches_its_recorded_render(record):
    import hashlib
    md = record.with_suffix("")           # strip .sha256, leaving <name>.md
    assert md.is_file(), f"{record.name} records a hash for a missing {md.name}"
    recorded = record.read_text(encoding="utf-8").split()[0]
    actual = hashlib.sha256(md.read_bytes()).hexdigest()
    assert recorded == actual, (
        f"{md.name} has changed since the guide was last rendered, so the "
        f"committed .tex and .pdf are older than the source of record. Re-sync "
        f"the .tex by hand, then run `python3 tools/render_guide.py {md}`. Do NOT "
        "edit the .sha256: that makes this pass and leaves the citation check to "
        "fail on the next real render. Note this catches only that the source "
        "MOVED -- nothing verifies the re-sync was faithful.")
