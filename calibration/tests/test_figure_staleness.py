"""A rendered figure must be paired with the SVG currently in the tree.

SCOPE: the SVG/PDF/sha256 triples under preprint/figures/, and nothing else.

This exists because the generator-and-artifact pair has already failed here once:
a figure generator was corrected two weeks before its committed images were, and
no check could open a figure to notice. FIG-01 through FIG-07 record it.

THE RECORD IS WRITTEN BY tools/render_figure_svg.py and names source, output and
their framed joint digest. Independent hashes describe two files; the pair digest
describes the one combination the renderer saw. Exact text binding then enters
through the PDF itself, including the one token the layout extraction clips.
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


def pair_sha256(svg: Path, pdf: Path) -> str:
    """Digest the ordered source/output pair with explicit framing."""
    digest = hashlib.sha256(b"biofilms-svg-pdf-pair-v1\0")
    for path in (svg, pdf):
        body = path.read_bytes()
        digest.update(len(body).to_bytes(8, "big"))
        digest.update(body)
    return digest.hexdigest()


def render_record_mismatches(svg: Path, pdf: Path, record_text: str) -> dict:
    """Return every field where a render record differs from its exact pair."""
    lines = record_text.splitlines()
    if len(lines) != 3:
        return {"format": f"expected three fields, got {len(lines)}"}
    try:
        fields = dict(line.split("=", 1) for line in lines)
    except ValueError:
        return {"format": "every line must be key=value"}
    if set(fields) != {"source", "pdf", "pair"}:
        return {"format": f"expected source/pdf/pair, got {sorted(fields)}"}
    expected = {
        "source": hashlib.sha256(svg.read_bytes()).hexdigest(),
        "pdf": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "pair": pair_sha256(svg, pdf),
    }
    return {key: (fields[key], value) for key, value in expected.items()
            if fields[key] != value}


@pytest.mark.parametrize("svg", svg_triples(), ids=lambda p: p.name)
def test_the_rendered_pdf_matches_the_svg_in_the_tree(svg):
    record = Path(str(svg) + ".sha256")
    pdf = svg.with_suffix(".pdf")
    assert record.is_file(), f"{svg.name} has no render record; run tools/render_figure_svg.py"
    assert pdf.is_file(), f"{svg.name} has no rendered PDF"
    mismatches = render_record_mismatches(
        svg, pdf, record.read_text(encoding="utf-8"))
    assert not mismatches, (
        f"{svg.name} and {pdf.name} are not the pair the renderer recorded "
        f"({mismatches}). Re-run "
        f"tools/render_figure_svg.py {svg.relative_to(REPO)} -- do NOT edit the "
        ".sha256 file: independent component hashes do not establish a pair.")


def test_the_render_record_rejects_a_stale_source_output_combination(tmp_path):
    """CONTROL: independently valid component hashes are not a valid pair.

    The source and PDF come from the production artifact. The source is mutated,
    its component hash is updated as the reported escape did, and the unchanged
    PDF retains its valid component hash. This calls the production comparison;
    a helper-only digest assertion would not prove that path consumes the pair.
    """
    original_svg = FIGDIR / "phase2_diffusion_cell.svg"
    original_pdf = original_svg.with_suffix(".pdf")
    svg = tmp_path / original_svg.name
    pdf = tmp_path / original_pdf.name
    source = original_svg.read_bytes()
    victim = b"self-diffusivity"
    assert victim in source, "artifact-derived mutation victim is absent"
    svg.write_bytes(source.replace(victim, b"self-diffusivNOTRENDERED", 1))
    pdf.write_bytes(original_pdf.read_bytes())
    stale_pair = pair_sha256(original_svg, original_pdf)
    forged_components = (
        f"source={hashlib.sha256(svg.read_bytes()).hexdigest()}\n"
        f"pdf={hashlib.sha256(pdf.read_bytes()).hexdigest()}\n"
        f"pair={stale_pair}\n")
    mismatches = render_record_mismatches(svg, pdf, forged_components)
    assert set(mismatches) == {"pair"}, mismatches


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
# WHAT THE COMPONENT HASHES CANNOT DO, AND THE FILE THEY GUARD SAID SO FIRST.
# Edit the SVG, update its component hash, and leave the old PDF with its still
# valid component hash: two independently true statements admit a stale pair.
# The framed pair field above makes that exact combination fail unless somebody
# deliberately rewrites the relation too. Raised by Codex on pull request #23.
#
# A HASH IS NOT PROVENANCE. The record is writable and a deliberate forgery is
# possible; its job is to prevent independent component updates from posing as a
# render receipt, not to authenticate an author.
#
# So the binding is content, not digest: every text run in the SVG must appear in
# the committed `.txt`, which the sidecar test independently pins to the PDF's
# actual bytes. The one layout-clipped suffix is verified through pdftohtml's
# hidden-text extraction from the PDF itself.
#
# SCOPE, STATED BECAUSE IT IS NARROWER THAN "THE PDF IS FRESH": this catches
# TEXT drift only. An SVG edit that moves a rectangle and touches no string
# passes, and no check in this repository would catch it. That is the honest
# reach -- and text is the failure mode the figure guards exist for, since
# FIG-01 through FIG-10 are every one of them a claim made in words inside an
# image.
#
# MATCHING COVERS THE WHOLE RUN, NOT ITS OPENING WORDS, AND THAT CHANGE IS THE
# POINT. The first version compared the first five words. Everything after word
# five could then diverge between SVG and PDF while this passed -- and the worst
# instance is the one that matters: changing `43` to `99` in the phase2 caption
# escaped entirely, and 43 is THE NUMBER THE FIGURE WAS CORRECTED FOR, sitting at
# word 24 of a 27-word run. Raised by Codex on pull request #23.
#
# THE RULE WAS MEASURED, NOT REASONED. Four candidates were run against all 39
# runs of the committed figure: full-run whitespace-collapsed (1 false failure),
# full-run whitespace-removed (1), first-five-words (0), head-5 + tail-5 (1).
# Every full-coverage candidate failed on the SAME run -- the one containing 43 --
# and the cause is that 26 of its 27 words match exactly while the final token is
# clipped at the column boundary: `self-diffusivity).` extracts as
# `self-diffusiv`. Hence: all words but the last must appear contiguously, and the
# last may be clipped.
#
# THAT RULE RESTS ON A PROPERTY OF THIS LAYOUT AND A FUTURE FAILURE SHOULD BE
# DIAGNOSED ACCORDINGLY. Clipping happens at the column boundary, which is why it
# lands on the final token. Widen a column, change the page, re-flow the figure,
# and a MIDDLE token could clip -- at which point this reports a false failure
# that looks exactly like a real one, and the obvious reading would be that
# somebody edited the SVG without re-rendering. A failure whose missing word is
# mid-run is more likely layout drift than tampering; check the render before
# concluding staleness.


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


CLIPPED_FINAL_TOKENS = {
    # LOCATION-SCOPED, NOT VALUE-KEYED OR GLOBALLY INFERRED. Exactly one run in
    # the measured figure crosses the right column boundary. The layout sidecar
    # ends at this prefix, while `pdftohtml -hidden` reads the complete token
    # from the PDF and is checked below. A changed suffix creates a different run
    # identity and does not inherit this exception.
    "phase2_diffusion_cell.svg": {
        "The film-scale effective diffusivity measured here is not the reactor-scale "
        "dispersion coefficient that currently carries the same name in radial_parms "
        "(1×10⁻³ cm² s⁻¹, about 43× water’s self-diffusivity).": "self-diffusiv",
    },
}


def unrendered_runs(svg_source: str, txt: str,
                    allowed_clips: dict[str, str] | None = None) -> list[str]:
    """Runs absent from the extraction, with exact location-earned clips only.

    Shared by the check and its controls, so a control cannot pass against a
    regressed check.
    """
    allowed_clips = allowed_clips or {}
    flat = " ".join(txt.split())
    out = []
    for s in svg_text_runs(svg_source):
        # THE CLIPPING ALLOWANCE IS EARNED PER RUN, NOT GRANTED TO ALL OF THEM.
        # The first version dropped the final token of EVERY run, because one run
        # in 39 has its last word clipped at the column boundary. That left the
        # last token of all 39 unchecked -- and one of them is a citation year:
        # `volume bromide did not (Hay et al. 2011,`. Changing 2011 to 2019 left
        # the guard green, forging a citation inside the figure whose entire
        # history is withdrawn claims surviving in images.
        #
        # So: try the whole run first. Only if that fails may the final token be
        # dropped, and only then. A run that does not clip is checked in full.
        if s in flat:
            continue
        w = s.split()
        allowed_prefix = allowed_clips.get(s)
        if allowed_prefix and len(w) > 1:
            head = " ".join(w[:-1])
            i = flat.find(head)
            if i >= 0:
                # The allowance belongs to this exact run at this exact figure,
                # and the observed fragment must still equal its recorded prefix.
                after = flat[i + len(head):].split()
                tail = w[-1]
                frag = after[0] if after else ""
                if (frag == allowed_prefix and tail.startswith(allowed_prefix)
                        and allowed_prefix != tail):
                    continue
        out.append(s)
    return out


def _run_extent(svg_source: str, run: str):
    """(start, end) of the <text> element whose normalised content is `run`.

    Returns the extent of THAT element, so a control can assert its mutation
    landed inside the run it selected rather than inside some run.
    """
    for m in re.finditer(r"<text[^>]*>(.*?)</text>", svg_source, re.S):
        inner = html.unescape(re.sub(r"<[^>]+>", "", m.group(1)))
        if " ".join(inner.split()) == run:
            return m.start(), m.end()
    return None


@pytest.mark.parametrize("svg", svg_triples(), ids=lambda p: p.name)
def test_the_svg_text_reaches_the_committed_pdf_text(svg):
    txt_path = svg.with_suffix(".txt")
    assert txt_path.is_file(), f"{svg.name} has no .txt sidecar; run tools/render_figure_svg.py"
    missing = unrendered_runs(
        svg.read_text(encoding="utf-8"),
        txt_path.read_text(encoding="utf-8"),
        CLIPPED_FINAL_TOKENS.get(svg.name))
    assert not missing, (
        f"{svg.name} carries text that is not in the rendered artifact, so the "
        f"committed PDF is older than the SVG: {missing[:3]}. Re-run "
        f"tools/render_figure_svg.py {svg.name}. Editing the .sha256 will make the "
        "hash checks pass and will not make this one. If the missing word is "
        "MID-run rather than at the end, suspect layout drift before staleness -- "
        "see the note on the clipped final token above.")


@pytest.mark.parametrize("svg", svg_triples(), ids=lambda p: p.name)
def test_every_declared_clip_is_complete_in_the_pdf(svg):
    """The layout exception is checked through a second, complete PDF path.

    `pdftotext -layout` clips at the visible column boundary. Poppler's hidden
    XML extraction retains the complete text object, so this verifies the suffix
    the layout sidecar cannot observe instead of granting it globally.
    """
    clips = CLIPPED_FINAL_TOKENS.get(svg.name, {})
    if not clips:
        return
    import shutil
    import subprocess
    if not shutil.which("pdftohtml"):
        pytest.skip("pdftohtml unavailable; the clipped suffix needs its PDF path")
    got = subprocess.run(
        ["pdftohtml", "-xml", "-hidden", "-stdout", str(svg.with_suffix(".pdf"))],
        capture_output=True, text=True)
    assert got.returncode == 0, f"pdftohtml failed for {svg.name}: {got.stderr}"
    hidden = html.unescape(re.sub(r"<[^>]+>", " ", got.stdout))
    for run in clips:
        suffix = " ".join(run.split()[-3:])
        assert suffix in hidden, (
            f"{svg.name}'s declared clipped suffix {suffix!r} is not complete in "
            "the PDF; the location-scoped layout exception is unsupported")


def _mutate_inside(source: str, run: str, old: str, new: str):
    """Replace `old` with `new` INSIDE the element whose content is `run`.

    Two properties the previous control had neither of. It did
    `source.replace(victim.split()[0], "Notarealword", 1)` -- the first occurrence
    of that word ANYWHERE in the file, which happened to fall inside the victim
    only because the longest run contains the file's first "The". Correct by
    coincidence. And the replacement operated on the NORMALISED run text, while
    the source stores entities (`&#215;`, `&#8217;`), so a mutation of any run
    containing one silently did nothing at all.
    """
    span = _run_extent(source, run)
    assert span, "victim run not locatable in the source"
    lo, hi = span
    seg = source[lo:hi]
    assert old in seg, f"{old!r} is not inside the selected run"
    return source[:lo] + seg.replace(old, new, 1) + source[hi:], (lo, hi)


def test_the_source_to_artifact_binding_detects_an_edited_svg():
    """CONTROL: an SVG edited without re-rendering must be caught, anywhere in a run.

    SCOPE: the phase2 figure's longest run for head/middle/tail, plus ONE
    DIFFERENT run so the rule is not tested only on the instance it was derived
    from. The known-bad is drawn from the artifact path rather than written by
    hand.

    AND THE MUTATION IS ASSERTED TO HAVE APPLIED, AND WHERE. The first attempt at
    reproducing the Codex finding was a no-op -- it replaced normalised text that
    the entity-encoded source does not contain -- and the guard's truthful
    "nothing missing" over an unmodified file read exactly like confirmation.
    That is a distinct failure from a vacuous assertion: the assertion was fine,
    THE INPUT NEVER CHANGED. So every mutation below asserts it applied and that
    it landed inside the run selected, by identity rather than by membership in
    the set of runs -- an offset off by enough lands in a neighbour and would
    otherwise pass.
    """
    svg = FIGDIR / "phase2_diffusion_cell.svg"
    source = svg.read_text(encoding="utf-8")
    txt = svg.with_suffix(".txt").read_text(encoding="utf-8")
    allowed_clips = CLIPPED_FINAL_TOKENS[svg.name]

    def missing(candidate):
        return unrendered_runs(candidate, txt, allowed_clips)

    assert not missing(source), "baseline is not clean; no control verdict"
    runs = svg_text_runs(source)
    assert runs, "no text runs found -- the extractor is broken, not the figure"

    # TRAINING SET: the run the clipped-final-token rule was derived from. Head,
    # middle and tail, because the defect was that only the head was ever read.
    victim = max(runs, key=len)
    lo0, hi0 = _run_extent(source, victim)
    for label, old, new in (("head", "The film-scale", "The notreal"),
                            ("middle", "reactor-scale", "notreal-scale"),
                            ("tail", "about 43&#215;", "about 99&#215;")):
        mutated, (lo, hi) = _mutate_inside(source, victim, old, new)
        assert mutated != source, f"{label}: mutation did not apply"
        # AN INVARIANT RESTATED, NOT A CHECK, AND SAYING SO IS THE POINT.
        # `_mutate_inside` locates the element by content identity and replaces
        # only within that slice, so the mutation CANNOT land in a neighbouring
        # run while that function is correct. Weakening this line to any
        # tautology leaves the suite green -- measured. It is kept because it
        # documents the invariant the construction provides, and it would catch a
        # future rewrite of `_mutate_inside` that located by offset arithmetic
        # instead. A membership test ("falls inside SOME run") would be the weak
        # version and is deliberately not what this compares.
        assert (lo, hi) == (lo0, hi0), (
            f"{label}: mutation landed in a different element than the one "
            "selected -- _mutate_inside no longer bounds the replace to the "
            "located run")
        assert missing(mutated), (
            f"{label} edit was not caught -- text after the opening words can "
            "diverge from the artifact unnoticed")

    # THE REPORTED ESCAPE: preserve the fragment visible in pdftotext's clipped
    # layout output and alter only the unobservable suffix. The exception is
    # attached to the original run identity, so the changed run cannot inherit
    # it. This is the control that distinguishes location scope from prefix
    # inference.
    mutated, _ = _mutate_inside(
        source, victim, "self-diffusivity).", "self-diffusivNOTRENDERED")
    assert mutated != source, "clipped-suffix control: mutation did not apply"
    assert missing(mutated), (
        "an edit beyond the visible clipped prefix inherited the exception")

    # THE CASE THE FOUR CONTROLS BELOW STRUCTURALLY CANNOT REACH: the final token
    # of a run that does NOT clip. Three of them mutate the one run that clips and
    # the fourth mutates a first word, so none could see that the clipping
    # allowance was being granted to all 39 runs rather than to the one that earns
    # it. That hole let `2011,` become `2019,` in a citation with the suite green.
    #
    # 35 of 39 runs are present in the extraction in full. This mutates one of
    # them, and asserts the mutation applied -- a no-op here would report ESCAPES
    # for the wrong reason, which happened while writing this: the first candidate
    # token came from a phrase withdrawn from the figure months ago and no longer
    # in the source at all.
    nonclipping = [r for r in runs if r in " ".join(txt.split()) and len(r.split()) > 3]
    assert nonclipping, "no full-present run to draw the control from"
    tail = nonclipping[0].split()[-1]
    mutated, _ = _mutate_inside(source, nonclipping[0], " " + tail, " Notarealword")
    assert mutated != source, "final-token control: mutation did not apply"
    assert missing(mutated), (
        "the final token of a NON-clipping run was edited and not caught -- the "
        "clipping allowance is being granted to runs that do not clip")

    # TEST SET, and it is labelled as such the way tools/absence_gate.py labels
    # its own: a DIFFERENT run, whose final token is not clipped. Deriving the
    # rule from one run and only ever testing it there is the training-set
    # problem, and this is the item that shows the rule generalises.
    other = next(r for r in runs
                 if r != victim and len(r.split()) >= 4
                 and " ".join(r.split()) in " ".join(txt.split()))
    first_word = other.split()[0]
    mutated, _ = _mutate_inside(source, other, first_word, "Notarealword")
    assert mutated != source, "test-set mutation did not apply"
    assert missing(mutated), (
        "a mutation in a run OTHER than the one the rule was derived from was "
        "not caught, so the rule is overfit to that run")
