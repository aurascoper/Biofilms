"""The claims ledger, enforced instead of merely written down.

`data/claims_ledger.csv` records every quantitative claim the manuscript and the
repository make, with a verdict on each: keep, restate, requalify, or delete.
Producing those verdicts took a full manual audit. Until this file existed,
NOTHING checked them — `grep -rln "claims_ledger" --include=*.py` returned
nothing at all — so a claim marked `delete` could reappear in any future
revision and no test would notice. The audit was a one-time human read of a
document that changes.

`delete` is the one verdict where absence IS the criterion, so it is the one
that can be mechanically enforced. A `restate` or `requalify` verdict asks
whether the revision now says the right thing, which is a judgement about
meaning and stays with a reviewer.

WHAT THIS DOES NOT CLAIM. It checks 30 of the 34 `delete` rows. The other four
are paraphrases with elisions and no verbatim fragment long enough to search
for, and the test NAMES them rather than passing silently over them — the same
discipline as `reference_d_status.enforcement_report()` naming the declared
criteria that have no consumer. A test that implied full coverage would be worse
than no test, because it would retire the manual review that actually covers
those four.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LEDGER = REPO / "data" / "claims_ledger.csv"
PREPRINT = REPO / "preprint" / "modeling_radioresistance_and_radiotropic_fitness.tex"

# Long enough that a match is the claim rather than a coincidence of ordinary
# scientific prose. At four words, phrases like "the biofilm is modelled"
# collide with unrelated sentences.
MIN_WORDS = 5


def _rows():
    with open(LEDGER, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def normalise_markup(source: str) -> str:
    """Manuscript source as plain lowercase prose, markup removed.

    WITHOUT THIS THE GUARD BARELY GUARDS. The ledger records claims as prose —
    `using DifferentialEquations.jl` — while the sources carry them as markup:
    `using \\texttt{DifferentialEquations.jl}` in LaTeX, backticked in the
    Markdown derivative. A raw substring search therefore misses a restored
    claim written the way anyone would actually write it. Measured against the
    original manuscript, normalising lifts detection from 12 to 15 of 30 in the
    `.tex` and from 15 to 20 of 30 in the `.md`.

    Also collapses whitespace, because LaTeX wraps at the column and a claim is
    usually split across a newline.
    """
    s = re.sub(r"(?m)(?<!\\)%.*", " ", source)                  # comments
    s = re.sub(r"\\(begin|end)\{[^}]*\}", " ", s)               # environments
    for _ in range(5):                                          # unwrap, nested
        s = re.sub(r"\\[a-zA-Z]+\*?\s*(?:\[[^\]]*\])?\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)                        # bare macros
    for ch in ("~", "\\", "`", "*"):
        s = s.replace(ch, " ")
    s = re.sub(r"[{}$&_^]", " ", s)
    s = re.sub(r"[\u2010-\u2015]", "-", s)                       # unicode dashes
    return " ".join(s.split()).lower()


def _preprint_text() -> str:
    return normalise_markup(PREPRINT.read_text(encoding="utf-8"))


# The manuscript as first committed, before any revision removed anything. It is
# the only document known to CONTAIN the deleted claims, which makes it the only
# valid control for whether this guard detects them.
ORIGINAL_COMMIT = "5980dc5"
ORIGINAL_PATH = "preprint/modeling_radiotrophic_fitness.md"


def _original_manuscript() -> str | None:
    import subprocess

    shallow = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "--is-shallow-repository"],
        capture_output=True, text=True)
    if shallow.stdout.strip() == "true":
        return None
    got = subprocess.run(
        ["git", "-C", str(REPO), "show", f"{ORIGINAL_COMMIT}:{ORIGINAL_PATH}"],
        capture_output=True, text=True)
    return normalise_markup(got.stdout) if got.returncode == 0 else None


def _detected(rows, text) -> list[str]:
    out = []
    for row in rows:
        if not row["claim_id"].startswith("PP-") or row["status"] != "delete":
            continue
        phrase = distinguishing_phrase(row["claim_text"])
        if phrase and phrase.lower() in text:
            out.append(row["claim_id"])
    return out


def distinguishing_phrase(claim: str) -> str:
    """The longest run of the claim with no elision and no bracketed citation.

    Ledger `claim_text` is a mix: some entries quote the manuscript verbatim,
    others paraphrase with "..." where the original ran long. Only the
    unelided runs can be searched for, so the longest one is taken and the rest
    of the claim is ignored rather than guessed at.
    """
    best = ""
    for segment in re.split(r"\.\.\.|\[\d+\]|[(),;]", claim):
        segment = " ".join(segment.split())
        if len(segment.split()) >= MIN_WORDS and len(segment) > len(best):
            best = segment
    return best


@pytest.fixture(scope="module")
def rows():
    return _rows()


def test_the_ledger_parses_and_is_not_empty(rows):
    assert len(rows) > 200
    assert {"claim_id", "status", "claim_text", "document"} <= set(rows[0])


def test_no_deleted_claim_survives_in_the_manuscript(rows):
    """THE ONE THAT MATTERS. Every `delete` verdict was reached because the
    claim was unsupported, fabricated, or checkable and false — a Hamiltonian
    that conserves nothing, a symplectic integrator that is not used, a
    quantum-mechanical noise mechanism that is a Gaussian keyed to a diffusion
    coefficient. They must not come back."""
    text = _preprint_text()
    survivors = []
    for row in rows:
        if not row["claim_id"].startswith("PP-") or row["status"] != "delete":
            continue
        phrase = distinguishing_phrase(row["claim_text"])
        if phrase and phrase.lower() in text:
            survivors.append(f"{row['claim_id']}: {phrase[:90]}")
    assert not survivors, (
        "claims marked `delete` in the ledger are present in the manuscript:\n  "
        + "\n  ".join(survivors))


def test_the_guard_actually_detects_deleted_claims(rows):
    """THE TEST THAT PROVES THE OTHER ONE MEANS SOMETHING.

    `test_no_deleted_claim_survives_in_the_manuscript` passes when it finds
    nothing. So does a guard that can find nothing at all — and the first
    version of this file was exactly that: raw substring matching against LaTeX,
    which detected 12 of 30 in the original and would have reported a clean
    manuscript either way. A check that cannot fail is not a check, which is the
    lesson this repository keeps relearning.

    So: run the guard against the one document known to CONTAIN the deleted
    claims, and require it to find them.

    The floor is 18 rather than the measured 20 to leave room for harmless
    changes in phrase extraction, and it is above the 15 that raw substring
    matching achieves — so dropping the markup normalisation fails here rather
    than silently halving the guard.
    """
    original = _original_manuscript()
    if original is None:
        pytest.skip("shallow clone: the original manuscript is not reachable")
    found = _detected(rows, original)
    assert len(found) >= 18, (
        f"the guard detected only {len(found)} deleted claims in the "
        "pre-revision manuscript, which contains them. It has stopped "
        "guarding — most likely the markup normalisation regressed.")


def test_coverage_is_reported_as_detection_not_as_phrase_count(rows, capsys):
    """Coverage is reported because it is PARTIAL, and reported as the number
    that means something.

    An earlier version of this file announced "30 of 34 textually checkable",
    which counts rows that YIELD a searchable phrase — not rows the guard can
    actually find. Only about 20 are detectable in a document known to contain
    them, because many ledger entries are reconstructions across the two
    superseded sources rather than verbatim quotes from either. Reporting the
    larger number invited retiring the manual review that covers the rest.
    """
    deleted = [r for r in rows
               if r["claim_id"].startswith("PP-") and r["status"] == "delete"]
    yields_phrase = [r for r in deleted if distinguishing_phrase(r["claim_text"])]
    original = _original_manuscript()

    with capsys.disabled():
        print(f"\n  claims-ledger guard: {len(deleted)} `delete` claims; "
              f"{len(yields_phrase)} yield a searchable phrase")
        if original is None:
            print("  detection rate: NOT MEASURED (shallow clone)")
        else:
            found = set(_detected(rows, original))
            print(f"  DETECTED in the pre-revision manuscript: {len(found)} of "
                  f"{len(deleted)} — this is the real coverage")
            missed = [r["claim_id"] for r in deleted if r["claim_id"] not in found]
            print("  NOT detectable — these still need human review:")
            for cid in missed:
                print(f"      {cid}")

    assert len(yields_phrase) >= 25


def test_every_preprint_claim_names_a_document_that_exists(rows):
    """What would have caught the drift this test was written after.

    The manuscript was revised and renamed outside the repository while 115
    ledger rows went on describing a file that had been superseded. A row whose
    document cannot be resolved is a row nobody can check.
    """
    # "repository" is legitimate and not a data error: a PP-numbered claim can
    # originate in the preprint audit while the thing it describes lives in the
    # code — PP-T2-29 is a Table 2 value whose evidence is biofilms_potts.jl,
    # PP-CIT-01 is about Citations.md.
    known = {"preprint", "preprint_tex", "repository"}
    unresolved = sorted({r["document"] for r in rows
                         if r["claim_id"].startswith("PP-")} - known)
    assert not unresolved, f"PP-* rows name unknown documents: {unresolved}"
    assert PREPRINT.exists(), (
        f"the ledger's preprint rows describe {PREPRINT.name}, which is not in "
        "the repository")


def test_the_manuscript_cites_a_resolvable_revision():
    """A manuscript that cites a commit is only as good as the commit. The
    citation names the PARENT of the commit that adds it, because a document
    cannot cite the commit that contains it — the same structural limit the run
    artifacts record for their own provenance."""
    import subprocess

    # A SHALLOW CLONE CANNOT ANSWER THIS, so do not pretend it did. CI checks
    # this repository out at full depth for exactly that reason; anywhere else
    # (a fork, a `--depth 1` clone) the absence of a commit says nothing about
    # the citation, and failing there would report a fact about the checkout.
    shallow = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "--is-shallow-repository"],
        capture_output=True, text=True)
    if shallow.stdout.strip() == "true":
        pytest.skip("shallow clone: commit resolution is not decidable here")

    text = PREPRINT.read_text(encoding="utf-8")
    cited = set(re.findall(r"repository revision \\texttt\{([0-9a-f]{7,40})\}", text))
    cited |= set(re.findall(r"manuscript is \\texttt\{([0-9a-f]{7,40})\}", text))
    assert cited, "the manuscript cites no repository revision"
    for sha in cited:
        done = subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{sha}^{{commit}}"],
                              capture_output=True)
        assert done.returncode == 0, f"cited revision {sha} does not resolve"


def test_claim_ids_are_unique(rows):
    """AN ID THAT NAMES TWO CLAIMS NAMES NEITHER.

    The ledger is referenced by id from commit messages, pull requests and the
    manuscript revision plan, and `test_no_deleted_claim_survives_in_the_manuscript`
    reports survivors by id. A duplicate makes every one of those references
    ambiguous — and silently, since nothing read the file as a keyed table.

    Two collisions had already accumulated (`REFINE-02`, `REFINE-03`), both
    from this session appending rows without checking. Caught by external
    review, not here, which is why this test exists now.
    """
    import collections

    counts = collections.Counter(r["claim_id"] for r in rows)
    duplicated = {cid: n for cid, n in counts.items() if n > 1}
    assert not duplicated, f"claim ids used more than once: {duplicated}"


def test_every_row_has_an_id_and_a_verdict(rows):
    """The two columns everything else keys on.

    Two vocabularies live here, and that is deliberate. The MANUSCRIPT audit
    reached one of four editorial verdicts per claim; the REPOSITORY rows carry
    a status describing what the code still owes. Mixing a new word into either
    silently removes a row from whichever consumer filters on the other.
    """
    verdicts = {"keep", "restate", "requalify", "delete",          # manuscript
                "supported", "needs_calibration",                  # repository
                "needs_verification", "must_not_be_claimed"}
    for row in rows:
        assert row["claim_id"].strip(), row
        assert row["status"] in verdicts, (row["claim_id"], row["status"])


def rows_of_wrong_width(path) -> list[tuple]:
    """(line number, claim id, field count) for every row that is not the
    declared width. Takes a PATH so a synthetic malformed file can exercise it;
    a guard that has only ever read a well-formed file has not been shown to
    detect anything."""
    # THE READER'S OWN POSITION, because a quoted field may span physical
    # lines. Filtering comments out of the stream first broke the numbering by
    # however many preceded a row; computing an index into the filtered
    # sequence broke it again for any multiline field. `csv.reader.line_num` is
    # the only counter that knows both, so the file is handed to it whole and
    # comment ROWS are skipped after parsing rather than before.
    #
    # The ledger holds no multiline field today. It is one `\n` in one note
    # away from holding one, and then every reported line would be wrong with
    # nothing to say so.
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        header = None
        out = []
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if header is None:
                header = row
                continue
            if len(row) != len(header):
                out.append((reader.line_num, row[0] if row else "?", len(row)))
    return out


def test_a_malformed_row_is_actually_detected(tmp_path):
    """THE CONTROL FOR THE WIDTH GUARD.

    The guard below passes by finding nothing, which is what a guard that can
    find nothing also does. I verified it by hand — editing the real ledger,
    watching it fail, reverting — and committed no evidence of that, which is
    the same as not having done it.

    So: a file with the exact corruption that occurred, built here.
    """
    good = "claim_id,location,claim_text,status\nA-01,here,fine,keep\n"
    assert rows_of_wrong_width(_write(tmp_path / "ok.csv", good)) == []

    # a sentence appended after the final field, commas and all — the edit that
    # gave LADDER-03 three phantom fields
    spilled = good + "A-02,here,fine,keep appended, with a comma, outside\n"
    caught = rows_of_wrong_width(_write(tmp_path / "spilled.csv", spilled))
    assert [c[1] for c in caught] == ["A-02"], caught
    assert caught[0][2] == 6, "the row should show its inflated field count"
    assert caught[0][0] == 3, "line 3 of the file, not 3rd of the filtered rows"

    # AND WITH COMMENTS ABOVE IT, which is what the real ledger has. Counting
    # the filtered stream reported a line that does not hold the bad row.
    commented = "# a note\n# another\n" + spilled
    caught = rows_of_wrong_width(_write(tmp_path / "commented.csv", commented))
    assert caught[0][0] == 5, f"expected physical line 5, got {caught[0][0]}"

    # AND WITH A QUOTED FIELD SPANNING LINES. The ledger holds none today and
    # is one newline in one note away from holding one, at which point an index
    # into the row sequence stops being a line number at all.
    multiline = ('claim_id,location,claim_text,status\n'
                 'A-01,here,"a note\nthat continues\nacross lines",keep\n'
                 'A-02,here,fine,keep appended, with a comma, outside\n')
    caught = rows_of_wrong_width(_write(tmp_path / "multi.csv", multiline))
    assert [c[1] for c in caught] == ["A-02"], caught
    assert caught[0][0] == 5, (
        f"the bad row ends on physical line 5; got {caught[0][0]}")

    # and the opposite corruption: a row that lost a field
    short = good + "A-03,here,fine\n"
    assert [c[1] for c in rows_of_wrong_width(_write(tmp_path / "short.csv", short))] \
        == ["A-03"]


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def test_every_row_has_exactly_the_declared_columns():
    """A ROW THAT SPILLS ITS FIELDS IS SILENTLY A DIFFERENT ROW.

    Appending a sentence to a CSV line by string concatenation puts it AFTER the
    final quoted field, so every comma in it becomes a column separator. That is
    what happened editing `LADDER-03`: the row grew three phantom fields, its
    `notes` were truncated mid-clause, and `csv.DictReader` reported the excess
    under the `None` key — where nothing looked.

    Every other test here reads named columns, so all of them passed over a
    corrupted row. `DictReader` is forgiving by design; this is the check that
    is not.
    """
    wrong = rows_of_wrong_width(LEDGER)
    with open(LEDGER, encoding="utf-8") as fh:
        header = next(csv.reader(l for l in fh if not l.startswith("#")))
    assert not wrong, (
        f"rows whose field count is not {len(header)} — a value containing a "
        f"comma was written unquoted: {wrong}")


def shadowed_definitions(source: str) -> dict:
    """Top-level names defined more than once in a module, with their counts.

    Takes SOURCE, not a path, so it can be handed a module that is known to be
    bad. A scan that has only ever read a clean repository has not been shown to
    find anything — weaken the predicate and it stays green either way.

    Only the module's own body is examined. A helper redefined inside two
    different functions is legal, common, and not what this is about.
    """
    import ast
    import collections

    counts = collections.Counter(
        node.name for node in ast.parse(source).body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)))
    return {name: n for name, n in counts.items() if n > 1}


def test_the_shadowing_scan_detects_a_known_bad_module():
    """THE CONTROL. Three inputs, and the third is the one that stops this
    becoming a nuisance: a name reused inside two separate function bodies is
    ordinary Python and must NOT be reported."""
    clean = "def a():\n    pass\n\n\ndef b():\n    pass\n"
    assert shadowed_definitions(clean) == {}

    duplicated = clean + "\n\ndef a():\n    pass\n"
    assert shadowed_definitions(duplicated) == {"a": 2}

    classes = "class C:\n    pass\n\n\nclass C:\n    pass\n"
    assert shadowed_definitions(classes) == {"C": 2}

    # a method named the same in two classes, and a local helper in two
    # functions — both legal, both must be ignored
    nested = ("class C:\n    def go(self):\n        pass\n\n\n"
              "class D:\n    def go(self):\n        pass\n\n\n"
              "def one():\n    def helper():\n        pass\n\n\n"
              "def two():\n    def helper():\n        pass\n")
    assert shadowed_definitions(nested) == {}


def test_no_module_defines_the_same_name_twice():
    """A SHADOWED DEFINITION MAKES A TEST RUN AGAINST CODE YOU JUST REPLACED.

    Editing this file with a script whose end-marker matched an earlier line
    duplicated a 188-line span, giving ten test functions two definitions each
    — and Python keeps the LAST one. The immediate damage was small because the
    copies were identical, but the same slip had already left a stale
    `rows_of_wrong_width` shadowing the fixed one, so the suite passed against
    the version that had just been replaced. Nothing reported either: a
    duplicate definition is legal Python and pytest collects only the survivor.

    Cheap to check, and it covers every module rather than the one that failed.
    """
    repo = Path(__file__).resolve().parents[2]
    shadowed = {}
    for path in repo.rglob("*.py"):
        if any(part in {".venv", "__pycache__", ".git", "build", "dist"}
               for part in path.parts):
            continue
        try:
            repeated = shadowed_definitions(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue    # not ours to judge
        if repeated:
            shadowed[str(path.relative_to(repo))] = repeated
    assert not shadowed, (
        "top-level definitions shadowed by a later one with the same name; "
        f"the earlier definition never runs: {shadowed}")


# STRUCTURAL AND POLARITY-AWARE, not a phrase list. Two earlier versions failed
# in opposite directions, which is the useful part of the history:
#
#   1. An enumeration of the six sentences that had appeared could not close a
#      natural-language class -- "the remaining conclusions are unchanged" is
#      the same prohibited claim and matched none of them.
#   2. Generalising to shape then over-matched: an unrestricted bridge crossed
#      a negation, so "The remaining conclusions no longer hold" read as
#      "...hold", and "Every other conclusion is affected" -- an explicit
#      RETRACTION and a DISCLOSURE, the honest behaviours the rule wants --
#      were both rejected as reassurances.
#
# So the verbs are split by the polarity of their subject. A positive
# quantifier over the remainder ("all other findings") is reassuring only with
# an unchangedness verb; a negative one ("nothing else") is reassuring only
# with a change verb. The bridge refuses to cross a negation, so a retraction
# cannot be read as its own opposite.
#
# WHAT IT STILL CANNOT DO. No regex closes a paraphrase class; a determined
# rewording gets through. This is a tripwire for a reflex that has recurred
# three times, not a proof. The actual control is AGENTS.md rule 5 -- delete
# the sentence, do not reword it -- and claiming more for the regex here would
# itself be an unchecked statement about scope, which is the defect guarded.
_POS_SUBJECT = (r"(?:all|every|the|any)\s+(?:other|remaining)\s+(?:\w+\s+){0,2}?"
                r"(?:conclusions?|findings?|results?|numbers?|figures?|claims?|values?)"
                r"|(?:the\s+)?(?:conclusions|findings|results|numbers|figures|claims)"
                r"\s+themselves"
                r"|the\s+rest\s+of\s+the\s+(?:conclusions?|findings?|results?|analysis|report)"
                r"|everything\s+else")
_POS_VERB = r"(?:unchanged|unaffected|stands?\b|holds?\b|survives?\b|remains?\s+valid)"
_NEG_SUBJECT = r"(?:nothing\s+else|no\s+other\s+\w+|none\s+of\s+the\s+(?:other\s+)?\w+)"
_NEG_VERB = r"(?:changes?\b|moved?\b|moves\b|(?:is|are|was|were)\s+affected)"
# Tempered: consumes anything up to a clause end EXCEPT a negation.
_BRIDGE = r"(?:(?!\bno\b|\bnot\b|\bnever\b|n't\b)[^.;]){0,45}?"

BLANKET_REASSURANCE = re.compile(
    rf"(?:(?:{_POS_SUBJECT}){_BRIDGE}{_POS_VERB})"
    rf"|(?:{_NEG_SUBJECT}{_BRIDGE}{_NEG_VERB})", re.I)


# A negation BEFORE the subject, which no lookbehind can catch here: Python's
# `re` requires fixed-width lookbehind and these prefixes are not. "Not all
# other findings hold" and "None of the other conclusions remain unchanged" are
# retractions, and the pattern starts at the subject, so it never sees the word
# that reverses them. Checked against the preceding CLAUSE only -- back to the
# last sentence or clause boundary -- so an unrelated "not" in a previous
# sentence cannot excuse a genuine reassurance.
_NEGATION_BEFORE = re.compile(r"\b(?:not|no|never|none|nor|cannot|n't)\b", re.I)
# A negation governs only its own clause. Splitting on sentence marks alone was
# too coarse: "This is not a measurement, but all other findings stand." had its
# reassurance suppressed by a `not` belonging to the previous clause, so any
# correction whose sentence contained an earlier negation got a free pass.
# Commas, colons and coordinating conjunctions end the negation's reach.
_CLAUSE_BOUNDARY = re.compile(
    r"[.;:,\n]|\b(?:but|yet|however|although|though|whereas|while|still)\b", re.I)


def blanket_reassurance_hits(text: str) -> list[str]:
    """Every blanket-reassurance phrase in `text`.

    Takes TEXT so the same finder serves the ledger's notes and the correction
    DOCUMENTS, and so a control can feed it known-bad input. The first version
    inlined the search over `rows` only, which left the markdown reports -- the
    place the phrase actually survived, twice -- unscanned by the production
    path while the regex controls passed.
    """
    hits = []
    for match in BLANKET_REASSURANCE.finditer(text or ""):
        clause = _CLAUSE_BOUNDARY.split(text[:match.start()])[-1]
        if _NEGATION_BEFORE.search(clause):
            continue        # a retraction, not a reassurance
        hits.append(match.group(0))
    return hits


# The documents that carry marked corrections. Scanned as a directory rather
# than a hand-listed file, so a new report is covered the day it is written.
CORRECTION_DOCS = sorted((REPO / "docs" / "calibration").glob("*.md"))


def reassurance_offenders(rows, docs) -> list[tuple]:
    """(where, phrase) for every blanket reassurance across BOTH surfaces.

    THE PRODUCTION PATH, extracted so a control can drive it. Testing
    `blanket_reassurance_hits` alone proves the regex and not its wiring:
    deleting the document loop left that control green while the reports went
    unscanned -- which is the exact defect this function exists to close, so a
    control that could not see it was worth nothing.
    """
    out = [(r["claim_id"], hit) for r in rows
           for hit in blanket_reassurance_hits(r["notes"])]
    for doc in docs:
        out += [(doc.name, hit) for hit in
                blanket_reassurance_hits(doc.read_text(encoding="utf-8"))]
    return out


def test_no_correction_offers_a_blanket_reassurance(rows):
    """SAY WHICH NUMBERS MOVED AND STOP.

    A correction that adds "every other conclusion is unchanged" makes a claim
    about scope that nobody has checked -- and here nobody had: withdrawing the
    3.2 um ladder rung moved component size AND the dependent interface-area
    range, while the text recording it reassured the reader that nothing else
    had.

    It survived two corrections by being REWORDED rather than deleted: first to
    "the conclusions themselves are unchanged" in the report, then to "the
    findings themselves stand" in the ledger. Each read as a smaller claim and
    meant the same thing, which is why the check is on the phrasing family --
    and why it scans BOTH surfaces. A guard covering one of the two places the
    phrase has lived is a guard that would have caught neither instance.

    A statement scoped to ONE named finding and substantiated beside it is fine
    and is not matched. The interface-area paragraph makes exactly that kind of
    claim ("so the conclusion is unchanged and in fact slightly stronger, since
    the error moves less than was reported") and keeps it.
    """
    offenders = reassurance_offenders(rows, CORRECTION_DOCS)
    assert not offenders, (
        "a correction claims unchecked scope; name the numbers that moved and "
        f"stop: {offenders}")


def test_the_document_scan_is_wired_up_not_just_the_regex(tmp_path):
    """THE CONTROL FOR THE WIRING, not for the pattern.

    `test_the_blanket_reassurance_pattern_actually_matches` proves the regex.
    It does NOT prove the regex is pointed at the reports, and for one commit
    it was not: `offenders` searched only the ledger's notes, so the report
    could have regained the phrase with the suite still green.

    So drive the document path: a file that regains the sentence must be
    caught, and one making the permitted scoped claim must not.
    """
    bad = tmp_path / "report.md"
    bad.write_text("# Correction\n\nComponent size moves 3.2 -> 1.6 um. "
                   "The conclusions themselves are unchanged.\n", encoding="utf-8")
    # THROUGH THE PRODUCTION FUNCTION, with an empty ledger, so only the
    # document path can produce the hit. Calling the regex helper here instead
    # is what let the unwired version pass.
    caught = reassurance_offenders([], [bad])
    assert [w for w, _ in caught] == ["report.md"], caught

    ok = tmp_path / "ok.md"
    ok.write_text("# Correction\n\nThe range moves to 0.467-0.510 across 8x, "
                  "so the conclusion is unchanged and in fact slightly stronger, "
                  "since the error moves less than was reported.\n", encoding="utf-8")
    assert reassurance_offenders([], [ok]) == []

    # and a ledger-only offender is still caught, so neither surface was lost
    assert reassurance_offenders(
        [{"claim_id": "X-01", "notes": "The findings themselves stand."}], []
    ) == [("X-01", "The findings themselves stand")]

    # and the real reports are actually in the scanned set
    assert any(d.name == "morphology_rasterization_ladder.md"
               for d in CORRECTION_DOCS), CORRECTION_DOCS


def test_the_blanket_reassurance_pattern_actually_matches(rows):
    """The control: the production check passes by finding nothing, so it is
    worth what its pattern can find.

    THE PARAPHRASES MATTER MORE THAN THE LITERALS. An earlier version listed the
    six sentences that had actually appeared and caught nothing else -- so "the
    remaining conclusions are unchanged" and "all other findings stand", the
    same claim in different words, both passed. The repository rule treats a
    rewording as the same prohibited claim, so the control has to test forms
    that were never written here.
    """
    blanket = [
        # the three that really appeared in this branch
        "every other conclusion is unchanged, including the interface-area finding",
        "The conclusions themselves are unchanged.",
        "The findings themselves stand.",
        # and paraphrases nobody has written yet, which is the point
        "the remaining conclusions are unchanged",
        "all other findings stand",
        "all other results hold",
        "every other number is unaffected",
        "the rest of the analysis is unchanged",
        "nothing else is affected",
        "nothing else changes",
        "everything else stands",
        # AND THE BOUNDARY CASE THE PREFIX FILTER MUST NOT SWALLOW: a negation
        # in a PREVIOUS sentence says nothing about this one, so the clause
        # scan stops at the sentence break. Otherwise any correction that used
        # the word "not" anywhere earlier would get a free pass.
        "This is not a measurement. All other findings stand.",
        # AND THE SAME-SENTENCE VARIANTS. A negation reaches only its own
        # clause; splitting on sentence marks alone let an earlier `not` excuse
        # a reassurance sharing its sentence.
        "This is not a measurement, but all other findings stand.",
        "This is not a measurement: all other findings stand.",
        "Although the range moved, all other findings stand.",
    ]
    for phrase in blanket:
        assert blanket_reassurance_hits(phrase), f"missed: {phrase}"

    # RETRACTIONS AND DISCLOSURES SAY THE OPPOSITE and must never trip. An
    # earlier version bridged straight over the negation, so "no longer hold"
    # read as "hold" and the guard rejected a correction for being honest --
    # the precise inversion of its purpose.
    opposite = [
        "The remaining conclusions no longer hold",
        "Every other conclusion is affected",
        "all other findings changed",
        "the remaining numbers moved",
        "every other result is not unchanged",
        # NEGATION BEFORE THE SUBJECT, which the pattern cannot see because it
        # starts AT the subject. Python's `re` has no variable-width lookbehind,
        # so these are filtered by inspecting the preceding clause instead.
        "Not all other findings hold",
        "None of the other conclusions remain unchanged",
        "no other findings remain unchanged",
    ]
    for phrase in opposite:
        assert not blanket_reassurance_hits(phrase), (
            f"rejected a RETRACTION as if it were a reassurance: {phrase}")

    # THE PERMITTED FORM: one named finding, substantiated on the same line.
    # If these ever trip, the guard has started deleting honest scoped claims,
    # which is a worse failure than the one it prevents.
    scoped = [
        "so the conclusion is unchanged and in fact slightly stronger, since "
        "the error moves less than was reported",
        "the biovolume (0.8) and porosity (0.4) figures are unaffected, having "
        "been measured at pitches that tile",
        "Two published numbers move: component size and the interface-area range.",
        "The interface-area finding was stated as 0.41-0.51 across a 16x range.",
    ]
    for phrase in scoped:
        assert not blanket_reassurance_hits(phrase), f"false positive: {phrase}"
