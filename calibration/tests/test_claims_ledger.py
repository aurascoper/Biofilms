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

EVERY `delete` ROW IS READ IN THE DOCUMENT IT NAMES. The first version filtered
`claim_id.startswith("PP-")` and read only the manuscript, which meant the 25
`delete` verdicts reached on anything else were carried by a test searching a
file those claims had never been in — a check that could not fail. Scope now
comes from the `document` column.

WHAT THIS DOES NOT CLAIM. Detection is partial. Of the 34 manuscript rows, 20
are findable in a document known to contain them; the rest are paraphrases with
elisions and no verbatim fragment long enough to search for. HANDOUT-01 is not
findable either, for a different reason: its claim is an equation, and
`normalise_markup` reduces LaTeX maths to nothing searchable. All of them are
NAMED in the coverage output rather than passed silently over — the same
discipline as `reference_d_status.enforcement_report()` naming the declared
criteria that have no consumer. A test that implied full coverage would be worse
than no test, because it would retire the manual review that actually covers
them.
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


def normalise_plaintext(source: str) -> str:
    """Extracted figure text, whitespace-collapsed and lowercased.

    NOT `normalise_markup`. That function strips `%` to end of line as a LaTeX
    comment, and figure text carries real percent signs — `c_mean = 1.2% c_ext`
    — so running it here would silently delete the rest of any line holding one.
    A normaliser aimed at the wrong language is a guard reading a document it
    cannot see.
    """
    return " ".join(source.split()).lower()


def _normalised(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    return (normalise_plaintext if path.suffix == ".txt"
            else normalise_markup)(text)


def _preprint_text() -> str:
    return normalise_markup(PREPRINT.read_text(encoding="utf-8"))


# THE ROW SAYS WHICH DOCUMENT IT IS ABOUT; the guard must not decide that from
# the id. Filtering on `claim_id.startswith("PP-")` meant every `delete` verdict
# reached on anything but the manuscript — the handout's two among them — was
# carried by a test that read a file those claims were never in, and so could
# not fail. Two aliases and one pseudo-document are all that stand between the
# `document` column and a path.
DOCUMENT_ALIASES = {
    "preprint": PREPRINT,
    "preprint_tex": PREPRINT,
    "README": REPO / "README.md",
}
# Not files: a claim about the repository as a whole, or about what was said to
# a collaborator. A row naming one cannot be searched, so it is NAMED rather
# than skipped in silence — an absent row and an unsearchable one otherwise
# print the same nothing.
PSEUDO_DOCUMENTS = {"repository", "correspondence"}

# NOT SOURCE, AND NOT IN THE REPOSITORY. `artifacts/` is gitignored
# (`.gitignore:76`) and nothing under it is tracked, so a row naming a pilot
# artifact resolves to a file only on a machine where that pilot has been run.
# Both document-resolution tests below were therefore passing on local state and
# would have failed on CI or any clean checkout — a check whose outcome depends
# on untracked files is rule 1 wearing rule 2's clothes. They are DECLARED here
# rather than removed from the ledger: the rows are real verdicts on real
# numbers, they are simply unenforceable from a clean tree, and saying that is
# the point. When the artifact happens to be present it is still read.
GENERATED_ARTIFACTS = {
    "artifacts/pilot/openmc_nested_pilot_budget.json",
    "artifacts/pilot/openmc_nested_pilot_verdict.json",
}


FIXTURES = Path(__file__).resolve().parent / "fixtures"

# The figures are the only PUBLISHED artifacts that are not source. A committed
# PDF carries claim text that no suite could read: `tests/runtests.jl` splits the
# monolith at `#  13. Figure export` and loads only what is above the cut, so
# `export_figures` is unreachable from Julia by construction. The generator was
# corrected on 2026-08-14 and the ledger recorded the verdict three times
# (RM-G08-01, RM-G10-01, PP-65-08, the last saying in as many words that "the
# committed PNGs still carry the old reversed zone labels"); the artifact went on
# saying the retracted thing for two weeks because nothing could fail.
FIGURES = REPO / "preprint" / "figures"


def document_path(document: str) -> Path | None:
    if document in PSEUDO_DOCUMENTS:
        return None
    path = DOCUMENT_ALIASES.get(document, REPO / document)
    # A row may name the artifact a reader actually opens. PDF text is
    # compressed, so searching the bytes finds nothing; the committed `.txt`
    # sidecar is what gets read, and the hash test below is what stops a stale
    # sidecar from standing in for a fresh PDF.
    if path.suffix == ".pdf":
        return path.with_suffix(".txt")
    return path


def deleted_rows(rows):
    return [r for r in rows if r["status"] == "delete"]


# THE CONTROLS ARE COMMITTED FILES, NOT GIT OBJECTS. Each of these was once
# recovered with `git show <sha>:<path>` — 5980dc5, 9319d43, and a third that
# would have been e24dbec — and each time the same exposure was written down and
# left in place: a squash-merge makes the pinned commit unreachable, `git show`
# fails, and the control degrades into a skip. Rule 2 says to read a skip as
# uncovered surface, and a control that can quietly stop controlling is rule 1
# again one level up. A committed file cannot become unreachable. Provenance is
# in `fixtures/README.md`.
ORIGINAL_FIXTURE = "modeling_radiotrophic_fitness_prerevision.md"   # was 5980dc5


def _fixture(name: str) -> str:
    """A known-bad document, read from disk. Never skips: if a control file is
    missing that is a failure, not an environment quirk."""
    path = FIXTURES / name
    assert path.is_file(), (
        f"negative-control fixture {name} is missing. The guard it feeds cannot "
        "fail without it, which makes it not a guard. See fixtures/README.md.")
    return normalise_markup(path.read_text(encoding="utf-8", errors="replace"))


def _original_manuscript() -> str:
    return _fixture(ORIGINAL_FIXTURE)


def _detected(rows, text, documents=("preprint", "preprint_tex")) -> list[str]:
    """Which `delete` rows for `documents` are findable in `text`.

    Takes the document set explicitly because it is used on CONTROLS — a file
    known to contain the claims — and a control is only valid for the rows whose
    document it actually is.
    """
    out = []
    for row in deleted_rows(rows):
        if row["document"] not in documents:
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


def test_no_deleted_claim_survives_in_the_document_it_names(rows):
    """THE ONE THAT MATTERS. Every `delete` verdict was reached because the
    claim was unsupported, fabricated, or checkable and false — a Hamiltonian
    that conserves nothing, a symplectic integrator that is not used, a
    quantum-mechanical noise mechanism that is a Gaussian keyed to a diffusion
    coefficient. They must not come back.

    Each row is searched in the document it names, not in the manuscript alone.
    """
    cache: dict[Path, str] = {}
    survivors = []
    for row in deleted_rows(rows):
        path = document_path(row["document"])
        if path is None or not path.is_file():
            continue                                    # named by the test below
        if path not in cache:
            cache[path] = _normalised(path)
        phrase = distinguishing_phrase(row["claim_text"])
        if phrase and phrase.lower() in cache[path]:
            survivors.append(f"{row['claim_id']} in {row['document']}: {phrase[:90]}")
    assert not survivors, (
        "claims marked `delete` in the ledger are present in the document that "
        "carried them:\n  " + "\n  ".join(survivors))


def test_the_guard_actually_detects_deleted_claims(rows):
    """THE TEST THAT PROVES THE OTHER ONE MEANS SOMETHING.

    `test_no_deleted_claim_survives_in_the_document_it_names` passes when it finds
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
    found = _detected(rows, _original_manuscript())
    assert len(found) >= 18, (
        f"the guard detected only {len(found)} deleted claims in the "
        "pre-revision manuscript, which contains them. It has stopped "
        "guarding — most likely the markup normalisation regressed.")


# The handout as first committed, before the azide line and the sorption box
# were corrected. Same role as ORIGINAL_FIXTURE, for the other document the
# ledger carries `delete` verdicts on.
HANDOUT_FIXTURE = "wan_meeting_handout_prefix.tex"          # was 9319d43
HANDOUT_PATH = "preprint/wan_meeting_handout.tex"


def test_the_guard_detects_the_handout_claims(rows):
    r"""THE SAME CONTROL, FOR THE OTHER DOCUMENT.

    Widening the filter from `PP-` to the `document` column is only worth
    something if the guard can actually bite on what it newly reads. Run it
    against the pre-correction handout, which contains both deleted claims.

    HANDOUT-01 is NOT expected here and is not counted: its claim is an
    equation, and `normalise_markup` reduces LaTeX maths to nothing searchable
    (`X_{\max}` becomes a bare `x`). A prose-phrase guard cannot check a
    formula, so that row stays with a human — see the coverage report.
    """
    found = _detected(rows, _fixture(HANDOUT_FIXTURE), documents=(HANDOUT_PATH,))
    assert "HANDOUT-02" in found, (
        "the guard found no deleted handout claim in the pre-correction "
        f"handout, which contains them (found: {found}). It is not guarding "
        "that document.")


def test_every_deleted_claim_names_a_document_the_guard_can_read(rows):
    """A ROW THE GUARD SILENTLY SKIPS IS NOT COVERED, AND MUST SAY SO.

    `test_no_deleted_claim_survives_in_the_document_it_names` passes over any
    row whose `document` does not resolve to a file. That is legitimate for the
    declared pseudo-documents and for nothing else: a typo, a renamed file or a
    new alias would otherwise retire a `delete` verdict without a word.
    """
    unreadable = sorted({r["document"] for r in deleted_rows(rows)
                         if (lambda q: q is None or not q.is_file())(
                             document_path(r["document"]))}
                        - PSEUDO_DOCUMENTS - GENERATED_ARTIFACTS)
    assert not unreadable, (
        "`delete` rows name documents that cannot be read, so those verdicts "
        f"are unguarded: {unreadable}. Add a path alias, fix the row, or "
        "declare it in PSEUDO_DOCUMENTS / GENERATED_ARTIFACTS.")


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
    deleted = [r for r in deleted_rows(rows)
               if r["document"] in ("preprint", "preprint_tex")]
    yields_phrase = [r for r in deleted if distinguishing_phrase(r["claim_text"])]
    original = _original_manuscript()

    with capsys.disabled():
        print(f"\n  claims-ledger guard: {len(deleted)} `delete` claims; "
              f"{len(yields_phrase)} yield a searchable phrase")
        found = set(_detected(rows, original))
        print(f"  DETECTED in the pre-revision manuscript: {len(found)} of "
              f"{len(deleted)} — this is the real coverage")
        missed = [r["claim_id"] for r in deleted if r["claim_id"] not in found]
        print("  NOT detectable — these still need human review:")
        for cid in missed:
            print(f"      {cid}")

        # THE OTHER DOCUMENTS ARE READ NOW, so their coverage is reported too —
        # a row the widened filter reaches but the prose guard cannot match is
        # still uncovered, and saying "34 delete claims" over it would hide that.
        others = [r for r in deleted_rows(rows)
                  if r["document"] not in ("preprint", "preprint_tex")]
        print(f"  outside the manuscript: {len(others)} `delete` claims in "
              f"{len({r['document'] for r in others})} documents")
        print("  NOT textually detectable there — equations and code, not prose:")
        print("      HANDOUT-01 — the claim is a formula; normalise_markup "
              "reduces LaTeX maths to nothing searchable")
        print("      FIG-01, FIG-02, FIG-05 — in-plot labels; no run of "
              "MIN_WORDS survives the comma and paren split. Covered by "
              "RETRACTED_IN_FIGURES, not by the phrase guard")
        print("      FIG-09, FIG-10 — full sentences that DO yield a phrase, and "
              "the phrase is still unreachable: pdftotext -layout reads the "
              "three columns across, so the run is contiguous nowhere. Also "
              "covered by RETRACTED_IN_FIGURES. Listing them beside the three "
              "above is the point — 'yields a phrase' is not 'is covered'")

    assert len(yields_phrase) >= 25


def test_every_claim_names_a_document_that_exists(rows):
    """What would have caught the drift this test was written after.

    The manuscript was revised and renamed outside the repository while 115
    ledger rows went on describing a file that had been superseded. A row whose
    document cannot be resolved is a row nobody can check.
    """
    # `PSEUDO_DOCUMENTS` is legitimate and not a data error: a PP-numbered claim
    # can originate in the preprint audit while the thing it describes lives in
    # the code — PP-T2-29 is a Table 2 value whose evidence is biofilms_potts.jl,
    # PP-CIT-01 is about Citations.md. `document_path` returns None for those,
    # and only a resolvable-but-missing path is a fault here.
    unresolved = sorted({r["document"] for r in rows
                         if (lambda q: q is not None and not q.is_file())(
                             document_path(r["document"]))}
                        - GENERATED_ARTIFACTS)
    assert not unresolved, f"rows name documents that do not exist: {unresolved}"

    # AND THE DECLARATION MUST STAY TRUE. Listing a document here excuses it
    # from the check above, so an entry that is no longer a gitignored artifact
    # — committed since, or renamed — would silently excuse nothing while
    # looking like it still did.
    misdeclared = sorted(d for d in GENERATED_ARTIFACTS
                         if d not in {r["document"] for r in rows})
    assert not misdeclared, (
        f"GENERATED_ARTIFACTS names documents no ledger row uses: {misdeclared}")
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
    manuscript revision plan, and `test_no_deleted_claim_survives_in_the_document_it_names`
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
# Commas, colons and conjunctions end the negation's reach. The first version
# said "coordinating conjunctions" and then listed mostly SUBORDINATING ones --
# although, though, whereas, while -- so "not a measurement and all other
# findings stand" kept its free pass. All seven coordinators are here now
# (FANBOYS), which is what the comment claimed all along.
_CLAUSE_BOUNDARY = re.compile(
    r"[.;:,\n]"
    r"|\b(?:and|but|or|nor|for|yet|so)\b"                    # coordinating
    r"|\b(?:however|although|though|whereas|while|still)\b",  # and the rest
    re.I)


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
        "This is not a measurement and all other findings stand.",
        "This is not a measurement or all other findings stand.",
        "This is not a measurement so all other findings stand.",
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


# ---------------------------------------------------------------------------
# THE FIGURES. A published artifact nothing could read.
# ---------------------------------------------------------------------------

# §2.6 of the manuscript: "Radiotrophy is not established for any of the seven
# species modelled." The word is retracted prose, so it must not survive inside
# an image either. This floor is coarse on purpose — it covers `fig1`, whose
# in-plot labels were "radiotrophic niche" and "radiosensitive core", two words
# each and so under MIN_WORDS: `distinguishing_phrase` cannot search them and
# the phrase guard never will. Coverage of that figure is this list, not the
# ledger row.
#
# "gy cumulative" is here for FIG-05, and it is a UNITS claim rather than a
# phenotype one — the list is "terms that must not survive inside an image", and
# RM-KR-01's verdict is NEVER PRINT Gy HERE, because D_cum = Ddot_R*t is
# dimensionless model time and no calibration converts it. fig3 carried
# "(50 Gy cumulative)" for fifteen days after d53f236 removed it from the
# generator, because main_coupled() hits RADIODIALYSIS: BLOCKED and nothing
# could rebuild the artifact. Three words, so the phrase guard could never have
# reached it either: without this term FIG-05 would have been covered by
# nothing, while sitting in a list named for rows the vocabulary covers.
#
# "charge-excluded" and "same trap" are FIG-09 and FIG-10, and they are here for
# a reason the first three do not share: THEIR CLAIM TEXT DOES YIELD A SEARCHABLE
# PHRASE AND THE PHRASE IS STILL NOT FINDABLE IN THE ARTIFACT. `pdftotext
# -layout` reads a three-column figure across, so "Bromide is size- and
# charge-excluded from fine porosity." comes back as "...charge-excluded
# breakthrough lag gives a lumped run alongside every measurement. from fine
# porosity." — the phrase interleaved with two other columns and contiguous
# nowhere. Verified against the real pre-correction sidecar, where `_detected`
# returns [] for both rows.
#
# THAT DISTINCTION IS THE WHOLE DEFECT THIS LIST NOW CLOSES. "Yields a phrase"
# was read as "is covered", and it is not the same predicate: coverage is a
# property of the ARTIFACT, not of the claim text. A hand-run control that
# appends the sentences as clean contiguous lines confirms a guard against an
# input shape a column layout cannot produce — a control passing on a mechanism
# that never ran, which is what the fig1/fig2 control below was already written
# to avoid. Between c2219a2 and this commit, FIG-09 and FIG-10 were covered by
# nothing at all.
RETRACTED_IN_FIGURES = ("radiotroph", "gy cumulative",
                        "charge-excluded", "same trap")


def figure_sidecars() -> list[Path]:
    return sorted(FIGURES.glob("*.txt"))


def retracted_vocabulary_hits(paths) -> list[tuple]:
    """(figure, term) for every retracted term found. Takes PATHS so the
    control can drive it with the pre-correction artifacts."""
    return [(p.name, term) for p in paths
            for term in RETRACTED_IN_FIGURES
            if term in normalise_plaintext(p.read_text(encoding="utf-8",
                                                       errors="replace"))]


def test_no_committed_figure_asserts_a_retracted_phenotype():
    """WHAT THE PROSE AUDIT DID NOT REACH.

    `fig2_melanin_accumulation.pdf` said, inside the image, "C. neoformans,
    C. sphaerospermum are radiotrophic (melanin-mediated energy gain)" — under a
    caption disowning "radiation-derived energy production" and against §2.6.
    `fig1` said "radiotrophic niche".

    Neither was unnoticed. The generator was corrected on 2026-08-14 and the
    ledger recorded the verdict three times. The artifact went on saying it
    because no test could open a figure: the claims guard read one `.tex`, and
    `tests/runtests.jl` splits the monolith above `#  13. Figure export`, so the
    Julia suite cannot reach `export_figures` at all. Fixing a generator does
    not fix what was already committed, and nothing said so.
    """
    hits = retracted_vocabulary_hits(figure_sidecars())
    assert not hits, (
        "committed figures assert a phenotype §2.6 retracts; the generator was "
        f"fixed but the artifact was not regenerated: {hits}")


def test_the_figure_guard_detects_the_pre_correction_artifacts():
    """THE CONTROL: run the floor against the artifacts known to carry the claim.

    NEITHER FIGURE IS REACHABLE BY THE PHRASE GUARD, and that is why this floor
    exists rather than duplicating it. `distinguishing_phrase` splits on commas
    and parentheses and needs MIN_WORDS=5 from one run. Fig 2's annotation —
    "C. neoformans, C. sphaerospermum are radiotrophic (melanin-mediated energy
    gain)" — yields runs of 2, 4 and 3 words, and fig 1's labels are two words
    each. So `_detected` returns nothing for FIG-01 and FIG-02 no matter what,
    exactly as it does for HANDOUT-01's equation, and the vocabulary floor is
    the whole of their coverage. Asserting a phrase hit here would have made the
    control pass on a mechanism that never ran.

    EXTENDED TO THE FOURTH ARTIFACT, and the glob had to widen to reach it.
    `phase2_diffusion_cell_prefix.txt` does not match `fig*`, so while this
    control globbed `fig*_prefix.txt` a committed known-bad sat in the fixtures
    directory covered by nothing — the count assertion pinned three and passed,
    which is the shape of every other defect on this page: a bound checked
    against the set it was derived from. FIG-09 and FIG-10 differ from the first
    three in that their claim text DOES yield a phrase; it is simply unfindable
    in a column layout. See RETRACTED_IN_FIGURES.
    """
    prefix = sorted(FIXTURES.glob("*_prefix.txt"))
    assert len(prefix) == 4, f"expected four pre-correction figures, got {prefix}"

    hits = retracted_vocabulary_hits(prefix)
    assert {p for p, _ in hits} == {f.name for f in prefix}, (
        "the vocabulary floor missed a figure that contains the retracted "
        f"phenotype: found {hits}")

    # and it is not matching everything it is handed
    assert retracted_vocabulary_hits(figure_sidecars()) == []


def test_no_figure_row_is_reachable_by_the_phrase_guard(rows):
    """COVERAGE IS A PROPERTY OF THE ARTIFACT, NOT OF THE CLAIM TEXT.

    The first version of this test split the figure rows into a vocabulary tier
    and a phrase tier by asking whether `distinguishing_phrase` returned
    anything, and put FIG-09 and FIG-10 in the phrase tier because their claim
    text yields a full sentence. THAT PREDICATE IS WRONG AND THE ROWS WERE
    COVERED BY NOTHING. `pdftotext -layout` reads a three-column figure across,
    so the sentence comes back interleaved with the other two columns and is
    contiguous nowhere; `_detected` returns [] for both against the real
    pre-correction sidecar. A phrase that exists in the ledger and cannot occur
    in the artifact is not a guard, and it looked exactly like one.

    The hand-run that "confirmed" the phrase guard bites appended both sentences
    as clean lines, which is a shape a column layout cannot produce -- a control
    passing on a mechanism that never ran, which is the same defect the fig1/fig2
    control below was written to avoid. Hence this test asks the artifacts.

    SCOPE: the `delete` rows whose `document` is under `preprint/figures/`,
    evaluated against the four committed `*_prefix.txt` sidecars and nothing
    else. It says nothing about rows in any other document.
    """
    known_bad = sorted(FIXTURES.glob("*_prefix.txt"))
    assert len(known_bad) == 4, (
        f"expected four pre-correction sidecars, got {[p.name for p in known_bad]}")
    joined = "\n".join(
        normalise_plaintext(p.read_text(encoding="utf-8", errors="replace"))
        for p in known_bad)

    figure_rows = [r for r in deleted_rows(rows)
                   if r["document"].startswith("preprint/figures/")]
    assert {r["claim_id"] for r in figure_rows} == {
        "FIG-01", "FIG-02", "FIG-05", "FIG-09", "FIG-10"}, (
        f"the figure row set moved: {sorted(r['claim_id'] for r in figure_rows)}. "
        "A new figure `delete` row needs a known-bad sidecar before it can be "
        "said to be covered by anything.")

    reachable = []
    for row in figure_rows:
        phrase = distinguishing_phrase(row["claim_text"])
        if phrase and phrase.lower() in joined:
            reachable.append(row["claim_id"])
    assert not reachable, (
        f"{reachable} are now findable as contiguous phrases in a known-bad "
        "figure sidecar, so the phrase guard genuinely covers them and this "
        "test's premise has changed. That is an improvement -- record which "
        "rows moved and why the layout now preserves the run.")

    # THE OTHER HALF, OR "COVERED BY NOTHING" AND "COVERED BY THE FLOOR" ARE THE
    # SAME RESULT HERE. Asserting only that no phrase reaches the artifacts is
    # satisfied by a word list that fires on none of them.
    covered = {name for name, _ in retracted_vocabulary_hits(known_bad)}
    assert covered == {p.name for p in known_bad}, (
        "the vocabulary floor is the whole of these rows' coverage and it does "
        f"not fire on every known-bad artifact: {sorted(covered)}")


def test_every_committed_figure_matches_its_hash():
    """THE HALF THAT MAKES THE SIDECAR WORTH ANYTHING.

    Everything above reads `<figure>.txt`, not the PDF a reader opens. Without
    this, regenerating a figure and forgetting its sidecar leaves the suite
    green on stale-but-clean text while the PDF still carries the claim — the
    artifact-versus-source split reproduced inside the guard written to close
    it. So the PDF's own bytes are pinned, in pure stdlib, and a figure that
    moves without its sidecar fails here.
    """
    import hashlib

    stale = []
    for pdf in sorted(FIGURES.glob("*.pdf")):
        recorded = pdf.with_suffix(".sha256")
        assert recorded.is_file(), f"{pdf.name} has no committed hash"
        actual = hashlib.sha256(pdf.read_bytes()).hexdigest()
        if actual != recorded.read_text(encoding="utf-8").strip():
            stale.append(pdf.name)
        assert pdf.with_suffix(".txt").is_file(), f"{pdf.name} has no text sidecar"
    assert not stale, (
        "these figures changed without their sidecar and hash being "
        f"regenerated, so every text check above read the OLD figure: {stale}. "
        "Re-run the export, then `pdftotext -layout` and `sha256sum` for each.")


def test_the_sidecars_are_what_pdftotext_actually_produces():
    """TIER 2: the sidecar is the PDF's text, not merely a file beside it.

    The hash test proves the PDF has not moved since the sidecar was written; it
    cannot prove the sidecar was ever a faithful extraction. This can, where
    `pdftotext` exists. It is not a bare skip when absent — the hash test and
    the vocabulary floor both still run — which is the difference rule 2 asks
    for between a gap and a hidden one.
    """
    import shutil
    import subprocess

    if shutil.which("pdftotext") is None:
        pytest.skip("pdftotext absent: hash + vocabulary tiers still cover this")

    drifted = []
    for pdf in sorted(FIGURES.glob("*.pdf")):
        got = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                             capture_output=True, text=True)
        assert got.returncode == 0, f"pdftotext failed on {pdf.name}"
        if normalise_plaintext(got.stdout) != normalise_plaintext(
                pdf.with_suffix(".txt").read_text(encoding="utf-8")):
            drifted.append(pdf.name)
    assert not drifted, (
        f"committed sidecars are not this PDF's text: {drifted}")


# ---------------------------------------------------------------------------
# A COINED TERM MUST NOT APPEAR BEFORE THE PARAGRAPH THAT DISOWNS IT
# ---------------------------------------------------------------------------

COINAGE = "radiodialysis"
DISCLAIMER = "coined for this three-equation system"

# A filename and a cross-reference label are not the coined term being used as
# though established: `biofilms_radiodialysis.R` is what the file is called, and
# `\label{}` never reaches a reader at all.
_EXEMPT = re.compile(r"\\(?:texttt|label|ref|eqref)\{[^}]*\}")


def coinage_uses_before_disclaimer(source: str) -> list[tuple[int, str]]:
    """(line number, line) for every reader-visible use of the coined term
    above the paragraph that declares it coined.

    Takes SOURCE rather than reading the file, so the control can drive it with
    a document known to violate the rule.
    """
    lines = source.split("\n")
    cut = next((i for i, l in enumerate(lines) if DISCLAIMER in l), None)
    if cut is None:
        return [(0, "the disclaimer paragraph is missing entirely")]
    return [(i + 1, l) for i, l in enumerate(lines[:cut])
            if COINAGE in _EXEMPT.sub("", l).lower()]


def test_the_coinage_never_precedes_its_disclaimer():
    """THE FIGURE DEFECT'S SHAPE, IN PROSE.

    §3.11 disowns `radiodialysis` in a full paragraph — "a targeted search finds
    no scientific usage of it outside this work". The term then appeared ten
    times, and the abstract called it "a radial radiodialysis solver" with no
    signal that the word was coined here. Caveat in the body, claim in the part
    people scan: the same split that left a retracted phenotype inside a figure
    while the prose around it was clean.

    The rule is positional and so is the check. Anywhere after the disclaimer
    the term is introduced and qualified; anywhere before it, it reads as
    established vocabulary.
    """
    hits = coinage_uses_before_disclaimer(PREPRINT.read_text(encoding="utf-8"))
    assert not hits, (
        f"'{COINAGE}' is used as established vocabulary before the paragraph "
        f"that declares it coined: {hits}")


def test_the_coinage_check_detects_a_use_above_the_disclaimer():
    """THE CONTROL. Without it this test passes on any document that never says
    the word, which is every document in the repository but one."""
    planted = (f"We solve the {COINAGE} system on a radial mesh.\n"
               f"A term {DISCLAIMER} and is not established.\n"
               f"The {COINAGE} solver is described in Section 6.\n")
    hits = coinage_uses_before_disclaimer(planted)
    assert [h[0] for h in hits] == [1], (
        "the check must flag the line above the disclaimer and only that one; "
        f"got {hits}")
    assert not coinage_uses_before_disclaimer(
        f"A term {DISCLAIMER}.\nThe {COINAGE} solver runs.\n")
    assert not coinage_uses_before_disclaimer(
        f"See \\texttt{{biofilms_{COINAGE}.R}} and \\label{{sec:{COINAGE}}}.\n"
        f"A term {DISCLAIMER}.\n"), "a filename and a label are not a use"


# ---------------------------------------------------------------------------
# SOME SENTENCES HAVE TO STAY
# ---------------------------------------------------------------------------

# The guard's other document tests all assert ABSENCE — a `delete` verdict is
# enforceable precisely because absence is its criterion. That leaves the
# opposite failure uncovered: a sentence carried deliberately, whose removal
# nothing notices. These rows name sentences that are load-bearing because of
# what they PREVENT the document from implying, which is exactly the kind that
# reads like a hedge on an editing pass and gets trimmed.
LOAD_BEARING = ("HOFFMAN-10",)


def load_bearing_absences(rows, read=None) -> list[tuple[str, str]]:
    """(claim_id, phrase) for each row whose sentence is no longer in its
    document. `read` overrides document lookup so the control can supply text."""
    out = []
    for row in rows:
        if row["claim_id"] not in LOAD_BEARING:
            continue
        phrase = distinguishing_phrase(row["claim_text"])
        assert phrase, (
            f"{row['claim_id']} is declared load-bearing but its claim_text has "
            f"no run of {MIN_WORDS}+ words to search for. Rephrase the row — "
            "the word floor is what stops a match being a coincidence.")
        path = document_path(row["document"])
        assert path is not None and path.is_file(), (
            f"{row['claim_id']} names {row['document']}, which does not resolve")
        text = read(path) if read else _normalised(path)
        if phrase.lower() not in text:
            out.append((row["claim_id"], phrase))
    return out


def test_a_load_bearing_sentence_survives_in_its_document(rows):
    """The Hoffman memo's geometry caveat.

    The implemented domain is a cylinder in water with a Robin boundary — no
    metal, no planar film, no interface. Without that sentence page 3 offers a
    dose field beneath a film on steel as something the code does, which is a
    geometry it does not have. It is also the sentence that reads like a hedge
    and goes first when someone tightens the page.
    """
    missing = load_bearing_absences(rows)
    assert not missing, (
        "a sentence the ledger records as load-bearing is no longer in the "
        f"document that needs it: {missing}")


def test_the_load_bearing_check_detects_a_deleted_sentence(rows):
    """THE CONTROL: strike the phrase from the document and require a report."""
    def struck(path):
        text = _normalised(path)
        for row in rows:
            if row["claim_id"] in LOAD_BEARING:
                text = text.replace(distinguishing_phrase(row["claim_text"]).lower(), "")
        return text

    assert [c for c, _ in load_bearing_absences(rows, read=struck)] == list(LOAD_BEARING), (
        "removing the phrase must make the check fail; if it does not, the "
        "check is reading something other than the sentence it names")


# --------------------------------------------------------- retracted citations
#
# THE THIRD TIER, AND IT EXISTS BECAUSE A CITATION CANNOT REACH THE FIRST TWO.
# `distinguishing_phrase` splits on commas and needs MIN_WORDS from one unbroken
# run. A bibliographic entry -- "Diele, Marangi, Ragni (2015), Math. Comput.
# Simul. 110, 40-52, doi ..." -- is all commas, so it yields '' and PP-REF-01 and
# PP-REF-02 were `delete` verdicts covered by nothing. Both were prescribed in the
# ledger and sat unapplied in the manuscript until 2026-08-31.
#
# NOTE THIS IS A SECOND CAUSE, NOT THE SAME ONE AS THE SOURCE-COMMENT GAP. A
# verdict is also unenforced when the claim lives in a .jl comment rather than in
# the file the row's `document` column names -- that is RM-G04-01, and it needs a
# RETRACTED_IN_SOURCES scan over SIM_FILES. Widening which files are scanned would
# not reach the two reference rows, because their invisibility is a property of the
# phrase EXTRACTOR. One symptom, two causes; this closes exactly one of them.
# THIS GUARD REACHES PP-REF-01 AND PP-REF-02 AND NOTHING ELSE.
RETRACTED_CITATIONS = {
    "10.1002/mma.3237": "PP-REF-01: resolves to an unrelated 3D MHD regularity paper",
    "10.1016/j.matcom.2014.02.006": "PP-REF-02: resolves to an unrelated QP paper",
}

# USE VERSUS MENTION, AND THE SET IS EXPLICIT NAMES RATHER THAN A GLOB.
# Recording that a DOI was withdrawn requires naming it, so the ledger and the
# red-team documents necessarily contain these strings. A whole-tree scan fails on
# its first run against the very files that record the withdrawal -- the shape of a
# test whose own control fixture contains the string it forbids.
#
# The resolution is an ALLOW-LIST, inverted: scan everything, and require every file
# carrying a withdrawn DOI to be declared here. A new document quoting one FAILS
# until it is consciously added, rather than passing because it fell outside a
# curated scan list. `docs/research/*_redteam.md` was rejected as a glob for exactly
# that reason: it would auto-admit every future red-team file, which is the property
# the allow-list was chosen to avoid.
#
# BOTH ENTRIES BELOW THE FIRST THREE WERE FOUND BY THIS GUARD ON ITS FIRST RUN, which
# is the allow-list earning its place: neither was anticipated, and a curated scan
# list would have silently omitted both instead of demanding a decision.
CITATION_RECORDING_FILES = {
    "data/claims_ledger.csv",
    "docs/research/session_claims_2026-08-24_redteam.md",
    "docs/research/external_reviews_2026-08-31_redteam.md",
    "calibration/tests/test_claims_ledger.py",
    # The pre-revision manuscript, which is ORIGINAL_FIXTURE -- the committed
    # known-bad this file already uses as its detection floor. It carries both
    # withdrawn DOIs BECAUSE it is the artifact from before the correction, and
    # fixtures/README.md says of these files: "Do not regenerate or 'fix' these
    # files. Their value is that they are wrong."
    "calibration/tests/fixtures/modeling_radiotrophic_fitness_prerevision.md",
}


def files_carrying_retracted_citations() -> dict:
    """{repo-relative path: [dois]} for every file containing a withdrawn DOI.

    SCOPE: every readable file in the repository except `.git`. Shared by the check
    and its control so the control cannot pass against a regressed check.
    """
    out: dict[str, list[str]] = {}
    for path in REPO.rglob("*"):
        # THE EXCLUSION IS "DERIVED FROM A SOURCE ALREADY SCANNED", NOT "IS A BUILD
        # ARTIFACT", AND THE TWO RULES HAVE DIFFERENT FUTURE COVERAGE. The first
        # version of this skipped `path.suffix in {.pyc, .pdf, .png, .so, .o}`,
        # which would silently admit any future artifact of those types carrying a
        # forbidden string -- and `preprint/*.pdf` is the manuscript's own build
        # output, so the rule excluded the one artifact a reader actually receives.
        # Measured: with no type exclusion at all, 14 PDFs in the tree produce zero
        # hits, because PDF text is compressed. The suffix rule bought nothing and
        # cost coverage.
        #
        # A `.pyc` is excludable for a narrower and checkable reason: it is a
        # compiled copy of a `.py` that is itself in scope, so its content is
        # already being read at its source. That justification is ASSERTED below
        # rather than assumed -- if the sibling source is missing, the file is
        # scanned rather than skipped.
        if not path.is_file() or ".git" in path.parts:
            continue
        if "__pycache__" in path.parts:
            stem = path.name.split(".")[0]
            if (path.parent.parent / f"{stem}.py").is_file():
                continue        # derived from a source this scan already reads
            # orphaned cache with no source in scope: scan it
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        hits = [d for d in RETRACTED_CITATIONS if d in text]
        if hits:
            out[str(path.relative_to(REPO)).replace("\\", "/")] = sorted(hits)
    return out


def test_no_withdrawn_citation_survives_outside_the_recording_layer():
    carrying = files_carrying_retracted_citations()
    undeclared = {p: d for p, d in carrying.items()
                  if p not in CITATION_RECORDING_FILES}
    assert not undeclared, (
        "these files carry a DOI the ledger withdrew, and are not declared as "
        f"recording documents: {undeclared}. If this is the manuscript or the "
        "README, the retraction was never applied. If it is a new document that "
        "records the withdrawal, add it to CITATION_RECORDING_FILES deliberately.")


def test_the_replacements_the_ledger_prescribed_are_present():
    """The other half: absence of the wrong DOI is not presence of the right one.

    SCOPE: the manuscript only. A guard that checked only for the withdrawn string
    would pass against a bibliography with the entry simply deleted.
    """
    tex = PREPRINT.read_text(encoding="utf-8")
    for doi, row in (("10.1007/s11538-015-0117-1", "PP-REF-01"),
                     ("10.3390/math8010025", "PP-REF-02")):
        assert doi in tex, f"{row}'s prescribed replacement {doi} is not in the manuscript"


def test_the_retracted_citation_guard_detects_the_pre_fix_manuscript():
    """CONTROL: the manuscript as it stood before the fix must be caught.

    SCOPE: one file, recovered from git rather than written here. The known-bad is
    drawn from the artifact path -- a hand-written string would test the scanner
    against my idea of the citation, which is the control-that-never-met-the-
    pipeline defect AGENTS.md rule 1 records from the figure-sidecar case.
    """
    import subprocess
    before = subprocess.run(
        ["git", "show",
         "a8ba2b0:preprint/modeling_radioresistance_and_radiotropic_fitness.tex"],
        cwd=REPO, capture_output=True, text=True)
    if before.returncode != 0:
        pytest.skip("pre-fix manuscript not reachable from git")
    found = [d for d in RETRACTED_CITATIONS if d in before.stdout]
    assert sorted(found) == sorted(RETRACTED_CITATIONS), (
        "the pre-fix manuscript should contain both withdrawn DOIs; the scanner "
        f"found {found}. It is matching nothing, not finding nothing.")
    assert not [d for d in RETRACTED_CITATIONS
                if d in PREPRINT.read_text(encoding="utf-8")], \
        "and the current manuscript must contain neither"


# ------------------------------------------------------ retracted in sources
#
# CAUSE 1 OF THE UNAPPLIED-VERDICT PATTERN, AND THE LAST OF THE THREE TIERS.
# A ledger verdict is enforced only in the file its `document` column names. Two
# ways that leaves a claim unguarded, and they are different:
#
#   (a) the column names ANOTHER file. RM-G04-01 was restated in README while
#       biofilms_potts.jl:16 went on listing a five-term Hamiltonian, because no
#       guard reads source comments.
#   (b) the column names NO file. `document = repository` is a PSEUDO_DOCUMENT,
#       so document_path returns None and the row is never searched at all.
#       Measured: 37 rows sit in that class, 5 of them carrying an unresolved
#       verdict that names a source file. PP-DEFF-01 is one, and its prescribed
#       comment fix sat unapplied at biofilms_radiodialysis.R:242 until
#       2026-08-31 -- the fifth instance of the pattern found this session.
#
# WHY A DECLARED VOCABULARY RATHER THAN A PHRASE SWEEP, MEASURED BEFORE CHOOSING.
# Running distinguishing_phrase over every `delete` row against all 144 source
# files returns ZERO hits. The phrase tier cannot reach these claims: the ones
# that live in source comments are short -- a formula, a range, a symbol -- and
# split into runs under MIN_WORDS, exactly as figure labels and citations do.
# So this is the third vocabulary tier, beside RETRACTED_IN_FIGURES and
# RETRACTED_CITATIONS, for the same reason each of those exists.
#
# SEEDED FROM ACTUAL RETRACTIONS, NOT FROM INVENTED VOCABULARY, AND THE SEED IS
# SMALL ON PURPOSE. Two entries is what the ledger actually supports today. A
# guard whose value is that the NEXT one is caught does not need a long list, and
# a long list of invented strings would be the can't-fail shape in a new costume.
# `restate` verdicts are NOT swept generally -- AGENTS.md is explicit that a
# restate asks whether the revision says the right thing, which is a judgement
# and stays with a reviewer. Only the specific withdrawn STRING is mechanical.
RETRACTED_IN_SOURCES = {
    "Table 2 range 1e-4..1e-2":
        "PP-DEFF-01: the range is refuted at its upper end and unresolved at its "
        "lower, and Table 2 (tab:params) carries no D_eff row at all",
    "H_radiation + H_pairwise":
        "RM-G04-01: the Hamiltonian has four terms; total_pairwise_energy is real "
        "but is called from take_snapshot alone and never enters acceptance",
}

# The same use-versus-mention resolution as RETRACTED_CITATIONS, and it is needed
# for the same reason: recording a withdrawal requires naming it. EXPLICIT NAMES,
# NOT A GLOB -- a glob would auto-admit every future file matching a pattern,
# which is the property the allow-list exists to prevent.
SOURCE_RECORDING_FILES = {
    "data/claims_ledger.csv",
    "docs/research/external_reviews_2026-08-31_redteam.md",
    "calibration/tests/test_claims_ledger.py",
}

SOURCE_SUFFIXES = {".jl", ".R", ".py"}


def sources_carrying_retracted_strings() -> dict:
    """{repo-relative path: [strings]} for every source carrying a withdrawn string.

    SCOPE: files with a SOURCE_SUFFIXES extension anywhere in the repository,
    excluding `.git` and any file derived from a source already scanned. Shared by
    the check and its control.
    """
    out: dict[str, list[str]] = {}
    for path in REPO.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        if "__pycache__" in path.parts or path.suffix not in SOURCE_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        hits = [s for s in RETRACTED_IN_SOURCES if s in text]
        if hits:
            out[str(path.relative_to(REPO)).replace("\\", "/")] = sorted(hits)
    return out


def test_no_retracted_string_survives_in_a_source_file():
    carrying = sources_carrying_retracted_strings()
    undeclared = {p: s for p, s in carrying.items()
                  if p not in SOURCE_RECORDING_FILES}
    assert not undeclared, (
        "these sources carry a string the ledger withdrew, and are not declared "
        f"recording files: {undeclared}. The ledger holds the withdrawn wording; "
        "a correction comment must point at the row rather than repeat the "
        "string, or the source someone copies from still carries it.")


def test_the_source_retraction_guard_detects_the_pre_fix_file():
    """CONTROL: the R file as it stood before the fix must be caught.

    SCOPE: one file, recovered from git. Drawn from the artifact path rather than
    written here -- a hand-typed string would test the scanner against my idea of
    the comment, which is the control-that-never-met-the-pipeline defect.
    """
    import subprocess
    before = subprocess.run(
        ["git", "show", "59b7d95:biofilms_radiodialysis.R"],
        cwd=REPO, capture_output=True, text=True)
    if before.returncode != 0:
        pytest.skip("pre-fix source not reachable from git")
    found = [s for s in RETRACTED_IN_SOURCES if s in before.stdout]
    assert "Table 2 range 1e-4..1e-2" in found, (
        "the pre-fix R file should carry the withdrawn range; the scanner found "
        f"{found}. It is matching nothing, not finding nothing.")
    assert not [s for s in RETRACTED_IN_SOURCES
                if s in (REPO / "biofilms_radiodialysis.R").read_text(encoding="utf-8")], \
        "and the current file must carry none"
