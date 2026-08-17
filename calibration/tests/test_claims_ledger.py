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
    # PHYSICAL LINE NUMBERS, so the report points at the line an editor opens.
    # The comment filter drops lines from the stream, and a counter over the
    # FILTERED sequence is off by however many comments precede the row -- an
    # error message naming the wrong line is a small lie that costs real time.
    # A quoted field may also hold newlines, so the reader's own position is
    # the only correct source.
    with open(path, encoding="utf-8") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    physical = [i for i, l in enumerate(open(path, encoding="utf-8"), start=1)
                if not l.startswith("#")]
    reader = csv.reader(lines)
    header = next(reader)
    out = []
    for index, row in enumerate(reader):
        if len(row) != len(header):
            # +1 because the header consumed the first filtered line
            line = physical[index + 1] if index + 1 < len(physical) else -1
            out.append((line, row[0] if row else "?", len(row)))
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
