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


def _preprint_text() -> str:
    """The manuscript as one normalised lowercase line.

    LaTeX wraps at the column, so a claim that survives revision is usually
    split across a newline. Collapsing whitespace is what makes the search find
    it; without this the test passes for the wrong reason.
    """
    return " ".join(PREPRINT.read_text(encoding="utf-8").split()).lower()


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


def test_the_uncheckable_claims_are_named_not_hidden(rows, capsys):
    """Coverage is reported because it is partial. A reader who believes this
    suite covers all 34 would retire the manual review that covers the rest."""
    deleted = [r for r in rows
               if r["claim_id"].startswith("PP-") and r["status"] == "delete"]
    uncheckable = [r["claim_id"] for r in deleted
                   if not distinguishing_phrase(r["claim_text"])]

    with capsys.disabled():
        print(f"\n  claims-ledger guard: {len(deleted) - len(uncheckable)} of "
              f"{len(deleted)} `delete` claims are textually checkable")
        if uncheckable:
            print("  NOT textually checkable — these need human review:")
            for cid in uncheckable:
                print(f"      {cid}")

    assert len(deleted) - len(uncheckable) >= 25, (
        "coverage collapsed: most `delete` claims no longer yield a searchable "
        "phrase, so this guard has quietly stopped guarding")


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

    text = PREPRINT.read_text(encoding="utf-8")
    cited = set(re.findall(r"repository revision \\texttt\{([0-9a-f]{7,40})\}", text))
    cited |= set(re.findall(r"manuscript is \\texttt\{([0-9a-f]{7,40})\}", text))
    assert cited, "the manuscript cites no repository revision"
    for sha in cited:
        done = subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{sha}^{{commit}}"],
                              capture_output=True)
        assert done.returncode == 0, f"cited revision {sha} does not resolve"
