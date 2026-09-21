"""Every code citation in the claims ledger resolves, and resolves to ONE place.

The ledger cited code as `file:line`. No test read a line number, and on 2026-09-21 only
49 of 153 still named the line the row was written about. A line number goes stale on any
insertion above it and fails nowhere.

A citation is now a file and a snippet:

    biofilms.R@`num_steps <- 1000`              the snippet occurs exactly once in the file
    README.md@9708a52@`four orders of`          exactly once in that file AT that commit,
                                                for a row about text that has since changed

Whitespace runs are flattened on both sides, so a snippet may span a line break and an
indentation change does not break it. Pattern from AGENTS.md rule 1: the control enters by
the production door, `citation_problems(root)`, on a scratch repository.

THREE OUTCOMES, NOT TWO. A snippet that occurs twice still "is present", and it points at
the wrong occurrence as easily as the right one. `ambiguous` is its own verdict.

SCOPE. This reads data/claims_ledger.csv alone. Line references into the superseded
preprint ("lines 711-721") name no file and are out of scope, as the ledger header says.
"""
import csv
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LEDGER = Path("data") / "claims_ledger.csv"

CITATION = re.compile(r"(?<![\w/.-])([\w./-]+\.[A-Za-z]\w*)@(?:([0-9a-f]{7,40})@)?`([^`]+)`")
# The retired form, with its `,:571` and ` and :24` tails.
BARE = re.compile(r"(?<![\w/.-])[\w./-]+\.[A-Za-z][A-Za-z0-9]*:\d+(?:[-–]\d+)?(?:(?:,|\s+and)\s*:\d+(?:[-–]\d+)?)*")

# A bare citation that may stay, keyed on the row AND the exact text. One entry: the deck was
# never tracked in this repository, and the note quotes the deck's sentence in full.
BARE_ALLOWED = {("HANDOUT-02", "wan-deck.tex:153-155")}

# A tripwire, not a target. A new citation changes it; update it in the same commit. It is
# here so that a pattern which silently stops matching cannot pass on zero citations.
N_CITATIONS = 152


def _flat(text: str) -> str:
    return " ".join(text.split())


def _rows(root: Path) -> list:
    with open(root / LEDGER, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def _text(root: Path, rel: str, commit):
    """The file's text today, or at a commit. None when it cannot be read."""
    if commit is None:
        path = root / rel
        return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else None
    shown = subprocess.run(["git", "-C", str(root), "show", f"{commit}:{rel}"],
                           capture_output=True, text=True)
    return shown.stdout if shown.returncode == 0 else None


def citations(root: Path = REPO) -> list:
    """(claim_id, column, file, commit or None, snippet) for every citation in the ledger."""
    return [(r["claim_id"], col, m.group(1), m.group(2), m.group(3))
            for r in _rows(root) for col, cell in r.items() if col
            for m in CITATION.finditer(cell or "")]


def resolve(root: Path, rel: str, commit, snippet: str) -> str:
    """`resolves`, `absent` or `ambiguous`."""
    text = _text(root, rel, commit)
    count = 0 if text is None else _flat(text).count(_flat(snippet))
    return "resolves" if count == 1 else "absent" if count == 0 else "ambiguous"


def citation_problems(root: Path = REPO) -> list:
    """One line per citation that does not resolve, and per bare `file:line` citation."""
    out = [f"{cid} [{col}] {rel}@`{snippet}`: {verdict}"
           for cid, col, rel, commit, snippet in citations(root)
           for verdict in [resolve(root, rel, commit, snippet)] if verdict != "resolves"]
    bare = {(r["claim_id"], m.group(0)) for r in _rows(root) for col, cell in r.items() if col
            for m in BARE.finditer(cell or "")}
    out += [f"{cid}: bare citation {lit}; write it as file@`snippet`" for cid, lit in sorted(bare - BARE_ALLOWED)]
    out += [f"{cid}: BARE_ALLOWED names {lit}, and the row no longer has it" for cid, lit in sorted(BARE_ALLOWED - bare)]
    return out


def test_every_citation_in_the_ledger_resolves_to_one_place():
    assert citation_problems() == []


def test_the_citation_count_is_the_pinned_one():
    found = citations()
    assert len(found) == N_CITATIONS, (
        f"{len(found)} citations, pinned {N_CITATIONS}. Read citation_problems() first, then "
        "update N_CITATIONS in the commit that adds or removes the citation.")


# --- controls, on a scratch repository ------------------------------------------------------

HEADER = "claim_id,code_location,notes\n"


def _scratch(tmp_path: Path, ledger_body: str, files: dict) -> Path:
    (tmp_path / "data").mkdir()
    (tmp_path / LEDGER).write_text("# a comment line\n" + HEADER + ledger_body, encoding="utf-8")
    for rel, text in files.items():
        (tmp_path / rel).write_text(text, encoding="utf-8")
    return tmp_path


STEPPER = "function step!(s)\n    s.t += 1\nend\n\nfunction other!(s)\n    s.u += 1\nend\n"


@pytest.mark.parametrize("label,cell,files,want", [
    ("the snippet occurs once", "stepper.jl@`s.t += 1`", {"stepper.jl": STEPPER}, []),
    ("it spans a line break", "stepper.jl@`function step!(s) s.t += 1`", {"stepper.jl": STEPPER}, []),
    ("a DUPLICATE of the cited line is planted", "stepper.jl@`s.t += 1`",
     {"stepper.jl": STEPPER + "function third!(s)\n    s.t += 1\nend\n"}, ["ambiguous"]),
    ("the cited line is removed", "stepper.jl@`s.t += 1`",
     {"stepper.jl": STEPPER.replace("s.t += 1", "s.t += 2")}, ["absent"]),
    ("the file is not there", "stepper.jl@`s.t += 1`", {}, ["absent"]),
    ("a bare line citation", "stepper.jl:2", {"stepper.jl": STEPPER}, ["bare citation stepper.jl:2"]),
    ("a bare citation with a tail", "stepper.jl:2,:6", {"stepper.jl": STEPPER}, ["bare citation stepper.jl:2,:6"]),
])
def test_the_checker_gives_each_verdict(tmp_path, monkeypatch, label, cell, files, want):
    monkeypatch.setitem(globals(), "BARE_ALLOWED", set())
    root = _scratch(tmp_path, f'X-01,"{cell}",\n', files)
    assert (root / LEDGER).read_text(encoding="utf-8").count(cell) == 1      # the plant is there
    got = citation_problems(root)
    assert len(got) == len(want) and all(w in g for w, g in zip(want, got)), (label, got)


def test_a_citation_at_a_commit_reads_that_commit(tmp_path, monkeypatch):
    """The text is gone today and present at the commit. Both directions are checked."""
    monkeypatch.setitem(globals(), "BARE_ALLOWED", set())
    root = _scratch(tmp_path, "", {"stepper.jl": STEPPER})

    def git(*args):
        done = subprocess.run(["git", "-C", str(root), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                              capture_output=True, text=True)
        if done.returncode:
            raise RuntimeError(done.stderr)          # plumbing, never the defect under test
        return done.stdout.strip()

    git("init", "-q")
    git("add", "stepper.jl")
    git("commit", "-q", "-m", "then")
    then = git("rev-parse", "--short=7", "HEAD")
    (root / "stepper.jl").write_text(STEPPER.replace("s.t += 1", "s.t += 2"), encoding="utf-8")
    (root / LEDGER).write_text(HEADER + f'X-01,"stepper.jl@{then}@`s.t += 1`",\n'
                                        f'X-02,"stepper.jl@`s.t += 1`",\n'
                                        f'X-03,"stepper.jl@{then}@`s.t += 2`",\n', encoding="utf-8")
    got = citation_problems(root)
    assert len(got) == 2 and "X-02" in got[0] and "absent" in got[0] and "X-03" in got[1], got


def test_a_stale_allowance_is_reported(tmp_path, monkeypatch):
    monkeypatch.setitem(globals(), "BARE_ALLOWED", {("X-01", "deck.tex:3")})
    root = _scratch(tmp_path, 'X-01,"deck.tex:3",\n', {})
    assert citation_problems(root) == []
    (root / LEDGER).write_text(HEADER + 'X-01,"no citation",\n', encoding="utf-8")
    assert citation_problems(root) == ["X-01: BARE_ALLOWED names deck.tex:3, and the row no longer has it"]
