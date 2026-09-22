"""Every code citation in the claims ledger resolves, and resolves to ONE place.

The ledger cited code as `file:line`. No test read a line number, and on 2026-09-21 only
41 of the 126 that a rule could place still named their line. A line number goes stale on any
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

WHAT THIS DOES NOT PROVE. A snippet that is unique and WRONG resolves. Codex found two on
pull request #49: a range whose first line was cited, where the row is about a later line.
That a citation names the right place is checked by reading the row beside the cited text,
and by nothing here.

READ BY EYE, 2026-09-21. Every citation of the conversion was read beside its row, the
audit-time text and today's line. 126 were placed by rule: 52 that a token comparison
flagged, with 7 wrong and fixed in pull request #49; and 74 unflagged, with 0 wrong (20 in
a seeded sample before the merge, then the other 54). A citation added later is read by its
author, and by nothing here.

LIVE MARKDOWN. Every tracked .md file whose first line is not `SUPERSEDED:` is read the same
way, with fenced blocks left out, and with two more verdicts:

- `revision missing`: a pinned commit that git cannot resolve, or that is not an ancestor of
  ANCESTOR_OF. A squash or a pruned branch can strand a pin that still reads fine locally.
- `context changed`: each citation also pins the digest of the element around it, so a
  snippet cannot stay true under a claim that was rewritten. The element is a table row, a
  list item, or else a blank-line block.

A dated record opens with `RECORD: line numbers refer to <commit>` and keeps its text; each
`file:line` in it must name a line that existed at that commit.

SCOPE. This reads data/claims_ledger.csv, its rows and its comment lines, and the live
Markdown. Line references into the superseded preprint ("lines 711-721") name no file and
are out of scope, as the ledger header says. A reference with no file name at all, such as
`.tex:580`, matches no pattern here and is read by nothing.
"""
import csv
import hashlib
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LEDGER = Path("data") / "claims_ledger.csv"

CITATION = re.compile(r"(?<![\w/.-])([\w./-]+\.[A-Za-z]\w*)@(?:([0-9a-f]{7,40})@)?`([^`]+)`")
# The retired form, with its `,:571` and ` and :24` tails.
BARE = re.compile(r"(?<![\w/.-])[\w./-]+\.[A-Za-z][A-Za-z0-9]*:\d+(?:[-–]\d+)?(?:(?:,|\s+and)\s*:\d+(?:[-–]\d+)?)*")

# A line number left behind a converted citation. Three forms were in the ledger after the
# first conversion: `snippet`,86,132 and "; :593 dH_mel" and "(:266)". Snippets and bare
# citations are taken out before this is applied, so neither is read twice. A journal's
# 34(6):433 has a bracket before its colon, and a time has a digit, so neither matches.
DETACHED = re.compile(r"(?<=[\s(;,]):\d+\b|`\s*,\s*\d+")

# A bare citation that may stay, keyed on the row AND the exact text. One entry: the deck was
# never tracked in this repository, and the note quotes the deck's sentence in full.
BARE_ALLOWED = {("HANDOUT-02", "wan-deck.tex:153-155")}

# A tripwire, not a target. A new citation changes it; update it in the same commit. It is
# here so that a pattern which silently stops matching cannot pass on zero citations.
N_CITATIONS = 160

# A pinned commit must be an ancestor of this ref. CI fetches with depth 0, so it is there.
ANCESTOR_OF = "origin/master"

# --- live Markdown ---------------------------------------------------------------------------
FENCE = re.compile(r"^(```|~~~).*?^\1[^\n]*$", re.M | re.S)
RECORD = re.compile(r"RECORD: line numbers refer to ([0-9a-f]{7,40})\s*$")
TAIL = re.compile(r":(\d+)(?:[-–](\d+))?")
# In Markdown only the colon form is a detached number: "`:1455`", "(:439, :450)". The
# backtick-comma form of the ledger gave 36 false alarms in prose.
MD_DETACHED = re.compile(r"(?<=[\s(;,`]):\d+\b")
LIST_ITEM = re.compile(r"(?:[-*+]|\d+[.)])\s")

# Bare citations that may stay in a live document or a record, keyed on the document AND the
# exact text. Both are unrecovered: the first names a file outside this repository, and the
# second matches no version of README.md at its record's commit or before it.
MD_BARE_ALLOWED = {("docs/metal_feasibility_macos_arm64.md", "JACC/src/threads/threads.jl:12-20"),
                   ("docs/audit_biofilms_potts.md", "README.md:101")}

# (document, citation, sha256 of the enclosing element's flattened text, first 16 hex digits).
# The element is the table row, the list item, or else the blank-line block around the
# citation. An edit to that element fails the gate until someone rereads the claim beside
# the cited text and updates one digest here. The digest reads raw text, markup and comments
# included, so a comment cannot hide words from it (issue #53 was that gap in another pin).
CONTEXT_PINS: set = set()


def _flat(text: str) -> str:
    return " ".join(text.split())


def _rows(root: Path) -> list:
    with open(root / LEDGER, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def _cells(root: Path) -> list:
    """(where, column, text) for every cell of every row AND for every comment line. The
    header holds citations too, and a reader that drops "#" lines never sees them."""
    with open(root / LEDGER, encoding="utf-8") as fh:
        comments = [(f"comment line {n}", "comment", line) for n, line in enumerate(fh, 1)
                    if line.startswith("#")]
    return comments + [(r["claim_id"], col, cell or "") for r in _rows(root)
                       for col, cell in r.items() if col]


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
    return [(where, col, m.group(1), m.group(2), m.group(3))
            for where, col, cell in _cells(root) for m in CITATION.finditer(cell)]


def _revision_ok(root: Path, commit: str) -> bool:
    """The commit exists AND is an ancestor of ANCESTOR_OF. A commit on a branch that was
    squashed or pruned can still be in a local object store, and would read fine there."""
    git = lambda *a: subprocess.run(["git", "-C", str(root), *a], capture_output=True, text=True).returncode
    if git("rev-parse", "--verify", "--quiet", f"{ANCESTOR_OF}^{{commit}}"):
        raise RuntimeError(f"{ANCESTOR_OF} is not a ref here; fetch it")    # plumbing, not a verdict
    return git("cat-file", "-e", f"{commit}^{{commit}}") == 0 and \
        git("merge-base", "--is-ancestor", commit, ANCESTOR_OF) == 0


def resolve(root: Path, rel: str, commit, snippet: str) -> str:
    """`resolves`, `absent`, `ambiguous` or `revision missing`. The fourth is its own verdict
    so that nobody "fixes" a pin that lost its commit by repointing it at today's text."""
    if commit is not None and not _revision_ok(root, commit):
        return "revision missing"
    text = _text(root, rel, commit)
    count = 0 if text is None else _flat(text).count(_flat(snippet))
    return "resolves" if count == 1 else "absent" if count == 0 else "ambiguous"


def citation_problems(root: Path = REPO) -> list:
    """One line per citation that does not resolve, per bare `file:line` citation and per
    detached line number, in the rows and in the comment lines."""
    out = [f"{cid} [{col}] {rel}@`{snippet}`: {verdict}"
           for cid, col, rel, commit, snippet in citations(root)
           for verdict in [resolve(root, rel, commit, snippet)] if verdict != "resolves"]
    bare = {(where, m.group(0)) for where, _col, cell in _cells(root) for m in BARE.finditer(cell)}
    out += [f"{cid}: bare citation {lit}; write it as file@`snippet`" for cid, lit in sorted(bare - BARE_ALLOWED)]
    out += [f"{where}: detached line number {m.group(0)!r}; give it a file and a snippet"
            for where, _col, cell in _cells(root)
            for m in DETACHED.finditer(BARE.sub(" ", re.sub(r"`[^`]*`", "`", cell)))]
    out += [f"{cid}: BARE_ALLOWED names {lit}, and the row no longer has it" for cid, lit in sorted(BARE_ALLOWED - bare)]
    return out


def _md_docs(root: Path) -> list:
    """(path, text) for every tracked .md file whose first line is not `SUPERSEDED:`."""
    listed = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "*.md"], capture_output=True, text=True)
    if listed.returncode:
        raise RuntimeError(listed.stderr)                   # plumbing, not a verdict
    docs = [(rel, (root / rel).read_text(encoding="utf-8", errors="replace"))
            for rel in listed.stdout.split("\0") if rel and (root / rel).is_file()]
    return [(rel, text) for rel, text in docs if not text.startswith("SUPERSEDED:")]


def _elements(text: str) -> list:
    """(first line number, raw text) per element: a table row, a list item with its
    continuation lines, a heading, or else a blank-line block. Fenced blocks are blanked
    first, keeping the line numbers. A whole table or list as one element would let an edit
    to one row void the context of every citation in the others."""
    text = FENCE.sub(lambda m: "\n" * m.group(0).count("\n"), text)
    out, cur, start = [], [], 0
    for n, line in enumerate(text.split("\n"), 1):
        s = line.lstrip()
        new = s.startswith(("|", "#")) or bool(LIST_ITEM.match(s))
        if cur and (not s or new):
            out.append((start, "\n".join(cur)))
            cur = []
        if s:
            start = n if not cur else start
            cur.append(line)
            if s.startswith(("|", "#")):
                out.append((start, line))
                cur = []
    if cur:
        out.append((start, "\n".join(cur)))
    return out


def context_digest(element: str) -> str:
    return hashlib.sha256(_flat(element).encode("utf-8")).hexdigest()[:16]


def _record_problems(root: Path, doc: str, commit: str, lit: str) -> list:
    """A record's `file:line` must name a line that existed at the record's commit. The path
    is read as written, then beside the record, then as the one file with that name."""
    name = lit.split(":", 1)[0]
    listed = subprocess.run(["git", "-C", str(root), "ls-tree", "-r", "--name-only", commit],
                            capture_output=True, text=True).stdout.split("\n")
    beside = str(Path(doc).parent / name)
    same = [p for p in listed if p.rsplit("/", 1)[-1] == name.rsplit("/", 1)[-1]]
    path = name if name in listed else beside if beside in listed else same[0] if len(same) == 1 else None
    if path is None:
        return [f"{doc}: {lit} names no file at {commit}"]
    n_lines = len((_text(root, path, commit) or "").split("\n"))
    ends = [int(t.group(2) or t.group(1)) for t in TAIL.finditer(lit[len(name):])]
    return [f"{doc}: {lit} is past the end of {path} at {commit} ({n_lines} lines)" for e in ends if e > n_lines][:1]


def markdown_problems(root: Path = REPO) -> list:
    """One line per problem in the live Markdown: a citation that does not resolve, an
    element whose text changed since its digest was pinned, a citation with no pin, a bare
    `file:line`, a detached number, and a stale pin or allowance. A document whose first line
    is `RECORD: line numbers refer to <commit>` keeps its bare citations, and each must
    name a line that existed at that commit."""
    out, seen, allowed = [], set(), set()
    for doc, text in _md_docs(root):
        record = RECORD.match(text.split("\n", 1)[0])
        record_ok = record is not None and _revision_ok(root, record.group(1))
        if record and not record_ok:
            out.append(f"{doc}: RECORD names {record.group(1)}: revision missing")
        for n, el in _elements(text):
            where = f"{doc}:{n}"
            for m in CITATION.finditer(el):
                verdict = resolve(root, m.group(1), m.group(2), m.group(3))
                if verdict != "resolves":
                    out.append(f"{where} {m.group(0)}: {verdict}")
                if record is None:
                    key = (doc, m.group(0), context_digest(el))
                    seen.add(key[:2])
                    if key not in CONTEXT_PINS:
                        known = any(p[:2] == key[:2] for p in CONTEXT_PINS)
                        out.append(f"{where} {m.group(0)}: " + (
                            f"context changed, digest now {key[2]}; reread the claim, then update its pin" if known
                            else f"no context pin; read the claim beside the cited text, then pin {key[2]}"))
            rest = CITATION.sub(" ", el)
            for m in BARE.finditer(rest):
                if (doc, m.group(0)) in MD_BARE_ALLOWED:
                    allowed.add((doc, m.group(0)))
                elif record:
                    out += _record_problems(root, doc, record.group(1), m.group(0)) if record_ok else []
                else:
                    out.append(f"{where}: bare citation {m.group(0)}; write it as file@`snippet`")
            if record is None:
                out += [f"{where}: detached line number {m.group(0)!r}; give it a file and a snippet"
                        for m in MD_DETACHED.finditer(BARE.sub(" ", rest))]
    out += [f"CONTEXT_PINS names {lit} in {doc}, and the document no longer has it"
            for doc, lit, _ in sorted(CONTEXT_PINS) if (doc, lit) not in seen]
    out += [f"MD_BARE_ALLOWED names {lit} in {doc}, and the document no longer has it"
            for doc, lit in sorted(MD_BARE_ALLOWED - allowed)]
    return out


@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="evidence: the live Markdown still cites by line number; the next commit converts it")
def test_every_citation_in_live_markdown_resolves_and_keeps_its_context():
    assert markdown_problems() == []


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
    # The three forms below are the ones the first conversion left in the real ledger.
    ("numbers left after the snippet", "stepper.jl@`s.t += 1`,6", {"stepper.jl": STEPPER}, ["detached line number"]),
    ("a tail after prose", "stepper.jl@`s.t += 1` T_cpm; :6 acceptance", {"stepper.jl": STEPPER},
     ["detached line number ':6'"]),
    ("a tail in brackets", "stepper.jl@`s.t += 1` and the other (:6)", {"stepper.jl": STEPPER},
     ["detached line number ':6'"]),
    ("a journal volume and page is not a line number", "stepper.jl@`s.t += 1`; J Ind Microbiol 34(6):433-441, at 16:48",
     {"stepper.jl": STEPPER}, []),
])
def test_the_checker_gives_each_verdict(tmp_path, monkeypatch, label, cell, files, want):
    monkeypatch.setitem(globals(), "BARE_ALLOWED", set())
    root = _scratch(tmp_path, f'X-01,"{cell}",\n', files)
    assert (root / LEDGER).read_text(encoding="utf-8").count(cell) == 1      # the plant is there
    got = citation_problems(root)
    assert len(got) == len(want) and all(w in g for w, g in zip(want, got)), (label, got)


# THE LEDGER'S COMMENT LINES. Raised by Codex on pull request #49 against 641f278, and
# reproduced first: header line 176 cited a file that is gone, as `file:157-182`, and the
# gate read data rows alone. The first form below is that header line's own.
@pytest.mark.parametrize("comment,want", [
    ("#   Counted directly at stepper.jl:2-6: two rows.", "bare citation stepper.jl:2-6"),
    ("#   The update is stepper.jl@`s.t += 9`.", "absent"),
    ("#   ...and the second update, at (:6).", "detached line number"),
])
def test_a_comment_line_of_the_ledger_is_read_too(tmp_path, monkeypatch, comment, want):
    monkeypatch.setitem(globals(), "BARE_ALLOWED", set())
    root = _scratch(tmp_path, 'X-01,"stepper.jl@`s.t += 1`",\n', {"stepper.jl": STEPPER})
    path = root / LEDGER
    path.write_text(comment + "\n" + path.read_text(encoding="utf-8"), encoding="utf-8")
    if path.read_text(encoding="utf-8").count(comment) != 1:
        raise RuntimeError("the plant did not apply")         # plumbing, never the defect
    got = citation_problems(root)
    assert len(got) == 1 and want in got[0], got


def test_a_citation_at_a_commit_reads_that_commit(tmp_path, monkeypatch):
    """The text is gone today and present at the commit. Both directions are checked."""
    monkeypatch.setitem(globals(), "BARE_ALLOWED", set())
    monkeypatch.setitem(globals(), "ANCESTOR_OF", "HEAD")
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


# --- controls for the Markdown scope, on a scratch repository ----------------------------------

def _git(root: Path, *args) -> str:
    done = subprocess.run(["git", "-C", str(root), "-c", "user.name=t", "-c", "user.email=t@t", *args],
                          capture_output=True, text=True)
    if done.returncode:
        raise RuntimeError(done.stderr)                  # plumbing, never the defect under test
    return done.stdout.strip()


def _md_scratch(tmp_path: Path, monkeypatch, files: dict) -> Path:
    """A git repository holding `files`, committed, with ANCESTOR_OF at its HEAD."""
    for rel, text in files.items():
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text(text, encoding="utf-8")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "scratch")
    monkeypatch.setitem(globals(), "ANCESTOR_OF", "HEAD")
    monkeypatch.setitem(globals(), "MD_BARE_ALLOWED", set())
    monkeypatch.setitem(globals(), "CONTEXT_PINS", set())
    return tmp_path


CLAIM = "The step adds one: stepper.jl@`s.t += 1`."
TABLE = "| step | stepper.jl@`s.t += 1` |\n| other | stepper.jl@`s.u += 1` |\n"


@pytest.mark.parametrize("label,doc,want", [
    ("a bare line citation", "The step is at `stepper.jl:2`.\n", ["bare citation stepper.jl:2"]),
    ("a tail left in its own code span", "The step is at stepper.jl@`s.t += 1` and `:6`.\n", ["detached line number ':6'"]),
    ("a superseded document is not read", "SUPERSEDED: 2026-09-21 — an old plan\n\n`stepper.jl:2`\n", []),
    ("a fenced block is not read", "Before.\n\n```\nstepper.jl:2\n```\n", []),
    ("a citation with no context pin", CLAIM + "\n", ["no context pin"]),
    ("a pinned commit that does not resolve", "Then: stepper.jl@0123456@`s.t += 1`.\n", ["revision missing"]),
])
def test_the_markdown_checker_gives_each_verdict(tmp_path, monkeypatch, label, doc, want):
    root = _md_scratch(tmp_path, monkeypatch, {"stepper.jl": STEPPER, "notes.md": doc})
    got = [g for g in markdown_problems(root) if "no context pin" not in g or "no context pin" in want]
    assert len(got) == len(want) and all(w in g for w, g in zip(want, got)), (label, got)


def test_a_pinned_element_edited_fails_and_its_neighbour_does_not(tmp_path, monkeypatch):
    """The plan's plant: an edit inside a cited table row must fail, and an edit to the next
    row must not, which is why a row and not the table is the element."""
    root = _md_scratch(tmp_path, monkeypatch, {"stepper.jl": STEPPER, "notes.md": TABLE})
    rows = TABLE.splitlines()
    monkeypatch.setitem(globals(), "CONTEXT_PINS", {("notes.md", "stepper.jl@`s.t += 1`", context_digest(rows[0])),
                                                    ("notes.md", "stepper.jl@`s.u += 1`", context_digest(rows[1]))})
    assert markdown_problems(root) == []
    for n, word in ((1, "another"), (0, "the step")):
        edited = TABLE.replace(rows[n], rows[n].replace("| other |" if n else "| step |", f"| {word} |"))
        if edited == TABLE:
            raise RuntimeError("the plant did not apply")         # plumbing, never the defect
        (root / "notes.md").write_text(edited, encoding="utf-8")
        got = markdown_problems(root)
        cited = "s.u += 1" if n else "s.t += 1"
        assert len(got) == 1 and "context changed" in got[0] and cited in got[0], (n, got)


def test_a_pin_on_a_commit_off_the_trunk_is_a_missing_revision(tmp_path, monkeypatch):
    """The commit exists in the object store, on a branch the trunk never merged."""
    root = _md_scratch(tmp_path, monkeypatch, {"stepper.jl": STEPPER, "notes.md": "none\n"})
    _git(root, "checkout", "-q", "-b", "side")
    (root / "stepper.jl").write_text(STEPPER + "# side\n", encoding="utf-8")
    _git(root, "commit", "-q", "-am", "side")
    side = _git(root, "rev-parse", "--short=7", "HEAD")
    _git(root, "checkout", "-q", "-")
    monkeypatch.setitem(globals(), "ANCESTOR_OF", _git(root, "rev-parse", "HEAD"))
    (root / "notes.md").write_text(f"Then: stepper.jl@{side}@`# side`.\n", encoding="utf-8")
    got = markdown_problems(root)
    assert len(got) == 2 and "revision missing" in got[0] and "no context pin" in got[1], got


def test_a_record_keeps_its_line_numbers_and_each_must_have_existed(tmp_path, monkeypatch):
    root = _md_scratch(tmp_path, monkeypatch, {"stepper.jl": STEPPER, "notes.md": "none\n"})
    then = _git(root, "rev-parse", "--short=7", "HEAD")
    for body, want in ((f"RECORD: line numbers refer to {then}\n\nAt `stepper.jl:2` and `:6`.\n", []),
                       (f"RECORD: line numbers refer to {then}\n\nAt `stepper.jl:99`.\n", ["past the end"]),
                       ("RECORD: line numbers refer to 0123456\n\nAt `stepper.jl:2`.\n", ["revision missing"])):
        (root / "notes.md").write_text(body, encoding="utf-8")
        got = markdown_problems(root)
        assert len(got) == len(want) and all(w in g for w, g in zip(want, got)), (body, got)


def test_a_stale_pin_and_a_stale_allowance_are_reported(tmp_path, monkeypatch):
    root = _md_scratch(tmp_path, monkeypatch, {"stepper.jl": STEPPER, "notes.md": "no citation\n"})
    monkeypatch.setitem(globals(), "CONTEXT_PINS", {("notes.md", "stepper.jl@`s.t += 1`", "0" * 16)})
    monkeypatch.setitem(globals(), "MD_BARE_ALLOWED", {("notes.md", "deck.tex:3")})
    got = markdown_problems(root)
    assert len(got) == 2 and "CONTEXT_PINS names" in got[0] and "MD_BARE_ALLOWED names" in got[1], got
