"""The withdrawn melanin/β_ion comparison must not survive anywhere in the tree.

WHY A WHOLE-TREE SCAN AND NOT A GREP. RM-PROP-01 recorded that the v1.2
correction of this comparison reached six working branches and reached neither
master nor five other remotes, and it stated its own method honestly: the
withdrawn phrase, in README.md, over every remote branch. The bound it declared
was real. The bound it did NOT declare is the one that bit: the sweep matched
LINE BY LINE, and README wraps at the column, so `four orders of` ended one line
and `magnitude more` began the next. A third instance sat three hundred lines
above the two the sweep fixed, in the same file it had just called clean, and an
external review (Codex, PR #19) found it along with four more files the
one-file scope never reached.

So: normalise whitespace before looking, and look in every tracked text file.

WHAT THIS GUARD CONSIDERS A DEFECT, and it is narrower than the phrase. "Four
orders of magnitude" is legitimate prose in this repository -- REFINE-08 uses it
for sub-voxel tally drift, and the dose audits in docs/research use it for
Bland 2022 sitting that far below every other ionizing record. The withdrawn claim is
the COMPARISON: the melanin acceptance bias set against the radiation term and
the gap called four orders. An occurrence counts only when both subjects are in
the window, and it is forgiven when the window also says the comparison is
withdrawn -- which is how the corrected sites, and this docstring, stay legal.

WHAT STILL BOUNDS IT. A restatement carrying the superseded picture in other
words -- "dwarfs", "negligible beside" -- with neither magnitude spelled out is
invisible here, exactly as it was to the grep. The term-set sweep in PP-62-13
covers that vocabulary for the manuscript; nothing covers it tree-wide.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

SUFFIXES = {".md", ".jl", ".csv", ".py", ".tex", ".txt", ".yml", ".yaml",
            ".json", ".R", ".r", ".toml", ".sh"}

# The withdrawn magnitude, either as the ratio in words or as the number the
# comparison against 1.000010 produces.
WITHDRAWN_MAGNITUDE = re.compile(r"(?:four|4) orders of magnitude|15,?500\s*(?:x|×|-?fold)")
MELANIN_TERM = re.compile(r"melanin|[dδ]h_mel|h_mel")
RADIATION_TERM = re.compile(r"beta_ion|β_ion|[dδ]h_rad|h_rad|1\.000010|one part in 1")
# Naming the withdrawal is the point of the corrected sites.
DECLARED = re.compile(r"withdrawn|withdraw|retracted")

WINDOW = 300


# A wrapped phrase inside a comment block does not merely break across a
# newline -- the next line's COMMENT MARKER lands in the middle of it. The
# compute_delta_H comment in biofilms_potts.jl read "Four orders of\n# magnitude"
# and flattening whitespace alone produced "four orders of # magnitude", which
# matches nothing. Collapsing whitespace was necessary and was not sufficient,
# and this guard shipped with that hole in it: the first whole-tree scan over
# every remote branch reported zero for a file that carried the claim on all of
# them. Strip the marker first.
COMMENT_PREFIX = re.compile(r"(?m)^[ \t]*(?:#+|//+|%+|;+|--)[ \t]?")


def flatten(text: str) -> str:
    """Comment markers dropped, whitespace collapsed, lowercased. THE WHOLE
    POINT: a phrase that wraps at the column is one phrase, and neither a
    line-oriented search nor a naive whitespace collapse can see it."""
    return " ".join(COMMENT_PREFIX.sub("", text).split()).lower()


def survivors_in(text: str) -> list[str]:
    flat = flatten(text)
    out = []
    for match in WITHDRAWN_MAGNITUDE.finditer(flat):
        window = flat[max(0, match.start() - WINDOW):match.end() + WINDOW]
        if not (MELANIN_TERM.search(window) and RADIATION_TERM.search(window)):
            continue                       # some other four orders of magnitude
        if DECLARED.search(window):
            continue                       # the correction, naming what it replaced
        out.append(window)
    return out


def tracked_text_files() -> list[Path]:
    listing = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z"],
                             capture_output=True, text=True, check=True).stdout
    return [REPO / name for name in listing.split("\0")
            if name and Path(name).suffix in SUFFIXES]


def test_the_withdrawn_comparison_survives_nowhere_in_the_tree():
    """THE ONE THAT MATTERS."""
    survivors = []
    for path in tracked_text_files():
        if path == Path(__file__):
            continue
        for window in survivors_in(path.read_text(encoding="utf-8", errors="replace")):
            rel = path.relative_to(REPO)
            survivors.append(f"{rel}: ...{window[WINDOW - 90:WINDOW + 120]}...")
    assert not survivors, (
        "the withdrawn melanin/β_ion comparison is asserted without being named "
        "as withdrawn:\n  " + "\n  ".join(survivors))


def test_the_scan_detects_the_instance_the_grep_missed():
    """THE CONTROL, AND IT IS THE REAL SENTENCE.

    README.md as it stood at c31917c, wrapped where README wrapped it. A
    line-oriented search for the withdrawn phrase returns nothing here; this
    scan must return one.
    """
    as_it_wrapped = (
        "`β_ion` — the one parameter Table 2 tabulates per species, and the one the entire sign convention is\n"
        "written around — biases acceptance by **one part in 10⁵** for exactly the species whose radial\n"
        "stratification is the headline result. The melanin term biases it by **15.5%**: four orders of\n"
        "magnitude more, through a coefficient of `0.5` hard-coded at its call site in `compute_delta_H`,\n"
        "appearing in no table and in no configuration file.\n")
    assert "four orders of magnitude" not in as_it_wrapped, (
        "the control no longer wraps, so it no longer reproduces the defect")
    assert not [l for l in as_it_wrapped.splitlines()
                if "four orders of magnitude" in l], "same, line by line"
    assert len(survivors_in(as_it_wrapped)) == 1


def test_the_scan_detects_the_instance_the_comment_marker_hid():
    """THE SECOND CONTROL, AND THE ONE THIS GUARD FIRST FAILED.

    compute_delta_H in biofilms_potts.jl as it stood on origin/master. The
    phrase wraps AND the next line opens with a comment marker, so flattening
    whitespace alone yields "four orders of # magnitude" and matches nothing.
    The first version of this scan reported this file clean on every remote
    branch that carried the claim -- a whole-tree walk is not a wider scan if
    its normaliser is narrower than the text.
    """
    as_it_was = (
        "    # T_cpm = 5.0, the radiation term for a radiotropic species is\n"
        "    # \u03b2_ion\u00b7I = -5e-5, an acceptance bias of 1.000010 \u2014 one part in 1e5. This\n"
        "    # term at the reported M = 1.44 is -0.72, a bias of 1.155. Four orders of\n"
        "    # magnitude. The radial stratification is therefore melanin-mediated, not\n"
        "    # \u03b2_ion-mediated.\n")
    assert not [l for l in as_it_was.splitlines()
                if "four orders of magnitude" in l.lower()], "line by line: invisible"
    assert "four orders of magnitude" not in " ".join(as_it_was.split()).lower(), (
        "whitespace collapse alone: still invisible, because the marker is in the way")
    assert len(survivors_in(as_it_was)) == 1


def test_the_scan_passes_the_correction_and_the_unrelated_uses():
    """THE OTHER HALF OF THE CONTROL. A guard that flags everything is no more
    use than one that flags nothing, and this one is aimed at a comparison, not
    at a phrase."""
    corrected = ("The melanin term biases it by 15.5%, about an order of magnitude more than "
                 "`β_ion` at the term's reach. Comparing it against the 1.000010 instead and "
                 "calling the gap four orders of magnitude is the withdrawn comparison.")
    unrelated = ("Float accumulation would predict the 8-term and 64-term sums to differ by "
                 "~8x, not by four orders of magnitude, so it is not accumulation.")
    dose = ("the total dose is 0.5-1 Gy over seven days, three to four orders of magnitude "
            "below every other ionizing record here. On melanin: pigmentation INCREASED "
            "under UV and DECREASED under gamma.")
    assert survivors_in(corrected) == []
    assert survivors_in(unrelated) == []
    assert survivors_in(dose) == []


def test_the_scan_reaches_more_than_one_file_type():
    """The sweep this replaces was bounded to README.md, and four of the five
    surfaces the review found were not README."""
    suffixes = {p.suffix for p in tracked_text_files()}
    assert {".md", ".jl", ".csv"} <= suffixes
    assert len(tracked_text_files()) > 100
