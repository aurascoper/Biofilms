#!/usr/bin/env python3
"""Mechanical candidate list for the absence-sentence gate.

Grep proposes; a human disposes. The point is to stop reading one's own prose for
a shape one is primed to miss, so this errs toward flagging and expects a human
'no' on some of what it returns.

THIS SURFACES CANDIDATES AND DOES NOT PROVE COVERAGE. It matches surface forms
for a semantic property, and English has unbounded ways to say a thing is absent.
KNOWN_UNHANDLED below carries the forms it is known to miss, asserted NOT to flag
and printed as a standing gap count on every run, so a clean pass over a document
can never be read as an all-clear. That list exists because the first version
matched `not found` and missed `did not find`, and an absence went through
unflagged on the same diff that introduced the tool.
"""
import re, sys

HAZARD = [
    (r"\bno\b",            "bare negation"),
    (r"\bnone\b",          "bare negation"),
    (r"\bnobody\b",        "bare negation"),
    (r"\bnot found\b",     "explicit absence"),
    (r"\bdid not find\b", "explicit absence"),
    (r"\bfound none\b",   "explicit absence"),
    (r"\bno study\b|\bno such\b", "explicit absence"),
    (r"\bnever\b",         "universal"),
    (r"\balways\b",        "universal"),
    (r"\bevery\b",         "universal"),
    (r"\bonly\b",          "universal in POSITIVE grammar"),
    (r"\bunstudied\b",     "explicit absence"),
    (r"\babundant\b",      "comparative corpus claim"),
    (r"\bacross\b",        "possible venue scope"),
    (r"\bruns\b|\buses\b", "corpus habit stated as assertion"),
]

def scan(text):
    out = []
    for sent in re.split(r"(?<=[.?!])\s+", text.strip()):
        if not sent:
            continue
        hits = sorted({why for pat, why in HAZARD if re.search(pat, sent, re.I)})
        if hits:
            out.append((sent, hits))
    return out

# --- controls -------------------------------------------------------------
# SHOULD FLAG. Four are the known-bads from this planning round -- and they are a
# TRAINING SET, since the gate was derived from them. #5 is the test-set item: a
# universal wearing positive grammar, a shape none of the four has.
SHOULD_FLAG = [
    "Not found across ScienceDirect, OSTI, NACE and AMPP.",
    "Not found in a search on 30 August 2026 using the term set below.",
    "There is abundant MIC literature on 304L, 316 and high-entropy alloys.",
    "Screening runs unirradiated coupons in sterile electrolyte.",
    "The only work in this area is on chromia-passivated substrates.",   # test set
]

# SHOULD NOT FLAG. Ordinary page-2 prose, asserting nothing about existence or
# universality. A gate that flags everything scores perfectly on SHOULD_FLAG and
# is useless, which is the serial_seed42 over-count control in another costume.
SHOULD_PASS = [
    "Biofilm developed on both materials.",
    "They measured its radioactivity by gamma-ray spectrometry and found it had "
    "retained radionuclides from the pool water.",
    "One of those isolates is B. subtilis, which is already one of the seven "
    "species in the model.",
    "What dose the biofilm itself receives, and whether a hydrated layer 50 to 200 "
    "micrometres thick perturbs the local field at the metal surface.",
]

# KNOWN GAPS, KEPT VISIBLE. These are absence claims the matcher does NOT catch.
# They are asserted not to flag, so incompleteness prints on every run instead of
# living in a comment nobody reads. If one starts flagging, that is an UNEXPECTED
# PASS: a gap closed, and this list needs updating rather than celebrating.
KNOWN_UNHANDLED = [
    "This question has yet to be characterised for the alloy in question.",
    "The literature is silent on the behaviour of that surface.",
    "The topic remains uncharacterized at pool temperature.",
    "Work in this area is confined to chromia-passivated substrates.",
    "It is an open question whether the film changes.",
]


# ---------------------------------------------------------------------------
# SECOND SURFACE: TEST CONTRACTS THAT OVERSTATE THEIR ASSERTIONS
# ---------------------------------------------------------------------------
# Three defects found in one review round shared a signature: the guard bounded
# its own claim correctly while the DOCSTRING ABOVE IT overstated the reach. A
# sweep over one file under a ledger note claiming the repository; `match`
# reading one hit under a docstring claiming every location; a set intersection
# under a comment claiming the registry. That is not three bugs, it is one seam
# -- between what a check does and what its description claims -- and it is the
# seam this project's discipline is organised around, unguarded in the one place
# that describes the guards.
#
# THE CHECK IS NOT "NO UNIVERSAL QUANTIFIER". Universals fire legitimately here:
# "asserts every registered absence still holds" is accurate. The check is that a
# universal's SCOPE IS STATED, in a form a reader can compare against the
# assertion below it. So the gate does not infer scope -- it requires the
# docstring to declare one, because inferring it is the same guess that produced
# the defect.
#
# AND IT DECLARES ITS OWN SCOPE, or it reproduces the defect one level up: it
# reads the files it is given and claims nothing about any file it was not.

UNIVERSAL = re.compile(r"\b(every|all|each|any|no)\b", re.I)
SCOPE_LINE = re.compile(r"^\s*(#\s*)?SCOPE:\s*\S", re.M)


def _doc_blocks(text: str):
    """(line number, block) for each docstring or run of >=2 comment lines."""
    blocks, buf, start = [], [], None
    for i, line in enumerate(text.split("\n"), 1):
        st = line.strip()
        if st.startswith("#"):
            if start is None:
                start = i
            buf.append(st.lstrip("# ").rstrip())
        else:
            if buf and len(buf) >= 2:
                blocks.append((start, "\n".join(buf)))
            buf, start = [], None
    if buf and len(buf) >= 2:
        blocks.append((start, "\n".join(buf)))
    for m in re.finditer(r'"""(.*?)"""', text, re.S):
        blocks.append((text[:m.start()].count("\n") + 1, m.group(1)))
    return blocks


# A CONTRACT IS A BLOCK THAT DESCRIBES A CHECK, not any prose containing "every".
# The first version flagged 19 blocks across three files, most of them ordinary
# commentary -- which is the noise failure that gets a gate muted within a week.
# A block counts only if a test definition sits within DEF_WINDOW lines of it, so
# what is judged is the boundary between what a check does and what its
# description claims, which is the seam this exists for.
DEF_WINDOW = 4
TEST_DEF = re.compile(r"^\s*(def\s+test_|@testset|function\s+test_)", re.M)


def scan_contracts(text: str):
    """Blocks that describe a check and assert a universal without declaring a scope."""
    lines = text.split("\n")
    def_lines = {i for i, l in enumerate(lines, 1) if TEST_DEF.match(l)}
    out = []
    for lineno, block in _doc_blocks(text):
        if not (UNIVERSAL.search(block) and not SCOPE_LINE.search(block)):
            continue
        span = range(lineno - DEF_WINDOW, lineno + block.count("\n") + DEF_WINDOW + 1)
        if any(d in span for d in def_lines):
            out.append((lineno, block.strip()[:180]))
    return out


# Controls in DOCSTRING FORM, because the scanner reads comment runs and
# triple-quoted blocks -- a bare sentence exercises nothing, which the first
# version of these controls did and reported as a miss.
CONTRACT_MUST_FLAG = (
    "# Every measured requirement has an SOP or a named gap.\n"
    "# The check compares the index against the register.\n"
    "def test_coverage():\n")
CONTRACT_MUST_PASS = (
    "# SCOPE: the 14 measured rows of reference_d_requirements.csv, no other row.\n"
    "# Every measured requirement has an SOP or a named gap.\n"
    "def test_coverage():\n")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--contracts":
        # SCOPE: the files named on this command line, and nothing else.
        print(f"contract scope: {len(sys.argv) - 2} file(s) named on the command line")
        n = 0
        for path in sys.argv[2:]:
            for lineno, block in scan_contracts(open(path).read()):
                n += 1
                print(f"\nUNSCOPED CONTRACT  {path}:{lineno}\n     {block}")
        print(f"\n{n} universal(s) without a declared SCOPE. "
              "A clean run is not an all-clear: this reads only the files it was given.")
        sys.exit(0)

    if len(sys.argv) > 1:
        for sent, hits in scan(open(sys.argv[1]).read()):
            print(f"FLAG {hits}\n     {sent}\n")
        sys.exit(0)

    ok = True
    print("=== should flag (5th is the test-set item) ===")
    for s in SHOULD_FLAG:
        hits = scan(s)
        mark = "ok  " if hits else "MISS"
        ok &= bool(hits)
        print(f"  {mark} {s[:70]:<70} {hits[0][1] if hits else '-'}")
    print("\n=== should pass clean ===")
    for s in SHOULD_PASS:
        hits = scan(s)
        mark = "ok  " if not hits else "FALSE POSITIVE"
        ok &= not hits
        print(f"  {mark} {s[:70]:<70} {hits[0][1] if hits else ''}")
    print("\n=== known gaps: absence claims this matcher does NOT catch ===")
    surprises = 0
    for s_ in KNOWN_UNHANDLED:
        hits = scan(s_)
        if hits:
            surprises += 1
            print(f"  UNEXPECTED PASS  {s_[:66]:<66} {hits[0][1]}")
        else:
            print(f"  gap (expected)   {s_[:66]}")
    print(f"\n{len(KNOWN_UNHANDLED) - surprises} known gaps still open. "
          "A clean run over a document is NOT an all-clear.")
    if surprises:
        print(f"{surprises} entr{'y' if surprises == 1 else 'ies'} now flag -- "
              "a gap closed; move them out of KNOWN_UNHANDLED.")

    print("\n=== contract mode: a universal must declare a scope ===")
    c_bad = scan_contracts(CONTRACT_MUST_FLAG)
    c_ok = scan_contracts(CONTRACT_MUST_PASS)
    print(f"  {'ok  ' if c_bad else 'MISS'} universal without SCOPE flags")
    print(f"  {'ok  ' if not c_ok else 'FALSE POSITIVE'} the same universal WITH a SCOPE line passes")
    ok &= bool(c_bad) and not c_ok

    print("\nGATE USABLE" if ok else "\nGATE NOT USABLE -- fix before running it on new prose")
    sys.exit(0 if ok else 1)
