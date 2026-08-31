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
# reads the files it is given and claims nothing about any file it was not,
# and it reports TWO TIERS because a proximity window alone has a
# false-negative mode -- see scan_contracts. As of 2026-08-31 the attached
# tier is triaged and the module-level tier is a BACKLOG, not a judged set:
# 17 blocks across six test files remain unscoped and unjudged. Saying so is
# the difference between a backlog and a silent gap.

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


# THE WINDOW'S TWO ERROR MODES ARE NOT SYMMETRIC, and the second is the one this
# tool exists to catch. A false positive costs a judgement call, which is cheap
# and gets recorded. A false NEGATIVE is a contract the gate silently does not
# cover -- and enumerating the drops showed the set is NOT empty: a module-level
# block stating "ALL FOUR BRANCHES MUST BE EXERCISED", and one stating which
# run_by values are exempt from a check, both sit further than DEF_WINDOW from
# any test definition. A gate claiming to guard the check-versus-description
# boundary while missing module-level contracts overstates its own reach, which
# is the subset defect inside the tool built to detect it.
#
# So there are two tiers and the caller is told which is which. Attached
# contracts are near a test definition. Module-level contracts are further away
# and are admitted only when they speak NORMATIVELY -- must, never, cannot,
# only, exempt -- because that is what separates a rule about a check from
# commentary about the physics.
NORMATIVE = re.compile(r"\b(must|never|cannot|only|exempt|required)\b", re.I)


def scan_contracts(text: str, *, module_level: bool = False):
    """Blocks that describe a check and assert a universal without declaring a scope.

    SCOPE: the text given. `module_level=False` returns contracts within
    DEF_WINDOW lines of a test definition; `module_level=True` returns the
    normatively-worded blocks outside that window. Neither set is the other's
    superset and both are reported.
    """
    lines = text.split("\n")
    def_lines = {i for i, l in enumerate(lines, 1) if TEST_DEF.match(l)}
    out = []
    for lineno, block in _doc_blocks(text):
        if not (UNIVERSAL.search(block) and not SCOPE_LINE.search(block)):
            continue
        # ASSOCIATE WITH THE ENCLOSING TEST, NOT ONLY A NEARBY DEFINITION. The
        # proximity window missed the defect that motivated this gate: the
        # "EVERY LOCATION" contract in prose_bounds.jl sits inside a @testset
        # whose header is 25 lines above it, so no definition was within four
        # lines and attached-mode returned nothing for the very file it was
        # built from.
        #
        # BUT "IS THERE A DEFINITION ABOVE ME" IS NOT THE SAME QUESTION AS "AM I
        # INSIDE ONE", AND USING IT MADE THE MODULE-LEVEL TIER UNREACHABLE. Once
        # a file's first test definition is passed, every later block satisfied
        # `preceding` -- so a block that had returned to module scope was
        # reported as attached, and module-level mode returned nothing for any
        # file with a test in it. Measured on calibration/tests/test_sop_index.py:
        # 8 attached, 0 module-level, while the no-signal registry at line 164
        # and three further module blocks sat there being misfiled.
        #
        # THAT IS RULE 3 INSIDE THE TOOL BUILT FOR RULE 3. The header advertises
        # the module-level count as a standing backlog printed every run so a
        # clean pass is never an all-clear; it printed 0, and an unreachable
        # tier and an empty tier print the same 0. Raised by Codex on pull
        # request #23.
        #
        # Membership is now INDENTATION, which is what "inside a body" actually
        # means in both languages this reads: a block at column 0 is module
        # scope, an indented block belongs to the definition enclosing it.
        raw = lines[lineno - 1] if lineno - 1 < len(lines) else ""
        indent = len(raw) - len(raw.lstrip())
        preceding = [d for d in def_lines if d <= lineno]
        attached = indent > 0 and bool(preceding)
        if module_level:
            if indent == 0 and def_lines and NORMATIVE.search(block):
                out.append((lineno, block.strip()[:180]))
        elif attached:
            out.append((lineno, block.strip()[:180]))
    return out


# Controls in DOCSTRING FORM, because the scanner reads comment runs and
# triple-quoted blocks -- a bare sentence exercises nothing, which the first
# version of these controls did and reported as a miss.
CONTRACT_MUST_FLAG = (
    "def test_coverage():\n"
    "    # Every measured requirement has an SOP or a named gap.\n"
    "    # The check compares the index against the register.\n")
# IN-BODY CONTROL: the motivating defect's own shape -- a contract documented
# beside assertions well inside a test body, far from its header. The proximity
# window returned nothing for this and reported the file clean.
CONTRACT_IN_BODY_MUST_FLAG = (
    "@testset \"ratios\" begin\n"
    "    x = 1\n    y = 2\n    z = 3\n    w = 4\n    v = 5\n"
    "    # EVERY LOCATION, NOT THE FIRST MATCH. `match` returns one hit.\n"
    "    # The docstring claimed every location and covered one.\n"
    "    @test x == 1\nend\n")
CONTRACT_MUST_PASS = (
    "def test_coverage():\n"
    "    # SCOPE: the 14 measured rows of reference_d_requirements.csv, no other row.\n"
    "    # Every measured requirement has an SOP or a named gap.\n")
# THE SHAPE THAT WAS INVISIBLE: a module-level contract sitting AFTER a completed
# test. Under the "is there a def above me" rule this was filed as attached, so
# the module-level tier reported nothing and its backlog count read as zero. It
# must appear in module-level mode and must NOT appear in attached mode -- both
# directions, because classifying it into the wrong tier is the actual defect and
# a one-sided check cannot see it.
CONTRACT_MODULE_AFTER_TEST_MUST_FLAG = (
    "def test_one():\n"
    "    assert True\n"
    "\n"
    "# EVERY registered absence must back an index row.\n"
    "# Nothing here declares which rows those are.\n")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--contracts":
        # SCOPE: the files named on this command line, and nothing else.
        print(f"contract scope: {len(sys.argv) - 2} file(s) named on the command line")
        n = 0
        m = 0
        for path in sys.argv[2:]:
            text = open(path).read()
            for lineno, block in scan_contracts(text):
                n += 1
                print(f"\nUNSCOPED CONTRACT (attached)  {path}:{lineno}\n     {block}")
            for lineno, block in scan_contracts(text, module_level=True):
                m += 1
                print(f"\nUNSCOPED CONTRACT (module-level)  {path}:{lineno}\n     {block}")
        print(f"\n{n} attached and {m} module-level universal(s) without a declared SCOPE. "
              "A clean run is not an all-clear: this reads only the files it was given, "
              "and module-level blocks are admitted only when worded normatively.")
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
    c_body = scan_contracts(CONTRACT_IN_BODY_MUST_FLAG)
    print(f"  {'ok  ' if c_body else 'MISS'} a contract INSIDE a test body, far from its header, flags")
    # BOTH DIRECTIONS: the defect was misclassification, not absence, so a check
    # that only asserted "module-level finds it" would pass while attached mode
    # also claimed it and the counts stayed wrong.
    m_after = scan_contracts(CONTRACT_MODULE_AFTER_TEST_MUST_FLAG, module_level=True)
    a_after = scan_contracts(CONTRACT_MODULE_AFTER_TEST_MUST_FLAG)
    print(f"  {'ok  ' if m_after else 'MISS'} a MODULE-LEVEL contract after a completed test "
          "reaches the module tier")
    print(f"  {'ok  ' if not a_after else 'MISFILED'} ...and is not also reported as attached")
    ok &= (bool(c_bad) and not c_ok and bool(c_body)
           and bool(m_after) and not a_after)

    print("\nGATE USABLE" if ok else "\nGATE NOT USABLE -- fix before running it on new prose")
    sys.exit(0 if ok else 1)
