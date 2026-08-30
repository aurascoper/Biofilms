#!/usr/bin/env python3
"""Mechanical candidate list for the absence-sentence gate.

Grep proposes; a human disposes. The point is to stop reading one's own prose for
a shape one is primed to miss, so this errs toward flagging and expects a human
'no' on some of what it returns.
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

if __name__ == "__main__":
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
    print("\nGATE USABLE" if ok else "\nGATE NOT USABLE -- fix before running it on new prose")
    sys.exit(0 if ok else 1)
