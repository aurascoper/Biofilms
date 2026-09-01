"""Section 3.4's overdamped paragraph states numbers; this binds them to the producer.

SCOPE: the numeric literals in the paragraph at
preprint/modeling_radioresistance_and_radiotropic_fitness.tex lines 571-584, and
the receipt analysis/overdamped_regime.py emits. Nothing else in section 3.4.

THE PARAGRAPH ENDS BY CLAIMING ITS OWN VERIFIABILITY -- "analysis/overdamped_regime.py
reproduces every number in this paragraph" -- and until this file existed nothing
checked that sentence. Worse, the producer was executed by NOTHING: no test, no
workflow, no script invoked it, so its own four assertions (Re < 1, the 1e8 ratio
floor, the 1 cm control, the Re>1 counter-control) had never run in CI. Running it
is the larger half of what this file does; the comparison below is the smaller half.

WHY THE ACCOUNTING IS COMPLETE RATHER THAN A LIST OF FIVE. The receipt used to emit
three scalars while the paragraph states five numbers, so a guard written from the
receipt would have covered two of five and looked finished. The defence is that
EVERY literal in the paragraph must be accounted for: extraction is mechanical, so a
number cannot be missed by an author who did not think to list it, and set equality
against ACCOUNTED is what fails when one appears that nobody mapped.

WHY THE FORMULA EXEMPTION IS SCOPED BY LOCATION AND NOT BY VALUE. The 2, 2 and 9 of
tau_p = 2*rho*a^2/(9*eta) are structure, not results, and must be exempt. A
value-keyed exemption set would be a hole with a comment on it: add a number, add it
to the set, suite green. Membership is therefore decided by whether the occurrence
falls inside FORMULA_SPAN. A stated result cannot be exempted because it does not
occur there -- see test_a_stated_number_cannot_be_moved_into_the_exemption.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
TEX = REPO / "preprint" / "modeling_radioresistance_and_radiotropic_fitness.tex"
PRODUCER = REPO / "analysis" / "overdamped_regime.py"

FIRST_LINE, LAST_LINE = 571, 584
UM = 1e-6

# The one span whose literals are structural. Matched as text, so it moves with the
# paragraph; if the formula is rewritten this raises rather than silently exempting
# nothing (or, worse, exempting a result that drifted into the old offsets).
FORMULA = r"\tau_p = 2\rho a^{2}/(9\eta)"

# Every non-structural literal in the paragraph, as written, with the SI quantity it
# denotes and the receipt path that must reproduce it. `written` is what appears in
# the .tex and fixes the rounding window; `si` is that value in the receipt's units.
ACCOUNTED = [
    # written,  si,        receipt path
    (1e3,       1e3,       ("inputs", "rho")),
    (1e-3,      1e-3,      ("inputs", "eta")),
    (1.0,       1 * UM,    ("lengths", "cell, 1 um at 20 um/s")),
    (10.0,      10 * UM,   ("lengths", "cell, 10 um at 20 um/s")),
    (20.0,      20 * UM,   ("speeds", "cell, 1 um at 20 um/s")),
    (2e-5,      2e-5,      ("reynolds", "cell, 1 um at 20 um/s")),
    (2e-4,      2e-4,      ("reynolds", "cell, 10 um at 20 um/s")),
    (100.0,     100 * UM,  ("lengths", "biofilm feature, 100 um at 20 um/s")),
    (2e-3,      2e-3,      ("reynolds", "biofilm feature, 100 um at 20 um/s")),
    (5.6e-8,    5.6e-8,    ("tau_p", "0.5 um")),
    (0.5,       0.5 * UM,  ("radii", "0.5 um")),
    (2.2e-7,    2.2e-7,    ("tau_p", "1 um")),
    (1.0,       1 * UM,    ("radii", "1 um")),
]


def paragraph() -> str:
    lines = TEX.read_text(encoding="utf-8").splitlines()
    return "\n".join(lines[FIRST_LINE - 1:LAST_LINE])


def _blank_unit_groups(text: str) -> str:
    """Replace \\mathrm{...} with spaces, preserving every character position.

    Unit exponents are not quantities: \\mathrm{kg\\,m^{-3}} carries a -3 that means
    "per cubic metre", not a number the producer computes. Blanking rather than
    deleting keeps offsets valid for the FORMULA_SPAN containment test.
    """
    out = list(text)
    i = 0
    while (j := text.find(r"\mathrm{", i)) != -1:
        depth, k = 0, j + len(r"\mathrm{") - 1
        while k < len(text):
            if text[k] == "{":
                depth += 1
            elif text[k] == "}":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        for m in range(j, min(k + 1, len(text))):
            out[m] = " "
        i = k + 1
    return "".join(out)


TOKEN = re.compile(
    r"(?P<sci>(?P<mant>\d+(?:\.\d+)?)\\times10\^\{(?P<sexp>-?\d+)\})"
    r"|(?P<pow>10\^\{(?P<pexp>-?\d+)\})"
    r"|(?P<bare>\d+(?:\.\d+)?)"
)


def literals(text: str):
    """Every numeric quantity in the paragraph, as (value, start offset)."""
    for m in TOKEN.finditer(_blank_unit_groups(text)):
        if m.group("sci"):
            yield float(m.group("mant")) * 10.0 ** int(m.group("sexp")), m.start()
        elif m.group("pow"):
            yield 10.0 ** int(m.group("pexp")), m.start()
        else:
            yield float(m.group("bare")), m.start()


def half_ulp(written: float) -> float:
    """Half a unit in the last place the .tex actually writes.

    The paragraph states 5.6e-8 for a computed 5.555...e-8, so exact comparison
    fails on correctly-rounded prose. Rounding the computed value and comparing
    strings would pass, but it is a binary verdict on the window and hides where
    inside it the producer sits -- a value at the far edge reads identically to one
    at the centre. A tolerance reports the margin, so drift becomes visible before
    it crosses.
    """
    s = repr(float(written))
    if "e" in s or "E" in s:
        mant, exp = s.lower().split("e")
        decimals = len(mant.split(".")[1].rstrip("0")) if "." in mant else 0
        return 0.5 * 10.0 ** (int(exp) - decimals)
    decimals = len(s.split(".")[1].rstrip("0")) if "." in s else 0
    return 0.5 * 10.0 ** (-decimals)


@pytest.fixture(scope="module")
def receipt(tmp_path_factory):
    out = tmp_path_factory.mktemp("overdamped") / "receipt.json"
    proc = subprocess.run(
        [sys.executable, str(PRODUCER), "--report", str(out)],
        cwd=REPO, capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"producer failed:\n{proc.stdout}\n{proc.stderr}"
    return json.loads(out.read_text())


def test_the_producer_runs_and_its_own_checks_pass(receipt):
    # Pinned so that adding a check to the producer forces someone to look here.
    assert receipt["checks_run"] == 4
    assert receipt["failures"] == 0
    assert receipt["skipped"] == 0 and receipt["skips"] == []
    assert receipt["complete"] is True


def test_the_formula_span_is_still_in_the_paragraph():
    assert paragraph().count(FORMULA) == 1, (
        "the tau_p formula moved or was rewritten; the location-scoped exemption "
        "below is anchored to it and must be re-anchored deliberately.")


def test_every_literal_in_the_paragraph_is_accounted_for():
    text = paragraph()
    start = text.index(FORMULA)
    span = range(start, start + len(FORMULA))

    found = sorted(v for v, off in literals(text) if off not in span)
    structural = sorted(v for v, off in literals(text) if off in span)

    assert structural == [2.0, 2.0, 9.0], (
        f"the formula span holds {structural}, not the 2, 2 and 9 of "
        r"tau_p = 2*rho*a^2/(9*eta)")

    expected = sorted(w for w, _, _ in ACCOUNTED)
    assert found == expected, (
        "section 3.4's paragraph states a number that no entry in ACCOUNTED maps to "
        "a receipt value (or maps one that is no longer there).\n"
        f"  in the paragraph, unaccounted: {sorted(set(found) - set(expected))}\n"
        f"  in ACCOUNTED, not in the text: {sorted(set(expected) - set(found))}")


@pytest.mark.parametrize("written,si,path", ACCOUNTED,
                         ids=[f"{p[0]}:{p[1]}" for _, _, p in ACCOUNTED])
def test_each_stated_number_is_within_half_a_unit_of_the_last_place(
        written, si, path, receipt):
    computed = receipt[path[0]][path[1]]
    tol = half_ulp(written) * (si / written)      # the window, in the receipt's units
    gap = abs(computed - si)
    assert gap <= tol, (
        f"{'/'.join(path)} = {computed!r} is outside the rounding window of the "
        f"{written!r} stated in section 3.4: gap {gap:.3e} > half-ULP {tol:.3e}")
    # Report the margin even on success, so drift inside the window is visible in -rs
    # output rather than silent until it crosses.
    print(f"    {'/'.join(path):45} margin {gap / tol:5.1%} of the rounding window")


def test_the_ten_orders_of_magnitude_claim(receipt):
    # "roughly ten orders of magnitude below the timescales this model addresses"
    import math
    for label in ("0.5 um", "1 um"):
        orders = math.log10(receipt["ratio"][label])
        assert round(orders) == 10, (
            f"T_bio/tau_p for a = {label} is 10^{orders:.2f}, which does not round "
            "to the ten orders of magnitude section 3.4 claims")
