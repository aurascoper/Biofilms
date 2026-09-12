#!/usr/bin/env python3
"""Assert the numpy depth-stencil receipt describes the run it claims.

WHY THIS IS NOT A COUNT.  The first version checked only that the receipt held
three cases, so three copies of ``default_regime`` passed and silently losing
``strong_gradient`` or the cylindrical negative control stayed green (Codex P2
on #19).  A gate that counts its inputs cannot tell which ones it got.

So the expected case names and the expected check names WITHIN each case are
both pinned.  A case whose internal checks were deleted has an empty ``checks``
dict, `all([])` is True, and it reports ``passed`` -- which is the same
can't-fail shape one level down.

``--self-test`` runs the gate against known-bad receipts and requires each to
be rejected, because a gate nobody can test is a gate nobody can trust.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

EXPECTED: dict[str, set[str]] = {
    "default_regime": {"accuracy", "second_order",
                       "negative_control_fails_slab_bound"},
    "strong_gradient": {"accuracy", "second_order",
                        "negative_control_fails_slab_bound"},
    "cylindrical_refactor_equivalence": {"matches_pre_change_stencil",
                                         "swap_control_differs"},
}


def check(report: dict) -> list[str]:
    """Return a list of problems; empty means the receipt is trustworthy."""
    bad: list[str] = []
    cases = report.get("cases")
    if not isinstance(cases, list):
        return ["receipt has no 'cases' list"]

    names = [c.get("case") for c in cases]
    if len(set(names)) != len(names):
        bad.append(f"duplicate case names: {names}")
    if set(names) != set(EXPECTED):
        missing = sorted(set(EXPECTED) - set(names))
        extra = sorted(set(names) - set(EXPECTED))
        bad.append(f"case identities differ -- missing {missing}, unexpected {extra}"
                   " (update EXPECTED if this is deliberate)")

    for case in cases:
        name = case.get("case")
        if name not in EXPECTED:
            continue
        got = set(case.get("checks", {}))
        if got != EXPECTED[name]:
            bad.append(f"{name}: checks differ -- missing "
                       f"{sorted(EXPECTED[name] - got)}, unexpected {sorted(got - EXPECTED[name])}")
        # `is not True`, not `not ok`: JSON "false" is a non-empty string and
        # therefore truthy, so a producer schema regression emitting string
        # verdicts was certified as passing (Codex on #19).  A verdict must be
        # an actual boolean true.
        for cname, ok in case.get("checks", {}).items():
            if ok is not True:
                bad.append(f"{name}: check {cname} is {ok!r}, not boolean true")
        if case.get("passed") is not True:
            bad.append(f"{name}: case verdict is {case.get('passed')!r}, not boolean true")

    if report.get("passed") is not True:
        bad.append(f"overall verdict is {report.get('passed')!r}, not boolean true")
    return bad


def self_test() -> int:
    """Every doctored receipt below must be REJECTED."""
    good = {
        "passed": True,
        "cases": [{"case": n, "passed": True, "checks": {c: True for c in ch}}
                  for n, ch in EXPECTED.items()],
    }
    if check(good):
        print(f"SELF-TEST BROKEN: a valid receipt was rejected: {check(good)}")
        return 1

    def drop_case(r):
        r["cases"] = r["cases"][:2]; return r

    def duplicate_case(r):
        r["cases"][1] = json.loads(json.dumps(r["cases"][0])); return r

    def gut_checks(r):
        r["cases"][0]["checks"] = {}; return r

    def flip_control(r):
        r["cases"][2]["checks"]["swap_control_differs"] = False; return r

    def rename_case(r):
        r["cases"][1]["case"] = "strong_gradient_v2"; return r

    def overall_false(r):
        r["passed"] = False; return r

    def string_verdicts(r):
        """The truthiness trap: JSON "false" is a non-empty string."""
        r["cases"][0]["passed"] = "false"
        r["cases"][0]["checks"] = {c: "false" for c in r["cases"][0]["checks"]}
        return r

    def overall_string_true(r):
        r["passed"] = "true"; return r

    failures = 0
    for name, mutate in [("a case dropped", drop_case),
                         ("a case duplicated", duplicate_case),
                         ("a case's checks deleted", gut_checks),
                         ("the swap control flipped", flip_control),
                         ("a case renamed", rename_case),
                         ("overall pass false", overall_false),
                         ("string verdicts \"false\"", string_verdicts),
                         ("overall verdict a string", overall_string_true)]:
        problems = check(mutate(json.loads(json.dumps(good))))
        status = "REJECTED" if problems else "ACCEPTED -- GATE IS BLIND"
        print(f"  {name:<28} {status}")
        if not problems:
            failures += 1
    print("self-test: gate rejects every doctored receipt."
          if not failures else f"self-test: {failures} doctored receipt(s) slipped through.")
    return 1 if failures else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("receipt", nargs="?", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.receipt is None:
        return ap.error("a receipt path is required unless --self-test is given")
    if not args.receipt.exists():
        print(f"no receipt at {args.receipt}: the verifier did not run to completion")
        return 1

    report = json.loads(args.receipt.read_text())  # invalid JSON raises, which is the point
    problems = check(report)
    names = [c.get("case") for c in report.get("cases", [])]
    print(f"receipt: {len(names)} cases -- {', '.join(map(str, names))}")
    for p in problems:
        print(f"  REJECTED: {p}")
    if problems:
        return 1
    print("numpy verifier ran complete, with the expected cases and checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
