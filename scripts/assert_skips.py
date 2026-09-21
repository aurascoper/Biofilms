#!/usr/bin/env python3
"""Fail when a pytest run skipped a test it is not declared to skip.

    python3 scripts/assert_skips.py <junit.xml> <expected file>
    python3 scripts/assert_skips.py --self-test

AGENTS.md rule 2: a skipped test is uncovered surface. `pytest -rs` makes a skip VISIBLE
in a log, and the workflow comment said "the only skips that should remain here are the
openmc ones". Nothing enforced that sentence. A runner that half-installs its packages
turns a failing test into a skipped one, the count stays green, and nobody reads the log.

The expected file has one skip per line: the test id, a tab, and text the reason must
contain. A `#` line is a comment. The id is "<classname>::<name>", or "<name>" alone for
a module that was skipped at collection, which is how pytest reports `importorskip`.

The comparison is on the SET of ids, in both directions. An undeclared skip fails. A
declared skip that did not happen fails too: the test now runs, which is good news, and
the list is then stale and must say so.
"""
import sys
import xml.etree.ElementTree as ET


def skipped(report: str) -> dict:
    """id -> reason text, for every skipped testcase. Refuses a report with no testcase."""
    root = ET.parse(report).getroot()
    cases = list(root.iter("testcase"))
    if not cases:
        raise SystemExit(f"{report}: no testcase in the report; nothing was checked")
    out = {}
    for tc in cases:
        mark = tc.find("skipped")
        if mark is None:
            continue
        cls, name = tc.attrib.get("classname", ""), tc.attrib["name"]
        out[f"{cls}::{name}" if cls else name] = f"{mark.attrib.get('message', '')} {mark.text or ''}"
    return out


def expected(path: str) -> dict:
    """id -> text the reason must contain."""
    out = {}
    for n, line in enumerate(open(path, encoding="utf-8").read().split("\n"), 1):
        if not line.strip() or line.startswith("#"):
            continue
        if "\t" not in line:
            raise SystemExit(f"{path}:{n}: want '<test id><TAB><reason text>', got {line!r}")
        ident, reason = line.split("\t", 1)
        out[ident.strip()] = reason.strip()
    return out


def problems(got: dict, want: dict) -> list:
    out = [f"UNDECLARED SKIP  {i}: {got[i].strip()[:160]}" for i in sorted(set(got) - set(want))]
    out += [f"DECLARED, BUT IT RAN  {i}: remove it from the list" for i in sorted(set(want) - set(got))]
    out += [f"SKIPPED FOR ANOTHER REASON  {i}: want {want[i]!r} in {got[i].strip()[:160]!r}"
            for i in sorted(set(got) & set(want)) if want[i] not in got[i]]
    return out


def self_test() -> int:
    """Each case is a report and a list, and the exit status THIS SCRIPT must give.

    Every case runs the script as CI runs it, as a command, and reads its exit status.
    The first version called problems() directly and made up its own verdict. With the
    last line of the script changed to `sys.exit(0)` it stayed green, while every CI
    assertion would have accepted an undeclared skip. Raised by Codex on pull request
    #48. AGENTS.md rule 1: a control enters by the production door."""
    import os
    import subprocess
    import tempfile

    def cli(report: str, wanted: str) -> int:
        return subprocess.run([sys.executable, os.path.abspath(__file__), report, wanted],
                              capture_output=True, text=True).returncode
    # The two shapes pytest writes, copied from real reports of this repository:
    # a module skipped at collection, and a test skipped by a skipif mark.
    module = ('<testcase classname="" name="tests.test_model_build" time="0.0">'
              '<skipped message="collection skipped">(\'/x/test_model_build.py\', 10, '
              '"Skipped: could not import \'openmc\': No module named \'openmc\'")</skipped></testcase>')
    marked = ('<testcase classname="tests.test_pilot" name="test_three_channels" time="0.0">'
              '<skipped type="pytest.skip" message="Dataset_3.nd2 not present; Dryad serves it">'
              '/x/test_pilot.py:243: Dataset_3.nd2 not present</skipped></testcase>')
    passed = '<testcase classname="tests.test_a" name="test_ok" time="0.1"/>'
    render = ('<testcase classname="" name="tests.test_viewer" time="0.0">'
              '<skipped message="collection skipped">(\'/x/test_viewer.py\', 3, '
              '"Skipped: could not import \'pyvista\': libGL.so.1: cannot open")</skipped></testcase>')
    want_both = "tests.test_model_build\tcould not import 'openmc'\ntests.test_pilot::test_three_channels\tnot present; Dryad\n"
    cases = [
        ("the declared skips and nothing else", passed + module + marked, want_both, 0),
        ("a seventh skip: the renderer did not import", passed + module + marked + render, want_both, 1),
        ("a declared skip that ran", passed + module, want_both, 1),
        ("a declared test skipped for another reason", passed + module + marked,
         want_both.replace("not present; Dryad", "could not import 'nd2'"), 1),
        ("no skip, and none declared", passed, "# none\n", 0),
        ("a skip, and none declared", passed + module, "# none\n", 1),
    ]
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for label, body, want, verdict in cases:
            rep, exp = os.path.join(d, "r.xml"), os.path.join(d, "e.txt")
            open(rep, "w").write(f'<testsuites><testsuite name="pytest">{body}</testsuite></testsuites>')
            open(exp, "w").write(want)
            got = cli(rep, exp)
            print(("ok     " if got == verdict else "WRONG  ") + f"{label}: verdict {got}, wanted {verdict}")
            bad |= got != verdict
        # A report with no testcase must be refused, never read as "no skips". So must a
        # list line with no tab, and a report that is not there.
        exp = os.path.join(d, "none.txt")
        open(exp, "w").write("# none\n")
        rep = os.path.join(d, "empty.xml")
        open(rep, "w").write('<testsuites><testsuite name="pytest"/></testsuites>')
        full = os.path.join(d, "full.xml")
        open(full, "w").write(f'<testsuites><testsuite name="pytest">{passed}</testsuite></testsuites>')
        untabbed = os.path.join(d, "untabbed.txt")
        open(untabbed, "w").write("tests.test_model_build could not import\n")
        for label, r, e in (("an empty report", rep, exp),
                            ("a list line with no tab", full, untabbed),
                            ("a report that is not there", os.path.join(d, "absent.xml"), exp)):
            got = cli(r, e)
            print(("ok     " if got else "WRONG  ") + f"{label} is refused: exit {got}")
            bad |= got == 0
    print("all skip pin controls gave the expected verdict" if not bad else "AT LEAST ONE CONTROL IS WRONG")
    return int(bad)


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        sys.exit(self_test())
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    found = problems(skipped(sys.argv[1]), expected(sys.argv[2]))
    print(f"{sys.argv[1]}: {len(skipped(sys.argv[1]))} skipped, {len(expected(sys.argv[2]))} declared")
    for line in found:
        print(line)
    sys.exit(1 if found else 0)
