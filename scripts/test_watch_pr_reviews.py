#!/usr/bin/env python3
"""Negative controls for scripts/watch_pr_reviews.py, against a fake `gh` on PATH.

Each control is a known-bad input the watcher must classify as what it is: an error is
an acquisition failure (never "no reviews"), a Codex pass without a resolvable commit is
`unrecognized` (never `none`, never `current`), a finding on page two is found only
through the cursor, a moving head is discarded, a hung `gh` is bounded by its timeout,
the deadline is monotonic, and a previous observation survives a later failure verbatim.

    python3 scripts/test_watch_pr_reviews.py
"""
import inspect
import json
import os
import stat
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import watch_pr_reviews as w  # noqa: E402

HEAD, OTHER, BASE = "a" * 40, "b" * 40, "c" * 40

# The fake decides by scenario (FAKE_GH_SCENARIO) and by what the real script actually
# passed on argv. Page two exists only behind `after=C1` as its OWN argv element, which is
# what `-f after=C1` produces; a query text containing `after:$after` never matches it.
FAKE_GH = r'''#!/usr/bin/env python3
import json, os, sys, time
argv = sys.argv[1:]; joined = " ".join(argv)
S = os.environ.get("FAKE_GH_SCENARIO", "current")
HEAD, OTHER, BASE = "a"*40, "b"*40, "c"*40
CODEX = "chatgpt-codex-connector"
if os.environ.get("FAKE_GH_SLEEP"):
    time.sleep(float(os.environ["FAKE_GH_SLEEP"]))
def out(obj):
    print(json.dumps(obj)); sys.exit(0)
def pr(inner):
    return {"data": {"repository": {"pullRequest": inner}}}
def conn(name, nodes, more=False):
    return pr({name: {"pageInfo": {"hasNextPage": more, "endCursor": "C1" if more else None},
                      "nodes": nodes}})
if argv[:2] == ["repo", "view"]:
    out({"nameWithOwner": "o/r"})
if S == "gh_fail":
    sys.stderr.write("boom\n"); sys.exit(1)
if S == "malformed":
    print("{not json"); sys.exit(0)
if S == "gql_errors":
    out({"data": None, "errors": [{"type": "RATE_LIMITED", "message": "API rate limit exceeded"}]})
page2 = "after=C1" in argv
if "headRefOid" in joined:
    head = HEAD
    if S == "moved":
        f = os.environ["FAKE_GH_STATE"]
        n = (int(open(f).read() or 0) if os.path.exists(f) else 0) + 1
        open(f, "w").write(str(n))
        head = HEAD if n == 1 else OTHER
    out(pr({"headRefOid": head, "baseRefOid": BASE, "baseRefName": "master"}))
if "object(expression" in joined:
    expr = [a for a in argv if a.startswith("expr=")][0][5:]
    oid = None
    if S in ("stale", "stale_comment") and OTHER.startswith(expr): oid = OTHER
    if S in ("comment_finding",) and HEAD.startswith(expr): oid = HEAD
    out({"data": {"repository": {"object": {"oid": oid} if oid else None}}})
if S == "cursor_fail" and page2:
    sys.stderr.write("cursor exploded\n"); sys.exit(1)
badge = lambda p: "**<sub><sub>![%s Badge](https://img.shields.io/badge/%s-orange?style=flat)</sub></sub>  " % (p, p)
review = {"author": {"login": CODEX}, "submittedAt": "2026-09-14T02:51:05Z", "state": "COMMENTED",
          "body": "### Codex Review\n\n**Reviewed commit:** `%s`" % HEAD[:10], "commit": {"oid": HEAD}}
thread = {"id": "PRRT_1", "isResolved": False, "isOutdated": False, "path": "a.py", "line": 7,
          "comments": {"nodes": [{"author": {"login": CODEX},
                                  "body": badge("P1") + "Keep unrelated initialization outside the health await**\n\nWhen the cache is absent..."}]}}
resolved = dict(thread, id="PRRT_0", isResolved=True, path="done.py")
if "reviews(" in joined:
    if S in ("empty", "unrecognized_comment", "stale_comment", "comment_finding"):
        out(conn("reviews", []))
    if S == "stale":
        out(conn("reviews", [dict(review, commit={"oid": OTHER}, body="**Reviewed commit:** `%s`" % OTHER[:10])]))
    if S == "unrecognized":
        out(conn("reviews", [dict(review, commit=None, body="Codex could not review this pull request.")]))
    out(conn("reviews", [review]))
if "reviewThreads(" in joined:
    if S == "empty":
        out(conn("reviewThreads", []))
    if S in ("page2", "cursor_fail"):
        if page2:
            out(conn("reviewThreads", [dict(thread, id="PRRT_2", path="page2.py")]))
        out(conn("reviewThreads", [resolved], more=True))
    out(conn("reviewThreads", [resolved, thread]))
if "comments(" in joined:
    if S == "comment_finding":
        out(conn("comments", [{"author": {"login": CODEX}, "createdAt": "2026-09-14T03:00:00Z",
            "body": "https://github.com/o/r/blob/%s/f.py#L9-L9\n" % HEAD[:12] + badge("P2") + "Comment-form finding**"}]))
    if S == "stale_comment":
        out(conn("comments", [{"author": {"login": CODEX}, "createdAt": "2026-09-14T03:00:00Z",
            "body": "https://github.com/o/r/blob/%s/f.py#L9-L9\n" % OTHER[:12] + badge("P2") + "Old finding**"}]))
    if S == "unrecognized_comment":
        out(conn("comments", [{"author": {"login": CODEX}, "createdAt": "2026-09-14T03:00:00Z",
            "body": "Codex usage limit reached; review not performed."}]))
    out(conn("comments", []))
sys.stderr.write("fake gh: unhandled call " + joined[:100] + "\n"); sys.exit(9)
'''


class WatcherControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        fake = Path(cls.tmp.name) / "bin" / "gh"
        fake.parent.mkdir()
        fake.write_text(FAKE_GH)
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        os.environ["PATH"] = str(fake.parent) + os.pathsep + os.environ["PATH"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        os.environ.pop("FAKE_GH_SLEEP", None)
        os.environ.pop("FAKE_GH_STATE", None)
        self.log = Path(self.tmp.name) / f"{self._testMethodName}.jsonl"

    def run_once(self, scenario, extra=(), log=None):
        os.environ["FAKE_GH_SCENARIO"] = scenario
        log = log or self.log
        rc = w.main(["--repo", "o/r", "--pr", "1", "--log", str(log), "--once", *extra])
        recs = [json.loads(line) for line in log.read_text().splitlines()]
        return rc, recs[-1], recs

    def assert_failed(self, rec, prefix):
        self.assertTrue(rec["acquisition"].startswith("failed:" + prefix), rec["acquisition"])
        self.assertEqual(rec["coverage"], "unknown")
        self.assertIsNone(rec["threads"])
        self.assertIsNone(rec["head"])

    # -- coverage classes --------------------------------------------------------------
    def test_a_current_head_finding_is_reported(self):
        rc, rec, _ = self.run_once("current")
        self.assertEqual(rc, 0)
        self.assertEqual(rec["acquisition"], "ok")
        self.assertEqual(rec["coverage"], "current")
        self.assertEqual(rec["reviewed_sha"], HEAD)
        self.assertEqual(rec["head"], HEAD)
        self.assertEqual([t["path"] for t in rec["threads"]], ["a.py"])   # the resolved one is not a finding
        self.assertEqual(rec["threads"][0]["severity"], "P1")
        self.assertTrue(rec["threads"][0]["head"].startswith("Keep unrelated initialization"))

    def test_valid_empty_collections_are_none_not_an_error(self):
        rc, rec, _ = self.run_once("empty")
        self.assertEqual(rc, 0)
        self.assertEqual(rec["acquisition"], "ok")
        self.assertEqual(rec["coverage"], "none")
        self.assertEqual(rec["threads"], [])
        self.assertIsNone(rec["reviewed_sha"])

    def test_a_stale_review_is_stale(self):
        rc, rec, _ = self.run_once("stale")
        self.assertEqual(rc, 0)
        self.assertEqual(rec["coverage"], "stale")
        self.assertEqual(rec["reviewed_sha"], OTHER)

    def test_a_stale_comment_form_review_resolves_through_the_object_query(self):
        _, rec, _ = self.run_once("stale_comment")
        self.assertEqual(rec["coverage"], "stale")
        self.assertEqual(rec["reviewed_kind"], "comment")
        self.assertEqual(rec["comment_findings"], [])          # its permalink is not the head

    def test_a_comment_form_finding_on_the_head_is_current_and_counted(self):
        _, rec, _ = self.run_once("comment_finding")
        self.assertEqual(rec["coverage"], "current")
        self.assertEqual(len(rec["comment_findings"]), 1)
        self.assertEqual(rec["comment_findings"][0]["severity"], "P2")

    def test_a_codex_review_without_a_resolvable_commit_is_unrecognized(self):
        rc, rec, _ = self.run_once("unrecognized")
        self.assertEqual(rc, 0)                                 # observed; it is the coverage that is unknown
        self.assertEqual(rec["coverage"], "unrecognized")
        self.assertIsNone(rec["reviewed_sha"])

    def test_a_codex_comment_without_a_commit_is_unrecognized_not_none(self):
        _, rec, _ = self.run_once("unrecognized_comment")
        self.assertEqual(rec["coverage"], "unrecognized")
        self.assertEqual(rec["reviewed_kind"], "comment")

    def test_a_usage_limit_notice_is_service_unavailable_not_none_and_not_current(self):
        rc, rec, _ = self.run_once("usage_limit")           # a current review exists, the notice is newer
        self.assertEqual(rc, 4)
        self.assertEqual(rec["acquisition"], "ok")
        self.assertEqual(rec["coverage"], "service_unavailable")
        self.assertIsNone(rec["reviewed_sha"])
        self.assertEqual([t["path"] for t in rec["threads"]], ["a.py"])   # threads are still observed
        self.assertIn("SERVICE UNAVAILABLE", w.one_line(rec))

    def test_a_declined_review_keeps_the_previously_established_coverage(self):
        _, good, _ = self.run_once("current")
        rc, rec, _ = self.run_once("usage_limit", log=self.log)
        self.assertEqual(rec["last_known"], good)
        _, later, _ = self.run_once("gh_fail", log=self.log)   # and a later failure still points at it
        self.assertEqual(later["last_known"], good)

    def test_the_loop_exits_4_when_the_service_declines_instead_of_polling(self):
        os.environ["FAKE_GH_SCENARIO"] = "usage_limit"
        rc = w.main(["--repo", "o/r", "--pr", "1", "--log", str(self.log)],
                    clock=lambda: 0.0, sleep=lambda s: self.fail("polled although the service declined"))
        self.assertEqual(rc, 4)

    # -- acquisition failures ----------------------------------------------------------
    def test_malformed_json_is_an_acquisition_failure(self):
        rc, rec, _ = self.run_once("malformed")
        self.assertEqual(rc, 2)
        self.assert_failed(rec, "non-json")

    def test_a_failing_gh_is_an_acquisition_failure(self):
        rc, rec, _ = self.run_once("gh_fail")
        self.assertEqual(rc, 2)
        self.assert_failed(rec, "exit:1")

    def test_a_graphql_errors_array_with_exit_zero_is_an_acquisition_failure(self):
        rc, rec, _ = self.run_once("gql_errors")
        self.assertEqual(rc, 2)
        self.assert_failed(rec, "graphql:RATE_LIMITED")

    # -- pagination --------------------------------------------------------------------
    def test_the_fake_hides_the_page_two_finding_from_a_one_page_reader(self):
        """The control on the control: page one alone carries no unresolved thread."""
        import subprocess
        os.environ["FAKE_GH_SCENARIO"] = "page2"
        one = json.loads(subprocess.run(["gh", "api", "graphql", "-f", "query=reviewThreads("],
                                        capture_output=True, text=True).stdout)
        conn = one["data"]["repository"]["pullRequest"]["reviewThreads"]
        self.assertTrue(conn["pageInfo"]["hasNextPage"])
        self.assertTrue(all(n["isResolved"] for n in conn["nodes"]))

    def test_a_finding_on_page_two_is_reached_through_the_cursor(self):
        rc, rec, _ = self.run_once("page2")
        self.assertEqual(rc, 0)
        self.assertEqual([t["path"] for t in rec["threads"]], ["page2.py"])

    def test_a_cursor_failure_yields_no_partial_findings(self):
        rc, rec, _ = self.run_once("cursor_fail")
        self.assertEqual(rc, 2)
        self.assert_failed(rec, "exit:1")

    # -- the head moving, hangs, deadlines ---------------------------------------------
    def test_an_observation_during_which_the_head_moved_is_discarded_and_retried(self):
        state = Path(self.tmp.name) / "moved.count"
        os.environ["FAKE_GH_STATE"] = str(state)
        rc, rec, _ = self.run_once("moved")
        self.assertEqual(rc, 0)
        self.assertEqual(rec["head"], OTHER)                    # the settled head, not the first read
        self.assertEqual(state.read_text(), "4")                # 2 reads discarded + 2 reads that agreed

    def test_a_hung_gh_is_bounded_by_the_per_call_timeout(self):
        os.environ["FAKE_GH_SLEEP"] = "5"
        t0 = time.monotonic()
        rc, rec, _ = self.run_once("current", extra=["--gh-timeout", "1"])
        self.assertLess(time.monotonic() - t0, 4.0)
        self.assertEqual(rc, 2)
        self.assert_failed(rec, "timeout")

    def test_the_deadline_is_monotonic_and_ends_the_loop_with_exit_3(self):
        self.assertIs(inspect.signature(w.main).parameters["clock"].default, time.monotonic)
        os.environ["FAKE_GH_SCENARIO"] = "stale"
        now = [0.0]
        rc = w.main(["--repo", "o/r", "--pr", "1", "--log", str(self.log),
                     "--deadline", "1", "--interval", "30"],
                    clock=lambda: now[0], sleep=lambda s: now.__setitem__(0, now[0] + s))
        self.assertEqual(rc, 3)
        self.assertGreaterEqual(len(self.log.read_text().splitlines()), 2)

    def test_the_loop_stops_at_zero_once_every_pull_request_is_current(self):
        os.environ["FAKE_GH_SCENARIO"] = "current"
        rc = w.main(["--repo", "o/r", "--pr", "1", "--log", str(self.log)],
                    clock=lambda: 0.0, sleep=lambda s: self.fail("slept although current"))
        self.assertEqual(rc, 0)
        self.assertEqual(len(self.log.read_text().splitlines()), 1)

    # -- persistence -------------------------------------------------------------------
    def test_a_previous_observation_survives_a_later_failure_verbatim_across_a_restart(self):
        _, good, _ = self.run_once("current")
        rc, bad, recs = self.run_once("gh_fail", log=self.log)   # a new process would reload the same log
        self.assertEqual(rc, 2)
        self.assertEqual(len(recs), 2)
        self.assert_failed(bad, "exit:1")
        self.assertEqual(bad["last_known"], good)               # original head, findings, observation time
        self.assertIsNone(good["last_known"])

    def test_nothing_previously_observed_is_said_so(self):
        _, bad, _ = self.run_once("gh_fail")
        self.assertIsNone(bad["last_known"])
        self.assertIn("nothing previously observed", w.one_line(bad))


if __name__ == "__main__":
    unittest.main(verbosity=2)
