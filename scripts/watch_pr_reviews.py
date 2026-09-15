#!/usr/bin/env python3
"""Observe Codex review coverage on pull requests. It REPORTS; it authorizes nothing.

    scripts/watch_pr_reviews.py --pr 12 --pr 18 --log watch.jsonl [--once] \
        [--repo OWNER/NAME] [--interval 60] [--deadline 90] [--gh-timeout 60]

Every pass writes one JSON line per pull request to --log and prints one line. Two
things are kept apart in every record, because they are different facts:

  acquisition  ok | failed:<reason>   -- could the service be observed at all?
  coverage     current | stale | none | unrecognized | unknown
               -- what the newest Codex review says about the CURRENT head.

A failed acquisition never reads as "no reviews": the record carries `last_known`, the
previous successful observation of that pull request verbatim (its head, its coverage,
its findings, its own observation time), reloaded from --log on every start.

The previous watcher passed `--arg` to `gh --jq`, discarded the error and printed "no
reviews" while seven findings arrived. Hence: every `gh` call is an argument list with a
timeout; every response is parsed as JSON and checked for shape and for a GraphQL
`errors` array; each of reviews, comments and reviewThreads is paginated on its own
cursor; the head and base are read before and after collection and an observation
during which either moved is discarded. The deadline is monotonic.

Exit codes. --once: 0 every pull request observed (whatever its coverage), 2 any
acquisition failed. Loop: 0 once every pull request is `current` (coverage has arrived
on the head; open findings are for a person to triage, and a fix moves the head again),
3 at the deadline. Merge authorization stays with scripts/preflight_merge.sh and triage.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

CODEX = "chatgpt-codex-connector"
FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
# The two ways Codex names the commit it read: a formal review says
# "Reviewed commit: `sha`"; a comment-form review carries a blob permalink.
SHA_RE = re.compile(r"Reviewed commit:\**\s*`?([0-9a-f]{7,40})`?|/blob/([0-9a-f]{7,40})/")
PERMALINK_RE = re.compile(r"/blob/([0-9a-f]{7,40})/")
BADGE_RE = re.compile(r"badge/(P[0-9])-")
MARKUP_RE = re.compile(r"</?sub>|!\[[^\]]*\]\([^)]*\)|\*\*|\s+")

_PR = "repository(owner:$owner,name:$name){pullRequest(number:$pr){%s}}"
Q_REFS = "query($owner:String!,$name:String!,$pr:Int!){%s}" % (
    _PR % "headRefOid baseRefOid baseRefName")
_PAGED = "query($owner:String!,$name:String!,$pr:Int!,$after:String){%s}"
Q_REVIEWS = _PAGED % (_PR % (
    "reviews(first:100,after:$after){pageInfo{hasNextPage endCursor} "
    "nodes{author{login} submittedAt state body commit{oid}}}"))
Q_COMMENTS = _PAGED % (_PR % (
    "comments(first:100,after:$after){pageInfo{hasNextPage endCursor} "
    "nodes{author{login} createdAt body}}"))
Q_THREADS = _PAGED % (_PR % (
    "reviewThreads(first:100,after:$after){pageInfo{hasNextPage endCursor} "
    "nodes{id isResolved isOutdated path line comments(first:1){nodes{author{login} body}}}}"))
Q_OBJECT = ("query($owner:String!,$name:String!,$expr:String!)"
            "{repository(owner:$owner,name:$name){object(expression:$expr){oid}}}")


class Unavailable(Exception):
    """The service could not be observed. str(self) is the reason."""


def gh(args: list[str], timeout: float) -> dict:
    """Run `gh <args>` and return its stdout as one JSON object, or raise Unavailable."""
    try:
        done = subprocess.run(["gh", *args], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise Unavailable("timeout")
    except FileNotFoundError:
        raise Unavailable("gh-missing")
    if done.returncode != 0:
        raise Unavailable(f"exit:{done.returncode}:{done.stderr.strip()[:160]}")
    try:
        doc = json.loads(done.stdout)
    except json.JSONDecodeError:
        raise Unavailable("non-json")
    if not isinstance(doc, dict):
        raise Unavailable("shape:document")
    if doc.get("errors"):
        kinds = sorted({str(e.get("type") or e.get("message") or "?")
                        for e in doc["errors"] if isinstance(e, dict)})
        raise Unavailable("graphql:" + ",".join(kinds or ["?"]))
    return doc


def gql(query: str, timeout: float, **variables) -> dict:
    args = ["api", "graphql", "-f", f"query={query}"]
    for key, value in variables.items():
        if value is None:            # page one: the nullable $after is simply not bound
            continue
        args += ["-F" if isinstance(value, int) else "-f", f"{key}={value}"]
    doc = gh(args, timeout)
    if not isinstance(doc.get("data"), dict):
        raise Unavailable("shape:data")
    return doc["data"]


def _pull_request(data: dict, what: str) -> dict:
    try:
        node = data["repository"]["pullRequest"]
    except (KeyError, TypeError):
        raise Unavailable(f"shape:{what}")
    if not isinstance(node, dict):
        raise Unavailable(f"shape:{what}")
    return node


def read_refs(repo: str, pr: int, timeout: float) -> tuple[str, str, str]:
    owner, name = repo.split("/")
    node = _pull_request(gql(Q_REFS, timeout, owner=owner, name=name, pr=pr), "refs")
    head, base, base_ref = node.get("headRefOid"), node.get("baseRefOid"), node.get("baseRefName")
    if not (isinstance(head, str) and FULL_SHA.match(head)
            and isinstance(base, str) and FULL_SHA.match(base) and isinstance(base_ref, str)):
        raise Unavailable("shape:refs")
    return head, base, base_ref


def paginate(query: str, connection: str, repo: str, pr: int, timeout: float) -> list:
    """Walk one connection on its own cursor. A failure on any page is a failure of the
    whole walk: no partial list is ever returned."""
    owner, name = repo.split("/")
    nodes, after = [], None
    while True:
        node = _pull_request(gql(query, timeout, owner=owner, name=name, pr=pr, after=after),
                             connection)
        conn = node.get(connection)
        try:
            page, got = conn["pageInfo"], conn["nodes"]
            more, cursor = page["hasNextPage"], page.get("endCursor")
        except (KeyError, TypeError):
            raise Unavailable(f"shape:{connection}")
        if not isinstance(got, list) or not isinstance(more, bool):
            raise Unavailable(f"shape:{connection}")
        nodes.extend(n for n in got if isinstance(n, dict))
        if not more:
            return nodes
        if not isinstance(cursor, str) or not cursor:
            raise Unavailable(f"shape:{connection}.endCursor")
        after = cursor


def resolve_sha(repo: str, prefix: str, timeout: float, cache: dict) -> str | None:
    """A full SHA for an abbreviated one, or None when the object is unknown or the
    prefix is ambiguous (GitHub answers null for both)."""
    if FULL_SHA.match(prefix):
        return prefix
    if prefix not in cache:
        owner, name = repo.split("/")
        data = gql(Q_OBJECT, timeout, owner=owner, name=name, expr=prefix)
        try:
            obj = data["repository"]["object"]
        except (KeyError, TypeError):
            raise Unavailable("shape:object")
        oid = obj.get("oid") if isinstance(obj, dict) else None
        cache[prefix] = oid if isinstance(oid, str) and FULL_SHA.match(oid) else None
    return cache[prefix]


def _login(node: dict) -> str | None:
    author = node.get("author")
    return author.get("login") if isinstance(author, dict) else None


def newest_codex(reviews: list, comments: list) -> dict | None:
    items = []
    for r in reviews:
        if _login(r) == CODEX:
            commit = r.get("commit")
            items.append({"at": r.get("submittedAt") or "", "body": r.get("body") or "",
                          "kind": "review",
                          "oid": commit.get("oid") if isinstance(commit, dict) else None})
    for c in comments:
        if _login(c) == CODEX:
            items.append({"at": c.get("createdAt") or "", "body": c.get("body") or "",
                          "kind": "comment", "oid": None})
    return max(items, key=lambda i: i["at"]) if items else None


def classify(newest: dict | None, head: str, resolve) -> tuple[str, str | None]:
    """current / stale / none / unrecognized, and the full SHA the newest Codex pass read.
    `unrecognized` is a Codex pass whose commit cannot be established; it is never
    reported as `none` and never as `current`."""
    if newest is None:
        return "none", None
    full = newest["oid"] if newest["kind"] == "review" else None
    if not (isinstance(full, str) and FULL_SHA.match(full)):
        m = SHA_RE.search(newest["body"])
        prefix = (m.group(1) or m.group(2)) if m else None
        if not prefix:
            return "unrecognized", None
        full = resolve(prefix)
        if full is None:
            return "unrecognized", None
    return ("current" if full == head else "stale"), full


def _head_of(body: str) -> str:
    return MARKUP_RE.sub(lambda m: " " if m.group(0).isspace() else "", body).strip()[:96]


def thread_findings(threads: list) -> list:
    out = []
    for t in threads:
        if t.get("isResolved"):
            continue
        comments = t.get("comments") or {}
        first = (comments.get("nodes") or [{}])[0] if isinstance(comments, dict) else {}
        body = first.get("body") or ""
        badge = BADGE_RE.search(body)
        out.append({"id": t.get("id"), "path": t.get("path"), "line": t.get("line"),
                    "outdated": bool(t.get("isOutdated")),
                    "severity": badge.group(1) if badge else "UNRANKED",
                    "author": _login(first), "head": _head_of(body)})
    return out


def comment_findings(comments: list, head: str, resolve) -> list:
    """Comment-form Codex findings whose permalink resolves to the current head."""
    out = []
    for c in comments:
        if _login(c) != CODEX:
            continue
        body = c.get("body") or ""
        badges = BADGE_RE.findall(body)
        if not badges:
            continue
        if any(resolve(m.group(1)) == head for m in PERMALINK_RE.finditer(body)):
            out.append({"severity": badges[0], "badges": len(badges),
                        "created_at": c.get("createdAt"), "head": _head_of(body)})
    return out


def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def observe(repo: str, pr: int, timeout: float, max_moves: int = 3) -> dict:
    cache: dict = {}
    resolve = lambda prefix: resolve_sha(repo, prefix, timeout, cache)  # noqa: E731
    for _ in range(max_moves):
        head0, base0, base_ref = read_refs(repo, pr, timeout)
        reviews = paginate(Q_REVIEWS, "reviews", repo, pr, timeout)
        comments = paginate(Q_COMMENTS, "comments", repo, pr, timeout)
        threads = paginate(Q_THREADS, "reviewThreads", repo, pr, timeout)
        head1, base1, _ = read_refs(repo, pr, timeout)
        if (head0, base0) != (head1, base1):
            print(f"pr {pr}: head or base moved during collection "
                  f"({head0[:10]}/{base0[:10]} -> {head1[:10]}/{base1[:10]}); discarded, retrying",
                  file=sys.stderr)
            continue
        newest = newest_codex(reviews, comments)
        coverage, reviewed = classify(newest, head0, resolve)
        return {"observed_at": now_iso(), "repo": repo, "pr": pr,
                "acquisition": "ok", "coverage": coverage,
                "head": head0, "base": base0, "base_ref": base_ref,
                "reviewed_sha": reviewed,
                "reviewed_kind": newest["kind"] if newest else None,
                "reviewed_at": newest["at"] if newest else None,
                "threads": thread_findings(threads),
                "comment_findings": comment_findings(comments, head0, resolve),
                "last_known": None}
    raise Unavailable("moved")


def failed_record(repo: str, pr: int, reason: str, last_known: dict | None) -> dict:
    return {"observed_at": now_iso(), "repo": repo, "pr": pr,
            "acquisition": f"failed:{reason}", "coverage": "unknown",
            "head": None, "base": None, "base_ref": None,
            "reviewed_sha": None, "reviewed_kind": None, "reviewed_at": None,
            "threads": None, "comment_findings": None,
            "last_known": last_known}


def load_last_known(log: Path) -> dict[int, dict]:
    """The newest successful observation per pull request already in the log."""
    last: dict[int, dict] = {}
    if not log.exists():
        return last
    for n, line in enumerate(log.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            print(f"{log}:{n}: unreadable record skipped", file=sys.stderr)
            continue
        if isinstance(rec, dict) and rec.get("acquisition") == "ok" and isinstance(rec.get("pr"), int):
            last[rec["pr"]] = rec
    return last


def one_line(rec: dict) -> str:
    if rec["acquisition"] != "ok":
        lk = rec.get("last_known")
        known = (f"last known {lk['coverage']} at {lk['head'][:10]} observed {lk['observed_at']}"
                 if lk else "nothing previously observed")
        return f"pr {rec['pr']}  ACQUISITION {rec['acquisition']}  coverage unknown  ({known})"
    sev = ",".join(t["severity"] for t in rec["threads"]) or "-"
    reviewed = f"{rec['reviewed_sha'][:10]} ({rec['reviewed_kind']} {rec['reviewed_at']})" \
        if rec["reviewed_sha"] else "none"
    return (f"pr {rec['pr']}  ok  {rec['coverage']:<12} head {rec['head'][:10]}  reviewed {reviewed}"
            f"  unresolved threads {len(rec['threads'])} [{sev}]"
            f"  comment findings {len(rec['comment_findings'])}")


def default_repo(timeout: float) -> str:
    doc = gh(["repo", "view", "--json", "nameWithOwner"], timeout)
    repo = doc.get("nameWithOwner")
    if not (isinstance(repo, str) and repo.count("/") == 1):
        raise Unavailable("shape:repo")
    return repo


def main(argv=None, clock=time.monotonic, sleep=time.sleep) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", help="OWNER/NAME (default: gh repo view)")
    ap.add_argument("--pr", type=int, action="append", required=True, help="repeatable")
    ap.add_argument("--log", type=Path, required=True, help="JSONL observation log")
    ap.add_argument("--once", action="store_true", help="one pass, then exit")
    ap.add_argument("--interval", type=float, default=60.0, help="seconds between passes")
    ap.add_argument("--deadline", type=float, default=90.0, help="minutes, monotonic")
    ap.add_argument("--gh-timeout", type=float, default=60.0, help="seconds per gh call")
    a = ap.parse_args(argv)

    try:
        repo = a.repo or default_repo(a.gh_timeout)
    except Unavailable as e:
        print(f"repository could not be determined: {e}", file=sys.stderr)
        return 2
    last = load_last_known(a.log)
    deadline = clock() + a.deadline * 60
    while True:
        recs = []
        for pr in a.pr:
            try:
                rec = observe(repo, pr, a.gh_timeout)
                last[pr] = rec
            except Unavailable as e:
                rec = failed_record(repo, pr, str(e), last.get(pr))
            with a.log.open("a") as fh:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
            print(one_line(rec), flush=True)
            recs.append(rec)
        if a.once:
            return 0 if all(r["acquisition"] == "ok" for r in recs) else 2
        if all(r["coverage"] == "current" for r in recs):
            return 0
        if clock() >= deadline:
            print(f"deadline reached after {a.deadline:g} minutes; coverage is not current on "
                  "every pull request", file=sys.stderr)
            return 3
        sleep(max(0.0, min(a.interval, deadline - clock())))


if __name__ == "__main__":
    raise SystemExit(main())
