# Review coverage: what the merge gate accepts

**Date:** 2026-09-15 · **Gate:** `scripts/preflight_merge.sh` · **Controls:**
`scripts/test_preflight_merge.sh` · **Self-review source expires:** 2026-10-15, enforced by the gate

## Before

Coverage meant one thing. The newest review by `chatgpt-codex-connector` had to name the head,
either as "Reviewed commit:" in a review body or through a blob permalink in a comment. Anything
else was "no review".

## Now

There are three sources. The newest one that names a commit decides whether the head is covered
or stale.

| source | what counts | its findings |
|---|---|---|
| Codex | unchanged | a badge in a body that names the head blocks; threads block |
| Copilot | a review by `copilot-pull-request-reviewer` whose commit is the head | arrive as review threads, which block |
| Self-review | a comment or review body by an owner, member or collaborator, with a `review-coverage:` line and its findings | `open` blocks; `fixed <sha>` and `deferred <where>` do not |

A self-review coverage comment has this shape:

```
review-coverage: <sha>
scope: <files or diff range>; checklist: <name>
findings:
- P1 fixed 71a5c90: <title>
- P2 deferred <where it is tracked>: <title>
- P2 open: <title>
```

`findings: none` is accepted only beside a `scope:` line that names both files and a checklist.
When a coverage comment names the head, the gate refuses it by name if it:
- has no findings section;
- lists nothing;
- says `none` without a scope;
- carries a finding line the gate cannot parse;
- fixes or defers a finding without saying where.

A marker from anyone who is not an owner, member or collaborator is reported and not counted.
Every coverage comment that names the head is read, so a disposition is changed by editing the
comment.

## Why

Codex began declining on its usage limit on 2026-09-15 at 05:40 UTC. Copilot's automatic review,
turned on by the `base` ruleset, has posted nothing since 2026-09-13 at 23:43 UTC. A gate that
counted only Codex refused every head, for a reason that says nothing about whether anyone read
the code. The cause is a subscription lapse, not a judgment of Codex's value.

A substitute review had already run: `/code-review` reviewed #12 and posted its findings. The gate
had no way to count it.

## What did not change

- Every unresolved review thread blocks, whoever opened it.
- A review of a commit that is no longer the head is stale, whichever source posted it.
- `scripts/watch_pr_reviews.py` observes Codex only and authorizes nothing. It keeps reporting
  `service_unavailable` while the gate accepts a substitute.

## Expiry: 2026-10-15

The self-review source lowers the bar, because the review is written by the same side that wrote
the code. After 2026-10-15 the gate stops counting self-review comments and prints this file's
path, because a lowering with no date becomes permanent by default. On that date, choose one:
- restore Codex-only coverage;
- extend the self-review source, with a new date and a reason;
- keep Copilot as the only substitute.

The Copilot source does not expire.
