# CLAUDE.md

The standards for this repository live in one file, shared with every other
agent that works here, so that two reviewers cannot be held to two different
bars:

@AGENTS.md

Everything above applies. What follows is only what is specific to working
through Claude Code.

## Reviews are not optional here

Codex reviews every pull request in this repository and has raised a P1 on four
consecutive ones. **Read its comments before merging, not after.**

```sh
gh pr view <pr> --json reviews -q '.reviews[].body'
scripts/preflight_merge.sh <pr>
```

Ten issues reached `master` across PRs #7–#10 because those comments were never
read — several P1, in code that had already been reported as finished. Reporting
work as done is a claim about the code, not about the effort spent on it.

## Verify a finding before acting on it, and before dismissing it

Reviews from other agents are evidence, not verdicts. Reproduce the failure —
run the function, print the value, check the import actually resolves — then fix
what you confirmed. Say plainly which findings you verified and which you did
not.

This cuts the other way too. A measurement of your own can be wrong: the damage
from the claims-ledger bug was first measured against an already-revised file
and read as 1/30, which nearly produced a second false conclusion on top of the
first. Check what the control actually contains before trusting what it reports.

## Do not narrate certainty you do not have

State what was run and what it returned. "Conservation to nine significant
figures" was true of one check and false about what that check verified; the
masses it compared differed by 3%. If a suite skipped, say it skipped.

## Work in a worktree, not the shared checkout

Hunter works in these trees at the same time you do. A shared checkout means a
shared branch pointer, and a branch can change under a running session between
one command and the next — so a commit lands wherever `HEAD` happens to point,
not where the work belongs.

This happened **twice in one session** on 2026-08-28, during PR #19:

- `0af4505`, Hunter's claims-ledger work, appeared on `feat/slab-depth-geometry`
  mid-session. Recovered onto `fix/claims-ledger-delete-verdicts`.
- `fix/ledger-guard-document-scope` was checked out under the running session,
  and census commit `5b100d7` landed on it instead of the feature branch.
  Recovered as `f49fe94`.

Both times the confirming command lied. `git push origin feat/slab-depth-geometry`
printed **`Everything up-to-date`** while local `HEAD` was a commit ahead — because
the branch *named* in the push was not the branch *checked out*, and the named one
genuinely had nothing new. **After committing, treat `Everything up-to-date` as
evidence the commit went somewhere else**, not as confirmation of anything.

Start in a worktree: `EnterWorktree`, or `git worktree add`. Note that
`worktree.baseRef` defaults to `fresh`, which branches from `origin/master`;
stacked work needs `head`. That a worktree already exists is not the rule — one
did (`Biofilms-v11`, on `fix/figure-artifacts-and-v11`) while both collisions
happened. The tooling was in use and the practice was not, because the mitigation
had been written into one plan rather than adopted as standing practice. A
mitigation scoped to a single piece of work is not in force for the next session.

Read the branch and the commit in **one invocation**:

    git status --porcelain=v2 --branch

which reports `branch.oid` and `branch.head` together, from a single read. Two
git commands on one shell line are not one snapshot — `&&` sequences them, it
does not make them atomic, so `git status -sb` can report the old branch while
`git log -1` reports the new commit. That is this section's own failure, and an
earlier draft of this paragraph prescribed exactly that two-command check.
A push's exit status is not evidence of where a commit landed either.
