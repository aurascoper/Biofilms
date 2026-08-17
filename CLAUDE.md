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
