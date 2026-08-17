#!/usr/bin/env bash
# Refuse to merge a pull request while a review finding is unresolved.
#
# WHY THIS EXISTS. Codex reviews every pull request in this repository and has
# raised a P1 on four consecutive ones. The reviews were never the missing
# piece -- nothing stopped a merge while their findings sat unread. Ten issues
# reached master across PRs #7-#10 that way, several P1, in code that had
# already been reported as done; PR #12 was one command away from merging over
# a NameError and a biosafety milestone that read AUTHORIZED while its own gate
# refused.
#
# CI answers "do the tests pass". It cannot answer "did anyone read the
# review", because a finding lives in a review thread and no status check
# knows the difference between a thread that was addressed and one that was
# ignored. This does.
#
# It also refuses a STALE review: if the newest Codex review names a commit
# that is no longer the head, the findings describe code that has since been
# rewritten, and their silence means nothing. Push, wait for the re-review,
# run this again.
#
# Usage:  scripts/preflight_merge.sh <pr-number>
# Exit:   0 nothing outstanding · 1 unresolved findings, stale review, or no
#         review at all
set -euo pipefail

PR="${1:?usage: preflight_merge.sh <pr-number>}"

# A gate nobody can test is a gate nobody can trust, and this one shipped with a
# branch that let P3 threads through while claiming to refuse every unresolved
# thread. PREFLIGHT_FIXTURE substitutes a canned API response so the decision
# logic can be exercised against synthetic threads. See test_preflight_merge.sh.
if [[ -n "${PREFLIGHT_FIXTURE:-}" ]]; then
  DATA="$(cat "$PREFLIGHT_FIXTURE")"
else
REPO="$(gh repo view --json nameWithOwner -q .nameWithOwner)"
OWNER="${REPO%%/*}"
NAME="${REPO##*/}"

# --paginate with a $endCursor in the query walks every page and concatenates
# the responses. A FIXED `first:100` WOULD MAKE THIS GATE FAIL OPEN: past a
# hundred threads the unseen ones are indistinguishable from resolved ones, and
# a gate that silently truncates its input is the same defect it exists to
# catch. `jq -s` folds the pages back into one document.
DATA="$(gh api graphql --paginate \
  -f owner="$OWNER" -f name="$NAME" -F pr="$PR" -f query='
query($owner:String!, $name:String!, $pr:Int!, $endCursor:String) {
  repository(owner:$owner, name:$name) {
    pullRequest(number:$pr) {
      title
      headRefOid
      reviews(last:20) {
        nodes { author { login } submittedAt body }
      }
      reviewThreads(first:100, after:$endCursor) {
        pageInfo { hasNextPage endCursor }
        nodes {
          isResolved isOutdated path line
          comments(first:1) { nodes { author { login } body url } }
        }
      }
    }
  }
}' | jq -s '{data:{repository:{pullRequest:
      (.[0].data.repository.pullRequest
       + {reviewThreads:{nodes:
           [.[].data.repository.pullRequest.reviewThreads.nodes[]]}})}}}')"
fi

PRJ="$(jq -r '.data.repository.pullRequest' <<<"$DATA")"
HEAD="$(jq -r '.headRefOid' <<<"$PRJ")"

printf '\n  PR #%s — %s\n' "$PR" "$(jq -r '.title' <<<"$PRJ")"
printf '  head %s\n' "${HEAD:0:10}"
printf '  %s\n' "────────────────────────────────────────────────────────────"

fail=0

# ---- 1. has the current head actually been reviewed? -----------------------
# Codex states "Reviewed commit: `<sha>`" in its review body. A review of an
# older commit is not evidence about this one.
reviewed="$(jq -r '
  [ .reviews.nodes[]
    | select(.author.login == "chatgpt-codex-connector")
    | .body ] | last // ""
  | capture("Reviewed commit:\\*{0,2}\\s*`?(?<sha>[0-9a-f]{7,40})`?").sha // ""
' <<<"$PRJ" 2>/dev/null || echo "")"

if [[ -z "$reviewed" ]]; then
  printf '  NO CODEX REVIEW FOUND on this pull request.\n'
  printf '  A merge here is unreviewed, not approved.\n'
  fail=1
elif [[ "${HEAD:0:${#reviewed}}" != "$reviewed" ]]; then
  printf '  STALE REVIEW: newest Codex review covers %s, head is %s.\n' \
         "${reviewed:0:10}" "${HEAD:0:10}"
  printf '  Its findings describe code that has since changed.\n'
  fail=1
else
  printf '  Codex reviewed %s — current.\n' "${reviewed:0:10}"
fi

# ---- 2. unresolved threads, by severity ------------------------------------
# Severity comes from the badge Codex puts in the comment: ![P1 Badge](...).
# An unresolved thread with no badge is still reported, at UNRANKED, because a
# human comment nobody answered is exactly as unread as a P1 nobody answered.
printf '\n'
threads="$(jq -r '
  .reviewThreads.nodes[]
  | select(.isResolved | not)
  | . as $t
  | ($t.comments.nodes[0] // {})   as $c
  | ($c.body // "")                as $b
  | ( if   ($b | test("P1 Badge")) then "P1"
      elif ($b | test("P2 Badge")) then "P2"
      elif ($b | test("P3 Badge")) then "P3"
      else "UNRANKED" end )        as $sev
  | [ $sev,
      ($t.path // "?"),
      ($t.line // 0 | tostring),
      ($t.isOutdated | tostring),
      ($c.author.login // "?"),
      ( $b | gsub("\\*\\*|<sub>|</sub>|!\\[[^]]*\\]\\([^)]*\\)"; "")
           | gsub("\\s+"; " ") | ltrimstr(" ") | .[0:96] )
    ] | @tsv
' <<<"$PRJ")"

if [[ -z "$threads" ]]; then
  printf '  No unresolved review threads.\n'
else
  while IFS=$'\t' read -r sev path line outdated who summary; do
    [[ -z "$sev" ]] && continue
    mark=""
    [[ "$outdated" == "true" ]] && mark=" (on an outdated diff)"
    printf '  OPEN  %-8s %s:%s%s\n        %s\n        — %s\n\n' \
           "$sev" "$path" "$line" "$mark" "$summary" "$who"
    # EVERY severity blocks, including P3 and unranked. The first version
    # listed P1|P2|UNRANKED and let a P3-only pull request print "Clear to
    # merge" while contradicting this script's own stated contract -- a check
    # that could not fail, in the check written to stop checks that cannot
    # fail. Severity ranks what to fix first; it does not rank what may be
    # ignored, and an unread P3 is exactly as unread as an unread P1.
    fail=1
  done <<<"$threads"
fi

printf '  %s\n' "────────────────────────────────────────────────────────────"
if (( fail )); then
  printf '  REFUSING MERGE.\n'
  printf '  Resolve each thread on GitHub once it is genuinely addressed —\n'
  printf '  resolving is the record that someone read it.\n\n'
  exit 1
fi
printf '  Clear to merge.\n\n'
