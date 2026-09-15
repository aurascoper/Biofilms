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
# It also refuses a STALE review: if the newest review coverage names a commit
# that is no longer the head, the findings describe code that has since been
# rewritten, and their silence means nothing. Push, wait for the re-review,
# run this again. What counts as coverage, and until when, is block 1 below and
# docs/review_coverage_decision.md.
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
        nodes { author { login } authorAssociation submittedAt state body commit { oid } }
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

# COMMENTS GET THEIR OWN PAGINATED CALL. `gh api --paginate` advances exactly
# one connection -- the one carrying `after:$endCursor` and `pageInfo` -- and
# that is reviewThreads above. A bare `comments(last:30)` alongside it is not
# paginated at all, so a comment-form finding followed by more than thirty
# comments simply vanished and the gate cleared without a push. That is the
# same fail-open as the unpaginated reviewThreads, reproduced in the surface
# added to fix a different fail-open.
COMMENTS="$(gh api graphql --paginate \
  -f owner="$OWNER" -f name="$NAME" -F pr="$PR" -f query='
query($owner:String!, $name:String!, $pr:Int!, $endCursor:String) {
  repository(owner:$owner, name:$name) {
    pullRequest(number:$pr) {
      comments(first:100, after:$endCursor) {
        pageInfo { hasNextPage endCursor }
        nodes { author { login } authorAssociation createdAt body }
      }
    }
  }
}' | jq -s '[.[].data.repository.pullRequest.comments.nodes[]]')"
fi

PRJ="$(jq -r '.data.repository.pullRequest' <<<"$DATA")"
# Fold the separately-paginated comments in, so everything below sees one shape.
if [[ -n "${PREFLIGHT_FIXTURE:-}" ]]; then
  PRJ="$(jq '.comments.nodes //= []' <<<"$PRJ")"
else
  PRJ="$(jq --argjson c "${COMMENTS:-[]}" '.comments = {nodes: $c}' <<<"$PRJ")"
fi
HEAD="$(jq -r '.headRefOid' <<<"$PRJ")"

printf '\n  PR #%s — %s\n' "$PR" "$(jq -r '.title' <<<"$PRJ")"
printf '  head %s\n' "${HEAD:0:10}"
printf '  %s\n' "────────────────────────────────────────────────────────────"

fail=0

# ---- 1. has the current head actually been reviewed? -----------------------
# Codex states "Reviewed commit: `<sha>`" in its review body. A review of an
# older commit is not evidence about this one.
# CODEX POSTS FINDINGS TWO WAYS and this gate only read one of them. A formal
# review states "Reviewed commit: `<sha>`"; a COMMENT-form review carries the
# sha inside a blob permalink instead. When the newest pass arrived as a
# comment, the gate reported a stale review and "no unresolved threads" while a
# live P2 finding sat in that comment -- a fail-open, and the fifth one found
# in this script. Take whichever surface is NEWER.
# A DISMISSED review was withdrawn from the record by a maintainer. It is not
# coverage, and the finding scan below does not read it either.
#
# THREE SOURCES, AND ONE OF THEM EXPIRES. Codex declined on its usage limit from
# 2026-09-15 and Copilot's automatic review went silent on 2026-09-13, so a gate
# that counted Codex alone refused every head for a reason that says nothing
# about whether anyone read the code. docs/review_coverage_decision.md records
# what was accepted before, what is accepted now, and why.
#   Codex        unchanged.
#   Copilot      a review by copilot-pull-request-reviewer; its commit is the
#                review's own.
#   self-review  a comment or review body by an owner, member or collaborator
#                carrying a `review-coverage: <sha>` line and its findings, each
#                with a disposition. Counted until SELF_REVIEW_EXPIRES, refused
#                by name after it: a lowering with no date is permanent.
# A candidate that names no commit is coverage from no source. The newest one
# that names a commit decides.
CODEX_LOGIN="chatgpt-codex-connector"
COPILOT_LOGIN="copilot-pull-request-reviewer"
SELF_REVIEW_EXPIRES="2026-10-15"
TODAY="${PREFLIGHT_TODAY:-$(date -u +%F)}"

JQ_DEFS='
def entries:
  [ (.reviews.nodes // [])[] | select((.state // "") != "DISMISSED")
    | {kind: "review body", at: .submittedAt, body: (.body // ""),
       login: (.author.login // ""), assoc: (.authorAssociation // ""),
       oid: (.commit.oid // "")} ]
  + [ (.comments.nodes // [])[]
    | {kind: "comment", at: .createdAt, body: (.body // ""),
       login: (.author.login // ""), assoc: (.authorAssociation // ""), oid: ""} ];

def codex_sha:
  (.body | capture("Reviewed commit:\\*{0,2}\\s*`?(?<s>[0-9a-f]{7,40})`?").s)
  // (.body | capture("/blob/(?<s>[0-9a-f]{7,40})/").s)
  // "";

def may_self_review:
  .login != $codex and .login != $copilot
  and (.assoc | IN("OWNER", "MEMBER", "COLLABORATOR"));

# null when the body carries no review-coverage line; otherwise the commit it
# names, its findings, and everything wrong with it.
def coverage:
  (.body | split("\n") | map(sub("\r$"; ""))) as $l
  | ([ $l | to_entries[] | select(.value | test("^review-coverage:")) | .key ] | first) as $m
  | if $m == null then null else
      (($l[$m] | capture("^review-coverage:\\s*(?<s>[0-9a-f]{7,40})\\s*$").s) // "") as $sha
    | ([ $l | to_entries[] | select(.value | test("^findings:")) | .key ] | first) as $h
    | (if $h == null then [] else [ $l[($h + 1):][] | select(startswith("- ")) ] end) as $items
    | ($h != null and ($l[$h] | test("^findings:\\s*none\\s*$"))) as $none
    | ([ $l[] | select(test("^scope:\\s*[^;\\s][^;]*;\\s*checklist:\\s*\\S")) ] | length > 0) as $scope
    | [ $items[]
        | (capture("^- (?<sev>P[0-9]) (?<disp>open|fixed|deferred)(?<ref>.*?):\\s+(?<title>\\S.*)$")
           | .ref |= ((. // "") | gsub("^\\s+|\\s+$"; "")))
          // {bad: .} ] as $parsed
    | { sha: $sha,
        findings: [ $parsed[] | select(has("bad") | not) ],
        errors: [
          (if $sha == "" then "its review-coverage line names no commit (7 to 40 hex digits)" else empty end),
          (if $h == null then "it has no findings: section" else empty end),
          (if $h != null and ($none | not) and ($items | length) == 0
             then "its findings: section lists nothing; write findings: none beside a scope: line"
             else empty end),
          (if $none and ($items | length) > 0 then "it says findings: none and then lists findings" else empty end),
          (if $none and ($scope | not)
             then "findings: none needs a scope: line naming files and a checklist" else empty end),
          ($parsed[] | select(has("bad")) | "a finding line does not parse: \(.bad)"),
          ($parsed[] | select(has("bad") | not) | select(.disp != "open" and .ref == "")
             | "\(.sev) \(.disp) does not say where: \(.title)")
        ] }
    end;

def self_reviews:
  [ entries[] | select(may_self_review) | . as $e | coverage as $c
    | select($c != null) | $e + {cov: $c} ];
'
jqc() {  # every program here shares JQ_DEFS and its arguments
  jq -r --arg codex "$CODEX_LOGIN" --arg copilot "$COPILOT_LOGIN" \
        --arg today "$TODAY" --arg expires "$SELF_REVIEW_EXPIRES" \
        --arg head "$HEAD" "$JQ_DEFS$1" <<<"$PRJ"
}

# NO 2>/dev/null. A coverage query that cannot run is an unknown, not "no
# review" and not "current"; the same reasoning as block 1b.
if ! cov="$(jqc '
  ( [ entries[] | select(.login == $codex) | {source: "Codex", at, sha: codex_sha} ]
  + [ entries[] | select(.login == $copilot and .kind == "review body")
      | {source: "Copilot", at, sha: .oid} ]
  + (if $today <= $expires then
       [ self_reviews[] | select((.cov.errors | length) == 0)
         | {source: "self-review by \(.login)", at, sha: .cov.sha} ]
     else [] end) )
  | map(select(.sha != "")) | sort_by(.at) | (last // {source: "", sha: ""})
  | "\(.source)\t\(.sha)"
')"; then
  printf '  COULD NOT READ REVIEW COVERAGE: the query failed.\n'
  printf '  Treating that as unknown, not as clean.\n'
  exit 1
fi
source="${cov%%$'\t'*}"
reviewed="${cov#*$'\t'}"

# What was NOT counted, and why, is printed rather than silently dropped. A
# malformed coverage comment that names the head fails the gate: it is a claim of
# review that cannot be read, and an unreadable claim is not a clean one.
if ! notices="$(jqc '
  ( entries[] | select(.login != $codex and .login != $copilot)
    | select(may_self_review | not) | select(coverage != null)
    | "NOT COUNTED\t\(.kind) by \(.login) (author association \(if .assoc == "" then "unknown" else .assoc end))" ),
  ( if $today > $expires and (self_reviews | length) > 0
      then "EXPIRED\t\(self_reviews | length)" else empty end ),
  ( if $today <= $expires then
      self_reviews[] | select((.cov.errors | length) > 0)
      | .cov.sha as $s | select($s == "" or ($head | startswith($s)))
      | "MALFORMED\t\(.kind) by \(.login)\t\(.cov.errors | join("; "))"
    else empty end )
')"; then
  printf '  COULD NOT READ REVIEW COVERAGE COMMENTS: the query failed.\n'
  printf '  Treating that as unknown, not as clean.\n'
  exit 1
fi
while IFS=$'\t' read -r what who detail; do
  case "$what" in
    "NOT COUNTED")
      printf '  NOT COUNTED: a review-coverage %s.\n' "$who"
      printf '  Only an owner, member or collaborator can cover a head.\n' ;;
    EXPIRED)
      printf '  SELF-REVIEW COVERAGE EXPIRED %s: %s review-coverage comment(s) not counted.\n' \
             "$SELF_REVIEW_EXPIRES" "$who"
      printf '  See docs/review_coverage_decision.md.\n' ;;
    MALFORMED)
      printf '  MALFORMED COVERAGE COMMENT (%s): %s\n' "$who" "$detail"
      fail=1 ;;
  esac
done <<<"$notices"

if [[ -z "$reviewed" ]]; then
  printf '  NO REVIEW COVERAGE FOUND on this pull request: no Codex review, Copilot\n'
  printf '  review or self-review coverage comment names a commit.\n'
  printf '  A merge here is unreviewed, not approved.\n'
  fail=1
elif [[ "${HEAD:0:${#reviewed}}" != "$reviewed" ]]; then
  printf '  STALE REVIEW: the newest coverage (%s) covers %s, head is %s.\n' \
         "$source" "${reviewed:0:10}" "${HEAD:0:10}"
  printf '  Its findings describe code that has since changed.\n'
  fail=1
else
  printf '  %s reviewed %s — current.\n' "$source" "${reviewed:0:10}"
fi

# ---- 1b. findings outside review threads ------------------------------------
# These carry no resolved state, so they block only while they name the CURRENT
# head: push a fix and they describe older code, exactly as a stale review does.
# That terminates, and it cannot be cleared by ignoring it.
#
# A BADGE IS A FINDING; THE TITLE IS ONLY ITS LABEL. This block used to emit one
# row per `</sub></sub>` title, required a `#L` line anchor, and read comments
# only. A badge without that markup, a file-level permalink, or a finding in a
# formal REVIEW BODY therefore produced no row, and the gate printed "Clear to
# merge" over it. The review surface is not hypothetical: Codex put a P2 in the
# body of its review of 675d332 on #12, with no thread for it. Every badge on a
# Codex body that names the head is now a row, and a title that cannot be parsed
# is reported as unparsed rather than dropped.
#
# NO 2>/dev/null HERE. Swallowing a jq error turns a broken query into a clean
# bill of health: the first version of this block mis-scoped `.` inside
# startswith(), jq failed on every comment, the error went to /dev/null and the
# gate printed "Clear to merge" over a live P2. A query that cannot run is an
# unknown, not an all-clear, so the failure is fatal here.
if ! outside_findings="$(jq -r --arg head "$HEAD" '
  ( [ .reviews.nodes[]
      | select(.author.login == "chatgpt-codex-connector")
      | select((.state // "") != "DISMISSED")
      | {kind: "review body", body: (.body // "")} ]
    + [ (.comments.nodes // [])[]
      | select(.author.login == "chatgpt-codex-connector")
      | {kind: "comment", body: (.body // "")} ] )[]
  | .kind as $k
  | .body as $b
  | ( ($b | capture("Reviewed commit:\\*{0,2}\\s*`?(?<s>[0-9a-f]{7,40})`?").s)
      // ($b | capture("/blob/(?<s>[0-9a-f]{7,40})/").s)
      // "" ) as $sha
  | select($sha != "" and ($head | startswith($sha)))
  | ([ $b | scan("badge/(P[0-9])-") ] | flatten) as $sev
  | ([ $b | scan("</sub></sub>\\s*([^*\n]+)") ] | flatten) as $titles
  | range(0; $sev | length) as $i
  | $sev[$i] + "\t" + $k + "\t"
    + ( ($titles[$i] // "") | sub("^\\s+";"") | sub("\\s+$";"")
        | if . == "" then "(title not parsed; read the \($k) itself)" else . end )
' <<<"$PRJ")"; then
  printf '  COULD NOT READ FINDINGS OUTSIDE THREADS: the query failed.\n'
  printf '  Treating that as unknown, not as clean.\n'
  exit 1
fi

if [[ -n "$outside_findings" ]]; then
  printf '\n  CODEX POSTED FINDINGS OUTSIDE REVIEW THREADS on this head, in a\n'
  printf '  review body or a comment. They have no resolve button; fix and push.\n\n'
  while IFS=$'\t' read -r sev kind title; do
    [[ -z "$sev" ]] && continue
    printf '  OPEN  %-8s (%s, no thread) %s\n' "$sev" "$kind" "$title"
    fail=1
  done <<<"$outside_findings"
fi

# ---- 1c. findings a self-review coverage comment lists for this head --------
# Every finding carries its disposition. `open` blocks until a push moves the
# head, exactly like a Codex finding outside a thread; `fixed` and `deferred`
# are the reviewer's record and are printed, not enforced. Every well-formed
# coverage comment that names the head is read, so a disposition is changed by
# editing the comment.
if [[ ! "$TODAY" > "$SELF_REVIEW_EXPIRES" ]]; then
  if ! listed="$(jqc '
    self_reviews[] | select((.cov.errors | length) == 0)
    | .cov.sha as $s | select($head | startswith($s))
    | .cov.findings[] | [.sev, .disp, .title, .ref] | @tsv
  ')"; then
    printf '  COULD NOT READ REVIEW-COVERAGE FINDINGS: the query failed.\n'
    printf '  Treating that as unknown, not as clean.\n'
    exit 1
  fi
  if [[ -n "$listed" ]]; then
    printf '\n  REVIEW-COVERAGE FINDINGS on this head:\n\n'
    # ref is LAST: tab is whitespace to `read`, so an empty field in the middle
    # would collapse and shift the title into it.
    while IFS=$'\t' read -r sev disp title ref; do
      [[ -z "$sev" ]] && continue
      if [[ "$disp" == open ]]; then
        printf '  OPEN  %-8s (self-review, no thread) %s\n' "$sev" "$title"
        fail=1
      else
        printf '  %-8s %-3s %s — %s\n' "$disp" "$sev" "$ref" "$title"
      fi
    done <<<"$listed"
  fi
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
