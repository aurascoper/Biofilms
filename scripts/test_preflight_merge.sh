#!/usr/bin/env bash
# Negative controls for the merge gate.
#
# The gate exists because a check that cannot fail is not a check -- and its
# first version was one. It listed `P1|P2|UNRANKED` as the failing cases and let
# a P3-only pull request print "Clear to merge" while its own documentation said
# it refuses while ANY thread is unresolved. Caught by external review.
#
# So the gate gets what it demands of everything else: inputs that must be
# rejected, and one that must be accepted so the whole thing cannot pass by
# refusing everything.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GATE="$HERE/preflight_merge.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

HEAD_SHA="abc1234def5678"
fail=0

# $1 name, $2 expected exit, $3 threads JSON array, $4 reviewed-commit sha
check() {
  local name="$1" want="$2" threads="$3" reviewed="$4"
  cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"fixture","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
    "submittedAt":"2026-08-16T00:00:00Z",
    "body":"### Codex Review\n\n**Reviewed commit:** \`$reviewed\`\n"}]},
  "comments":{"nodes":[]},
  "reviewThreads":{"nodes":$threads}
}}}}
JSON
  PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1
  local got=$?
  if [[ "$got" == "$want" ]]; then
    printf '  ok    %s (exit %s)\n' "$name" "$got"
  else
    printf '  FAIL  %s: expected exit %s, got %s\n' "$name" "$want" "$got"
    sed 's/^/        | /' "$TMP/out"
    fail=1
  fi
}

thread() {  # $1 severity-badge text, $2 isResolved
  printf '[{"isResolved":%s,"isOutdated":false,"path":"a/b.py","line":7,
           "comments":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
           "body":"![%s Badge](x) something","url":"u"}]}}]' "$2" "$1"
}

printf '\n  merge gate — negative controls\n  %s\n' \
       "──────────────────────────────────────────"

# THE ONE THAT WAS BROKEN. A P3 is a finding; nobody read it either.
check "an unresolved P3 blocks"        1 "$(thread P3 false)" "$HEAD_SHA"
check "an unresolved P1 blocks"        1 "$(thread P1 false)" "$HEAD_SHA"
check "an unresolved P2 blocks"        1 "$(thread P2 false)" "$HEAD_SHA"
# A human comment with no badge is exactly as unread as a badged one.
check "an unranked thread blocks"      1 "$(thread XX false)" "$HEAD_SHA"
# A review of a commit that is no longer the head says nothing about the head.
check "a stale review blocks"          1 "[]"                 "0000000000"
# And no review at all is not the same as a clean review.
check "no codex review blocks"         1 "[]"                 ""

# THE CONTROL THAT KEEPS THE OTHERS HONEST: a gate that refuses everything
# would pass every test above and be useless.
check "resolved threads clear"         0 "$(thread P1 true)"  "$HEAD_SHA"
check "no threads at all clears"       0 "[]"                 "$HEAD_SHA"

# THE COMMENT-FORM REVIEW. Codex sometimes posts its pass as a comment whose
# only commit reference is a blob permalink, with no "Reviewed commit:" text.
# jq's `capture` yields empty on no match, so the `//` fallback is evaluated
# and the permalink sha is read -- a review claimed otherwise, and this pins
# the behaviour it questioned rather than leaving it to the reader of the jq.
# $1 name, $2 expected exit, $3 the sha in the permalink
check_comment() {
  local name="$1" want="$2" sha="$3"
  cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"fixture","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":[]},
  "comments":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
    "createdAt":"2026-08-16T00:00:00Z",
    "body":"Codex Review: Didn't find any major issues. https://github.com/o/r/blob/$sha/a/b.py#L7"}]},
  "reviewThreads":{"nodes":[]}
}}}}
JSON
  PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1
  local got=$?
  if [[ "$got" == "$want" ]]; then
    printf '  ok    %s (exit %s)\n' "$name" "$got"
  else
    printf '  FAIL  %s: expected exit %s, got %s\n' "$name" "$want" "$got"
    sed 's/^/        | /' "$TMP/out"
    fail=1
  fi
}
check_comment "a permalink-only comment review of the head clears" 0 "$HEAD_SHA"
check_comment "a permalink-only comment review of another commit blocks" 1 "0000000000"

# ---------------------------------------------------------------------------
# The pagination path, which PREFLIGHT_FIXTURE bypasses entirely.
#
# Every control above supplies an already-flattened response, so removing
# --paginate, $endCursor, pageInfo or the jq fold left them all green -- the
# fixture hook that made the gate testable also hid the newest thing added to
# it. So: a fake `gh` on PATH that emits TWO pages, with the only unresolved
# thread on page two. A gate that reads one page reports "Clear to merge".
printf '\n  pagination — the fixture cannot reach this\n  %s\n' \
       "──────────────────────────────────────────"

FAKEBIN="$TMP/bin"; mkdir -p "$FAKEBIN"
cat >"$FAKEBIN/gh" <<'FAKE'
#!/usr/bin/env bash
# `gh repo view` -> the repo; `gh api graphql --paginate` -> two pages, exactly
# as the real command emits them: one JSON document per page, concatenated.
if [[ "$1" == "repo" ]]; then echo "aurascoper/Biofilms"; exit 0; fi
# The gate makes TWO paginated calls -- reviewThreads and comments. The
# comments one gets TWO pages with the finding only on page two, and page two
# is emitted only when --paginate is present AND the query carries the cursor
# contract. A fixture holding all the comments at once proves nothing about
# pagination; this is the only path that does.
if [[ "$*" == *"comments(first:100"* ]]; then
  cpage() { cat <<CJSON
{"data":{"repository":{"pullRequest":{"comments":{
  "pageInfo":{"hasNextPage":$1,"endCursor":"CC"},"nodes":$2}}}}}
CJSON
  }
  cpage true '[{"author":{"login":"aurascoper"},"createdAt":"2026-08-17T19:00:00Z","body":"page one chatter"}]'
  if [[ "$*" == *--paginate* ]] \
     && [[ "$*" == *'$endCursor: String'* || "$*" == *'$endCursor:String'* ]] \
     && [[ "$*" == *'after: $endCursor'* || "$*" == *'after:$endCursor'* ]] \
     && [[ "$*" == *hasNextPage* ]]; then
    cpage false '[{"author":{"login":"chatgpt-codex-connector"},"createdAt":"2026-08-17T19:05:00Z","body":"### Codex Review\n\nhttps://github.com/o/r/blob/abc1234def5678/f.py#L9-L9\n**<sub><sub>![P1 Badge](https://img.shields.io/badge/P1-orange?style=flat)</sub></sub>  Finding only on comment page two**\n"}]'
  fi
  exit 0
fi
page() {  # $1 hasNextPage, $2 threads-json
  cat <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"paged fixture","headRefOid":"abc1234def5678",
  "reviews":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
    "submittedAt":"2026-08-16T00:00:00Z",
    "body":"**Reviewed commit:** \`abc1234def5678\`"}]},
  "reviewThreads":{"pageInfo":{"hasNextPage":$1,"endCursor":"CUR"},
                   "nodes":$2}}}}}
JSON
}
# Page one: nothing but a RESOLVED thread. Page two: the open P1.
page true  '[{"isResolved":true,"isOutdated":false,"path":"p1.py","line":1,
             "comments":{"nodes":[{"author":{"login":"x"},"body":"resolved","url":"u"}]}}]'
# THE FAKE ENFORCES THE REAL CONTRACT, or it tests nothing. `gh api --paginate`
# only walks a GraphQL connection when the query declares an `$endCursor:
# String` variable AND selects `pageInfo { hasNextPage endCursor }`. A fake
# that emitted page two on the flag alone would stay green after either element
# was deleted, while real pagination silently stopped -- which is the same
# fail-open shape the pagination fix existed to close.
# `after:$endCursor` ON THE CONNECTION, not merely the variable declaration.
# Checking that `$endCursor` appears somewhere is satisfied by the declaration
# alone: delete the binding and cursor updates cannot advance reviewThreads,
# so the real gate sticks on page one while this control stays green. That was
# verified by removing only the binding -- all nine controls passed.
if [[ "$*" == *--paginate* ]] \
   && [[ "$*" == *'$endCursor: String'* || "$*" == *'$endCursor:String'* ]] \
   && [[ "$*" == *'after: $endCursor'* || "$*" == *'after:$endCursor'* ]] \
   && [[ "$*" == *hasNextPage* ]] && [[ "$*" == *endCursor* ]]; then
  page false '[{"isResolved":false,"isOutdated":false,"path":"page2.py","line":9,
                "comments":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
                "body":"![P1 Badge](x) only visible on page two","url":"u"}]}}]'
fi
FAKE
chmod +x "$FAKEBIN/gh"

out="$(PATH="$FAKEBIN:$PATH" "$GATE" 1 2>&1)"; got=$?
if [[ "$got" == 1 ]] && grep -q "page2.py" <<<"$out"; then
  printf '  ok    a thread on page two still blocks (exit 1)\n'
else
  printf '  FAIL  page-two thread not seen: exit %s\n' "$got"
  sed 's/^/        | /' <<<"$out"
  fail=1
fi

# ---------------------------------------------------------------------------
# COMMENT-FORM REVIEWS. Codex posts findings two ways: a formal review, and an
# issue comment carrying the sha in a blob permalink. The gate read only
# reviews, so when the newest pass arrived as a comment it reported a stale
# review and "no unresolved threads" -- while a live P2 sat in that comment.
printf '\n  comment-form findings — the surface the gate was blind to\n  %s\n' \
       "──────────────────────────────────────────"

BODY='### Codex Review\n\nhttps://github.com/o/r/blob/'"$HEAD_SHA"'/f.py#L505-L507\n**<sub><sub>![P2 Badge](https://img.shields.io/badge/P2-yellow?style=flat)</sub></sub>  Account for negation before the positive subject**\n\nDetail.'

comment_fixture() {  # $1 = comment body JSON string, $2 = reviews nodes
  cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"fixture","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":$2},
  "comments":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
    "createdAt":"2026-08-17T18:00:00Z","body":$1}]},
  "reviewThreads":{"nodes":[]}
}}}}
JSON
  PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1
}

# A comment naming the CURRENT head, with a finding -> must block, and must
# name the finding rather than merely refusing.
comment_fixture "$(jq -Rn --arg b "$(printf '%b' "$BODY")" '$b')" "[]"
got=$?
if [[ "$got" == 1 ]] && grep -q "Account for negation" "$TMP/out"; then
  printf '  ok    a comment-form finding on the head blocks, and is named\n'
else
  printf '  FAIL  comment-form finding not surfaced (exit %s)\n' "$got"
  sed 's/^/        | /' "$TMP/out"; fail=1
fi

# THE SAME COMMENT ON AN OLD SHA must NOT block: pushing a fix is what clears
# it, exactly as with a stale review. Otherwise the gate never terminates.
OLDBODY=${BODY//$HEAD_SHA/0000000000}
comment_fixture "$(jq -Rn --arg b "$(printf '%b' "$OLDBODY")" '$b')" \
  '[{"author":{"login":"chatgpt-codex-connector"},"submittedAt":"2026-08-17T18:01:00Z","body":"**Reviewed commit:** `'"$HEAD_SHA"'`"}]'
got=$?
if [[ "$got" == 0 ]]; then
  printf '  ok    the same finding on an older sha no longer blocks\n'
else
  printf '  FAIL  a superseded comment still blocks (exit %s)\n' "$got"
  sed 's/^/        | /' "$TMP/out"; fail=1
fi

# STALENESS MUST READ THE COMMENT TOO. A comment-only pass on the current head
# is not a stale review; before the fix the gate called it one.
comment_fixture "$(jq -Rn --arg b "$(printf '%b' "$BODY")" '$b')" "[]"
if grep -q "STALE REVIEW" "$TMP/out"; then
  printf '  FAIL  a current comment-form review was reported as stale\n'
  sed 's/^/        | /' "$TMP/out"; fail=1
else
  printf '  ok    a comment-form review counts as covering the head\n'
fi

# THE FINDING BURIED PAST THE FIRST PAGE. `comments(last:30)` in the main query
# was not paginated at all -- one --paginate call advances one connection, and
# that one is reviewThreads. A comment-form finding followed by thirty later
# comments vanished, and the gate cleared without a push.
printf '\n  comment pagination — the finding buried under later chatter\n  %s\n' \
       "──────────────────────────────────────────"

# 40 innocuous comments AFTER the finding, so anything reading only a recent
# window misses it entirely.
FILLER="$(python3 - <<'PY'
import json
print(",".join(json.dumps({"author":{"login":"aurascoper"},
                           "createdAt":"2026-08-17T19:%02d:00Z" % i,
                           "body":"routine follow-up %d" % i}) for i in range(40)))
PY
)"
FINDING="$(python3 - "$HEAD_SHA" <<'PY'
import json,sys
sha=sys.argv[1]
print(json.dumps({"author":{"login":"chatgpt-codex-connector"},
  "createdAt":"2026-08-17T18:00:00Z",
  "body":"### Codex Review\n\nhttps://github.com/o/r/blob/%s/f.py#L1-L2\n"
         "**<sub><sub>![P2 Badge](https://img.shields.io/badge/P2-yellow?style=flat)"
         "</sub></sub>  Buried finding past the first page**\n" % sha}))
PY
)"
cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"paged comments","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":[{"author":{"login":"chatgpt-codex-connector"},
    "submittedAt":"2026-08-17T18:30:00Z",
    "body":"**Reviewed commit:** \`$HEAD_SHA\`"}]},
  "comments":{"nodes":[$FINDING,$FILLER]},
  "reviewThreads":{"nodes":[]}
}}}}
JSON
PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1; got=$?
if [[ "$got" == 1 ]] && grep -q "Buried finding" "$TMP/out"; then
  printf '  ok    a finding behind 40 later comments still blocks\n'
else
  printf '  FAIL  buried comment-form finding was missed (exit %s)\n' "$got"
  sed 's/^/        | /' "$TMP/out"; fail=1
fi

# THROUGH THE REAL QUERY, not a grep over it. The fake serves the finding only
# on comment page TWO, and only when --paginate and the cursor contract are
# both present -- so removing either leaves the gate seeing page one alone.
out="$(PATH="$FAKEBIN:$PATH" "$GATE" 1 2>&1)"; got=$?
if [[ "$got" == 1 ]] && grep -q "Finding only on comment page two" <<<"$out"; then
  printf '  ok    a finding on comment page two blocks (real --paginate path)\n'
else
  printf '  FAIL  comment page two never fetched: exit %s\n' "$got"
  sed 's/^/        | /' <<<"$out"; fail=1
fi

# ---------------------------------------------------------------------------
# FINDINGS OUTSIDE THREADS, IN EVERY SHAPE CODEX USES. The scanner keyed on the
# `</sub></sub>` title markup and a `#L` anchor, and never read review bodies:
# the first three shapes below each printed "Clear to merge" before the fix.
# The review-body fixture is the shape of Codex's review of 675d332 on #12,
# which carried a P2 in its body with no thread, pointed at the head.
printf '\n  findings outside threads — every shape Codex uses\n  %s\n' \
       "──────────────────────────────────────────"

json_str() { python3 -c 'import json,sys; print(json.dumps(sys.argv[1]))' "$1"; }

verdict() {  # $1 name, $2 expected exit, $3 actual exit, $4 text the output must carry ("" for none)
  if [[ "$3" == "$2" ]] && { [[ -z "$4" ]] || grep -qF -- "$4" "$TMP/out"; }; then
    printf '  ok    %s (exit %s)\n' "$1" "$3"
  else
    printf '  FAIL  %s: expected exit %s%s, got %s\n' "$1" "$2" "${4:+ naming \"$4\"}" "$3"
    sed 's/^/        | /' "$TMP/out"; fail=1
  fi
}

CODEX='"chatgpt-codex-connector"'
QUIET='"routine chatter, no commit named"'

# A badge with no `<sub>` markup around it: zero titles used to mean zero rows.
BODY_NOSUB="$(json_str "$(printf 'https://github.com/o/r/blob/%s/f.py#L5-L6\n![P1 Badge](https://img.shields.io/badge/P1-orange?style=flat) Live P1 without the sub markup' "$HEAD_SHA")")"
comment_fixture "$BODY_NOSUB" "[]"
verdict "a badge without the title markup blocks, reported as unparsed" 1 $? "title not parsed"

# A file-level permalink with no `#L` anchor: the select() used to drop it.
BODY_NOL="$(json_str "$(printf 'https://github.com/o/r/blob/%s/f.py\n**<sub><sub>![P2 Badge](https://img.shields.io/badge/P2-yellow?style=flat)</sub></sub>  File-level finding with no line anchor**' "$HEAD_SHA")")"
comment_fixture "$BODY_NOL" "[]"
verdict "a finding on a file-level permalink blocks" 1 $? "File-level finding with no line anchor"

# A finding in a formal review body, with no thread and no comment finding.
BODY_REVIEW="$(json_str "$(printf '### Codex Review\n\nhttps://github.com/o/r/blob/%s/viewer.py#L315\n**<sub><sub>![P2 Badge](https://img.shields.io/badge/P2-yellow?style=flat)</sub></sub>  Require reductions to be coarser on every axis**\n\nDetail.' "$HEAD_SHA")")"
comment_fixture "$QUIET" "[{\"author\":{\"login\":$CODEX},\"submittedAt\":\"2026-08-17T18:01:00Z\",\"state\":\"COMMENTED\",\"body\":$BODY_REVIEW}]"
verdict "a finding in a formal review body on the head blocks" 1 $? "Require reductions to be coarser on every axis"

# ...and it terminates: the same review body on an older sha, beside a current pass, clears.
BODY_REVIEW_OLD="${BODY_REVIEW//$HEAD_SHA/0000000000}"
comment_fixture "$(json_str "https://github.com/o/r/blob/$HEAD_SHA/a.py#L1 no findings")" "[{\"author\":{\"login\":$CODEX},\"submittedAt\":\"2026-08-17T17:00:00Z\",\"state\":\"COMMENTED\",\"body\":$BODY_REVIEW_OLD}]"
verdict "the same review-body finding on an older sha does not block" 0 $? ""

# A DISMISSED review was withdrawn from the record: it is not coverage.
REVIEW_OF_HEAD="$(json_str "**Reviewed commit:** \`$HEAD_SHA\`")"
comment_fixture "$QUIET" "[{\"author\":{\"login\":$CODEX},\"submittedAt\":\"2026-08-17T18:01:00Z\",\"state\":\"DISMISSED\",\"body\":$REVIEW_OF_HEAD}]"
verdict "a dismissed review is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"
# The control: the identical review, not dismissed, covers the head.
comment_fixture "$QUIET" "[{\"author\":{\"login\":$CODEX},\"submittedAt\":\"2026-08-17T18:01:00Z\",\"state\":\"COMMENTED\",\"body\":$REVIEW_OF_HEAD}]"
verdict "the same review, not dismissed, is coverage" 0 $? ""

# ---------------------------------------------------------------------------
# COVERAGE FROM SOURCES OTHER THAN CODEX (docs/review_coverage_decision.md).
# A self-review coverage comment is a deliberate lowering of this gate, so each
# way it can be vacuous, stale, unauthorised or past its date gets the input
# that must be refused, beside the input that must be accepted.
printf '\n  coverage other than Codex — self-review and Copilot\n  %s\n' \
       "──────────────────────────────────────────"

# $1 comment body (plain text), $2 author association, $3 reviews nodes (default [])
self_fixture() {
  local body; body="$(json_str "$1")"
  cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"fixture","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":${3:-[]}},
  "comments":{"nodes":[{"author":{"login":"aurascoper"},"authorAssociation":"$2",
    "createdAt":"2026-09-15T20:00:00Z","body":$body}]},
  "reviewThreads":{"nodes":[]}
}}}}
JSON
  PREFLIGHT_TODAY="${TODAY_OVERRIDE:-2026-09-15}" PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1
}

SCOPE='scope: scripts/preflight_merge.sh; checklist: AGENTS.md, the six rules'
COVERED="review-coverage: $HEAD_SHA
$SCOPE
findings:
- P1 fixed 71a5c90: the gate cleared a comment finding it could not parse
- P2 deferred https://github.com/o/r/pull/12#issuecomment-1: the pilot publication guard"

# The three the decision names: stale, vacuous, and the one that must pass.
self_fixture "${COVERED//$HEAD_SHA/0000000000}" OWNER
verdict "a coverage comment on a stale sha blocks" 1 $? "STALE REVIEW"
self_fixture "review-coverage: $HEAD_SHA
$SCOPE" OWNER
verdict "the marker with no findings blocks, by name" 1 $? "MALFORMED COVERAGE COMMENT"
self_fixture "$COVERED" OWNER
verdict "fixed and deferred findings on the head clear" 0 $? "self-review by aurascoper reviewed"

# An open finding blocks and is named; the disposition is what the gate enforces.
self_fixture "$COVERED
- P2 open: the finding this control names" OWNER
verdict "an open finding in a coverage comment blocks, named" 1 $? "the finding this control names"

# "No findings" is accepted only as a statement with its scope.
self_fixture "review-coverage: $HEAD_SHA
$SCOPE
findings: none" OWNER
verdict "findings: none beside a scope line clears" 0 $? ""
self_fixture "review-coverage: $HEAD_SHA
findings: none" OWNER
verdict "findings: none with no scope line blocks" 1 $? "needs a scope: line"

# A deferral that does not say where is a finding nobody will find again.
self_fixture "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P2 deferred: nothing says where this went" OWNER
verdict "a deferral that names no place blocks" 1 $? "does not say where"

# On a public repository anyone can type the marker.
self_fixture "$COVERED" NONE
verdict "a marker from outside the repository is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"

# The expiry is enforced, not recorded: the last day counts, the next does not.
TODAY_OVERRIDE=2026-10-15 self_fixture "$COVERED" OWNER
verdict "self-review coverage still counts on its expiry date" 0 $? ""
TODAY_OVERRIDE=2026-10-16 self_fixture "$COVERED" OWNER
verdict "the day after expiry it is refused by name" 1 $? "SELF-REVIEW COVERAGE EXPIRED"

# Copilot: the review's own commit is what it covers.
COPILOT_REVIEW() { printf '[{"author":{"login":"copilot-pull-request-reviewer"},"authorAssociation":"NONE","submittedAt":"2026-09-15T21:00:00Z","state":"COMMENTED","body":"## Pull request overview","commit":{"oid":"%s"}}]' "$1"; }
self_fixture "routine chatter, no marker" OWNER "$(COPILOT_REVIEW "$HEAD_SHA")"
verdict "a Copilot review of the head is coverage" 0 $? "Copilot reviewed"
self_fixture "routine chatter, no marker" OWNER "$(COPILOT_REVIEW 0000000000000000)"
verdict "a Copilot review of an older commit is stale" 1 $? "STALE REVIEW"

# ---------------------------------------------------------------------------
# THE PATHS A FIRST PASS LEFT UNCONTROLLED. Two reviews of the commit above
# found six defects inside the new code and nine paths with no control at all --
# including the `fail=1` that makes a malformed coverage comment block, which
# could be deleted with every control still green. Each one below is the input
# that fails when its line is removed.
printf '\n  coverage comments — every path, and every way to write one wrong\n  %s\n' \
       "──────────────────────────────────────────"

c_node() {  # $1 body, $2 association (OWNER), $3 login (aurascoper), $4 createdAt
  printf '{"author":{"login":"%s"},"authorAssociation":"%s","createdAt":"%s","body":%s}' \
         "${3:-aurascoper}" "${2:-OWNER}" "${4:-2026-09-15T20:00:00Z}" "$(json_str "$1")"
}
r_codex() {  # $1 sha, $2 submittedAt
  printf '{"author":{"login":"chatgpt-codex-connector"},"authorAssociation":"NONE","submittedAt":"%s","state":"COMMENTED","body":%s}' \
         "${2:-2026-09-15T18:00:00Z}" "$(json_str "**Reviewed commit:** \`$1\`")"
}
r_copilot() {  # $1 commit oid, $2 submittedAt, $3 state
  printf '{"author":{"login":"copilot-pull-request-reviewer"},"authorAssociation":"NONE","submittedAt":"%s","state":"%s","body":"## Pull request overview","commit":{"oid":"%s"}}' \
         "${2:-2026-09-15T10:00:00Z}" "${3:-COMMENTED}" "$1"
}
r_self() {  # $1 body, $2 state, $3 submittedAt as JSON (null for a draft)
  printf '{"author":{"login":"aurascoper"},"authorAssociation":"OWNER","submittedAt":%s,"state":"%s","body":%s}' \
         "${3:-\"2026-09-15T20:00:00Z\"}" "${2:-COMMENTED}" "$(json_str "$1")"
}
nodes() { local IFS=,; printf '[%s]' "$*"; }
run_gate() {  # $1 comments array, $2 reviews array
  cat >"$TMP/fixture.json" <<JSON
{"data":{"repository":{"pullRequest":{
  "title":"fixture","headRefOid":"$HEAD_SHA",
  "reviews":{"nodes":${2:-[]}},"comments":{"nodes":${1:-[]}},
  "reviewThreads":{"nodes":[]}
}}}}
JSON
  PREFLIGHT_TODAY="${TODAY_OVERRIDE:-2026-09-15}" PREFLIGHT_FIXTURE="$TMP/fixture.json" "$GATE" 1 >"$TMP/out" 2>&1
}

# THE ONE THAT PINNED NOTHING. Deleting the fail=1 behind MALFORMED left all 34
# controls green, because each took its exit from the ABSENCE of coverage. Beside
# a current Codex review there is no absence, and an unreadable claim of review
# has to refuse on its own.
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE")")" "$(nodes "$(r_codex "$HEAD_SHA")")"
verdict "a malformed comment blocks beside current Codex coverage" 1 $? "MALFORMED COVERAGE COMMENT"

# A finding line the parser cannot read is reported, never dropped.
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P1: forgot the disposition")")"
verdict "a finding line that does not parse blocks, named" 1 $? "does not parse"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P1 fixedd: a misspelling is not a disposition")")"
verdict "a misspelled disposition blocks" 1 $? "does not parse"

# Inside the section everything is a finding; below it, prose is prose.
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
see the other pull request
- P2 open: something")")"
verdict "prose inside the findings section blocks" 1 $? "is not a finding"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings: none

- see also the other pull request")")"
verdict "a bullet below the section is prose, not a finding" 0 $? ""

# EVERY MARKDOWN BULLET. Keying on "- " dropped an indented sub-item and a star
# bullet: the reviewer wrote the finding down and the gate cleared over it.
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P1 fixed 71a5c90: a real one
  - P1 open: the indented finding this control names")")"
verdict "an indented open finding blocks" 1 $? "the indented finding this control names"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
* P2 open: the starred finding this control names")")"
verdict "a starred open finding blocks" 1 $? "the starred finding this control names"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
+ P3 deferred docs/x.md: a plus bullet is a bullet")")"
verdict "a plus bullet parses as a finding" 0 $? "deferred"

# One marker per body. The first winning silently let a quoted example decide.
run_gate "$(nodes "$(c_node "For reference:
review-coverage: ffffffffffff

And the real one:
review-coverage: $HEAD_SHA
$SCOPE
findings: none")")"
verdict "two markers in one body block, by name" 1 $? "review-coverage lines"

# Quoting the format, or somebody else's comment, certifies nothing.
run_gate "$(nodes "$(c_node "The shape is:

\`\`\`
review-coverage: $HEAD_SHA
$SCOPE
findings: none
\`\`\`")")"
verdict "a marker inside a fence is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"
# ...and that has to be the MARKER's own fence check. With the findings section
# outside the fence, the body is otherwise complete, so only the marker line's
# fence awareness stands between a quoted sha and a certified head.
run_gate "$(nodes "$(c_node "Here is the marker I will use:

\`\`\`
review-coverage: $HEAD_SHA
\`\`\`

$SCOPE
findings: none")")"
verdict "a fenced marker beside an unfenced findings section is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"
run_gate "$(nodes "$(c_node "> review-coverage: $HEAD_SHA
> findings: none")")"
verdict "a quoted marker is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"

# A draft review is a review nobody submitted.
run_gate "[]" "$(nodes "$(r_self "review-coverage: $HEAD_SHA
$SCOPE
findings: none" PENDING null)")"
verdict "a PENDING draft review is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"

# Evidence about the head beats newer evidence about something else.
run_gate "[]" "$(nodes "$(r_copilot "$HEAD_SHA" 2026-09-15T10:00:00Z)" "$(r_codex 0000000000 2026-09-15T12:00:00Z)")"
verdict "coverage of the head is not masked by a newer stale review" 0 $? "Copilot reviewed"

# The deliberate change, pinned: a body naming no commit is not the newest word.
run_gate "$(nodes "$(c_node "Codex has hit its usage limit." NONE chatgpt-codex-connector 2026-09-15T19:00:00Z)")" \
         "$(nodes "$(r_codex "$HEAD_SHA" 2026-09-15T18:00:00Z)")"
verdict "a Codex comment naming no commit does not mask its own review" 0 $? "Codex reviewed"

# Block 1c terminates, exactly as block 1b does: push, and the finding goes.
run_gate "$(nodes "$(c_node "review-coverage: 0000000000
$SCOPE
findings:
- P2 open: an open finding on a commit that is gone")")" "$(nodes "$(r_codex "$HEAD_SHA")")"
verdict "an open finding on an older sha stops blocking" 0 $? ""

# The error clauses that nothing covered.
run_gate "$(nodes "$(c_node "review-coverage: nothexatall
$SCOPE
findings: none")")"
verdict "a marker naming no commit blocks, by name" 1 $? "names no commit"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:")")"
verdict "an empty findings section blocks, by name" 1 $? "lists nothing"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings: none
- P2 open: contradicting the none above")")"
verdict "none followed by findings blocks, by name" 1 $? "and then lists findings"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
scope: ; checklist: AGENTS.md
findings: none")")"
verdict "a scope line naming no files blocks" 1 $? "needs a scope: line"

# The notice is printed, not merely implied by a refusal elsewhere.
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings: none" NONE passer-by)")" "$(nodes "$(r_codex "$HEAD_SHA")")"
verdict "an outsider marker is reported as not counted" 0 $? "NOT COUNTED"

# Past the date, a self-review covers nothing and enforces nothing.
TODAY_OVERRIDE=2026-10-16 run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P2 open: an expired open finding")")" "$(nodes "$(r_codex "$HEAD_SHA")")"
verdict "after expiry a self-review finding is not enforced" 0 $? "EXPIRED"

# Correct behaviours that nothing pinned. A review withdrawn from the record is
# not coverage from either source, and a later clean comment does not cancel an
# earlier open finding: every coverage comment naming the head is read, so a
# disposition changes by editing the comment that carries it.
run_gate "[]" "$(nodes "$(r_self "review-coverage: $HEAD_SHA
$SCOPE
findings: none" DISMISSED)")"
verdict "a dismissed review carrying a marker is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"
run_gate "[]" "$(nodes "$(r_copilot "$HEAD_SHA" 2026-09-15T10:00:00Z DISMISSED)")"
verdict "a dismissed Copilot review is not coverage" 1 $? "NO REVIEW COVERAGE FOUND"
run_gate "$(nodes "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings:
- P2 open: the earlier comment still says open" OWNER aurascoper 2026-09-15T20:00:00Z)" \
          "$(c_node "review-coverage: $HEAD_SHA
$SCOPE
findings: none" OWNER aurascoper 2026-09-15T21:00:00Z)")"
verdict "a later clean comment does not cancel an earlier open finding" 1 $? "the earlier comment still says open"

printf '  %s\n' "──────────────────────────────────────────"
if (( fail )); then printf '  FAILED\n\n'; exit 1; fi
printf '  all pass\n\n'
