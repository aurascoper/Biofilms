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
verdict "a dismissed review is not coverage" 1 $? "NO CODEX REVIEW FOUND"
# The control: the identical review, not dismissed, covers the head.
comment_fixture "$QUIET" "[{\"author\":{\"login\":$CODEX},\"submittedAt\":\"2026-08-17T18:01:00Z\",\"state\":\"COMMENTED\",\"body\":$REVIEW_OF_HEAD}]"
verdict "the same review, not dismissed, is coverage" 0 $? ""

printf '  %s\n' "──────────────────────────────────────────"
if (( fail )); then printf '  FAILED\n\n'; exit 1; fi
printf '  all pass\n\n'
