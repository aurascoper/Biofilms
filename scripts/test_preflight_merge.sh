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

printf '  %s\n' "──────────────────────────────────────────"
if (( fail )); then printf '  FAILED\n\n'; exit 1; fi
printf '  all pass\n\n'
