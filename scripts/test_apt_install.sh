#!/usr/bin/env bash
# Negative controls for scripts/apt_install.sh. No apt, no network, no root: a stub
# `sudo` on PATH records what it was asked to run and exits as each case tells it to.
#
#   scripts/test_apt_install.sh
#
# The cases, each with the result it must get:
#   update fails, install succeeds     exit 0, a warning on stderr, the install was run
#   update succeeds, install fails     non-zero
#   update fails, install fails        non-zero
#   no arguments                       exit 2, and nothing was run
#   a mutant that forgives the install  exit 0 on a failed install: the case above can
#                                      tell a broken script from a sound one
set -u -o pipefail
here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
bad=0

mkdir "$tmp/bin"
cat > "$tmp/bin/sudo" <<'EOF'
#!/usr/bin/env bash
echo "$*" >> "$STUB_LOG"
case "$2" in
    update)  exit "$STUB_UPDATE_RC" ;;
    install) exit "$STUB_INSTALL_RC" ;;
    *)       echo "stub sudo: unexpected call: $*" >&2; exit 99 ;;
esac
EOF
chmod +x "$tmp/bin/sudo"

run() {  # run <script> <update rc> <install rc> [arguments...]; sets rc, err, log
    local script=$1 u=$2 i=$3; shift 3
    : > "$tmp/log"
    err=$(PATH="$tmp/bin:$PATH" STUB_LOG="$tmp/log" STUB_UPDATE_RC=$u STUB_INSTALL_RC=$i \
          bash "$script" "$@" 2>&1 >/dev/null)
    rc=$?
    log=$(cat "$tmp/log")
}
ok()    { echo "ok     $1"; }
wrong() { echo "WRONG  $1"; bad=1; }

run "$here/apt_install.sh" 100 0 -y xvfb libgl1
if [ "$rc" = 0 ] && [[ "$err" == *warning* ]] && [[ "$log" == *"apt-get install -y xvfb libgl1"* ]]; then
    ok "a failed update warns, and the install still runs with its arguments"
else
    wrong "failed update: rc=$rc err='$err' log='$log'"
fi

run "$here/apt_install.sh" 0 100 -y no-such-package
[ "$rc" != 0 ] && ok "a failed install fails the script (rc=$rc)" || wrong "a failed install exited 0"

run "$here/apt_install.sh" 100 100 -y xvfb
[ "$rc" != 0 ] && ok "both failing fails the script (rc=$rc)" || wrong "both failing exited 0"

run "$here/apt_install.sh" 0 0
[ "$rc" = 2 ] && [ -z "$log" ] && ok "no arguments is refused, and nothing is run" \
    || wrong "no arguments: rc=$rc log='$log'"

# The control on the control. A script that forgives the install must be told apart.
sed 's|^sudo apt-get install "\$@"$|sudo apt-get install "$@" \|\| true|' "$here/apt_install.sh" > "$tmp/mutant.sh"
if ! grep -q 'apt-get install "\$@" || true' "$tmp/mutant.sh"; then
    wrong "the mutation did not apply"
else
    run "$tmp/mutant.sh" 0 100 -y no-such-package
    [ "$rc" = 0 ] && ok "the forgiving mutant exits 0 on a failed install, so that case can fail" \
        || wrong "the mutant was not forgiving: rc=$rc"
fi

[ "$bad" = 0 ] && echo "all apt gate controls gave the expected result" || echo "AT LEAST ONE CONTROL IS WRONG"
exit "$bad"
