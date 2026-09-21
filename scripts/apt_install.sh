#!/usr/bin/env bash
# Install apt packages on a CI runner. `apt-get update` is NOT the gate. The install is.
#
#   scripts/apt_install.sh -y --no-install-recommends poppler-utils
#
# WHY. GitHub's Ubuntu images list third-party sources, packages.microsoft.com among
# them, and those answer 403 from time to time. `apt-get update` then exits 100, although
# the Ubuntu archive, where every package named in these workflows lives, refreshed fine.
# On 2026-09-21 that failed python-unit on a60be28 before one test had run.
#
# WHAT THIS DOES NOT DO, AND WHY. `apt-get update || true` on its own turns a loud failure
# into a silent one. Here a failed update is printed, and the install must still succeed:
# a package that cannot be resolved fails this script by name.
# scripts/test_apt_install.sh proves both halves against a stub `sudo`.
#
# The arguments go to `apt-get install` unchanged, so each caller keeps its own flags.
# xvfb needs its recommended packages for xvfb-run, and the LaTeX step must not have them.
set -euo pipefail

[ "$#" -gt 0 ] || { echo "usage: apt_install.sh <apt-get install arguments>" >&2; exit 2; }

if ! sudo apt-get update; then
    echo "warning: apt-get update reported errors; the install below is the gate" >&2
fi
sudo apt-get install "$@"
