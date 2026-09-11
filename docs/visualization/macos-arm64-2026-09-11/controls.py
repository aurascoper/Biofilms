#!/usr/bin/env python3
"""Reproducible negative controls for the ParaView evidence tier.

viewing/README.md documents these controls firing. It documents them in PROSE --
they were run ad hoc in a shell and written up, with no script anyone can re-run.
In a repository whose own rule is "a control that cannot fire is not a control",
a control that cannot be re-run is a claim, not a receipt. This is the receipt.

Every control runs on an APFS copy-on-write clone (`cp -c -R`). The evidence is
opened read-only and never written.

Every mutation asserts it changed EXACTLY the number of sites it intended. That is
not defensive padding -- controls 1 and 2 originally reported "did not fire" and
made the verifier look weak, when in fact the regex targeted timestep="50" while
the file writes timestep="50.0", so the control silently no-opped and the verifier
correctly passed an unmodified file. A control that cannot modify its target is
indistinguishable from a check that cannot fail.

The two state controls need biofilm_signal_views.pvsm, which is 550 KB of saved
ParaView state and is deliberately NOT committed. Pass --state <pvsm> to point at it;
without one they report BLOCKED rather than crashing, because a traceback reads as
"the harness is broken" when the truth is "you did not give it the state file."

usage: controls.py <inert_signal_dir> [--receipt <render_manifest.json>]
                   [--with-state] [--state <pvsm>] [--out <dir>]
"""
import sys, os, json, shutil, hashlib, subprocess, tempfile, re

HERE = os.path.dirname(os.path.abspath(__file__))
PVPYTHON = "/Applications/ParaView-6.1.1.app/Contents/bin/pvpython"
TARGET_VTI = "paraview/signal_mcs000050.vti"
STATE = None  # set from --state; see main()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clone(src, dst):
    """APFS copy-on-write clone. Falls back to a real copy off APFS."""
    if subprocess.run(["cp", "-c", "-R", src, dst], capture_output=True).returncode != 0:
        shutil.copytree(src, dst)
    return dst


def sub_once(path, pattern, repl, expect=1):
    """Rewrite `path`, asserting the edit landed exactly `expect` times."""
    text = open(path).read()
    new, n = re.subn(pattern, repl, text)
    assert n == expect, "control mutation touched %d sites, expected %d: %s in %s" % (n, expect, pattern, path)
    open(path, "w").write(new)
    return n


def flip_one_bit(path, offset=None):
    """Flip the low bit of one byte inside the appended binary payload."""
    data = bytearray(open(path, "rb").read())
    if offset is None:
        offset = len(data) // 2
    before = data[offset]
    data[offset] ^= 0x01
    open(path, "wb").write(bytes(data))
    assert open(path, "rb").read()[offset] != before, "bit flip did not take"
    return offset


def regenerate_manifest(root):
    """What a careful tamperer does for free: rewrite the manifest and its receipt
    so the bundle is internally consistent again."""
    man = os.path.join(root, "derived_manifest.json")
    d = json.load(open(man))
    for rel in d["artifacts"]:
        p = os.path.join(root, rel)
        if os.path.isfile(p):
            d["artifacts"][rel] = sha256(p)
    json.dump(d, open(man, "w"), indent=2)
    open(os.path.join(root, "derived_manifest.sha256"), "w").write(
        "%s  derived_manifest.json\n" % sha256(man))


def run_artifacts(root, receipt, receipt_dir):
    cmd = [sys.executable, os.path.join(HERE, "verify_artifacts.py"), root]
    if receipt:
        cmd += ["--receipt", receipt]
    env = dict(os.environ, VERIFY_RECEIPT_DIR=receipt_dir)
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    nfail = None
    m = re.search(r"(\d+) checks, (\d+) failures", p.stdout)
    if m:
        nfail = int(m.group(2))
    return p.returncode, nfail, p.stdout + p.stderr


def run_state(statefile, receipt, receipt_dir):
    if not os.path.exists(PVPYTHON):
        return None, None, "pvpython not found at " + PVPYTHON
    cmd = [PVPYTHON, os.path.join(HERE, "verify_state.py"), statefile]
    if receipt:
        cmd += ["--receipt", receipt]
    p = subprocess.run(cmd, capture_output=True, text=True,
                       env=dict(os.environ, VERIFY_RECEIPT_DIR=receipt_dir))
    m = re.search(r"(\d+) checks, (\d+) failures", p.stdout)
    return p.returncode, (int(m.group(2)) if m else None), p.stdout + p.stderr


# ---------------------------------------------------------------- the controls
# Each returns (label, expectation, observed, verdict) where verdict is one of
# FIRES / DID-NOT-FIRE / BLOCKED. A control that does not fire is the finding.

def control_baseline(root, receipt, wd):
    c = clone(root, os.path.join(wd, "baseline"))
    rc, nf, log = run_artifacts(c, receipt, wd)
    return ("baseline: unmutated clone", "0 failures", "%s failures" % nf,
            "FIRES" if nf == 0 else "DID-NOT-FIRE", log)


def control_missing_timestep(root, receipt, wd):
    c = clone(root, os.path.join(wd, "missing_ts"))
    pvd = os.path.join(c, "paraview", "signal_trajectory.pvd")
    sub_once(pvd, r'\s*<DataSet timestep="50\.0"[^/]*/>\n', "\n", expect=1)
    rc, nf, log = run_artifacts(c, receipt, wd)
    return ("1: timestep 50.0 removed from .pvd", ">=1 failure", "%s failures" % nf,
            "FIRES" if nf else "DID-NOT-FIRE", log)


def control_duplicate_timestep(root, receipt, wd):
    c = clone(root, os.path.join(wd, "dup_ts"))
    pvd = os.path.join(c, "paraview", "signal_trajectory.pvd")
    text = open(pvd).read()
    m = re.search(r'[ \t]*<DataSet timestep="50\.0"[^/]*/>\n', text)
    assert m, "control could not locate the timestep 50.0 DataSet line"
    sub_once(pvd, re.escape(m.group(0)), (m.group(0) + m.group(0)).replace("\\", "\\\\"), expect=1)
    rc, nf, log = run_artifacts(c, receipt, wd)
    return ("2: timestep 50.0 duplicated", ">=1 failure", "%s failures" % nf,
            "FIRES" if nf else "DID-NOT-FIRE", log)


def control_vti_bitflip(root, receipt, wd):
    c = clone(root, os.path.join(wd, "bitflip"))
    off = flip_one_bit(os.path.join(c, TARGET_VTI))
    rc, nf, log = run_artifacts(c, receipt, wd)
    return ("3: one bit flipped in %s (byte %d)" % (os.path.basename(TARGET_VTI), off),
            ">=1 failure", "%s failures" % nf,
            "FIRES" if nf else "DID-NOT-FIRE", log)


def control_consistent_tamper(root, receipt, wd):
    """THE control this harness exists for -- correction #8.

    Alter the data, then regenerate the manifest AND its .sha256 receipt so the
    bundle self-certifies. Unpinned, every check passes: the manifest the verifier
    trusts is the one the tamperer rewrote. Pinned against a frozen receipt that
    lives beside the renderers, it must fail.
    """
    c = clone(root, os.path.join(wd, "consistent"))
    flip_one_bit(os.path.join(c, TARGET_VTI))
    regenerate_manifest(c)
    rc_u, nf_u, log_u = run_artifacts(c, None, wd)          # unpinned: expected to PASS
    rc_p, nf_p, log_p = run_artifacts(c, receipt, wd)       # pinned:   must FAIL
    verdict = "FIRES" if nf_p else "DID-NOT-FIRE"
    return ("5: consistent tamper (data + manifest + receipt all regenerated)",
            "unpinned PASSES (the hole), pinned FAILS (the fix)",
            "unpinned %s failures / pinned %s failures" % (nf_u, nf_p),
            verdict, log_u + "\n--- pinned ---\n" + log_p)


def control_state_decoy(root, receipt, wd):
    state = os.path.join(wd, "DECOY_state.pvsm")
    if not (STATE and os.path.isfile(STATE)):
        raise FileNotFoundError("no .pvsm: pass --state <biofilm_signal_views.pvsm>")
    shutil.copy(STATE, state)
    real = os.path.join(root, "paraview", "signal_trajectory.pvd")
    decoy = os.path.join(wd, "DECOY_trajectory.pvd")
    open(decoy, "w").write('<?xml version="1.0"?>\n<VTKFile type="Collection"><Collection/></VTKFile>\n')
    n = sub_once(state, re.escape(real), decoy.replace("\\", "\\\\"), expect=1)
    rc, nf, log = run_state(state, receipt, wd)
    return ("0: .pvsm repointed at a decoy .pvd", ">=1 failure", "%s failures" % nf,
            "BLOCKED" if nf is None else ("FIRES" if nf else "DID-NOT-FIRE"), log)


def control_state_lut(root, receipt, wd):
    state = os.path.join(wd, "LUT_state.pvsm")
    if not (STATE and os.path.isfile(STATE)):
        raise FileNotFoundError("no .pvsm: pass --state <biofilm_signal_views.pvsm>")
    shutil.copy(STATE, state)
    # Narrow the signal LUT 0..9 -> 0..5 without re-rendering. The rendered frames
    # would then no longer correspond to the state that claims to have produced them.
    text = open(state).read()
    m = re.search(r'(<Property name="RGBPoints".*?)(9)(\b)', text, re.S)
    assert m, "control could not locate the signal LUT upper bound in the state"
    sub_once(state, re.escape(m.group(0)), (m.group(1) + "5" + m.group(3)).replace("\\", "\\\\"), expect=1)
    rc, nf, log = run_state(state, receipt, wd)
    return ("4a: signal LUT narrowed 0..9 -> 0..5 without re-render", ">=1 failure",
            "%s failures" % nf,
            "BLOCKED" if nf is None else ("FIRES" if nf else "DID-NOT-FIRE"), log)


DATA_CONTROLS = [control_baseline, control_missing_timestep, control_duplicate_timestep,
                 control_vti_bitflip, control_consistent_tamper]
STATE_CONTROLS = [control_state_decoy, control_state_lut]


def main(argv):
    global STATE
    root = os.path.abspath(argv[1])
    STATE = os.path.abspath(argv[argv.index("--state") + 1]) if "--state" in argv \
        else os.path.join(HERE, "biofilm_signal_views.pvsm")
    # cwd, not HERE: once these scripts are committed, defaulting to HERE means a
    # verification run dirties the tree it is verifying.
    out_dir = os.path.abspath(argv[argv.index("--out") + 1]) if "--out" in argv else os.getcwd()
    receipt = argv[argv.index("--receipt") + 1] if "--receipt" in argv \
        else os.path.join(HERE, "render_manifest.json")
    receipt = os.path.abspath(receipt) if os.path.isfile(receipt) else None
    controls = DATA_CONTROLS + (STATE_CONTROLS if "--with-state" in argv else [])

    wd = tempfile.mkdtemp(prefix="biofilm_controls_", dir="/private/tmp")
    rows, logs = [], {}
    try:
        for fn in controls:
            try:
                label, expect, obs, verdict, log = fn(root, receipt, wd)
            except AssertionError as e:
                label, expect, obs, verdict, log = fn.__name__, "mutation applies", str(e), "BLOCKED", str(e)
            except FileNotFoundError as e:
                label, expect, obs, verdict, log = fn.__name__, "state file available", str(e), "BLOCKED", str(e)
            except Exception as e:
                label, expect, obs, verdict, log = fn.__name__, "runs", "%s: %s" % (type(e).__name__, e), "BLOCKED", str(e)
            rows.append({"control": label, "expectation": expect, "observed": obs, "verdict": verdict})
            logs[label] = log[-4000:]
            print("[%-12s] %-62s %s" % (verdict, label, obs))
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    bad = [r for r in rows if r["verdict"] == "DID-NOT-FIRE"]
    blocked = [r for r in rows if r["verdict"] == "BLOCKED"]
    out = {"evidence_root": root, "receipt": receipt, "state": STATE,
           "controls": rows, "all_fired": not bad and not blocked, "logs": logs}
    json.dump(out, open(os.path.join(out_dir, "control_verification.json"), "w"), indent=2)
    print("\n%d controls, %d did not fire, %d blocked" % (len(rows), len(bad), len(blocked)))
    if blocked:
        print("BLOCKED is not a pass: a control that could not run has established nothing.")
    return 1 if (bad or blocked) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
