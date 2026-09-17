#!/usr/bin/env python3
"""Artifact-integrity tier: the .pvd's timestep set and every file it references.

Separate from verify_state.py, which checks the saved ParaView state. This checks the
DATA the state points at. Exits non-zero on any failure.

With --receipt, the expected manifest hash comes from a FROZEN render_manifest.json that
lives beside the renderers, not from the derived_manifest.sha256 sitting beside the data.
That is correction #8: a clone tampered consistently -- data altered, manifest regenerated,
receipt regenerated -- passes every unpinned check, because the manifest it is checked
against is the one the tamperer rewrote. A mutable manifest must not be its own authority.

usage: verify_artifacts.py <inert_signal_dir> [--receipt <render_manifest.json>]
"""
import sys, os, json, hashlib, xml.etree.ElementTree as ET

EXPECTED_TIMESTEPS = 101


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt_policy(receipt, receipt_arg_given):
    """Decide what running WITHOUT a frozen receipt means.

    `receipt` is the loaded render_manifest.json dict, or None when --receipt was
    not supplied (or the file was unreadable). `receipt_arg_given` distinguishes
    "user asked for a receipt and it failed to load" from "user never asked".

    Return (proceed, severity, note):
      proceed  -- False aborts the run before any check
      severity -- "ok" | "warn" | "fail"; "fail" records a failure in the receipt
                  JSON so the exit code is non-zero even if every hash matches
      note     -- one line, printed and stored

    The trade-off: fail-closed is honest -- an unpinned run cannot distinguish real
    evidence from a consistently tampered clone, which is the exact hole this flag
    exists to close -- but it breaks every existing invocation, including the
    baseline in viewing/README.md and anything that shells this script. Degrading
    with a loud warning keeps those working at the cost of a mode that can quietly
    self-certify if nobody reads stderr.

    Chosen policy: fail-closed only on a BAD --receipt. A missing pin is almost
    always a typo or a moved file, and a typo that silently downgrades the trust
    model is exactly the bug class this change exists to close. A deliberate
    unpinned run stays possible, loudly, for the README baseline.
    """
    if receipt:
        return True, "ok", "pinned against frozen receipt"
    if receipt_arg_given:
        return False, "fail", ("--receipt given but unreadable; refusing to fall back "
                               "to the manifest beside the data")
    return True, "warn", ("UNPINNED: a consistently tampered clone passes every check "
                          "in this mode")


def main(argv):
    root = argv[1]
    receipt_arg_given = "--receipt" in argv
    receipt = None
    if receipt_arg_given:
        rp = argv[argv.index("--receipt") + 1]
        if os.path.isfile(rp):
            receipt = json.load(open(rp))

    proceed, severity, note = receipt_policy(receipt, receipt_arg_given)
    print("receipt policy: [%s] %s" % (severity, note))
    if not proceed:
        return 2

    fails, checks = [], []

    def ck(name, cond, got):
        checks.append({"check": name, "pass": bool(cond), "observed": str(got)})
        if not cond:
            fails.append(name)
        print(("  PASS  " if cond else "  FAIL  ") + name + "  ->  " + str(got))

    if severity == "fail":
        ck("frozen receipt supplied", False, note)

    pvd = os.path.join(root, "paraview", "signal_trajectory.pvd")
    ck("pvd exists", os.path.exists(pvd), pvd)
    if not os.path.exists(pvd):
        return 1

    coll = ET.parse(pvd).getroot().find("Collection")
    ck("pvd has a <Collection> element", coll is not None, coll is not None)
    if coll is None:
        return 1
    datasets = coll.findall("DataSet")
    times = [float(d.get("timestep")) for d in datasets]
    files = [d.get("file") for d in datasets]

    ck(f"{EXPECTED_TIMESTEPS} DataSet entries", len(datasets) == EXPECTED_TIMESTEPS, len(datasets))
    ck("timesteps unique", len(set(times)) == len(times),
       f"{len(set(times))} unique of {len(times)}")
    ck("timesteps ascending", times == sorted(times), "ascending" if times == sorted(times) else "OUT OF ORDER")
    ck("timesteps are 0..100 complete", set(times) == set(float(i) for i in range(EXPECTED_TIMESTEPS)),
       f"{min(times)}..{max(times)}" if times else "EMPTY")

    missing = [f for f in files if not os.path.exists(os.path.join(root, "paraview", f))]
    ck("every referenced file exists", not missing, f"{len(missing)} missing" + (f": {missing[:3]}" if missing else ""))

    man = os.path.join(root, "derived_manifest.json")
    ck("derived manifest exists", os.path.exists(man), man)
    if os.path.exists(man):
        man_hash = sha256(man)
        d = json.load(open(man))
        rpath = os.path.join(root, "derived_manifest.sha256")
        rtext = open(rpath).read().split() if os.path.exists(rpath) else []
        ck("receipt file exists and is non-empty", bool(rtext),
           "missing" if not os.path.exists(rpath) else ("empty" if not rtext else "ok"))
        if rtext:
            # Self-consistency only. Passes on a consistently tampered clone by design;
            # the pinned checks below are what actually bite.
            ck("manifest matches its own local receipt (self-consistency)",
               rtext[0] == man_hash, rtext[0][:16] + "...")

        # --- pinned tier: the mutable manifest is checked against the frozen receipt ---
        if receipt:
            pin = receipt["pinned"]
            ck("derived manifest matches FROZEN receipt hash", pin["derived_manifest_sha256"] == man_hash,
               "pinned %s / got %s" % (pin["derived_manifest_sha256"][:16], man_hash[:16]))
            ck("parent manifest hash matches FROZEN receipt",
               pin["parent_manifest_sha256"] == d.get("parent_manifest_sha256"),
               "pinned %s / got %s" % (pin["parent_manifest_sha256"][:16],
                                       str(d.get("parent_manifest_sha256"))[:16]))
            pvd_hash = sha256(pvd)
            ck("pvd matches FROZEN receipt hash", pin["pvd_sha256"] == pvd_hash,
               "pinned %s / got %s" % (pin["pvd_sha256"][:16], pvd_hash[:16]))
            ck("artifact count matches FROZEN receipt",
               pin["artifact_count"] == len(d["artifacts"]),
               "pinned %d / got %d" % (pin["artifact_count"], len(d["artifacts"])))

        arts = d["artifacts"]
        bad, absent = [], []
        for rel, want in arts.items():
            p = os.path.join(root, rel)
            if not os.path.exists(p):
                absent.append(rel)
            elif sha256(p) != want:
                bad.append(rel)
        ck(f"all {len(arts)} artifact hashes match", not bad and not absent,
           f"{len(arts) - len(bad) - len(absent)}/{len(arts)} ok, {len(bad)} altered, {len(absent)} missing"
           + (f"; first altered: {bad[0]}" if bad else ""))

        # every .vti the pvd references must be a hashed artifact
        unregistered = [f for f in files if os.path.join("paraview", f) not in arts]
        ck("every pvd-referenced vti is registered in the manifest", not unregistered,
           f"{len(unregistered)} unregistered" + (f": {unregistered[:3]}" if unregistered else ""))

    print("\n%d checks, %d failures" % (len(checks), len(fails)))
    # Write to the CURRENT DIRECTORY. Two earlier defaults were both wrong in the same
    # way -- each wrote into a tree that is read-only by intent. Beside the data wrote a
    # receipt into the evidence bundle on every run; beside this script wrote into the
    # committed tree once these files were checked in, so verifying the repo dirtied it.
    # cwd surprises nobody: run from this directory and you get the old behaviour, run
    # from anywhere else and the receipt lands where you are. The APFS-clone clobber that
    # motivated the data-adjacent default is handled where it arises -- controls.py sets
    # VERIFY_RECEIPT_DIR to its own scratch dir.
    out = os.environ.get("VERIFY_RECEIPT_DIR") or os.getcwd()
    json.dump({"root": os.path.abspath(root), "pinned_by": (rp if receipt else None),
               "receipt_policy": {"severity": severity, "note": note},
               "checks": checks, "failures": fails},
              open(os.path.join(out, "artifact_verification.json"), "w"), indent=2)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
