#!/usr/bin/env python3
"""Freeze the render receipt: render_manifest.json.

The point of this file is correction #8. verify_artifacts.py used to read
derived_manifest.sha256 from beside the data, so a clone tampered *consistently*
-- bit flipped, manifest regenerated, receipt regenerated -- certified itself.
A mutable manifest sitting next to its own data is not an authority.

So the expected manifest hash is frozen HERE, beside the renderers, and
verify_artifacts.py --receipt reads it from here. Regenerating this file is a
deliberate act with its own timestamp; regenerating the one beside the data is
what an attacker (or a careless rsync) does for free.

Records nothing it did not measure. No rendering is performed and no evidence
file is written.

usage: make_render_manifest.py <inert_signal_dir> [--repo <Biofilms checkout>]
"""
import sys, os, json, hashlib, subprocess, platform, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
# Declared, not inferred: these are read off build_views.py / animate_4d.py and
# restated in one place so a figure caption can cite a single line.
RENDERING = {
    "camera": "ResetCamera, then Azimuth 35, Elevation 25, ResetCamera (build_views.py:66)",
    "slice_plane_origin": [20.0, 20.0, 28.5],
    "slice_plane_normal": [0.0, 0.0, 1.0],
    "slice_note": "z=28.5 is the CENTRE of zero-based cell layer 28. Cells span [k,k+1], "
                  "so an integer z is a face, not a layer. Chosen after inspecting the "
                  "trajectory because it carries the MCS 8 onset and the MCS 100 peak; "
                  "NOT canonical, and NOT every frame's maximum (MCS 30 peaks at z=29).",
    "species_threshold": "CELLS/species in [1,7] inclusive",
    "signal_view_threshold": "CELLS/interior_mask in [1,1]",
    "quorum_threshold_predicate": "occupied_above_threshold != 0, i.e. (occupied & signal >= 5.0). "
                                  "Thresholding `signal >= 5.0` alone disagrees with metrics.json on "
                                  "16 of 101 frames because it also catches UNOCCUPIED interior voxels.",
    "adjacency_rule": "Connectivity filter = 26-connectivity (point-sharing). Face-sharing "
                      "(6-connectivity) gives 7 at MCS 8 and 19 at the MCS 11 peak, not 6 and 18. "
                      "Every connectivity number must carry its adjacency rule.",
    "signal_lut": {"preset": "Viridis", "range": [0.0, 9.0], "AutomaticRescaleRangeMode": "Never"},
    "species_lut": "categorical, 7 annotated categories, palette verbatim from viewer/paraview_species.py",
    "raytracing": "off at render time in both views; see build_report.json.raytracing_props_zeroed. "
                  "NOTE: ParaView 6.1.1 does not serialise EnableRayTracing into a .pvsm at all, "
                  "so the reloaded-state assertion documents an invariant and is NOT a control.",
    "movie": "ffmpeg, 12 fps, fixed camera and fixed 0-9 scale across all 101 frames; see encode_movies.sh",
}
SCOPE = ("Registers artifacts and rendering settings. Confers NO evidential status on any "
         "binding or sorbate assumption, and verifies a view of a declared model -- nothing "
         "biological and nothing isotope-specific. The .pvd time key is `mcs`, a step count, "
         "not seconds; spacing is 1.0 lattice units/site and the D-PITCH refusal is intact.")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_tree(base, names):
    out = {}
    for n in names:
        p = os.path.join(base, n)
        out[n] = sha256(p) if os.path.isfile(p) else None
    return out


def protected_csvs(repo):
    """Prove the protected CSVs did not move -- against HEAD, not against a prior run.

    A before/after pair only proves this pipeline did not touch them. Comparing the
    working file's git blob id to HEAD's proves nothing touched them since the commit,
    which is the claim actually worth making.
    """
    out = {}
    for rel in ("data/claims_ledger.csv", "data/parameter_provenance.csv"):
        p = os.path.join(repo, rel)
        if not os.path.isfile(p):
            out[rel] = {"status": "ABSENT", "path": p}
            continue
        g = lambda *a: subprocess.run(a, cwd=repo, capture_output=True, text=True).stdout.strip()
        work, head = g("git", "hash-object", rel), g("git", "rev-parse", "HEAD:" + rel)
        out[rel] = {"sha256": sha256(p), "git_blob_worktree": work, "git_blob_HEAD": head,
                    "unmodified_vs_HEAD": bool(head) and work == head}
    return out


def main(argv):
    root = os.path.abspath(argv[1])
    repo = os.path.abspath(argv[argv.index("--repo") + 1]) if "--repo" in argv \
        else os.path.expanduser("~/Developer/Biofilms-lattice-viewer-mac")

    man = os.path.join(root, "derived_manifest.json")
    d = json.load(open(man))
    pvd = os.path.join(root, "paraview", "signal_trajectory.pvd")
    rcpt = os.path.join(root, "derived_manifest.sha256")

    out = {
        "schema_version": 1,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "purpose": "Frozen expected-hash receipt for the ParaView macOS/arm64 render. "
                   "verify_artifacts.py --receipt pins against THIS file, not against the "
                   "mutable manifest beside the data.",
        "scope": SCOPE,
        "evidence_root": root,
        # --- the freeze ---
        "pinned": {
            "derived_manifest_sha256": sha256(man),
            "derived_manifest_receipt_text": open(rcpt).read().strip() if os.path.isfile(rcpt) else None,
            "parent_manifest_sha256": d["parent_manifest_sha256"],
            "pvd_sha256": sha256(pvd),
            "artifact_count": len(d["artifacts"]),
            "expected_timesteps": 101,
        },
        "parent": {"run_id": d["run_id"], "parent_run_id": d["parent_run_id"],
                   "source_commit": d["source_commit"], "julia_version": d["julia_version"],
                   "julia_executable_sha256": d["julia_executable_sha256"],
                   "hash_contract": d["hash_contract"], "created_utc": d["created_utc"],
                   "configuration": d["configuration"], "scope": d["scope"]},
        "environment": {
            "paraview": "6.1.1", "paraview_app": "/Applications/ParaView-6.1.1.app",
            "paraview_arch": "native arm64 (brew --cask)",
            "python": platform.python_version(),
            "platform": platform.platform(), "machine": platform.machine(),
            "note": "libopenvkl_module_cpu_device.dylib is absent from the macOS bundle; "
                    "every pvpython run prints [openvkl] INITIALIZATION ERROR. That is the "
                    "OSPRay volume path, unused here.",
        },
        "renderer_sources": hash_tree(HERE, [
            "build_views.py", "animate_4d.py", "verify_state.py", "verify_artifacts.py",
            "controls.py", "encode_movies.sh", "make_render_manifest.py"]),
        "exporter_sources": hash_tree(repo, ["export_vti.jl", "viewer/paraview_species.py"]),
        "rendering": RENDERING,
        "outputs": hash_tree(HERE, sorted(
            f for f in os.listdir(HERE)
            if f.endswith((".png", ".mp4", ".pvsm")))),
        "protected_csv": protected_csvs(repo),
    }
    dest = os.path.join(HERE, "render_manifest.json")
    json.dump(out, open(dest, "w"), indent=2)
    print("wrote %s" % dest)
    print("  pinned derived_manifest_sha256 = %s" % out["pinned"]["derived_manifest_sha256"])
    print("  pinned pvd_sha256              = %s" % out["pinned"]["pvd_sha256"])
    print("  artifacts pinned               = %d" % out["pinned"]["artifact_count"])
    for k, v in out["protected_csv"].items():
        print("  %-34s unmodified_vs_HEAD=%s" % (k, v.get("unmodified_vs_HEAD")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
