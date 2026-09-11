# Reopen the saved state in a FRESH ParaView process and confirm data paths, timesteps,
# and the settings the recipe specifies. Exits non-zero on any failure.
from paraview.simple import *
import sys, os, json, hashlib
state = sys.argv[1]
# Optional frozen receipt. verify_state.py had the same hole verify_artifacts.py did:
# it bound the .pvd to "the manifest that governs the bundle it sits in", which a
# consistently tampered clone rewrites for free. --receipt pins to render_manifest.json,
# which lives beside the renderers, not beside the data.
# Same contract as verify_artifacts.py:receipt_policy -- fail-closed on a BAD --receipt.
# This used to be `if os.path.isfile(_rp): RECEIPT = ...` with no else, so a typo'd path
# silently ran UNPINNED and, because nothing was printed, produced output byte-identical
# to a deliberate unpinned run. A typo that downgrades the trust model is the exact bug
# class the pinning exists to close.
_REQUIRED = ("derived_manifest_sha256", "parent_manifest_sha256", "pvd_sha256", "artifact_count")
RECEIPT = None
if "--receipt" in sys.argv:
    _rp = sys.argv[sys.argv.index("--receipt") + 1]
    if not os.path.isfile(_rp):
        print("receipt policy: [fail] --receipt given but not readable: %s" % _rp)
        print("  refusing to fall back to the manifest beside the data")
        sys.exit(2)
    try:
        _pin = json.load(open(_rp))["pinned"]
    except Exception as _e:
        print("receipt policy: [fail] --receipt could not be parsed: %s: %s"
              % (type(_e).__name__, _e))
        sys.exit(2)
    _missing = [k for k in _REQUIRED if k not in _pin]
    if _missing:
        print("receipt policy: [fail] receipt is missing required keys: %s" % _missing)
        sys.exit(2)
    RECEIPT = _pin
    print("receipt policy: [ok] pinned against frozen receipt (trusted-identity mode)")
else:
    print("receipt policy: [warn] UNPINNED (self-consistency mode): the .pvd is checked "
          "against the manifest beside the data, which a consistently tampered clone rewrites")
LoadState(state)
fails, checks = [], []
def ck(name, cond, got):
    checks.append({"check": name, "pass": bool(cond), "observed": str(got)})
    if not cond: fails.append(name)
    print(("  PASS  " if cond else "  FAIL  ") + name + "  ->  " + str(got))

srcs = GetSources()
by_kind = {k[0]: v for k, v in srcs.items()}
reader = next(v for k, v in srcs.items() if v.GetXMLName() in ("PVDReader",))
fn = reader.FileName
ck("reader points at a signal_trajectory.pvd", fn.endswith("inert_signal/paraview/signal_trajectory.pvd"), fn)
ck("referenced .pvd exists on disk", os.path.exists(fn), os.path.exists(fn))
# A path suffix identifies a filename, not a file. Bind to content: the .pvd is itself a
# hashed artifact, so check it against the manifest that governs the bundle it sits in.
if os.path.exists(fn):
    root = os.path.dirname(os.path.dirname(fn))
    manp = os.path.join(root, "derived_manifest.json")
    if os.path.exists(manp):
        want = json.load(open(manp))["artifacts"].get("paraview/signal_trajectory.pvd")
        got = hashlib.sha256(open(fn, "rb").read()).hexdigest()
        ck("the .pvd content matches its manifest hash (self-consistency)", want == got,
           (got[:16] + "...") if want else "not registered in manifest")
        if RECEIPT:
            ck("the .pvd matches the FROZEN receipt hash", RECEIPT["pvd_sha256"] == got,
               "pinned %s / got %s" % (RECEIPT["pvd_sha256"][:16], got[:16]))
            mh = hashlib.sha256(open(manp, "rb").read()).hexdigest()
            ck("derived manifest matches the FROZEN receipt hash",
               RECEIPT["derived_manifest_sha256"] == mh,
               "pinned %s / got %s" % (RECEIPT["derived_manifest_sha256"][:16], mh[:16]))
    else:
        ck("a derived_manifest.json governs the referenced .pvd", False, "absent: " + manp)
reader.UpdatePipelineInformation()
ts = list(reader.TimestepValues)
ck("101 timesteps", len(ts) == 101, len(ts))
ck("timesteps are exactly {0..100}", set(ts) == set(float(i) for i in range(101)),
   f"{min(ts) if ts else 'EMPTY'}..{max(ts) if ts else ''}, {len(set(ts))} unique")

# Presence of every proxy is established BEFORE any property on it is read. The checks
# used to read thr_sp.LowerThreshold immediately after the presence check and only abort
# further down, so a missing filter raised AttributeError instead of reporting a failure.
thr_sp = by_kind.get("Species_1to7")
thr_in = by_kind.get("Interior_mask_1")
sl = by_kind.get("Signal_slice")
ck("Species_1to7 filter present", thr_sp is not None, thr_sp is not None)
ck("Interior_mask_1 filter present", thr_in is not None, thr_in is not None)
ck("Signal_slice filter present", sl is not None, sl is not None)
if None in (thr_sp, thr_in, sl):
    json.dump({"state": state, "mode": "pinned" if RECEIPT else "self-consistency",
               "checks": checks, "failures": fails},
              open(os.path.join(os.path.dirname(state) or ".", "state_verification.json"), "w"), indent=2)
    print("\n%d checks, %d failures (aborted: a required filter is missing)" % (len(checks), len(fails)))
    sys.exit(1)

ck("species threshold 1..7", thr_sp.LowerThreshold == 1 and thr_sp.UpperThreshold == 7,
   f"{thr_sp.LowerThreshold}..{thr_sp.UpperThreshold}")
ck("species threshold on CELLS/species", list(thr_sp.Scalars) == ["CELLS", "species"], list(thr_sp.Scalars))

ck("interior_mask threshold 1..1", thr_in.LowerThreshold == 1 and thr_in.UpperThreshold == 1,
   f"{thr_in.LowerThreshold}..{thr_in.UpperThreshold}")
ck("slice plane at z=28.5", abs(sl.SliceType.Origin[2] - 28.5) < 1e-9, list(sl.SliceType.Origin))

lut_sp = GetColorTransferFunction("species")
ck("species LUT categorical", lut_sp.InterpretValuesAsCategories == 1, lut_sp.InterpretValuesAsCategories)
ck("species LUT has 7 annotated categories", len(lut_sp.Annotations) == 14, len(lut_sp.Annotations) // 2)

lut_sg = GetColorTransferFunction("signal")
pts = list(lut_sg.RGBPoints)
lo, hi = pts[0], pts[-4]
ck("signal LUT range 0..9", abs(lo - 0.0) < 1e-9 and abs(hi - 9.0) < 1e-9, f"{lo}..{hi}")
ck("signal LUT rescale locked (Never)", lut_sg.AutomaticRescaleRangeMode == "Never", lut_sg.AutomaticRescaleRangeMode)

# NOT LOAD-BEARING, and kept only as documentation of the invariant.
# ParaView 6.1.1 does not serialise EnableRayTracing into a .pvsm at all -- verified by
# saving one state with it at the default 0 and another with it explicitly 1: neither file
# mentions the property. So a reloaded state always reports 0 regardless of how the frames
# were actually rendered, and this assertion would pass on a state saved with ray tracing on.
# The render-time value is evidence; it lives in build_report.json["raytracing_props_zeroed"].
for i, v in enumerate(GetRenderViews()):
    ck(f"view {i+1} reports ray tracing off on reload (invariant, not a control)",
       v.EnableRayTracing == 0, v.EnableRayTracing)

json.dump({"state": state, "mode": "pinned" if RECEIPT else "self-consistency",
           "checks": checks, "failures": fails},
          open(os.path.join(os.path.dirname(state) or ".", "state_verification.json"), "w"), indent=2)
print("\n%d checks, %d failures" % (len(checks), len(fails)))
sys.exit(1 if fails else 0)
