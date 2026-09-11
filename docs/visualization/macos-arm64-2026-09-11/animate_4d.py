#!/usr/bin/env python3
"""Render all 101 frames as 4D voxels: parcels, and the supra-threshold signal region.

Two independent 3D voxel views, one PNG each per MCS, fixed camera and fixed colour
scales so frames are comparable. Ray tracing off. Writes nothing outside <outdir>.

usage: pvpython animate_4d.py <signal_trajectory.pvd> <outdir>
"""
from paraview.simple import *
import sys, os, json

pvd, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
SIGNAL_RANGE = (0.0, 9.0)
THRESHOLD = 5.0                      # the declared threshold from the run's config.toml
RES = [1100, 900]

COLORS = {1: (0.902, 0.098, 0.294), 2: (0.235, 0.706, 0.294), 3: (0.263, 0.388, 0.847),
          4: (0.961, 0.510, 0.192), 5: (0.569, 0.118, 0.706), 6: (0.259, 0.831, 0.957),
          7: (0.941, 0.196, 0.902)}
LABELS = {1: "C. neoformans", 2: "D. radiodurans", 3: "C. sphaerospermum", 4: "B. subtilis",
          5: "A. niger", 6: "S. oneidensis", 7: "O. intermedium"}

reader = OpenDataFile(pvd)
reader.UpdatePipelineInformation()
ts = [int(t) for t in reader.TimestepValues]
print("frames:", len(ts))

di = reader.GetCellDataInformation()
arrays = [di.GetArray(i).GetName() for i in range(di.GetNumberOfArrays())]
print("cell arrays:", arrays)
for need in ("species", "occupied_above_threshold", "signal"):
    if need not in arrays:
        raise SystemExit("FATAL: required cell array %r absent from %s. ParaView's Threshold "
                         "returns empty output for a missing array rather than raising, so this "
                         "would render blank frames and write a clean all-zero receipt." % (need, pvd))

# THRESHOLD is used only for the burned-in frame label and the receipt; the voxel counts come
# from the data's own occupied_above_threshold array. A mismatch would mislabel every frame
# while the receipt stayed self-consistent, so reconcile against the run's config.
cfg = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(pvd))), "provenance", "config.toml")
if os.path.exists(cfg):
    import re as _re
    m = _re.search(r"^threshold\s*=\s*([0-9.]+)", open(cfg).read(), _re.M)
    if m and abs(float(m.group(1)) - THRESHOLD) > 1e-12:
        raise SystemExit("FATAL: THRESHOLD=%s but %s declares %s" % (THRESHOLD, cfg, m.group(1)))
    print("threshold reconciled against", cfg)
else:
    print("WARNING: no config.toml beside the data; THRESHOLD=%s is unreconciled" % THRESHOLD)

LoadPalette(paletteName='WhiteBackground')

# --- View A: parcel voxels ---
va = CreateRenderView(); va.ViewSize = RES; va.OrientationAxesVisibility = 1
va.EnableRayTracing = 0
thr_sp = Threshold(Input=reader, registrationName="Parcels")
thr_sp.Scalars = ["CELLS", "species"]; thr_sp.LowerThreshold = 1; thr_sp.UpperThreshold = 7
da = Show(thr_sp, va); da.SetRepresentationType("Surface With Edges"); da.EdgeColor = [0.25]*3
ColorBy(da, ("CELLS", "species"))
lut_sp = GetColorTransferFunction("species"); lut_sp.InterpretValuesAsCategories = 1
lut_sp.Annotations = [s for k in sorted(COLORS) for s in (str(k), LABELS[k])]
lut_sp.IndexedColors = [c for k in sorted(COLORS) for c in COLORS[k]]
da.SetScalarBarVisibility(va, True)
GetScalarBar(lut_sp, va).Title = "species"; GetScalarBar(lut_sp, va).ComponentTitle = ""

# --- View B: supra-threshold signal voxels ---
vb = CreateRenderView(); vb.ViewSize = RES; vb.OrientationAxesVisibility = 1
vb.EnableRayTracing = 0
# Threshold on the run's OWN endpoint array, not on `signal` directly.
# InertSignal.jl:82,86 define the endpoint as `occupied & (A >= threshold)` -- occupied
# interior voxels only. Thresholding `signal >= 5.0` instead also catches unoccupied voxels
# the field has diffused into, and disagrees with metrics.json on 16 of 101 frames.
# Summing `occupied_above_threshold` reproduces metrics.json exactly on 101/101 frames.
thr_sg = Threshold(Input=reader, registrationName="OccupiedAboveThreshold")
thr_sg.Scalars = ["CELLS", "occupied_above_threshold"]
thr_sg.LowerThreshold = 1; thr_sg.UpperThreshold = 1
db = Show(thr_sg, vb); db.SetRepresentationType("Surface With Edges"); db.EdgeColor = [0.25]*3
ColorBy(db, ("CELLS", "signal"))
lut_sg = GetColorTransferFunction("signal"); lut_sg.ApplyPreset("Viridis", True)
lut_sg.RescaleTransferFunction(*SIGNAL_RANGE); lut_sg.AutomaticRescaleRangeMode = "Never"
GetOpacityTransferFunction("signal").RescaleTransferFunction(*SIGNAL_RANGE)
db.SetScalarBarVisibility(vb, True)
GetScalarBar(lut_sg, vb).Title = "signal in quorum voxels"
GetScalarBar(lut_sg, vb).ComponentTitle = "(lattice units)"

# Fix BOTH cameras on the full lattice bounds so no frame rescales.
scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = float(ts[-1])
for v, src in ((va, thr_sp), (vb, thr_sg)):
    v.ResetCamera()
Render(va); Render(vb)
# widen to the whole 40^3 domain, identical for both, so growth reads as growth
for v in (va, vb):
    v.CameraFocalPoint = [20, 20, 20]
    v.CameraPosition = [20 + 95, 20 - 62, 20 + 68]
    v.CameraViewUp = [0, 0, 1]
    v.CameraParallelProjection = 0

ta = Text(registrationName="LabelA", Text=""); Show(ta, va).FontSize = 24
tb = Text(registrationName="LabelB", Text=""); Show(tb, vb).FontSize = 24

counts = {}
for t in ts:
    scene.AnimationTime = float(t)
    thr_sg.UpdatePipeline(float(t))
    n = thr_sg.GetDataInformation().GetNumberOfCells()
    counts[t] = n
    ta.Text = "MCS %3d   parcels" % t
    tb.Text = "MCS %3d   occupied & signal >= %.1f : %d voxels" % (t, THRESHOLD, n)
    Render(va); Render(vb)
    SaveScreenshot("%s/anim_parcels_%03d.png" % (outdir, t), va, ImageResolution=RES)
    SaveScreenshot("%s/anim_quorum_%03d.png" % (outdir, t), vb, ImageResolution=RES)
    if t % 20 == 0 or t in (7, 8):
        print("  MCS %3d  supra-threshold voxels: %d" % (t, n))

json.dump({"threshold": THRESHOLD, "signal_range": list(SIGNAL_RANGE),
           "frames": len(ts), "supra_threshold_voxels_by_mcs": counts,
           "first_nonzero_mcs": next((t for t in ts if counts[t] > 0), None)},
          open("%s/anim_report.json" % outdir, "w"), indent=2)
print("done")
