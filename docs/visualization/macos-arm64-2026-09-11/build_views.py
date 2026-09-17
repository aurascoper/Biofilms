# Build the species and signal views from signal_trajectory.pvd, screenshot the named MCS
# frames, and save one .pvsm. Palette is taken verbatim from Biofilms viewer/paraview_species.py.
# usage: pvpython build_views.py <signal_trajectory.pvd> <outdir>
from paraview.simple import *
import sys, os, json

pvd, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
FRAMES = [0, 7, 8, 30, 100]          # 7 = before, 8 = first threshold crossing
SIGNAL_RANGE = (0.0, 9.0)            # global max across 101 frames is 8.9034
# Slice plane chosen from the data, not from the geometric centre. The z=20 mid-plane
# peaks at 4.45 at MCS 8 and never crosses the 5.0 threshold, so it cannot show the onset.
# The 3D maximum sits near z=28 at MCS 8 (5.2971), 30 (8.4933) and 100 (8.9034).
SLICE_Z = 28.5

COLORS = {1: (0.902, 0.098, 0.294), 2: (0.235, 0.706, 0.294), 3: (0.263, 0.388, 0.847),
          4: (0.961, 0.510, 0.192), 5: (0.569, 0.118, 0.706), 6: (0.259, 0.831, 0.957),
          7: (0.941, 0.196, 0.902)}
LABELS = {1: "C. neoformans", 2: "D. radiodurans", 3: "C. sphaerospermum", 4: "B. subtilis",
          5: "A. niger", 6: "S. oneidensis", 7: "O. intermedium"}

def no_raytracing(v):
    # ParaView 5.7+ exposes only EnableRayTracing; probing the old EnableOSPRay name
    # raises NotSupportedException rather than returning False, so do not hasattr it.
    try:
        v.EnableRayTracing = 0
        return "EnableRayTracing=0"
    except Exception as e:
        return "raytracing property unavailable: %s" % e

reader = OpenDataFile(pvd)
reader.UpdatePipelineInformation()
ts = list(reader.TimestepValues)
print("TIMESTEPS:", len(ts), "first", ts[0], "last", ts[-1])
di = reader.GetCellDataInformation()
arrays = [di.GetArray(i).GetName() for i in range(di.GetNumberOfArrays())]
print("CELL ARRAYS:", arrays)

# ---- View 1: parcel geometry ----
# pvpython batch mode creates no implicit layout; make one so the .pvsm opens side-by-side.
v1 = CreateRenderView()
try:
    layout = CreateLayout(name="Biofilm signal views")
    AssignViewToLayout(view=v1, layout=layout, hint=0)
    layout.SplitHorizontal(0, 0.5)
    LAYOUT_OK = True
except Exception as e:
    print("LAYOUT FALLBACK (views still saved, GUI will arrange):", e)
    layout = None; LAYOUT_OK = False
LoadPalette(paletteName='WhiteBackground')
v1.ViewSize = [1200, 1000]; v1.OrientationAxesVisibility = 1
rt1 = no_raytracing(v1)

thr_sp = Threshold(Input=reader, registrationName="Species_1to7")
thr_sp.Scalars = ["CELLS", "species"]
thr_sp.LowerThreshold = 1; thr_sp.UpperThreshold = 7
d1 = Show(thr_sp, v1)
d1.SetRepresentationType("Surface With Edges"); d1.EdgeColor = [0.2, 0.2, 0.2]
ColorBy(d1, ("CELLS", "species"))
lut_sp = GetColorTransferFunction("species")
lut_sp.InterpretValuesAsCategories = 1
lut_sp.Annotations = [s for k in sorted(COLORS) for s in (str(k), LABELS[k])]
lut_sp.IndexedColors = [c for k in sorted(COLORS) for c in COLORS[k]]
d1.SetScalarBarVisibility(v1, True)
GetScalarBar(lut_sp, v1).Title = "species"; GetScalarBar(lut_sp, v1).ComponentTitle = ""
v1.ResetCamera(); GetActiveCamera().Azimuth(35); GetActiveCamera().Elevation(25); v1.ResetCamera()
t1 = Text(registrationName="Label_species", Text="MCS 0")
td1 = Show(t1, v1); td1.FontSize = 26; td1.Color = [0, 0, 0]

# ---- View 2: signal distribution ----
v2 = CreateRenderView()
if LAYOUT_OK:
    AssignViewToLayout(view=v2, layout=layout, hint=2)
v2.ViewSize = [1200, 1000]; v2.OrientationAxesVisibility = 1
rt2 = no_raytracing(v2)

thr_in = Threshold(Input=reader, registrationName="Interior_mask_1")
thr_in.Scalars = ["CELLS", "interior_mask"]
thr_in.LowerThreshold = 1; thr_in.UpperThreshold = 1
sl = Slice(Input=thr_in, registrationName="Signal_slice")
sl.SliceType = "Plane"; sl.SliceType.Origin = [20.0, 20.0, SLICE_Z]; sl.SliceType.Normal = [0.0, 0.0, 1.0]
d2 = Show(sl, v2)
d2.SetRepresentationType("Surface")
ColorBy(d2, ("CELLS", "signal"))
lut_sg = GetColorTransferFunction("signal")
lut_sg.ApplyPreset("Viridis", True)  # 6.x dropped the "(matplotlib)" suffix
lut_sg.RescaleTransferFunction(*SIGNAL_RANGE)
lut_sg.AutomaticRescaleRangeMode = "Never"
GetOpacityTransferFunction("signal").RescaleTransferFunction(*SIGNAL_RANGE)
d2.SetScalarBarVisibility(v2, True)
GetScalarBar(lut_sg, v2).Title = "signal (lattice units)"; GetScalarBar(lut_sg, v2).ComponentTitle = ""
v2.ResetCamera()
t2 = Text(registrationName="Label_signal", Text="MCS 0")
td2 = Show(t2, v2); td2.FontSize = 26; td2.Color = [0, 0, 0]

scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps()
report = {"pvd": os.path.abspath(pvd), "timesteps": len(ts),
          "raytracing_props_zeroed": {"view1": rt1, "view2": rt2},
          "signal_range": list(SIGNAL_RANGE), "cell_arrays": arrays,
          "layout_split": LAYOUT_OK, "slice_z": SLICE_Z, "paraview": GetParaViewVersion().GetVersion() if hasattr(GetParaViewVersion(),"GetVersion") else str(GetParaViewVersion()), "frames": {}}

for t in FRAMES:
    scene.AnimationTime = float(t)
    t1.Text = "MCS %d  -- parcels" % t
    t2.Text = "MCS %d  -- signal" % t
    Render(v1); Render(v2)
    p1 = "%s/species_mcs%03d.png" % (outdir, t)
    p2 = "%s/signal_z%s_mcs%03d.png" % (outdir, ("%g" % SLICE_Z).replace(".", "p"), t)
    SaveScreenshot(p1, v1, ImageResolution=[1200, 1000])
    SaveScreenshot(p2, v2, ImageResolution=[1200, 1000])
    sl.UpdatePipeline(float(t))
    _a = sl.GetCellDataInformation().GetArray("signal")
    rng = _a.GetRange() if _a is not None else (float("nan"), float("nan"))
    report["frames"][t] = {"species_png": os.path.basename(p1),
                           "signal_png": os.path.basename(p2),
                           "slice_signal_range": [round(rng[0], 6), round(rng[1], 6)]}
    print("MCS %3d  slice signal range %.4f .. %.4f" % (t, rng[0], rng[1]))

state = "%s/biofilm_signal_views.pvsm" % outdir
SaveState(state)
report["pvsm"] = os.path.abspath(state)
json.dump(report, open("%s/build_report.json" % outdir, "w"), indent=2)
print("WROTE STATE:", state)
