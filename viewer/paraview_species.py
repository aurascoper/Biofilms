# Render lattice.pvd: species voxels, categorical colours, air removed; save a screenshot per MCS
# and a .pvsm state the GUI can load.
from paraview.simple import *
import sys
pvd, outdir = sys.argv[1], sys.argv[2]   # usage: pvpython viewer/paraview_species.py lattice.pvd <outdir>
COLORS = {1: (0.902, 0.098, 0.294), 2: (0.235, 0.706, 0.294), 3: (0.263, 0.388, 0.847), 4: (0.961, 0.510, 0.192),
          5: (0.569, 0.118, 0.706), 6: (0.259, 0.831, 0.957), 7: (0.941, 0.196, 0.902)}
LABELS = {1: "C. neoformans", 2: "D. radiodurans", 3: "C. sphaerospermum", 4: "B. subtilis",
          5: "A. niger", 6: "S. oneidensis", 7: "O. intermedium"}
reader = OpenDataFile(pvd)
reader.UpdatePipelineInformation()
thr = Threshold(Input=reader)
thr.Scalars = ["CELLS", "species"]
thr.LowerThreshold = 1; thr.UpperThreshold = 7
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1400, 1000]; view.Background = [1, 1, 1]; view.OrientationAxesVisibility = 1
view.UseColorPaletteForBackground = 0   # else ParaView 6 paints its palette grey over Background
view.OrientationAxesLabelColor = [0, 0, 0]
disp = Show(thr, view)
disp.SetRepresentationType("Surface With Edges"); disp.EdgeColor = [0.2, 0.2, 0.2]
ColorBy(disp, ("CELLS", "species"))
lut = GetColorTransferFunction("species")
lut.InterpretValuesAsCategories = 1
lut.Annotations = [s for k in sorted(COLORS) for s in (str(k), LABELS[k])]
lut.IndexedColors = [c for k in sorted(COLORS) for c in COLORS[k]]
disp.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view); bar.Title = "species"; bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]   # palette text is white; invisible on the white background
view.ResetCamera(); cam = GetActiveCamera(); cam.Azimuth(35); cam.Elevation(25); view.ResetCamera()
text = Text(Text="MCS 30"); tdisp = Show(text, view); tdisp.FontSize = 28; tdisp.Color = [0, 0, 0]
scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps()
for t in reader.TimestepValues:
    scene.AnimationTime = t
    text.Text = "MCS %d" % t
    Render(view)
    SaveScreenshot("%s/species_mcs%03d.png" % (outdir, int(t)), view, ImageResolution=[1400, 1000])
SaveState("%s/species_view.pvsm" % outdir)
print("wrote", len(reader.TimestepValues), "frames and species_view.pvsm to", outdir)
