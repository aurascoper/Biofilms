# Render a scalar CPM field from lattice.pvd: a longitudinal slice through the cylinder
# axis beside a volume view, one PNG per MCS, a .pvsm state and a JSON render report.
# Mechanics follow docs/visualization/macos-arm64-2026-09-11/build_views.py (ParaView
# 6.1.1 arm64): CreateLayout, WhiteBackground, Viridis, ray tracing off, fixed LUT.
#
# usage: pvpython viewer/paraview_fields.py <lattice.pvd> <outdir> <melanin|nutrient>
#        pvpython viewer/paraview_fields.py verify <state.pvsm> <report.json>
#
# Coordinates are lattice sites and both fields are model-scale. Nothing here attaches
# a physical unit, and a changing nutrient image is not evidence that the uptake law is
# Monod: that is established by the model's own tests.
from paraview.simple import *
import sys, os, json, math

FIELDS = {"melanin": "melanin (model field units)",
          "nutrient": "nutrient (model concentration scale)"}
WIDTH, HEIGHT = 1200, 1000          # per view; the layout PNG is two views wide
# Viewport crop per view, as fractions, clear of the caption (top left) and the
# scalar bar (right edge). The image controls read only this region.
CROP = (0.15, 0.70, 0.20, 0.85)     # x0, x1, y0, y1 of each view


def no_raytracing(v):
    # Probe only EnableRayTracing: hasattr(view, "EnableOSPRay") raises on 6.1.1.
    try:
        v.EnableRayTracing = 0
        return "EnableRayTracing=0"
    except Exception as e:
        return "raytracing property unavailable: %s" % e


def view_crop_stats(png, n_views=2):
    """Per view: fraction of crop pixels that are not the white background, and a
    hash of the crop. Text and scalar bars are outside the crop by construction."""
    from vtkmodules.vtkIOImage import vtkPNGReader
    from vtkmodules.util.numpy_support import vtk_to_numpy
    import hashlib
    r = vtkPNGReader(); r.SetFileName(png); r.Update()
    img = r.GetOutput(); w, h, _ = img.GetDimensions()
    px = vtk_to_numpy(img.GetPointData().GetScalars()).reshape(h, w, -1)[:, :, :3]
    out = []
    for k in range(n_views):
        x0, x1 = k * w // n_views, (k + 1) * w // n_views
        vw = x1 - x0
        c = px[int(CROP[2] * h):int(CROP[3] * h),
               x0 + int(CROP[0] * vw):x0 + int(CROP[1] * vw)]
        nonbg = float((c.min(axis=2) < 245).mean())
        out.append({"nonbackground_fraction": round(nonbg, 6),
                    "crop_sha256": hashlib.sha256(c.tobytes()).hexdigest()})
    return out


def field_string(reader, name):
    d = servermanager.Fetch(reader)
    a = d.GetFieldData().GetAbstractArray(name) if d is not None else None
    return a.GetValue(0) if a is not None and a.GetNumberOfValues() else None


def render(pvd, outdir, array):
    os.makedirs(outdir, exist_ok=True)
    reader = OpenDataFile(pvd)
    reader.UpdatePipelineInformation()
    ts = list(reader.TimestepValues)
    if not ts:
        sys.exit("no timesteps in %s" % pvd)
    if array not in reader.CellData.keys():
        sys.exit("%s is not CellData in %s (cell arrays: %s; point arrays: %s)"
                 % (array, pvd, list(reader.CellData.keys()), list(reader.PointData.keys())))

    # One range for the whole trajectory, measured before anything is drawn, so a
    # colour change between frames is a change in the field and not a rescale.
    lo, hi = math.inf, -math.inf
    for t in ts:
        reader.UpdatePipeline(t)
        a, b = reader.CellData[array].GetRange()
        lo, hi = min(lo, a), max(hi, b)
    if not (math.isfinite(lo) and math.isfinite(hi)):
        sys.exit("nonfinite global range for %s: %r..%r" % (array, lo, hi))
    constant = hi <= lo
    rng = (lo - 0.5, lo + 0.5) if constant else (lo, hi)

    reader.UpdatePipeline(ts[0])
    b = reader.GetDataInformation().GetBounds()
    centre = [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2]
    extent = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])

    LoadPalette(paletteName="WhiteBackground")
    v1 = CreateRenderView(); v2 = CreateRenderView()
    layout = CreateLayout(name="%s field views" % array)
    AssignViewToLayout(view=v1, layout=layout, hint=0)
    layout.SplitHorizontal(0, 0.5)
    AssignViewToLayout(view=v2, layout=layout, hint=2)
    rt = {}
    for name, v in (("slice", v1), ("volume", v2)):
        v.ViewSize = [WIDTH, HEIGHT]; v.OrientationAxesVisibility = 1
        rt[name] = no_raytracing(v)

    lut = GetColorTransferFunction(array)
    lut.ApplyPreset("Viridis", True)
    otf = GetOpacityTransferFunction(array)

    def lock_range():
        # ColorBy and representation changes rescale to the data in view, so the
        # global range is (re)applied after every one of them.
        lut.AutomaticRescaleRangeMode = "Never"
        lut.RescaleTransferFunction(*rng)
        otf.RescaleTransferFunction(*rng)
        try:
            otf.AutomaticRescaleRangeMode = "Never"
        except Exception:
            pass

    # Slice: the cylinder's axis is parallel to z, so a plane containing it has its
    # normal in x/y. Geometric, through the domain centre, not chosen from the data.
    thr = Threshold(Input=reader, registrationName="Interior_mask_1")
    thr.Scalars = ["CELLS", "interior_mask"]
    thr.LowerThreshold = 1; thr.UpperThreshold = 1
    sl = Slice(Input=thr, registrationName="Axial_slice")
    sl.SliceType = "Plane"; sl.SliceType.Origin = centre; sl.SliceType.Normal = [1.0, 0.0, 0.0]
    d1 = Show(sl, v1); d1.SetRepresentationType("Surface")
    ColorBy(d1, ("CELLS", array))
    d1.SetScalarBarVisibility(v1, True)

    # Volume: the original ImageData and the original cell array first. A GPU mode is
    # a request; what the representation accepted is what gets reported.
    d2 = Show(reader, v2)
    ColorBy(d2, ("CELLS", array))
    d2.SetRepresentationType("Volume")
    # Accepted means it reads back as set; an unknown value may only warn.
    tried = []
    for mode in ("GPU Based", "Smart"):
        try:
            d2.VolumeRenderingMode = mode
        except Exception as e:
            tried.append("%s: refused (%s)" % (mode, e)); continue
        got = str(d2.GetProperty("VolumeRenderingMode").GetData())
        if got == mode:
            tried.append("%s: accepted" % mode); break
        tried.append("%s: read back as %s" % (mode, got))
    accepted = str(d2.GetProperty("VolumeRenderingMode").GetData())
    # Sites outside the cylinder are pinned near C_wall, the field's maximum, and
    # drawn as an opaque box that hides the interior. Opacity therefore comes from
    # interior_mask (0 outside, uniform inside) while colour stays on the field.
    opacity = "field (global range)"
    try:
        d2.UseSeparateOpacityArray = 1
        d2.OpacityArrayName = ["CELLS", "interior_mask"]
        mask_otf = d2.ScalarOpacityFunction
        # (x, opacity, midpoint, sharpness) per point: a step from 0 to 0.06 at 0.5.
        mask_otf.Points = [0.0, 0.0, 0.5, 0.0, 0.49, 0.0, 0.5, 0.0,
                           0.51, 0.06, 0.5, 0.0, 1.0, 0.06, 0.5, 0.0]
        opacity = "separate array interior_mask: 0 outside, 0.06 inside"
    except Exception as e:
        opacity = "field (global range); separate opacity array refused: %s" % e
    d2.SetScalarBarVisibility(v2, True)
    lock_range()
    for v in (v1, v2):
        bar = GetScalarBar(lut, v); bar.Title = FIELDS[array]; bar.ComponentTitle = ""

    cam1 = v1.GetActiveCamera()
    cam1.SetFocalPoint(*centre); cam1.SetPosition(centre[0] + 3 * extent, centre[1], centre[2])
    cam1.SetViewUp(0, 0, 1); v1.ResetCamera(False)
    cam2 = v2.GetActiveCamera()
    cam2.SetFocalPoint(*centre); cam2.SetPosition(centre[0] + 2 * extent, centre[1] - 2 * extent,
                                                  centre[2] + 1.5 * extent)
    cam2.SetViewUp(0, 0, 1); v2.ResetCamera(False)

    labels = []
    for v in (v1, v2):
        t = Text(Text=""); td = Show(t, v); td.FontSize = 24; td.Color = [0, 0, 0]
        td.WindowLocation = "Upper Left Corner"; labels.append(t)

    scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps()

    def shoot(t, path):
        scene.AnimationTime = t
        for lab in labels:
            lab.Text = "MCS %d — %s\nspatial coordinates: lattice sites" % (int(t), array)
        Render(v1); Render(v2)
        pts = list(lut.RGBPoints)
        if abs(pts[0] - rng[0]) > 1e-12 or abs(pts[-4] - rng[1]) > 1e-12:
            sys.exit("LUT drifted from the global range at MCS %d: %r..%r" % (int(t), pts[0], pts[-4]))
        SaveScreenshot(path, layout, ImageResolution=[2 * WIDTH, HEIGHT])
        return view_crop_stats(path)

    # Fall back to point data for the VOLUME only, and only if the cell-data volume
    # actually drew nothing; the slice stays on the original cell data either way.
    association, fallback = "CELLS", None
    probe = os.path.join(outdir, ".%s_volume_probe.png" % array)
    if shoot(ts[-1], probe)[1]["nonbackground_fraction"] == 0.0:
        Hide(reader, v2)
        c2p = CellDatatoPointData(Input=reader, registrationName="Volume_cell_to_point")
        d2 = Show(c2p, v2); ColorBy(d2, ("POINTS", array)); d2.SetRepresentationType("Volume")
        d2.SetScalarBarVisibility(v2, True)
        lock_range()
        association, fallback = "POINTS", "CellDatatoPointData"
        accepted = str(d2.GetProperty("VolumeRenderingMode").GetData())
        opacity = "field (global range); mask not applied after fallback"
    os.remove(probe)

    frames = {}
    for t in ts:
        png = os.path.join(outdir, "%s_mcs%06d.png" % (array, int(t)))
        frames[int(t)] = {"png": os.path.basename(png), "views": shoot(t, png)}
        print("MCS %d  slice %.4f  volume %.4f nonbackground" % (
            int(t), frames[int(t)]["views"][0]["nonbackground_fraction"],
            frames[int(t)]["views"][1]["nonbackground_fraction"]))

    state = os.path.join(outdir, "%s_view.pvsm" % array)
    SaveState(state)
    cam = lambda c: {"position": list(c.GetPosition()), "focal_point": list(c.GetFocalPoint()),
                     "view_up": list(c.GetViewUp())}
    report = {
        "pvd": os.path.abspath(pvd), "git_sha": field_string(reader, "git_sha"),
        "paraview": GetParaViewSourceVersion(), "array": array, "timesteps": [int(t) for t in ts],
        "slice_association": "CELLS", "volume_association": association,
        "volume_fallback": fallback, "global_range": [lo, hi], "constant_field": constant,
        "lut_range": list(rng), "lut_and_opacity_rescale": "Never",
        "slice": {"origin": centre, "normal": [1.0, 0.0, 0.0], "input": "interior_mask == 1"},
        "cameras": {"slice": cam(cam1), "volume": cam(cam2)},
        "raytracing_at_render": rt,
        "volume_rendering_mode": {"attempts": tried, "accepted": accepted},
        "volume_opacity": opacity,
        "crop": {"x0": CROP[0], "x1": CROP[1], "y0": CROP[2], "y1": CROP[3]},
        "frames": frames, "pvsm": os.path.basename(state),
    }
    path = os.path.join(outdir, "%s_render_report.json" % array)
    json.dump(report, open(path, "w"), indent=2)
    print("wrote %d frames, %s and %s" % (len(ts), os.path.basename(state), os.path.basename(path)))


def verify(state, report_path):
    """Reopen a saved state in this (fresh) process and check it against its report."""
    rep = json.load(open(report_path))
    LoadState(state)
    fails = []
    def ck(name, ok, got):
        print(("  PASS  " if ok else "  FAIL  ") + name + "  ->  " + str(got))
        ok or fails.append(name)
    srcs = GetSources()
    reader = next((v for v in srcs.values() if v.GetXMLName() == "PVDReader"), None)
    ck("a PVD reader is present", reader is not None, reader is not None)
    if reader is None:
        sys.exit(1)
    ck("reader points at the reported .pvd", reader.FileName == rep["pvd"], reader.FileName)
    reader.UpdatePipelineInformation()
    ts = [int(t) for t in reader.TimestepValues]
    ck("timesteps match the report", ts == rep["timesteps"], ts)
    ck("the array is still CellData at the source", rep["array"] in reader.CellData.keys(),
       list(reader.CellData.keys()))
    ck("two render views", len(GetRenderViews()) == 2, len(GetRenderViews()))
    by = {k[0]: v for k, v in srcs.items()}
    sl = by.get("Axial_slice")
    ck("slice present", sl is not None, sl is not None)
    if sl is not None:
        ck("slice origin", all(abs(a - b) < 1e-9 for a, b in zip(sl.SliceType.Origin, rep["slice"]["origin"])),
           list(sl.SliceType.Origin))
        ck("slice normal", list(sl.SliceType.Normal) == rep["slice"]["normal"], list(sl.SliceType.Normal))
    lut = GetColorTransferFunction(rep["array"])
    pts = list(lut.RGBPoints)
    ck("LUT range survives", abs(pts[0] - rep["lut_range"][0]) < 1e-12 and abs(pts[-4] - rep["lut_range"][1]) < 1e-12,
       (pts[0], pts[-4]))
    ck("LUT rescale locked", lut.AutomaticRescaleRangeMode == "Never", lut.AutomaticRescaleRangeMode)
    has_c2p = "Volume_cell_to_point" in by
    ck("volume association matches the report", has_c2p == (rep["volume_association"] == "POINTS"),
       "cell-to-point filter present" if has_c2p else "original cell data")
    print("\n%d failures" % len(fails))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "verify":
        verify(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[3] in FIELDS:
        render(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        sys.exit("usage: pvpython paraview_fields.py <lattice.pvd> <outdir> <melanin|nutrient>\n"
                 "       pvpython paraview_fields.py verify <state.pvsm> <report.json>")
