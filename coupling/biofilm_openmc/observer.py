"""Read-only display of a viewer bundle. Renders; never computes.

THE SEPARATION IS THE POINT. This module turns a bundle into things a renderer
can draw and does nothing else: it evaluates no metric, applies no threshold,
returns no verdict, and writes no file. `feedback_uq` and `synthetic_gate` are
not imported, and a test asserts it — a display layer that can compute is a
display layer that will eventually be quoted.

AND ALMOST NONE OF IT IS RENDERING. Choosing which layers may be drawn, what
each must be labelled with, how a categorical field maps to colour, and what the
provenance panel says are all decisions about DATA. They are pure numpy here and
testable in the bare tier. Only `to_image_data` and `plot_layer` touch PyVista,
and they import it function-locally, so the whole contract is exercised on a
runner with no renderer installed.

WHAT IT REFUSES TO DO SILENTLY. A bundle marks every layer `native` and
`authoritative_for_quantitation`. An upsampled dose field is legitimate to draw
— overlaying dose on the label grid is exactly what a viewer is for — and is not
legitimate to read numbers off. `display_plan` therefore attaches a REQUIRED
banner to every non-quotable layer, and `plot_layer` refuses to render one
unless the caller has passed the banner through to a title. The alternative is a
picture that looks like per-parcel dose and is one coarse value repeated.

Units travel with the layer for the same reason: the transport fields are Gy per
SOURCE PARTICLE, and a colour bar labelled Gy/s has invented an emission rate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .viewer import (BOOLEAN, CATEGORICAL, OUT_OF_DOMAIN, UNDECLARED,
                     is_undeclared, read_layer, read_manifest)

# Same palette and labels as the serial model's figure section, carried across
# from the CairoMakie prototype on `feat/visualize-3d`. That branch is otherwise
# unusable — it predates the whole coupling stack and deletes it — but the
# species colours are a real convention and should not be reinvented per viewer.
SPECIES_COLOURS = ("#e6194b", "#3cb44b", "#4363d8", "#f58231",
                   "#911eb4", "#42d4f4", "#f032e6")
SPECIES_LABELS = ("C. neoformans", "D. radiodurans", "C. sphaerospermum",
                  "B. subtilis", "A. niger", "S. oneidensis", "O. intermedium")

# Orbit parameters from the same prototype, kept so successive viewers of this
# model look like one another.
DEFAULT_CAMERA = {"azimuth": 0.4, "elevation": 0.5}

BACKGROUND = 0          # cell_id background, per the exchange schema
WALL = OUT_OF_DOMAIN    # outside the biological domain; NOT a membrane material


@dataclass(frozen=True)
class DisplayLayer:
    """One drawable layer, and everything a viewer must say about it."""
    name: str
    grid_id: str
    unit: str
    semantic_kind: str
    quotable: bool
    banner: str             # empty when quotable; REQUIRED on screen otherwise
    derivation_note: str    # how a derived layer was made; may be quotable
    colour_by: str          # "species" | "scalar" | "mask"
    background: float | None | str = UNDECLARED   # producer's, not renderer's
    occupancy_from: str | None = None  # the layer that says where biomass is

    def as_dict(self) -> dict:
        return {"name": self.name, "grid_id": self.grid_id, "unit": self.unit,
                "semantic_kind": self.semantic_kind, "quotable": self.quotable,
                "banner": self.banner, "derivation_note": self.derivation_note,
                "colour_by": self.colour_by}


def _banner(layer: dict) -> str:
    """What a non-quotable layer must carry on screen.

    Named for the mechanism rather than 'derived', because the reader needs to
    know WHY the numbers are not readable, and the two reasons differ: an
    upsampled field was never measured at this resolution, whereas a display
    resampling may be exact and is still not the grid it was computed on.
    """
    if layer["authoritative_for_quantitation"]:
        return ""
    derivation = layer.get("derivation") or "derived"
    source = layer.get("source_grid_id") or "another grid"
    if derivation == "upsampled_coarse_dose":
        return (f"NOT QUANTITATIVE — upsampled from {source}. Every value here "
                "is one coarse-bin value repeated; the detail is the label "
                "grid's, not the dose's.")
    if derivation == "display_resampling":
        return (f"NOT QUANTITATIVE — resampled from {source} for display only.")
    return f"NOT QUANTITATIVE — {derivation} from {source}."


def _derivation_note(layer: dict) -> str:
    """Descriptive, not a warning. A layer can be derived AND quotable — an
    exact mass-weighted reduction is both — and stamping that with the same
    alarming banner an upsampling gets would teach a reader to ignore banners.
    """
    if layer["native"]:
        return ""
    source = layer.get("source_grid_id") or "another grid"
    return f"{layer.get('derivation') or 'derived'} from {source}"


def _colour_by(layer: dict) -> str:
    if layer["semantic_kind"] == "boolean":
        return "mask"
    if layer["semantic_kind"] != "categorical":
        return "scalar"
    # Only a species field maps onto the palette; cell and lineage ids are
    # unbounded labels and get a categorical colormap, not seven fixed colours.
    return "species" if "species" in layer["name"] else "scalar"


def display_plan(manifest: dict) -> list[DisplayLayer]:
    """Every layer, with its banner and colour rule resolved.

    Returns ALL layers including the non-quotable ones: hiding them would push a
    viewer into recomputing the overlay itself, which is how the display layer
    grows a metric.
    """
    return [DisplayLayer(
        name=l["name"], grid_id=l["grid_id"], unit=l["unit"],
        semantic_kind=l["semantic_kind"],
        quotable=bool(l["authoritative_for_quantitation"]),
        banner=_banner(l), derivation_note=_derivation_note(l),
        colour_by=_colour_by(l),
        background=(l.get("background", UNDECLARED)
                    if l.get("background", UNDECLARED) is None
                    or is_undeclared(l.get("background", UNDECLARED))
                    else float(l["background"])),
        occupancy_from=l.get("occupancy_from")) for l in manifest["layers"]]


def grid_geometry(manifest: dict, grid_id: str) -> dict:
    """Origin, spacing and cell counts, in cm, for one grid.

    A renderer needs POINT dimensions (cells + 1) while the bundle stores CELL
    counts, and confusing them shifts every field by half a voxel — visible only
    as a misregistration between the label grid and a 4x dose grid, which is
    exactly the overlay this bundle exists to support.
    """
    for grid in manifest["grids"]:
        if grid["grid_id"] == grid_id:
            shape = tuple(int(n) for n in grid["shape_xyz"])
            return {"grid_id": grid_id,
                    "cell_shape_xyz": shape,
                    "point_shape_xyz": tuple(n + 1 for n in shape),
                    "origin_cm": tuple(float(v) for v in grid["origin_cm"]),
                    "spacing_cm": tuple(float(v) for v in grid["spacing_cm"]),
                    "extent_cm": tuple(float(v) for v in grid["extent_cm"]),
                    "material_resolution_grid":
                        grid.get("material_resolution_grid") or None}
    raise KeyError(f"grid {grid_id!r} is not in this bundle")


def provenance_panel(manifest: dict) -> list[tuple]:
    """What the viewer must show about where the data came from.

    Ordered, and `target_calibration` is first: a reader who takes one glance
    should learn whether these numbers are allowed to mean anything about a real
    system before they learn anything else about them.
    """
    prov = manifest.get("provenance", {})
    rows = [("target_calibration", str(prov.get("target_calibration", "unset"))),
            ("reference_system_id", str(prov.get("reference_system_id", "unset"))),
            ("evidence_policy", str(prov.get("evidence_policy", "unset"))),
            ("system_provenance", str(prov.get("system_provenance", "unset")))]
    for key in ("study", "openmc_version", "git_commit", "particles", "batches"):
        if key in prov:
            rows.append((key, str(prov[key])))
    rows.append(("grids", ", ".join(g["grid_id"] for g in manifest["grids"])))
    quotable = sum(1 for l in manifest["layers"]
                   if l["authoritative_for_quantitation"])
    rows.append(("layers", f"{len(manifest['layers'])} ({quotable} quantitative)"))
    return rows


def species_legend(present_ids=None) -> list[tuple]:
    """(id, label, colour) for the species actually present, or all seven.

    Filtered by default because a legend listing seven organisms for a snapshot
    containing two is a legend that misdescribes the picture.
    """
    ids = (range(1, len(SPECIES_LABELS) + 1) if present_ids is None
           else sorted({int(v) for v in np.asarray(present_ids).ravel()
                        if 1 <= int(v) <= len(SPECIES_LABELS)}))
    return [(i, SPECIES_LABELS[i - 1], SPECIES_COLOURS[i - 1]) for i in ids]


def species_palette(legend) -> tuple[list[str], tuple[float, float]]:
    """Colormap and colour limits giving each species id its OWN palette slot.

    ONE SLOT PER ID ACROSS THE SPANNED RANGE, not one per species present.
    Handing a renderer the three colours of {1, 2, 7} with limits spanning 1..7
    makes it distribute three colours over six units, so species 2 draws in
    whatever colour falls at that fraction while the legend labels it with the
    second palette entry. The picture and its key then disagree -- which is
    precisely the failure a fixed palette exists to prevent, arriving through
    the code that implements the fixed palette.

    Half-unit padding puts each integer id in the middle of its own band, so a
    species keeps its colour no matter which others are present.
    """
    if not legend:
        # AN EMPTY SNAPSHOT IS A LEGITIMATE FIELD, not a broken one: a species
        # layer of pure background has no ids present, `species_legend` returns
        # nothing, and min() over nothing raised ValueError before the layer
        # could draw. Hand back the full fixed palette -- there is nothing to
        # colour, so nothing can be miscoloured, and the caller stays on one
        # code path.
        return (list(SPECIES_COLOURS), (0.5, len(SPECIES_COLOURS) + 0.5))
    lo = min(i for i, _, _ in legend)
    hi = max(i for i, _, _ in legend)
    return ([SPECIES_COLOURS[i - 1] for i in range(lo, hi + 1)],
            (lo - 0.5, hi + 0.5))


def occupied_mask(cell_id) -> np.ndarray:
    """Voxels holding biomass.

    `> 0` and not `!= 0`: the schema uses -1 for OUTSIDE the biological domain
    and 0 for background, so a `!= 0` test draws the wall as if it were biomass.
    """
    return np.asarray(cell_id) > BACKGROUND


# ---------------------------------------------------------------- rendering

def to_image_data(path, layer_name, manifest=None):
    """One layer as a `pyvista.ImageData`, positioned in centimetres.

    The array is attached as CELL data, never point data. A voxel field is
    piecewise constant over its cell; interpolating it to points would smooth
    the material boundaries the model asserts are sharp, and would make a
    categorical field meaningless outright.
    """
    import pyvista as pv

    doc = manifest or read_manifest(path)
    data, attrs = read_layer(path, layer_name)
    geometry = grid_geometry(doc, attrs["grid_id"])

    image = pv.ImageData(dimensions=geometry["point_shape_xyz"],
                         spacing=geometry["spacing_cm"],
                         origin=geometry["origin_cm"])
    # Fortran order: the bundle stores logical (x,y,z) and VTK wants x fastest.
    image.cell_data[layer_name] = np.asarray(data).ravel(order="F")
    return image


def plot_layer(path, layer_name, *, plotter=None, show_banner=True,
               manifest=None):
    """Draw one layer so its INTERIOR is visible, with its declared colour rule.

    `add_mesh` on an `ImageData` renders only the exterior surface of the block.
    For every bundle this repository produces that surface is background — the
    cube circumscribing the cylinder is 21.5% corner void — so a naive call
    yields an opaque shell with the biofilm hidden inside it. Two paths, chosen
    by the layer's own semantic kind:

    - CATEGORICAL and BOOLEAN layers are drawn as a surface, hiding what the
      PRODUCER declared to be absent — never what this function guessed. Three
      cases, because a value does not always know its own meaning: a layer that
      names an `occupancy_from` layer is masked by THAT field, one that declares
      a `background` value hides exactly that value, and one that declares
      neither is drawn in full.

      `generation` is why. It is 0 for a founder and 0 for empty lattice sites,
      since `export_checkpoint.jl` zero-fills and skips unoccupied voxels — so
      thresholding at 0.5 deleted the founding cohort, and drawing everything
      rendered the void as generation-0 biomass. Neither is a background value;
      the answer is `cell_id`. A species layer gets the seven fixed colours
      rather than a continuous colormap, which is the whole point of having a
      palette: the same organism must not change colour because a different
      subset is present.
    - INTENSIVE and EXTENSIVE layers are volume-rendered, because thresholding a
      dose field at an arbitrary value would hide exactly the structure the
      viewer exists to show.

    Refuses to suppress a non-quotable layer's banner, and shows the neutral
    derivation note for a layer that is derived but still quotable — otherwise
    an exact reduction is indistinguishable on screen from a native measurement
    and silently loses its source grid.
    """
    import pyvista as pv

    doc = manifest or read_manifest(path)
    plan = {l.name: l for l in display_plan(doc)}
    if layer_name not in plan:
        raise KeyError(f"layer {layer_name!r} is not in this bundle")
    layer = plan[layer_name]
    if layer.banner and not show_banner:
        raise ValueError(
            f"layer {layer_name!r} is not quantitative and its banner may not "
            f"be suppressed: {layer.banner}")

    image = to_image_data(path, layer_name, manifest=doc)
    p = plotter or pv.Plotter()
    title = f"{layer_name} [{layer.unit}]"

    if layer.semantic_kind in (CATEGORICAL, BOOLEAN):
        if layer.occupancy_from:
            # THE VALUE CANNOT DISAMBIGUATE ITSELF, so do not ask it to.
            # `generation` is 0 for a founder AND for empty space, so mask by
            # the layer the producer named as occupancy instead of by value.
            occ, _ = read_layer(path, layer.occupancy_from)
            # `occupied_mask`, not `!= 0`: the schema uses -1 for outside the
            # biological domain, so a truthiness test draws the wall as biomass.
            image.cell_data["_occupied"] = (
                occupied_mask(occ).ravel(order="F")).astype(np.uint8)
            occupied = image.threshold(0.5, scalars="_occupied")
            occupied.set_active_scalars(layer_name)
        elif layer.background is None:
            occupied = image.threshold(scalars=layer_name)   # keeps everything
        elif is_undeclared(layer.background):
            # write_bundle refuses this, so a bundle in hand cannot reach here.
            raise ValueError(
                f"layer {layer_name!r} declares no absence semantics")
        else:
            # ONE BRANCH FOR EVERY DECLARED BACKGROUND. This was split by
            # `background == 0` and the two halves disagreed about the
            # sentinel: the zero half excluded OUT_OF_DOMAIN, the other
            # dropped only the declared background via an inverted equality
            # band, so a layer declaring background 5 and carrying -1 drew
            # the out-of-domain wall as label data. Nothing refuses that
            # encoding either -- `_has_unknown_sentinel` runs in
            # `bundle_problems` ONLY under `if layer.occupancy_from is not
            # None`, so a directly-rendered layer is never validated against
            # the domain contract at all.
            #
            # NOT `occupied_mask` (`> 0`): that assumes the same contract,
            # and would silently discard a legitimate negative label (a real
            # -7, not a sentinel) on a layer nothing validated. Exclude ONLY
            # the two values that ARE sentinels, by equality -- the declared
            # background and the schema's reserved OUT_OF_DOMAIN. Anything
            # else, negative or not, is data. Background -1 collapses the two
            # terms, which is correct: it is one value either way.
            raw = np.asarray(image.cell_data[layer_name])
            keep = (raw != layer.background) & (raw != OUT_OF_DOMAIN)
            image.cell_data["_occupied"] = keep.astype(np.uint8)
            occupied = image.threshold(0.5, scalars="_occupied")
            occupied.set_active_scalars(layer_name)
        if occupied.n_cells == 0:
            # AN EMPTY FIELD MUST STILL DRAW, because "nothing is here" is a
            # result a reader needs to see. Thresholding an all-background
            # layer leaves a dataset with no cells and `add_mesh` rejects it,
            # so making `species_palette` empty-safe only moved the crash one
            # line down -- the bare-tier test exercised the palette, not this.
            # Say so on the canvas and fall through to the banner.
            p.add_text(f"{layer_name}: no occupied cells",
                       position="lower_left", font_size=9)
        elif layer.colour_by == "species":
            # THE LEGEND IS BUILT FROM THE CELLS THAT SURVIVE, not from the raw
            # field. Reading the array again gave two independent derivations
            # of one fact, and they disagreed the moment a producer declared an
            # in-range background: with values {1, 9} and background 1, the raw
            # field still offered species 1 to the legend while `occupied` had
            # already dropped it, so the bounds selected nothing and `add_mesh`
            # was handed an empty dataset. With {1, 2} and background 1 the
            # same gap instead LISTED species 1 in a picture that does not
            # contain it.
            #
            # One source. Whatever is in `occupied` is what gets named.
            present = np.asarray(occupied.cell_data[layer_name])
            legend = species_legend(present)
            if not legend:
                # Occupied cells, none of them a species this palette names.
                p.add_text(f"{layer_name}: no valid species ids",
                           position="lower_left", font_size=9)
            else:
                # WHAT IS DRAWN MUST BE WHAT THE LEGEND DESCRIBES.
                # `species_legend` drops an out-of-range label -- 9, when there
                # are seven species -- but the occupied dataset still held that
                # cell, so it rendered in whatever colour its value landed on
                # with nothing in the key to explain it. A viewer showing an
                # organism it cannot name is worse than one showing fewer
                # cells: the reader counts it.
                #
                # The legend's own span is the valid range by construction, so
                # thresholding to it drops exactly the labels the legend
                # dropped, and nothing else.
                lo = min(i for i, _, _ in legend)
                hi = max(i for i, _, _ in legend)
                drawn = occupied.threshold((lo, hi), scalars=layer_name)
                cmap, clim = species_palette(legend)
                p.add_mesh(drawn, scalars=layer_name, show_scalar_bar=False,
                           cmap=cmap, clim=clim)
                p.add_legend([(label, colour) for _, label, colour in legend])
        else:
            p.add_mesh(occupied, scalars=layer_name,
                       scalar_bar_args={"title": title})
    else:
        p.add_volume(image, scalars=layer_name,
                     scalar_bar_args={"title": title})

    if layer.banner and show_banner:
        p.add_text(layer.banner, position="upper_edge", font_size=8)
    elif layer.derivation_note:
        p.add_text(layer.derivation_note, position="upper_edge", font_size=7)
    return p


def plot_overlay(path, label_layer, dose_layer, *, plotter=None,
                 cylinder=None, manifest=None):
    """Draw a label layer and a dose layer in ONE scene, in physical space.

    THIS IS THE DIAGNOSTIC, NOT THE CHECK. `tests/test_grid_coregistration.py`
    is what fails when the two grids disagree; this is what a reader opens
    afterwards to see HOW they disagree. A picture nobody is obliged to look at
    cannot be a gate, and pretending otherwise would be the
    check-that-cannot-fail this repository keeps finding.

    It composes rather than reimplements: `plot_layer` already accepts a
    `plotter`, and both layers are positioned by `to_image_data` from their own
    grid's declared `origin_cm`/`spacing_cm`. So the two grids land in the same
    world frame BECAUSE the bundle says where they are -- if one is offset, the
    picture shows it offset, which is the entire point. Nothing here recomputes
    a coordinate.

    The dose layer keeps whatever banner `display_plan` gives it. An upsampled
    dose field drawn over labels is exactly the case `_banner` exists for, and
    `plot_layer` refuses to suppress it.

    `cylinder`, when given, is `(x0, y0, r, z_lo, z_hi)` in cm -- the CSG
    boundary the rectilinear lattice sits inside. Voxels at the curved wall are
    partially inside, which is why `mesh_material_volumes` raytraces them
    instead of counting bins; drawing the surface beside the lattice is how a
    reader sees which voxels those are. It is wireframe and unlit on purpose:
    it is geometry being quoted, not a field being read.
    """
    import pyvista as pv

    doc = manifest or read_manifest(path)
    p = plotter or pv.Plotter()

    plot_layer(path, label_layer, plotter=p, manifest=doc)
    plot_layer(path, dose_layer, plotter=p, manifest=doc)

    if cylinder is not None:
        x0, y0, r, z_lo, z_hi = (float(v) for v in cylinder)
        surface = pv.Cylinder(center=(x0, y0, 0.5 * (z_lo + z_hi)),
                              direction=(0.0, 0.0, 1.0),
                              radius=r, height=(z_hi - z_lo),
                              resolution=64, capping=False)
        p.add_mesh(surface, style="wireframe", opacity=0.35,
                   lighting=False, color="grey")
        p.add_text("grey wireframe: CSG boundary, not data",
                   position="lower_right", font_size=7)
    return p
