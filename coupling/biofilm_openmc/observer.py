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

from .viewer import read_layer, read_manifest

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
WALL = -1               # outside the biological domain; NOT a membrane material


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
        colour_by=_colour_by(l)) for l in manifest["layers"]]


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
    """Draw one layer. Refuses to hide a non-quotable layer's banner.

    `show_banner=False` on a derived layer raises rather than warns: the banner
    is the only thing distinguishing a real 4x dose field from a coarse field
    broadcast to look like one, and a caller who suppresses it has removed the
    distinction from the screen.
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
    p.add_mesh(image, scalars=layer_name,
               scalar_bar_args={"title": f"{layer_name} [{layer.unit}]"})
    if layer.banner and show_banner:
        p.add_text(layer.banner, position="upper_edge", font_size=8)
    return p
