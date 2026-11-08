"""The multi-grid viewer bundle: several native grids, one coordinate system.

WHY THIS EXISTS. A viewer that puts every field on the CPM lattice invents
detail it does not have. `synthetic_e2e.py` upsamples a coarse dose field to
lattice resolution BEFORE attribution, so every per-cell and per-lineage number
it produces is a broadcast coarse plane joined to fine labels. Displayed on the
lattice it looks like per-parcel dose; it is one coarse value repeated across
every parcel in a bin. Since the subvoxel refinement study, the opposite case
also exists — a genuinely 2x or 4x native OpenMC dose field, which the CPM grid
cannot represent at all.

So the bundle carries GRIDS, not one grid. Each layer names the grid it lives on
and whether it is NATIVE there. Resampling stays legal, because a viewer must be
able to overlay things, but the DIRECTION of the resampling decides whether the
result may still be quoted:

    upsampling  broadcasts one value across bins that were never separately
                measured. It INVENTS information, so it may be drawn and may
                never be quoted.
    reduction   block-sums an extensive quantity, or takes a mass-weighted mean
                of an intensive one. It DESTROYS information and invents none —
                the refinement study confirmed a 4x dose field summed back to
                the lattice reproduces the native lattice field to nine
                significant figures — so it may be quoted.

Both implications are enforced rather than documented, and a "reduction" onto a
grid that is not coarser is refused, because that label is the only thing
standing between upsampling and being quoted.

SEMANTIC KIND DECIDES HOW A LAYER MAY BE RESAMPLED, and getting it wrong is off
by factor**3 rather than visibly broken. The distinction is the same one
`calibration/biofilm_calibration/spatial/field.py` states from the other side:

    extensive    energy, mass          reduce by block SUM
    intensive    dose, volume fraction reduce by a WEIGHTED mean (see
                                       `mesh.coarsen_ratio`) — never a plain one
    categorical  cell_id, species_id   never resampled by arithmetic at all;
                                       a mean of two labels is a third label
    boolean      interior_mask, Omega_b  any/all, declared

UNITS ARE NOT DECORATION. The pilot's fields are Gy per SOURCE PARTICLE, not
Gy/s: they differ by the source activity, which for the synthetic system does
not exist. A bundle that labels them Gy/s has fabricated an emission rate.

Pure numpy + h5py, so this lives in the bare tier alongside everything else that
reads a snapshot. Nothing here imports a rendering library, and nothing here
computes a metric.

NO `viewer` EXTRA IS DECLARED YET, deliberately. The bundle needs nothing this
package does not already depend on, and an extra that installs nothing is a
promise about a renderer that does not exist. It arrives with the renderer.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from .snapshot import _h5py_to_logical, _logical_to_h5py

SCHEMA_VERSION = 1

# How a layer behaves under resampling. See the module docstring.
EXTENSIVE = "extensive"
INTENSIVE = "intensive"
# THE DEFAULT MUST NOT BE A STATEMENT. `background=None` says "every cell
# carries information, hide nothing" -- true for `generation`, and wrong for
# `omega_b`. If it were also the default, a producer who simply forgot would
# make that claim silently, which is how two layers came to render their empty
# space as data. UNDECLARED is the default and is refused at write time.
UNDECLARED = "undeclared"

# The exchange schema's out-of-domain sentinel: 0 is empty space INSIDE the
# biological domain, -1 is outside it. `observer.occupied_mask` tests `> 0`, so
# it treats both as absent -- which makes this constant half of the occupancy
# contract, not a detail of one renderer. It lives here, beside the Layer that
# has to validate against it.
OUT_OF_DOMAIN = -1


def is_undeclared(value) -> bool:
    """Whether an absence declaration is missing.

    BY VALUE, NEVER BY IDENTITY. The sentinel has to survive a round trip
    through the manifest, and `json.loads('"undeclared"')` is a DIFFERENT string
    object that compares equal -- so `value is UNDECLARED` accepted a categorical
    layer with no absence semantics, and `plot_layer` then sent the string into
    thresholding as though it were a number. Every representation must be
    refused the same way, so every caller asks this one question.
    """
    return isinstance(value, str) and value == UNDECLARED

CATEGORICAL = "categorical"
BOOLEAN = "boolean"
SEMANTIC_KINDS = frozenset({EXTENSIVE, INTENSIVE, CATEGORICAL, BOOLEAN})

# A layer that is not native must say how it was produced, and the split below
# is the one that decides whether it may still be quoted.
#
# UPSAMPLING INVENTS INFORMATION. Broadcasting one coarse value across bins that
# were never separately measured produces a field that LOOKS resolved and is
# not, which is exactly the `dose_attribution.h5` defect. Never authoritative,
# whatever the arithmetic.
#
# REDUCTION DESTROYS INFORMATION AND INVENTS NONE. A block sum of an extensive
# quantity, or a mass-weighted mean of an intensive one, is exact: the
# refinement study confirmed a 4x dose field summed back to the lattice
# reproduces the native lattice field to nine significant figures. So a
# correctly reduced layer MAY be quoted — the earlier blanket rule that
# "derived means unquotable" was blunter than the physics.
UPSAMPLING_DERIVATIONS = frozenset({"upsampled_coarse_dose", "display_resampling"})
REDUCING_DERIVATIONS = frozenset({"block_sum", "mass_weighted_mean"})
# Neither a resampling nor a measurement: Omega_b from CSG volumes, a difference
# field, a debiased effect map. Quotable if the producer says so.
COMPUTED = "computed"
DERIVATIONS = UPSAMPLING_DERIVATIONS | REDUCING_DERIVATIONS | {COMPUTED}


@dataclass(frozen=True)
class Grid:
    """One sampling of the shared physical volume.

    `origin_cm` and `spacing_cm` are required rather than derived: two grids are
    only comparable if both state where they start and how big their bins are,
    and a viewer that infers either will eventually infer it wrong.

    `material_resolution_grid` is the honesty field for a refined dose grid. A
    4x tally is native transport data, but the MATERIAL it integrates is still
    piecewise constant at the CPM pitch, so a viewer must not present its
    sub-voxel structure as biological structure.
    """
    grid_id: str
    shape_xyz: tuple
    origin_cm: tuple
    spacing_cm: tuple
    material_resolution_grid: str | None = None
    note: str = ""

    @property
    def extent_cm(self) -> tuple:
        return tuple(float(n) * float(s)
                     for n, s in zip(self.shape_xyz, self.spacing_cm))

    def as_dict(self) -> dict:
        return {"grid_id": self.grid_id,
                "shape_xyz": [int(n) for n in self.shape_xyz],
                "origin_cm": [float(v) for v in self.origin_cm],
                "spacing_cm": [float(v) for v in self.spacing_cm],
                "extent_cm": [float(v) for v in self.extent_cm],
                "material_resolution_grid": self.material_resolution_grid,
                "note": self.note}


@dataclass(frozen=True)
class Layer:
    """One field on one grid, with everything a reader needs to refuse it."""
    name: str
    grid_id: str
    unit: str
    semantic_kind: str
    data: np.ndarray                      # logical (x, y, z)
    native: bool = True
    authoritative_for_quantitation: bool = True
    source_grid_id: str | None = None     # where a derived layer came from
    derivation: str | None = None
    # THE PRODUCER DECLARES THIS; A RENDERER MAY NOT ASSUME IT. Which value, if
    # any, means "nothing is here". `cell_id` 0 is genuinely empty space, but
    # `generation` 0 is a founder and thresholding it away deletes the first
    # cohort from the picture. Only the code that built the field knows which it
    # is, so `None` means every cell carries information and none may be hidden.
    background: float | None | str = UNDECLARED
    # WHEN A VALUE CANNOT DISAMBIGUATE ITSELF. `generation` is 0 both for a
    # founder and for empty lattice sites, because `export_checkpoint.jl`
    # zero-fills the array and skips unoccupied voxels -- so no `background`
    # can be right: 0 deletes the founders, None draws the void as biomass.
    # Naming the layer that says WHERE BIOMASS IS resolves it, and only the
    # producer knows which layer that is.
    occupancy_from: str | None = None
    note: str = ""

    def as_dict(self) -> dict:
        return {"name": self.name, "grid_id": self.grid_id, "unit": self.unit,
                "semantic_kind": self.semantic_kind,
                "shape_xyz": [int(n) for n in self.data.shape],
                "dtype": str(self.data.dtype),
                "native": bool(self.native),
                "authoritative_for_quantitation":
                    bool(self.authoritative_for_quantitation),
                "source_grid_id": self.source_grid_id,
                "derivation": self.derivation,
                "background": (self.background
                               if self.background is None
                               or is_undeclared(self.background)
                               else float(self.background)),
                "occupancy_from": self.occupancy_from,
                "note": self.note}


@dataclass(frozen=True)
class Table:
    """A columnar product, for anything that is NOT a field.

    Per-cell and per-lineage dose attribution belongs here and not on a grid.
    It is a join of a dose field with labels, and when that field was upsampled
    first, the result has one value per coarse bin repeated across every label
    in it. Drawn as a voxel image it reads as per-parcel detail that was never
    computed; as a table with `basis` stamped on it, it reads as what it is.
    """
    name: str
    columns: dict
    basis: str
    units: dict = field(default_factory=dict)
    note: str = ""

    def as_dict(self) -> dict:
        return {"name": self.name, "basis": self.basis,
                "columns": {k: len(np.asarray(v)) for k, v in self.columns.items()},
                "units": dict(self.units), "note": self.note}


def _has_unknown_sentinel(data) -> bool:
    """Any negative value that is not the schema's out-of-domain marker."""
    arr = np.asarray(data)
    if not np.issubdtype(arr.dtype, np.number):
        return False
    negative = arr[arr < 0]
    return bool(negative.size and (negative != OUT_OF_DOMAIN).any())


def bundle_problems(grids, layers, tables=()) -> list[str]:
    """Everything wrong with a bundle, collected — never the first thing wrong.

    The same reporting style as the config loader: a caller assembling a bundle
    wants the whole list, not to fix one field per run.
    """
    out: list[str] = []
    by_id = {}
    for grid in grids:
        if grid.grid_id in by_id:
            out.append(f"duplicate grid_id {grid.grid_id!r}")
        by_id[grid.grid_id] = grid
        if len(grid.shape_xyz) != 3 or any(int(n) < 1 for n in grid.shape_xyz):
            out.append(f"grid {grid.grid_id!r} shape_xyz must be three positive "
                       f"integers, got {list(grid.shape_xyz)}")
        if len(grid.spacing_cm) != 3 or any(float(v) <= 0 for v in grid.spacing_cm):
            out.append(f"grid {grid.grid_id!r} spacing_cm must be three positive "
                       f"lengths, got {list(grid.spacing_cm)}")
        if len(grid.origin_cm) != 3:
            out.append(f"grid {grid.grid_id!r} origin_cm must be three coordinates")

    for grid in grids:
        ref = grid.material_resolution_grid
        if ref is not None and ref not in by_id:
            out.append(f"grid {grid.grid_id!r} names material_resolution_grid "
                       f"{ref!r}, which is not declared in this bundle")

    # ONE COORDINATE SYSTEM IS THE PREMISE. Grids may sample the volume
    # differently, but a bundle whose grids describe different volumes is not a
    # bundle — it is two datasets in one file, and any overlay drawn from it is
    # a lie about where things are.
    if len(by_id) > 1:
        reference = grids[0]
        for grid in grids[1:]:
            for axis, (a, b) in enumerate(zip(grid.extent_cm, reference.extent_cm)):
                if abs(a - b) > 1e-9 * max(1.0, abs(b)):
                    out.append(
                        f"grid {grid.grid_id!r} spans {a:.6g} cm on axis {axis} "
                        f"but {reference.grid_id!r} spans {b:.6g} cm — grids in "
                        "one bundle must describe the same physical volume")
            if tuple(float(v) for v in grid.origin_cm) != \
                    tuple(float(v) for v in reference.origin_cm):
                out.append(f"grid {grid.grid_id!r} origin differs from "
                           f"{reference.grid_id!r}; one coordinate system means "
                           "one origin")

    seen: set = set()
    for layer in layers:
        if layer.name in seen:
            out.append(f"duplicate layer name {layer.name!r}")
        seen.add(layer.name)
        grid = by_id.get(layer.grid_id)
        if grid is None:
            out.append(f"layer {layer.name!r} names grid {layer.grid_id!r}, "
                       "which is not declared in this bundle")
        elif tuple(int(n) for n in layer.data.shape) != \
                tuple(int(n) for n in grid.shape_xyz):
            out.append(f"layer {layer.name!r} has shape {list(layer.data.shape)} "
                       f"but grid {layer.grid_id!r} is {list(grid.shape_xyz)}")
        if layer.semantic_kind not in SEMANTIC_KINDS:
            out.append(f"layer {layer.name!r} semantic_kind "
                       f"{layer.semantic_kind!r} is not one of "
                       f"{sorted(SEMANTIC_KINDS)}")
        if not layer.unit:
            out.append(f"layer {layer.name!r} has no unit; a blank unit is never "
                       "read as dimensionless, declare 'dimensionless'")

        # THE RULE THE FORMAT EXISTS FOR: upsampling invents information, so an
        # upsampled layer may be drawn and may not be quoted. Reduction destroys
        # information and invents none, so it may be both.
        if layer.derivation in UPSAMPLING_DERIVATIONS and \
                layer.authoritative_for_quantitation:
            out.append(
                f"layer {layer.name!r} is {layer.derivation} and claims "
                "authoritative_for_quantitation; upsampling broadcasts one "
                "value across bins that were never separately measured, so it "
                "may be drawn and may not be quoted")
        # A reduction must actually reduce. A layer labelled `block_sum` onto a
        # FINER grid is mislabelled upsampling, and the label is the only thing
        # standing between it and being quoted.
        if layer.derivation in REDUCING_DERIVATIONS and grid is not None \
                and layer.source_grid_id in by_id:
            source = by_id[layer.source_grid_id]
            if int(np.prod(grid.shape_xyz)) >= int(np.prod(source.shape_xyz)):
                out.append(
                    f"layer {layer.name!r} claims {layer.derivation} from "
                    f"{layer.source_grid_id!r} onto a grid that is not coarser; "
                    "a reduction that does not reduce is upsampling wearing the "
                    "wrong label")
        if not layer.native and layer.derivation is None:
            out.append(f"layer {layer.name!r} is derived but does not say how")
        if layer.native and layer.derivation is not None:
            out.append(f"layer {layer.name!r} is native but carries derivation "
                       f"{layer.derivation!r}; native means measured on its own "
                       "grid, so the two cannot both be true")
        if layer.derivation is not None and layer.derivation not in DERIVATIONS:
            out.append(f"layer {layer.name!r} derivation {layer.derivation!r} "
                       f"is not one of {sorted(DERIVATIONS)}")
        if layer.source_grid_id is not None and layer.source_grid_id not in by_id:
            out.append(f"layer {layer.name!r} names source_grid_id "
                       f"{layer.source_grid_id!r}, which is not declared")
        if not layer.native and layer.source_grid_id is None:
            out.append(f"layer {layer.name!r} is derived but does not say from "
                       "which grid")

        # A DANGLING OCCUPANCY REFERENCE FAILS OPEN. `occupancy_from` exists
        # because `generation` cannot state its own absence -- 0 is a founder
        # and 0 is empty space -- so a misspelt or cross-grid reference does
        # not merely mislabel the layer, it silently returns the renderer to
        # the ambiguity the field was added to resolve. Same treatment as
        # source_grid_id above: name it at write time, not at draw time.
        # OMISSION MUST NOT ACQUIRE A SEMANTICS. `background=None` means "every
        # cell carries information, hide nothing" -- a real and sometimes
        # correct statement, but a terrible DEFAULT, because a producer that
        # simply forgot silently claims it. Two layers did: `omega_b` drew its
        # false cells as part of the region, and the upsampled `cell_id_on_r*`
        # overlay grew a shell of empty space. Both were found by review, not
        # here. A categorical or boolean layer must now SAY which it means.
        if layer.semantic_kind in (CATEGORICAL, BOOLEAN) \
                and is_undeclared(layer.background) \
                and layer.occupancy_from is None:
            out.append(
                f"layer {layer.name!r} is {layer.semantic_kind} but declares "
                "neither `background` nor `occupancy_from`, so a renderer "
                "cannot tell absence from data. Say which value means "
                "'nothing here', name the layer that does, or pass "
                "`background=None` to state that every cell is informative")

        if layer.occupancy_from is not None:
            other = next((l for l in layers if l.name == layer.occupancy_from),
                         None)
            if other is None:
                out.append(f"layer {layer.name!r} names occupancy_from "
                           f"{layer.occupancy_from!r}, which is not a layer in "
                           "this bundle")
            elif other.grid_id != layer.grid_id:
                out.append(f"layer {layer.name!r} takes occupancy from "
                           f"{layer.occupancy_from!r}, which is on grid "
                           f"{other.grid_id!r} and not {layer.grid_id!r}; a "
                           "mask must be cell-for-cell with what it masks")
            elif layer.background is not None \
                    and not is_undeclared(layer.background):
                out.append(f"layer {layer.name!r} declares both a background "
                           "value and an occupancy layer; they are two answers "
                           "to one question and the renderer would pick one")
            elif other.background == 0 and _has_unknown_sentinel(other.data):
                # THE OTHER HALF OF THE CONTRACT. Declaring background 0 says
                # which value is empty; it does NOT say that every negative is
                # a sentinel, and `Layer.background` cannot -- it names one
                # value. But the mask tests `> 0`, so a source using, say, -7
                # as legitimate data loses those cells silently. Only the
                # schema's own -1 is a known sentinel; anything else is an
                # encoding this mask does not implement.
                out.append(
                    f"layer {layer.name!r} takes occupancy from "
                    f"{layer.occupancy_from!r}, which declares background 0 but "
                    f"holds negative values other than {OUT_OF_DOMAIN} "
                    "(out-of-domain). The occupancy mask tests `> 0`, so those "
                    "cells would be dropped as if they were empty")
            elif other.background != 0:
                # THE MASK IMPLEMENTS ONE ENCODING, so only that one is
                # accepted. `observer.occupied_mask` tests `> 0`: background 0
                # is empty and any NEGATIVE value is a sentinel outside the
                # biological domain, which is the exchange schema's convention
                # (cell_id 0 background, -1 wall). A source declaring
                # `background = -1` would mean label 0 is valid data, and the
                # mask would silently delete every zero-valued cell.
                #
                # Refusing is honest; quietly rendering the wrong cells is not.
                out.append(
                    f"layer {layer.name!r} takes occupancy from "
                    f"{layer.occupancy_from!r}, which declares background "
                    f"{other.background!r}. The occupancy mask tests `> 0`, so "
                    "only a source with background 0 (and negative values as "
                    "out-of-domain sentinels) is supported")

        # A label is not a number. Averaging two cell ids gives a third cell id,
        # which is why categorical layers may never be resampled arithmetically
        # and why a float dtype here is a defect rather than a style choice.
        if layer.semantic_kind == CATEGORICAL and \
                not np.issubdtype(layer.data.dtype, np.integer):
            out.append(f"layer {layer.name!r} is categorical but has dtype "
                       f"{layer.data.dtype}; labels must be integers")
        if layer.semantic_kind == BOOLEAN and layer.data.dtype != np.bool_:
            out.append(f"layer {layer.name!r} is boolean but has dtype "
                       f"{layer.data.dtype}")
        if layer.semantic_kind == CATEGORICAL and \
                layer.derivation in {"block_sum", "mass_weighted_mean"}:
            out.append(f"layer {layer.name!r} is categorical and was resampled "
                       f"by {layer.derivation!r}; arithmetic on labels produces "
                       "a label that was never in the data")

    for table in tables:
        if not table.basis:
            out.append(f"table {table.name!r} has no declared basis")
        lengths = {k: len(np.asarray(v)) for k, v in table.columns.items()}
        if len(set(lengths.values())) > 1:
            out.append(f"table {table.name!r} columns have differing lengths: "
                       f"{lengths}")
    return out


def manifest(grids, layers, tables=(), provenance=None) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "logical_axis_order": "xyz",
        "dataset_axis_order_h5py": "zyx",
        "grids": [g.as_dict() for g in grids],
        "layers": [l.as_dict() for l in layers],
        "tables": [t.as_dict() for t in tables],
        "provenance": dict(provenance or {}),
    }


def write_bundle(path, grids, layers, tables=(), provenance=None) -> dict:
    """Write the bundle, or refuse with every problem at once.

    Returns the manifest so a caller can record or assert on it without
    reopening the file.
    """
    import h5py

    problems = bundle_problems(grids, layers, tables)
    if problems:
        raise ValueError("refusing to write a viewer bundle:\n  - "
                         + "\n  - ".join(problems))
    doc = manifest(grids, layers, tables, provenance)
    with h5py.File(path, "w") as f:
        f.attrs["schema_version"] = SCHEMA_VERSION
        f.attrs["logical_axis_order"] = "xyz"
        f.attrs["dataset_axis_order_h5py"] = "zyx"
        # The manifest travels as one JSON blob AND as per-dataset attrs. The
        # blob is what a reader parses; the attrs are what survives someone
        # opening the file in a generic HDF5 browser and looking at one array.
        f.attrs["manifest_json"] = json.dumps(doc, sort_keys=True)
        for key, value in (provenance or {}).items():
            f.attrs[key] = "" if value is None else value
        for grid in grids:
            group = f.create_group(f"grids/{grid.grid_id}")
            for key, value in grid.as_dict().items():
                group.attrs[key] = "" if value is None else value
        for layer in layers:
            dset = f.create_dataset(f"layers/{layer.name}",
                                    data=_logical_to_h5py(layer.data))
            for key, value in layer.as_dict().items():
                dset.attrs[key] = "" if value is None else value
        for table in tables:
            group = f.create_group(f"tables/{table.name}")
            group.attrs["basis"] = table.basis
            group.attrs["note"] = table.note
            for column, values in table.columns.items():
                group[column] = np.asarray(values)
                if column in table.units:
                    group[column].attrs["unit"] = table.units[column]
    return doc


def read_manifest(path) -> dict:
    import h5py

    with h5py.File(path, "r") as f:
        return json.loads(f.attrs["manifest_json"])


def read_layer(path, name) -> tuple:
    """The array in logical (x,y,z) order, and its declared attributes.

    Returned together on purpose: an array whose `authoritative_for_quantitation`
    flag was left behind in the file is an array a caller will quote.
    """
    import h5py

    with h5py.File(path, "r") as f:
        dset = f[f"layers/{name}"]
        return _h5py_to_logical(dset[()]), dict(dset.attrs)


def quantitative_layers(manifest_doc) -> list[str]:
    """The layers a reader may compute with. Everything else is for looking at.

    KEYED ON `authoritative_for_quantitation` ALONE, because that is the field
    named for the decision and `native` is descriptive. Requiring both was a
    leftover from the blunter first rule, and it silently demoted the case the
    sharpened rule exists to permit: an exact mass-weighted reduction is derived
    (native = false) and IS quotable, since it reproduces the native field at
    that resolution to 5.6e-16. `bundle_problems` already guarantees the flag
    cannot be set on an upsampling, so the flag alone is sufficient.
    """
    return [l["name"] for l in manifest_doc["layers"]
            if l["authoritative_for_quantitation"]]


MELANIZED_FRAC_KEY = "melanized_fraction"


@dataclass(frozen=True)
class LayerSpec:
    name: str
    units: str
    dtype: str = "float64"
    quotable: bool = True
    description: str = ""
    cmap: str = "viridis"
    clim: tuple[float, float] = (0.0, 1.0)
    categorical: bool = False
    derived: bool = False
    native_openmc_field: bool = False
    grid_representation: str = "resampled_extruded_display_grid"
    synthetic_surrogate: bool = True
    measured_physical_dose: bool = False


LAYER_CONTRACTS: dict[str, LayerSpec] = {
    "cell_id": LayerSpec(
        "cell_id", "id", "int32", True,
        "Discrete CPM spin ID (0=water, -1=void/steam)",
        cmap="tab20", clim=(-1.0, 64.0), categorical=True, derived=False,
    ),
    "species_id": LayerSpec(
        "species_id", "index", "int32", True,
        "Species index 1..5 (0=water, -1=void)",
        cmap="tab10", clim=(1.0, 7.0), categorical=True, derived=False,
    ),
    "local_attenuation_kappa": LayerSpec(
        "local_attenuation_kappa", "1/lattice", "float64", True,
        "State-dependent macro attenuation coefficient kappa(sigma)",
        cmap="viridis", clim=(0.0, 6.0), categorical=False, derived=True,
    ),
    "fast_flux": LayerSpec(
        "fast_flux", "rel_flux", "float64", True,
        "Primary unmoderated fast radiation group I_fast (synthetic surrogate, not measured Gy)",
        cmap="inferno", clim=(0.0, 1.0), categorical=False, derived=False,
    ),
    "thermal_flux": LayerSpec(
        "thermal_flux", "rel_flux", "float64", True,
        "Moderated secondary radiolytic group I_therm (synthetic surrogate, not measured Gy)",
        cmap="inferno", clim=(0.0, 1.0), categorical=False, derived=False,
    ),
    "total_radiation": LayerSpec(
        "total_radiation", "rel_flux", "float64", True,
        "Total effective ionizing field I_fast + I_therm (synthetic surrogate, not measured Gy)",
        cmap="inferno", clim=(0.0, 1.0), categorical=False, derived=True,
    ),
    "precursor_damage_M": LayerSpec(
        "precursor_damage_M", "conc_au", "float64", True,
        "Accumulated precursor damage (Iodine-135 analogue)",
        cmap="magma", clim=(0.0, 1.5), categorical=False, derived=False,
    ),
    "xenon_inhibitor_C": LayerSpec(
        "xenon_inhibitor_C", "conc_au", "float64", True,
        "Delayed radiolytic inhibitor (Xenon-135 poison analogue)",
        cmap="magma", clim=(0.0, 1.5), categorical=False, derived=True,
    ),
    "bilinear_burnout_sink": LayerSpec(
        "bilinear_burnout_sink", "conc_au/s", "float64", True,
        "Active flux-scavenging burnout rate sigma_burn * I * rho_mel * C",
        cmap="plasma", clim=(0.0, 0.8), categorical=False, derived=True,
    ),
    "void_mask": LayerSpec(
        "void_mask", "bool", "uint8", False,
        "Binary mask of low-attenuation radiolytic/steam voids",
        cmap="gray", clim=(0.0, 1.0), categorical=True, derived=True,
    ),
    "melanin_shield_mask": LayerSpec(
        "melanin_shield_mask", "bool", "uint8", False,
        "Binary mask of high-attenuation melanized cells (Sp 4 & 5)",
        cmap="inferno", clim=(0.0, 1.0), categorical=True, derived=True,
    ),
    MELANIZED_FRAC_KEY: LayerSpec(
        MELANIZED_FRAC_KEY, "fraction", "float64", True,
        "Biomass fraction of melanized species (Sp 4 + Sp 5)",
        cmap="viridis", clim=(0.0, 1.0), categorical=False, derived=True,
    ),
    "shannon_diversity_H": LayerSpec(
        "shannon_diversity_H", "nats", "float64", True,
        "Community Shannon diversity index H",
        cmap="viridis", clim=(0.0, 2.0), categorical=False, derived=True,
    ),
    "void_reactivity_coeff": LayerSpec(
        "void_reactivity_coeff", "dI/df_void", "float64", True,
        "Effective void amplification diagnostic (positive void coefficient)",
        cmap="coolwarm", clim=(-1.0, 1.0), categorical=False, derived=True,
    ),
}


def compute_and_validate_layers(raw_state: dict, params: dict) -> dict:
    """Compute all derived spatial and scalar layers in viewer.py and verify
    that every LAYER_CONTRACTS entry is present and finite.
    Supports both the two-group state dictionary ('cell_id', 'species_id',
    'fast_flux', 'thermal_flux', 'precursor_damage_M', 'xenon_inhibitor_C')
    and the continuum-melanin snapshot dictionary ('spin', 'species',
    'melanin', 'radiation', 'poison').

    Scientific & Software Contract:
    - Never mutates `raw_state` or `params`.
    - All radiation/flux layers (`fast_flux`, `thermal_flux`, `total_radiation`)
      are dimensionless synthetic surrogates (`rel_flux`), never measured
      physical dose (`Gy` or `Gy/s`), and are defined on the CPM display grid
      rather than a native OpenMC transport tally mesh.
    """
    continuum_mode = "spin" in raw_state and "cell_id" not in raw_state

    if continuum_mode:
        sigma = np.array(raw_state["spin"], dtype=np.int32, copy=True)
        occupied = sigma > 0
        raw_sp = np.asarray(raw_state["species"], dtype=np.int32)
        if raw_sp.ndim == 1:
            species = np.zeros_like(sigma, dtype=np.int32)
            species[occupied] = raw_sp[sigma[occupied] - 1]
        else:
            species = np.array(raw_sp, dtype=np.int32, copy=True)
        melanin = np.array(raw_state["melanin"], dtype=np.float64, copy=True)
        I_tot = np.array(raw_state["radiation"], dtype=np.float64, copy=True)
        I_fast = I_tot.copy()
        I_therm = np.zeros_like(I_tot)
        C_inh = np.array(raw_state["poison"], dtype=np.float64, copy=True)
        M_dam = np.array(
            raw_state.get("poison_lag", melanin), dtype=np.float64, copy=True
        )
        if not (
            sigma.shape
            == species.shape
            == melanin.shape
            == I_tot.shape
            == C_inh.shape
        ):
            raise ValueError("Spatial input arrays must share identical dimensions")

        kappa_base = float(params.get("kappa", 1.5))
        kappa_void = float(params.get("kappa_void", 0.15))
        kappa_cell = float(params.get("kappa_cell", 1.0))
        kappa_mel = float(params.get("kappa_mel", 2.5))
        if kappa_base == 0.0:
            kappa_map = np.zeros_like(I_tot)
        else:
            kappa_map = np.where(
                occupied,
                kappa_base * (kappa_cell + kappa_mel * melanin),
                kappa_base * kappa_void,
            )
        sigma_burn = float(params.get("sigma_burn", 0.65))
        burnout_sink = sigma_burn * I_tot * (0.2 + melanin) * C_inh
        void_mask = (~occupied).astype(np.uint8)
        melanin_mask = ((species >= 4) | (melanin >= 0.4)).astype(np.uint8)
        kappa_water = kappa_base * kappa_cell
    else:
        sigma = np.array(raw_state["cell_id"], dtype=np.int32, copy=True)
        species = np.array(raw_state["species_id"], dtype=np.int32, copy=True)
        I_fast = np.array(raw_state["fast_flux"], dtype=np.float64, copy=True)
        I_therm = np.array(raw_state["thermal_flux"], dtype=np.float64, copy=True)
        M_dam = np.array(raw_state["precursor_damage_M"], dtype=np.float64, copy=True)
        C_inh = np.array(raw_state["xenon_inhibitor_C"], dtype=np.float64, copy=True)

        if not (
            sigma.shape
            == species.shape
            == I_fast.shape
            == I_therm.shape
            == M_dam.shape
            == C_inh.shape
        ):
            raise ValueError("Spatial input arrays must share identical dimensions")

        I_tot = I_fast + I_therm
        void_mask = (sigma == -1).astype(np.uint8)
        melanin_mask = np.isin(species, [4, 5]).astype(np.uint8)

        kappa_water = float(params.get("kappa_water", 0.35))
        kappa_void = float(params.get("kappa_void", 0.02))
        kappa_sp = list(params.get("kappa_species", [0.35, 0.45, 0.50, 1.40, 2.20]))

        kappa_map = np.full(sigma.shape, kappa_water, dtype=np.float64)
        kappa_map[sigma == -1] = kappa_void
        for sp_idx, k_val in enumerate(kappa_sp, start=1):
            kappa_map[species == sp_idx] = float(k_val)

        sigma_burn = float(params.get("sigma_burn", 0.45))
        burnout_sink = (
            sigma_burn * I_tot * melanin_mask.astype(np.float64) * C_inh
        )

    bio_mask = species > 0
    total_bio = int(np.sum(bio_mask))
    if total_bio > 0:
        mel_frac = float(np.sum(melanin_mask[bio_mask]) / total_bio)
        max_sp = max(5, int(np.max(species)))
        counts = np.array(
            [np.sum(species == s) for s in range(1, max_sp + 1)],
            dtype=np.float64,
        )
        probs = counts[counts > 0] / total_bio
        shannon_H = float(-np.sum(probs * np.log(probs)))
    else:
        mel_frac = 0.0
        shannon_H = 0.0

    mean_void = float(np.mean(void_mask)) if void_mask.size > 0 else 0.0
    mean_flux = float(np.mean(I_tot)) if I_tot.size > 0 else 0.0
    baseline_flux = float(params.get("I0", 1.0)) * float(
        np.exp(-0.5 * kappa_water)
    )
    void_reactivity = (mean_flux - baseline_flux) / max(mean_void, 1e-6)

    out = {
        "cell_id": sigma,
        "species_id": species,
        "local_attenuation_kappa": kappa_map,
        "fast_flux": I_fast,
        "thermal_flux": I_therm,
        "total_radiation": I_tot,
        "precursor_damage_M": M_dam,
        "xenon_inhibitor_C": C_inh,
        "bilinear_burnout_sink": burnout_sink,
        "void_mask": void_mask,
        "melanin_shield_mask": melanin_mask,
        MELANIZED_FRAC_KEY: mel_frac,
        "shannon_diversity_H": shannon_H,
        "void_reactivity_coeff": void_reactivity,
    }

    if continuum_mode:
        out.update(
            {
                "spin": sigma,
                "melanin": melanin,
                "radiation": I_tot,
                "local_attenuation": kappa_map,
                "poison": C_inh,
                "burnout_rate": burnout_sink,
            }
        )

    for key, spec in LAYER_CONTRACTS.items():
        if key not in out:
            raise KeyError(f"Contract violation: missing derived layer '{key}'")
        val = out[key]
        if spec.quotable and not np.all(np.isfinite(val)):
            raise ValueError(
                f"Contract violation: quotable layer '{key}' contains NaN/Inf"
            )

    return out
