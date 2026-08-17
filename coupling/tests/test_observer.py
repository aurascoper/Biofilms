"""The read-only observer: what it must say, and what it must not become.

Almost everything here runs in the BARE TIER, because almost nothing the
observer does is rendering. Which layers may be drawn, what each must be
labelled with, how a categorical field maps to colour, where a grid sits in
centimetres, what the provenance panel says — those are decisions about data,
and a runner with no renderer installed can check every one of them.

Only the two PyVista entry points are gated, and they are the thin part.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from biofilm_openmc import observer
from biofilm_openmc.observer import (SPECIES_COLOURS, SPECIES_LABELS,
                                     display_plan, grid_geometry,
                                     occupied_mask, provenance_panel,
                                     species_legend)
from biofilm_openmc.viewer import (UNDECLARED, Grid, Layer, manifest,
                                   write_bundle)

LATTICE = Grid("cpm_labels", (4, 4, 4), (0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
REFINED = Grid("dose_refinement_4", (16, 16, 16), (0.0, 0.0, 0.0),
               (0.0625, 0.0625, 0.0625), material_resolution_grid="cpm_labels")


def _bundle_manifest():
    labels = Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                   np.arange(64, dtype=np.int32).reshape((4, 4, 4)),
                   background=0)
    species = Layer("species_id", "cpm_labels", "dimensionless", "categorical",
                    np.full((4, 4, 4), 2, dtype=np.int32), background=0)
    omega = Layer("omega_b", "cpm_labels", "dimensionless", "boolean",
                  np.ones((4, 4, 4), dtype=bool), background=0)
    dose = Layer("dose_per_source_r4", "dose_refinement_4",
                 "Gy/source-particle", "intensive", np.zeros((16, 16, 16)))
    upsampled = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                      "intensive", np.zeros((4, 4, 4)),
                      native=False, authoritative_for_quantitation=False,
                      source_grid_id="dose_refinement_4",
                      derivation="upsampled_coarse_dose")
    return manifest([LATTICE, REFINED], [labels, species, omega, dose, upsampled],
                    provenance={"reference_system_id": "synthetic_biofilm_e2e",
                                "target_calibration": False,
                                "evidence_policy": "synthetic",
                                "openmc_version": "0.15.3"})


# ------------------------------------------------------- what it must not be

def test_the_observer_computes_nothing():
    """A display layer that can compute is a display layer that will eventually
    be quoted.

    Checked on the IMPORT GRAPH, not on the text: the module docstring names the
    modules it refuses to import, and a substring test would fail on the very
    sentence that documents the rule.
    """
    import ast
    from pathlib import Path

    tree = ast.parse(Path(observer.__file__).read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
            imported.update(f"{node.module or ''}.{a.name}" for a in node.names)
    for forbidden in ("feedback_uq", "synthetic_gate", "feedback_gate"):
        assert not any(forbidden in name for name in imported), (
            f"observer.py imports {forbidden}; rendering and deciding must not "
            "live in the same module")


def test_pyvista_is_imported_function_locally_only():
    """So the whole contract is exercised on a runner with no renderer."""
    from pathlib import Path
    src = (Path(observer.__file__)).read_text()
    for line in src.splitlines():
        if line.startswith("import pyvista") or line.startswith("from pyvista"):
            pytest.fail("pyvista is imported at module scope; the bare tier "
                        "must be able to import this module")
    assert "    import pyvista as pv" in src


def _unresolved_globals(fn):
    """Names `fn` loads as globals that resolve to nothing.

    LOAD_GLOBAL, not `co_names`. `co_names` mixes global lookups with attribute
    names -- `pv.Plotter` contributes "Plotter" -- so a check built on it has to
    guess which entries are real, and the first version of this test guessed by
    requiring the name be UPPERCASE *and* an actual export of viewer.py. That
    matched the one bug in front of it and nothing else: a misspelling like
    `CATGORICAL` is not a viewer export, `occupied_maks` is not uppercase, and
    both would have sailed through while `plot_layer` raised NameError on any
    host with a renderer. A check narrowed until it matches exactly one known
    defect has stopped being a check.

    The bytecode already draws the distinction, so ask it instead of inferring.
    """
    import builtins
    import dis

    return sorted({
        ins.argval for ins in dis.get_instructions(fn)
        if ins.opname == "LOAD_GLOBAL"
        and ins.argval not in vars(observer)
        and not hasattr(builtins, ins.argval)})


def test_every_global_the_render_path_uses_actually_resolves():
    """THE BARE TIER MUST STILL BE ABLE TO FAIL ON THE RENDER PATH.

    `plot_layer` shipped referencing CATEGORICAL and BOOLEAN, which live in
    viewer.py and were never imported here. Every call on a machine with a
    renderer raised NameError on the first line of the branch. Nothing caught
    it: pyvista is absent from CI and from the dev environment, so all eight
    render tests SKIPPED and the suite reported green over a dead function.

    "8 skipped" is not a neutral line in a test report -- it is the uncovered
    surface, and this is the check that covers the part of it needing no
    renderer. Name resolution is decidable statically, so decide it statically.
    """
    for name in dir(observer):
        fn = getattr(observer, name)
        if not callable(fn) or not hasattr(fn, "__code__"):
            continue
        if fn.__code__.co_filename != observer.__file__:
            continue
        assert not _unresolved_globals(fn), (
            f"observer.{name} loads {_unresolved_globals(fn)} as globals, and "
            "nothing defines them; this raises NameError wherever it runs")


def test_the_name_check_catches_a_misspelling():
    """THE CONTROL FOR THE CONTROL.

    The test above passes when it finds nothing, so it is worth exactly as much
    as its ability to find something. Its previous version could find only an
    uppercase name that viewer.py really exports -- it would have missed
    `CATGORICAL` and `occupied_maks`, which are the realistic ways this breaks.
    """
    src = ("def broken(x):\n"
           "    if x in (CATGORICAL, BOOLEAN):\n"
           "        return occupied_maks(x)\n"
           "    return np.zeros(3)\n")
    ns = {}
    exec(compile(src, "<broken>", "exec"), ns)

    found = _unresolved_globals(ns["broken"])
    assert "CATGORICAL" in found, "a misspelled constant must be caught"
    assert "occupied_maks" in found, "a misspelled lowercase helper too"
    # ...and the names that DO resolve in observer are not reported
    assert "BOOLEAN" not in found
    assert "np" not in found


# --------------------------------------------------------- what it must say

def test_a_non_quotable_layer_carries_a_banner_naming_the_mechanism():
    """The reader needs to know WHY the numbers are unreadable. 'Derived' does
    not distinguish a field that was never measured at this resolution from one
    that was resampled exactly."""
    plan = {l.name: l for l in display_plan(_bundle_manifest())}
    up = plan["dose_on_lattice"]
    assert not up.quotable
    assert "NOT QUANTITATIVE" in up.banner
    assert "upsampled" in up.banner
    assert "dose_refinement_4" in up.banner
    assert "repeated" in up.banner, (
        "the banner must say what is actually wrong: one coarse value repeated")


def test_every_quotable_layer_has_an_empty_banner_and_vice_versa():
    for layer in display_plan(_bundle_manifest()):
        assert bool(layer.banner) != layer.quotable


def test_a_derived_but_quotable_layer_gets_a_note_and_no_warning():
    """QUOTABILITY KEYS ON `authoritative_for_quantitation` ALONE. Requiring
    `native` too was a leftover from the blunter first rule and it demoted
    exactly the case the sharpened rule permits: an exact mass-weighted
    reduction is derived AND quotable, since it reproduces the native field at
    that resolution to 5.6e-16.

    Caught by running the observer against the real bundle, where it labelled
    `dose_on_lattice_from_r4` NOT QUANTITATIVE. Stamping that with the same
    alarming banner an upsampling gets would teach a reader to ignore banners.
    """
    reduced = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                    "intensive", np.zeros((4, 4, 4)),
                    native=False, authoritative_for_quantitation=True,
                    source_grid_id="dose_refinement_4",
                    derivation="mass_weighted_mean")
    doc = manifest([LATTICE, REFINED], [reduced])
    layer = display_plan(doc)[0]
    assert layer.quotable
    assert layer.banner == ""
    assert layer.derivation_note == "mass_weighted_mean from dose_refinement_4"

    from biofilm_openmc.viewer import quantitative_layers
    assert quantitative_layers(doc) == ["dose_on_lattice"]


def test_units_travel_with_the_layer():
    """A colour bar labelled Gy/s on a per-source field has invented an
    emission rate the synthetic system does not have."""
    plan = {l.name: l for l in display_plan(_bundle_manifest())}
    assert plan["dose_per_source_r4"].unit == "Gy/source-particle"
    assert plan["cell_id"].unit == "dimensionless"


def test_only_species_maps_onto_the_seven_colour_palette():
    """cell_id and lineage_id are unbounded labels; giving them seven fixed
    colours would assert a species identity they do not carry."""
    plan = {l.name: l for l in display_plan(_bundle_manifest())}
    assert plan["species_id"].colour_by == "species"
    assert plan["cell_id"].colour_by == "scalar"
    assert plan["omega_b"].colour_by == "mask"


def test_the_provenance_panel_leads_with_target_calibration():
    rows = provenance_panel(_bundle_manifest())
    assert rows[0][0] == "target_calibration"
    assert rows[0][1] == "False"
    keys = dict(rows)
    assert keys["reference_system_id"] == "synthetic_biofilm_e2e"
    assert "1 quantitative" in keys["layers"] or "quantitative" in keys["layers"]


# ------------------------------------------------------------ the geometry

def test_grid_geometry_returns_point_dimensions_not_cell_counts():
    """A renderer wants cells+1; the bundle stores cells. Confusing them shifts
    every field by half a voxel — visible only as a misregistration between the
    label grid and the 4x dose grid, which is the overlay this exists for."""
    g = grid_geometry(_bundle_manifest(), "dose_refinement_4")
    assert g["cell_shape_xyz"] == (16, 16, 16)
    assert g["point_shape_xyz"] == (17, 17, 17)
    assert g["spacing_cm"] == (0.0625, 0.0625, 0.0625)
    assert g["material_resolution_grid"] == "cpm_labels"


def test_both_grids_land_on_the_same_physical_cube():
    doc = _bundle_manifest()
    a, b = grid_geometry(doc, "cpm_labels"), grid_geometry(doc, "dose_refinement_4")
    assert a["origin_cm"] == b["origin_cm"]
    assert a["extent_cm"] == pytest.approx(b["extent_cm"])


def test_an_unknown_grid_is_refused():
    with pytest.raises(KeyError):
        grid_geometry(_bundle_manifest(), "nope")


# -------------------------------------------------------------- the labels

def test_occupancy_uses_greater_than_zero_because_minus_one_is_the_wall():
    """The schema uses -1 for OUTSIDE the biological domain and 0 for
    background, so `!= 0` would draw the wall as biomass."""
    cell_id = np.array([[[0, 1]], [[-1, 3]]], dtype=np.int32)
    assert np.array_equal(occupied_mask(cell_id),
                          np.array([[[False, True]], [[False, True]]]))


def test_the_legend_lists_only_the_species_present():
    """Seven organisms named for a snapshot containing two is a legend that
    misdescribes the picture."""
    legend = species_legend(np.array([2, 2, 6]))
    assert [i for i, _, _ in legend] == [2, 6]
    assert [lab for _, lab, _ in legend] == ["D. radiodurans", "S. oneidensis"]
    assert len(species_legend()) == 7


def test_out_of_range_species_ids_are_dropped_not_coloured():
    assert species_legend(np.array([0, -1, 9])) == []


def test_the_palette_matches_the_serial_figure_convention():
    """Carried across from the CairoMakie prototype so successive viewers of
    this model look like one another."""
    assert SPECIES_COLOURS[0] == "#e6194b"
    assert len(SPECIES_COLOURS) == len(SPECIES_LABELS) == 7
    assert SPECIES_LABELS[1] == "D. radiodurans"


# ------------------------------------------------------ the rendering tier
#
# GATED PER TEST, NOT PER MODULE. A module-scope `importorskip` here would skip
# every test above it too -- 15 bare-tier checks silently reported as SKIPPED
# because a renderer is absent, which is precisely the accounting failure the
# project's "no -m deselection" rule exists to prevent.

_needs_pyvista = pytest.mark.skipif(
    importlib.util.find_spec("pyvista") is None,
    reason="viewer extra not installed")


@pytest.mark.viewer
@_needs_pyvista
def test_a_layer_becomes_cell_data_positioned_in_centimetres(tmp_path):
    """Cell data, never point data: a voxel field is piecewise constant over its
    cell, and interpolating it to points would smooth the material boundaries
    the model asserts are sharp — and make a categorical field meaningless."""
    path = tmp_path / "b.h5"
    doc = _bundle_manifest()
    write_bundle(path, [LATTICE, REFINED],
                 [Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                        np.arange(64, dtype=np.int32).reshape((4, 4, 4)),
                        background=0)])
    image = observer.to_image_data(path, "cell_id")
    assert image.dimensions == (5, 5, 5)
    assert "cell_id" in image.cell_data
    assert "cell_id" not in image.point_data
    assert image.spacing == pytest.approx((0.25, 0.25, 0.25))


@pytest.mark.viewer
@_needs_pyvista
def test_suppressing_a_banner_is_refused_rather_than_warned(tmp_path):
    """The banner is the only thing on screen distinguishing a real 4x dose
    field from a coarse field broadcast to look like one."""
    path = tmp_path / "b.h5"
    write_bundle(path, [LATTICE, REFINED], [
        Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle", "intensive",
              np.zeros((4, 4, 4)), native=False,
              authoritative_for_quantitation=False,
              source_grid_id="dose_refinement_4",
              derivation="upsampled_coarse_dose")])
    with pytest.raises(ValueError, match="banner may not be suppressed"):
        observer.plot_layer(path, "dose_on_lattice", show_banner=False)


def test_absence_semantics_survive_the_round_trip_to_the_display_plan():
    """A VALUE DOES NOT ALWAYS KNOW ITS OWN MEANING, so the producer says.

    `generation` is 0 for a founder and 0 for empty lattice sites, because
    `export_checkpoint.jl` zero-fills the array and skips unoccupied voxels. No
    background VALUE can separate those: 0 deletes the founding cohort from
    every picture of it, and None draws the void as generation-0 biomass. Both
    of those shipped, in that order — the second was the fix for the first.

    The disambiguator is a different layer, so the producer names it. This
    checks the declaration actually reaches the renderer; a field that is
    written to the manifest and dropped on the way back is worse than no field,
    because the producer has been told it was heard.
    """
    doc = manifest(
        [Grid("cpm_labels", (2, 2, 2), (0.0, 0.0, 0.0), (1e-3, 1e-3, 1e-3))],
        [Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
               np.zeros((2, 2, 2), np.int32), background=0),
         Layer("generation", "cpm_labels", "dimensionless", "categorical",
               np.zeros((2, 2, 2), np.int32), occupancy_from="cell_id")],
        [], provenance={})
    plan = {l.name: l for l in display_plan(doc)}

    assert plan["cell_id"].background == 0
    assert plan["cell_id"].occupancy_from is None
    # The ambiguous layer must NOT claim a background value...
    assert plan["generation"].background is UNDECLARED
    # ...and must not be left claiming every voxel is informative either.
    assert plan["generation"].occupancy_from == "cell_id"


@pytest.mark.parametrize("kwargs,expect", [
    ({"occupancy_from": "cel_id"},              "not a layer in this bundle"),
    ({"occupancy_from": "dose_r4"},             "and not"),
    ({"occupancy_from": "cell_id",
      "background": 0},                         "two answers to one question"),
])
def test_a_dangling_occupancy_reference_is_refused_at_write_time(kwargs, expect):
    """A BROKEN REFERENCE HERE FAILS OPEN.

    `occupancy_from` exists because `generation` cannot state its own absence —
    0 is a founder and 0 is empty space. So a misspelt name, or one pointing at
    a layer on a different grid, does not merely mislabel the layer: it drops
    the renderer straight back into the ambiguity the field was added to
    resolve, and it does it silently.

    Same treatment `source_grid_id` already gets — named when the bundle is
    written, not when someone tries to draw it.
    """
    from biofilm_openmc.viewer import bundle_problems

    grids = [LATTICE, REFINED]
    layers = [
        Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
              np.zeros((4, 4, 4), np.int32), background=0),
        Layer("dose_r4", "dose_refinement_4", "Gy/source-particle", "intensive",
              np.zeros((16, 16, 16))),
        Layer("generation", "cpm_labels", "dimensionless", "categorical",
              np.zeros((4, 4, 4), np.int32), **kwargs),
    ]
    found = bundle_problems(grids, layers)
    assert any(expect in p for p in found), found

    # and the correct declaration raises nothing, so this is not a test that
    # every occupancy reference is rejected
    layers[-1] = Layer("generation", "cpm_labels", "dimensionless",
                       "categorical", np.zeros((4, 4, 4), np.int32),
                       occupancy_from="cell_id")
    assert not [p for p in bundle_problems(grids, layers)
                if "occupancy" in p], bundle_problems(grids, layers)


def test_a_species_keeps_its_colour_whatever_subset_is_present():
    """THE INVARIANT THE PALETTE EXISTS FOR, broken by the palette code.

    `species_legend` correctly returns only the species present, each with its
    fixed colour. Passing those colours straight to the renderer with limits
    spanning the id range then distributes them ACROSS that range: for ids
    {1, 2, 7} three colours are spread over six units, so species 2 draws in
    whatever falls at that fraction while the legend labels it with the second
    palette entry. The picture disagrees with its own key.

    So the check is not "the right colours appear" but "the same organism lands
    in the same slot under different subsets", which is the actual claim.
    """
    from biofilm_openmc.observer import species_palette

    def slot_of(species_id, present):
        cmap, (lo, hi) = species_palette(species_legend(np.array(present)))
        # where the renderer places this id: which band of a uniform colormap
        # the value falls into
        index = int((species_id - lo) / (hi - lo) * len(cmap))
        return cmap[min(index, len(cmap) - 1)]

    for sid in (1, 2, 7):
        alone = slot_of(sid, [sid])
        assert alone == SPECIES_COLOURS[sid - 1], sid

    # the uneven subset from the finding
    for sid in (1, 2, 7):
        assert slot_of(sid, [1, 2, 7]) == SPECIES_COLOURS[sid - 1], (
            f"species {sid} changes colour when {{1,2,7}} are present")

    # and a contiguous subset that does not start at 1
    for sid in (3, 4, 5):
        assert slot_of(sid, [3, 4, 5]) == SPECIES_COLOURS[sid - 1], sid


def test_an_empty_species_layer_still_produces_a_palette():
    """AN EMPTY SNAPSHOT IS A LEGITIMATE FIELD, not a broken one.

    A species layer of pure background has no ids present, so `species_legend`
    returns nothing — which the existing out-of-range test already permits —
    and `min()` over nothing raised ValueError before the layer could draw.
    A viewer that crashes on an empty snapshot cannot show that the snapshot is
    empty, which is a thing a reader needs to see.
    """
    from biofilm_openmc.observer import species_palette

    assert species_legend(np.zeros((2, 2, 2), np.int32)) == []
    cmap, (lo, hi) = species_palette([])
    assert len(cmap) == len(SPECIES_COLOURS)
    assert lo < 1 and hi > len(SPECIES_COLOURS)


@_needs_pyvista
def test_an_all_background_species_layer_draws_instead_of_crashing(tmp_path):
    """"NOTHING IS HERE" IS A RESULT A READER NEEDS TO SEE.

    Making `species_palette` empty-safe only moved the crash one line down:
    thresholding an all-background layer leaves a dataset with zero cells, and
    `add_mesh` rejects it. The bare-tier control exercised `species_palette([])`
    and so could not detect that — the fix was verified one level below where
    it failed.

    This drives `plot_layer` itself, which is the only place the two interact.
    """
    path = tmp_path / "empty.h5"
    write_bundle(path, [LATTICE],
                 [Layer("species_id", "cpm_labels", "dimensionless",
                        "categorical", np.zeros((4, 4, 4), np.int32),
                        background=0)],
                 [], provenance={"reference_system_id": "synthetic",
                                 "target_calibration": False,
                                 "evidence_policy": "synthetic",
                                 "openmc_version": "0.15.3"})
    plotter = observer.plot_layer(path, "species_id")
    assert plotter is not None


def test_a_categorical_layer_must_declare_what_absence_means():
    """OMISSION MUST NOT ACQUIRE A SEMANTICS.

    `background=None` says "every cell carries information, hide nothing". That
    is true for `generation` and wrong for `omega_b`, whose false cells are not
    part of the region. When None was also the DEFAULT, a producer who simply
    forgot made that claim silently — and two did: `omega_b` drew its false
    cells as region, and the upsampled `cell_id_on_r*` overlay grew a shell of
    empty space around itself. Both were found by review, not here.

    So the default is UNDECLARED and is refused at write time. Saying "every
    cell is informative" is still available; it just has to be said.
    """
    from biofilm_openmc.viewer import bundle_problems

    def problems(**kwargs):
        return bundle_problems([LATTICE], [Layer(
            "cell_id", "cpm_labels", "dimensionless", "categorical",
            np.zeros((4, 4, 4), np.int32), **kwargs)])

    assert any("declares neither" in p for p in problems()), problems()
    # each of the three ways to declare is accepted
    assert not [p for p in problems(background=0) if "declares neither" in p]
    assert not [p for p in problems(background=None) if "declares neither" in p]
    assert not [p for p in problems(occupancy_from="cell_id")
                if "declares neither" in p]


def test_the_real_refinement_bundle_declares_absence_on_every_label_layer():
    """COVER THE PRODUCED MANIFEST, not a fixture that resembles it.

    The two layers that omitted a declaration — `omega_b` and the upsampled
    `cell_id_on_r*` — are built by `subvoxel_refinement.py`, and every test here
    used its own hand-written layers, so nothing looked at what the producer
    actually emits.
    """
    import ast
    from pathlib import Path

    producer = (Path(observer.__file__).resolve().parents[1] / "scripts"
                / "subvoxel_refinement.py")
    assert producer.exists(), f"the producer moved: {producer}"
    src = producer.read_text()
    tree = ast.parse(src)
    undeclared = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", "") == "Layer"):
            continue
        kinds = [a.value for a in node.args if isinstance(a, ast.Constant)]
        if not ({"categorical", "boolean"} & set(kinds)):
            continue
        declared = {k.arg for k in node.keywords}
        if not (declared & {"background", "occupancy_from"}):
            undeclared.append(kinds[0] if kinds else "?")
    assert not undeclared, (
        f"{len(undeclared)} categorical/boolean layers in subvoxel_refinement.py "
        "declare neither background nor occupancy_from")
