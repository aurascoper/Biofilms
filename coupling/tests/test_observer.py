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
from biofilm_openmc.viewer import Grid, Layer, manifest, write_bundle

LATTICE = Grid("cpm_labels", (4, 4, 4), (0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
REFINED = Grid("dose_refinement_4", (16, 16, 16), (0.0, 0.0, 0.0),
               (0.0625, 0.0625, 0.0625), material_resolution_grid="cpm_labels")


def _bundle_manifest():
    labels = Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                   np.arange(64, dtype=np.int32).reshape((4, 4, 4)))
    species = Layer("species_id", "cpm_labels", "dimensionless", "categorical",
                    np.full((4, 4, 4), 2, dtype=np.int32))
    omega = Layer("omega_b", "cpm_labels", "dimensionless", "boolean",
                  np.ones((4, 4, 4), dtype=bool))
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


def test_every_global_the_render_path_uses_actually_resolves():
    """THE BARE TIER MUST STILL BE ABLE TO FAIL ON THE RENDER PATH.

    `plot_layer` shipped referencing CATEGORICAL and BOOLEAN, which live in
    viewer.py and were never imported here. Every call on a machine with a
    renderer raised NameError on the first line of the branch. Nothing caught
    it: pyvista is absent from CI and from the dev environment, so all eight
    render tests SKIPPED and the suite reported green over a dead function.

    "8 skipped" is not a neutral line in a test report — it is the uncovered
    surface, and this is the check that covers the part of it that needs no
    renderer. Name resolution is decidable statically, so it is checked
    statically.
    """
    import builtins

    for fn in (observer.plot_layer, observer.to_image_data,
               observer.species_legend, observer.display_plan,
               observer.grid_geometry, observer.provenance_panel):
        unresolved = [name for name in fn.__code__.co_names
                      if name not in vars(observer)
                      and not hasattr(builtins, name)
                      # attribute access shares co_names with global lookups;
                      # `pv.Plotter` puts "Plotter" here. Locals are the import
                      # site for those, so anything bound locally is fine.
                      and name not in fn.__code__.co_varnames]
        # Attribute names on objects are indistinguishable from globals in
        # co_names, so this cannot be a bare emptiness assert. What it CAN do is
        # require that no unresolved name is one the module also uses as a
        # module-level constant elsewhere -- which is exactly the CATEGORICAL
        # case. Cross-check against the names viewer.py exports.
        from biofilm_openmc import viewer
        leaked = [n for n in unresolved
                  if n.isupper() and hasattr(viewer, n)]
        assert not leaked, (
            f"{fn.__name__} uses {leaked} from viewer.py without importing "
            "them; this raises NameError wherever the code actually runs")


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
                        np.arange(64, dtype=np.int32).reshape((4, 4, 4)))])
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
