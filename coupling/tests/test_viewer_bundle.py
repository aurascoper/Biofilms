"""The viewer bundle contract: what may be drawn, and what may be quoted.

The format exists because of one measured defect. `synthetic_e2e.py` upsamples a
coarse dose field to lattice resolution BEFORE attribution, so its per-cell and
per-lineage numbers are one coarse value repeated across every parcel in a bin.
On a lattice-shaped grid that reads as per-parcel dose. Since the refinement
study the opposite case exists too — a native 4x OpenMC dose field carrying real
structure the CPM grid cannot hold.

Both are legitimate data. Neither is legitimate to display as the other, so the
bundle carries several grids and every layer says which one it is native on.

Bare tier: numpy + h5py, no OpenMC, no marker.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_openmc.viewer import (Grid, Layer, Table, bundle_problems,
                                   manifest, quantitative_layers, read_layer,
                                   read_manifest, write_bundle)

LATTICE = Grid("cpm_labels", (4, 4, 4), (0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
# Same volume, four times the bins per axis: 4 * 0.25 == 16 * 0.0625.
REFINED = Grid("dose_refinement_4", (16, 16, 16), (0.0, 0.0, 0.0),
               (0.0625, 0.0625, 0.0625),
               material_resolution_grid="cpm_labels")


def _labels(shape=(4, 4, 4)):
    return Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                 np.arange(int(np.prod(shape)), dtype=np.int32).reshape(shape),
                 background=0)


def _dose(grid="dose_refinement_4", shape=(16, 16, 16)):
    return Layer("dose_per_source", grid, "Gy/source-particle", "intensive",
                 np.random.default_rng(0).random(shape))


def test_a_valid_multi_grid_bundle_round_trips(tmp_path):
    path = tmp_path / "bundle.h5"
    doc = write_bundle(path, [LATTICE, REFINED], [_labels(), _dose()],
                       provenance={"reference_system_id": "synthetic_biofilm_e2e",
                                   "target_calibration": False})
    assert read_manifest(path) == doc

    data, attrs = read_layer(path, "cell_id")
    assert data.shape == (4, 4, 4)
    assert np.array_equal(data, _labels().data), "logical (x,y,z) order survives"
    assert attrs["semantic_kind"] == "categorical"
    assert attrs["grid_id"] == "cpm_labels"


def test_upsampling_may_be_drawn_and_may_not_be_quoted():
    """THE RULE THE WHOLE FORMAT EXISTS TO CARRY. Upsampling broadcasts one
    value across bins that were never separately measured, so it INVENTS
    information — enforced rather than left to whoever assembles the bundle."""
    upsampled = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                      "intensive", np.zeros((4, 4, 4)),
                      native=False, authoritative_for_quantitation=True,
                      source_grid_id="dose_refinement_4",
                      derivation="upsampled_coarse_dose")
    problems = bundle_problems([LATTICE, REFINED], [upsampled])
    assert any("may be drawn and may not be quoted" in p for p in problems)


def test_a_correct_reduction_may_still_be_quoted():
    """Reduction destroys information and invents none. The refinement study
    confirmed a 4x dose field summed back to the lattice reproduces the native
    lattice field to nine significant figures, so forbidding this would be
    blunter than the physics."""
    reduced = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                    "intensive", np.zeros((4, 4, 4)),
                    native=False, authoritative_for_quantitation=True,
                    source_grid_id="dose_refinement_4",
                    derivation="mass_weighted_mean")
    assert bundle_problems([LATTICE, REFINED], [reduced]) == []


def test_a_reduction_that_does_not_reduce_is_refused():
    """A `block_sum` onto a FINER grid is upsampling wearing the wrong label,
    and that label is the only thing standing between it and being quoted."""
    backwards = Layer("dose_fine", "dose_refinement_4", "Gy/source-particle",
                      "intensive", np.zeros((16, 16, 16)),
                      native=False, authoritative_for_quantitation=True,
                      source_grid_id="cpm_labels",
                      derivation="mass_weighted_mean")
    assert any("not coarser" in p
               for p in bundle_problems([LATTICE, REFINED], [backwards]))


def test_a_derived_layer_that_is_honest_is_accepted():
    honest = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                   "intensive", np.zeros((4, 4, 4)),
                   native=False, authoritative_for_quantitation=False,
                   source_grid_id="dose_refinement_4",
                   derivation="upsampled_coarse_dose")
    assert bundle_problems([LATTICE, REFINED], [honest]) == []
    doc = manifest([LATTICE, REFINED], [honest])
    assert quantitative_layers(doc) == [], (
        "an upsampled layer must not appear in the quotable set")


def test_a_derived_layer_must_say_how_and_from_where():
    for kwargs, expected in (
            ({"native": False, "authoritative_for_quantitation": False,
              "source_grid_id": "dose_refinement_4"}, "does not say how"),
            ({"native": False, "authoritative_for_quantitation": False,
              "derivation": "display_resampling"}, "from which grid"),
            ({"derivation": "display_resampling"}, "cannot both be true")):
        layer = Layer("x", "cpm_labels", "Gy/source-particle", "intensive",
                      np.zeros((4, 4, 4)), **kwargs)
        assert any(expected in p for p in bundle_problems([LATTICE, REFINED],
                                                          [layer]))


def test_labels_may_not_be_resampled_by_arithmetic():
    """A mean of two cell ids is a third cell id that was never in the data —
    a defect with no shape mismatch to give it away."""
    smeared = Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                    np.zeros((4, 4, 4), dtype=np.int32),
                    native=False, authoritative_for_quantitation=False,
                    source_grid_id="dose_refinement_4",
                    derivation="mass_weighted_mean")
    assert any("arithmetic on labels" in p
               for p in bundle_problems([LATTICE], [smeared]))


def test_a_categorical_layer_must_hold_integers():
    floaty = Layer("cell_id", "cpm_labels", "dimensionless", "categorical",
                   np.zeros((4, 4, 4), dtype=float))
    assert any("labels must be integers" in p
               for p in bundle_problems([LATTICE], [floaty]))


def test_grids_in_one_bundle_must_describe_the_same_volume():
    """Grids may sample the volume differently; a bundle whose grids describe
    different volumes is two datasets in one file, and any overlay drawn from it
    misplaces things."""
    elsewhere = Grid("other", (4, 4, 4), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    problems = bundle_problems([LATTICE, elsewhere], [])
    assert any("same physical volume" in p for p in problems)

    shifted = Grid("shifted", (4, 4, 4), (5.0, 0.0, 0.0), (0.25, 0.25, 0.25))
    assert any("one origin" in p for p in bundle_problems([LATTICE, shifted], []))


def test_shape_must_match_the_grid_it_names():
    wrong = Layer("dose", "cpm_labels", "Gy/source-particle", "intensive",
                  np.zeros((8, 8, 8)))
    assert any("but grid" in p for p in bundle_problems([LATTICE], [wrong]))


def test_unknown_grid_kind_and_blank_unit_are_refused():
    problems = bundle_problems([LATTICE], [
        Layer("a", "nope", "Gy", "intensive", np.zeros((4, 4, 4))),
        Layer("b", "cpm_labels", "Gy", "vibes", np.zeros((4, 4, 4))),
        Layer("c", "cpm_labels", "", "intensive", np.zeros((4, 4, 4))),
    ])
    assert any("not declared in this bundle" in p for p in problems)
    assert any("semantic_kind" in p for p in problems)
    assert any("blank unit is never" in p for p in problems)


def test_every_problem_is_reported_not_just_the_first():
    """Same reporting contract as the config loader: a caller assembling a
    bundle wants the whole list."""
    problems = bundle_problems([LATTICE], [
        Layer("a", "nope", "", "vibes", np.zeros((2, 2, 2))),
    ])
    # Grid, semantic kind and unit are all reported at once. The SHAPE is not,
    # and that is correct: there is no declared grid to compare it against, so
    # reporting a shape mismatch would be inventing the comparison.
    assert len(problems) == 3
    assert any("not declared" in p for p in problems)
    assert any("semantic_kind" in p for p in problems)
    assert any("blank unit" in p for p in problems)


def test_attribution_is_a_table_and_carries_its_basis(tmp_path):
    """Per-cell dose is a JOIN of a field with labels, and when that field was
    upsampled first it holds one coarse value repeated across every parcel in a
    bin. As a table with `basis` stamped on it, it reads as what it is."""
    table = Table("dose_attribution",
                  {"label": np.array([1, 2, 3]),
                   "dose_per_source": np.array([1.0, 2.0, 3.0])},
                  basis="upsampled_coarse_dose",
                  units={"dose_per_source": "Gy/source-particle"})
    path = tmp_path / "b.h5"
    doc = write_bundle(path, [LATTICE], [_labels()], [table])
    assert doc["tables"][0]["basis"] == "upsampled_coarse_dose"

    import h5py
    with h5py.File(path, "r") as f:
        assert f["tables/dose_attribution"].attrs["basis"] == "upsampled_coarse_dose"
        assert f["tables/dose_attribution/dose_per_source"].attrs["unit"] \
            == "Gy/source-particle"


def test_a_table_with_ragged_columns_is_refused():
    ragged = Table("t", {"a": np.arange(3), "b": np.arange(4)}, basis="computed")
    assert any("differing lengths" in p
               for p in bundle_problems([LATTICE], [], [ragged]))


def test_write_refuses_and_names_every_problem(tmp_path):
    with pytest.raises(ValueError, match="refusing to write a viewer bundle"):
        write_bundle(tmp_path / "b.h5", [LATTICE],
                     [Layer("a", "nope", "", "vibes", np.zeros((2, 2, 2)))])
    assert not (tmp_path / "b.h5").exists(), "a refused bundle leaves no file"


def test_a_refined_grid_declares_what_resolves_its_material():
    """A 4x tally is native transport data, but the material it integrates is
    still piecewise constant at the CPM pitch — so its sub-voxel structure is
    transport structure and a viewer must not read it as biology."""
    doc = manifest([LATTICE, REFINED], [_labels(), _dose()])
    refined = next(g for g in doc["grids"] if g["grid_id"] == "dose_refinement_4")
    assert refined["material_resolution_grid"] == "cpm_labels"

    dangling = Grid("g", (4, 4, 4), (0.0,) * 3, (0.25,) * 3,
                    material_resolution_grid="absent")
    assert any("not declared" in p for p in bundle_problems([dangling], []))


def test_quantitative_layers_selects_on_the_authoritative_flag_alone():
    doc = manifest([LATTICE, REFINED], [
        _labels(),
        _dose(),
        Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle", "intensive",
              np.zeros((4, 4, 4)), native=False,
              authoritative_for_quantitation=False,
              source_grid_id="dose_refinement_4",
              derivation="upsampled_coarse_dose"),
    ])
    assert set(quantitative_layers(doc)) == {"cell_id", "dose_per_source"}


def test_a_correct_reduction_appears_in_the_quotable_set():
    """`native` is descriptive; `authoritative_for_quantitation` decides. Keying
    on both demoted an exact reduction, which the sharpened rule permits."""
    reduced = Layer("dose_on_lattice", "cpm_labels", "Gy/source-particle",
                    "intensive", np.zeros((4, 4, 4)), native=False,
                    authoritative_for_quantitation=True,
                    source_grid_id="dose_refinement_4",
                    derivation="mass_weighted_mean")
    doc = manifest([LATTICE, REFINED], [reduced])
    assert quantitative_layers(doc) == ["dose_on_lattice"]


def _sample_two_way_raw_state(shape=(6, 6)):
    cell_id = np.zeros(shape, dtype=np.int32)
    cell_id[0, 0] = -1
    cell_id[1, 1] = 1
    cell_id[2, 2] = 4
    cell_id[3, 3] = 5
    species_id = np.zeros(shape, dtype=np.int32)
    species_id[0, 0] = -1
    species_id[1, 1] = 1
    species_id[2, 2] = 4
    species_id[3, 3] = 5
    fast_flux = np.full(shape, 0.8, dtype=np.float64)
    thermal_flux = np.full(shape, 0.2, dtype=np.float64)
    precursor_M = np.full(shape, 0.4, dtype=np.float64)
    inhibitor_C = np.full(shape, 0.5, dtype=np.float64)
    return {
        "cell_id": cell_id,
        "species_id": species_id,
        "fast_flux": fast_flux,
        "thermal_flux": thermal_flux,
        "precursor_damage_M": precursor_M,
        "xenon_inhibitor_C": inhibitor_C,
    }


def test_compute_and_validate_layers_satisfies_all_contracts():
    from biofilm_openmc.viewer import (
        LAYER_CONTRACTS,
        MELANIZED_FRAC_KEY,
        compute_and_validate_layers,
    )

    raw = _sample_two_way_raw_state()
    params = {
        "I0": 1.0,
        "kappa_water": 0.35,
        "kappa_void": 0.02,
        "kappa_species": [0.35, 0.45, 0.50, 1.40, 2.20],
        "sigma_burn": 0.45,
    }
    out = compute_and_validate_layers(raw, params)
    assert set(out.keys()) == set(LAYER_CONTRACTS.keys())
    assert np.isclose(out["local_attenuation_kappa"][0, 0], 0.02)
    assert np.isclose(out["local_attenuation_kappa"][2, 2], 1.40)
    assert np.isclose(out["local_attenuation_kappa"][3, 3], 2.20)
    assert np.isclose(out["bilinear_burnout_sink"][1, 1], 0.0)
    assert np.isclose(out["bilinear_burnout_sink"][2, 2], 0.45 * 1.0 * 0.5)
    assert np.isclose(out[MELANIZED_FRAC_KEY], 2.0 / 3.0)
    assert out["shannon_diversity_H"] > 0.0


def test_compute_and_validate_layers_edge_cases_and_negative_controls():
    from biofilm_openmc.observer import BiofilmOpenMCObserver
    from biofilm_openmc.viewer import (
        MELANIZED_FRAC_KEY,
        compute_and_validate_layers,
    )

    # Empty biomass (all water / voids): no division by zero
    raw_empty = _sample_two_way_raw_state()
    raw_empty["cell_id"][:] = 0
    raw_empty["species_id"][:] = 0
    out_empty = compute_and_validate_layers(raw_empty, {})
    assert out_empty[MELANIZED_FRAC_KEY] == 0.0
    assert out_empty["shannon_diversity_H"] == 0.0

    # Shape mismatch rejected
    raw_bad_shape = _sample_two_way_raw_state()
    raw_bad_shape["fast_flux"] = np.ones((4, 4), dtype=np.float64)
    with pytest.raises(ValueError, match="identical dimensions"):
        compute_and_validate_layers(raw_bad_shape, {})

    # NaN in quotable layer rejected
    raw_nan = _sample_two_way_raw_state()
    raw_nan["xenon_inhibitor_C"][1, 1] = np.nan
    with pytest.raises(ValueError, match="contains NaN/Inf"):
        compute_and_validate_layers(raw_nan, {})

    # Observer rejects unvalidated or incomplete payloads
    with pytest.raises(TypeError, match="validated layer dictionary"):
        BiofilmOpenMCObserver(["not_a_dict"])
    with pytest.raises(KeyError, match="missing"):
        BiofilmOpenMCObserver({"cell_id": raw_empty["cell_id"]})
    obs = BiofilmOpenMCObserver(out_empty)
    assert obs.layers is out_empty


def test_biofilm_openmc_observer_four_panel_audit_with_stub_pyvista(
    monkeypatch, tmp_path
):
    import sys
    import types

    from biofilm_openmc.observer import BiofilmOpenMCObserver
    from biofilm_openmc.viewer import compute_and_validate_layers

    class _StubContours:
        n_points = 4

    class _StubPointGrid:
        def contour(self, isosurfaces=4, scalars=None):
            return _StubContours()

    class _StubImageData:
        def __init__(self, dimensions, spacing, origin):
            self.dimensions = dimensions
            self.spacing = spacing
            self.origin = origin
            self.cell_data = {}

        def cell_data_to_point_data(self):
            return _StubPointGrid()

    class _StubPlotter:
        def __init__(self, shape=(2, 2), off_screen=True, window_size=(1600, 1200)):
            self.shape = shape
            self.off_screen = off_screen
            self.window_size = window_size
            self.subplots_visited = []
            self.titles = []
            self.meshes = []
            self.legends = []
            self.texts = []
            self.views_linked = False
            self.screenshots = []

        def subplot(self, r, c):
            self.subplots_visited.append((r, c))

        def add_title(self, title, font_size=11):
            self.titles.append(title)

        def add_mesh(self, mesh, **kwargs):
            self.meshes.append((mesh, kwargs))

        def add_legend(self, entries, **kwargs):
            self.legends.append(entries)

        def add_text(self, text, **kwargs):
            self.texts.append(text)

        def link_views(self):
            self.views_linked = True

        def screenshot(self, path):
            self.screenshots.append(str(path))

    stub_pv = types.ModuleType("pyvista")
    stub_pv.ImageData = _StubImageData
    stub_pv.Plotter = _StubPlotter
    monkeypatch.setitem(sys.modules, "pyvista", stub_pv)

    raw = _sample_two_way_raw_state()
    validated = compute_and_validate_layers(raw, {})
    obs = BiofilmOpenMCObserver(validated)
    shot = tmp_path / "four_panel.png"
    pl = obs.render_four_panel_audit(screenshot_path=str(shot), off_screen=True)
    assert pl.subplots_visited == [(0, 0), (0, 1), (1, 0), (1, 1)]
    assert len(pl.titles) == 4
    assert len(pl.legends) == 1
    assert pl.views_linked is True
    assert pl.screenshots == [str(shot)]
    assert "mean_dose_rate" in pl.telemetry


def test_continuum_melanin_and_threshold_clipped_observer_telemetry(
    monkeypatch, tmp_path
):
    import sys
    import types

    from biofilm_openmc.observer import (
        SPECIES_ANNOTATIONS,
        BiofilmOpenMCObserver,
    )
    from biofilm_openmc.viewer import compute_and_validate_layers

    # Exact 64x64 snapshot from the specification
    N = 64
    spin = np.zeros((N, N), dtype=np.int32)
    spin[26:38, 26:38] = 1
    spin[16:26, 28:36] = 2
    species_lut = np.array([4, 2], dtype=np.int32)

    melanin = np.where(spin == 1, 0.75, np.where(spin == 2, 0.15, 0.0))
    y, x = np.mgrid[0:N, 0:N]
    r = np.hypot(x - 31.5, y - 31.5)
    radiation = np.exp(-1.5 * r / N)
    poison = 0.4 * np.exp(-((r - 12.0) / 6.0) ** 2)

    raw_state = {
        "spin": spin,
        "species": species_lut,
        "melanin": melanin,
        "radiation": radiation,
        "poison": poison,
    }
    params = {
        "kappa": 1.5,
        "kappa_void": 0.15,
        "kappa_cell": 1.0,
        "kappa_mel": 2.5,
        "sigma_burn": 0.65,
    }

    layers = compute_and_validate_layers(raw_state, params)
    assert np.isclose(np.mean(layers["radiation"]), 0.5766, atol=1e-4)

    # Negative control: kappa = 0.0 -> local_attenuation == 0 everywhere
    layers_neg = compute_and_validate_layers(raw_state, {"kappa": 0.0})
    assert np.all(layers_neg["local_attenuation"] == 0.0)

    class _StubContours:
        n_points = 4

    class _StubPointGrid:
        def contour(self, isosurfaces=4, scalars=None):
            return _StubContours()

    class _StubFullMesh:
        def __init__(self, dimensions=(65, 65, 9), spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
            self.dimensions = dimensions
            self.spacing = spacing
            self.origin = origin
            self.cell_data = {}
            self.n_cells = 224
            self.threshold_calls = []
            self.clip_calls = []

        def threshold(self, value=0.5, scalars=None):
            self.threshold_calls.append((value, scalars))
            return self

        def clip(self, normal="-y", crinkle=False):
            self.clip_calls.append((normal, crinkle))
            return self

        def cell_data_to_point_data(self):
            return _StubPointGrid()

    created_plotters = []

    class _StubPlotter:
        def __init__(self, shape=(2, 2), off_screen=True, window_size=(1600, 1200)):
            self.meshes = []
            self.closed = False
            created_plotters.append(self)

        def subplot(self, r, c):
            pass

        def add_title(self, title, font_size=11):
            pass

        def add_mesh(self, mesh, **kwargs):
            self.meshes.append((mesh, kwargs))

        def add_legend(self, entries, **kwargs):
            pass

        def add_text(self, text, **kwargs):
            pass

        def link_views(self):
            pass

        def screenshot(self, path):
            pass

        def close(self):
            self.closed = True

    stub_pv = types.ModuleType("pyvista")
    stub_pv.ImageData = _StubFullMesh
    stub_pv.Plotter = _StubPlotter
    monkeypatch.setitem(sys.modules, "pyvista", stub_pv)

    obs = BiofilmOpenMCObserver(spacing=(1.0, 1.0, 1.0))
    metrics = obs.render_four_panel_audit(
        raw_state, params, screenshot_path=str(tmp_path / "audit.png")
    )
    assert round(metrics["mean_dose_rate"], 4) == 0.5766
    assert round(metrics["void_fraction"], 4) == 0.9453
    assert round(metrics["mean_kappa_occupied"], 4) == 3.5089
    assert round(metrics["mean_kappa_void"], 4) == 0.225
    assert round(metrics["poison_inventory"], 4) == 320.8142
    assert round(metrics["mean_burnout_rate"], 4) == 0.0094

    pl = created_plotters[-1]
    mesh0, kwargs0 = pl.meshes[0]
    assert mesh0.dimensions == (65, 65, 9)
    assert (0.5, "species_id") in mesh0.threshold_calls
    assert ("-y", True) in mesh0.clip_calls
    assert ("-y", False) in mesh0.clip_calls
    assert kwargs0.get("categories") is True
    assert kwargs0.get("annotations") == SPECIES_ANNOTATIONS


def test_two_way_attenuation_and_observer_scientific_review_contracts():
    """Verify the five scientific and software review contracts for the
    two-way attenuation cascade, poison-burnout surrogate, display vs native
    OpenMC grid distinction, read-only observer isolation, and synthetic
    surrogate dose labeling.
    """
    from biofilm_openmc.observer import BiofilmOpenMCObserver
    from biofilm_openmc.viewer import (
        LAYER_CONTRACTS,
        MELANIZED_FRAC_KEY,
        compute_and_validate_layers,
    )

    # 1. Two-group attenuation cascade conservation and normalization
    raw = _sample_two_way_raw_state(shape=(8, 8))
    raw_copy = {k: v.copy() for k, v in raw.items()}
    params = {
        "I0": 1.0,
        "kappa_water": 0.35,
        "kappa_void": 0.02,
        "kappa_species": [0.35, 0.45, 0.50, 1.40, 2.20],
        "sigma_burn": 0.45,
    }
    out = compute_and_validate_layers(raw, params)
    assert np.all(out["fast_flux"] >= 0.0) and np.all(out["fast_flux"] <= params["I0"])
    assert np.all(out["thermal_flux"] >= 0.0)
    assert np.allclose(out["total_radiation"], out["fast_flux"] + out["thermal_flux"])
    assert np.all(out["total_radiation"] <= params["I0"] + 1e-12)

    # 2. Energy/state accounting in the poison-burnout surrogate
    mel_mask = np.isin(raw["species_id"], [4, 5]).astype(np.float64)
    expected_burn = params["sigma_burn"] * out["total_radiation"] * mel_mask * raw["xenon_inhibitor_C"]
    assert np.allclose(out["bilinear_burnout_sink"], expected_burn)
    assert np.all(out["bilinear_burnout_sink"] >= 0.0)

    # Discrete zero-flux Neumann 5-point Laplacian sum is identically zero
    C = raw["xenon_inhibitor_C"].copy()
    C[2, 3] = 1.25
    C[0, 0] = 0.05
    C[7, 7] = 0.90
    nx, ny = C.shape
    lap = np.zeros_like(C)
    for i in range(nx):
        for j in range(ny):
            im1 = max(0, i - 1)
            ip1 = min(nx - 1, i + 1)
            jm1 = max(0, j - 1)
            jp1 = min(ny - 1, j + 1)
            lap[i, j] = C[im1, j] + C[ip1, j] + C[i, jm1] + C[i, jp1] - 4.0 * C[i, j]
    assert np.isclose(np.sum(lap), 0.0, atol=1e-14)

    # 3. Native OpenMC fields vs resampled/extruded display fields
    assert BiofilmOpenMCObserver.NATIVE_OPENMC_FIELD is False
    assert BiofilmOpenMCObserver.GRID_REPRESENTATION == "resampled_extruded_display_grid"
    for spec in LAYER_CONTRACTS.values():
        assert spec.native_openmc_field is False
        assert spec.grid_representation == "resampled_extruded_display_grid"

    # 4. Observer remains strictly read-only and non-authoritative for quantitation
    assert BiofilmOpenMCObserver.READ_ONLY_OBSERVER is True
    assert BiofilmOpenMCObserver.AUTHORITATIVE_FOR_QUANTITATION is False
    for k, arr in raw.items():
        assert np.array_equal(arr, raw_copy[k])
    out["total_radiation"][0, 0] = 999.0
    assert np.array_equal(raw["fast_flux"], raw_copy["fast_flux"])

    # 5. Synthetic neutron/dose outputs cannot be mistaken for measured physical dose
    assert BiofilmOpenMCObserver.SYNTHETIC_SURROGATE_ONLY is True
    assert BiofilmOpenMCObserver.MEASURED_PHYSICAL_DOSE is False
    for flux_key in ("fast_flux", "thermal_flux", "total_radiation"):
        spec = LAYER_CONTRACTS[flux_key]
        assert spec.synthetic_surrogate is True
        assert spec.measured_physical_dose is False
        assert spec.units == "rel_flux"
        assert "synthetic surrogate" in spec.description
    assert MELANIZED_FRAC_KEY == "melanized_fraction"

