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
                 np.arange(int(np.prod(shape)), dtype=np.int32).reshape(shape))


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
