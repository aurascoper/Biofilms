"""Mass closure, the export gate, seeded uncertainty, and the shipped tables."""

from __future__ import annotations

import math

import pytest

from biofilm_calibration.materials.basis_conversion import BasisError
from biofilm_calibration.materials import report
from biofilm_calibration.materials.export import (MaterialClass, problems_for,
                                                  to_config_fragment)
from biofilm_calibration.materials.mixture import (Component, blend,
                                                   ideal_mixture_density)
from biofilm_calibration.materials.schema import (BULK_MEASUREMENTS,
                                                  COMPONENT_DEFINITIONS,
                                                  ELEMENTAL_ANALYSIS,
                                                  METAL_LOADING,
                                                  SAMPLE_METADATA)
from biofilm_calibration.materials.uncertainty import (propagate,
                                                       relative_uncertainty)
from biofilm_calibration.schema import read_table

WATER = {"H": 0.111894, "O": 0.888106}
BIOMASS = {"H": 0.10, "C": 0.50, "N": 0.10, "O": 0.30}


def test_blend_reproduces_an_analytic_mixture():
    """90% water + 10% dry biomass, worked by hand."""
    comps = [Component("water", 0.9, WATER, density_g_cm3=1.0),
             Component("dry_biomass", 0.1, BIOMASS, density_g_cm3=1.4)]
    got = blend(comps)
    assert got["C"] == pytest.approx(0.1 * 0.50)
    assert got["H"] == pytest.approx(0.9 * 0.111894 + 0.1 * 0.10)
    assert got["O"] == pytest.approx(0.9 * 0.888106 + 0.1 * 0.30)
    assert sum(got.values()) == pytest.approx(1.0)


def test_component_fractions_must_close():
    """A forgotten component would otherwise be absorbed as a silent
    renormalization."""
    comps = [Component("water", 0.9, WATER),
             Component("dry_biomass", 0.05, BIOMASS)]      # 0.95, not 1
    with pytest.raises(BasisError, match="silent renormalization"):
        blend(comps)


def test_component_elemental_fractions_must_close():
    with pytest.raises(BasisError, match="sum to"):
        Component("bad", 1.0, {"C": 0.5, "H": 0.2})


def test_negative_fractions_are_rejected():
    with pytest.raises(BasisError, match="negative"):
        Component("bad", 1.0, {"C": 1.2, "H": -0.2})
    with pytest.raises(BasisError, match="outside"):
        Component("bad", 1.5, WATER)


def test_metal_double_counted_in_ash_and_assay_is_caught():
    """Easy to create: ash comes from combustion and the metal is also assayed
    directly, so it lands in both components."""
    comps = [Component("water", 0.8, WATER),
             Component("ash", 0.1, {"Th": 1.0},
                       contributes_metals=frozenset({"Th"})),
             Component("sorbed_metal", 0.1, {"Th": 1.0},
                       contributes_metals=frozenset({"Th"}))]
    with pytest.raises(BasisError, match="counted twice"):
        blend(comps)


def test_ideal_mixture_density_is_labelled_a_proxy():
    comps = [Component("water", 0.9, WATER, density_g_cm3=1.0),
             Component("dry_biomass", 0.1, BIOMASS, density_g_cm3=1.4)]
    value, evidence = ideal_mixture_density(comps)
    assert evidence == "derived_proxy"
    assert 1.0 < value < 1.05
    # volume additivity, by hand
    assert value == pytest.approx(1.0 / (0.9 / 1.0 + 0.1 / 1.4))


def test_proxy_material_cannot_be_exported_as_calibrated():
    """THE export guard: a plausible-looking proxy composition must not be
    laundered into a measurement by being shipped."""
    m = MaterialClass("hydrated biofilm", 1.02, blend(
        [Component("water", 0.9, WATER), Component("dry_biomass", 0.1, BIOMASS)]),
        evidence_basis="proxy", source_id="S1", system_conditions="30 C, pH 7")
    with pytest.raises(BasisError, match="cannot be exported"):
        to_config_fragment(m, "baseline_biomass")
    assert any("proxy" in p for p in problems_for(m))


def test_export_requires_attribution_and_conditions():
    els = blend([Component("water", 0.9, WATER),
                 Component("dry_biomass", 0.1, BIOMASS)])
    no_src = MaterialClass("m", 1.02, els, "direct_measurement", "", "30 C")
    assert any("attributable" in p for p in problems_for(no_src))
    no_cond = MaterialClass("m", 1.02, els, "direct_measurement", "S1", "")
    assert any("system_conditions" in p for p in problems_for(no_cond))


def test_explicit_components_kind_is_incompatible_with_the_mapper():
    els = blend([Component("water", 0.9, WATER),
                 Component("dry_biomass", 0.1, BIOMASS)])
    m = MaterialClass("m", 1.02, els, "direct_measurement", "S1", "30 C",
                      material_model_kind="explicit_components")
    assert any("sub-voxel homogenization" in p for p in problems_for(m))


def test_a_measured_material_exports_in_the_loader_s_shape():
    els = blend([Component("water", 0.9, WATER),
                 Component("dry_biomass", 0.1, BIOMASS)])
    m = MaterialClass("hydrated biofilm", 1.02, els, "direct_measurement",
                      "S1", "30 C, pH 7, aerobic")
    frag = to_config_fragment(m, "baseline_biomass")
    assert "[materials.baseline_biomass]" in frag
    assert "density_g_cm3 = 1.02" in frag
    assert "[materials.baseline_biomass.elements]" in frag
    # the coupling loader requires a closed composition
    assert sum(els.values()) == pytest.approx(1.0)


def test_uncertainty_is_seeded_and_reproducible():
    """An interval that moves between runs cannot be reviewed or cited."""
    inputs = {"m": (1.0, 0.05), "v": (10.0, 0.3)}
    fn = lambda m, v: m / v
    a = propagate(fn, inputs, seed=7, draws=4000)
    b = propagate(fn, inputs, seed=7, draws=4000)
    assert a == b
    c = propagate(fn, inputs, seed=8, draws=4000)
    assert c != a
    assert a["mean"] == pytest.approx(0.1, rel=0.02)
    assert 0.0 < relative_uncertainty(a) < 0.2


def test_exact_inputs_give_zero_spread():
    out = propagate(lambda m, v: m / v, {"m": (1.0, 0), "v": (10.0, None)},
                    seed=1, draws=100)
    assert out["sd"] == pytest.approx(0.0, abs=1e-15)   # float noise, not spread
    assert out["mean"] == pytest.approx(0.1)
    assert out["n_valid"] == 100


@pytest.fixture(scope="module")
def materials_dir(data_dir):
    return data_dir / "materials"


@pytest.mark.parametrize("name,schema", [
    ("sample_metadata", SAMPLE_METADATA),
    ("bulk_measurements", BULK_MEASUREMENTS),
    ("elemental_analysis", ELEMENTAL_ANALYSIS),
    ("metal_loading", METAL_LOADING),
    ("component_definitions", COMPONENT_DEFINITIONS),
])
def test_shipped_tables_validate_and_are_empty(materials_dir, name, schema):
    """Header-only is honest: nothing has been measured."""
    assert read_table(materials_dir / f"{name}.csv", schema) == []


def test_gates_are_provisional_and_blocked(materials_dir):
    """The expected, valid outcome: transport is merely short of data, while
    the radiodialysis mapping is blocked by a units error in the model."""
    gates = report.evaluate(materials_dir)
    assert gates.openmc == report.PROVISIONAL
    assert gates.radiodialysis == report.BLOCKED
    joined = " ".join(gates.radiodialysis_blockers)
    assert "site occupancy" in joined and "units error" in joined
