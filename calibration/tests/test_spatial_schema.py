"""The shipped spatial tables, and the gate over them."""

from __future__ import annotations

import pytest

from biofilm_calibration.schema import read_table
from biofilm_calibration.spatial import report
from biofilm_calibration.spatial.schema import (BIOFILM_STRUCTURE,
                                                ENTITY_SEMANTICS,
                                                OBJECT_MORPHOLOGY,
                                                SAMPLE_METADATA, SOURCES)

JULIA_SPECIES = ["CN", "DR", "CS", "BS", "AN", "SO", "OI"]


@pytest.fixture(scope="module")
def spatial_dir(data_dir):
    return data_dir / "spatial"


def test_entity_semantics_covers_the_authoritative_registry(spatial_dir):
    """The Julia seven-species registry is authoritative (audit §7); a missing
    species would leave a cell ID with no declared meaning."""
    rows = read_table(spatial_dir / "entity_semantics.csv", ENTITY_SEMANTICS)
    assert [r["species_id"] for r in rows] == JULIA_SPECIES


def test_no_species_claims_a_biological_generation(spatial_dir):
    """The model has no growth dynamics and no shape term, so a literal
    generation claim would be unsupportable — biofilms_potts.jl says as much at
    L1578-1583."""
    rows = read_table(spatial_dir / "entity_semantics.csv", ENTITY_SEMANTICS)
    for r in rows:
        assert r["literal_generation_claim_allowed"] == "false", r["species_id"]
        assert r["lineage_semantics"] == "computational", r["species_id"]
        assert r["entity_kind"] == "biomass_parcel", r["species_id"]


def test_every_entity_row_resolves_to_a_source(spatial_dir):
    rows = read_table(spatial_dir / "entity_semantics.csv", ENTITY_SEMANTICS)
    known = {r["source_id"] for r in read_table(spatial_dir / "sources.csv", SOURCES)}
    for r in rows:
        assert r["source_id"] in known, r["species_id"]


@pytest.mark.parametrize("name,schema", [
    ("sample_metadata", SAMPLE_METADATA),
    ("object_morphology", OBJECT_MORPHOLOGY),
    ("biofilm_structure", BIOFILM_STRUCTURE),
])
def test_measurement_tables_are_valid_and_empty(spatial_dir, name, schema):
    """Header-only is the honest state: no morphology has been measured. The
    schemas must still validate, so the day data arrives it is checked."""
    assert read_table(spatial_dir / f"{name}.csv", schema) == []


def test_gate_is_provisional_not_ready(spatial_dir, config_dir):
    """Semantics are declared but nothing is measured, so the gate must not
    claim readiness for time calibration."""
    gate = report.evaluate(
        spatial_dir, config_dir / "cpm_spatial_acceptance_template.toml")
    assert gate.verdict == report.PROVISIONAL
    joined = " ".join(gate.blockers)
    assert "object_morphology is empty" in joined
    assert "acceptance thresholds unset" in joined
    assert "domain semantics undeclared" in joined
    assert any("entity semantics declared" in r for r in gate.reasons)


def test_gate_returns_model_revision_required_for_a_literal_cell_claim(tmp_path):
    """If literal genealogy is ever asserted, the gate must stop rather than
    calibrate a representation that cannot carry the claim."""
    src = (tmp_path / "entity_semantics.csv")
    header = ",".join(ENTITY_SEMANTICS.names)
    src.write_text(
        header + "\n"
        "CN,C. neoformans,biological_cell,cell_division,biological,rod,true,"
        "declared,S1,ready,claims literal generations\n")
    gate = report.evaluate(tmp_path)
    assert gate.verdict == report.MODEL_REVISION_REQUIRED
    assert "no shape or connectivity term" in " ".join(gate.blockers)


def test_gate_is_not_evaluated_without_semantics(tmp_path):
    gate = report.evaluate(tmp_path)
    assert gate.verdict == report.NOT_EVALUATED
