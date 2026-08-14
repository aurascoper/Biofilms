"""Sample linkage, blanks, and the gates that now depend on them.

A pitch and a density that came from different biofilms would produce a
perfectly well-formed calibration describing nothing. The linkage columns and
the pairing blockers are what make that failure visible.
"""

from __future__ import annotations

import pytest

from biofilm_calibration.acquisition import (BASELINE_CONDITION, BLANKS,
                                             linkage_columns)
from biofilm_calibration.materials import report as material_report
from biofilm_calibration.materials.basis_conversion import (BasisError,
                                                            blank_corrected_mass)
from biofilm_calibration.materials.schema import BULK_MEASUREMENTS
from biofilm_calibration.materials.schema import SAMPLE_METADATA as MAT_SAMPLES
from biofilm_calibration.schema import read_table
from biofilm_calibration.spatial import report as spatial_report
from biofilm_calibration.spatial.schema import SAMPLE_METADATA as SPA_SAMPLES

LINKAGE = {"culture_batch_id", "paired_sample_group", "blank_sample_id",
           "measurement_order", "growth_condition_id", "medium_batch_id"}


def test_both_branches_carry_the_same_linkage_columns():
    """Defined once, so the two branches cannot drift apart on the names."""
    assert LINKAGE <= set(SPA_SAMPLES.names)
    assert LINKAGE <= set(MAT_SAMPLES.names)
    assert LINKAGE <= set(BULK_MEASUREMENTS.names)
    assert {c.name for c in linkage_columns()} == LINKAGE


def test_held_out_is_a_column_on_imaged_samples():
    assert "held_out" in SPA_SAMPLES.names


def test_surface_water_procedure_is_structural_not_free_text():
    """Six required columns, because 'wet mass' is not a measurement until the
    surface water is defined."""
    for c in ("drain_orientation", "drain_time_s", "blot_material",
              "blot_contact_time_s", "ambient_temperature_C",
              "time_to_weighing_s"):
        assert c in BULK_MEASUREMENTS.names, c
        assert BULK_MEASUREMENTS.column(c).required, c
    assert "surface_water_removal_protocol" not in BULK_MEASUREMENTS.names


def test_masses_are_recorded_as_weighed():
    """The raw reading includes the substrate; the biofilm mass is derived by
    blank subtraction, so the correction cannot be silently skipped."""
    assert "wet_mass_sample_plus_substrate_g" in BULK_MEASUREMENTS.names
    assert "dry_mass_sample_plus_substrate_g" in BULK_MEASUREMENTS.names
    assert "wet_mass_g" not in BULK_MEASUREMENTS.names


def test_blank_subtraction_is_required_not_defaulted():
    assert blank_corrected_mass(5.0, 2.0, label="wet") == pytest.approx(3.0)
    with pytest.raises(BasisError, match="no blank mass"):
        blank_corrected_mass(5.0, None, label="wet")
    with pytest.raises(BasisError, match="no biofilm"):
        blank_corrected_mass(5.0, 5.0, label="wet")
    with pytest.raises(BasisError, match="negative blank"):
        blank_corrected_mass(5.0, -1.0, label="wet")


def test_three_blank_kinds_are_declared():
    """The hydrated substrate blank is the one most easily skipped and the
    costliest to skip."""
    kinds = BLANKS.column("blank_kind").vocabulary
    assert kinds == {"dry_substrate", "hydrated_substrate",
                     "medium_exposed_abiotic"}


def test_baseline_condition_ships_unset(data_dir):
    """A protocol contract, not a claim about a specimen that exists."""
    assert read_table(data_dir / "baseline_condition.csv", BASELINE_CONDITION) == []


def test_baseline_condition_can_declare_a_surrogate(tmp_path):
    """A surrogate exercises the harness; is_target_system records that it
    cannot clear the gate."""
    p = tmp_path / "b.csv"
    p.write_text(",".join(BASELINE_CONDITION.names) + "\n"
                 "gc1,P. aeruginosa surrogate,PAO1,LB,30,7.0,aerobic,glass,48,"
                 "static,none,microvolume,false,declared,S1,provisional,pilot\n")
    rows = read_table(p, BASELINE_CONDITION)
    assert rows[0]["is_target_system"] == "false"
    assert rows[0]["domain_semantics"] == "microvolume"


def test_domain_semantics_vocabulary_admits_microvolume():
    assert "microvolume" in BASELINE_CONDITION.column("domain_semantics").vocabulary


def test_blanks_table_ships_empty(data_dir):
    assert read_table(data_dir / "materials" / "blanks.csv", BLANKS) == []


def test_spatial_gate_blocks_on_pairing_and_held_out(tmp_path):
    """An imaged stack with no paired sibling and no held-out designation must
    not be able to produce a pitch."""
    from biofilm_calibration.spatial.schema import (BIOFILM_STRUCTURE,
                                                    ENTITY_SEMANTICS,
                                                    OBJECT_MORPHOLOGY)
    from biofilm_calibration.schema import write_template

    src = tmp_path
    (src / "entity_semantics.csv").write_text(
        ",".join(ENTITY_SEMANTICS.names) + "\n"
        "CN,C. neoformans,biomass_parcel,parcel_split,computational,"
        "isotropic_blob,false,declared,S1,ready,\n")
    for schema, name in ((OBJECT_MORPHOLOGY, "object_morphology"),
                         (BIOFILM_STRUCTURE, "biofilm_structure")):
        write_template(src / f"{name}.csv", schema)
    # one imaged sample, unpaired and not held out
    (src / "sample_metadata.csv").write_text(
        ",".join(SPA_SAMPLES.names) + "\n"
        "s1,batch1,,blank1,1,gc1,med1,false,Pa,PAO1,LB,30,7.0,aerobic,"
        "stationary,48,confocal,0.1,0.1,0.5,otsu,v1,r1,S1\n")

    gate = spatial_report.evaluate(src)
    joined = " ".join(gate.blockers)
    assert "no paired_sample_group" in joined
    assert "held_out" in joined


def test_material_gate_blocks_on_unblanked_and_unpaired(data_dir):
    """Against the shipped (empty) tables the emptiness blockers fire; the new
    checks are exercised on populated fixtures above and here for wiring."""
    gates = material_report.evaluate(data_dir / "materials")
    assert gates.openmc == material_report.PROVISIONAL
    assert gates.radiodialysis == material_report.BLOCKED


def test_material_gate_names_unblanked_samples(tmp_path):
    from biofilm_calibration.materials.schema import (COMPONENT_DEFINITIONS,
                                                      ELEMENTAL_ANALYSIS,
                                                      METAL_LOADING)
    from biofilm_calibration.schema import write_template

    for schema, name in ((ELEMENTAL_ANALYSIS, "elemental_analysis"),
                         (METAL_LOADING, "metal_loading"),
                         (COMPONENT_DEFINITIONS, "component_definitions"),
                         (MAT_SAMPLES, "sample_metadata")):
        write_template(tmp_path / f"{name}.csv", schema)
    # a bulk measurement naming no blank
    (tmp_path / "bulk_measurements.csv").write_text(
        ",".join(BULK_MEASUREMENTS.names) + "\n"
        "s1,r1,batch1,pair1,,1,gc1,med1,5.0,1.0,,10.0,,,,imaging,vertical,60,"
        "lint-free,5,21.0,120,105C,constant mass,,ok,S1\n")
    gates = material_report.evaluate(tmp_path)
    assert any("name no blank" in b for b in gates.openmc_blockers)
