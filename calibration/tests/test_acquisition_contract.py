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


def row_for(schema, **values):
    """Build a CSV row from a schema, so adding a column cannot silently
    shift a hand-counted fixture into the wrong fields."""
    missing = set(values) - set(schema.names)
    assert not missing, f"unknown column(s) for {schema.name}: {missing}"
    return ",".join(str(values.get(c, "")) for c in schema.names)


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
    p.write_text(",".join(BASELINE_CONDITION.names) + "\n" + row_for(
        BASELINE_CONDITION, growth_condition_id="gc1",
        consortium_or_surrogate="P. aeruginosa surrogate",
        strain_identities="PAO1", biosafety_level_by_strain="PAO1:BSL2",
        institutional_approval_id="IBC-2026-001",
        containment_facility="BSL2 core", growth_medium="LB",
        temperature_C="30", pH="7.0", oxygen_condition="aerobic",
        substrate_or_membrane="glass", biofilm_age_h="48",
        flow_condition="static", irradiation="none",
        domain_semantics="microvolume", is_target_system="false",
        evidence_basis="declared", source_id="S1", status="provisional",
        notes="pilot") + "\n")
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
        ",".join(SPA_SAMPLES.names) + "\n" + row_for(
            SPA_SAMPLES, sample_id="s1", culture_batch_id="batch1",
            blank_sample_id="blank1", growth_condition_id="gc1",
            held_out="false", species="Pa", medium="LB", temperature_C="30",
            pH="7.0", oxygen_condition="aerobic", growth_phase="stationary",
            imaging_method="confocal", voxel_size_x_um="0.1",
            voxel_size_y_um="0.1", voxel_size_z_um="0.5",
            segmentation_method="otsu",
            segmentation_basis="whole_biofilm_envelope",
            source_id="S1") + "\n")

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
        ",".join(BULK_MEASUREMENTS.names) + "\n" + row_for(
            BULK_MEASUREMENTS, sample_id="s1", replicate_id="r1",
            culture_batch_id="batch1", paired_sample_group="pair1",
            growth_condition_id="gc1",
            wet_mass_sample_plus_substrate_g="5.0",
            dry_mass_sample_plus_substrate_g="1.0",
            hydrated_volume_cm3="10.0", volume_method="imaging",
            volume_basis="whole_biofilm_envelope",
            volume_support="full_coupon",
            drain_orientation="vertical", drain_time_s="60",
            blot_material="lint-free", blot_contact_time_s="5",
            ambient_temperature_C="21.0", time_to_weighing_s="120",
            drying_protocol="105C", drying_endpoint="constant mass",
            quality_flag="ok", source_id="S1") + "\n")
    gates = material_report.evaluate(tmp_path)
    assert any("name no blank" in b for b in gates.openmc_blockers)


# --- dynamic-observable contract -------------------------------------------

def test_exactly_one_dynamic_observable_is_selected(data_dir):
    """Selected 2026-08-14 rather than left open: the choice determines whether
    the campaign needs time-resolved stacks, and discovering that after a
    static campaign would mean repeating it."""
    from biofilm_calibration.spatial.time_observable import (TIME_OBSERVABLE,
                                                             problems, selectable)
    rows = read_table(data_dir / "spatial" / "time_observable.csv", TIME_OBSERVABLE)
    chosen = [r for r in rows if r["selected"] == "true"]
    assert len(chosen) == 1, "exactly one observable defines S_exp"
    assert selectable(chosen[0])
    assert problems(rows) == []


def test_growth_dominated_observables_are_unsupported_not_blocked(data_dir):
    """A blocked quantity awaits a measurement; an unsupported one awaits a
    model change. No amount of thickening data calibrates a model that cannot
    thicken."""
    from biofilm_calibration.spatial.time_observable import TIME_OBSERVABLE
    rows = {r["observable_id"]: r for r in
            read_table(data_dir / "spatial" / "time_observable.csv", TIME_OBSERVABLE)}
    for oid in ("biofilm_thickening_rate", "specific_growth_rate"):
        assert rows[oid]["represented_by_model"] == "false", oid
        assert rows[oid]["status"] == "unsupported_by_current_model", oid


def test_single_cell_msd_is_semantically_incompatible(data_dir):
    """The conventional choice, and a category mismatch under parcel
    semantics: a cell_id is not an observable bacterium."""
    from biofilm_calibration.spatial.time_observable import (TIME_OBSERVABLE,
                                                             selectable)
    rows = {r["observable_id"]: r for r in
            read_table(data_dir / "spatial" / "time_observable.csv", TIME_OBSERVABLE)}
    msd = rows["single_cell_msd"]
    assert msd["represented_by_model"] == "true"
    assert msd["measurable"] == "true"
    assert msd["semantically_compatible_with_parcel"] == "false"
    assert not selectable(msd)


def test_a_selected_but_unqualified_observable_is_rejected():
    """Selection must satisfy all three conditions, not merely be marked."""
    from biofilm_calibration.spatial.time_observable import problems
    bad = [{"observable_id": "x", "selected": "true",
            "represented_by_model": "false", "measurable": "true",
            "semantically_compatible_with_parcel": "true"}]
    assert any("is selected but fails" in p for p in problems(bad))


def test_spatial_gate_no_longer_blocks_on_the_observable(data_dir, config_dir):
    """That half is now settled; the gate stays PROVISIONAL for the rest."""
    from biofilm_calibration.spatial import report
    gate = report.evaluate(data_dir / "spatial",
                           config_dir / "cpm_spatial_acceptance_template.toml")
    assert gate.verdict == report.PROVISIONAL
    assert not any("dynamic observable" in b for b in gate.blockers)
    assert any("no baseline condition declared" in b for b in gate.blockers)


def test_unsupported_status_cannot_carry_a_value():
    """A new status must not fall through the no-value check the way an
    unlisted one silently would."""
    from biofilm_calibration.schema import (NO_VALUE_STATUSES, SchemaError,
                                            Column, TableSchema, read_table)
    import pytest as _pt
    assert "unsupported_by_current_model" in NO_VALUE_STATUSES
    S = TableSchema(name="t", columns=(Column("id"), Column("v", numeric=True,
                                                            required=False),
                                       Column("status")), key=("id",))
    import tempfile, pathlib as _pl
    with tempfile.TemporaryDirectory() as td:
        p = _pl.Path(td) / "t.csv"
        p.write_text("id,v,status\na,1.5,unsupported_by_current_model\n")
        with _pt.raises(SchemaError, match="has no value"):
            read_table(p, S)


# --- dataset candidate ledger ----------------------------------------------

def test_candidate_ledger_validates_and_records_rejections(data_dir):
    """A search that records only its accepts is not a provenance trail."""
    from biofilm_calibration.spatial.datasets import DATASET_CANDIDATES, problems
    rows = read_table(data_dir / "spatial" / "dataset_candidates.csv",
                      DATASET_CANDIDATES)
    assert problems(rows) == []
    rejected = [r for r in rows if r["verdict"] == "rejected"]
    assert len(rejected) >= 4
    for r in rejected:
        assert r["rationale"].strip(), r["candidate_id"]


def test_no_public_candidate_clears_the_target_gate(data_dir):
    """PUBLIC_COMPONENTS_FOUND_BUT_NOT_PAIRED: a surrogate cannot clear the
    target gate however good its data are."""
    from biofilm_calibration.spatial.datasets import DATASET_CANDIDATES
    rows = read_table(data_dir / "spatial" / "dataset_candidates.csv",
                      DATASET_CANDIDATES)
    assert not any(r["clears_target_gate"] == "true" for r in rows)


def test_a_surrogate_claiming_the_target_gate_is_rejected():
    """The rule is enforced, not conventional."""
    from biofilm_calibration.spatial.datasets import problems
    bad = [{"candidate_id": "x", "clears_target_gate": "true",
            "target_relationship": "surrogate", "verdict": "use_first",
            "rationale": "r", "has_3d_stacks": "true"}]
    assert any("only the target consortium" in p for p in problems(bad))


def test_the_verified_pilot_candidate_carries_its_evidence(data_dir):
    """Sizes were verified against the repository API rather than transcribed."""
    from biofilm_calibration.spatial.datasets import DATASET_CANDIDATES
    rows = {r["candidate_id"]: r for r in read_table(
        data_dir / "spatial" / "dataset_candidates.csv", DATASET_CANDIDATES)}
    vc = rows["vcholerae_dryad_zcrjdfnph"]
    assert vc["verdict"] == "use_first"
    assert vc["size_bytes"] == 1589998422
    assert vc["size_basis"] == "current_version"
    assert "Dryad API" in vc["verified_via"]
    assert vc["target_relationship"] == "surrogate"


def test_cumulative_size_trap_is_recorded(data_dir):
    """Dryad storageSize is cumulative across versions and overestimated the
    B. subtilis download by 4.3x — a download budget taken from it would be
    badly wrong."""
    from biofilm_calibration.spatial.datasets import DATASET_CANDIDATES
    rows = {r["candidate_id"]: r for r in read_table(
        data_dir / "spatial" / "dataset_candidates.csv", DATASET_CANDIDATES)}
    bs = rows["bsubtilis_dryad_f4qrfj71n"]
    assert bs["size_basis"] == "current_version"
    assert "215" in bs["rationale"] and "cumulative" in bs["rationale"].lower()


def test_holographic_motility_rejection_confirms_the_observable_choice(data_dir):
    """An excellent exact-species dataset that answers the wrong question."""
    from biofilm_calibration.spatial.datasets import DATASET_CANDIDATES
    rows = {r["candidate_id"]: r for r in read_table(
        data_dir / "spatial" / "dataset_candidates.csv", DATASET_CANDIDATES)}
    h = rows["shewanella_holographic_motility"]
    assert h["verdict"] == "rejected"
    assert h["target_relationship"] == "exact_species"
    assert "parcel" in h["rationale"]
