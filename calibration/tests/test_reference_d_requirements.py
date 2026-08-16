"""The frozen Reference D protocol, checked against the contract it describes.

A measurement protocol is a document, and documents drift. The failure is quiet
and expensive: the table keeps listing a field nothing requires any more, or
stops listing one that was added, and nobody finds out until a campaign comes
back missing something that cannot be measured again on the same coupons.

So the freeze is enforced rather than asserted. These tests fail if the table
and the live contract disagree in either direction.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "calibration" / "scripts"))

import reference_d_status as status  # noqa: E402


@pytest.fixture(scope="module")
def requirements():
    return status.load_requirements()


def test_the_table_is_well_formed(requirements):
    # load_requirements() raises on a bad vocabulary, a duplicate id or an
    # empty artifact/acceptance/unblocks cell, so reaching here is the check.
    assert requirements
    assert len({r["requirement_id"] for r in requirements}) == len(requirements)


def test_every_required_field_has_a_requirement(requirements):
    """The direction that costs a repeated campaign: the code refuses on a
    field and the protocol never asked anyone to measure it."""
    spatial, material = status.spatial_report, status.material_report
    s = spatial.evaluate(status.DATA / "spatial",
                         REPO / "config" / "cpm_spatial_acceptance_template.toml")
    m = material.evaluate(status.DATA / "materials")
    fields = status.required_fields(status.declared_binding(), status.target_spec(),
                                    (s.verdict, m.openmc))
    ledgered = {r["config_key"] for r in requirements}
    assert not (fields - ledgered), \
        f"the contract requires fields the table does not name: {sorted(fields - ledgered)}"


def test_the_check_mode_passes(requirements):
    """The same cross-check the script runs, as a test rather than a habit."""
    assert status.main(["--check"]) == 0


def test_every_requirement_lands_somewhere_nameable(requirements):
    """A requirement whose destination is vague is a requirement nobody can
    close. Each row lands in a file, or explicitly in the emitted config —
    never in a description of a file."""
    direct = "enters the config directly"
    for r in requirements:
        lands = r["lands_in"]
        assert ".csv" in lands or ".toml" in lands or direct in lands, \
            f"{r['requirement_id']} lands in {lands!r}, which is neither a file " \
            "nor the config"


def test_declared_requirements_are_not_laboratory_work(requirements):
    """`declared` is a MODELLING CHOICE. Sending one to a lab returns a
    measurement of something nobody asked a question about — and worse, one that
    then looks like evidence for a choice that was never evidential."""
    declared = [r for r in requirements if r["supplied_by"] == "declared"]
    assert declared
    for r in declared:
        assert "measure" not in r["artifact"].lower(), \
            f"{r['requirement_id']} is declared but its artifact asks for a measurement"


def test_model_blocked_requirements_are_not_measurements_in_waiting(requirements):
    """A row that awaits a MODEL change must say so, or it ends up in a campaign
    budget as a measurement that would not have unblocked anything."""
    blocked = [r for r in requirements
               if r["status"] == "unsupported_by_current_model"]
    assert blocked, "D-XRED at least should be blocked_on_model"
    for r in blocked:
        assert "not by missing data" in r["notes"] or "units error" in r["notes"], \
            f"{r['requirement_id']} is blocked_on_model but does not say why"


def test_the_two_central_blanks_are_both_covered(requirements):
    """Reference D's two central numbers, named in the reference-system doc as
    the whole remaining gap."""
    keys = {r["config_key"] for r in requirements}
    assert "lattice_pitch_um" in keys
    assert "density_g_cm3" in keys


def test_satisfied_rows_really_are_satisfied(requirements):
    """A row marked satisfied must not appear in the live gate's blockers."""
    s = status.spatial_report.evaluate(
        status.DATA / "spatial",
        REPO / "config" / "cpm_spatial_acceptance_template.toml")
    blockers = " ".join(s.blockers)
    satisfied = {r["requirement_id"] for r in requirements
                 if r["status"] == "satisfied"}
    # entity semantics and the time observable are the two that currently pass.
    assert "D-ENTITY" in satisfied and "D-TIMEOBS" in satisfied
    assert "entity semantics unresolved" not in blockers
    assert "time_observable" not in blockers


def test_the_protocol_document_exists_and_points_at_the_table():
    doc = REPO / "docs" / "calibration" / "reference_d_measurement_protocol.md"
    text = doc.read_text(encoding="utf-8")
    assert "reference_d_requirements.csv" in text
    assert "reference_d_status.py" in text
    # The two things a reader must not miss.
    assert "cells_only" in text          # the refused volume basis
    assert "A1" in text                  # the sealed capsule that cannot be modelled


def test_the_campaign_is_not_blocked_by_its_own_outputs(requirements):
    """D-NUMERICS and D-AXIALFACE need data the campaign produces. If they
    gated its start the programme would deadlock: no campaign until the sweep,
    no sweep until the pitch, no pitch until the campaign."""
    analysis = {r["requirement_id"] for r in requirements
                if r["status"] == "awaiting_analysis"}
    assert "D-NUMERICS" in analysis
    assert "D-AXIALFACE" in analysis
    assert not (analysis & set(status._BLOCKS_CAMPAIGN))


def test_every_modelling_decision_is_frozen(requirements):
    """The point of the branch: nothing is left `awaiting_declaration`, because
    a laboratory cannot answer one and measurements taken before the question
    was framed answer nothing."""
    undecided = [r["requirement_id"] for r in requirements
                 if r["status"] == "awaiting_declaration"]
    assert undecided == []


def test_campaign_readiness_is_separate_from_config_readiness(requirements):
    """Four verdicts, because they gate different things. Collapsing them into
    one flag is how "the software can do feedback" becomes "feedback matters".

    The fourth is institutional: it gates CULTURING, which is not a question
    about measurements at all, and it is the only one whose answer comes from
    outside this repository.
    """
    s = status.spatial_report.evaluate(
        status.DATA / "spatial",
        REPO / "config" / "reference_d_spatial_acceptance.toml")
    m = status.material_report.evaluate(
        status.DATA / "materials",
        REPO / "config" / "reference_d_material_acceptance.toml")
    verdicts = status.readiness(requirements, s.verdict, m.openmc,
                                status.declared_binding())
    assert set(verdicts) == {status.AUTHORIZED, status.CAMPAIGN_READY,
                             status.CONFIG_READY, status.SWEEP_READY}
    # Nothing is authorised, and every criterion says why.
    authorized_ok, authorized_blockers = verdicts[status.AUTHORIZED]
    assert not authorized_ok
    assert len(authorized_blockers) == len(status.AUTHORIZATION_CRITERIA) == 9

    # The campaign waits on the institutional approval and nothing measured.
    campaign_ok, campaign_blockers = verdicts[status.CAMPAIGN_READY]
    assert not campaign_ok
    assert any("awaiting_approval" in b for b in campaign_blockers)
    assert any("institutional authorization" in b for b in campaign_blockers)
    assert not any("awaiting_measurement" in b for b in campaign_blockers), (
        "the campaign must not be blocked by the measurements it exists to make")
    # The config waits on the measurements, which is a different threshold.
    assert not verdicts[status.CONFIG_READY][0]
    assert not verdicts[status.SWEEP_READY][0]


def test_authorization_cannot_be_reached_without_a_baseline_row():
    """CAMPAIGN_READY requires the institutional verdict, so the two can never
    disagree — which is the point of deriving one from the other rather than
    letting a reader reconcile them."""
    assert all(not met for _, met in status.authorization_criteria([]))


def test_criteria_with_no_consumer_are_named(requirements):
    """A gate must not close on a threshold nobody compares against."""
    report = status.enforcement_report()
    decorative = set(report["not_implemented"])
    assert "spatial:maximum_axis_ratio_error" in decorative
    assert "spatial:minimum_objects_per_species" in decorative
    assert "spatial:require_independent_validation_sample" in decorative
    assert "spatial:maximum_biovolume_fraction_error" in set(report["hard"])
    assert "spatial:maximum_porosity_error" in set(report["derived"])
