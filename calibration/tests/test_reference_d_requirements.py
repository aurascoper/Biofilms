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


def test_a_valid_approval_meets_every_criterion():
    """The suite would be worthless if the milestone could not be reached. It
    could not, for a while: `authorization_criteria` was called without
    `sources`, so criterion 8 was unmeetable by any input — indistinguishable
    from a correctly withheld authorization, because "not authorized" is the
    expected state."""
    from test_approval import SOURCES, TODAY, row

    unmet = [c for c, ok in status.authorization_criteria(
        [row()], SOURCES, today=TODAY) if not ok]
    assert unmet == [], unmet


@pytest.mark.parametrize("field,bad", [
    ("strain_identities", "unknown"),
    ("containment_facility", "TBD"),
    ("risk_assessment_reference", "pending"),
    ("institutional_approval_id", "d approved"),
    ("institutional_approval_authority", "n/a"),
    ("approved_protocol_version", ""),
    ("approval_source_id", "IBC_NOPE"),
    ("approval_effective_date", "June 2026"),
    ("approval_expiration_date", "2026-08-01"),
    ("culturing_start_date", "2026-05-01"),
    ("is_target_system", "false"),
    ("scope_hash", " "),
])
def test_every_refusal_reaches_a_criterion(field, bad):
    """THE NEGATIVE CONTROL THE MAPPING NEVER HAD.

    `authorization_criteria` translates `approval.problems` into nine criteria
    by substring. Any refusal matching no pattern used to be read as evidence
    for nothing, and so defaulted to MET — a false positive on an institutional
    biosafety milestone.

    That is exactly what happened: adding the `approved_protocol_version` check
    to the gate produced a refusal the milestone could not see. The gate said
    no; all nine criteria said yes. Nothing failed, because nothing compared
    them.

    So break one field at a time and require the milestone to notice every
    single one. A new check in `approval.problems` that lands nowhere now fails
    here instead of silently widening the gap.
    """
    from biofilm_calibration import approval
    from test_approval import SOURCES, TODAY, row

    r = row(**{field: bad})
    assert approval.problems([r], SOURCES, today=TODAY), (
        "the fixture is stale: this input no longer refuses at the gate, so "
        "the test proves nothing about the criteria")
    unmet = [c for c, ok in status.authorization_criteria(
        [r], SOURCES, today=TODAY) if not ok]
    assert unmet, (
        f"breaking {field!r} refuses at the gate but leaves all criteria met. "
        "The milestone reads AUTHORIZED over a live refusal.")


def test_criteria_with_no_consumer_are_named(requirements):
    """A gate must not close on a threshold nobody compares against."""
    report = status.enforcement_report()
    decorative = set(report["not_implemented"])
    assert "spatial:maximum_axis_ratio_error" in decorative
    assert "spatial:minimum_objects_per_species" in decorative
    assert "spatial:require_independent_validation_sample" in decorative
    assert "spatial:maximum_biovolume_fraction_error" in set(report["hard"])
    assert "spatial:maximum_porosity_error" in set(report["derived"])


@pytest.mark.parametrize("kwargs,expect_number,must_stay_met", [
    ({"approval_source_id": ""},          8, (7,)),  # provenance, NOT dates
    ({"scope_hash": " "},                 2, (7,)),  # binding, NOT dates
    ({"approval_effective_date": ""},     7, (8,)),  # an actual date
    ({"approval_expiration_date": ""},    7, (8,)),
    ({"culturing_start_date": ""},        7, (8,)),
    # AN IDENTIFIER MUST NOT BE ABLE TO MANUFACTURE A DATE BLOCKER. These
    # values are unresolvable source ids and nothing more; the refusal echoes
    # them, and a substring match on the concatenated text made criterion 7
    # fire on the word inside somebody's identifier.
    ({"approval_source_id": "expired"},   8, (7,)),
    ({"approval_source_id": "not-an-ISO-date"}, 8, (7,)),
    # AND NO VALUE MAY IMPERSONATE A CROSS-FIELD SENTENCE. The relation
    # phrases are checked only AFTER the field is identified, so an identifier
    # spelling one out is still just an unresolvable identifier.
    ({"approval_source_id": "does not match the conditions"}, 8, (6, 7)),
    ({"approval_source_id": "is not the target system"}, 8, (9, 7)),
    ({"strain_identities": "unknown"},    1, (7, 8)),
    ({"containment_facility": "TBD"},     3, (7, 8)),
    ({"is_target_system": "false"},       9, (7, 8)),
])
def test_an_unset_field_blames_the_criterion_it_actually_belongs_to(
        kwargs, expect_number, must_stay_met):
    """BEING TOLD THE WRONG THING IS WRONG IS WORSE THAN BEING TOLD NOTHING.

    Criterion 7 matched a bare `"is unset"`, and `approval.problems` emits that
    phrase for the three dates AND for `approval_source_id` and
    `approval_scope_hash`. So an approval with no registered document reported
    criterion 7 — about dates — as its blocker, while criterion 8, about
    provenance, stayed green. Whoever read that milestone would have gone to
    fix a date that was never the problem.

    The blocker list is the whole output of this function. If it names the
    wrong criterion it is not a partial answer, it is a misdirection.
    """
    from test_approval import SOURCES, TODAY, row

    unmet = [c for c, ok in status.authorization_criteria(
        [row(**kwargs)], SOURCES, today=TODAY) if not ok]
    assert any(c.startswith(f"{expect_number}.") for c in unmet), (
        f"expected criterion {expect_number} to be the blocker, got: {unmet}")
    # and nothing fell through to the unmapped catch, which would mean the
    # criteria still cannot see this refusal at all
    assert not any(c.startswith("UNMAPPED") for c in unmet), unmet

    # THE HALF THAT ACTUALLY CATCHES THE BUG. Asserting the right criterion is
    # present does NOT detect a misdirection -- with the bare `"is unset"`
    # pattern restored, criterion 8 still fires for an unset source id and this
    # test passed while criterion 7 was also, wrongly, being reported. A blocker
    # list that names an extra criterion sends someone to fix a field that was
    # never the problem, so the criteria that are NOT implicated must stay met.
    for number in must_stay_met:
        assert not any(c.startswith(f"{number}.") for c in unmet), (
            f"criterion {number} is reported unmet, but {kwargs} has nothing to "
            f"do with it; blockers were: {unmet}")


def test_the_two_ways_a_scope_hash_can_fail_are_different_criteria():
    """ONE FIELD, TWO FAILURES, and they are not the same problem.

    An UNSET `approval_scope_hash` means nothing binds the approval to this row
    — the conditions are not frozen, criterion 2. A hash that DOES NOT MATCH
    means someone edited a growth condition after approval, so the approval no
    longer describes what it covers — criterion 6. Reporting either as the
    other sends a reader to the wrong remedy.
    """
    from test_approval import SOURCES, TODAY, row

    def blockers(r):
        return [c.split(".")[0] for c, ok in status.authorization_criteria(
            [r], SOURCES, today=TODAY) if not ok]

    edited = row()
    edited["temperature_C"] = "37"          # approved at 30
    assert blockers(edited) == ["6"]
    assert blockers(row(scope_hash=" ")) == ["2"]


@pytest.mark.parametrize("gid", [
    "GC1", "O'Brien-1", 'a"b', "GC 1", "R2A/30C",
    # THE ONE NO DELIMITER PAIR CAN BRACKET: repr must escape one of the two
    # quote characters, so every "strip the quoted prefix" regex stops early.
    'O\'Brien "lab"',
    'both \' and " and: a colon',
])
def test_the_condition_id_cannot_change_which_criterion_is_blamed(gid):
    """THE IDENTIFIER IS DATA. It must not steer the checklist.

    `approval.problems` formats the id with `!r`, and Python switches to double
    quotes when the value contains an apostrophe — so an ordinary condition id
    like `O'Brien-1` produced a prefix the classifier's regex could not strip.
    The head became the word "condition", nothing matched, and an unset source
    id was reported as a scope failure plus an UNMAPPED refusal.

    Four fixes narrowed that defect without closing it, because the defect was
    never the regex — it was asking a sentence built from unrestricted CSV data
    to say which field it was about. `approval.classified` now carries the
    subject, so these ids cannot reach the decision at all.

    Whether a growth condition happens to be named after someone Irish is not a
    fact about its biosafety approval.
    """
    from test_approval import SOURCES, TODAY, row

    def blockers(**kw):
        return [c.split(".")[0] for c, ok in status.authorization_criteria(
            [row(growth_condition_id=gid, **kw)], SOURCES, today=TODAY)
            if not ok]

    assert blockers(approval_source_id="") == ["8"]
    assert blockers(approval_effective_date="") == ["7"]
    assert blockers(strain_identities="unknown") == ["1"]
    assert blockers() == []
