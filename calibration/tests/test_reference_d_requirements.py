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
    blocked = [r for r in requirements if r["status"] == "blocked_on_model"]
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
