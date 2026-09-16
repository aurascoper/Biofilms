"""Controls for the two named pilot outputs exempted from existence checks.

The exemption is an exact path list, not a directory-wide waiver. These controls
call the same guards as the calibration suite, including the content guard when
a generated output exists.
"""

import json

import pytest

import test_claims_ledger as ledger


@pytest.mark.parametrize("document", (
    "docs/missing-source.md",
    "artifacts/pilot/undeclared-output.json",
))
@pytest.mark.parametrize("guard", (
    ledger.test_every_deleted_claim_names_a_document_the_guard_can_read,
    ledger.test_every_claim_names_a_document_that_exists,
))
def test_missing_undeclared_documents_still_fail(document, guard, tmp_path, monkeypatch):
    monkeypatch.setattr(ledger, "REPO", tmp_path)
    row = dict(ledger._rows()[0], document=document, status="delete")
    assert not ledger.document_path(document).exists()
    with pytest.raises(AssertionError, match=document):
        guard([row])


def test_only_declared_pilot_outputs_may_be_absent(tmp_path, monkeypatch, capsys):
    rows = [r for r in ledger._rows() if r["document"] in ledger.GENERATED_ARTIFACTS]
    assert {r["document"] for r in rows} == ledger.GENERATED_ARTIFACTS
    monkeypatch.setattr(ledger, "REPO", tmp_path)
    for path in ledger.GENERATED_ARTIFACTS:
        assert not ledger.document_path(path).exists()
    ledger.test_every_deleted_claim_names_a_document_the_guard_can_read(rows)
    ledger.test_every_claim_names_a_document_that_exists(rows)
    # An absent output is uncovered content, even though its declared absence
    # is legitimate. Keep the affected ledger ids visible in pytest's output.
    with capsys.disabled():
        for row in rows:
            print(f"\n  generated output absent: {row['claim_id']} [{row['status']}] "
                  f"in {row['document']}; content not checked")


def test_present_generated_output_still_rejects_a_withdrawn_claim(tmp_path, monkeypatch):
    """THE CONTROL IS THE ARTIFACT'S REAL DELETED FIELD, NOT LEDGER PROSE.

    This test used to write `{"claim": <the ledger's claim_text>}` and pass. It
    proved that substring matching works, which was never in doubt, and it could
    not prove the thing it exists for. The producer never wrote that sentence:
    `coupling/scripts/openmc_nested_pilot.py` emitted
    `material_lever_sensitivity_rel_l2` into `budget_doc` as six hardcoded
    numbers, and the regression this guards against is a regenerated pilot
    putting them back. So rebuild THAT shape and run the guard on it.
    """
    row = next(r for r in ledger._rows() if r["claim_id"] == "PILOT-LEV-01")
    assert row["document"] in ledger.GENERATED_ARTIFACTS
    field = ledger.deleted_field(row)
    assert field == "material_lever_sensitivity_rel_l2", (
        "the ledger row no longer declares the removed field in `location`; "
        "without it this control is back to matching prose")

    monkeypatch.setattr(ledger, "REPO", tmp_path)
    path = ledger.document_path(row["document"])
    path.parent.mkdir(parents=True)

    # `budget_doc` as the producer wrote it before the field was dropped. The
    # six values are PILOT-LEV-01's own `reported_value`; the key is what the
    # guard keys on, and the values are here so the fixture is a real artifact
    # rather than a stub carrying one string.
    regressed = {
        "schema_version": 1, "tier": "S0", "target_calibration": False,
        "openmc_runs": 12, "histories": 200000,
        field: {"density_x1.35": 0.0137, "Fe_5pct": 0.0033, "Gd_5pct": 0.0376,
                "dehydration": 0.0459, "Gd_20pct": 0.1241, "Gd_40pct": 0.2242},
        "noise_floor_basis": "decorrelated_seeds_identical_material",
    }
    path.write_text(json.dumps(regressed, indent=2), encoding="utf-8")
    with pytest.raises(AssertionError, match=row["claim_id"]):
        ledger.test_no_deleted_claim_survives_in_the_document_it_names([row])

    # The same artifact as the producer writes it now: every other field intact,
    # so a pass here means the field was detected and not the whole document.
    del regressed[field]
    path.write_text(json.dumps(regressed, indent=2), encoding="utf-8")
    ledger.test_no_deleted_claim_survives_in_the_document_it_names([row])


def test_prose_matching_alone_could_not_see_that_regression(tmp_path, monkeypatch):
    """WHY THE FIXTURE ABOVE STOPPED USING `claim_text`.

    Pinned rather than described: the phrase the prose guard would search for is
    ledger wording, and it does not occur in the artifact the producer writes.
    If phrase extraction ever changes so that it does, this test says so instead
    of the coverage quietly moving.
    """
    row = next(r for r in ledger._rows() if r["claim_id"] == "PILOT-LEV-01")
    phrase = ledger.distinguishing_phrase(row["claim_text"])
    assert phrase == "emitted into the budget artifact."

    artifact = json.dumps({"schema_version": 1,
                           ledger.deleted_field(row): {"Gd_40pct": 0.2242}})
    assert phrase.lower() not in ledger.normalise_markup(artifact)
    assert ledger.normalise_markup(ledger.deleted_field(row)) in \
        ledger.normalise_markup(artifact)
