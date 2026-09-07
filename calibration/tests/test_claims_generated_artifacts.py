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
    row = next(r for r in ledger._rows() if r["claim_id"] == "PILOT-LEV-01")
    assert row["document"] in ledger.GENERATED_ARTIFACTS
    phrase = ledger.distinguishing_phrase(row["claim_text"])
    assert phrase, "the selected ledger claim must be detectable"
    monkeypatch.setattr(ledger, "REPO", tmp_path)
    path = ledger.document_path(row["document"])
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"claim": row["claim_text"]}), encoding="utf-8")
    assert phrase.lower() in ledger.normalise_markup(path.read_text())
    with pytest.raises(AssertionError, match=row["claim_id"]):
        ledger.test_no_deleted_claim_survives_in_the_document_it_names([row])
    path.write_text(json.dumps({"claim": "withdrawn; see ledger"}), encoding="utf-8")
    ledger.test_no_deleted_claim_survives_in_the_document_it_names([row])
