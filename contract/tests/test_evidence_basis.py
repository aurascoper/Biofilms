"""The evidence vocabulary, enforced as two relations instead of one word.

A ledger VALUE and a manuscript CLAIM are backed by different kinds of thing,
and until 2026-09-06 the repository held four disagreeing copies of "the"
evidence vocabulary. These tests pin the two sets in physical_contract, read
every file that carries an evidence_basis column (the parameter ledger, the
species table, which had ZERO consumers, and the claims ledger, whose column no
test read), and hold the claims ledger's legacy spellings to an allowlist that
may not grow. It does not shrink when a row is superseded either, which this
said until 2026-09-12: legacy_pairs() reads every row and the ledger is
additive, so a superseding row leaves the old pair in place. An entry goes only
when its stored pair leaves the ledger.

CONTROLS. Each guard below has a planted failure that must fire and a clean
case that must pass. The relation control is over CONSUMERS as well as sets:
the claims ledger fed through the value check must reject something, and the
parameter ledger fed through the claim check must reject something. Both fail
on ZERO hits, so a superseded row cannot make either pass vacuously.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from physical_contract import (CLAIM_EVIDENCE_BASIS, EVIDENCE_NULLS,
                               LEGACY_CLAIM_EVIDENCE, PARAMETER_EVIDENCE_BASIS,
                               PROVENANCE_SOURCES, canonical_claim_evidence,
                               claim_evidence_problems,
                               parameter_evidence_problems)

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
CLAIMS = DATA / "claims_ledger.csv"
PARAMETERS = DATA / "parameter_provenance.csv"
SPECIES = DATA / "species_parameter_provenance.csv"
ALLOWLIST = Path(__file__).resolve().parent / "fixtures" / "claims_ledger_legacy_evidence.csv"

VALUE_ONLY = PARAMETER_EVIDENCE_BASIS - CLAIM_EVIDENCE_BASIS
CLAIM_ONLY = CLAIM_EVIDENCE_BASIS - PARAMETER_EVIDENCE_BASIS


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    assert rows, f"{path.name} read back empty; nothing below would be checked"
    return rows


def legacy_pairs(rows) -> set:
    """(claim_id, stored_value) for every row still carrying a legacy value."""
    return {(r["claim_id"], r["evidence_basis"]) for r in rows
            if r["evidence_basis"] in LEGACY_CLAIM_EVIDENCE}


# --- the sets ------------------------------------------------------------

def test_the_two_relations_overlap_by_design_and_differ_by_design():
    assert PARAMETER_EVIDENCE_BASIS & CLAIM_EVIDENCE_BASIS == {
        "primary_literature", "derived", "declared", "synthetic"}
    assert "direct_measurement" in VALUE_ONLY   # a claim is not measured
    assert "prior_search" in CLAIM_ONLY         # a value is not searched for
    assert "computational" in CLAIM_ONLY
    # The source axis and the nulls are neither relation's vocabulary.
    assert not (PROVENANCE_SOURCES & (PARAMETER_EVIDENCE_BASIS | CLAIM_EVIDENCE_BASIS))
    assert not (EVIDENCE_NULLS & (PARAMETER_EVIDENCE_BASIS | CLAIM_EVIDENCE_BASIS))
    # `unresolved` is a status; it must not creep back in as a basis.
    assert "unresolved" not in PARAMETER_EVIDENCE_BASIS


# The calibration schema must read PARAMETER_EVIDENCE_BASIS as the same object.
# That assertion lives in calibration/tests/test_evidence_vocabulary.py, because
# CI's shared-contract tier installs contract and coupling only, and a test that
# cannot import its subject would skip or error rather than check.


# --- the canonicaliser ---------------------------------------------------

def test_canonicaliser_maps_legacy_reads_canonical_and_refuses_the_rest():
    assert canonical_claim_evidence("code") == ("computational", "code")
    assert canonical_claim_evidence("literature") == ("primary_literature", None)
    assert canonical_claim_evidence("none") == ("absent", None)
    assert canonical_claim_evidence("declared") == ("declared", None)
    assert canonical_claim_evidence("") == ("", None)
    for stored, (basis, source) in LEGACY_CLAIM_EVIDENCE.items():
        assert basis in CLAIM_EVIDENCE_BASIS | EVIDENCE_NULLS, stored
        assert source is None or source in PROVENANCE_SOURCES, stored
    with pytest.raises(ValueError, match="neither"):
        canonical_claim_evidence("measured")
    # A legacy spelling is readable but not writable.
    assert claim_evidence_problems("code") == [
        "'code' is a legacy spelling; new rows write basis 'computational' with source 'code'"]
    assert claim_evidence_problems("computational") == []


# --- the files -----------------------------------------------------------

def test_parameter_ledgers_use_the_value_vocabulary():
    for path in (PARAMETERS, SPECIES):
        bad = [(r.get("config_key") or r.get("claim_id"), r["evidence_basis"])
               for r in _read(path) if parameter_evidence_problems(r["evidence_basis"])]
        assert not bad, f"{path.name}: {bad}"


def test_claims_ledger_rows_are_canonical_or_allowlisted():
    rows = _read(CLAIMS)
    for r in rows:
        canonical_claim_evidence(r["evidence_basis"])   # raises on a new word
    have = legacy_pairs(rows)
    listed = {(r["claim_id"], r["stored_value"]) for r in _read(ALLOWLIST)}
    grew = have - listed
    assert not grew, f"legacy evidence values on rows not in the allowlist: {sorted(grew)}"
    stale = listed - have
    assert not stale, f"allowlist entries no longer carried by the ledger: {sorted(stale)}"


def test_allowlist_control_fires_in_both_directions():
    rows = _read(CLAIMS)
    listed = {(r["claim_id"], r["stored_value"]) for r in _read(ALLOWLIST)}
    assert legacy_pairs(rows) == listed          # clean case

    grown = rows + [{"claim_id": "NEW-01", "evidence_basis": "code"}]
    assert legacy_pairs(grown) - listed == {("NEW-01", "code")}

    drifted = [dict(r) for r in rows]
    target = next(r for r in drifted if r["evidence_basis"] == "code")
    target["evidence_basis"] = "simulation_output"     # same row, other legacy value
    assert (target["claim_id"], "code") in listed
    assert (target["claim_id"], "code") not in legacy_pairs(drifted)       # shrank
    assert (target["claim_id"], "simulation_output") not in listed         # and grew


# --- the relation, over consumers ----------------------------------------

def test_relation_control_over_planted_rows():
    claim_row = {"claim_id": "FIXTURE-CLAIM", "evidence_basis": "prior_search"}
    value_row = {"config_key": "fixture.value", "evidence_basis": "direct_measurement"}
    assert claim_evidence_problems(claim_row["evidence_basis"]) == []
    assert parameter_evidence_problems(claim_row["evidence_basis"]) == [
        "'prior_search' is a claim's evidence basis, not a value's"]
    assert parameter_evidence_problems(value_row["evidence_basis"]) == []
    assert claim_evidence_problems(value_row["evidence_basis"]) == [
        "'direct_measurement' is a value's evidence basis, not a claim's"]


def test_relation_control_over_live_data_fails_on_zero_hits():
    # Claims through the VALUE check: canonicalise first, so the rejection is
    # about the relation and not about a legacy spelling.
    # `absent` is a claim-side null the value check rightly refuses; it is
    # excluded here so the hits are about the relation, not the nulls.
    bases = {canonical_claim_evidence(r["evidence_basis"])[0] for r in _read(CLAIMS)}
    claim_hits = {b for b in bases - EVIDENCE_NULLS if parameter_evidence_problems(b)}
    assert claim_hits, "no claims-ledger row is rejected by the value check; the control is vacuous"
    assert claim_hits <= CLAIM_ONLY, claim_hits

    value_hits = {r["evidence_basis"] for r in _read(PARAMETERS)
                  if claim_evidence_problems(r["evidence_basis"])}
    assert value_hits, "no parameter-ledger row is rejected by the claim check; the control is vacuous"
    assert value_hits <= VALUE_ONLY, value_hits
