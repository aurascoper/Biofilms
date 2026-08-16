"""The institutional approval gate: everything a presence check let through.

`D-APPROVAL` was once asserted as the words "d approved". The refusal was human;
the code was a presence check, and since the schema already refuses a blank it
was very nearly a no-op — any non-empty string cleared it.

Every test here is a string or a state that USED TO PASS. None of them is
hypothetical.

WHAT THIS SUITE DOES NOT CLAIM. None of it verifies an approval is genuine. It
verifies that the fields a reviewer would ask for are present, resolve, agree
with each other, and still describe this row. A determined forger defeats all of
it; a careless edit does not, and the careless edit is the realistic failure.
"""

from __future__ import annotations

from datetime import date

import pytest

from biofilm_calibration import approval
from biofilm_calibration.approval import scope_digest

TODAY = date(2026, 8, 16)

SOURCES = [
    {"source_id": "IBC_2026_0147",
     "document_type": "institutional_biosafety_approval",
     "url": "https://example.edu/ibc/0147", "sha256": ""},
    {"source_id": "IBC_NO_LOCATOR",
     "document_type": "institutional_biosafety_approval",
     "url": "", "sha256": ""},
    {"source_id": "NIST_WATER", "document_type": "declaration",
     "url": "https://nist.gov", "sha256": ""},
]


def row(**overrides) -> dict:
    """A complete, internally consistent approval. Every test breaks ONE thing."""
    r = dict(
        growth_condition_id="GC1",
        strain_identities="D. radiodurans R1; C. neoformans H99",
        biosafety_level_by_strain="DR:BSL1;CN:BSL2",
        containment_facility="Building 4, BSL-2 suite",
        risk_assessment_reference="RA-2026-0031",
        is_target_system="true",
        institutional_approval_id="IBC-2026-0147",
        institutional_approval_authority="Institutional Biosafety Committee",
        approval_source_id="IBC_2026_0147",
        approval_effective_date="2026-06-01",
        approval_expiration_date="2027-06-01",
        culturing_start_date="2026-07-01",
        approved_protocol_version="v3",
        growth_medium="R2A", temperature_C="30", pH="7.0",
        oxygen_condition="aerobic", substrate_or_membrane="PES",
        biofilm_age_h="72", flow_condition="static", irradiation="none",
    )
    scope_override = overrides.pop("scope_hash", None)
    r.update(overrides)
    r["approval_scope_hash"] = scope_override or scope_digest(r)
    return r


def problems(r, sources=None, today=TODAY):
    return approval.problems([r], SOURCES if sources is None else sources,
                             today=today)


def test_a_complete_consistent_approval_is_accepted():
    """The suite would be worthless if nothing could pass it."""
    assert problems(row()) == []


# ------------------------------------------------------------- placeholders

@pytest.mark.parametrize("filler", [
    "d approved",       # the string that was actually offered
    "D Approved",       # case
    "  d approved  ",   # whitespace
    "approved", "pending", "TBD", "n/a", "none", "unknown", "placeholder",
    "submitted", "in review", "-",
])
def test_filler_text_is_refused_however_it_is_spelled(filler):
    """An approval identifier is free text by necessity — institutional formats
    differ, so constraining the identifier itself would be wrong. Naming the
    fillers is the only instrument left, and it is weaker than a vocabulary:
    it stops the careless case, not the determined one."""
    found = problems(row(institutional_approval_id=filler))
    assert any("names nothing" in p for p in found), found


def test_an_identifier_without_an_issuer_names_nothing():
    assert any("institutional_approval_authority" in p
               for p in problems(row(institutional_approval_authority="TBD")))


def test_the_strains_must_be_named_because_biosafety_follows_strains():
    """C. neoformans is BSL-2 where the others are BSL-1, so an approval that
    cannot say which organisms it covers covers none of them."""
    assert any("strain_identities" in p
               for p in problems(row(strain_identities="unknown")))


# ---------------------------------------------------------- the document

def test_an_approval_that_resolves_to_the_wrong_kind_of_document_is_refused():
    """The failure this closes: a source_id that is real, resolves, and is a
    cross-section table. Without document_type nothing in the data can see it."""
    found = problems(row(approval_source_id="NIST_WATER"))
    assert any("not an approval" in p for p in found), found


def test_an_unresolvable_source_is_refused():
    assert any("does not resolve" in p
               for p in problems(row(approval_source_id="IBC_NOPE")))


def test_an_approval_nobody_can_retrieve_is_not_evidence():
    found = problems(row(approval_source_id="IBC_NO_LOCATOR"))
    assert any("neither a locator nor a digest" in p for p in found), found


def test_an_identifier_with_no_registered_document_is_refused():
    assert any("cannot be checked against anything" in p
               for p in problems(row(approval_source_id="")))


# -------------------------------------------------------------- the dates

def test_culturing_may_not_precede_the_approval():
    """The requirement says the approval must be obtained BEFORE culturing.
    This is that sentence, enforced."""
    found = problems(row(culturing_start_date="2026-05-01"))
    assert any("must precede culturing" in p for p in found), found


def test_a_lapsed_approval_authorises_nothing():
    found = problems(row(approval_expiration_date="2026-08-01"))
    assert any("expired" in p for p in found), found


def test_expiry_is_judged_against_an_injected_date():
    """The verdict genuinely depends on when it is asked, which is correct for
    an approval — and means a test that did not pin the date would pass today
    and fail later for no reason anyone could find."""
    r = row(approval_expiration_date="2026-09-01")
    assert problems(r, today=date(2026, 8, 16)) == []
    assert any("expired" in p for p in problems(r, today=date(2026, 9, 2)))


def test_culturing_after_expiry_is_refused():
    found = problems(row(culturing_start_date="2027-07-01"))
    assert any("after the approval expired" in p for p in found), found


def test_an_approval_expiring_before_it_starts_is_refused():
    found = problems(row(approval_effective_date="2026-06-01",
                         approval_expiration_date="2026-05-01"))
    assert any("before it becomes effective" in p for p in found), found


@pytest.mark.parametrize("bad", ["June 2026", "01/06/2026", "2026-13-01", "soon"])
def test_dates_must_be_iso(bad):
    assert any("is not an ISO date" in p
               for p in problems(row(approval_effective_date=bad)))


# -------------------------------------------------------------- the scope

def test_editing_a_growth_condition_after_approval_invalidates_it():
    """THE REALISTIC FAILURE. Not forgery — someone changes the temperature and
    the approval silently comes to cover something it never covered."""
    r = row()
    r["temperature_C"] = "37"
    found = approval.problems([r], SOURCES, today=TODAY)
    assert any("no longer covers it" in p for p in found), found


def test_the_digest_covers_the_scope_and_not_the_bookkeeping():
    """A committee never saw `notes`, so changing it must not invalidate an
    approval — otherwise the check cries wolf and gets recomputed reflexively."""
    r = row()
    before = scope_digest(r)
    r["notes"] = "typo fixed"
    r["source_id"] = "SOMETHING_ELSE"
    assert scope_digest(r) == before


def test_the_digest_is_stable_under_reflowing_but_not_under_meaning():
    r = row()
    before = scope_digest(r)
    r["growth_medium"] = "  R2A  "
    assert scope_digest(r) == before, "whitespace is not meaning"
    r["growth_medium"] = "LB"
    assert scope_digest(r) != before, "medium is meaning"


def test_an_unset_scope_hash_binds_nothing():
    assert any("nothing binds this row" in p
               for p in problems(row(scope_hash=" ")))


# ------------------------------------------------------------ the surrogate

def test_a_surrogate_cannot_clear_the_target_gate():
    """A surrogate exercises the harness and can never clear D-APPROVAL, whether
    or not its own paperwork is in order."""
    found = problems(row(is_target_system="false"))
    assert any("not the target system" in p for p in found), found


# ----------------------------------------------------- the authoritative file

def test_the_gate_reads_only_the_authoritative_table(tmp_path):
    """The proposal file must never be a second, softer path to the same
    verdict. `approval.problems` takes ROWS, so it cannot be pointed at a
    different file by accident — the only caller is the spatial gate, which
    names baseline_condition.csv explicitly."""
    import inspect

    from biofilm_calibration.spatial import report

    src = inspect.getsource(report.evaluate)
    assert "baseline_condition.csv" in src
    assert "proposal" not in src, (
        "the gate must not read the non-authoritative proposal table")
