"""Institutional biosafety approval: the checks a presence test cannot make.

WHY THIS MODULE EXISTS. `D-APPROVAL` was asserted in conversation as the words
"d approved". The refusal was correct and it was human — the code was weaker than
the refusal. `institutional_approval_id` is already `required=True` at the
schema, so a BLANK was already refused; what passed was any non-empty string.
`approved`, `pending`, `TBD` and `d approved` all cleared the gate, and so would
a real-looking identifier belonging to another protocol, covering different
strains, at another facility, issued after culturing had started, or expired.

WHAT THIS CAN AND CANNOT DO, stated first so nothing downstream overstates it.

It can refuse filler text; require an issuing authority; require the approval to
resolve to a registered document of the right KIND with somewhere to find it;
require the dates to bracket the culturing and not to have lapsed; detect a
growth condition edited after approval; and keep a surrogate out of the target
gate.

It cannot verify that the document exists, that the committee issued it, or that
the identifier is real. No check local to this repository can. The value is that
it turns silent drift into a deliberate act and makes the external document the
anchor — not that it authenticates anything. An approval is evidence produced by
an institution, and the most software can do is refuse to pretend otherwise.

THE SCOPE DIGEST IS A DRIFT DETECTOR, NOT A FORGERY DETECTOR. Anyone who edits a
growth condition AND recomputes the hash defeats it. What it makes impossible is
editing a condition and silently inheriting an approval that never covered it —
which is the realistic failure, because it happens by accident.
"""

from __future__ import annotations

import hashlib
import json
from datetime import date

from physical_contract import APPROVAL_DOCUMENT_TYPES, is_placeholder

# Fields whose filler would make the approval unidentifiable. `strain_identities`
# is here because biosafety follows strains: an approval that cannot say which
# organisms it covers covers none of them.
_MUST_NAME_SOMETHING = (
    "institutional_approval_id",
    "institutional_approval_authority",
    "containment_facility",
    "strain_identities",
    "risk_assessment_reference",
    # WITHOUT THIS THE SCOPE DIGEST BINDS NOTHING USEFUL. `approved_protocol_version`
    # is one of the SCOPE_COLUMNS, so a blank value simply hashes as blank and the
    # digest still matches -- an approval could be recorded without identifying
    # WHICH protocol the institution reviewed, which is the binding the digest
    # exists to create.
    "approved_protocol_version",
)

# What the approval was issued AGAINST. Editing any of these after the fact
# changes what would have been reviewed, so the digest covers exactly them —
# not the whole row, which carries bookkeeping that no committee ever saw.
SCOPE_COLUMNS = (
    "strain_identities",
    "biosafety_level_by_strain",
    "containment_facility",
    "growth_medium",
    "temperature_C",
    "pH",
    "oxygen_condition",
    "substrate_or_membrane",
    "biofilm_age_h",
    "flow_condition",
    "irradiation",
    "approved_protocol_version",
)

_DATE_FIELDS = ("approval_effective_date", "approval_expiration_date",
                "culturing_start_date")


def scope_digest(row: dict) -> str:
    """SHA-256 over the approved scope, canonicalised.

    Canonical means sorted keys, collapsed whitespace and a stable separator, so
    a digest does not change because someone reflowed a cell. It is deliberately
    NOT case-folded: strain designations are case-significant.
    """
    scope = {k: " ".join(str(row.get(k) or "").split()) for k in SCOPE_COLUMNS}
    payload = json.dumps(scope, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _parse(value: str):
    try:
        return date.fromisoformat(" ".join(str(value or "").split()))
    except ValueError:
        return None


def problems(rows, sources=None, *, today: date | None = None) -> list[str]:
    """Everything wrong with the declared approvals, collected.

    `today` is INJECTED and defaults to the real date. An approval expires, so
    the verdict genuinely depends on when it is asked — that is correct, and it
    also means a test that does not pin the date is a time bomb. The caller
    passes the date it judged against and reports it.

    `sources` is the spatial branch registry, already read. Resolution happens
    against rows rather than a path so this module does no I/O and stays as
    testable as the rest of the gate.
    """
    out: list[str] = []
    today = today or date.today()
    by_id = {r["source_id"]: r for r in (sources or [])}

    for row in rows:
        gid = row.get("growth_condition_id", "<unnamed>")

        for field in _MUST_NAME_SOMETHING:
            if is_placeholder(row.get(field)):
                out.append(
                    f"growth condition {gid!r}: {field} = "
                    f"{(row.get(field) or '')!r} names nothing. An approval is "
                    "evidence produced by an institution; filler text has no "
                    "evidentiary force whatever it says")

        if row.get("is_target_system") != "true":
            out.append(
                f"growth condition {gid!r} is not the target system: a "
                "surrogate exercises the harness but cannot clear the gate for "
                "the seven-species consortium")

        # --- the approval document -------------------------------------
        source_id = (row.get("approval_source_id") or "").strip()
        if is_placeholder(source_id):
            out.append(
                f"growth condition {gid!r}: approval_source_id is unset — an "
                "approval identifier with no registered document behind it "
                "cannot be checked against anything")
        elif source_id not in by_id:
            out.append(
                f"growth condition {gid!r}: approval_source_id {source_id!r} "
                "does not resolve in the spatial source registry")
        else:
            src = by_id[source_id]
            kind = (src.get("document_type") or "").strip()
            if kind not in APPROVAL_DOCUMENT_TYPES:
                out.append(
                    f"growth condition {gid!r}: approval_source_id "
                    f"{source_id!r} is document_type {kind or 'unset'!r}, not "
                    "an approval — a source_id resolving to the wrong KIND of "
                    "document is how an approval comes to cite a datasheet")
            if is_placeholder(src.get("url")) and is_placeholder(src.get("sha256")):
                out.append(
                    f"growth condition {gid!r}: approval document "
                    f"{source_id!r} has neither a locator nor a digest — an "
                    "approval nobody can retrieve is not evidence")

        # --- the dates --------------------------------------------------
        parsed = {}
        for field in _DATE_FIELDS:
            raw = row.get(field)
            if is_placeholder(raw):
                out.append(f"growth condition {gid!r}: {field} is unset")
                continue
            value = _parse(raw)
            if value is None:
                out.append(
                    f"growth condition {gid!r}: {field} = {raw!r} is not an "
                    "ISO date (YYYY-MM-DD)")
            else:
                parsed[field] = value

        effective = parsed.get("approval_effective_date")
        expires = parsed.get("approval_expiration_date")
        started = parsed.get("culturing_start_date")

        if effective and started and started < effective:
            out.append(
                f"growth condition {gid!r}: culturing started {started} but the "
                f"approval was not effective until {effective} — the approval "
                "must precede culturing, which is the whole point of it")
        if expires and started and started >= expires:
            out.append(
                f"growth condition {gid!r}: culturing started {started}, on or "
                f"after the approval expired {expires}")
        if effective and expires and expires <= effective:
            out.append(
                f"growth condition {gid!r}: approval expires {expires}, on or "
                f"before it becomes effective {effective}")
        if expires and today >= expires:
            out.append(
                f"growth condition {gid!r}: the approval expired {expires} and "
                f"today is {today} — a lapsed approval authorises nothing")

        # --- the scope --------------------------------------------------
        declared = (row.get("approval_scope_hash") or "").strip()
        if is_placeholder(declared):
            out.append(
                f"growth condition {gid!r}: approval_scope_hash is unset, so "
                "nothing binds this row's conditions to what was approved")
        else:
            computed = scope_digest(row)
            if declared.casefold() != computed.casefold():
                out.append(
                    f"growth condition {gid!r}: approval_scope_hash does not "
                    f"match the conditions on this row (declared "
                    f"{declared[:16]}…, computed {computed[:16]}…). A growth "
                    "condition has been changed since the approval was "
                    "recorded, so the approval no longer covers it")
    return out
