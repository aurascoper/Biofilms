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
from typing import NamedTuple

from physical_contract import APPROVAL_DOCUMENT_TYPES, is_placeholder

# Fields whose filler would make the approval unidentifiable. `strain_identities`
# is here because biosafety follows strains: an approval that cannot say which
# organisms it covers covers none of them.
_MUST_NAME_SOMETHING = (
    "institutional_approval_id",
    "institutional_approval_authority",
    "containment_facility",
    "strain_identities",
    "biosafety_level_by_strain",
    "risk_assessment_reference",
    # WITHOUT THIS THE SCOPE DIGEST BINDS NOTHING USEFUL. `approved_protocol_version`
    # is one of the SCOPE_COLUMNS, so a blank value simply hashes as blank and the
    # digest still matches -- an approval could be recorded without identifying
    # WHICH protocol the institution reviewed, which is the binding the digest
    # exists to create.
    "approved_protocol_version",
)

def _entries(value) -> list[str]:
    return [p.strip() for p in str(value or "").split(";") if p.strip()]


def _match_key(value: str) -> str:
    """Normalise for comparison ONLY -- typography, never identity.

    Case and run-length of internal whitespace are how a hand-typed identifier
    varies between two cells that mean the same organism. Nothing here can turn
    one strain into another, which is the property that makes it safe: the
    comparison is looser about how a name is written and exactly as strict
    about which name it is.
    """
    return " ".join(str(value).split()).casefold()


def _biosafety_mapping_problems(row) -> list[str]:
    """Whether `biosafety_level_by_strain` is a per-strain MAPPING at all.

    ADDING THE FIELD TO `_MUST_NAME_SOMETHING` MADE IT REACHABLE, NOT VALID.
    That closed the fail-open on placeholder text and left a wider one: any
    non-filler string cleared the check. `"BSL2"` -- a single level for a
    mixed-BSL consortium, which is the precise error the schema warns about --
    passed. So did `"DR:BSL1"` against two declared strains, silently covering
    one organism and clearing the institutional milestone for both.

    THE CONVENTION IS THE PRODUCER'S, NOT INVENTED HERE. `acquisition.py`
    declares this field as keyed by the `strain_identities` entry VERBATIM.
    This checks that declaration and nothing beyond it.

    AND A COUNT IS NOT A BINDING. This check was structural once -- shape,
    count, distinctness -- and said so in a docstring that called key-to-strain
    binding impossible: nothing declared how `DR` related to `D. radiodurans
    R1`, so inferring it from initials would have been the consumer inventing
    semantics. True, and the wrong conclusion. `XX:BSL1;YY:BSL2` against two
    declared strains has the right shape, the right count and distinct keys, so
    it passed, and every authorization criterion then reported met with
    biosafety evidence for NEITHER declared organism. An approval that names no
    organism this row declares covers none of them, which is the same failure
    the count check was written to catch, one size larger.

    The fix was not to guess the convention but to make the producer state one,
    and the one it states removes the guesswork entirely: the key IS the
    identifier. Matching is case-insensitive and collapses whitespace --
    typography, not identity -- so a strain identifier containing ':' cannot be
    written as a key at all and is refused against `strain_identities` rather
    than misreported as a malformed pair.
    """
    value = row.get("biosafety_level_by_strain")
    if is_placeholder(value):
        return []          # already refused by the placeholder pass

    strains = _entries(row.get("strain_identities"))
    pairs = _entries(value)
    out: list[str] = []

    # REFUSED AGAINST THE FIELD THAT IS ACTUALLY WRONG. A strain identifier
    # carrying ':' cannot be a key, and letting it through would surface below
    # as "not a strain:level pair" -- pointing the reader at the mapping when
    # the mapping is the only correct thing about the row.
    unexpressible = [s for s in strains if ":" in s]
    if unexpressible:
        out.append(
            f"strain_identities entries {unexpressible} contain ':', so they "
            "cannot be used as biosafety_level_by_strain keys, which are the "
            "identifiers verbatim. Rename the strain or the approval cannot "
            "state a level for it")
        return out

    malformed = [p for p in pairs if p.count(":") != 1
                 or not p.split(":")[0].strip()
                 or not p.split(":")[1].strip()]
    if malformed:
        out.append(
            f"biosafety_level_by_strain has entries that are not "
            f"strain:level pairs: {malformed}. A single level for a "
            "mixed-BSL consortium is the error this field exists to prevent "
            "-- biosafety follows strains, not species")
        return out

    levels = [p.split(":")[1].strip().upper().replace("-", "")
              for p in pairs]
    unknown = [lv for lv in levels if lv not in _BIOSAFETY_LEVELS]
    if unknown:
        out.append(
            f"biosafety_level_by_strain names levels {unknown}, which are not "
            f"recognised biosafety levels {sorted(_BIOSAFETY_LEVELS)}")

    keys = [p.split(":")[0].strip() for p in pairs]
    if len(set(keys)) != len(keys):
        out.append(
            f"biosafety_level_by_strain repeats a strain key in {keys}; one "
            "entry per strain, or a level is silently overridden")

    # THE BINDING. This replaces a count comparison, which `XX:BSL1;YY:BSL2`
    # satisfied against two declared strains while naming neither of them. Sets,
    # not lengths -- and reported in both directions, because an uncovered
    # strain and a key naming nothing declared are different mistakes with the
    # same cause and different repairs.
    if strains:
        declared = {_match_key(s): s for s in strains}
        keyed = {_match_key(k): k for k in keys}
        uncovered = [declared[n] for n in declared if n not in keyed]
        unrecognised = [keyed[n] for n in keyed if n not in declared]
        if uncovered or unrecognised:
            detail = []
            if uncovered:
                detail.append(f"declares {uncovered} with no level")
            if unrecognised:
                detail.append(
                    f"names {unrecognised}, which strain_identities does not "
                    "declare")
            out.append(
                "biosafety_level_by_strain is keyed by the strain_identities "
                "entry verbatim, and this row " + " and ".join(detail) + ". An "
                "approval that omits a strain does not cover it, and one that "
                "names an organism this row never declared is evidence about "
                "something else -- neither omission may read as coverage")
    return out


# The levels an approval may name. Not a guess: BSL-1..4 is the standard
# containment scale, and anything outside it is a value nobody issuing an
# approval would write.
_BIOSAFETY_LEVELS = frozenset({"BSL1", "BSL2", "BSL3", "BSL4"})

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


class Refusal(NamedTuple):
    """One refusal, with the field it is ABOUT carried beside the prose.

    THE SUBJECT IS DATA, NOT SOMETHING TO RECOVER FROM THE SENTENCE. The
    institutional checklist used to classify these by parsing the message --
    strip the "growth condition <id>" prefix, take the first token as the field
    name -- and that arrangement failed four times running, each time because a
    VALUE reached a decision it had no part in: an echoed field value matching
    a criterion pattern; a longer echoed value impersonating a relation phrase;
    an apostrophe in a condition id flipping `repr` to double quotes; and then
    an id holding BOTH quote characters, which `repr` must escape and no
    delimiter pair can bracket.

    Condition ids are unrestricted CSV data. Prose built from them cannot be
    parsed back into its parts, so the parts are kept instead.

    `subject` is the column the refusal is about. It is None for the sentences
    describing a relation BETWEEN fields, which have no single subject, and
    `"approval_scope_hash:mismatch"` distinguishes an edited condition from an
    unset hash -- one field, two failures, two different remedies.
    """

    subject: str | None
    text: str


def classified(rows, sources=None, *,
               today: date | None = None) -> list[Refusal]:
    """Everything wrong with the declared approvals, each with its subject.

    `today` is INJECTED and defaults to the real date. An approval expires, so
    the verdict genuinely depends on when it is asked — that is correct, and it
    also means a test that does not pin the date is a time bomb. The caller
    passes the date it judged against and reports it.

    `sources` is the spatial branch registry, already read. Resolution happens
    against rows rather than a path so this module does no I/O and stays as
    testable as the rest of the gate.
    """
    out: list[Refusal] = []
    today = today or date.today()
    by_id = {r["source_id"]: r for r in (sources or [])}

    def add(subject, text):
        out.append(Refusal(subject, f"growth condition {gid!r}: {text}"))

    for row in rows:
        gid = row.get("growth_condition_id", "<unnamed>")

        for field in _MUST_NAME_SOMETHING:
            if is_placeholder(row.get(field)):
                add(field,
                    f"{field} = {(row.get(field) or '')!r} names nothing. An "
                    "approval is evidence produced by an institution; filler "
                    "text has no evidentiary force whatever it says")

        for text in _biosafety_mapping_problems(row):
            add("biosafety_level_by_strain", text)

        if row.get("is_target_system") != "true":
            add("is_target_system",
                "this is not the target system: a surrogate exercises the "
                "harness but cannot clear the gate for the seven-species "
                "consortium")

        # --- the approval document -------------------------------------
        source_id = (row.get("approval_source_id") or "").strip()
        if is_placeholder(source_id):
            add("approval_source_id",
                "approval_source_id is unset — an approval identifier with no "
                "registered document behind it cannot be checked against "
                "anything")
        elif source_id not in by_id:
            add("approval_source_id",
                f"approval_source_id {source_id!r} does not resolve in the "
                "spatial source registry")
        else:
            src = by_id[source_id]
            kind = (src.get("document_type") or "").strip()
            if kind not in APPROVAL_DOCUMENT_TYPES:
                add("approval_source_id",
                    f"approval_source_id {source_id!r} is document_type "
                    f"{kind or 'unset'!r}, not an approval — a source_id "
                    "resolving to the wrong KIND of document is how an "
                    "approval comes to cite a datasheet")
            if is_placeholder(src.get("url")) and is_placeholder(src.get("sha256")):
                add("approval_source_id",
                    f"approval document {source_id!r} has neither a locator "
                    "nor a digest — an approval nobody can retrieve is not "
                    "evidence")

        # --- the dates --------------------------------------------------
        parsed = {}
        for field in _DATE_FIELDS:
            raw = row.get(field)
            if is_placeholder(raw):
                add(field, f"{field} is unset")
                continue
            value = _parse(raw)
            if value is None:
                add(field, f"{field} = {raw!r} is not an ISO date (YYYY-MM-DD)")
            else:
                parsed[field] = value

        effective = parsed.get("approval_effective_date")
        expires = parsed.get("approval_expiration_date")
        started = parsed.get("culturing_start_date")

        # No single subject: each describes a relation BETWEEN dates.
        if effective and started and started < effective:
            add(None,
                f"culturing started {started} but the approval was not "
                f"effective until {effective} — the approval must precede "
                "culturing, which is the whole point of it")
        if expires and started and started >= expires:
            add(None, f"culturing started {started}, on or after the approval "
                      f"expired {expires}")
        if effective and expires and expires <= effective:
            add(None, f"the approval expires {expires}, on or before it "
                      f"becomes effective {effective}")
        if expires and today >= expires:
            add(None, f"the approval expired {expires} and today is {today} — "
                      "a lapsed approval authorises nothing")

        # --- the scope --------------------------------------------------
        declared = (row.get("approval_scope_hash") or "").strip()
        if is_placeholder(declared):
            add("approval_scope_hash",
                "approval_scope_hash is unset, so nothing binds this row's "
                "conditions to what was approved")
        else:
            computed = scope_digest(row)
            if declared.casefold() != computed.casefold():
                add("approval_scope_hash:mismatch",
                    "approval_scope_hash does not match the conditions on this "
                    f"row (declared {declared[:16]}…, computed "
                    f"{computed[:16]}…). A growth condition has been changed "
                    "since the approval was recorded, so the approval no "
                    "longer covers it")
    return out


def problems(rows, sources=None, *, today: date | None = None) -> list[str]:
    """The same refusals as prose, for display and for existing callers.

    Built FROM `classified` rather than beside it, so the two cannot disagree
    about what is wrong.
    """
    return [r.text for r in classified(rows, sources, today=today)]
