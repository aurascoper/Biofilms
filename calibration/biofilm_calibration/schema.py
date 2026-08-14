"""Shared CSV table contract for the calibration branches.

Both the spatial and material workstreams need the same thing: a table whose
columns are declared, whose vocabularies are closed, whose key is unique, and
whose numeric fields either parse or are explicitly empty. This is that, as a
library rather than as a test, so the harnesses can refuse bad input at read
time instead of discovering it downstream.

Two rules are baked in rather than left to each caller:

EMPTY IS NOT ZERO. A blank numeric cell means "not measured" and reads back as
None. Nothing here ever substitutes 0.0 for a missing measurement — that is the
single most dangerous default in a calibration ledger, because zero density and
unmeasured density are both plausible-looking numbers.

READY REQUIRES ATTRIBUTION. A row whose status is `ready` must carry a
source_id, and `unresolved`/`blocked` rows must not carry a value. The
provenance vocabulary is shared here so the two branches cannot drift apart on
what "ready" means.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

# Shared provenance vocabulary. Kept here, not duplicated per branch, so that
# "ready" means the same thing in the spatial and material ledgers.
# `unsupported_by_current_model` is deliberately distinct from `blocked`. A
# blocked quantity awaits a measurement; an unsupported one awaits a MODEL
# CHANGE, and no quantity of data will clear it. Filtering on status must be
# able to separate "go and measure it" from "the simulation cannot represent
# this yet".
STATUS = frozenset({"ready", "provisional", "blocked", "unresolved",
                    "not_applicable", "unsupported_by_current_model"})

# Statuses that assert a quantity has no value. Kept explicit rather than
# inferred, so a new status cannot silently fall through the check.
NO_VALUE_STATUSES = frozenset({"blocked", "unresolved",
                               "unsupported_by_current_model"})
EVIDENCE_BASIS = frozenset({
    "direct_measurement", "manufacturer_datasheet", "primary_literature",
    "derived", "declared", "proxy", "synthetic", "unresolved",
})


class SchemaError(ValueError):
    """Raised with every problem in a table, collected into one message."""


@dataclass(frozen=True)
class Column:
    name: str
    numeric: bool = False
    vocabulary: frozenset | None = None
    required: bool = True          # must be non-empty in every row
    doc: str = ""


@dataclass(frozen=True)
class TableSchema:
    name: str
    columns: tuple[Column, ...]
    key: tuple[str, ...] = ()      # uniqueness key; empty = no uniqueness rule
    doc: str = ""
    comments: tuple[str, ...] = field(default_factory=tuple)

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(c.name for c in self.columns)

    def column(self, name: str) -> Column:
        for c in self.columns:
            if c.name == name:
                return c
        raise KeyError(name)


def _parse_numeric(raw: str):
    """Blank -> None. Never zero: an unmeasured quantity is not a measured
    zero, and silently conflating them is how a calibration ledger lies."""
    text = raw.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{text!r} is not a number") from None


def read_table(path, schema: TableSchema) -> list[dict]:
    """Read and validate a CSV against `schema`. Collects every problem and
    raises once, in the repository's fail-loud style."""
    path = Path(path)
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(l for l in fh if not l.startswith("#"))
        header = tuple(reader.fieldnames or ())
        rows = list(reader)

    problems: list[str] = []
    if header != schema.names:
        missing = [c for c in schema.names if c not in header]
        extra = [c for c in header if c not in schema.names]
        if missing:
            problems.append(f"missing column(s): {', '.join(missing)}")
        if extra:
            problems.append(f"unknown column(s): {', '.join(extra)}")
        if not missing and not extra:
            problems.append("columns are out of declared order")

    parsed: list[dict] = []
    for i, row in enumerate(rows, start=2):          # header is line 1
        out = dict(row)
        for col in schema.columns:
            if col.name not in row:
                continue
            raw = (row[col.name] or "")
            if col.required and not raw.strip():
                problems.append(f"row {i}: {col.name} is empty (required)")
            if col.vocabulary is not None and raw.strip() \
                    and raw.strip() not in col.vocabulary:
                problems.append(
                    f"row {i}: {col.name} = {raw.strip()!r} not in "
                    f"{sorted(col.vocabulary)}")
            if col.numeric:
                try:
                    out[col.name] = _parse_numeric(raw)
                except ValueError as e:
                    problems.append(f"row {i}: {col.name} {e}")
        parsed.append(out)

    if schema.key:
        seen: dict[tuple, int] = {}
        for i, row in enumerate(parsed, start=2):
            k = tuple(str(row.get(c, "")) for c in schema.key)
            if k in seen:
                problems.append(
                    f"row {i}: duplicate {'+'.join(schema.key)} = {k} "
                    f"(first seen at row {seen[k]})")
            else:
                seen[k] = i

    problems.extend(_provenance_problems(parsed, schema))

    if problems:
        raise SchemaError(
            f"{path.name} does not satisfy the {schema.name} schema:\n  - "
            + "\n  - ".join(problems))
    return parsed


def _provenance_problems(rows: list[dict], schema: TableSchema) -> list[str]:
    """`ready` must be attributed; `blocked`/`unresolved` must not carry a
    value. Applied only to tables that declare the relevant columns."""
    names = schema.names
    if "status" not in names:
        return []
    problems = []
    value_cols = [c.name for c in schema.columns
                  if c.numeric and c.name not in schema.key]
    for i, row in enumerate(rows, start=2):
        status = (row.get("status") or "").strip()
        if status == "ready" and "source_id" in names \
                and not (row.get("source_id") or "").strip():
            problems.append(f"row {i}: status=ready with no source_id")
        if status in NO_VALUE_STATUSES:
            carried = [c for c in value_cols if row.get(c) is not None]
            if carried:
                problems.append(
                    f"row {i}: status={status} but carries a value in "
                    f"{', '.join(carried)} — a blocked quantity has no value")
    return problems


def write_template(path, schema: TableSchema) -> None:
    """Write a header-only CSV with the schema's comment block.

    Header-only is the honest state for a table whose measurements do not
    exist. It is not a placeholder for invented rows.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {schema.name} — {schema.doc}"] if schema.doc else []
    lines += [f"# {c}" for c in schema.comments]
    if schema.key:
        lines.append(f"# UNIQUENESS KEY: {' + '.join(schema.key)}")
    lines.append("# Blank numeric cells mean NOT MEASURED. They are never read as zero.")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
        csv.writer(fh, lineterminator="\n").writerow(schema.names)
