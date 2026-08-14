"""The shared table contract, including the two defaults that matter."""

from __future__ import annotations

import pytest

from biofilm_calibration.schema import (Column, SchemaError, TableSchema,
                                        read_table, write_template)

S = TableSchema(
    name="demo",
    columns=(
        Column("sample_id"),
        Column("value", numeric=True, required=False),
        Column("status", vocabulary=frozenset({"ready", "blocked"})),
        Column("source_id", required=False),
    ),
    key=("sample_id",),
    doc="test fixture",
)


def _write(tmp_path, rows, header=None):
    p = tmp_path / "t.csv"
    head = ",".join(header or S.names)
    p.write_text(head + "\n" + "\n".join(rows) + ("\n" if rows else ""))
    return p


def test_blank_numeric_is_none_never_zero(tmp_path):
    """An unmeasured quantity and a measured zero are not the same fact."""
    rows = read_table(_write(tmp_path, ["s1,,blocked,"]), S)
    assert rows[0]["value"] is None
    assert rows[0]["value"] != 0.0


def test_ready_requires_a_source(tmp_path):
    with pytest.raises(SchemaError, match="no source_id"):
        read_table(_write(tmp_path, ["s1,1.5,ready,"]), S)
    read_table(_write(tmp_path, ["s1,1.5,ready,SRC1"]), S)   # ok


def test_blocked_rows_must_not_carry_a_value(tmp_path):
    with pytest.raises(SchemaError, match="blocked quantity has no value"):
        read_table(_write(tmp_path, ["s1,0.05,blocked,"]), S)


def test_vocabulary_uniqueness_and_numeric_are_enforced(tmp_path):
    with pytest.raises(SchemaError, match="not in"):
        read_table(_write(tmp_path, ["s1,1.0,maybe,SRC1"]), S)
    with pytest.raises(SchemaError, match="duplicate sample_id"):
        read_table(_write(tmp_path, ["s1,1.0,ready,SRC1",
                                     "s1,2.0,ready,SRC1"]), S)
    with pytest.raises(SchemaError, match="not a number"):
        read_table(_write(tmp_path, ["s1,abc,ready,SRC1"]), S)


def test_unknown_and_missing_columns_are_named(tmp_path):
    with pytest.raises(SchemaError, match="unknown column"):
        read_table(_write(tmp_path, [], header=list(S.names) + ["extra"]), S)
    with pytest.raises(SchemaError, match="missing column"):
        read_table(_write(tmp_path, [], header=["sample_id"]), S)


def test_every_problem_is_reported_at_once(tmp_path):
    """Fail-loud, repository style: one error listing everything wrong."""
    with pytest.raises(SchemaError) as exc:
        read_table(_write(tmp_path, ["s1,abc,maybe,", "s1,,ready,"]), S)
    msg = str(exc.value)
    assert "not a number" in msg and "not in" in msg and "duplicate" in msg \
        and "no source_id" in msg


def test_template_round_trips_empty(tmp_path):
    p = tmp_path / "tmpl.csv"
    write_template(p, S)
    assert read_table(p, S) == []
    assert "NOT MEASURED" in p.read_text()
