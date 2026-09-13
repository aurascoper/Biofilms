import json
import math
import os
import sqlite3
from copy import deepcopy
from datetime import datetime, timedelta, timezone

import pytest

from geolocator.market import FUEL_TICKERS, MARKET_PROXY
from geolocator.research import record_bar, research_snapshot, validate_family

NOW = datetime(2026, 9, 13, tzinfo=timezone.utc)


def test_legacy_store_never_claims_historical_availability(tmp_path):
    path = tmp_path / "bars.sqlite"
    with sqlite3.connect(path) as con:
        con.execute("CREATE TABLE bars (symbol TEXT, ts INTEGER, close REAL)")
        con.execute("INSERT INTO bars VALUES('SECRET_ACCOUNT_BALANCE',1,9000)")
    before = path.read_bytes()
    doc = research_snapshot(path, NOW)
    assert (
        doc["status"] == "unavailable"
        and doc["reason"] == "availability_receipts_missing"
    )
    assert "SECRET" not in json.dumps(doc)
    assert path.read_bytes() == before
    os.utime(path, (NOW.timestamp(), NOW.timestamp()))
    assert research_snapshot(path, NOW) == doc
    assert (
        research_snapshot(tmp_path / "missing", NOW)["reason"] == "market_store_absent"
    )
    assert not (tmp_path / "missing").exists()


def test_all_pairs_returned_and_suppression_detected(tmp_path):
    path = tmp_path / "bars.sqlite"
    members = [s for values in FUEL_TICKERS.values() for s in values] + [MARKET_PROXY]
    for s, symbol in enumerate(members):
        price = 100.0
        for i in range(33):
            at = NOW - timedelta(days=32 - i)
            price *= math.exp(0.003 * math.sin(i * (s + 1) + 0.7 * s))
            record_bar(
                path,
                symbol=symbol,
                completed_at=at.isoformat(),
                close=str(price),
                source_url="https://example.test/synthetic",
                now=at + timedelta(seconds=1),
            )
    doc = research_snapshot(path, NOW + timedelta(seconds=2))
    assert doc["status"] == "available"
    assert len(doc["correlations"]) == doc["familySize"] == doc["computedCount"] == 6
    assert len(doc["baskets"]["Gas"]["returns"]) == 32
    validate_family(doc)
    bad = deepcopy(doc)
    bad["correlations"].pop()
    with pytest.raises(ValueError, match="incomplete_correlation_family"):
        validate_family(bad)
    bad = deepcopy(doc)
    bad["computedCount"] = 5
    with pytest.raises(ValueError, match="count_mismatch"):
        validate_family(bad)
    # Financial redaction holds by construction: no source accepts account or
    # lattice fields, including nested fields or ranks with unknown lineage.
    assert not {
        "balance",
        "positions",
        "equity",
        "cells",
        "signal",
        "has_position",
    } & set(doc)
    assert doc["latticeEvidence"]["reason"] == "strategy_lineage_not_declared"
    assert (
        research_snapshot(path, NOW + timedelta(days=4))["reason"]
        == "completed_bars_stale"
    )


def test_future_and_backfilled_bars_have_explicit_availability(tmp_path):
    path = tmp_path / "bars.sqlite"
    with pytest.raises(ValueError, match="bar_not_completed"):
        record_bar(
            path,
            symbol="NGAS_USDT",
            completed_at=(NOW + timedelta(seconds=1)).isoformat(),
            close="10",
            source_url="https://example.test",
            now=NOW,
        )
    past = NOW - timedelta(days=30)
    record_bar(
        path,
        symbol="NGAS_USDT",
        completed_at=past.isoformat(),
        close="10",
        source_url="https://example.test",
        now=NOW,
    )
    with sqlite3.connect(path) as con:
        assert (
            con.execute("SELECT available_at FROM research_bars").fetchone()[0]
            == NOW.isoformat()
        )
    assert research_snapshot(path, past)["reason"] == "future_or_noncausal_bar"
