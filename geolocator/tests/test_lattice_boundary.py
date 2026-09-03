"""
The /api/lattice redaction boundary, tested adversarially.

The point of these tests is not that the boundary works on today's upstream
payload -- it demonstrably does. The point is that it keeps working when
upstream changes underneath it. So the mocked lattice here deliberately ships
fields that must never appear, including ones nobody has written yet.
"""

import importlib
import json
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from geolocator import lattice as lat  # noqa: E402
from geolocator.api import app  # noqa: E402

client = TestClient(app)


# A payload shaped like the real one, salted with things that must not escape.
HOSTILE_UPSTREAM = {
    "now": 1787000000000,
    "state_mtime_ms": 1787000000000,
    "state": {
        "generated": "2099-01-01T00:00:00+00:00",  # far future -> nominal
        "mission": {
            "phase": "MANAGE",
            "phase_reason": "drawdown ok",
            "entries_allowed": True,
            "timestamp": "2099-01-01T00:00:00+00:00",
            # Deliberately synthetic. Real account figures must never be used
            # as fixtures: this file is the test that PROVES the boundary hides
            # financial state, so putting genuine numbers in it would disclose
            # exactly what the boundary exists to withhold.
            "equity": 1234.56,              # must never escape
            "peak_equity": 9876.54,         # must never escape
            "drawdown_pct": 87.5,           # must never escape
        },
        "buckets": {
            "USDT": {"equity": 1111.11, "available": 222.22,
                     "margin_used": 888.88, "positions": [{"symbol": "BTC_USDT", "im": 44.44}]},
        },
        "cells": [
            {
                "symbol": "BTC_USDT", "quote": "USDT", "class": "crypto", "rank": 1,
                "has_position": True,
                "position_size": 0.42,      # must never escape
                "entry_px": 61234.5,        # must never escape
                "signal": {
                    "direction": "LONG", "score": 2.5, "z_4h": -1.9, "price": 61000.0,
                    "liquidity": "1M", "asset": "crypto",
                    "api_key": "sk-live-DEADBEEF",   # must never escape
                    "leverage": 20,                  # must never escape
                },
            },
            {"symbol": "ETH_USDT", "quote": "USDT", "class": "crypto",
             "rank": 2, "signal": None},
            {"no_symbol": True},             # malformed -> dropped
        ],
        "econ": {"gridcoin": {"price_usd": 0.005, "mcap_usd": 2.7e6,
                              "vol24h_usd": 13.9, "chg24h_pct": -14.8,
                              "private_note": "leak"}},   # must never escape
        "secrets": {"hl_private_key": "0xdeadbeef"},      # must never escape
        "babysitter": {"model": "x", "parity": {}},
    },
}

# Two checks that mean different things. Conflating them produces a false
# positive the moment an allowlisted free-text VALUE happens to contain a
# forbidden WORD -- which is exactly what phase_reason="drawdown ok" does below,
# and why that string is deliberately left in the fixture.
FORBIDDEN_KEYS = [
    "equity", "available", "margin_used", "peak_equity", "drawdown_pct",
    "position_size", "entry_px", "api_key", "leverage", "secrets",
    "hl_private_key", "private_note", "babysitter", "buckets", "positions",
]
FORBIDDEN_VALUES = ["sk-live", "DEADBEEF", "0xdeadbeef", "SMUGGLED", "1234.56", "9876"]


def walk_keys(node):
    if isinstance(node, dict):
        for k, v in node.items():
            yield k
            yield from walk_keys(v)
    elif isinstance(node, list):
        for v in node:
            yield from walk_keys(v)


@pytest.fixture(autouse=True)
def _mock_upstream(monkeypatch):
    monkeypatch.setattr(lat, "_fetch_upstream", lambda: json.loads(json.dumps(HOSTILE_UPSTREAM)))
    lat._cache["at"] = None
    lat._cache["doc"] = None
    yield
    lat._cache["at"] = None
    lat._cache["doc"] = None


# ── the boundary holds under injection ───────────────────────────────────────
def test_no_forbidden_key_survives_the_boundary():
    keys = {k.lower() for k in walk_keys(client.get("/api/lattice").json())}
    leaked = sorted(keys & {k.lower() for k in FORBIDDEN_KEYS})
    assert not leaked, f"forbidden keys reached /api/lattice: {leaked}"


def test_no_secret_value_survives_the_boundary():
    body = client.get("/api/lattice").text
    leaked = [w for w in FORBIDDEN_VALUES if w.lower() in body.lower()]
    assert not leaked, f"secret values reached /api/lattice: {leaked}"


def test_free_text_is_capped_not_trusted():
    """phase_reason is upstream free text. It is bounded, not sanitised."""
    hostile = json.loads(json.dumps(HOSTILE_UPSTREAM))
    hostile["state"]["mission"]["phase_reason"] = "A" * 5000
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(lat, "_fetch_upstream", lambda: hostile)
        lat._cache["at"] = None
        reason = client.get("/api/lattice").json()["mission"]["phase_reason"]
    assert len(reason) == lat.FREE_TEXT_CAP + 1  # cap + ellipsis
    lat._cache["at"] = None


def test_only_allowlisted_keys_exist():
    d = client.get("/api/lattice").json()
    assert set(d["mission"]) <= set(lat.MISSION_FIELDS)
    for cell in d["cells"]:
        assert set(cell) <= set(lat.CELL_FIELDS) | {"signal", "has_position"}
        if cell["signal"] is not None:
            assert set(cell["signal"]) <= set(lat.SIGNAL_FIELDS)
    assert set(d["econ"]["gridcoin"]) <= set(lat.ECON_GRIDCOIN_FIELDS)


def test_malformed_cells_are_dropped_not_passed():
    d = client.get("/api/lattice").json()
    assert [c["symbol"] for c in d["cells"]] == ["BTC_USDT", "ETH_USDT"]


def test_allowlisted_key_cannot_smuggle_a_nested_object(monkeypatch):
    """A key allowlist with unbounded values is not a boundary."""
    hostile = json.loads(json.dumps(HOSTILE_UPSTREAM))
    hostile["state"]["cells"][0]["rank"] = {"nested": "sk-live-SMUGGLED"}
    hostile["state"]["mission"]["phase"] = {"nested": "sk-live-SMUGGLED"}
    monkeypatch.setattr(lat, "_fetch_upstream", lambda: hostile)
    lat._cache["at"] = None
    body = client.get("/api/lattice").text
    assert "SMUGGLED" not in body
    assert client.get("/api/lattice").json()["cells"][0]["rank"] is None


# ── the client cannot ask for more ───────────────────────────────────────────
@pytest.mark.parametrize("qs", ["", "?view=full", "?view=", "?view=stream",
                                "?view=stream&view=full", "?VIEW=full"])
def test_view_parameter_is_never_honoured(qs):
    baseline = client.get("/api/lattice").json()
    got = client.get(f"/api/lattice{qs}").json()
    baseline.pop("fetched_at", None), got.pop("fetched_at", None)
    assert got == baseline


def test_upstream_view_is_a_server_constant():
    assert lat.UPSTREAM_VIEW == "stream"


# ── position presence is opt-in ──────────────────────────────────────────────
def test_has_position_absent_by_default():
    d = client.get("/api/lattice").json()
    assert all("has_position" not in c for c in d["cells"])
    assert d["position_presence_exposed"] is False


def test_has_position_only_behind_the_explicit_flag(monkeypatch):
    monkeypatch.setattr(lat, "EXPOSE_POSITION_PRESENCE", True)
    lat._cache["at"] = None
    d = client.get("/api/lattice").json()
    assert d["cells"][0]["has_position"] is True
    assert d["cells"][1]["has_position"] is False


# ── the globe survives the board being down ──────────────────────────────────
def test_board_down_is_unavailable_not_a_crash(monkeypatch):
    import urllib.error

    def boom():
        raise urllib.error.URLError("Connection refused")

    monkeypatch.setattr(lat, "_fetch_upstream", boom)
    lat._cache["at"] = None
    r = client.get("/api/lattice")
    assert r.status_code == 200
    d = r.json()
    assert d["status"] == "unavailable" and d["cells"] == [] and d["error"]


def test_missing_state_is_unavailable(monkeypatch):
    monkeypatch.setattr(lat, "_fetch_upstream", lambda: {"state": None})
    lat._cache["at"] = None
    d = client.get("/api/lattice").json()
    assert d["status"] == "unavailable" and d["cells"] == []


# ── the allowlist is load-bearing for CA parity ──────────────────────────────
def test_allowlist_preserves_the_boards_ca_seed_fields():
    """The lattice panel reproduces the board's automaton bit-for-bit on any cell
    with no open position, because the board's 12-field seed string ends in five
    position fields that evaluate to "" when there is no position -- and the panel
    always joins five empty fields.

    Measured against the live board on 2026-08-26: 1,102 cells, 27 held,
    1,075 bit-identical, 0 unheld divergences.

    Dropping any field below silently breaks that parity, so it is pinned here.
    """
    assert {"symbol", "quote", "class", "rank"} <= set(lat.CELL_FIELDS)
    assert {"direction", "score", "z_4h"} <= set(lat.SIGNAL_FIELDS)
    # And the five that must NOT cross, or rule 30 would leak position side.
    assert not ({"side", "vol", "leverage", "entry", "liq"} & set(lat.SIGNAL_FIELDS))
    assert not ({"side", "vol", "leverage", "entry", "liq"} & set(lat.CELL_FIELDS))
