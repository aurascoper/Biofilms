"""
Rev 6A agri_overlay: hash verification, and unknown-capacity handling in the two
generic endpoints (/api/stats, /api/plants) that assume capacity_mw is always a float.

capacity_mw=None is a real value for this layer (no NDVI baseline for a county), not an
edge case. The load-bearing tests here are the ones that would have caught the bug that
actually shipped once already: a payload with a mix of known and unknown capacity_mw
must not crash either endpoint, must exclude unknowns from any sum, and must report how
many were excluded rather than let that count disappear silently.

The hash tests are a positive/negative pair on purpose. Corrupting `cells` proves the
hash mechanism works at all; corrupting `generated_at` proves the widening from a
cells-only hash to a whole-payload hash actually did something, since a cells-only
scheme would have let that corruption straight through.
"""

import copy
import hashlib
import json
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geolocator.api import SOURCES, _load_agri_overlay, app  # noqa: E402
from geolocator.freshness import REFERENCE, TrackedSource  # noqa: E402

client = TestClient(app)


def _payload(cells):
    # deepcopy: callers mutate the returned payload's cells in place (to corrupt it after
    # hashing) -- without this, that mutation would leak into the shared MIXED_CELLS module
    # constant and pollute every later test in the session.
    p = {
        "schema_version": 2,
        "source_repo": "agri_yield_pipeline",
        "source_git_sha": "deadbeef",
        "generated_at": "2026-08-27T00:00:00+00:00",
        "county_set_note": "test fixture",
        "provenance": {},
        "cells": copy.deepcopy(cells),
    }
    p["payload_sha256"] = hashlib.sha256(
        json.dumps(p, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return p


def _write(path, payload):
    path.write_text(json.dumps(payload))


MIXED_CELLS = [
    {"id": "county:mo:boone", "kind": "county", "name": "Boone", "state": "MO",
     "lat": 38.9, "lon": -92.3,
     "ndvi": {"value": 0.4, "baseline_mu": 0.6, "baseline_sigma": 0.05, "z": -4.0,
               "n_baseline_years": 23, "as_of": "2026-04-20"},
     "weather": None, "yield_sensitivity": None, "sar_vv_db": None},
    {"id": "county:mo:stlouiscity", "kind": "county", "name": "St. Louis City", "state": "MO",
     "lat": 38.6, "lon": -90.2,
     "ndvi": None, "ndvi_unavailable_reason": "no_cropland",
     "weather": None, "yield_sensitivity": None, "sar_vv_db": None},
]


# ── loader-level: hash verification ───────────────────────────────────────────


def test_loader_accepts_a_valid_payload(tmp_path):
    f = tmp_path / "overlay.json"
    _write(f, _payload(MIXED_CELLS))
    result = _load_agri_overlay(f)
    assert len(result["items"]) == 2
    by_name = {i["name"]: i for i in result["items"]}
    assert by_name["Boone"]["capacity_mw"] == pytest.approx(4.0)
    # The whole point: unknown is None, not 0.0 or any other stand-in number.
    assert by_name["St. Louis City"]["capacity_mw"] is None
    assert by_name["St. Louis City"]["extra"]["ndvi_unavailable_reason"] == "no_cropland"


def test_corrupted_cells_is_rejected(tmp_path):
    """Positive control: the hash mechanism catches corruption at all."""
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS)
    payload["cells"][0]["ndvi"]["z"] = 999.0  # mutate after the hash was computed
    _write(f, payload)
    with pytest.raises(ValueError, match="corruption"):
        _load_agri_overlay(f)


def test_corrupted_generated_at_is_rejected(tmp_path):
    """Negative control for the payload-wide hash widening.

    This is the case a cells-only hash (schema_version 1's cells_sha256) would have
    missed entirely -- generated_at lives outside `cells`. If this test passes on a
    cells-only scheme it proves nothing about the fix; it has to fail before the
    widening and pass after."""
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS)
    payload["generated_at"] = "2099-01-01T00:00:00+00:00"  # mutate after hashing
    _write(f, payload)
    with pytest.raises(ValueError, match="corruption"):
        _load_agri_overlay(f)


def test_corrupted_source_git_sha_is_rejected(tmp_path):
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS)
    payload["source_git_sha"] = "0000000"
    _write(f, payload)
    with pytest.raises(ValueError, match="corruption"):
        _load_agri_overlay(f)


def test_schema_v1_cells_only_hash_still_verifies(tmp_path):
    """Backward compatibility: an old export with only cells_sha256 (no payload_sha256)
    still verifies against cells-only corruption, via the fallback path."""
    f = tmp_path / "overlay_v1.json"
    cells_hash = hashlib.sha256(
        json.dumps(MIXED_CELLS, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    payload = {
        "schema_version": 1, "source_repo": "agri_yield_pipeline", "source_git_sha": "x",
        "generated_at": "2026-08-27T00:00:00+00:00", "cells_sha256": cells_hash,
        "provenance": {}, "cells": MIXED_CELLS,
    }
    _write(f, payload)
    result = _load_agri_overlay(f)
    assert len(result["items"]) == 2


# ── endpoint-level: unknown capacity_mw handling ──────────────────────────────


@pytest.fixture
def agri_source(tmp_path, monkeypatch):
    f = tmp_path / "overlay.json"
    _write(f, _payload(MIXED_CELLS))
    src = TrackedSource(
        id="agri_overlay", path=f, layer_class=REFERENCE,
        loader=_load_agri_overlay, empty={"items": [], "meta": {}},
        count_of=lambda d: len((d or {}).get("items", [])),
    )
    monkeypatch.setitem(SOURCES, "agri_overlay", src)
    return src


def test_stats_excludes_unknown_from_sum_and_reports_the_count(agri_source):
    r = client.get("/api/stats", params={"layer": "agri_overlay"})
    assert r.status_code == 200
    body = r.json()
    assert body["count"] == 2
    assert body["unknown_capacity_count"] == 1
    # total_capacity_mw is Boone's 4.0 only -- St. Louis City's unknown must not have
    # been coerced into the sum as 0.0 (which would be indistinguishable from this
    # correct value here, but the point is it was *excluded*, not coincidentally equal).
    assert body["total_capacity_mw"] == pytest.approx(4.0)


def test_plants_default_filter_reports_zero_excluded(agri_source):
    """min_capacity defaults to 0.0 -- an unknown capacity_mw must not be excluded just
    because 0.0 is a technically-passable floor. excluded_unknown_capacity must still be
    present (0), not merely absent, so a consumer can tell 'no filter active' apart from
    'filter active, nothing excluded'."""
    r = client.get("/api/plants", params={"layer": "agri_overlay"})
    assert r.status_code == 200
    body = r.json()
    assert len(body["features"]) == 2
    assert body["excluded_unknown_capacity"] == 0


def test_plants_active_filter_excludes_unknown_and_reports_the_count(agri_source):
    """The case that actually recurs: an explicit min_capacity can't be verified against
    an unknown value, so it's excluded -- and unlike the bug being fixed, that exclusion
    is now counted, not silent."""
    r = client.get("/api/plants", params={"layer": "agri_overlay", "min_capacity": 1.0})
    assert r.status_code == 200
    body = r.json()
    names = {f["properties"]["name"] for f in body["features"]}
    assert names == {"Boone"}
    assert body["excluded_unknown_capacity"] == 1
