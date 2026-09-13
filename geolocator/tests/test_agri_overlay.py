"""
Rev 6A agri_overlay: hash verification, and unknown-capacity handling in the two
generic endpoints (/api/stats, /api/plants) that assume capacity_mw is always a float.

capacity_mw is None for EVERY county in this layer: it has no MW figure, and |NDVI z|
was once stored there as a marker magnitude, which /api/stats summed as megawatts. The
endpoint tests use a synthetic layer with a mix of known and unknown capacity_mw, which
must not crash either endpoint, must exclude unknowns from any sum, and must report how
many were excluded rather than let that count disappear silently.

The hash tests are a positive/negative pair on purpose. Corrupting `cells` proves the
hash mechanism works at all; corrupting `generated_at` proves the widening from a
cells-only hash to a whole-payload hash actually did something, since a cells-only
scheme would have let that corruption straight through.
"""

import copy
import dataclasses
import hashlib
import json
import re
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from geolocator.api import LAYER_IDS, SOURCES, _load_agri_overlay, app  # noqa: E402
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
    # No county has a capacity: a z-score is not a megawatt figure, and it used to be
    # stored here as |z| (Boone would have read 4.0 MW). The anomaly stays in extra.
    assert by_name["Boone"]["capacity_mw"] is None
    assert by_name["Boone"]["extra"]["ndvi"]["z"] == pytest.approx(-4.0)
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


def _v1_payload():
    cells_hash = hashlib.sha256(
        json.dumps(MIXED_CELLS, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": 1, "source_repo": "agri_yield_pipeline", "source_git_sha": "x",
        "generated_at": "2026-08-27T00:00:00+00:00", "cells_sha256": cells_hash,
        "provenance": {}, "cells": copy.deepcopy(MIXED_CELLS),
    }


def test_schema_v1_cells_only_hash_still_verifies(tmp_path):
    """Backward compatibility: an old export with only cells_sha256 (no payload_sha256)
    still verifies against cells-only corruption, via the version-1 path."""
    f = tmp_path / "overlay_v1.json"
    _write(f, _v1_payload())
    result = _load_agri_overlay(f)
    assert len(result["items"]) == 2


@pytest.mark.parametrize("version, missing", [(2, "payload_sha256"), (1, "cells_sha256")])
def test_a_missing_hash_is_refused_not_skipped(tmp_path, version, missing):
    """Negative control for failing closed. The old code selected the cells-only path by
    the ABSENCE of payload_sha256 and skipped the check when cells_sha256 was absent too,
    so a v2 export with its hash deleted (or no hash at all) was served unverified."""
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS) if version == 2 else _v1_payload()
    del payload[missing]
    _write(f, payload)
    with pytest.raises(ValueError, match=f"carries no {missing}"):
        _load_agri_overlay(f)


@pytest.mark.parametrize("version", [None, 3, "2"])
def test_an_unsupported_schema_version_is_refused_before_hashing(tmp_path, version):
    """The else branch used to read every non-1 version, including a missing one, under
    the v2 rule, so an unversioned export with a valid payload hash was accepted."""
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS)
    del payload["payload_sha256"]
    if version is None:
        del payload["schema_version"]
    else:
        payload["schema_version"] = version
    payload["payload_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()  # a VALID hash for the export as written
    _write(f, payload)
    with pytest.raises(ValueError, match="unsupported schema_version"):
        _load_agri_overlay(f)


@pytest.mark.parametrize("version", [True, 1.0])
def test_a_version_that_merely_equals_1_is_not_version_1(tmp_path, version):
    """`True == 1` and `1.0 == 1`, so a malformed version with a valid cells hash used
    to enter the v1 path; the gate requires an exact int."""
    f = tmp_path / "overlay.json"
    payload = _v1_payload()
    payload["schema_version"] = version
    _write(f, payload)
    with pytest.raises(ValueError, match="unsupported schema_version"):
        _load_agri_overlay(f)


def test_a_v2_export_cannot_downgrade_to_the_cells_only_hash(tmp_path):
    """A v2 export carrying only cells_sha256 must not be verified by the v1 rule: that
    would let generated_at / source_git_sha be edited under a still-valid cells hash."""
    f = tmp_path / "overlay.json"
    payload = _payload(MIXED_CELLS)
    del payload["payload_sha256"]
    payload["cells_sha256"] = hashlib.sha256(
        json.dumps(payload["cells"], sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    _write(f, payload)
    with pytest.raises(ValueError, match="carries no payload_sha256"):
        _load_agri_overlay(f)


# ── endpoint-level: unknown capacity_mw handling ──────────────────────────────


@pytest.fixture
def agri_source(tmp_path, monkeypatch):
    """The REAL registry entry, pointed at a fixture file: dataclasses.replace keeps every
    other field exactly as registered, so these tests exercise the registration and not a
    copy of it."""
    f = tmp_path / "overlay.json"
    _write(f, _payload(MIXED_CELLS))
    src = dataclasses.replace(SOURCES["agri_overlay"], path=f)
    monkeypatch.setitem(SOURCES, "agri_overlay", src)
    return src


def test_stats_never_sums_a_z_score_as_megawatts(agri_source):
    """Boone's z of -4.0 used to reach total_capacity_mw as 4.0 MW."""
    r = client.get("/api/stats", params={"layer": "agri_overlay"})
    assert r.status_code == 200
    body = r.json()
    assert body["count"] == 2
    assert body["unknown_capacity_count"] == 2
    assert body["total_capacity_mw"] == 0.0
    assert body["capacity_by_fuel"] == {}


def test_plants_mw_filter_never_admits_a_z_score(agri_source):
    """min_capacity=1 used to pass Boone through on |z| = 4. Every county is unknown, so an
    explicit floor excludes them all and says so; the default floor keeps them all."""
    r = client.get("/api/plants", params={"layer": "agri_overlay", "min_capacity": 1.0})
    assert r.status_code == 200
    assert r.json()["features"] == []
    assert r.json()["excluded_unknown_capacity"] == 2
    r = client.get("/api/plants", params={"layer": "agri_overlay"})
    assert len(r.json()["features"]) == 2
    assert r.json()["excluded_unknown_capacity"] == 0
    assert all(f["properties"]["capacity_mw"] is None for f in r.json()["features"])


# ── endpoint-level: the unknown-capacity counter on a mixed layer ─────────────


@pytest.fixture
def mixed_layer(monkeypatch):
    """A layer with known and unknown capacities across two countries, registered under
    the overlay's id so the routes see it through the same registry."""
    items = [
        {"name": "A", "country": "MO", "latitude": 1, "longitude": 1, "color_key": "x",
         "capacity_mw": 4.0},
        {"name": "B", "country": "MO", "latitude": 1, "longitude": 1, "color_key": "x",
         "capacity_mw": None},
        {"name": "C", "country": "IL", "latitude": 1, "longitude": 1, "color_key": "x",
         "capacity_mw": None},
    ]
    src = TrackedSource(
        id="agri_overlay", path=Path(__file__), layer_class=REFERENCE,
        loader=lambda _p: {"items": items, "meta": {}}, empty={"items": [], "meta": {}},
        count_of=lambda d: len((d or {}).get("items", [])),
    )
    monkeypatch.setitem(SOURCES, "agri_overlay", src)


def test_stats_excludes_unknown_from_sum_and_reports_the_count(mixed_layer):
    body = client.get("/api/stats", params={"layer": "agri_overlay"}).json()
    assert body["count"] == 3
    assert body["unknown_capacity_count"] == 2
    assert body["total_capacity_mw"] == pytest.approx(4.0)


def test_plants_default_filter_reports_zero_excluded(mixed_layer):
    """min_capacity defaults to 0.0 -- an unknown capacity_mw must not be excluded just
    because 0.0 is a technically-passable floor. excluded_unknown_capacity must still be
    present (0), not merely absent."""
    body = client.get("/api/plants", params={"layer": "agri_overlay"}).json()
    assert len(body["features"]) == 3
    assert body["excluded_unknown_capacity"] == 0


def test_plants_active_filter_excludes_unknown_and_reports_the_count(mixed_layer):
    body = client.get("/api/plants", params={"layer": "agri_overlay", "min_capacity": 1.0}).json()
    assert {f["properties"]["name"] for f in body["features"]} == {"A"}
    assert body["excluded_unknown_capacity"] == 2


def test_layers_reports_the_vintage_and_no_retrieval_time(agri_source):
    """generated_at is when the dataset was made, not when these bytes arrived. The entry
    used to pass the same extractor to retrieved_of, so /api/layers reported both."""
    entry = next(l for l in client.get("/api/layers").json()["layers"] if l["id"] == "agri_overlay")
    assert entry["vintage"] == "2026-08-27T00:00:00+00:00"
    assert entry["retrieved_at"] is None


def test_stats_source_names_the_overlay_repo_not_wri(agri_source):
    r = client.get("/api/stats", params={"layer": "agri_overlay"})
    assert r.json()["source"] == "agri_yield_pipeline"
    # and the WRI label still belongs to the power layer
    assert client.get("/api/stats", params={"layer": "power"}).json()["source"] != "agri_yield_pipeline"


def test_plants_default_floor_still_applies_to_a_known_capacity(monkeypatch):
    """The unknown special case skipped the floor for every known value when min_capacity
    was the default 0.0, so a negative capacity_mw came back; before this layer existed
    `cap < min_capacity` excluded it. Unknowns stay in under the default floor."""
    items = [
        {"name": "neg", "country": "MO", "latitude": 1, "longitude": 1, "color_key": "x",
         "capacity_mw": -5.0},
        {"name": "unk", "country": "MO", "latitude": 1, "longitude": 1, "color_key": "x",
         "capacity_mw": None},
    ]
    src = TrackedSource(
        id="agri_overlay", path=Path(__file__), layer_class=REFERENCE,
        loader=lambda _p: {"items": items, "meta": {}}, empty={"items": [], "meta": {}},
        count_of=lambda d: len((d or {}).get("items", [])),
    )
    monkeypatch.setitem(SOURCES, "agri_overlay", src)
    body = client.get("/api/plants", params={"layer": "agri_overlay"}).json()
    assert {f["properties"]["name"] for f in body["features"]} == {"unk"}
    assert body["excluded_unknown_capacity"] == 0


def test_plants_country_filter_scopes_the_unknown_counter(mixed_layer):
    """The capacity check used to run before the country predicate, so a country-scoped
    request counted unknowns from every other country as excluded by capacity."""
    body = client.get("/api/plants", params={"layer": "agri_overlay", "min_capacity": 1.0,
                                             "country": "MO"}).json()
    assert {f["properties"]["name"] for f in body["features"]} == {"A"}
    assert body["excluded_unknown_capacity"] == 1


def test_plants_unknown_counter_is_not_truncated_by_limit(mixed_layer):
    """The loop used to break at `limit`, so an unknown row past the first page was never
    inspected and went uncounted."""
    body = client.get("/api/plants", params={"layer": "agri_overlay", "min_capacity": 1.0,
                                             "limit": 1}).json()
    assert len(body["features"]) == 1
    assert body["excluded_unknown_capacity"] == 2


def test_plants_holds_only_the_page_while_scanning_past_it(mixed_layer, monkeypatch):
    """Scanning past `limit` for the counter must not accumulate every match: the power
    layer is 34,936 rows against a default limit of 5,000. Measured, not inferred: a
    trace on the route's frame records the largest `out` it ever holds."""
    import sys
    from geolocator import api as api_mod
    peak = 0

    def tracer(frame, event, arg):
        nonlocal peak
        if frame.f_code is api_mod.plants.__code__:
            out = frame.f_locals.get("out")
            if out is not None:
                peak = max(peak, len(out))
            return tracer
        return None

    monkeypatch.setattr(api_mod, "_to_geojson", lambda items, layer="power": {"type": "FeatureCollection", "features": [None] * len(items)})
    sys.settrace(tracer)
    try:
        body = api_mod.plants(layer="agri_overlay", fuel=None, min_capacity=0.0,
                              max_capacity=None, country=None, limit=1)
    finally:
        sys.settrace(None)
    assert len(body["features"]) == 1
    assert peak == 1


# ── client: every server layer is selectable ──────────────────────────────────


def _client(name):
    return (Path(__file__).resolve().parents[1] / "static" / "app" / name).read_text()


def test_the_min_mw_filter_does_not_apply_to_a_layer_without_capacity():
    """Every overlay item has capacity_mw null, so `(capacity_mw || 0) >= minCap` hid the
    whole layer the moment the Min MW control went above zero. The layer declares mw:
    false and the predicate is gated on hasCapacity(id). Source-level: no JS runtime here."""
    js = _client("layers.js")
    block = js.split("const SITE_LAYERS = [", 1)[1].split("];", 1)[0]
    assert re.search(r"id:\s*'agri_overlay'.*mw:\s*false", block)
    predicate = js.split("const shown = activeFeatures().filter(", 1)[1].split(");", 1)[0]
    assert "!hasCapacity(id) ||" in predicate
    assert ">= minCap" in predicate


def test_the_hud_never_reports_a_non_capacity_layer_as_zero_mw():
    """main.js summed capacity_mw || 0 over every shown feature and index.html fixed the
    unit outside the readout, so selecting only the overlay displayed "0 MW". The sum is
    now over hasCapacity layers only and the unit travels with the number, so a
    non-capacity selection reads "n/a" rather than a measured zero."""
    stats = _client("main.js").split("function updateStats(", 1)[1].split("\n}\n", 1)[0]
    assert "hasCapacity(id)" in stats
    assert "withMw.length ?" in stats and "'n/a'" in stats
    html = (Path(__file__).resolve().parents[1] / "static" / "index.html").read_text()
    assert "</b> MW" not in html.split('id="stat-cap"', 1)[1].split("</div>", 1)[0]


def test_every_server_layer_is_selectable_in_the_client():
    """The browser's site-layer list is hard-coded in layers.js; a layer registered only
    on the server appears in /api/layers and can never be toggled or rendered. This
    reads the client's list from source because there is no JS runtime in this tier."""
    js = (Path(__file__).resolve().parents[1] / "static" / "app" / "layers.js").read_text()
    block = js.split("const SITE_LAYERS = [", 1)[1].split("];", 1)[0]
    client_ids = set(re.findall(r"id:\s*'([a-z_]+)'", block))
    assert client_ids == set(LAYER_IDS)
