"""
Blue Marble as a reference layer.

The whole point of this layer is that it must NOT weaken the freshness model it
plugs into. A 2004 composite retrieved today has to come back NOMINAL, while a
market bar cache last written on Aug 19 has to come back STALE, in the same
response. If imagery ever reports STALE, the distinction between "depicts an old
moment" and "is out of date" has collapsed.
"""

import hashlib
import json
import sys
import urllib.error
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from geolocator import imagery  # noqa: E402
from geolocator.api import app  # noqa: E402

client = TestClient(app)

JPEG = b"\xff\xd8\xff" + b"\x00" * 512


@pytest.fixture
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(imagery, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(imagery, "RASTER", tmp_path / "blue_marble.jpg")
    monkeypatch.setattr(imagery, "META", tmp_path / "blue_marble.meta.json")
    return tmp_path


# ── epistemic class ──────────────────────────────────────────────────────────
def test_blue_marble_is_a_reference_layer():
    st = client.get("/api/imagery/blue-marble/meta").json()
    assert st["layer_class"] == "reference"
    assert st["kind"] == "imagery"


def test_a_twenty_two_year_old_vintage_is_never_stale():
    st = client.get("/api/imagery/blue-marble/meta").json()
    assert st["vintage"] == "2004"
    assert st["status"] == "nominal"
    assert st["status"] != "stale"
    assert st["age_s"] is None, "a vintage has no age to measure against a clock"


def test_retrieved_at_does_not_masquerade_as_a_source_timestamp():
    st = client.get("/api/imagery/blue-marble/meta").json()
    assert st["source_timestamp"] is None
    assert st["retrieved_at"], "transport time is still recorded, just not as authority"
    assert st["authority"] == "dataset-vintage"


def test_provider_outage_is_unavailable_not_stale(isolated_cache, monkeypatch):
    def boom():
        raise urllib.error.URLError("Network is unreachable")
    monkeypatch.setattr(imagery, "_fetch", boom)
    st = imagery.state()
    assert st["status"] == "unavailable"
    assert st["status"] != "stale", "not having the bytes is not the same as having old ones"
    assert st["availability"] == "unavailable"
    assert st["vintage"] is None and st["authority"] is None


def test_outage_with_a_cache_serves_the_cache(isolated_cache, monkeypatch):
    monkeypatch.setattr(imagery, "_fetch", lambda: JPEG)
    assert imagery.state()["availability"] == "live"

    def boom():
        raise urllib.error.URLError("down")
    monkeypatch.setattr(imagery, "_fetch", boom)
    monkeypatch.setattr(imagery, "CACHE_TTL", imagery.timedelta(seconds=-1))  # force refetch
    st = imagery.state()
    assert st["availability"] == "cached"
    assert st["status"] == "nominal", "cached reference imagery is still valid imagery"
    assert "unreachable" in (st["note"] or "")


# ── provenance ───────────────────────────────────────────────────────────────
def test_sha256_covers_the_exact_served_bytes():
    meta = client.get("/api/imagery/blue-marble/meta").json()
    raw = client.get("/api/imagery/blue-marble").content
    assert hashlib.sha256(raw).hexdigest() == meta["content_sha256"]
    assert raw.startswith(b"\xff\xd8\xff")


def test_attribution_travels_with_the_bytes():
    r = client.get("/api/imagery/blue-marble")
    assert r.headers["x-imagery-attribution"] == "NASA Earth Observatory"
    assert r.headers["x-imagery-vintage"] == "2004"
    assert client.get("/api/imagery/blue-marble/meta").json()["attribution"] == "NASA Earth Observatory"


def test_request_descriptor_is_recorded(isolated_cache, monkeypatch):
    monkeypatch.setattr(imagery, "_fetch", lambda: JPEG)
    st = imagery.state()
    assert st["request_descriptor"].startswith(imagery.GIBS_ENDPOINT)
    assert "BlueMarble_NextGeneration" in st["request_descriptor"]


# ── SSRF surface ─────────────────────────────────────────────────────────────
@pytest.mark.parametrize("qs", [
    "?url=http://169.254.169.254/latest/meta-data/",
    "?url=file:///etc/passwd",
    "?LAYERS=SomethingElse",
    "?WIDTH=30000&HEIGHT=30000",
    "?endpoint=http://evil.example",
    "?BBOX=0,0,1,1",
])
def test_no_client_parameter_can_steer_the_request(qs):
    """The GetMap request is a server constant. Extra params must be inert."""
    base = client.get("/api/imagery/blue-marble/meta").json()
    got = client.get(f"/api/imagery/blue-marble/meta{qs}").json()
    assert got["request_descriptor"] == base["request_descriptor"]
    assert got["content_sha256"] == base["content_sha256"]
    assert got["width"] == 2048 and got["height"] == 1024


def test_there_is_no_generic_imagery_proxy_route():
    for path in ("/api/imagery", "/api/imagery/", "/api/imagery/fetch"):
        assert client.get(path, params={"url": "http://evil.example"}).status_code == 404


def test_oversized_response_is_refused(isolated_cache, monkeypatch):
    monkeypatch.setattr(imagery, "MAX_BYTES", 128)
    monkeypatch.setattr(imagery, "_fetch", imagery._fetch)
    class FakeResp:
        headers = {"Content-Type": "image/jpeg"}
        def read(self, n): return b"\xff\xd8\xff" + b"\x00" * n
        def __enter__(self): return self
        def __exit__(self, *a): return False
    monkeypatch.setattr(imagery.urllib.request, "urlopen", lambda *a, **k: FakeResp())
    st = imagery.state()
    assert st["status"] == "unavailable"
    assert "exceeds" in (st["error"] or "")


def test_non_image_response_is_refused(isolated_cache, monkeypatch):
    """GIBS reports errors as XML with HTTP 200, so content-type is the real check."""
    class FakeResp:
        headers = {"Content-Type": "text/xml"}
        def read(self, n): return b"<ServiceException/>"
        def __enter__(self): return self
        def __exit__(self, *a): return False
    monkeypatch.setattr(imagery.urllib.request, "urlopen", lambda *a, **k: FakeResp())
    st = imagery.state()
    assert st["status"] == "unavailable"
    assert "content-type" in (st["error"] or "")


# ── the layer must not weaken what it plugs into ─────────────────────────────
def test_reference_vintage_and_measured_staleness_coexist():
    """The whole justification for the layer, in one assertion.

    Blue Marble depicts 2004 and is NOMINAL. market_bars was last observed
    2026-08-19 and is STALE. Both in the same response, from the same model.
    """
    f = client.get("/api/health").json()["freshness"]
    bm, mb = f["blue_marble"], f["market_bars"]

    assert bm["layer_class"] == "reference" and bm["status"] == "nominal"
    assert bm["vintage"] == "2004"

    if mb["status"] != "unavailable":          # no bar cache on a fresh clone
        assert mb["layer_class"] == "live"
        assert mb["status"] == "stale"
        assert mb["source_timestamp"].startswith("2026-08-19")
        # The older dataset is the healthy one. That is the point.
        assert bm["vintage"] < mb["source_timestamp"][:4]


def test_lattice_boundary_is_untouched_by_imagery(isolated_lattice_cache, monkeypatch):
    """This test owns both memo state and upstream input.

    It used to clear lattice._cache in place and then hit whatever happened to
    be listening on the developer's loopback port. That made module order and
    local process state part of the assertion. A boundary test should instead
    supply the exact upstream document it is proving the imagery work leaves
    untouched.
    """
    from geolocator import lattice as _lat

    # The upstream deliberately CARRIES the fields the boundary must drop, which
    # is what makes this a negative control rather than a shape check. With a
    # clean fixture the assertions below cannot fail: widen MISSION_FIELDS to
    # admit "equity" and a sanitised upstream still passes, because there was
    # never anything present for the boundary to have leaked.
    upstream = {
        "state": {
            "generated": "2099-01-01T00:00:00+00:00",
            "mission": {
                "phase": "TEST",
                "entries_allowed": False,
                "equity": 1234.56,
                "margin_used": 9876.54,
            },
            "cells": [
                {"symbol": "AAA_USDT", "quote": "USDT", "class": "crypto",
                 "rank": 1, "has_position": True, "position_size": 4.2,
                 "signal": None},
            ],
            "econ": None,
        }
    }
    monkeypatch.setattr(_lat, "_fetch_upstream", lambda: upstream)

    d = client.get("/api/lattice").json()
    assert set(d["mission"]) <= set(_lat.MISSION_FIELDS)
    assert d["position_presence_exposed"] is False
    body = json.dumps(d).lower()
    for forbidden in ("equity", "margin_used", "position_size", "has_position"):
        assert forbidden not in body, forbidden


def test_legacy_page_is_still_byte_identical():
    served = client.get("/legacy").content
    assert hashlib.sha256(served).hexdigest() == (
        "d8e7d8400ee76cd66490f673132dbd33fc31a4e0ffaebaf06f5d5d10d73ec722")
