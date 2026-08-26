"""Route-level behaviour: layer validation, freshness surfacing, rollback path."""

import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geolocator.api import LAYER_IDS, app  # noqa: E402

client = TestClient(app)


# ── unknown layers are an error, not an empty result ─────────────────────────
@pytest.mark.parametrize("route", ["/api/plants", "/api/fuels", "/api/stats"])
def test_unknown_layer_404s(route):
    r = client.get(route, params={"layer": "does_not_exist"})
    assert r.status_code == 404
    assert "unknown layer" in r.json()["detail"]


@pytest.mark.parametrize("route", ["/api/plants", "/api/fuels", "/api/stats"])
@pytest.mark.parametrize("layer", LAYER_IDS)
def test_known_layers_still_work(route, layer):
    assert client.get(route, params={"layer": layer}).status_code == 200


def test_empty_result_is_distinguishable_from_unknown_layer():
    """The bug this replaces: a typo and a genuinely empty filter looked alike."""
    typo = client.get("/api/plants", params={"layer": "powr"})
    empty = client.get("/api/plants", params={"layer": "power", "fuel": "Unobtainium"})
    assert typo.status_code == 404
    assert empty.status_code == 200 and empty.json()["features"] == []


# ── freshness is surfaced, per layer ─────────────────────────────────────────
def test_health_reports_class_and_authority_per_layer():
    d = client.get("/api/health").json()
    assert set(d["freshness"]) == set(LAYER_IDS) | {
        "market_bars", "blue_marble", "market_correlation"}
    for name, st in d["freshness"].items():
        assert st["layer_class"] in ("live", "reference", "authored")
        assert st["status"] in (
            "nominal", "loading", "degraded", "stale", "fallback", "unavailable")


def test_authored_layers_never_claim_a_source_timestamp():
    fresh = client.get("/api/health").json()["freshness"]
    for name in ("nuclear", "battery", "gridcoin", "stars"):
        st = fresh[name]
        assert st["layer_class"] == "authored"
        assert st["status"] == "fallback"
        assert st["authority"] == "filesystem-mtime"


def test_measured_band_carries_its_own_freshness():
    """Its authority is MAX(bar_date), so it can be stale while everything else is not.

    market_bars.sqlite is gitignored (regenerable via `python3 -m geolocator.market
    backfill`), so a fresh clone legitimately has no cache. `unavailable` with no
    authority is the correct answer there, not a failure.
    """
    d = client.get("/api/links").json()
    st = d["measured_freshness"]
    assert st["layer_class"] == "live"
    if st["status"] == "unavailable":
        assert st.get("authority") is None       # no cache: nothing to be current about
    else:
        assert st["authority"] == "source-timestamp"


def test_layers_route_colours_worldgrid_by_fuel():
    layers = {x["id"]: x for x in client.get("/api/layers").json()["layers"]}
    assert layers["worldgrid"]["color_field"] == "primary_fuel"
    assert layers["power"]["color_field"] == "primary_fuel"
    assert layers["nuclear"]["color_field"] == "stage"


# ── the rollback path is real ────────────────────────────────────────────────
def test_legacy_route_serves_the_preserved_globe():
    r = client.get("/legacy")
    assert r.status_code == 200 and "<title>" in r.text


def test_links_still_ship_the_epistemic_statuses():
    feats = client.get("/api/links").json()["features"]
    statuses = {f["properties"]["status"] for f in feats}
    assert statuses <= {"structural", "measured", "speculative"}
    assert "structural" in statuses
    # Bands with no geography must keep shipping geometry: null, not a fake point.
    for f in feats:
        if f["properties"]["a"]["frame"] == "fuel":
            assert f["geometry"] is None


# ── versionless module graph must always revalidate ──────────────────────────
@pytest.mark.parametrize("mod", [
    "/app/main.js", "/app/layerManager.js", "/app/imagery.js", "/app/styles/thermal.js",
])
def test_app_modules_are_never_served_without_revalidation(mod):
    """Every module lives at a stable path with no content hash.

    Without an explicit Cache-Control, StaticFiles sends only Last-Modified and
    ETag, which lets a browser apply heuristic freshness and reuse a module
    WITHOUT revalidating. That was observed in the field: after an edit,
    layerManager.js came back with transferSize 0 while its siblings took 304s,
    so the page ran a mixed old/new module graph and the symptom presented as an
    application bug. `no-cache` means revalidate-before-use, not do-not-store --
    unchanged modules still answer 304.
    """
    r = client.get(mod)
    assert r.status_code == 200
    assert r.headers.get("cache-control") == "no-cache", mod


def test_index_and_legacy_are_not_accidentally_cacheable_forever():
    for path in ("/", "/legacy"):
        cc = client.get(path).headers.get("cache-control", "")
        assert "max-age" not in cc or "no-cache" in cc


# ── the measured band's producer is optional, and its absence is explicit ────
def test_absent_market_adapter_is_declared_not_silent():
    """An uninstalled producer must not read as "no correlations found".

    Absence of the producer and a null result from the producer are different
    facts. Collapsing them is the same error this codebase refuses everywhere
    else, and a swallowed ImportError collapses them by default.
    """
    from geolocator.api import market_adapter_installed

    st = client.get("/api/health").json()["freshness"]["market_correlation"]
    assert st["layer_class"] == "live"
    if market_adapter_installed():
        assert st["installed"] is True
    else:
        assert st["installed"] is False
        assert st["status"] == "unavailable"
        assert "not installed" in st["reason"]
        assert "geolocator.market" in st["reason"]


def test_adapter_contract_is_published_so_it_can_be_supplied():
    """A future user must be able to add a producer without editing api.py."""
    from geolocator.api import market_adapter_installed

    st = client.get("/api/links").json()["measured_freshness"]
    assert st["id"] == "market_correlation"
    if not market_adapter_installed():
        assert "fuel_correlations()" in st["contract"]


def test_missing_adapter_does_not_take_the_other_bands_with_it():
    feats = client.get("/api/links").json()["features"]
    statuses = {f["properties"]["status"] for f in feats}
    assert "structural" in statuses
    assert len(feats) > 0
