"""
The vendored WRI snapshot, and the vintage/retrieval distinction.

A provenance file that can drift from the bytes it describes is worse than no
provenance file, because it converts "unknown" into "confidently wrong". So the
declared SHA-256 is checked against the actual committed CSV here, and that test
is the reason the declaration is trustworthy at all.
"""

import hashlib
import json
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geolocator.api import REPO, SOURCES, app  # noqa: E402

client = TestClient(app)
PROV_PATH = REPO / "power_plant_database_global.provenance.json"
CSV_PATH = REPO / "power_plant_database_global.csv"


@pytest.fixture(scope="module")
def prov():
    return json.loads(PROV_PATH.read_text())


def test_declared_sha256_matches_the_vendored_bytes(prov):
    """The whole reason the declaration is worth anything."""
    if not CSV_PATH.exists():
        pytest.skip("full WRI CSV not present (sample fallback in use)")
    actual = hashlib.sha256(CSV_PATH.read_bytes()).hexdigest()
    assert actual == prov["sha256"], (
        "provenance has drifted from the file it describes; if the dataset was "
        "intentionally replaced, that is a new dataset-version commit")
    assert CSV_PATH.stat().st_size == prov["bytes"]


def test_declaration_carries_publisher_licence_and_attribution(prov):
    assert prov["publisher"] == "World Resources Institute"
    assert prov["license"] == "CC BY 4.0"
    assert "World Resources Institute" in prov["attribution"]
    assert prov["source"].startswith("https://")
    assert prov["semantic_class"] == "reference"


def test_declaration_separates_vintage_from_retrieval(prov):
    """1.3.0 is what the data depicts. 2026-08-17 is when the bytes arrived."""
    assert prov["vintage"] == "1.3.0"
    assert prov["retrieved"] != prov["vintage"]
    assert prov["data_coverage_through"] == 2019


def test_snapshot_is_declared_immutable_not_a_cache(prov):
    joined = " ".join(prov["notes"]).lower()
    assert "not a cache" in joined
    assert "new explicit dataset-version commit" in joined


# ── the layers must report the declaration, not a retrieval date ─────────────
@pytest.mark.parametrize("layer", ["power", "worldgrid"])
def test_reference_layers_report_the_declared_vintage(layer):
    """`worldgrid` derives from a file OUTSIDE this repository.

    Its source is the lattice board's energy_world.json, which a public clone
    legitimately does not have -- CI found this within minutes of existing,
    because every local run passed only on the machine that happened to hold
    that file. Absent, the layer must report `unavailable` with no vintage;
    present, it must report the declared release. Both are correct answers, and
    asserting only the first would be asserting the author's filesystem.
    """
    st = SOURCES[layer].state()
    assert st.layer_class == "reference"
    if st.status == "unavailable":
        assert layer == "worldgrid", "only the cross-repo source may be absent"
        # UNAVAILABLE has two causes: the file is absent, or its loader threw.
        # Only the first is a portability fact; the second is a parsing
        # regression, and accepting both would make this allowance hide the very
        # class of bug the suite exists to catch.
        assert st.error and st.error.startswith("missing:"), st.error
        assert st.vintage is None and st.authority is None
        return
    assert st.status == "nominal"
    assert st.authority == "dataset-vintage"
    assert st.vintage == "1.3.0", "a WRI release version, not a date"


@pytest.mark.parametrize("layer", ["power", "worldgrid"])
def test_retrieval_date_is_never_the_vintage(layer):
    """The bug this replaced: `retrieved` was reported under `dataset-vintage`.

    Skipped for a layer whose source is absent -- there is no retrieval date to
    confuse with a vintage when nothing was retrieved.
    """
    st = SOURCES[layer].state()
    if st.status == "unavailable":
        assert layer == "worldgrid"
        # Absent is a skip. Broken is a failure. See the note above.
        assert st.error and st.error.startswith("missing:"), st.error
        pytest.skip(f"{layer} source not present in this checkout")
    assert st.retrieved_at == "2026-08-17"
    assert st.vintage != st.retrieved_at
    assert st.source_timestamp is None


def test_reference_layers_never_go_stale_by_clock():
    """A reference layer may be absent, but it must never be STALE.

    That is the whole distinction: not having the data and having old data are
    different states, and only one of them is a clock measurement.
    """
    f = client.get("/api/health").json()["freshness"]
    for name in ("power", "worldgrid", "blue_marble"):
        assert f[name]["status"] != "stale", name
        assert f[name]["layer_class"] == "reference"


# ── cache policy is an invariant, not a one-off ──────────────────────────────
@pytest.mark.parametrize("path", ["/api/health", "/api/layers", "/api/lattice", "/api/links"])
def test_api_responses_are_never_stored(path):
    """/api/lattice carries a redacted view of a trading book. Not on disk."""
    assert client.get(path).headers.get("cache-control") == "no-store"


@pytest.mark.parametrize("path", ["/", "/legacy", "/app/main.js", "/app/styles/noir.js"])
def test_versionless_assets_must_revalidate(path):
    assert client.get(path).headers.get("cache-control") == "no-cache"


def test_public_imagery_revalidates_rather_than_never_storing():
    """Large, unchanging, not sensitive -- a 304 on 209 KB is worth having."""
    assert client.get("/api/imagery/blue-marble").headers.get("cache-control") == "no-cache"


def test_nothing_claims_a_long_immutable_cache():
    """Nothing here is content-hashed, so nothing here earns max-age."""
    for path in ("/", "/legacy", "/app/main.js", "/api/health", "/api/imagery/blue-marble"):
        cc = client.get(path).headers.get("cache-control", "")
        assert "max-age" not in cc and "immutable" not in cc, path


def test_layers_route_exposes_vintage_and_retrieval_separately():
    """The provenance footer renders both; conflating them printed
    'retrieved 1.3.0' on the live page."""
    layers = {x["id"]: x for x in client.get("/api/layers").json()["layers"]}
    power = layers["power"]
    assert power["vintage"] == "1.3.0"
    assert power["retrieved_at"] == "2026-08-17"
    assert power["vintage"] != power["retrieved_at"]
