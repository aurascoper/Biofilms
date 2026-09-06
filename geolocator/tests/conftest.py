"""Shared fixtures for the geolocator test suite."""

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from geolocator import imagery, lattice as lat  # noqa: E402


@pytest.fixture
def isolated_lattice_cache(monkeypatch):
    """Give one test its own lattice memo object.

    Replacing the cache object, rather than clearing the process-global dict in
    place, makes ownership explicit: pytest restores the original object when
    the test ends, and no cached document can leak in from test order or a
    previous module.
    """
    cache = {"at": None, "doc": None}
    monkeypatch.setattr(lat, "_cache", cache)
    return cache


# NASA's uptime is not a property of this repository.
#
# Every test in this module runs against a synthetic upstream by default, so the
# blocking CI tier is deterministic and offline. The stub is installed at the
# urlopen seam rather than by replacing imagery._fetch, which means the real
# validation path still executes -- content-type check, JPEG magic bytes, and
# the MAX_BYTES ceiling -- and the tests below that install their own urlopen
# stub simply override this one.
#
# Set GEOLOCATOR_LIVE_GIBS=1 to run the identical suite against the real
# endpoint. That is interoperability, a separate claim from "our contracts
# hold", and it runs in a non-blocking scheduled job.
LIVE_GIBS = os.environ.get("GEOLOCATOR_LIVE_GIBS") == "1"

SYNTHETIC_JPEG = b"\xff\xd8\xff\xe0" + b"\x00" * 8192 + b"\xff\xd9"


class _StubGibsResponse:
    """Shaped like the object urlopen hands to imagery._fetch."""

    headers = {"Content-Type": "image/jpeg"}

    def __init__(self, data):
        self._data = data

    def read(self, n=-1):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture(autouse=True)
def offline_gibs(tmp_path, monkeypatch):
    """Isolate the imagery cache and stub the provider, unless running live."""
    if LIVE_GIBS:
        return
    cache = tmp_path / "imagery"
    monkeypatch.setattr(imagery, "CACHE_DIR", cache)
    monkeypatch.setattr(imagery, "RASTER", cache / "blue_marble.jpg")
    monkeypatch.setattr(imagery, "META", cache / "blue_marble.meta.json")
    monkeypatch.setattr(
        imagery.urllib.request, "urlopen", lambda *a, **k: _StubGibsResponse(SYNTHETIC_JPEG))
