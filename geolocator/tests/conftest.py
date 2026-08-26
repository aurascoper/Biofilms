"""Shared fixtures for the geolocator test suite."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from geolocator import lattice as lat  # noqa: E402


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
