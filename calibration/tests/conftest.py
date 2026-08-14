"""Shared fixtures for the calibration tiers."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT_MARKER = "biofilms_potts.jl"


@pytest.fixture(scope="session")
def repo_root():
    for parent in Path(__file__).resolve().parents:
        if (parent / REPO_ROOT_MARKER).exists():
            return parent
    raise RuntimeError("repository root not found")


@pytest.fixture(scope="session")
def data_dir(repo_root):
    return repo_root / "data" / "calibration"


@pytest.fixture(scope="session")
def config_dir(repo_root):
    return repo_root / "config"
