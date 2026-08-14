"""SCHEMA-INTEGRITY tests for data/parameter_provenance.csv.

These prove the ledger is internally consistent and mechanically complete —
that it covers the loader's required keys, uses the declared vocabularies, and
never claims readiness without a resolvable, located source. They CANNOT prove
a value is honest: a fabricated number with a plausible citation passes every
assertion here. Verifying a value against its primary source is a review
process, not something a test can do offline.

The point of importing the loader's PUBLIC registry (FIELD_SPECS,
required_keys) rather than restating the keys is that the ledger cannot drift
from the config contract without a test failing.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from biofilm_openmc.config import STAGES, FIELD_SPECS, required_keys

DATA = Path(__file__).resolve().parents[2] / "data"
LEDGER = DATA / "parameter_provenance.csv"
SOURCES = DATA / "sources.csv"

STATUS = {"ready", "provisional", "blocked", "not_applicable"}
MAPPING = {"mapped", "unmapped", "unsupported"}
COMPAT = {"direct", "requires_transform", "unsupported"}
PROVENANCE = {"published_replica", "certified_component", "engineered_composite",
              "declared"}
EVIDENCE = {"direct_measurement", "assay_certificate", "manufacturer_datasheet",
            "primary_literature", "evaluated_nuclear_data", "derived",
            "declared", "synthetic"}
DISTRIBUTION = {"point", "discrete_lines", "uniform", "lognormal", "n/a", ""}
RANK = {"high", "medium", "low", "unknown"}
STAGE_VALUES = set(STAGES) | {"none"}


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


@pytest.fixture(scope="module")
def rows():
    return _read(LEDGER)


@pytest.fixture(scope="module")
def source_ids():
    return {r["source_id"] for r in _read(SOURCES)}


def test_uniqueness_key_holds(rows):
    """A ledger keyed on config_key alone silently merges reference systems."""
    keys = [(r["reference_system_id"], r["required_for_stage"], r["config_key"])
            for r in rows]
    dupes = {k for k in keys if keys.count(k) > 1}
    assert not dupes, f"duplicate (system, stage, key) rows: {sorted(dupes)}"


def test_controlled_vocabularies(rows):
    for r in rows:
        where = f"{r['reference_system_id']}/{r['config_key']}"
        assert r["status"] in STATUS, f"{where}: status {r['status']!r}"
        assert r["mapping_status"] in MAPPING, f"{where}: {r['mapping_status']!r}"
        assert r["model_compatibility"] in COMPAT, f"{where}: {r['model_compatibility']!r}"
        assert r["system_provenance"] in PROVENANCE, f"{where}: {r['system_provenance']!r}"
        assert r["evidence_basis"] in EVIDENCE, f"{where}: {r['evidence_basis']!r}"
        assert r["distribution"] in DISTRIBUTION, f"{where}: {r['distribution']!r}"
        assert r["sensitivity_rank"] in RANK, f"{where}: {r['sensitivity_rank']!r}"
        assert r["required_for_stage"] in STAGE_VALUES, \
            f"{where}: stage {r['required_for_stage']!r}"


def test_mapped_rows_name_a_real_config_key(rows):
    """Every mapped key must exist in the loader's registry — no phantom keys.
    Material-class rows are exempt: the class table is open by design."""
    known = {f.dotted for f in FIELD_SPECS}
    for r in rows:
        key = r["config_key"]
        if r["mapping_status"] != "mapped" or key.startswith("materials."):
            continue
        assert key in known, f"{key} is not a config key the loader knows about"


def test_mapped_rows_sit_at_the_stage_the_loader_says(rows):
    known = {f.dotted: f.min_stage for f in FIELD_SPECS}
    for r in rows:
        key = r["config_key"]
        if r["mapping_status"] != "mapped" or key not in known:
            continue
        assert r["required_for_stage"] == known[key], (
            f"{key}: ledger says {r['required_for_stage']}, "
            f"loader says {known[key]}")


def test_each_system_covers_the_strongest_stage_it_claims(rows):
    """A system that claims a stage must carry a row for every key that stage
    requires — otherwise its coverage is unknowable rather than merely
    incomplete. `not_applicable` rows are excluded: a transport benchmark is
    not 'missing' a CPM clock it never needs."""
    by_system: dict[str, list[dict]] = {}
    for r in rows:
        if r["mapping_status"] == "mapped" and r["status"] != "not_applicable":
            by_system.setdefault(r["reference_system_id"], []).append(r)

    for system, srows in by_system.items():
        claimed = {r["required_for_stage"] for r in srows} - {"none"}
        if not claimed:
            continue
        strongest = max(claimed, key=STAGES.index)
        have = {r["config_key"] for r in srows}
        missing = required_keys(strongest) - have
        assert not missing, (
            f"{system} claims stage {strongest} but has no row for: "
            f"{sorted(missing)}")
        assert any(k.startswith("materials.") for k in have), \
            f"{system} claims stage {strongest} but declares no material class"


def test_ready_rows_are_fully_attributed(rows):
    """Readiness requires a value AND a resolvable, located source. This is the
    weakest useful check: it catches an unsourced claim, not a false one."""
    for r in rows:
        if r["status"] != "ready":
            continue
        where = f"{r['reference_system_id']}/{r['config_key']}"
        assert r["value"].strip(), f"{where}: ready with no value"
        assert r["source_id"].strip(), f"{where}: ready with no source_id"
        assert r["source_locator"].strip(), f"{where}: ready with no source_locator"
        assert r["source_locator"].strip() != "NOT LOCATED", \
            f"{where}: ready but its source was never located"


def test_every_source_id_resolves(rows, source_ids):
    for r in rows:
        sid = r["source_id"].strip()
        if sid:
            assert sid in source_ids, \
                f"{r['config_key']}: source_id {sid!r} not in data/sources.csv"


def test_derived_rows_show_their_arithmetic(rows):
    for r in rows:
        if r["evidence_basis"] == "derived":
            assert r["derivation"].strip(), \
                f"{r['config_key']}: derived with no derivation recorded"


def test_uncertain_distributions_carry_parameters(rows):
    """`uniform` with no interval says nothing."""
    for r in rows:
        if r["distribution"] in {"uniform", "lognormal"}:
            assert r["distribution_parameters"].strip(), \
                f"{r['config_key']}: {r['distribution']} with no parameters"


def test_unsupported_rows_are_never_ready(rows):
    """A value the model cannot represent must not read as usable."""
    for r in rows:
        if r["model_compatibility"] == "unsupported":
            assert r["status"] != "ready", (
                f"{r['config_key']}: marked ready but model_compatibility is "
                "unsupported")
            assert r["applicability"].strip(), \
                f"{r['config_key']}: unsupported with no explanation"
