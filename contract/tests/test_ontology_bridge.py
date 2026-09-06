"""data/ontology_bridge.csv, enforced offline.

The bridge maps every vocabulary value, unit string and phenomenon word to an
IRI. An OLS-fetching validator would make this a network test, and a skipped
test is uncovered surface, so the label and the date it was verified travel in
the file and these checks are arithmetic and set algebra only.

What is checked: both directions between the contract sets and the file; every
unit string in the four ledgers that carry a unit column has a row; a unit's
exponents equal its quantity kind's; mirrored and minted IRIs live in disjoint
namespaces and every minted term names a mirrored parent; every phenomenon row
states its endpoint (the manuscript's rule as a column constraint); the only
basis without an IRI is `derived`; blank verified_on rows equal an explicit
list. Each rule has a planted failure below.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from physical_contract import (BRIDGE_AXES, CLAIM_EVIDENCE_BASIS,
                               EVIDENCE_NULLS, PARAMETER_EVIDENCE_BASIS,
                               PARENT_RELATIONS, PHENOMENA,
                               PROJECT_NAMESPACE, PROVENANCE_SOURCES)

REPO = Path(__file__).resolve().parents[2]
BRIDGE = REPO / "data" / "ontology_bridge.csv"
UNIT_LEDGERS = {  # file -> its unit column
    REPO / "data" / "parameter_provenance.csv": "unit",
    REPO / "data" / "species_parameter_provenance.csv": "units",
    REPO / "data" / "uncertainty" / "feedback_parameter_distributions.csv": "unit",
    REPO / "data" / "calibration" / "suspended_isotherm_proposal.csv": "units",
}
# Unit rows that no ledger carries: manuscript strings and the parent of one.
MANUSCRIPT_UNITS = frozenset({"Gy", "mGy h^-1", "R", "R h^-1"})
# Rows allowed a blank verified_on. Empty: every mirrored term was re-verified
# on 2026-09-06 and minted rows carry their mint date. A new entry here is a
# term someone wrote down without looking it up, and the list may not grow.
VERIFIED_ON_EXCEPTIONS: frozenset = frozenset()
# `derived` on either basis axis is the only row allowed no IRI.
UNMAPPED = frozenset({("derived", "value_basis"), ("derived", "claim_basis")})
# Units whose `count_of` must be filled: a count is dimensionless, so this is
# the only place the counted thing survives.
COUNTED = {"cells mm^-3": "cell", "ug cell^-1 Gy^-1": "cell", "per_s": "photon",
           "species": "species", "rows": "rows"}
DIMS = ("L", "M", "T", "I", "Theta", "N", "J")
QUDT_VERSION = "QUDT v3.5.1"
OBO = "http://purl.obolibrary.org/obo/"
# A model term is an equation or a process, never a quantity.
MODEL_TERM_PARENTS = frozenset({OBO + "IAO_0000030", OBO + "BFO_0000015"})


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    assert rows, f"{path.name} read back empty"
    return rows


@pytest.fixture(scope="module")
def rows():
    return _read(BRIDGE)


def by_axis(rows, axis):
    return {r["local_term"] for r in rows if r["axis"] == axis}


def exponents(r):
    return tuple(int(r[d]) for d in DIMS)


def bridge_problems(rows) -> list:
    """Every structural rule, as a list so planted rows can be checked."""
    out = []
    mirrored = {r["iri"] for r in rows if r["namespace"] == "mirrored"}
    by_iri = {r["iri"]: r for r in rows if r["iri"]}
    kinds = {r["local_term"]: r for r in rows if r["axis"] == "quantity_kind"}
    for r in rows:
        key = (r["local_term"], r["axis"])
        if r["axis"] not in BRIDGE_AXES:
            out.append(f"{key}: unknown axis")
        ns, iri = r["namespace"], r["iri"]
        if ns == "mirrored":
            if not iri.startswith("http") or iri.startswith(PROJECT_NAMESPACE):
                out.append(f"{key}: mirrored row carries a project IRI")
        elif ns == "minted":
            if not iri.startswith(PROJECT_NAMESPACE):
                out.append(f"{key}: minted row carries an external IRI")
            if r["axis"] == "null":
                if r["nearest_parent"]:
                    out.append(f"{key}: a null has no external parent")
            elif r["nearest_parent"] not in mirrored:
                out.append(f"{key}: nearest_parent is not a mirrored row")
            elif r["parent_relation"] not in PARENT_RELATIONS:
                out.append(f"{key}: parent_relation {r['parent_relation']!r}")
            elif r["parent_relation"] == "same_dimension":
                par = by_iri[r["nearest_parent"]]
                if par["axis"] not in ("unit", "quantity_kind") or exponents(par) != exponents(r):
                    out.append(f"{key}: same_dimension parent {par['local_term']} has other exponents")
            if r["axis"] == "model_term" and r["nearest_parent"] not in MODEL_TERM_PARENTS:
                out.append(f"{key}: a model term's parent is an information entity or a process")
        elif ns == "unmapped":
            if key not in UNMAPPED or iri:
                out.append(f"{key}: only derived may be unmapped")
        else:
            out.append(f"{key}: namespace {ns!r}")
        if ns != "minted" and r["parent_relation"]:
            out.append(f"{key}: parent_relation on a row that is not minted")
        if r["xref"] and r["xref"] not in mirrored:
            out.append(f"{key}: xref is not a mirrored row")
        needs_definition = (r["axis"] == "model_term"
                            or (r["axis"] == "quantity_kind" and (ns == "minted" or exponents(r) == (0,) * 7)))
        if needs_definition and not r["definition"].strip():
            out.append(f"{key}: no definition; exponents alone cannot identify it")
        if not iri and key not in UNMAPPED:
            out.append(f"{key}: no IRI")
        if not r["verified_on"] and key not in UNMAPPED and key[0] not in VERIFIED_ON_EXCEPTIONS:
            out.append(f"{key}: blank verified_on outside the exception list")
        if r["axis"] == "phenomenon" and not r["endpoint_assay"].strip():
            out.append(f"{key}: phenomenon with no statable endpoint")
        if r["axis"] == "unit":
            kind = r["quantity_kind"]
            if kind == "none":
                if r["local_term"] != "n/a":
                    out.append(f"{key}: only n/a may have no quantity kind")
            elif kind not in kinds:
                out.append(f"{key}: quantity kind {kind!r} has no row")
            elif exponents(r) != exponents(kinds[kind]):
                out.append(f"{key}: exponents {exponents(r)} differ from {kind} {exponents(kinds[kind])}")
            if COUNTED.get(r["local_term"], "") != r["count_of"]:
                out.append(f"{key}: count_of {r['count_of']!r}")
    return out


# --- clean cases ----------------------------------------------------------

def test_header_names_the_qudt_version():
    head = BRIDGE.read_text(encoding="utf-8").split("\n", 12)[:12]
    assert any(QUDT_VERSION in l for l in head)


def test_bridge_is_structurally_clean(rows):
    assert bridge_problems(rows) == []
    keys = [(r["local_term"], r["axis"]) for r in rows]
    assert len(keys) == len(set(keys)), "duplicate (term, axis)"


def test_contract_sets_and_bridge_agree_in_both_directions(rows):
    assert by_axis(rows, "value_basis") == PARAMETER_EVIDENCE_BASIS
    assert by_axis(rows, "claim_basis") == CLAIM_EVIDENCE_BASIS
    assert by_axis(rows, "source") == PROVENANCE_SOURCES
    assert by_axis(rows, "null") == {"blank" if v == "" else v for v in EVIDENCE_NULLS}
    assert by_axis(rows, "phenomenon") == PHENOMENA
    assert by_axis(rows, "status") == set(), "status rows are not in scope yet"


def test_every_ledger_unit_string_has_a_row(rows):
    units = by_axis(rows, "unit")
    seen = set()
    for path, col in UNIT_LEDGERS.items():
        strings = {r[col] for r in _read(path)}
        assert strings, f"{path.name}: no unit strings"
        missing = strings - units
        assert not missing, f"{path.name}: {sorted(missing)}"
        seen |= strings
    assert units - seen == MANUSCRIPT_UNITS, units - seen


def test_the_two_rate_strings_are_not_dose_rates(rows):
    r = {x["local_term"]: x for x in rows if x["axis"] == "unit"}
    assert r["rad hr^-1"]["quantity_kind"] == "AngularVelocity" and exponents(r["rad hr^-1"]) == (0, 0, -1, 0, 0, 0, 0)
    assert r["R h^-1"]["quantity_kind"] == "ExposureRate" and exponents(r["R h^-1"]) == (0, -1, 0, 1, 0, 0, 0)
    assert exponents(r["mGy h^-1"]) == (2, 0, -3, 0, 0, 0, 0)


# --- planted failures ------------------------------------------------------

def _plant(rows, term, axis, **changes):
    copy = [dict(r) for r in rows]
    target = next(r for r in copy if r["local_term"] == term and r["axis"] == axis)
    before = {k: target[k] for k in changes}
    target.update(changes)
    assert {k: target[k] for k in changes} != before, "mutation did not change the row"
    return copy


def test_rad_per_hour_mapped_to_the_dose_rate_fails_on_arithmetic(rows):
    planted = _plant(rows, "rad hr^-1", "unit", quantity_kind="AbsorbedDoseRate")
    assert any("exponents (0, 0, -1, 0, 0, 0, 0) differ from AbsorbedDoseRate (2, 0, -3" in p
               for p in bridge_problems(planted))


def test_roentgen_per_hour_mapped_to_the_dose_rate_fails_on_arithmetic(rows):
    planted = _plant(rows, "R h^-1", "unit", quantity_kind="AbsorbedDoseRate")
    assert any("exponents (0, -1, 0, 1, 0, 0, 0) differ from AbsorbedDoseRate" in p
               for p in bridge_problems(planted))


def test_namespace_rules_fire_both_ways(rows):
    minted_external = _plant(rows, "radiotropic", "phenomenon", iri="http://purl.obolibrary.org/obo/GO_0009606")
    assert any("minted row carries an external IRI" in p for p in bridge_problems(minted_external))
    mirrored_project = _plant(rows, "Gy", "unit", iri=PROJECT_NAMESPACE + "gray")
    assert any("mirrored row carries a project IRI" in p for p in bridge_problems(mirrored_project))
    orphan = _plant(rows, "declared", "value_basis", nearest_parent=PROJECT_NAMESPACE + "nothing")
    assert any("nearest_parent is not a mirrored row" in p for p in bridge_problems(orphan))


def test_endpoint_rule_derived_rule_and_verified_on_rule_fire(rows):
    assert any("no statable endpoint" in p
               for p in bridge_problems(_plant(rows, "radiotrophic", "phenomenon", endpoint_assay="  ")))
    given_an_iri = _plant(rows, "derived", "value_basis", iri="http://purl.obolibrary.org/obo/ECO_0000000")
    assert any("only derived may be unmapped" in p for p in bridge_problems(given_an_iri))
    another_unmapped = _plant(rows, "proxy", "value_basis", namespace="unmapped", iri="")
    assert any("only derived may be unmapped" in p for p in bridge_problems(another_unmapped))
    assert any("blank verified_on" in p
               for p in bridge_problems(_plant(rows, "cm", "unit", verified_on="")))
    assert any("count_of" in p
               for p in bridge_problems(_plant(rows, "cells mm^-3", "unit", count_of="")))


def test_step3_controls_fire(rows):
    orphan = _plant(rows, "henry_isotherm", "model_term", nearest_parent=PROJECT_NAMESPACE + "nothing")
    assert any("nearest_parent is not a mirrored row" in p for p in bridge_problems(orphan))
    quantity_parent = _plant(rows, "donnan_dialysis", "model_term", nearest_parent="http://qudt.org/vocab/quantitykind/Velocity")
    assert any("information entity or a process" in p for p in bridge_problems(quantity_parent))
    nameless = _plant(rows, "thiele_modulus", "quantity_kind", definition="")
    assert any("no definition" in p for p in bridge_problems(nameless))
    wrong_parent = _plant(rows, "porosity", "quantity_kind", nearest_parent="http://qudt.org/vocab/quantitykind/Length")
    assert any("same_dimension parent Length has other exponents" in p for p in bridge_problems(wrong_parent))
    bad_xref = _plant(rows, "radiolysis_of_water", "model_term", xref=PROJECT_NAMESPACE + "hydroxyl")
    assert any("xref is not a mirrored row" in p for p in bridge_problems(bad_xref))
