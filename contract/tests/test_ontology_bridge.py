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
                               CONVERSION_STATUSES, EVIDENCE_NULLS,
                               PARAMETER_EVIDENCE_BASIS, PARENT_RELATIONS,
                               PHENOMENA, PROJECT_NAMESPACE,
                               PROVENANCE_SOURCES, SUBSTITUTION_RELATIONS,
                               UNIT_SYSTEMS)

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
           "species": "species", "rows": "rows", "sites": "site"}
DIMS = ("L", "M", "T", "I", "Theta", "N", "J", "H")
REQUIREMENTS = REPO / "data" / "calibration" / "reference_d_requirements.csv"
# Unit rows that carry no quantity kind: a marker for "no value" and a marker for
# "no unit yet". Anything else with kind `none` is a number hiding its unit.
NO_KIND = frozenset({"n/a", "placeholder"})
TEX = REPO / "preprint" / "modeling_radioresistance_and_radiotropic_fitness.tex"
SPECIES_TABLE = REPO / "data" / "species_parameter_provenance.csv"
PARAMETERS = REPO / "data" / "parameter_provenance.csv"
POTTS = REPO / "biofilms_potts.jl"
# Code species tag -> species table name; OI has no Table 2 row at all.
SPECIES = {"CN": "C. neoformans", "DR": "D. radiodurans", "CS": "C. sphaerospermum",
           "BS": "B. subtilis", "AN": "A. niger", "SO": "S. oneidensis", "OI": None}
# A code coefficient with no tabulated row must be zero, except these, which
# the source comments as estimated.
UNTABULATED_NONZERO = frozenset({("beta_s_ion", "OI")})
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


def requirement_status(path=REQUIREMENTS) -> dict:
    return {r["requirement_id"]: r["status"] for r in _read(path)}


def bridge_problems(rows, requirements=None) -> list:
    """Every structural rule, as a list so planted rows can be checked."""
    out = []
    requirements = requirement_status() if requirements is None else requirements
    by_term = {(r["local_term"], r["axis"]): r for r in rows}
    ledger_basis = {r["claim_id"]: r["evidence_basis"] for r in _read(SPECIES_TABLE)}
    ledger_basis.update({r["config_key"]: r["evidence_basis"] for r in _read(PARAMETERS)})
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
            if r["axis"] in ("null", "conversion", "coefficient"):
                if r["nearest_parent"]:
                    out.append(f"{key}: a {r['axis']} row has no external parent")
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
        if r["axis"] == "conversion":
            st, factor, req = r["conversion_status"], r["conversion_factor"].strip(), r["requirement_id"].strip()
            if st not in CONVERSION_STATUSES:
                out.append(f"{key}: conversion_status {st!r}")
            elif st == "blocked":
                if factor:
                    out.append(f"{key}: blocked with a factor")
                if requirements.get(req) != "awaiting_measurement":
                    out.append(f"{key}: blocked on {req!r}, whose status is {requirements.get(req)!r}; fill the factor or re-block")
            elif st == "declared":
                if not factor or r["basis"] not in ("synthetic", "declared") or r["unit_system"] != "synthetic_reference":
                    out.append(f"{key}: a declared conversion needs a factor, a synthetic or declared basis and the synthetic_reference system")
            elif st == "ready":
                if not factor:
                    out.append(f"{key}: ready with no factor")
                if requirements.get(req) not in ("satisfied", "measured"):
                    out.append(f"{key}: ready while {req!r} is {requirements.get(req)!r}")
                dep = r["depends_on"].strip()
                if dep and by_term.get((dep, "conversion"), {}).get("conversion_status") != "ready":
                    out.append(f"{key}: ready before {dep!r} is")
        if r["axis"] == "coefficient":
            unit = by_term.get((r["unit"], "unit"))
            if unit is None:
                out.append(f"{key}: unit {r['unit']!r} is not a unit row")
            elif unit["quantity_kind"] != r["quantity_kind"] or unit["unit_system"] != r["unit_system"]:
                out.append(f"{key}: kind or system disagrees with its unit row")
            if r["unit_system"] == "SI":
                if r["substitution_of"] or r["relation"]:
                    out.append(f"{key}: a tabulated prior substitutes nothing")
            else:
                prior = by_term.get((r["substitution_of"], "coefficient"))
                if prior is None:
                    out.append(f"{key}: substitution_of names no coefficient row")
                elif prior["unit_system"] != "SI" or prior["unit"] == r["unit"]:
                    out.append(f"{key}: a shipped coefficient must substitute a tabulated SI prior in another unit")
                if r["relation"] not in SUBSTITUTION_RELATIONS:
                    out.append(f"{key}: relation {r['relation']!r}")
                if not r["code_location"].strip():
                    out.append(f"{key}: a shipped coefficient names where it ships")
            if not r["ledger_rows"].strip():
                out.append(f"{key}: no ledger rows")
            # The row's basis is the ledger's, not the bridge's: it must agree
            # with every ledger row it names, and a hard-coded literal is a
            # declared choice by definition (PP-T2-29).
            if r["basis"] not in PARAMETER_EVIDENCE_BASIS:
                out.append(f"{key}: basis {r['basis']!r}")
            else:
                for lid in r["ledger_rows"].split(";"):
                    lb = ledger_basis.get(lid)
                    if lb is not None and lb != r["basis"]:
                        out.append(f"{key}: basis {r['basis']!r} disagrees with {lid} ({lb!r})")
            if r["relation"] == "hard_coded_replacement" and r["basis"] != "declared":
                out.append(f"{key}: a hard-coded literal has no basis but declared")
        elif r["substitution_of"] or r["relation"] or r["unit"]:
            out.append(f"{key}: coefficient columns on a non-coefficient row")
        if r["axis"] == "coefficient" or r["conversion_status"] or r["conversion_factor"] or r["requirement_id"]:
            pass
        if r["axis"] != "conversion" and (r["conversion_status"] or r["conversion_factor"] or r["requirement_id"]):
            out.append(f"{key}: conversion columns on a non-conversion row")
        needs_definition = (r["axis"] in ("model_term", "conversion", "coefficient")
                            or (r["axis"] == "quantity_kind" and (ns == "minted" or exponents(r) == (0,) * 8)))
        if needs_definition and not r["definition"].strip():
            out.append(f"{key}: no definition; exponents alone cannot identify it")
        if not iri and key not in UNMAPPED:
            out.append(f"{key}: no IRI")
        if not r["verified_on"] and key not in UNMAPPED and key[0] not in VERIFIED_ON_EXCEPTIONS:
            out.append(f"{key}: blank verified_on outside the exception list")
        if r["axis"] == "phenomenon" and not r["endpoint_assay"].strip():
            out.append(f"{key}: phenomenon with no statable endpoint")
        if r["axis"] == "unit":
            if r["unit_system"] not in UNIT_SYSTEMS:
                out.append(f"{key}: unit_system {r['unit_system']!r}")
            elif r["unit_system"] == "SI" and r["H"] != "0":
                out.append(f"{key}: an SI unit has no Hamiltonian exponent")
            kind = r["quantity_kind"]
            if r["local_term"] in NO_KIND and kind != "none":
                out.append(f"{key}: a placeholder may not acquire a quantity kind")
            if kind == "none":
                if r["local_term"] not in NO_KIND:
                    out.append(f"{key}: only {sorted(NO_KIND)} may have no quantity kind")
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
    # The reverse direction is over SI rows: a lattice unit comes from the
    # code, not a ledger, and must say which field in its definition.
    units = by_axis(rows, "unit")
    lattice = [r for r in rows if r["axis"] == "unit" and r["unit_system"] != "SI"]
    assert lattice and all(r["definition"].strip() for r in lattice), "a lattice unit with no definition"
    units -= {r["local_term"] for r in lattice}
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
    assert r["rad hr^-1"]["quantity_kind"] == "AngularVelocity" and exponents(r["rad hr^-1"]) == (0, 0, -1, 0, 0, 0, 0, 0)
    assert r["R h^-1"]["quantity_kind"] == "ExposureRate" and exponents(r["R h^-1"]) == (0, -1, 0, 1, 0, 0, 0, 0)
    assert exponents(r["mGy h^-1"]) == (2, 0, -3, 0, 0, 0, 0, 0)


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
    assert any("exponents (0, 0, -1, 0, 0, 0, 0, 0) differ from AbsorbedDoseRate (2, 0, -3" in p
               for p in bridge_problems(planted))


def test_roentgen_per_hour_mapped_to_the_dose_rate_fails_on_arithmetic(rows):
    planted = _plant(rows, "R h^-1", "unit", quantity_kind="AbsorbedDoseRate")
    assert any("exponents (0, -1, 0, 1, 0, 0, 0, 0) differ from AbsorbedDoseRate" in p
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


def test_lattice_system_and_the_gate(rows):
    lattice = {r["local_term"] for r in rows if r["axis"] == "unit" and r["unit_system"] == "lattice"}
    assert {"latt", "MCS", "H", "latt^2/field_step", "placeholder", "normalised"} <= lattice
    conv = {r["local_term"]: r for r in rows if r["axis"] == "conversion"}
    assert conv["latt->m"]["conversion_status"] == "blocked" and conv["latt->m"]["requirement_id"] == "D-PITCH"
    assert conv["MCS->s"]["depends_on"] == "latt->m"
    assert conv["latt->cm (synthetic dosimetry)"]["unit_system"] == "synthetic_reference"
    assert requirement_status()["D-PITCH"] == "awaiting_measurement"


def test_step4_controls_fire(rows, tmp_path):
    ready_no_factor = _plant(rows, "latt->m", "conversion", conversion_status="ready")
    problems = bridge_problems(ready_no_factor)
    assert any("ready with no factor" in p for p in problems)
    assert any("ready while 'D-PITCH' is 'awaiting_measurement'" in p for p in problems)

    dosed = _plant(rows, "placeholder", "unit", quantity_kind="AbsorbedDoseRate")
    assert any("a placeholder may not acquire a quantity kind" in p for p in bridge_problems(dosed))

    # The register flips in a temporary copy: D-PITCH measured, the bridge unchanged.
    text = REQUIREMENTS.read_text(encoding="utf-8")
    assert text.count("D-PITCH,lattice_pitch_um,") == 1
    flipped = tmp_path / "reference_d_requirements.csv"
    body = [l for l in text.split("\n") if l.startswith("D-PITCH,")][0]
    assert ",awaiting_measurement," in body
    flipped.write_text(text.replace(body, body.replace(",awaiting_measurement,", ",measured,", 1)), encoding="utf-8")
    reg = requirement_status(flipped)
    assert reg["D-PITCH"] == "measured", "mutation did not reach the parsed column"
    assert any("blocked on 'D-PITCH', whose status is 'measured'" in p for p in bridge_problems(rows, reg))
    assert bridge_problems(rows) == []

    clock_first = _plant(rows, "MCS->s", "conversion", conversion_status="ready", conversion_factor="1 s")
    assert any("ready before 'latt->m' is" in p for p in bridge_problems(clock_first))


# --- step 5: the tabulated prior and the shipped coefficient -----------------

def code_vectors(text=None) -> dict:
    """{symbol: {species tag: value}} parsed from CPMParams' commented vectors."""
    import re
    text = POTTS.read_text(encoding="utf-8") if text is None else text
    out = {}
    for field, symbol in (("β_ion", "beta_s_ion"), ("α_M_species", "alpha_M")):
        block = re.search(field + r"::Vector\{Float64\} = \[(.*?)\n\s*\]", text, re.S)
        assert block, f"{field} vector not found"
        vals = {tag: v for v, tag in re.findall(r"^\s*(-?[0-9.e+-]+),?\s*#\s*([A-Z]{2})\b", block.group(1), re.M)}
        assert set(vals) == set(SPECIES), (field, sorted(vals))
        out[symbol] = {k: float(v) for k, v in vals.items()}
    return out


def table_ranges() -> dict:
    """{(symbol, species name): (lo, hi)} from the species table."""
    import re
    out = {}
    for r in _read(SPECIES_TABLE):
        if r["symbol"] in ("beta_s_ion", "alpha_M"):
            lo, hi = re.split(r"\s+to\s+|(?<=\d)-(?=\d)", r["range"])
            out[(r["symbol"], r["species"])] = (float(lo), float(hi))
    assert len(out) == 9, sorted(out)
    return out


def substitution_problems(vectors, ranges) -> list:
    out = []
    for symbol, values in vectors.items():
        for tag, value in values.items():
            name = SPECIES[tag]
            rng = ranges.get((symbol, name)) if name else None
            if rng is None:
                if value != 0.0 and (symbol, tag) not in UNTABULATED_NONZERO:
                    out.append(f"{symbol}[{tag}] = {value} with no Table 2 row")
            elif not rng[0] <= abs(value) <= rng[1]:
                out.append(f"{symbol}[{tag}] = {value} outside Table 2's {rng}")
    return out


def test_coefficient_rows_pair_prior_with_shipped(rows):
    coef = {r["local_term"]: r for r in rows if r["axis"] == "coefficient"}
    assert set(coef) == {"beta_ion_prior", "beta_ion_cpm", "alpha_M_prior", "alpha_M_cpm", "melanin_coupling_cpm"}
    for r in coef.values():
        assert f"\\label{{{r['declared_in']}}}" in TEX.read_text(encoding="utf-8"), r["declared_in"]
    ids = {r["claim_id"] for r in _read(SPECIES_TABLE)} | {r["config_key"] for r in _read(PARAMETERS)}
    for r in coef.values():
        missing = set(r["ledger_rows"].split(";")) - ids
        assert not missing, (r["local_term"], missing)
    # The distinguishing fact: same number, different unit system.
    assert coef["beta_ion_cpm"]["unit_system"] != coef["beta_ion_prior"]["unit_system"]
    assert coef["melanin_coupling_cpm"]["relation"] == "hard_coded_replacement"
    text = POTTS.read_text(encoding="utf-8")
    assert text.count("0.5 * M_local") == 2, "the hard-coded coupling moved"


def test_shipped_numbers_are_the_tabulated_priors():
    vectors, ranges = code_vectors(), table_ranges()
    assert substitution_problems(vectors, ranges) == []
    assert vectors["beta_s_ion"]["CN"] < 0 and vectors["alpha_M"]["DR"] == 0.0   # sign by role; no melanin


def test_step5_controls_fire(rows):
    vectors, ranges = code_vectors(), table_ranges()
    drifted = {k: dict(v) for k, v in vectors.items()}; drifted["beta_s_ion"]["CN"] = -5e-3
    assert any("outside Table 2" in p for p in substitution_problems(drifted, ranges))
    invented = {k: dict(v) for k, v in vectors.items()}; invented["alpha_M"]["DR"] = 0.2
    assert any("no Table 2 row" in p for p in substitution_problems(invented, ranges))
    assert any("relation ''" in p for p in bridge_problems(_plant(rows, "beta_ion_cpm", "coefficient", relation="")))
    same_system = _plant(rows, "beta_ion_cpm", "coefficient", substitution_of="alpha_M_cpm")
    assert any("must substitute a tabulated SI prior" in p for p in bridge_problems(same_system))
    prior_substituting = _plant(rows, "beta_ion_prior", "coefficient", substitution_of="alpha_M_prior")
    assert any("substitutes nothing" in p for p in bridge_problems(prior_substituting))
    wrong_unit = _plant(rows, "beta_ion_cpm", "coefficient", unit="Gy^-1")
    assert any("kind or system disagrees" in p for p in bridge_problems(wrong_unit))


def test_the_hard_coded_coupling_can_only_be_declared(rows):
    coef = {r["local_term"]: r for r in rows if r["axis"] == "coefficient"}
    assert coef["melanin_coupling_cpm"]["basis"] == "declared"
    assert coef["beta_ion_prior"]["basis"] == "derived"     # the ledger's word for these six rows
    cited = _plant(rows, "melanin_coupling_cpm", "coefficient", basis="primary_literature")
    problems = bridge_problems(cited)
    assert any("has no basis but declared" in p for p in problems)
    assert any("disagrees with PP-T2-29" in p for p in problems)
    relabelled = _plant(rows, "beta_ion_prior", "coefficient", basis="declared")
    assert any("disagrees with PP-T2-12" in p for p in bridge_problems(relabelled))
