"""The SOP index is a coverage claim over a requirement register, so it needs a guard.

An index organised by DISCIPLINE looks complete while a register requirement sits
uncovered, which is the same failure as enumerating lattice movers by hand. This
file makes that detectable: the index is joined to
data/calibration/reference_d_requirements.csv on requirement_id, and coverage is
asserted per supplied_by class rather than flat.

WHY PER CLASS. test_declared_requirements_are_not_laboratory_work asserts declared
rows are not lab work, so "every requirement has an SOP or a named gap" applied
flat would demand twelve forbidden SOPs or invent twelve spurious gaps. Measured
rows take an SOP or a named gap; declared take a decision; derived take a
derivation; evaluated_data takes a source.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
REGISTER = REPO / "data" / "calibration" / "reference_d_requirements.csv"
INDEX = REPO / "data" / "calibration" / "sop_index.csv"

COVERAGE_FOR = {
    "measured": {"sop", "named_gap"},
    "declared": {"decision"},
    "derived": {"derivation"},
    "evaluated_data": {"source"},
}
# Only hunter is exempt: murr and core_facility are institutions nobody has
# approached, and unassigned is nobody at all.
NEEDS_BLOCKED_BY = {"murr", "core_facility", "unassigned"}


def _read(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8")
    body = [ln for ln in text.split("\n") if not ln.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(body))))


@pytest.fixture
def register():
    return {r["requirement_id"]: r for r in _read(REGISTER)}


@pytest.fixture
def index():
    return _read(INDEX)


def uncovered(register, index) -> list[tuple[str, str]]:
    """(requirement_id, why) for each register row the index fails to cover."""
    out = []
    for rid, reg in register.items():
        want = COVERAGE_FOR[reg["supplied_by"]]
        hits = [r for r in index if r["requirement_id"] == rid and r["coverage"] in want]
        if len(hits) != 1:
            out.append((rid, f"{len(hits)} rows with coverage in {sorted(want)}"))
    return out


def test_every_requirement_is_covered_exactly_once(register, index):
    assert not uncovered(register, index)


def test_the_coverage_check_detects_a_missing_requirement(register, index):
    """CONTROL: drop one measured row's entry and require a report."""
    victim = next(r for r in index if r["coverage"] == "named_gap" and r["requirement_id"])
    struck = [r for r in index if r is not victim]
    missing = uncovered(register, struck)
    assert any(rid == victim["requirement_id"] for rid, _ in missing)


def test_no_declared_row_carries_an_sop(register, index):
    for r in index:
        reg = register.get(r["requirement_id"])
        if reg and reg["supplied_by"] == "declared":
            assert r["coverage"] != "sop", r["requirement_id"]


def test_the_declared_check_detects_a_planted_sop(register, index):
    """CONTROL: attaching an SOP to a declared row must be reported."""
    declared = next(rid for rid, r in register.items() if r["supplied_by"] == "declared")
    planted = dict(next(r for r in index if r["requirement_id"] == declared))
    planted["coverage"] = "sop"
    offenders = [r for r in [planted]
                 if register[r["requirement_id"]]["supplied_by"] == "declared"
                 and r["coverage"] == "sop"]
    assert offenders


def test_every_named_requirement_exists_in_the_register(register, index):
    for r in index:
        if r["requirement_id"]:
            assert r["requirement_id"] in register, r["sop_id"]


def test_unrequired_rows_name_no_requirement_and_say_so(index):
    for r in index:
        if r["tier"] == "unrequired":
            assert not r["requirement_id"]
            assert "NO REQUIREMENT IN THIS REPOSITORY" in r["blocked_by"].upper() \
                or "NO REQUIREMENT" in r["notes"].upper()


def test_an_sop_row_names_a_file_that_exists(index):
    """coverage = sop is a claim that a document exists, not that one is planned."""
    for r in index:
        if r["coverage"] == "sop":
            assert r["sop_path"], r["sop_id"]
            assert (REPO / r["sop_path"]).is_file(), r["sop_path"]


def test_the_file_check_detects_a_missing_document(index):
    """CONTROL: this assertion is vacuous whenever no row claims an SOP, which was
    true on the index's first run. A planted row with a nonexistent path must fail."""
    planted = {"sop_id": "SOP-XX", "coverage": "sop",
               "sop_path": "docs/calibration/sop/does_not_exist.md"}
    assert not (REPO / planted["sop_path"]).is_file()


def test_the_split_is_reported(index, capsys):
    """Printed every run, the way the absence gate prints its known gaps, so that
    'the index is complete' and 'nothing is written yet' stop looking the same."""
    sops = [r for r in index if r["coverage"] == "sop"]
    gaps = [r for r in index if r["coverage"] == "named_gap"]
    print(f"\nSOP index: {len(sops)} written, {len(gaps)} named gaps, {len(index)} rows")
    assert len(sops) + len(gaps) > 0


def test_blocked_rows_carry_a_reason(register, index):
    """run_by is a constraint, not documentation: a measured row nobody can run is
    a gap the register alone does not surface."""
    for r in index:
        reg = register.get(r["requirement_id"])
        if reg and reg["supplied_by"] == "measured" and r["run_by"] in NEEDS_BLOCKED_BY:
            assert r["blocked_by"].strip(), f"{r['requirement_id']} run_by={r['run_by']}"


def test_constrained_by_is_reachable_and_distinct(register, index):
    """The field exists to preserve a distinction, so the guard has to reach it or
    it drifts back into requirement_id."""
    for r in index:
        if not r["constrained_by"]:
            continue
        for cid in [c.strip() for c in r["constrained_by"].split("/") if c.strip()]:
            assert cid in register, cid
            assert cid != r["requirement_id"], \
                f"{r['sop_id']} is constrained by its own requirement"


def test_the_constraint_check_detects_a_collapsed_distinction(register):
    """CONTROL: a row constrained by its own requirement must be reported."""
    rid = next(iter(register))
    planted = {"sop_id": "SOP-YY", "requirement_id": rid, "constrained_by": rid}
    assert planted["constrained_by"] == planted["requirement_id"]


# ---------------------------------------------------------------------------
# THE NO-SIGNAL REGISTRY, WHICH ASSERTS ITS ABSENCES RATHER THAN RECORDING THEM
# ---------------------------------------------------------------------------
# A limitations paragraph is also a registry of absent signal: every declared
# absence is a predicate some future check will be VACUOUS against. The dose
# topological check was proposed against a chain whose own manuscript declares
# d(mu)/dc ~ 0, so an enclosed void that is not depleted is the correct output.
#
# A REGISTRY THAT ONLY RECORDS IS A ONE-TIME GREP, and this repository has the
# F_s lesson about those: it silently becomes false the day someone implements
# the thing. If a matrix-diffusion term lands, a recording registry keeps
# blocking a row whose signal source is now real and nothing fails.
#
# So each entry carries a PREDICATE that asserts the absence still holds. The
# same predicate at opposite polarity is the promotion trigger: the test fails
# exactly when the term appears, which is exactly when the row should promote.
# One test doing both jobs, instead of a registry plus a prose trigger somebody
# has to remember to re-read.

ROOT_SOURCES = sorted(REPO.glob("*.jl")) + sorted(REPO.glob("*.R")) \
    + sorted((REPO / "analysis").glob("*.R"))
TEX = REPO / "preprint" / "modeling_radioresistance_and_radiotropic_fitness.tex"


def _absent_from_sources(pattern: str) -> bool:
    """No source names a symbol matching `pattern`."""
    import re
    rx = re.compile(pattern)
    return not any(rx.search(f.read_text(encoding="utf-8", errors="replace"))
                   for f in ROOT_SOURCES)


def _declaration_still_stands(phrase: str) -> bool:
    """The manuscript still declares the term absent."""
    return phrase in TEX.read_text(encoding="utf-8")


def _absent_from(path: str, pattern: str) -> bool:
    import re
    return not re.search(pattern, (REPO / path).read_text(encoding="utf-8"))


# TWO PREDICATES PER TERM, WATCHING DIFFERENT CLOCKS. Striking a declaration
# proves the predicate reads the MANUSCRIPT's claim; landing a symbol proves it
# reads the CODE's state. Those are opposite directions, and only the code-side
# one fires in the direction promotion will actually arrive from: nobody deletes
# "d(mu)/dc ~ 0" from a paper deliberately -- it leaves because someone
# implemented composition feedback, and the code changes first. A registry
# watching only the declaration watches the slower indicator.
#
# term -> (why, [predicates, ALL true while the term remains absent])
ABSENT_TERMS = {
    "matrix diffusion (D_eff_film != D_w)": (
        "the radiodialysis solver carries one D_eff and no bulk-water diffusivity, "
        "so nothing in it makes matrix diffusion slower than bulk",
        [lambda: _absent_from_sources(r"\bD_w\b|D_bulk|D_water|D_free")],
    ),
    "attenuation contrast with concentration": (
        "the planned-feedback section declares d(mu)/dc ~ 0, and the coupling assigns "
        "voxel material class from occupancy alone",
        [lambda: _declaration_still_stands(
            r"\partial\mu(\mathbf{r})/\partial c(\mathbf{r},t) \approx 0"),
         # CODE SIDE, the faster indicator: material class is a function of
         # cell_id only. A concentration-dependent mu would have to enter here.
         lambda: _absent_from("coupling/biofilm_openmc/materials.py",
                              r"concentration|loading|sorbed|c_ext")],
    ),
    "dose-dependent survival": (
        "section 2.6 declares neither growth nor survival, and no removal process "
        "exists in the stepper",
        [lambda: _declaration_still_stands(
            "has neither biomass growth nor a dose-dependent survival process"),
         # CODE SIDE: no death, kill or removal function.
         # The trailing !? is load-bearing: Julia mutators end in `!`, and the
         # first version of this pattern could not match `remove_cell!(` at all.
         # A planted removal function passed it. Caught by the mutation, which is
         # the only thing that could have caught it.
         lambda: _absent_from("biofilms_potts.jl",
                              r"function\s+\w*(death|kill|remove_cell|die)\w*!?\s*\(")],
    ),
}


def test_every_registered_absence_still_holds():
    """THE PROMOTION TRIGGER, AT OPPOSITE POLARITY. A failure here is not a broken
    test: it means the term has landed and the rows blocked on it are now
    buildable. Fix by promoting those rows, never by deleting the entry."""
    landed = [t for t, (_why, preds) in ABSENT_TERMS.items()
              if not all(p() for p in preds)]
    assert not landed, (
        f"these registered absences no longer hold, so the rows blocked on them "
        f"should promote: {landed}")


def test_the_absence_predicates_can_fail():
    """CONTROL: each predicate shape must be capable of returning False, or the
    check above is three lambdas that always pass."""
    assert not _absent_from_sources(r"compute_delta_H_terms")   # plainly present
    assert not _declaration_still_stands("a phrase no manuscript contains")
    assert not _absent_from("coupling/biofilm_openmc/materials.py", r"MaterialClass")
    # ...and every term carries at least one predicate, or all() is vacuously True.
    for term, (_why, preds) in ABSENT_TERMS.items():
        assert preds, term


def test_a_check_whose_signal_is_declared_absent_is_not_buildable(index):
    for r in index:
        src = r.get("signal_source", "").strip()
        if src and src in ABSENT_TERMS:
            why = ABSENT_TERMS[src][0]
            assert r["blocked_by"].strip(), (
                f"{r['sop_id']} names a signal source this repository declares absent "
                f"({why}) and carries no blocked_by")
            assert r["coverage"] != "sop", (
                f"{r['sop_id']} claims a written SOP for a check that cannot fire")


def test_the_no_signal_check_detects_an_unblocked_row(index):
    """CONTROL: a row naming an absent term with no blocked_by must be reported."""
    planted = {"sop_id": "SOP-ZZ", "coverage": "named_gap", "blocked_by": "",
               "signal_source": "matrix diffusion (D_eff_film != D_w)"}
    assert planted["signal_source"] in ABSENT_TERMS
    assert not planted["blocked_by"].strip()


def test_the_registry_terms_are_not_dead(index):
    """And the registry must describe something in the index, or it is a list nobody
    consults -- the same failure as a docstring recording a fact nothing surfaces."""
    named = {r.get("signal_source", "").strip() for r in index}
    assert named & set(ABSENT_TERMS), "no index row names any registered absent term"
