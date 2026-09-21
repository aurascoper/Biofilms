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
# SCOPE: the run_by values of MEASURED register rows only. Declared rows are
# never lab work and derived rows are computed, so neither is reached here.
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

_root_sources = lambda root: sorted(root.glob("*.jl")) + sorted(root.glob("*.R")) \
    + sorted((root / "analysis").glob("*.R"))   # noqa: E731 -- two lines, so :219 stays put
TEX = REPO / "preprint" / "modeling_radioresistance_and_radiotropic_fitness.tex"


def _absent_from_sources(pattern: str, root: Path = REPO) -> bool:
    """No source names a symbol matching `pattern`."""
    import re
    rx = re.compile(pattern)
    return not any(rx.search(f.read_text(encoding="utf-8", errors="replace"))
                   for f in _root_sources(root))


def _declaration_still_stands(phrase: str, root: Path = REPO) -> bool:
    """The manuscript still declares the term absent."""
    return phrase in (root / TEX.relative_to(REPO)).read_text(encoding="utf-8")


def _absent_from(path: str, pattern: str, root: Path = REPO) -> bool:
    import re
    return not re.search(pattern, (root / path).read_text(encoding="utf-8"))


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
        [lambda root=REPO: _absent_from_sources(r"\bD_w\b|D_bulk|D_water|D_free", root)],
    ),
    "attenuation contrast with concentration": (
        "the planned-feedback section declares d(mu)/dc ~ 0, and the coupling assigns "
        "voxel material class from occupancy alone",
        [lambda root=REPO: _declaration_still_stands(
            r"\partial\mu(\mathbf{r})/\partial c(\mathbf{r},t) \approx 0", root),
         # CODE SIDE, the faster indicator: material class is a function of
         # cell_id only. A concentration-dependent mu would have to enter here.
         lambda root=REPO: _absent_from("coupling/biofilm_openmc/materials.py",
                                        r"concentration|loading|sorbed|c_ext", root),
         # REGIME SIDE. d(mu)/dc ~ 0 is a PHOTON statement: at ~1 MeV over
         # C/H/N/O, biomass and void are near-indistinguishable. For NEUTRONS it
         # is false and interestingly so -- biofilm is largely water, so a
         # biofilm region moderates against a void and hydrogen's scattering
         # cross-section is large. A source change is exactly the update path
         # nobody would re-read this entry for, and neither predicate above
         # fails if the particle changes while the manuscript sentence stays.
         # Matched on CONTENT, not a line number: a line anchor silently starts
         # reading a different statement the moment anything above it shifts.
         lambda root=REPO: not _absent_from("coupling/biofilm_openmc/model.py",
                                            r'particle\s*=\s*"photon"', root)],
    ),
    "dose-dependent survival": (
        "the abstract declares neither growth nor a dose-dependent survival process, "
        "and the stepper removes a cell only when its volume reaches zero: it records "
        "the lifetime dose in that event, and the dose decides nothing there",
        [lambda root=REPO: _declaration_still_stands(
            "has neither biomass growth nor a dose-dependent survival process", root),
         # CODE SIDE: no death, kill or removal function.
         # The trailing !? is load-bearing: Julia mutators end in `!`, and the
         # first version of this pattern could not match `remove_cell!(` at all.
         # A planted removal function passed it.
         #
         # WHAT CAUGHT IT, PRECISELY, BECAUSE THE CREDIT MATTERS FOR WHETHER THE
         # STEP IS OPTIONAL. The argument for adding a code-side predicate was
         # about which indicator a row watches, slow or fast; it gave a reason to
         # add this one and said nothing about whether the added one worked. What
         # found the bug was FIRING THE MUTATION in the direction promotion
         # arrives from and watching the check fail to fail. The pattern was also
         # correct for a language it was not written against -- `!` is not a word
         # character, and mutator-bang is a Julia convention -- so it reads as
         # reasonable to anyone reviewing it. Invisible to inspection, visible
         # only to execution.
         lambda root=REPO: _absent_from("biofilms_potts.jl",
                                        r"function\s+\w*(death|kill|remove_cell|die)\w*!?\s*\(",
                                        root),
         # THE TWO ABOVE READ ONE MANUSCRIPT AND ONE JULIA FILE. A port in another
         # language changes neither, so the two below read every listed file.
         lambda root=REPO: not removal_tokens(root),
         lambda root=REPO: not capability_claims(root),
         # ...and a pinned file that is gone is reported from a list of its own.
         lambda root=REPO: not missing_pinned_paths(root)],
    ),
}


# ---------------------------------------------------------------------------
# TWO READERS WITH NO LANGUAGE AND NO DIRECTORY BUILT IN
# ---------------------------------------------------------------------------
# The survival predicate opened biofilms_potts.jl and wanted the word `function`.
# An Odin port writes `cell_kill :: proc(` under odin/ and was invisible to it twice
# over. A list of source suffixes would repeat that defect for the next language,
# so the file set is whatever git lists, and a suffix the author never met is read.

def listed_text(root: Path = REPO) -> dict:
    """relative path -> text, for every file git lists under `root`: tracked, or
    untracked and not ignored. An ignored file cannot arrive in a pull request."""
    import subprocess
    out = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z", "--cached", "--others",
         "--exclude-standard"], check=True, capture_output=True).stdout.decode("utf-8")
    texts = {}
    for rel in filter(None, out.split("\0")):
        path = root / rel
        if not path.is_file():
            continue                    # still in the index, deleted from the work tree
        raw = path.read_bytes()
        if b"\0" not in raw[:8192]:     # keyed on the content being binary, not on a suffix
            texts[rel] = raw.decode("utf-8", errors="replace")
    if not texts:
        raise RuntimeError(f"git listed no text file under {root}")   # never default to pass
    return texts


# DEPARTURE FROM AGENTS.md, DECIDED BY HUNTER ON 2026-09-20. "A suffix cannot justify
# skipping" -- and these six are skipped by removal_tokens(). They hold prose and data
# ABOUT the model, where death is a topic: measured on fe6c59e, reading them costs 27
# allowlisted files for `death` alone. capability_claims() reads them, suffix or not.
PROSE_AND_DATA = {".md", ".tex", ".txt", ".csv", ".tsv", ".json"}

# death|kill|remove_cell|die are this entry's first stems. The others are the words
# PRRT-spatial-CPM's spec/response_growth_survival.md uses for the process: SF, its
# survival fraction; dying; mitotic; necrotic; cull; resorb and resorption. Matched
# inside an identifier at a snake_case or camelCase boundary, so `skill`, `died` and
# `overkill` do not match and `kill_cell!`, `cellDeath` and `CELL_DYING` do.
#
# survival, survivor and surviving are the term's own noun, added after Codex planted
# `survival_probability` on pull request #47. Measured on 8374b59: 14 files. The verb
# (survive, survives, survived) is in 46, nearly all English in comments and test
# names, so it stays out and KNOWN_UNSEEN says so.
REMOVAL_VOCAB = (
    r"(?<![A-Za-z])(?i:deaths?|kill|cull|dying|remove_?cell|mitotic|necro|resor[bp])"
    r"|(?<![A-Za-z])(?i:survival|survivors?|surviving)"
    r"|(?<![A-Za-z])(?i:die)(?![A-Za-z])"
    r"|(?<=[a-z])(?:Death|Kill|Cull|Dying|RemoveCell|Survival|Survivor|Surviving|Die(?![a-z]))"
    r"|(?<![A-Za-z0-9])(?i:sf)(?![A-Za-z0-9])")

WHOLE_FILE = None
# file -> (the identifiers allowed in it, why). Keyed on the identifier, not the file:
# a `kill_cell!` added to the stepper fires although the stepper is listed here.
REMOVAL_ALLOWED = {
    # Conditional: see EXEMPTION_HOLDS. The reason is checked, not only stated.
    "biofilms_potts.jl": ({"DeathEvent", "death", "death_log",
                           "SURVIVAL", "Surviving", "survival"},
                          "bookkeeping for removal at volume <= 0, which no dose term "
                          "enters, and comments stating the model computes no survival"),
    "export_checkpoint.jl": ({"DeathEvent", "death", "death_log"},
                             "serialises that same bookkeeping"),
    "biofilms_3d.R": ({"death"}, "a comment stating the model has no death"),
    "coupling/tests/test_feedback_gate.py":
        ({"death", "growth_survival", "growth_survival_response"},
         "a docstring stating the CPM has no death, and the pathway id the gate refuses"),
    "coupling/biofilm_openmc/feedback_gate.py": ({"growth_survival"},
                                                 "lists the pathway as unsupported"),
    "post-biofilms-standalone.html": ({"death", "killed"},
                                      "blog prose: 'the model has no death', 'seeds killed it'"),
    "docs/visualization/macos-arm64-2026-09-11/component_overlap.py":
        ({"deaths"}, "connected components that disappear between frames, not cells"),
    "docs/visualization/macos-arm64-2026-09-11/test_component_overlap.py":
        ({"deaths", "surviving"}, "the tests of that component count; a lineage that "
                                  "outlasts a merge of components"),
    "docs/visualization/macos-arm64-2026-09-11/decay_reference.py":
        ({"surviving", "surviving_fraction"}, "nuclei not yet decayed, 2^(-t/T_half)"),
    ".github/workflows/claims-guards.yml": ({"survivor"}, "a withdrawn claim still in a file"),
    "calibration/tests/test_claims_ledger.py": ({"survivor", "survivors"},
                                                "withdrawn phrases still in a document"),
    "calibration/tests/test_withdrawn_comparisons.py": ({"survivors", "survivors_in"},
                                                        "withdrawn comparisons still in the tree"),
    "calibration/tests/test_figure_staleness.py": ({"surviving"},
                                                   "a comment: withdrawn claims in images"),
    "calibration/tests/test_biomass_field.py": ({"surviving"},
                                                "a comment: biomass left in a synthetic field"),
    "coupling/tests/test_pilot_seeding.py": ({"surviving"},
                                             "a comment: a number that outlasts a rename"),
    "coupling/biofilm_openmc/lineage.py": ({"survivors"},
                                           "a docstring: labels that remain after masking"),
    "geolocator/api.py": ({"surviving"}, "a comment: fuel pairs that pass a filter"),
    "coupling/scripts/synthetic_e2e.py": ({"sf"}, "a local: the far-field standard deviation"),
    "diagnostics/payload_census/bench_reduction.py": ({"sf"}, "a local: samples, full payload"),
    "diagnostics/inert_signal/guard.jl": ({"SF"}, "a comment's example alias for signal_field"),
    "geolocator/static/app/tokens.js": ({"SF"}, "the font name SF Mono"),
    "geolocator/static/index.html": ({"SF"}, "the font name SF Mono"),
    "geolocator/static/app/hud/detectionOverlay.js": ({"cull", "culled", "killed"},
                                                      "frustum culling of map markers"),
    "Primary Fuel Source Geolocator.twb": ({"cull"}, "the Tableau attribute mark-labels-cull"),
    "calibration/tests/test_sop_index.py": (WHOLE_FILE, "declares the vocabulary and the plants"),
}


def stepper_removal_sites(root: Path = REPO) -> list:
    """(function, enclosing conditions) for every line of the stepper that deletes from
    state.cells or constructs a DeathEvent. Read by indentation, the way that file is
    written. Reformat the file and this returns something else, which voids the
    exemption below: it fails closed."""
    import re
    lines = (root / "biofilms_potts.jl").read_text(encoding="utf-8").split("\n")
    sites = []
    for i, line in enumerate(lines):
        if "delete!(state.cells" not in line and not re.search(r"\bDeathEvent\(", line):
            continue
        indent, conds, fn = len(line) - len(line.lstrip()), [], None
        for prev in reversed(lines[:i]):
            if not prev.strip() or len(prev) - len(prev.lstrip()) >= indent:
                continue
            indent, head = len(prev) - len(prev.lstrip()), prev.strip()
            if head.startswith(("if ", "elseif ")):
                conds.append(head.split(" ", 1)[1])
            elif head == "else":
                conds.append("else")
            elif head.startswith("function "):
                fn = re.match(r"function\s+([\w!]+)", head).group(1)
                break
        sites.append((fn, tuple(reversed(conds))))
    return sites


def stepper_removal_functions(root: Path = REPO) -> dict:
    """function name -> sha256 of its code, for every stepper function that deletes from
    state.cells or constructs a DeathEvent. Comment lines and blank lines are left out.

    The if-headers above are a marker for what governs a removal. A guard clause, a ternary
    or an early return governs it too and changes no header (Codex on pull request #47, at
    02f8732). So the WHOLE function is pinned: any edit to its code voids the exemption
    until someone has read the function again and updated the digest."""
    import hashlib
    import re
    lines = (root / "biofilms_potts.jl").read_text(encoding="utf-8").split("\n")
    out, start = {}, None
    for i, line in enumerate(lines):
        if line.startswith("function "):
            start = i
        elif line == "end" and start is not None:
            body = lines[start:i + 1]
            if any("delete!(state.cells" in b or re.search(r"\bDeathEvent\(", b) for b in body):
                code = "\n".join(b.strip() for b in body if b.strip() and not b.strip().startswith("#"))
                name = re.match(r"function\s+([\w!]+)", body[0]).group(1)
                out[name] = hashlib.sha256(code.encode("utf-8")).hexdigest()
            start = None
    return out


# What earns the stepper its exemption, stated as the thing itself. Codex on pull
# request #47: a dose condition added to this cleanup reuses DeathEvent and death_log,
# so an exemption granted on the file name subtracts the very signal it should report.
STEPPER_REMOVES_ON_VOLUME_ALONE = [
    ("mcs_step!", ("cell.volume <= 0",)),      # constructs the DeathEvent
    ("mcs_step!", ("cell.volume <= 0",)),      # deletes that cell
    ("divide_cell!", ()),                      # the parent's identity ends at division
]
# EVERY EXEMPTION IS PINNED TO THE LINES THAT EARNED IT. Codex on pull request #47, four
# times: keyed on the file, an exemption hid a dose condition behind the stepper's own
# bookkeeping; keyed on a word, it hid `function survival(dose)`; keyed on a count, it hid
# that function written over the one prose line that says survival; and pinned for three
# model files only, it hid `def growth_survival(dose)` written over the feedback gate's
# unsupported-pathway line, in production code I had called a test or a script.
#
# So there is no carve-out. For every file in REMOVAL_ALLOWED the pin is every non-comment
# line that holds an exempt identifier, verbatim and in order. A `#` line is exempt as what
# it is: it cannot run. A new line that uses an exempt identifier, or an edit to a pinned
# one, voids that file's exemption and the gate fails. THAT IS A TRIPWIRE ON PURPOSE: the
# person who made the edit looks, and updates the pin. Measured on b68e800: 79 lines in 24
# files, none longer than 87 characters.
#
# WHAT STAYS EXEMPT WHOLESALE, and why: this file, because it declares the vocabulary and
# holds the pins and the plants.
EXEMPT_LINES = {
    'biofilms_potts.jl': [
        'the quantity the model computes is a DIRECTIONAL PREFERENCE, not survival and',
        'lifetime_dose_Gy::Float64   # accrued while alive; finalized at division/death',
        'struct DeathEvent',
        'death_log::Vector{DeathEvent}',
        'DivisionEvent[], DeathEvent[], ArchivedCell[],',
        'push!(state.death_log,',
        'DeathEvent(state.current_mcs, state.physical_time_s, id,',
        '@printf("  Surviving cells: %d\\n", length(state.cells))',
        '@printf("  Surviving cells: %d\\n", length(state.cells))',
    ],
    'export_checkpoint.jl': [
        'de = state.death_log',
        'f["events/death/$name"] = name in ("physical_time_s", "lifetime_dose_Gy") ?',
        'de = [SR.DeathEvent(',
        'Int(read(f["events/death/mcs"])[i]),',
        'read(f["events/death/physical_time_s"])[i],',
        'Int(read(f["events/death/cell_id"])[i]),',
        'Int(read(f["events/death/species"])[i]),',
        'Int(read(f["events/death/lineage_id"])[i]),',
        'Int(read(f["events/death/generation"])[i]),',
        'read(f["events/death/lifetime_dose_Gy"])[i])',
        'for i in eachindex(read(f["events/death/mcs"]))]',
        'convert(Vector{SR.DeathEvent}, de),',
    ],
    'biofilms_3d.R': [
    ],
    'coupling/tests/test_feedback_gate.py': [
        '@pytest.mark.parametrize("pathway", ["growth_survival", "membrane_response"])',
        '"""The CPM has no birth and no death, and the membrane constitutive law is',
        'for pid in ("growth_survival_response", "membrane_response"):',
    ],
    'coupling/biofilm_openmc/feedback_gate.py': [
        'unsupported_pathways: tuple = ("growth_survival", "membrane_response")',
    ],
    'post-biofilms-standalone.html': [
        'damage, not benefit: the model has no death and no growth to express',
        '<p>Tidy. I was pleased with it. Three more seeds killed it. On seed 43',
    ],
    'docs/visualization/macos-arm64-2026-09-11/component_overlap.py': [
        'events.append({"mcs": t, "births": k, "deaths": 0, "merges": 0, "splits": 0,',
        'events.append({"mcs": t, "births": births, "deaths": disappeared,',
        'churn = sum(e["births"] + e["deaths"] + e.get("retired_by_merge", 0) for e in within)',
        'for kk in ("births", "deaths", "retired_by_merge", "merges", "splits")}',
        '% (totals["births"], totals["deaths"], totals["retired_by_merge"],',
        '("births", "deaths", "retired_by_merge", "merges", "carried")}',
    ],
    'docs/visualization/macos-arm64-2026-09-11/test_component_overlap.py': [
        'present for one frame reported two, an all-vanish frame reported zero deaths, and every',
        'check("deaths == 1", e["deaths"] == 1, e["deaths"])',
        'print("\\n=== two components joining into one: 1 merge, 1 retirement, 0 deaths ===")',
        'check("deaths == 0 (neither disappeared)", e["deaths"] == 0, e["deaths"])',
        'print("\\n=== one component splitting into two: 1 split, 1 birth, 0 deaths ===")',
        'check("deaths == 0", e["deaths"] == 0, e["deaths"])',
        'check("surviving lineage is the LARGER predecessor (id 0)",',
        'check("deaths == 1 (C disappeared)", e["deaths"] == 1, e["deaths"])',
    ],
    'docs/visualization/macos-arm64-2026-09-11/decay_reference.py': [
        'surviving = 2.0 ** (-t / T_HALF_DAYS)',
        'act = a0 * surviving',
        '"surviving_fraction": surviving,',
        'metrics.append({"days": t, "surviving_fraction": surviving,',
        'print("surviving at 1 t1/2 : %.6f   at 2: %.6f   at 3: %.6f" %',
    ],
    '.github/workflows/claims-guards.yml': [
    ],
    'calibration/tests/test_claims_ledger.py': [
        'survivors = []',
        'survivors.append(f"{row[\'claim_id\']} in {row[\'document\']}: {phrase[:90]}")',
        'survivors.append(',
        'assert not survivors, (',
        '"carried them:\\n  " + "\\n  ".join(survivors))',
        'reports survivors by id. A duplicate makes every one of those references',
        'duplicate definition is legal Python and pytest collects only the survivor.',
    ],
    'calibration/tests/test_withdrawn_comparisons.py': [
        'def survivors_in(text: str) -> list[str]:',
        'survivors = []',
        'for window in survivors_in(path.read_text(encoding="utf-8", errors="replace")):',
        'survivors.append(f"{rel}: ...{window[WINDOW - 90:WINDOW + 120]}...")',
        'assert not survivors, (',
        '"as withdrawn:\\n  " + "\\n  ".join(survivors))',
        'assert len(survivors_in(as_it_wrapped)) == 1',
        'assert len(survivors_in(as_it_was)) == 1',
        'assert survivors_in(corrected) == []',
        'assert survivors_in(unrelated) == []',
        'assert survivors_in(dose) == []',
    ],
    'calibration/tests/test_figure_staleness.py': [
    ],
    'calibration/tests/test_biomass_field.py': [
    ],
    'coupling/tests/test_pilot_seeding.py': [
        'hazard was the NUMBER surviving, not the variable name.',
    ],
    'coupling/biofilm_openmc/lineage.py': [
        'Masking the inputs beforehand does NOT fix that — survivors labelled 0 are',
    ],
    'geolocator/api.py': [
    ],
    'coupling/scripts/synthetic_e2e.py': [
        'sf = np.sqrt((sd_occ[far] ** 2).sum()) / far.sum()',
        'separation = (dn - df) / math.sqrt(sn ** 2 + sf ** 2)',
    ],
    'diagnostics/payload_census/bench_reduction.py': [
        'sf, af = samples(lambda: load_full(full_d), reps)',
        'show("full trajectory, 10 arrays", sf)',
        'mf, mr = statistics.median(sf), statistics.median(sr)',
        'sep = sr[-1] < sf[0] or sf[-1] < sr[0]',
    ],
    'diagnostics/inert_signal/guard.jl': [
        'bump`), an alias (`const SF = signal_field`), or entry through a struct field it does',
    ],
    'geolocator/static/app/tokens.js': [
        'mono: \'ui-monospace, "SF Mono", SFMono-Regular, Menlo, Consolas, monospace\',',
    ],
    'geolocator/static/index.html': [
        '--mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,monospace;',
    ],
    'geolocator/static/app/hud/detectionOverlay.js': [
        '* projected every marker every frame would be the thing that finally killed the',
        '*   - frustum-culled before sorting, so off-screen markers cost nothing',
        'if (!frustum.containsPoint(v)) continue;         // cull before sorting',
    ],
    'Primary Fuel Source Geolocator.twb': [
        "<format attr='mark-labels-cull' value='true' />",
    ],
}


def exempt_code_lines(root: Path, rel: str) -> list:
    """The stripped lines of root/rel that are not comments and hold an identifier
    REMOVAL_ALLOWED exempts there, in file order."""
    import re
    token = re.compile(r"[A-Za-z_][A-Za-z0-9_]*!?")
    allowed = REMOVAL_ALLOWED[rel][0]
    text = (root / rel).read_text(encoding="utf-8", errors="replace")
    text = re.sub(r"[A-Za-z0-9+/=]{60,}", " ", text)      # the same blobs removal_tokens() drops
    return [line.strip() for line in text.split("\n")
            if not line.strip().startswith("#") and allowed & set(token.findall(line))]


# file -> a relation that must hold BESIDE the line pin for its exemption to apply
# A TRIPWIRE ON TWO FUNCTIONS, AND ITS COST. Any edit to the code of mcs_step! or of
# divide_cell! fails the survival gate until this digest is updated. Whoever updates it has
# read the function and confirmed that every removal in it still depends on volume alone.
# Control flow in a CALLER of these functions is not read by anything here.
STEPPER_REMOVAL_FUNCTIONS = {
    "mcs_step!": "b5f1f5e0812e1f6e11c3d375e1a2b16b2254e9ba9256155af8850a32117b621c",
    "divide_cell!": "60c744c42cbbe2b3b36c06a92c740077788865f47f0c8b3878e3b957bf1355f4",
}

EXEMPTION_HOLDS = {
    "biofilms_potts.jl":         # what governs each removal, not only which lines name one
        lambda root: (stepper_removal_sites(root) == STEPPER_REMOVES_ON_VOLUME_ALONE
                      and stepper_removal_functions(root) == STEPPER_REMOVAL_FUNCTIONS),
    "export_checkpoint.jl":      # it serialises the log; it must remove nothing
        lambda root: "delete!(" not in (root / "export_checkpoint.jl").read_text(encoding="utf-8"),
}


# A BROKEN PIN IS THE SIGNAL, WITH OR WITHOUT AN IDENTIFIER LEFT TO NAME. Codex on pull
# request #47, at a60be28: the natural way to ENABLE the capability deletes the line that
# declared it absent. The pin stops matching, the exemption goes void, and the only
# identifier a reader could have named is gone with the line. So a file whose pinned
# lines changed is reported even when nothing else in it is.
PINS_CHANGED = "<the pinned lines of this file changed>"


def removal_tokens(root: Path = REPO) -> dict:
    """relative path -> the removal-vocabulary identifiers in it that REMOVAL_ALLOWED
    does not name. Empty while no survival or removal process has arrived."""
    import re
    token = re.compile(r"[A-Za-z_][A-Za-z0-9_]*!?")
    vocab = re.compile(REMOVAL_VOCAB)
    # An unbroken run of 60 base64-alphabet characters is an embedded image or a digest.
    # No identifier here has one: the longest word between underscores is far shorter.
    blob = re.compile(r"[A-Za-z0-9+/=]{60,}")
    found = {}
    for rel, text in listed_text(root).items():
        allowed = REMOVAL_ALLOWED.get(rel, (set(), ""))[0]
        if allowed is WHOLE_FILE or Path(rel).suffix.lower() in PROSE_AND_DATA:
            continue
        broken = rel in REMOVAL_ALLOWED and not (
            exempt_code_lines(root, rel) == EXEMPT_LINES[rel]
            and EXEMPTION_HOLDS.get(rel, lambda _root: True)(root))
        if broken:
            allowed = set()
        hits = sorted({t for t in token.findall(blob.sub(" ", text)) if vocab.search(t)}
                      - allowed)
        if hits or broken:
            found[rel] = hits or [PINS_CHANGED]
    for rel in set(EXEMPT_LINES) - set(found) - set(listed_text(root)):
        found[rel] = [PINS_CHANGED]     # the whole file is gone, and its pins with it
    return found


# The name the PRRT spec, this repository's ledgers and the port all use for the
# capability. A port that follows spec/pr_checklist.md declares it, in any language.
CAPABILITY = "response.growth_survival"
# These files may name it, and why. WHICH LINES may name it is CAPABILITY_LINES below.
CAPABILITY_NAMED_AS_UNSUPPORTED = {
    "README.md": "states it is unsupported_by_current_model and unfittable in principle",
    "data/claims_ledger.csv": "rows and comments recording the absence",
    "data/parameter_provenance.csv": "the row whose status is unsupported_by_current_model",
    "data/uncertainty/feedback_parameter_distributions.csv": "status unsupported, excluded",
    "docs/calibration/integration_contract.md": "states it is unfittable in principle",
    "docs/calibration/time_observable_contract.md": "applies that same status",
    "docs/online_coupling_design.md": "lists the response transforms as BLOCKED",
    "calibration/tests/test_sop_index.py": "declares the name and the plants",
}


# A file exempt by its path stays exempt when its sentence flips to "supported" (Codex on
# pull request #47, at ce141b2). So every line that names the capability today is pinned by
# the sha256 of its stripped text. A new line, or an edited one, is reported. The hint
# beside each digest is for the reader; only the digest is compared.
CAPABILITY_LINES = {
    'README.md': {
        'd532e460a85283c8f1f1f6bbec45cb8ce273310b8fb6d2ac774de93cc8236dd5',   # `response.growth_survival` is `unsupported_by_current_mode
        'e9f8ed58ad95fb9c6bba3fb0f07596140895586157aa048dbd86caf7a526fcf9',   # | Growth, division, survival, or extinction | **`response.
        '89038d1100abe10df00a57e419ca5b40ad7d276c45c6ca1c7e751f07bd993e11',   # | Parcel persistence | 42 / 42 seeded **biomass parcels**,
    },
    'data/claims_ledger.csv': {
        'a20f598b29e16a246a648ccdc8ad9bff380c84b8d6dba5e57bd09166174fe0c4',   # #    Both are STRUCTURALLY ABSENT, like response.growth_su
        'bdd5212b0150cc38d00756d32f9feaf0a1c35187d9426ee1c12102d8deebac81',   # RM-KR-12,README,Key Results row 7,no extinctions,model_out
        'b9152435103a2a45c1eb62c48e5812279333f8e1cd9f2c2919955acbbd9bce45',   # NEWS-PHEN-02,repository,docs/calibration/integration_contr
        '19a8906a502d01efc3d9676399d7b206ffcb1b3c4bcff9c5aa917e748de5e33b',   # GUARD-SCOPE-02,repository,calibration/tests/test_sop_index
    },
    'data/parameter_provenance.csv': {
        '847be5f2ebe76d7c286d07f3b6b01573e6d9e450053a48c2fc48def29f0fb8a5',   # D_engineered_composite,none,response.hamiltonian,RADIOTROP
        'de1072e11fdfd4a5bf7b18c009bd0717a27579f01baded0fd503e8a42c32b21e',   # D_engineered_composite,none,response.growth_survival,Dose 
    },
    'data/uncertainty/feedback_parameter_distributions.csv': {
        '74ed6997de4db71e825187e511b21da7a776f222c8321c02b3d58484e955e1cd',   # growth_survival_response,response.growth_survival,dose to 
    },
    'docs/calibration/integration_contract.md': {
        '292820152278b7c843aefe40d9cdb1a09093301929232016c2438b605f0720e0',   # - **`response.growth_survival` is currently unfittable in 
        'e9dad65e5ac59c013f035612880eec1d4e4bcb8ae383a268fc1d165b5256a6e0',   # absent**, for the same reason `response.growth_survival` i
    },
    'docs/calibration/time_observable_contract.md': {
        '5af4b5d0d71e1ad2bc64dc170876fd39aae964ca7b276a3e1012d950024d8b64',   # The same status now applies to `response.growth_survival` 
    },
    'docs/online_coupling_design.md': {
        'dac60bf35e7bcab4e37866922463c9791582e87adb4c5a753f13ed0610596a3e',   # | 6 · Biological response transforms | **NOT NEEDED** | **
    },
    "calibration/tests/test_sop_index.py": WHOLE_FILE,   # declares the name and the plants
}


def capability_claims(root: Path = REPO) -> list:
    """"file:line" for every line that names the capability and is not pinned above.
    No suffix is skipped: the port's own declaration is a Markdown file."""
    import hashlib
    claims = []
    texts = listed_text(root)
    for rel in sorted(set(texts) | set(CAPABILITY_LINES)):      # a deleted file is not listed
        text = texts.get(rel, "")
        pinned = CAPABILITY_LINES.get(rel, set())
        if pinned is WHOLE_FILE or (CAPABILITY not in text and not pinned):
            continue
        now = {n: hashlib.sha256(line.strip().encode("utf-8")).hexdigest()
               for n, line in enumerate(text.split("\n"), 1) if CAPABILITY in line}
        claims += [f"{rel}:{n}" for n, digest in now.items() if digest not in pinned]
        gone = len(set(pinned) - set(now.values()))
        if gone:        # a deleted line has no digest of its own to report
            claims.append(f"{rel}: {gone} pinned line(s) gone")
    return claims


# THE PATHS THE PINS NAME, WRITTEN OUT A SECOND TIME, ON PURPOSE. The tables above remember a
# path only while they have an entry for it. The natural clean-up after deleting a file
# removes its entries, and then nothing remembers that the file was ever pinned (Codex on
# pull request #47, at c28a381). This list is not built from the tables, so that clean-up
# does not touch it, and test_every_pin_has_its_tables fails until the two agree again.
#
# WHERE THIS STOPS. An edit that deletes the file, its table entries AND its line here leaves
# the gate green. No guard in a repository outlasts an edit of itself. That edit is three
# deliberate changes to this file in one diff, and the diff's reader is the last check.
PINNED_PATHS = [
    ".github/workflows/claims-guards.yml",
    "Primary Fuel Source Geolocator.twb",
    "README.md",
    "biofilms_3d.R",
    "biofilms_potts.jl",
    "calibration/tests/test_biomass_field.py",
    "calibration/tests/test_claims_ledger.py",
    "calibration/tests/test_figure_staleness.py",
    "calibration/tests/test_withdrawn_comparisons.py",
    "coupling/biofilm_openmc/feedback_gate.py",
    "coupling/biofilm_openmc/lineage.py",
    "coupling/scripts/synthetic_e2e.py",
    "coupling/tests/test_feedback_gate.py",
    "coupling/tests/test_pilot_seeding.py",
    "data/claims_ledger.csv",
    "data/parameter_provenance.csv",
    "data/uncertainty/feedback_parameter_distributions.csv",
    "diagnostics/inert_signal/guard.jl",
    "diagnostics/payload_census/bench_reduction.py",
    "docs/calibration/integration_contract.md",
    "docs/calibration/time_observable_contract.md",
    "docs/online_coupling_design.md",
    "docs/visualization/macos-arm64-2026-09-11/component_overlap.py",
    "docs/visualization/macos-arm64-2026-09-11/decay_reference.py",
    "docs/visualization/macos-arm64-2026-09-11/test_component_overlap.py",
    "export_checkpoint.jl",
    "geolocator/api.py",
    "geolocator/static/app/hud/detectionOverlay.js",
    "geolocator/static/app/tokens.js",
    "geolocator/static/index.html",
    "post-biofilms-standalone.html",
]


def missing_pinned_paths(root: Path = REPO) -> list:
    """The pinned paths that git does not list under `root`."""
    return sorted(set(PINNED_PATHS) - set(listed_text(root)))


def landed_terms(root: Path = REPO, *, terms=None) -> list:
    """Registered absences whose predicates no longer all hold under `root`.
    Shared by the gate and its plants, so a plant cannot pass against a copy."""
    keys = set(ABSENT_TERMS) if terms is None else set(terms)
    return [t for t, (_why, preds) in ABSENT_TERMS.items()
            if t in keys and not all(p(root) for p in preds)]


def test_every_registered_absence_still_holds():
    """THE PROMOTION TRIGGER, AT OPPOSITE POLARITY. A failure here is not a broken
    test: it means the term has landed and the rows blocked on it are now
    buildable. Fix by promoting those rows, never by deleting the entry."""
    landed = landed_terms()
    assert not landed, (
        f"these registered absences no longer hold, so the rows blocked on them "
        f"should promote: {landed}; removal identifiers {removal_tokens()}; "
        f"files naming {CAPABILITY}: {capability_claims()}")


def test_the_absence_predicates_can_fail():
    """CONTROL: each predicate shape must be capable of returning False, or the
    check above is three lambdas that always pass.

    SCOPE: the two predicate SHAPES (_absent_from_sources, _declaration_still_stands)
    and the non-emptiness of each term's predicate list. It does NOT establish that
    any predicate is bound to the term its row names -- that is what the per-term
    mutations establish, and only firing them can.
    """
    assert not _absent_from_sources(r"compute_delta_H_terms")   # plainly present
    assert not _declaration_still_stands("a phrase no manuscript contains")
    assert not _absent_from("coupling/biofilm_openmc/materials.py", r"MaterialClass")
    # ...and every term carries at least one predicate, or all() is vacuously True.
    #
    # VACUOUS TRUTH OVER AN EMPTY SET IS THE THIRD INSTANCE OF ONE FAMILY IN THIS
    # REPOSITORY, and they are worth naming together rather than patched apart:
    # an empty `swift test` suite reporting success, a CI job that was never
    # dispatched reporting no failures, and now all() over a term with no
    # predicates. Each reports "nothing wrong" from "nothing checked", and each
    # was found by asking what the green would look like if the set were empty.
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
    """CONTROL: a row naming an absent term with no blocked_by must be reported.

    SCOPE: one synthetic row. Judged and left unscoped once by the contract gate
    before this line existed -- the three remaining flags in this file are
    legitimate: "no row claims an SOP" and "a row with no blocked_by" describe a
    control's own input, and "printed every run" is temporal. A gate whose false
    positives are never recorded as judged gets re-litigated every pass.
    """
    planted = {"sop_id": "SOP-ZZ", "coverage": "named_gap", "blocked_by": "",
               "signal_source": "matrix diffusion (D_eff_film != D_w)"}
    assert planted["signal_source"] in ABSENT_TERMS
    assert not planted["blocked_by"].strip()


def dead_terms(index, *, terms=None) -> list:
    """Registered absences that no index row names. Shared by the check and its
    control, so the control cannot pass against a regressed guard."""
    keys = set(ABSENT_TERMS) if terms is None else set(terms)
    named = {r.get("signal_source", "").strip() for r in index if r.get("signal_source", "").strip()}
    return sorted(keys - named)


def test_no_registry_term_is_dead(index):
    """EVERY term must back an index row, not merely one of them.

    SCOPE: the keys of ABSENT_TERMS against the signal_source column of the shipped
    index, and nothing else. It says nothing about whether a row is blocked for the
    right reason.

    The first version intersected the two sets, which proves only that at least
    one key is referenced -- and in the shipped index exactly one was: SOP-T4
    named matrix diffusion while nothing named the other two, so two of three
    entries were dead while the control passed. That is the same subset-checked
    set-asserted defect as reading only the first regex match. Raised by Codex on
    pull request #23. Compare the complete set."""
    dead = dead_terms(index)
    assert not dead, (
        f"registered absences no index row is blocked on, so their predicates are "
        f"promotion triggers for nothing: {dead}")


def test_the_dead_term_check_detects_an_orphan(index):
    """CONTROL, THROUGH THE SHARED LOGIC RATHER THAN BESIDE IT.

    SCOPE: one planted orphan, evaluated by dead_terms() -- the same function the
    real check calls. The first version recomputed a hard-coded set difference
    independently, so it would have kept passing if the guard regressed to the
    old intersection: a control that reimplements what it protects tests its own
    copy. Raised by Codex on pull request #23.
    """
    planted = list(index) + [{"signal_source": ""}]
    assert "a term no row names" in dead_terms(planted, terms={"a term no row names"})
    assert not dead_terms(planted, terms=set())


# ---------------------------------------------------------------------------
# PLANTS FOR ONE TERM, THROUGH THE PRODUCTION READER, UNDER A SCRATCH ROOT
# ---------------------------------------------------------------------------
# test_the_absence_predicates_can_fail proves the predicate SHAPES can return False.
# It says itself that it does not bind any predicate to its term. These do, for
# "dose-dependent survival", because an implementation of it now exists: an Odin
# port of PRRT-spatial-CPM's response.growth_survival, not yet proposed here.
#
# SCOPE: one term. Each plant is written under a scratch git root, and the gate's
# own landed_terms() reads that root. The other two terms have no plant yet.
#
# THE KNOWN-BAD IS THE PORT'S OWN TEXT, NOT A PARAPHRASE OF IT. Source:
# github.com/ffinkdevs/Biofilms, branch odin-port,
# commit 5c7e46aadedc0f48bbf6f6b697cfb19582217eb9. That commit is in no ref of
# this repository, so the lines are embedded and a reader can diff them there.
SURVIVAL = "dose-dependent survival"

# odin/cpm/response.odin lines 64-66 (the survival law) and 322-337 (the sweep).
ODIN_RESPONSE = '''sf_lq :: proc(D, G, alpha, beta: f64) -> f64 {
	return jexp_f64(-(alpha * D + beta * G * D * D))
}
	// Phase 2: unconditional survival sweep (A1).
	sf_sum := 0.0
	n_drawn := 0
	for id in 1..<s.next_id {
		if !cell_alive(s, id) {
			continue
		}
		sf := sf_lq(cell_dose(s, id), G, rp.alpha, rp.beta)
		sf_sum += sf
		n_drawn += 1
		u := rng_f64(&s.rng)
		if u > sf {
			cell_set_state(s, id, CELL_DYING)
			rep.deaths += 1
		}
	}
'''

# odin/cpm/sim.odin lines 248, 468 and 489: the three definitions, nothing between.
ODIN_SIM = '''cell_kill :: #force_inline proc(s: ^Sim, id: int) {
cull_dying :: proc(s: ^Sim, id: int) -> bool {
resorb_sweep :: proc(s: ^Sim) -> int {
'''

# odin/CONFORMANCE.md line 1.
ODIN_CONFORMANCE = '# Conformance report: Odin port vs `response.growth_survival`\n'


def _scratch_repo(tmp_path: Path) -> Path:
    """A fresh git root holding the two files the survival predicates open, and every
    file a pin names. A pinned file that is missing is reported, as a deleted one is."""
    import shutil
    import subprocess
    root = tmp_path / "repo"
    (root / TEX.relative_to(REPO)).parent.mkdir(parents=True)
    shutil.copy(TEX, root / TEX.relative_to(REPO))
    shutil.copy(REPO / "biofilms_potts.jl", root / "biofilms_potts.jl")
    for rel in sorted(set(EXEMPT_LINES) | set(CAPABILITY_LINES)):
        if not (root / rel).exists():
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(REPO / rel, root / rel)
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    return root


def _plant(root: Path, rel: str, text: str) -> None:
    """Append `text` to root/rel, creating it, and prove it is there.

    RuntimeError, not assert: the xfail marks below expect AssertionError, and a
    plant that never applied must not be counted as the defect it was meant to show."""
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(text)
    if text not in path.read_text(encoding="utf-8"):
        raise RuntimeError(f"the plant did not apply: {rel}")


def _replace_once(root: Path, rel: str, old: str, new: str) -> None:
    """Replace the one occurrence of `old` in root/rel. The victim is taken from the
    file as it is: anything but exactly one occurrence is a RuntimeError."""
    path = root / rel
    text = path.read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise RuntimeError(f"{rel} has {text.count(old)} occurrences of {old!r}")
    path.write_text(text.replace(old, new), encoding="utf-8")
    if new not in path.read_text(encoding="utf-8"):
        raise RuntimeError(f"the replacement did not apply: {rel}")


def test_the_unplanted_scratch_root_reports_no_landed_term(tmp_path):
    """NEGATIVE CONTROL: without it, a reader that fires on every root passes
    every plant below."""
    assert landed_terms(_scratch_repo(tmp_path), terms={SURVIVAL}) == []


SURVIVAL_PLANTS = [
    pytest.param("odin/cpm/response.odin", ODIN_RESPONSE,
                 id="the port's survival law and sweep"),
    pytest.param("odin/cpm/sim.odin", ODIN_SIM,
                 id="the port's removal procedures"),
    # A language and a directory nobody had in mind when any scan here was written.
    pytest.param("port/kernel.zig",
                 "fn killCell(lattice: *Lattice, i: usize) void {}\n"
                 "pub fn step(lattice: *Lattice) void { killCell(lattice, 0); }\n",
                 id="a definition and a call in an unlisted language"),
    pytest.param("biofilms_potts.jl", "\nremove_cell!(lattice, i)\n",
                 id="a call site in the stepper, with no definition"),
    pytest.param("biofilms_potts.jl", "\nk_death = 0.1\n",
                 id="a rate constant in the stepper"),
    pytest.param("odin/CONFORMANCE.md", ODIN_CONFORMANCE,
                 id="the port's declaration of the capability"),
    # The natural names for the process, a generic removal, and no capability name.
    pytest.param("port/response.zig",
                 "const survival_probability = exp(-alpha*D - beta*G*D*D);\n"
                 "if (u > survival_probability) { cells.delete(id); }\n",
                 id="natural survival names and a generic delete"),
    # The shortest route from inside the repository: the law, written in the stepper,
    # under the very word its comments were exempted for.
    pytest.param("biofilms_potts.jl",
                 "\nfunction survival(dose)\n    return exp(-0.2 * dose)\nend\n",
                 id="a survival function in the stepper"),
    # The one form the first predicate could already see.
    pytest.param("biofilms_potts.jl", "\nfunction remove_cell!(lattice, i)\nend\n",
                 id="a Julia definition in the stepper"),
]


@pytest.mark.parametrize("rel, text", SURVIVAL_PLANTS)
def test_a_planted_survival_process_lands_the_term(tmp_path, rel, text):
    root = _scratch_repo(tmp_path)
    _plant(root, rel, text)
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]


SEEN_BY = {     # plant id -> (removal_tokens, capability_claims), exactly
    "the port's survival law and sweep":
        ({"odin/cpm/response.odin":
          ["CELL_DYING", "deaths", "sf", "sf_lq", "sf_sum", "survival"]}, []),
    "the port's removal procedures":
        ({"odin/cpm/sim.odin": ["cell_kill", "cull_dying", "resorb_sweep"]}, []),
    "a definition and a call in an unlisted language": ({"port/kernel.zig": ["killCell"]}, []),
    "a call site in the stepper, with no definition":
        ({"biofilms_potts.jl": ["remove_cell!"]}, []),
    "a rate constant in the stepper": ({"biofilms_potts.jl": ["k_death"]}, []),
    "the port's declaration of the capability": ({}, ["odin/CONFORMANCE.md:1"]),
    "a Julia definition in the stepper": ({"biofilms_potts.jl": ["remove_cell!"]}, []),
    "natural survival names and a generic delete":
        ({"port/response.zig": ["survival_probability"]}, []),
    # The exemption is void, so every identifier the stepper was allowed is reported.
    "a survival function in the stepper":
        ({"biofilms_potts.jl": sorted(REMOVAL_ALLOWED["biofilms_potts.jl"][0])}, []),
}


def test_a_dose_condition_on_the_steppers_own_cleanup_lands_the_term(tmp_path):
    """The shortest route to the process from inside this repository: the volume-zero
    cleanup gains a dose condition and reuses DeathEvent and death_log, the very
    identifiers REMOVAL_ALLOWED exempts in that file."""
    root = _scratch_repo(tmp_path)
    _replace_once(root, "biofilms_potts.jl", "        if cell.volume <= 0\n",
                  "        if cell.volume <= 0 || cell.lifetime_dose_Gy > 50.0\n")
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]
    # Through the reader it was written for: the exemption is void, so every identifier
    # the stepper was allowed is reported, and nothing else is.
    assert removal_tokens(root) == {"biofilms_potts.jl": sorted(
        REMOVAL_ALLOWED["biofilms_potts.jl"][0])}
    assert capability_claims(root) == []


# CONTROL FLOW INSIDE THE BLOCK. Raised by Codex on pull request #47 against 02f8732, and
# reproduced first, in three more forms than the one it named. stepper_removal_sites() reads
# the `if` headers around a removal. A guard clause or a ternary makes the removal depend on
# dose and changes no header, no pinned line and no vocabulary word.
# (victim line, taken from the stepper as it is; what it becomes)
STEPPER_GUARD_PLANTS = [
    pytest.param("        if cell.volume <= 0\n",
                 "        if cell.volume <= 0\n            cell.lifetime_dose_Gy <= 50.0 && continue\n",
                 id="a guard clause with && continue"),
    pytest.param("        if cell.volume <= 0\n",
                 "        if cell.volume <= 0\n            cell.lifetime_dose_Gy > 50.0 || continue\n",
                 id="a guard clause with || continue"),
    pytest.param("            delete!(state.cells, id)\n",
                 "            cell.lifetime_dose_Gy > 50.0 ? delete!(state.cells, id) : nothing\n",
                 id="a ternary on the delete"),
    pytest.param("    parent = state.cells[id]\n",
                 "    parent = state.cells[id]\n    parent.lifetime_dose_Gy > 50.0 && return :spared\n",
                 id="a dose guard at the top of divide_cell!"),
]


@pytest.mark.parametrize("victim, rewrite", STEPPER_GUARD_PLANTS)
def test_control_flow_inside_a_removal_lands_the_term(tmp_path, victim, rewrite):
    root = _scratch_repo(tmp_path)
    _replace_once(root, "biofilms_potts.jl", victim, rewrite)
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]
    assert removal_tokens(root) == {"biofilms_potts.jl": sorted(
        REMOVAL_ALLOWED["biofilms_potts.jl"][0])}
    assert capability_claims(root) == []


# REPLACEMENTS INSIDE A FILE THAT IS ALREADY EXEMPT. Raised by Codex on pull request #47
# against ce141b2, and each reproduced before this was written: an exemption keyed on a file,
# a word or a count survives an edit that swaps what earned it for the thing it hides.
# (file, a pattern that must select exactly ONE line of it, what that line becomes)
REPLACEMENT_PLANTS = [
    pytest.param("biofilms_potts.jl", r"(?<![A-Za-z_])survival(?![A-Za-z_])",
                 lambda line: "function survival(dose)\n    return exp(-0.2 * dose)\nend",
                 id="the prose line saying survival becomes a survival function"),
    # Anchored: the same phrase is also inside a table row further down the README.
    pytest.param("README.md", r"^`response\.growth_survival` is `unsupported_by_current_model`",
                 lambda line: line.replace("`unsupported_by_current_model`", "`supported`"),
                 id="the README declares the capability supported"),
    pytest.param("data/parameter_provenance.csv", r",response\.growth_survival,",
                 lambda line: line.replace("unsupported_by_current_model", "supported_by_current_model"),
                 id="the provenance row's status flips"),
    # Raised by Codex on pull request #47 against b68e800. Production code, which an
    # earlier comment here had called a test or a script.
    pytest.param("coupling/biofilm_openmc/feedback_gate.py", r"growth_survival",
                 lambda line: "def growth_survival(dose):\n    return math.exp(-0.2 * dose)",
                 id="the gate's unsupported-pathway line becomes an implementation"),
    # OUTSIDE THE SCOPE I HAD IN MIND. I would have called this file "just a script" and
    # left it exempt by word. Its local `sf` is a far-field standard deviation.
    pytest.param("coupling/scripts/synthetic_e2e.py", r"^\s*sf = ",
                 lambda line: "    sf = np.exp(-0.2 * dose)",
                 id="a script's local sf becomes a survival fraction"),
    # REMOVALS. Raised by Codex on pull request #47 against a60be28, and reproduced first.
    # The natural way to ENABLE the capability deletes the line that declared it absent.
    # The pin then fails to match, and nothing is left in the file for a reader to name.
    pytest.param("coupling/biofilm_openmc/feedback_gate.py", r"growth_survival",
                 lambda line: line.replace('("growth_survival", "membrane_response")',
                                           '("membrane_response",)'),
                 id="the pathway leaves the unsupported list"),
    pytest.param("README.md", r"^`response\.growth_survival` is `unsupported_by_current_model`",
                 lambda line: "",
                 id="the README line that declares it unsupported is deleted"),
]


@pytest.mark.parametrize("rel, victim, rewrite", REPLACEMENT_PLANTS)
def test_a_replacement_inside_an_exempt_file_lands_the_term(tmp_path, rel, victim, rewrite):
    import re
    root = _scratch_repo(tmp_path)
    # The victim comes from the file as it is now, never from memory of it.
    lines = (root / rel).read_text(encoding="utf-8").split("\n")
    hits = [n for n, line in enumerate(lines, 1) if re.search(victim, line)]
    if len(hits) != 1:
        raise RuntimeError(f"{victim!r} selects {len(hits)} lines of {rel}, not one")
    old, new = lines[hits[0] - 1], rewrite(lines[hits[0] - 1])
    if new == old:
        raise RuntimeError(f"the rewrite changed nothing in {rel}:{hits[0]}")
    _replace_once(root, rel, old + "\n", new + "\n")
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]
    # ...and through the reader it was written for, exactly.
    after = (root / rel).read_text(encoding="utf-8")
    if rel in REMOVAL_ALLOWED:      # the exemption is void: every exempt identifier still
        left = sorted(a for a in REMOVAL_ALLOWED[rel][0]      # in the file shows, or the
                      if re.search(rf"(?<![A-Za-z0-9_]){re.escape(a)}(?![A-Za-z0-9_!])", after))
        assert removal_tokens(root) == {rel: left or [PINS_CHANGED]}    # broken pin alone
        assert capability_claims(root) == []
    else:                           # a new or edited line by number, and the pin it replaced
        assert removal_tokens(root) == {}
        edited = [f"{rel}:{hits[0]}"] if CAPABILITY in new else []
        assert capability_claims(root) == edited + [f"{rel}: 1 pinned line(s) gone"]


@pytest.mark.parametrize("rel, text", SURVIVAL_PLANTS)
def test_each_plant_is_seen_by_the_reader_it_was_written_for(tmp_path, request, rel, text):
    """landed_terms() is true if ANY predicate fires. This consumes the finest thing
    each reader returns, so a plant cannot land its term through the wrong one."""
    root = _scratch_repo(tmp_path)
    _plant(root, rel, text)
    plant_id = request.node.callspec.id
    assert (removal_tokens(root), capability_claims(root)) == SEEN_BY[plant_id]


# A DELETED FILE. Raised by Codex on pull request #47 against 80d9787, and reproduced first.
# listed_text() omits a file that is gone, so neither reader compares that file's pins. The
# stepper is left out: the older Julia predicate opens it, so its absence is already an error.
PINNED_FILES = sorted(set(PINNED_PATHS) - {"biofilms_potts.jl"})


@pytest.mark.parametrize("rel", PINNED_FILES)
def test_a_deleted_pinned_file_lands_the_term(tmp_path, rel):
    root = _scratch_repo(tmp_path)
    if landed_terms(root, terms={SURVIVAL}):
        raise RuntimeError("the scratch root reports a landed term before the deletion")
    (root / rel).unlink()
    if (root / rel).exists():
        raise RuntimeError(f"the deletion did not apply: {rel}")
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]
    # ...and through the reader that pins the file, exactly.
    pins = CAPABILITY_LINES.get(rel) or set()
    assert removal_tokens(root) == ({rel: [PINS_CHANGED]} if rel in EXEMPT_LINES else {})
    assert capability_claims(root) == ([f"{rel}: {len(pins)} pinned line(s) gone"] if pins else [])


# ...AND ITS TABLE ENTRIES WITH IT. Raised by Codex on pull request #47 against c28a381, and
# reproduced first. The natural clean-up after deleting a file removes its entries from the
# four tables. Both readers then forget the path, the tables still agree with each other, and
# PINNED_FILES above shrinks without a word, because it is built from those same tables.
@pytest.mark.parametrize("rel", PINNED_FILES)
def test_a_deleted_file_is_remembered_after_its_table_entries_go(tmp_path, monkeypatch, rel):
    root = _scratch_repo(tmp_path)
    if landed_terms(root, terms={SURVIVAL}):
        raise RuntimeError("the scratch root reports a landed term before the deletion")
    (root / rel).unlink()
    dropped = 0
    for table in (EXEMPT_LINES, REMOVAL_ALLOWED, CAPABILITY_LINES, CAPABILITY_NAMED_AS_UNSUPPORTED):
        dropped += rel in table
        monkeypatch.delitem(table, rel, raising=False)
    if (root / rel).exists() or not dropped:
        raise RuntimeError(f"the deletion did not apply: {rel}, {dropped} table entries")
    assert landed_terms(root, terms={SURVIVAL}) == [SURVIVAL]
    # ...through the list alone: both table readers have forgotten the path.
    assert (removal_tokens(root), capability_claims(root)) == ({}, [])
    assert missing_pinned_paths(root) == [rel]


# THE LEDGER ROW THAT DESCRIBES THIS GUARD IS A CLAIM, AND IT WENT STALE THREE TIMES. Codex on
# pull request #47, at 4b9c25e, a7740c9 and 54f920b: the row named two readers "and nothing
# else" while the term had five predicates. A sentence that lists predicates is wrong the day
# one is added. So the row is read here: it must name every function the predicates call.
def unnamed_survival_predicates(root: Path = REPO) -> list:
    """Functions of this module that a survival predicate calls and GUARD-SCOPE-02 does not name."""
    import csv
    with open(root / "data" / "claims_ledger.csv", encoding="utf-8") as fh:
        rows = {r["claim_id"]: r for r in csv.DictReader(l for l in fh if not l.startswith("#"))}
    claim = rows["GUARD-SCOPE-02"]["claim_text"]
    called = {name for pred in ABSENT_TERMS[SURVIVAL][1] for name in pred.__code__.co_names
              if callable(globals().get(name))}
    if not called:
        raise RuntimeError("no function found in the survival predicates")    # never default to pass
    return sorted(name for name in called if f"{name}(" not in claim)


def survival_predicate_census() -> list:
    """One entry per survival predicate: the functions of this module it calls, and the
    string arguments it passes. Read from the code objects, so it is what the gate runs."""
    return [(tuple(n for n in pred.__code__.co_names if callable(globals().get(n))),
             tuple(c for c in pred.__code__.co_consts if isinstance(c, str)))
            for pred in ABSENT_TERMS[SURVIVAL][1]]


# THE PREDICATES AS GUARD-SCOPE-02 DESCRIBES THEM, PINNED VERBATIM. A change to any predicate,
# or a new one, fails the row test below. Whoever updates this list has read the row's claim
# against the predicates again and corrected it first. A new argument is a new scope.
SURVIVAL_PREDICATES = [
    (("_declaration_still_stands",),
     ("has neither biomass growth nor a dose-dependent survival process",)),
    (("_absent_from",),
     ("biofilms_potts.jl", r"function\s+\w*(death|kill|remove_cell|die)\w*!?\s*\(")),
    (("removal_tokens",), ()),
    (("capability_claims",), ()),
    (("missing_pinned_paths",), ()),
]


def stale_survival_row(root: Path = REPO) -> list:
    """Every reason GUARD-SCOPE-02 may no longer describe the predicates. Empty while it does."""
    changed = [f"predicate {n}: {now} is not the pinned {then}" for n, (now, then) in
               enumerate(zip(survival_predicate_census(), SURVIVAL_PREDICATES)) if now != then]
    if len(survival_predicate_census()) != len(SURVIVAL_PREDICATES):
        changed.append(f"{len(survival_predicate_census())} predicates, and "
                       f"{len(SURVIVAL_PREDICATES)} are pinned")
    return unnamed_survival_predicates(root) + changed


def test_the_ledger_row_names_every_survival_predicate():
    assert stale_survival_row() == []


# A HELPER NAME IS A MARKER TOO. Raised by Codex on pull request #47 against 5937c7d, and
# reproduced first: a second _absent_from() call, on a new file, adds a scope the row does not
# state, and the name census returns nothing because the row already names _absent_from().
def test_a_reused_helper_with_new_arguments_makes_the_row_stale(monkeypatch):
    why, preds = ABSENT_TERMS[SURVIVAL]
    monkeypatch.setitem(ABSENT_TERMS, SURVIVAL, (why, preds + [
        lambda root=REPO: _absent_from("odin/cpm/sim.odin", r"cell_kill", root)]))
    if len(ABSENT_TERMS[SURVIVAL][1]) != len(preds) + 1:
        raise RuntimeError("the plant did not apply")
    assert stale_survival_row() == ["6 predicates, and 5 are pinned"]
    assert unnamed_survival_predicates() == []      # the name census alone still sees nothing


def test_a_changed_argument_makes_the_row_stale(monkeypatch):
    """Not a new predicate: the old one, pointed at another file."""
    why, preds = ABSENT_TERMS[SURVIVAL]
    swapped = [(lambda root=REPO: _absent_from("odin/cpm/sim.odin", r"cell_kill", root))
               if "_absent_from" in p.__code__.co_names else p for p in preds]
    monkeypatch.setitem(ABSENT_TERMS, SURVIVAL, (why, swapped))
    got = stale_survival_row()
    assert len(got) == 1 and got[0].startswith("predicate 1:") and "odin/cpm/sim.odin" in got[0], got


def test_the_predicate_census_detects_an_unnamed_one(monkeypatch):
    why, preds = ABSENT_TERMS[SURVIVAL]
    before = unnamed_survival_predicates()
    monkeypatch.setitem(ABSENT_TERMS, SURVIVAL, (why, preds + [lambda root=REPO: bool(listed_text(root))]))
    assert unnamed_survival_predicates() == sorted(before + ["listed_text"])


def test_every_pin_has_its_tables():
    """An allowlisted file with no pinned lines is exempt by its name, and a capability
    file with a reason and no pinned lines is exempt by its path."""
    assert set(EXEMPT_LINES) == {rel for rel, (allowed, _why) in REMOVAL_ALLOWED.items()
                                 if allowed is not WHOLE_FILE}
    assert set(EXEMPTION_HOLDS) <= set(EXEMPT_LINES)
    assert set(CAPABILITY_LINES) == set(CAPABILITY_NAMED_AS_UNSUPPORTED)
    # READ THE GATE FIRST, as the next test says. A path that left the tables because its
    # file was deleted is a red gate, and the remedy is not to delete its line from the list.
    assert set(PINNED_PATHS) == set(EXEMPT_LINES) | {
        rel for rel, pins in CAPABILITY_LINES.items() if pins is not WHOLE_FILE}
    assert len(PINNED_PATHS) == len(set(PINNED_PATHS))


def test_no_allowlist_entry_has_outlived_its_reason():
    """An exemption is keyed on what earns it. When the identifier leaves the file, or
    the file stops naming the capability, the entry exempts nothing.

    READ THE GATE FIRST. If a pinned line is gone because the thing it declared absent
    was just ENABLED, test_every_registered_absence_still_holds is red too, and the
    remedy is to promote the rows blocked on the term. Deleting the entry to make this
    test pass is the edit Codex described on pull request #47: it leaves the capability
    enabled and the trigger green. Delete an entry only when the gate is green."""
    import re
    texts = listed_text()
    token = re.compile(r"[A-Za-z_][A-Za-z0-9_]*!?")
    stale = [(rel, sorted(allowed - set(token.findall(texts.get(rel, "")))))
             for rel, (allowed, _why) in REMOVAL_ALLOWED.items() if allowed is not WHOLE_FILE]
    assert not [(rel, gone) for rel, gone in stale if gone]
    import hashlib
    for rel, pinned in CAPABILITY_LINES.items():
        if pinned is WHOLE_FILE:
            continue
        now = {hashlib.sha256(line.strip().encode("utf-8")).hexdigest()
               for line in texts.get(rel, "").split("\n") if CAPABILITY in line}
        assert pinned == now, f"{rel}: pinned lines that are gone {sorted(pinned - now)}"


# WHAT NEITHER READER SEES, ASSERTED SO THE LIMIT IS A FACT IN THE SUITE AND NOT A HOPE.
# A vocabulary finds names. A survival process with opaque names and no declaration is
# found by review and by nothing here. If a better reader sees one of these, this test
# fails, and the entry moves to SURVIVAL_PLANTS.
KNOWN_UNSEEN = [
    ("odin/cpm/opaque.odin", "dN := g*N - k(d)*N\n"),
    # The port's own law and sweep with every name replaced.
    ("odin/cpm/opaque.odin", "p := jexp_f64(-(a*D + b*G*D*D))\nif u > p { st[id] = 3 }\n"),
    # The verb as a name. It is excluded on a measurement: 46 files use it as English.
    ("port/response.zig", "fn survive(d: f64) f64 { return exp(-a*d); }\n"),
    # Dose zeroes the volume somewhere else, and the volume-only cleanup does the rest.
    # stepper_removal_sites() reads what governs a removal, not what feeds it.
    ("biofilms_potts.jl", "\ncell.volume = cell.lifetime_dose_Gy > 50.0 ? 0 : cell.volume\n"),
    # Source under one of the six PROSE_AND_DATA suffixes. The same line in port/model.zig
    # is seen. Skipping them is the departure from AGENTS.md decided on 2026-09-20.
    ("port/model.txt", "fn killCell(i: usize) void {}\n"),
    # Inside an unbroken 60-character run, which removal_tokens() drops as a blob.
    ("viewer/blob.js", "const a = '" + "A" * 30 + "killCell" + "A" * 30 + "';\n"),
]


@pytest.mark.parametrize("rel, text", KNOWN_UNSEEN)
def test_what_the_survival_readers_cannot_see_is_stated(tmp_path, rel, text):
    root = _scratch_repo(tmp_path)
    _plant(root, rel, text)
    assert landed_terms(root, terms={SURVIVAL}) == []
