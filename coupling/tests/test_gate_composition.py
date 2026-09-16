"""The composition seam nothing else in this repo tests: dose output feeding
an actual gate decision.

`coupling/scripts/synthetic_e2e.py` already composes the physics chain
end-to-end with real OpenMC transport (snapshot -> model -> transport ->
exact-CSG mass -> dose -> mesh -> lineage attribution -> results), and stops.
`synthetic_gate.decide()` and `feedback_gate`'s gates are tested elsewhere,
exclusively against hand-fabricated numpy arrays (`test_synthetic_gate_fixtures.py`,
`test_feedback_gate.py`, `test_uq_estimator.py`) -- never against anything
`dose.py`/`feedback_uq.py` actually derived from a tally. The one place in the
repo that wires real transport -> dose -> a gate decision is
`openmc_nested_pilot.py`'s `_report()` (the `decide()` call), and it has no
test coverage for that composition at all.

This file closes that seam using a REAL pinned transport result, not a mock:
`coupling/tests/fixtures/golden_tally_water_phantom.json` is 12 actual OpenMC
runs (2 outer draws x 3 replicates x {baseline, feedback}, density x1.35 --
the same DENSITY_SCALE lever `openmc_nested_pilot.py` already uses), generated
by `coupling/scripts/regenerate_golden_tally.py` under the verified
openmc-biofilms env. No OpenMC is needed HERE: this test only replays the
pinned heating tally through the real, live repo code -- specific_energy_per_source
-> debiased_squared_effect -> decide() -- the same compare-only-in-CI split
`tests/contract_csv.jl` uses for the serial fixture.
"""

from __future__ import annotations

import ast
import pytest
import sys
import json
from pathlib import Path

import numpy as np

from biofilm_openmc.dose import specific_energy_per_source
from biofilm_openmc.feedback_uq import debiased_squared_effect
from biofilm_openmc.synthetic_gate import (PASS_SYNTHETIC_GATE, ThresholdPolicy,
                                           VarianceBudget, decide)

# PRODUCTION'S POLICY, NOT A BARE ONE. The draws below are `e_squared`, which
# feedback_uq declares as `debiased_relative_l2_squared`; a default
# ThresholdPolicy() is denominated in `relative_l2` and carries 0.10 / 0.02.
# Comparing one against the other is the error decide() exists to refuse, and it
# refuses only when the caller declares the metric. The first version of this
# file passed neither the production policy nor the metric, and pinned
# EFFECT_BELOW_THRESHOLD -- a verdict production's own policy does not produce.
# The paths mirror what openmc_nested_pilot.py does for itself, so the bare
# no-OpenMC tier can import it without calibration installed.
_COUPLING = Path(__file__).resolve().parents[1]
_ROOT = _COUPLING.parent
for _p in (_COUPLING / "scripts", _ROOT / "calibration", _ROOT / "contract"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from openmc_nested_pilot import POLICY  # noqa: E402

_FIXTURE = (Path(__file__).parent / "fixtures"
           / "golden_tally_water_phantom.json")

with open(_FIXTURE) as f:
    _DATA = json.load(f)

_MASS_B = np.array(_DATA["mass_kg_baseline"])
_MASS_F = np.array(_DATA["mass_kg_feedback"])
_N_OUTER = _DATA["n_outer"]
_N_REP = _DATA["n_replicates"]


def _per_source_fields(label: str, mass: np.ndarray, heating_scale: float = 1.0):
    """specific_energy_per_source for every replicate run of one condition,
    from the pinned raw heating tally -- REAL repo code, not a re-derivation."""
    out = []
    for run in _DATA["runs"][label]:
        heating = np.array(run["heating_mean_eV_per_src"]) * heating_scale
        out.append(specific_energy_per_source(heating, mass))
    return out


def _effect_draws(feedback_heating_scale: float = 1.0):
    """One e_squared per outer draw, with the metric the producer declares.

    THE DECLARATION TRAVELS WITH THE NUMBERS. Returning the draws alone is what
    let this file hand `e_squared` to a policy denominated in E; the consumer
    has to read what the producer said it computed.
    """
    baseline_fields = _per_source_fields("baseline", _MASS_B)
    feedback_fields = _per_source_fields("feedback", _MASS_F,
                                         feedback_heating_scale)
    draws, metrics = [], set()
    for o in range(_N_OUTER):
        lo, hi = o * _N_REP, (o + 1) * _N_REP
        effect = debiased_squared_effect(baseline_fields[lo:hi],
                                         feedback_fields[lo:hi],
                                         _MASS_B.ravel())
        draws.append(effect.e_squared)
        metrics.add(effect.metric_id)
    assert len(metrics) == 1, f"the producer changed metric mid-run: {metrics}"
    return np.array(draws), metrics.pop()


def test_pinned_transport_composes_into_a_verdict():
    """The real seam: pinned heating -> dose -> debiased effect -> decide().

    The measured effect is small (~0.0011 in E^2): a x1.35 density increase
    raises heating by ~26% but raises mass by 35%, so specific energy per kg
    barely moves -- a real, unforced result, not tuned to land anywhere.
    Recorded here as the pinned expectation, the same way `serial_seed42.csv`
    pins what `validate_serial.jl` actually produces.

    Under production's POLICY that is EFFECT_DETECTED_BUT_NOT_PRACTICALLY_
    IMPORTANT: the draws resolve above the practical floor (0.02^2 = 0.0004)
    and sit entirely below the threshold (0.10^2 = 0.01). The verdict pinned
    here until 2026-09-16 was EFFECT_BELOW_THRESHOLD, which this same fixture
    produces only against a policy denominated in E rather than E^2.
    """
    draws, metric = _effect_draws()
    assert draws.shape == (_N_OUTER,)
    assert np.all(np.isfinite(draws))
    assert np.all(draws < 0.01), (
        f"effect_draws {draws} drifted far from the pinned expectation; "
        "regenerate the fixture only if OpenMC or the nuclear data changed")

    verdict = decide(draws, VarianceBudget(transport=1e-4), POLICY,
                     metric_id=metric)
    assert verdict.verdict == "EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT", \
        verdict.reason


def test_scaling_the_dose_changes_the_verdict():
    """THE SECOND NEGATIVE CONTROL. Proving decide() gets CALLED with real
    data doesn't prove the VALUE it's called with matters -- a decide() stub
    returning a constant would still be reachable and still get called, and
    a call-site check alone would not catch it. So: scale the pinned
    feedback heating by a real, meaningful factor and require the verdict to
    actually change. If it didn't, the composition would be running but
    inert -- the gate-that-never-fires failure one level up from the seam
    this file exists to close.
    """
    base_draws, metric = _effect_draws()
    scaled_draws, scaled_metric = _effect_draws(feedback_heating_scale=5.0)
    baseline_verdict = decide(base_draws, VarianceBudget(transport=1e-4),
                              POLICY, metric_id=metric)
    scaled_verdict = decide(scaled_draws, VarianceBudget(transport=1e-4),
                            POLICY, metric_id=scaled_metric)

    assert scaled_verdict.verdict != baseline_verdict.verdict, (
        "scaling the feedback dose 5x did not change the verdict -- the "
        "pipeline runs but the gate's output does not depend on its input")
    assert scaled_verdict.verdict == PASS_SYNTHETIC_GATE


def test_a_policy_in_the_wrong_denomination_is_refused():
    """THE CONTROL THE PINNED TEST WAS MISSING, and the reason it was wrong.

    decide() refuses a threshold that is not denominated in the draws' own
    metric -- but only when the caller declares that metric. Passing nothing is
    silent, and silence is what produced the old pinned verdict. Both halves are
    asserted here: declared, the mismatch is NOT_EVALUATED; undeclared, the same
    numbers and the same policy produce a verdict that reads as a result.
    """
    draws, metric = _effect_draws()
    assert POLICY.metric_id == metric, (
        "production's policy is no longer denominated in the metric the "
        "producer declares; the pinned verdict below is then meaningless")

    refused = decide(draws, VarianceBudget(transport=1e-4), ThresholdPolicy(),
                     metric_id=metric)
    assert refused.verdict == "NOT_EVALUATED", refused.verdict
    assert "not transferable" in refused.reason, refused.reason

    silent = decide(draws, VarianceBudget(transport=1e-4), ThresholdPolicy())
    assert silent.verdict == "EFFECT_BELOW_THRESHOLD", (
        "the undeclared call no longer produces the old verdict, so this "
        "control no longer reproduces the defect it exists for")


_REPO = Path(__file__).resolve().parents[2]
_REGEN = _REPO / "coupling" / "scripts" / "regenerate_golden_tally.py"
_WORKFLOW = (_REPO / ".github" / "workflows"
             / "golden-tally-verification.yml")


# Every directory a first-party absolute import can resolve against. The pilot
# puts `calibration` on sys.path and imports `biofilm_calibration` from it, so
# the calibration root is a base like the others; leaving it out hid
# joint_uncertainty.py from the closure.
_FIRST_PARTY_BASES = ("coupling", "coupling/scripts", "contract", "calibration")


def _first_party_imports(path: Path) -> set[Path]:
    """The first-party modules one file imports, as repo-relative paths.
    Walks the real AST, function-local imports included. A dotted name
    resolves to `name.py` or to a package's `__init__.py`."""
    tree = ast.parse(path.read_text())
    out: set[Path] = set()

    def add(rel: Path, roots) -> None:
        for root in roots:
            for candidate in (root / (rel.as_posix() + ".py"),
                              root / rel / "__init__.py"):
                if candidate.exists():
                    out.add(candidate.relative_to(_REPO))
                    # IMPLICIT PACKAGE INITIALISERS EXECUTE TOO. Importing
                    # `biofilm_openmc.config` runs `biofilm_openmc/__init__.py`
                    # first, whether or not anything names it; every
                    # `__init__.py` between the root and the module is part of
                    # the generating run and belongs in the closure.
                    for parent in candidate.relative_to(root).parents:
                        init = root / parent / "__init__.py"
                        if init.exists():
                            out.add(init.relative_to(_REPO))

    absolute_roots = [_REPO / b for b in _FIRST_PARTY_BASES]
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.level:
                # RELATIVE IMPORTS RESOLVE AGAINST THE IMPORTING PACKAGE.
                # `from .snapshot import ...` has `module == "snapshot"` and
                # `level == 1`; the first version skipped every node whose
                # module was None and resolved the rest as absolute, so a
                # relative dependency was invisible to the closure and could
                # be absent from both workflow path lists while this passed.
                pkg = path.parent
                for _ in range(node.level - 1):
                    pkg = pkg.parent
                base, roots = Path(*node.module.split(".")) if node.module else Path(), [pkg]
            else:
                base, roots = Path(*node.module.split(".")), absolute_roots
            if node.module:
                add(base, roots)
            # AN ALIAS MAY BE A SUBMODULE. `from biofilm_openmc import x`
            # loads `biofilm_openmc.x` when x is a module; the walker used to
            # record the package initialiser and stop, so a producer pulled
            # in by name ran without being in either workflow filter. `add`
            # keeps a name only when a file exists for it, so an attribute
            # resolves to nothing and is ignored.
            for alias in node.names:
                add(base / alias.name, roots)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                add(Path(*alias.name.split(".")), absolute_roots)
    return out


def _modules_that_produce_the_fixture() -> set[Path]:
    """Every first-party module the regeneration script reaches, TRANSITIVELY.

    The first version walked the script's own imports and nothing further, so
    `physical_contract` -- imported by `biofilm_openmc.config` and `.model`,
    which the script imports inside `_run_one` -- was never in the inventory,
    and the trigger test could not know the workflow omitted it. A change to
    the shared vocabularies or validation there alters the generating run
    while neither verification job fires. The closure over first-party
    imports is the inventory; a base directory is added to
    `_FIRST_PARTY_BASES`, not a module name."""
    # THE ENTRY POINT IS ITSELF A PRODUCER. Seeding only the frontier returned the
    # script's imports and not the script, so an edit to the generator's own body
    # could omit it from both filters while this closure reported nothing missing.
    seen: set[Path] = {_REGEN.relative_to(_REPO)}
    frontier = [_REGEN.relative_to(_REPO)]
    while frontier:
        here = frontier.pop()
        for dep in _first_party_imports(_REPO / here):
            if dep not in seen:
                seen.add(dep)
                frontier.append(dep)
    return seen


# THE BUILD SPECIFICATIONS ARE INPUTS TOO. The verification job installs `contract` and
# `coupling[dev]` editably, so a dependency pin or a package-data rule in either manifest
# changes the generating run without touching a single producer module.
_BUILD_SPECS = ("coupling/pyproject.toml", "contract/pyproject.toml")


def _fixture_inputs() -> set[str]:
    """Every repo-relative path whose change alters the generating run: the transitive
    producer closure, the regeneration entry point included, plus the build specs."""
    return {p.as_posix() for p in _modules_that_produce_the_fixture()} | set(_BUILD_SPECS)


def _paths_missing_from_triggers(required: set[str], triggers: dict) -> dict[str, list[str]]:
    """Per path-filtered event, the required inputs its `paths:` list does not name."""
    return {event: sorted(p for p in required
                          if p not in set(triggers.get(event, {}).get("paths", [])))
            for event in ("push", "pull_request")}


def _workflow_triggers() -> dict:
    """The workflow's `on:` block, PARSED. Not grepped: a substring search over
    the file text finds a path listed under `push:` alone and reports the
    `pull_request:` filter covered too, which is exactly the hole this file
    had. YAML also parses the bare key `on` as the boolean True."""
    import yaml
    doc = yaml.safe_load(_WORKFLOW.read_text())
    return doc[True] if True in doc else doc["on"]


def test_every_fixture_producing_module_triggers_verification():
    """THE CONTROL FOR THE WORKFLOW'S `paths:` FILTERS, which are hand-
    maintained lists and therefore things that drift.

    The filter's first version covered only the OpenMC pin and the
    nuclear-data pin, on the stated premise that nothing else can change what
    real transport produces. False: `_run_one` imports config, model, mesh,
    materials and dose to build and tally the phantom, and edits to any of
    them move the real output while `test_pinned_transport_composes_into_a_
    verdict` above keeps replaying the committed JSON and passing. A stale
    pin that nothing can invalidate is a check that cannot fail.

    So derive the list from the script's actual imports and require EVERY
    path-filtered trigger to name each one. Add an import to `_run_one`
    without touching the filters and this fails, naming the file to add.
    """
    triggers = _workflow_triggers()
    producers = _modules_that_produce_the_fixture()

    # The walk itself must find something, or this test passes vacuously the
    # moment the AST parse or the path resolution breaks.
    assert len(producers) >= 6, (
        f"only found {sorted(map(str, producers))}; the import walk is not "
        "resolving the modules it is supposed to check")

    for event, missing in _paths_missing_from_triggers(_fixture_inputs(), triggers).items():
        assert not missing, (
            f"the generating run reads {missing}, so a change to any of them "
            f"changes the real tally -- but {_WORKFLOW.name}'s `{event}:` "
            "paths filter does not list them, so verification would not run "
            "and the committed fixture would go stale while CI stayed green.")


def test_the_entry_point_is_in_its_own_closure():
    """The closure used to seed only its frontier with the generator, so the returned
    set held everything the generator imports and not the generator. Seed the frontier
    alone again and this fails."""
    producers = {p.as_posix() for p in _modules_that_produce_the_fixture()}
    assert _REGEN.relative_to(_REPO).as_posix() in producers, sorted(producers)


@pytest.mark.parametrize("event", ["push", "pull_request"])
@pytest.mark.parametrize("path", ["coupling/scripts/regenerate_golden_tally.py", *_BUILD_SPECS])
def test_removing_a_required_input_from_one_trigger_is_caught(event, path):
    """KNOWN-BAD: the real parsed workflow with ONE required path dropped from ONE event.
    The guard must name exactly that event and that path -- a closure without its entry
    point, or an input list without the manifests, could not. The two lists are copied
    separately on purpose: in the parsed YAML they are one anchored object, and a
    deepcopy keeps them shared, so removing from one would remove from both."""
    real = _workflow_triggers()
    triggers = {e: {"paths": list(real[e]["paths"])} for e in ("push", "pull_request")}
    assert path in triggers[event]["paths"], "the real workflow must list it: that is the premise"
    triggers[event]["paths"].remove(path)
    missing = _paths_missing_from_triggers(_fixture_inputs(), triggers)
    assert missing[event] == [path], missing
    other = "pull_request" if event == "push" else "push"
    assert missing[other] == [], missing


def test_relative_imports_are_part_of_the_closure():
    """`biofilm_openmc.model` reaches `snapshot` as `from .snapshot import`,
    which has no absolute module name. A walk that resolved only absolute
    imports could not see it, so a relative-only dependency that changes the
    fixture would be missing from both path lists while the trigger test
    stayed green. Drop the `node.level` branch and this fails."""
    deps = {p.as_posix() for p in _first_party_imports(
        _REPO / "coupling" / "biofilm_openmc" / "model.py")}
    assert "coupling/biofilm_openmc/snapshot.py" in deps, sorted(deps)


def test_from_import_aliases_that_name_submodules_are_resolved(tmp_path, monkeypatch):
    """`from biofilm_openmc import producer` LOADS `biofilm_openmc.producer`
    when that name is a submodule, and so does `from .sub import leaf`. The
    walker recorded the package initialiser and stopped, so a producer pulled
    in by name entered the generating run without entering either workflow
    filter. Each alias is now tried as a child module and kept when the file
    exists; a name that resolves to no file is an attribute and is ignored.

    KNOWN-BAD, on a throwaway tree so nothing is written into the repository:
    a package with a child module and a nested package with a leaf, imported
    only through from-import aliases. Disable the alias branch and this fails.
    """
    root = tmp_path / "coupling"
    (root / "pkg" / "sub").mkdir(parents=True)
    (root / "pkg" / "__init__.py").write_text("VALUE = 1\n")
    (root / "pkg" / "child.py").write_text("")
    (root / "pkg" / "sub" / "__init__.py").write_text("")
    (root / "pkg" / "sub" / "leaf.py").write_text("")
    (root / "pkg" / "user.py").write_text("from .sub import leaf\n")
    (root / "main.py").write_text("from pkg import child, VALUE\n")
    monkeypatch.setattr(sys.modules[__name__], "_REPO", tmp_path)
    monkeypatch.setattr(sys.modules[__name__], "_FIRST_PARTY_BASES", ("coupling",))

    absolute = {p.as_posix() for p in _first_party_imports(root / "main.py")}
    assert "coupling/pkg/child.py" in absolute, sorted(absolute)
    assert "coupling/pkg/__init__.py" in absolute, sorted(absolute)
    assert not [p for p in absolute if "VALUE" in p], sorted(absolute)

    relative = {p.as_posix() for p in _first_party_imports(root / "pkg" / "user.py")}
    assert "coupling/pkg/sub/leaf.py" in relative, sorted(relative)
    assert "coupling/pkg/sub/__init__.py" in relative, sorted(relative)


def test_the_calibration_root_is_part_of_the_closure():
    """`openmc_nested_pilot.py` puts `calibration` on sys.path and imports
    `biofilm_calibration.joint_uncertainty`, the correlated-draw sampler the
    pilot's outer draws come from. With the calibration root missing from
    `_FIRST_PARTY_BASES` the module resolved nowhere and was silently absent
    from the closure and the workflow. Drop "calibration" from the bases and
    this fails."""
    producers = {p.as_posix() for p in _modules_that_produce_the_fixture()}
    assert "calibration/biofilm_calibration/joint_uncertainty.py" in producers, \
        sorted(producers)


def test_implicit_package_initialisers_are_part_of_the_closure():
    """Importing any `biofilm_openmc.*` module executes
    `coupling/biofilm_openmc/__init__.py`; nothing names it, so a walk over
    named imports never returned it. Disable the initialiser rule in `add`
    and this fails."""
    producers = {p.as_posix() for p in _modules_that_produce_the_fixture()}
    for init in ("coupling/biofilm_openmc/__init__.py",
                 "calibration/biofilm_calibration/__init__.py"):
        assert init in producers, sorted(producers)


def test_the_shared_contract_package_is_a_fixture_producer():
    """THE INVENTORY MUST SEE THROUGH ONE IMPORT. `regenerate_golden_tally.py`
    never names `physical_contract`; `biofilm_openmc.config` and `.model` do,
    and both are imported inside `_run_one`. A walk that stopped at the
    script's own imports listed neither the package nor, therefore, the
    workflow's omission of it. Remove "contract" from `_FIRST_PARTY_BASES`
    and this fails; remove the path from the workflow and the trigger test
    above fails."""
    producers = {p.as_posix() for p in _modules_that_produce_the_fixture()}
    assert "contract/physical_contract/__init__.py" in producers, sorted(producers)
    assert "coupling/biofilm_openmc/config.py" in producers, sorted(producers)


def test_a_pull_request_can_reach_this_workflow_at_all():
    """THE TRIGGER-LEVEL CONTROL. The test above proves the paths are right;
    it cannot notice that an entire EVENT is absent, and that was the real
    defect: `push` is branch-filtered to master/feat/ci/research, so a PR from
    `fix/**` -- or from any fork, which matches no branch filter -- never
    reached this job, while coupling-tests.yml passed by replaying the
    committed tally. A fixture-producing change could merge with the pin never
    revalidated.

    Asserting the paths without asserting the event is the hole that let it
    ship, so assert the event, and assert the two lists are the SAME list --
    not merely equal today. YAML anchors give one source; two copies drift,
    and the copy that drifts is the pre-merge one nobody watches.
    """
    triggers = _workflow_triggers()

    assert "pull_request" in triggers, (
        "no pull_request trigger: this job cannot run before a merge from any "
        "branch outside the push filter, which includes every fork")
    assert triggers["pull_request"].get("paths"), (
        "the pull_request trigger has no paths filter, so it either never "
        "fires or fires on everything -- both are wrong for a 45-minute job")
    assert triggers["pull_request"]["paths"] == triggers["push"]["paths"]
    # Same list object, from the YAML anchor -- not two lists that happen to
    # match on the day someone last synchronised them.
    assert triggers["pull_request"]["paths"] is triggers["push"]["paths"], (
        "the two path lists are equal but separate, so they can drift; use "
        "the `&fixture_inputs` / `*fixture_inputs` anchor")


_ENV_FILE = _REPO / "environment.yml"
_STACK_DOC = _REPO / "docs" / "openmc_stack.md"


def test_the_transport_environment_has_exactly_one_spec():
    """THE FOURTH COPY IS THE ONE THAT DRIFTS.

    The package list lived in three places -- this doc and the inline
    `create-args` of both workflows -- and two bugs came straight out of that:
    a cross-sections path written twice and a `paths:` filter that covered
    less than its comment claimed. Collapsing to `environment.yml` fixes it
    only if nothing quietly restates the contents, so:

      * both workflows must consume the file, never inline args, and
      * the doc must NAME the file rather than reproduce its dependency list.

    A doc that quotes the list is a fourth copy wearing prose.

    AND THE TRIGGER HAS TO WATCH IT. Collapsing the copies moved the package
    list out from under this workflow's `paths:` filter, which still named
    docs/openmc_stack.md -- the file that used to hold it. The check above
    cannot see that: it derives its list from the regeneration script's
    IMPORTS, and environment.yml is data, not a module. So assert it here,
    where the single spec is the subject.
    """
    import re
    import yaml

    spec = yaml.safe_load(_ENV_FILE.read_text())
    # Split on the first comparison character, not just `=`. `str(d).split("=")`
    # leaves `vtk<9.7` whole, so a version-bounded entry would fail the
    # membership test below for the wrong reason -- reading as "the package is
    # missing" when it is present and merely pinned.
    packages = {re.split(r"[=<>!~ ]", str(d))[0].strip()
                for d in spec["dependencies"]}
    # pyvista, NOT vtk, and deliberately: vtk arrives as pyvista's dependency
    # so that ONE resolver owns the pair. Naming vtk here let conda pick 9.7.0
    # while `pip install -e coupling[dev]` needed pyvista's `<9.7.0`, and pip
    # cannot uninstall a conda-installed package -- see environment.yml, and
    # `test_the_installed_vtk_satisfies_pyvistas_own_requirement` below.
    assert {"openmc", "pyvista"} <= packages, packages

    triggers = _workflow_triggers()
    for event in ("push", "pull_request"):
        listed = set(triggers.get(event, {}).get("paths", []))
        assert "environment.yml" in listed, (
            f"{_WORKFLOW.name}'s `{event}:` paths filter does not list "
            "environment.yml, but the job builds its transport env from that "
            "file -- so changing the OpenMC pin alone would regenerate no "
            "fixture, and the committed tally would stay pinned to a stack "
            "nothing reran")

    # PARSED, not grepped -- the same lesson as the paths-filter check above.
    # A substring scan for "create-args" fires on the COMMENT explaining that
    # create-args was removed, which is a check failing for the opposite of
    # its reason.
    for wf in ("golden-tally-verification.yml", "coupling-tests.yml"):
        doc = yaml.safe_load((_REPO / ".github" / "workflows" / wf).read_text())
        setups = [s for job in doc["jobs"].values() for s in job["steps"]
                  if "setup-micromamba" in str(s.get("uses", ""))]
        assert setups, f"{wf} no longer sets up the transport env"
        for step in setups:
            args = step.get("with", {})
            assert args.get("environment-file") == "environment.yml", (
                f"{wf} does not build the env from environment.yml: {args}")
            assert "create-args" not in args, (
                f"{wf} still builds the env from inline args, so it is a "
                "second spec that can drift from environment.yml")

    doc = _STACK_DOC.read_text()
    assert "environment.yml" in doc, "the doc must name the single spec"
    # THE DOC'S RESOLVED TABLE MUST AGREE WITH THE PIN. This asserted only
    # that the literal `openmc=0.15.3` was absent, while the table two
    # paragraphs down still said `| openmc | 0.15.3 |`: bump the pin in
    # environment.yml and the doc goes stale with the check green. Read both
    # and require the same version, whatever it is.
    pinned = [str(d) for d in spec["dependencies"]
              if re.split(r"[=<>!~ ]", str(d))[0].strip() == "openmc"]
    assert len(pinned) == 1 and "=" in pinned[0], pinned
    pin = pinned[0].split("=", 1)[1].strip()
    table = re.findall(r"^\|\s*openmc\s*\|\s*([^|]+?)\s*\|", doc, re.M)
    assert table, "docs/openmc_stack.md no longer tabulates the resolved openmc"
    assert table == [pin], (
        f"docs/openmc_stack.md tabulates openmc {table} while environment.yml "
        f"pins {pin}; the doc is stale against the single spec")


def test_the_installed_vtk_satisfies_pyvistas_own_requirement():
    """PIP CANNOT UNINSTALL A CONDA PACKAGE.

    Both workflows run `pip install -e "coupling[dev]"` inside the conda env
    built from environment.yml. If the vtk conda resolved falls outside the
    range pyvista declares, pip has to replace it and cannot -- a conda-built
    distribution ships no RECORD, so pip refuses to delete files it has no
    manifest for. The install dies before any test runs:

        error: uninstall-no-record-file
        x Cannot uninstall vtk 9.7.0

    That is not hypothetical. environment.yml named `vtk` unpinned, conda-forge
    served 9.7.0, pyvista 0.48.4 requires `<9.7.0`, and the golden-tally job
    failed on every run from the minute that line landed. Naming `pyvista`
    instead hands the pair to one resolver -- this asserts it stayed handed
    over.

    ASSERTED AGAINST THE LIVE ENVIRONMENT, because neither half of the conflict
    is written down here: conda decides the version, pyvista's own metadata
    decides the range, and this repository states neither. A check that read
    environment.yml would be reading the file that has nothing to say about it.
    """
    import pytest
    from importlib.metadata import PackageNotFoundError, requires, version
    from packaging.requirements import Requirement

    # DISTRIBUTIONS, NOT IMPORTS. `pytest.importorskip("pyvista")` was the
    # first version of this line and it is the wrong gate twice over: pyvista
    # imports a dozen optional things, so one missing transitive dependency
    # turns the guard into a skip -- and a render stack too broken to import
    # is precisely the state worth failing on, not the state worth excusing.
    # The check is metadata-only; it never needs the module. Skip only when
    # the distribution is genuinely absent.
    try:
        installed, pyvista_version = version("vtk"), version("pyvista")
    except PackageNotFoundError as exc:
        pytest.skip(f"{exc.name} is not installed, so the pip/conda seam this "
                    "guards does not exist here -- it exists in every tier "
                    "that installs [dev]")
    wanted = [Requirement(r) for r in (requires("pyvista") or [])]
    wanted = [r for r in wanted
              if r.name == "vtk" and r.marker is None]

    # The control: if pyvista ever stops constraining vtk, this test would pass
    # against anything, which is the shape it exists to refuse.
    assert wanted, (
        "pyvista declares no unconditional vtk requirement, so this check "
        f"proves nothing; its metadata now reads {requires('pyvista')}")

    for req in wanted:
        assert req.specifier.contains(installed, prereleases=True), (
            f"conda installed vtk {installed}, but pyvista "
            f"{pyvista_version} requires {req}. pip will try to replace it "
            "on `pip install -e coupling[dev]` and cannot, so the install "
            "fails before a single test runs. Do not pin vtk in "
            "environment.yml to patch this -- that copies pyvista's ceiling "
            "into a second file by hand. Name only pyvista and let conda "
            "resolve both.")
