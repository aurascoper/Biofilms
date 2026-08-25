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
import json
from pathlib import Path

import numpy as np

from biofilm_openmc.dose import specific_energy_per_source
from biofilm_openmc.feedback_uq import debiased_squared_effect
from biofilm_openmc.synthetic_gate import (PASS_SYNTHETIC_GATE, ThresholdPolicy,
                                           VarianceBudget, decide)

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


def _effect_draws(feedback_heating_scale: float = 1.0) -> np.ndarray:
    """One e_squared per outer draw, exactly the shape decide() consumes."""
    baseline_fields = _per_source_fields("baseline", _MASS_B)
    feedback_fields = _per_source_fields("feedback", _MASS_F,
                                         feedback_heating_scale)
    draws = []
    for o in range(_N_OUTER):
        lo, hi = o * _N_REP, (o + 1) * _N_REP
        effect = debiased_squared_effect(baseline_fields[lo:hi],
                                         feedback_fields[lo:hi],
                                         _MASS_B.ravel())
        draws.append(effect.e_squared)
    return np.array(draws)


def test_pinned_transport_composes_into_a_verdict():
    """The real seam: pinned heating -> dose -> debiased effect -> decide().

    The measured effect is small (~0.0011): a x1.35 density increase raises
    heating by ~26% but raises mass by 35%, so specific energy per kg barely
    moves -- a real, unforced result, not tuned to land anywhere. Recorded
    here as the pinned expectation, the same way `serial_seed42.csv` pins
    what `validate_serial.jl` actually produces.
    """
    draws = _effect_draws()
    assert draws.shape == (_N_OUTER,)
    assert np.all(np.isfinite(draws))
    assert np.all(draws < 0.01), (
        f"effect_draws {draws} drifted far from the pinned expectation; "
        "regenerate the fixture only if OpenMC or the nuclear data changed")

    verdict = decide(draws, VarianceBudget(transport=1e-4), ThresholdPolicy())
    assert verdict.verdict == "EFFECT_BELOW_THRESHOLD", verdict.reason


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
    baseline_verdict = decide(_effect_draws(),
                              VarianceBudget(transport=1e-4),
                              ThresholdPolicy())
    scaled_verdict = decide(_effect_draws(feedback_heating_scale=5.0),
                            VarianceBudget(transport=1e-4),
                            ThresholdPolicy())

    assert scaled_verdict.verdict != baseline_verdict.verdict, (
        "scaling the feedback dose 5x did not change the verdict -- the "
        "pipeline runs but the gate's output does not depend on its input")
    assert scaled_verdict.verdict == PASS_SYNTHETIC_GATE


_REPO = Path(__file__).resolve().parents[2]
_REGEN = _REPO / "coupling" / "scripts" / "regenerate_golden_tally.py"
_WORKFLOW = (_REPO / ".github" / "workflows"
             / "golden-tally-verification.yml")


def _modules_that_produce_the_fixture() -> set[Path]:
    """Every first-party module `regenerate_golden_tally.py` imports, as
    repo-relative paths. Walks the real AST -- including the function-local
    imports inside `_run_one`, which is where all of them live -- so an
    import added later is picked up without anyone remembering to."""
    tree = ast.parse(_REGEN.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
        elif isinstance(node, ast.Import):
            names.update(a.name for a in node.names)

    out: set[Path] = set()
    for name in names:
        parts = name.split(".")
        for base in ("coupling", "coupling/scripts"):
            candidate = _REPO / base / (Path(*parts).as_posix() + ".py")
            if candidate.exists():
                out.add(candidate.relative_to(_REPO))
    return out


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

    for event in ("push", "pull_request"):
        listed = set(triggers.get(event, {}).get("paths", []))
        missing = sorted(p.as_posix() for p in producers
                         if p.as_posix() not in listed)
        assert not missing, (
            f"{_REGEN.name} imports {missing}, so a change to any of them "
            f"changes the real tally -- but {_WORKFLOW.name}'s `{event}:` "
            "paths filter does not list them, so verification would not run "
            "and the committed fixture would go stale while CI stayed green.")


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
    # The version pin is the thing worth restating nowhere: if it appears in
    # the doc, the doc can disagree with the file about which OpenMC is pinned.
    assert "openmc=0.15.3" not in doc, (
        "docs/openmc_stack.md restates the pinned version, which is a fourth "
        "copy of environment.yml's contents -- name the file instead")


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
