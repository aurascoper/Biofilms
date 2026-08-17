"""The pilot's seeding and Omega_b rules, which are what its numbers mean.

BOTH DEFECTS GUARDED HERE WERE SHIPPED AND PUBLISHED. They are not hypothetical
failure modes: they produced a report whose false-positive control was vacuous
and whose two mesh factors were compared over different regions. A unit test is
cheap; re-deriving the finding a second time was not.

No OpenMC needed — the seeding rule and the mask are pure functions, which is
the point. A rule that decides what a measurement MEANS should be checkable
without running the measurement.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "calibration"))
_spec = importlib.util.spec_from_file_location(
    "openmc_nested_pilot", _ROOT / "coupling" / "scripts" / "openmc_nested_pilot.py")
pilot = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pilot)

_rspec = importlib.util.spec_from_file_location(
    "subvoxel_refinement", _ROOT / "coupling" / "scripts" / "subvoxel_refinement.py")
refinement = importlib.util.module_from_spec(_rspec)
_rspec.loader.exec_module(refinement)


def test_unpaired_scenarios_give_the_two_states_different_seeds():
    """THE DEFECT THAT MADE THE FALSE-POSITIVE CONTROL VACUOUS.

    The seed used to depend only on (draw, replicate). With an unchanged
    material the two states were then bit-identical, so the metric returned
    ~1e-16 by construction — at any energy, for any physics. The control could
    not have failed, so it licensed nothing.
    """
    a = pilot._seed(draw=0, rep=0, state="baseline", paired=False)
    b = pilot._seed(draw=0, rep=0, state="feedback", paired=False)
    assert a != b, ("the noise-floor scenario must decorrelate the two states, "
                    "or it measures determinism rather than the noise floor")


def test_paired_scenarios_keep_common_random_numbers():
    """Matched seeds are legitimate variance reduction and must be preserved:
    the levers are measurable only because pairing cancels most of the noise."""
    for state in ("baseline", "feedback"):
        assert (pilot._seed(draw=2, rep=1, state=state, paired=True)
                == pilot._seed(draw=2, rep=1, state="baseline", paired=True))


def test_seeds_stay_distinct_across_draws_and_replicates():
    """The unpaired offset must not collide with another (draw, rep) cell."""
    seeds = [pilot._seed(d, r, s, p)
             for d in range(8) for r in range(4)
             for s in ("baseline", "feedback") for p in (True, False)]
    # baseline/feedback coincide by design when paired, so dedupe that pairing.
    distinct = {pilot._seed(d, r, s, p)
                for d in range(8) for r in range(4)
                for s in ("baseline", "feedback") for p in (True, False)}
    assert len(distinct) == 8 * 4 * 2, "unpaired offset collides with a draw"
    assert len(seeds) > len(distinct)  # the paired duplicates are intentional


def test_omega_b_takes_no_field_at_all():
    """THE DEFECT THAT MADE THE MESH FACTORS INCOMPARABLE, AND ITS SEQUEL.

    Omega_b first required rel_err < 0.25 — a STATISTICAL cut on a REGION
    definition, which removed about a third of the bins at factor 1 and none at
    factor 2, so the factors covered regions differing ~2x in coverage and the
    region moved with the history count.

    It then still required a dose floor taken from REPLICATE 0's field. Same
    error, subtler: the region became a random set drawn from one replicate, so
    of the R(R-1) ordered pairs the U-statistic averages, 2(R-1) involve a row
    that helped choose the weights — four of six at R = 3 — and the
    independence the debiasing rests on does not hold for them.

    The strongest form of the guarantee is that no field is passed in at all: a
    criterion that cannot be reached cannot be applied.
    """
    params = pilot.omega_b.__code__.co_varnames[:pilot.omega_b.__code__.co_argcount]
    assert params == ("mass", "biomass_fraction")
    assert "rel_err" not in pilot.omega_b.__code__.co_varnames
    assert not hasattr(pilot, "DOSE_FLOOR_FRACTION")

    rng = np.random.default_rng(0)
    shape = (4, 4, 4)
    mass = np.full(shape, 0.5)
    frac = rng.uniform(0.0, 1.0, shape)
    mask = pilot.omega_b(mass, frac)
    assert np.array_equal(mask, frac >= pilot.MIN_BIOMASS_VOLUME_FRACTION)


def test_omega_b_still_refuses_massless_bins():
    """The geometric criteria must survive the removal of the statistical ones,
    or the metric starts including bins that hold nothing."""
    shape = (2, 2, 2)
    frac = np.ones(shape)
    mass = np.ones(shape)
    mass[1, 1, 1] = 0.0
    mask = pilot.omega_b(mass, frac)
    assert not mask[1, 1, 1]
    assert mask.sum() == 7


def test_the_lever_table_is_not_hardcoded_anywhere():
    """The six constants that used to be written into the budget artifact
    unconditionally, with no run behind them. Guarded by value, because the
    hazard was the NUMBER surviving, not the variable name.

    Comments are excluded: the module explains what those numbers were and why
    they were removed, and that prose is the record. Only executable code can
    put a transcribed value back into an artifact.
    """
    src = (_ROOT / "coupling" / "scripts" / "openmc_nested_pilot.py").read_text()
    code = "\n".join(ln for ln in src.splitlines()
                     if not ln.lstrip().startswith("#"))
    for stale in ("0.1241", "0.2242", "0.0459", "0.0376", "0.0137", "0.0033"):
        assert stale not in code, (
            f"{stale} is a transcribed lever value; measure it with --levers "
            "so it arrives with a sample row instead")


# ------------------------------------------- publishing an unfinished run

@pytest.mark.parametrize("publish,stopped,expected", [
    (True,  "budget exhausted before scenario 3", "refuse"),
    (True,  None,                                 "allow"),
    (False, "budget exhausted before scenario 3", "allow"),
    (False, None,                                 "allow"),
])
def test_partial_results_may_not_replace_the_canonical_tables(
        publish, stopped, expected):
    """THE ONLY THING STANDING BETWEEN A TIMED-OUT RUN AND THE PUBLISHED
    EVIDENCE, and until now it had no test.

    `--publish` writes to `data/calibration/`, replacing the canonical tables.
    The scenario-completeness check upstream confirms the full set was
    REQUESTED; it cannot know whether the run finished. A budget exhausted
    partway sets `stopped_early`, and publishing then swaps complete evidence
    for partial rows.

    The guard lived inside `_report`, which imports openmc on its first line and
    therefore cannot run anywhere in this suite — so the protection existed and
    was never once exercised. All four combinations are checked, not just the
    refusing one: a guard that raised unconditionally would also stop the bad
    case, and would make the pilot unpublishable.
    """
    if expected == "refuse":
        with pytest.raises(SystemExit, match="stopped early"):
            pilot.refuse_partial_publish(publish, stopped)
    else:
        assert pilot.refuse_partial_publish(publish, stopped) is None


# ------------------------------------ declaring significance with no test

@pytest.mark.parametrize("df,expected", [
    (0, None),      # a single outer draw: nothing to test with
    (-1, None),
    (1, 318.3),     # and the value that makes the refusal matter
    (2, 22.33),
    (100, 3.232),   # untabulated: the STRICTEST df not exceeding it
])
def test_significance_is_refused_where_there_are_no_degrees_of_freedom(
        df, expected):
    """A ONE-DRAW SCENARIO MUST NOT BE CALLED DISTINGUISHABLE FROM ZERO.

    With `m` outer draws the test has `m - 1` degrees of freedom, so a single
    draw has none. Falling back to the normal critical value 3.09 understates
    the truth at every finite df — at df = 1 it is 318.3, a factor of 103 — so
    a scenario with no degrees of freedom would be declared significant on a
    test that cannot be performed.

    This lived nested inside `_report`, whose first statement imports openmc.
    The refusal could have been deleted and no suite in this repository would
    have noticed.
    """
    assert pilot.t_critical_999(df) == expected


def test_the_low_df_refusal_is_not_merely_a_smaller_number():
    """The value returned for df = 1 must be enormous, or the refusal is
    decoration: a critical value near 3.09 would let a one-draw scenario pass
    on noise. 318.3 is the actual one-sided t at 0.999."""
    assert pilot.t_critical_999(1) > 100
    # and it must fall monotonically toward the normal limit, never below it
    values = [pilot.t_critical_999(d) for d in sorted(pilot.T_CRIT_999)]
    assert values == sorted(values, reverse=True)
    assert min(values) > 3.09


# ------------------------------------------- ratios that cannot share a mass

@pytest.mark.parametrize("ratios,accepted", [
    ([1, 2, 4], True),
    ([1, 2], True),
    ([1, 3, 4], False),     # 3 does not divide 4
    ([2, 3, 6], True),
    ([1, 4, 6], False),     # 4 does not divide 6
])
def test_a_ratio_that_cannot_share_the_mass_denominator_is_refused(
        ratios, accepted):
    """THE ONLY THING STANDING BETWEEN A BAD RATIO AND A PAID-FOR RAYTRACE.

    The mass denominator is raytraced once at the finest ratio and coarsened by
    `finest // ratio`. Integer division silently yields step 1 for a ratio of 3
    under a finest of 4, pairing a 4x mass array with a 3x tally — and the
    mismatch surfaces only when dose conversion combines them, after the
    raytrace and the histories are spent.

    The guard lived in `main()` below `import openmc`, so it was unreachable
    from every suite that runs here. Both directions are checked: a guard that
    refused everything would also stop `1,2,4`, which is the ladder the study
    actually runs.
    """
    if accepted:
        assert refinement.refuse_non_divisor_ratios(ratios) is None
    else:
        with pytest.raises(SystemExit, match="do not divide the finest"):
            refinement.refuse_non_divisor_ratios(ratios)


# ------------------------- the guards, reached the way production reaches them

def test_the_canonical_tables_are_unreachable_after_a_partial_run(tmp_path):
    """THE REFUSAL IS THE DOOR, not a sign beside it.

    `refuse_partial_publish` alone proves the predicate, not the wiring: delete
    its call and every helper-level test stays green while `_report` reopens
    data/calibration/ in write mode after a budget-exhausted run. So ask the
    question production asks — "where do the tables go?" — and require that the
    canonical answer cannot be obtained without passing the guard.
    """
    from types import SimpleNamespace

    canonical = pilot.REPO / "data" / "calibration"

    with pytest.raises(SystemExit, match="stopped early"):
        pilot.resolve_output_dir(
            SimpleNamespace(publish=True, outdir=tmp_path),
            "budget exhausted before scenario 3")

    # AND THE OTHER ROUTE TO THE SAME DIRECTORY. `--publish` is one way in;
    # `--outdir data/calibration` is another, and it reached the canonical
    # tables without passing any check at all. The guard is on the destination
    # now, so every path to that directory goes through it.
    for outdir in (canonical, canonical / "sub"):
        with pytest.raises(SystemExit, match="stopped early"):
            pilot.resolve_output_dir(
                SimpleNamespace(publish=False, outdir=outdir),
                "budget exhausted before scenario 3")

    # the same run pointed somewhere harmless is fine
    assert pilot.resolve_output_dir(
        SimpleNamespace(publish=False, outdir=tmp_path),
        "budget exhausted before scenario 3") == tmp_path.resolve()
    # and a COMPLETE run may publish, or the guard would make the pilot useless
    assert pilot.resolve_output_dir(
        SimpleNamespace(publish=True, outdir=tmp_path), None) == canonical


@pytest.mark.parametrize("m,e2,se,expected", [
    (1, 10.0, 0.001, False),   # one outer draw: no degrees of freedom at all
    (0, 10.0, 0.001, False),
    (2, 10.0, 0.001, True),    # df 1, crit 318.3, and 10/0.001 clears it
    (2, 10.0, 1.0,   False),   # df 1: 10 does NOT clear 318.3
    (31, 10.0, 1.0,  True),    # df 30, crit 3.385
    (31, 10.0, 0.0,  False),   # zero standard error is not infinite confidence
])
def test_the_verdict_itself_refuses_when_there_are_no_degrees_of_freedom(
        m, e2, se, expected):
    """`t_critical_999(0) is None` proves the table. It does NOT prove the None
    reaches the verdict — swapping the call site back to the 3.09 fallback left
    every table-level test green while a one-draw scenario became reportable as
    significant again.

    The row at m=2, se=1.0 is the one that matters most: it FAILS, because at
    one degree of freedom the critical value is 318.3 rather than 3.09. If the
    normal fallback ever returns, that row flips to True.
    """
    assert pilot.distinguishable_from_zero(e2, se, m) is expected


@pytest.mark.parametrize("ratios,message", [
    ("1,3,4", "do not divide the finest"),
    ("1,4,6", "do not divide the finest"),
    ("2,4",   "ratio 1 is the reference grid"),
])
def test_main_refuses_a_bad_ratio_before_it_imports_openmc(
        ratios, message, tmp_path):
    """THE WIRING AND THE ORDERING, not just the predicate.

    `refuse_non_divisor_ratios` alone proves the rule; delete its call from
    `main()` and the suite stays green while the run proceeds to the finest-grid
    raytrace and only then combines incompatible mass and tally shapes.

    The check sat BELOW `import openmc` and below `load_snapshot`, so it fired
    once the heavy stack was already up — and never at all in this tier, since
    the import fails first. Moving it above the import is what makes the
    refusal cheap, and this test is what proves it moved: openmc is not
    installed here, so reaching the import at all raises ModuleNotFoundError
    rather than SystemExit. The paths below do not exist either, for the same
    reason — nothing may be opened before the ratios are judged.
    """
    with pytest.raises(SystemExit, match=message):
        refinement.main(["--snapshot", str(tmp_path / "nope.h5"),
                         "--config", str(tmp_path / "nope.toml"),
                         "--outdir", str(tmp_path / "out"),
                         "--ratios", ratios])
    assert not (tmp_path / "out").exists(), (
        "main() created its output directory before judging the ratios")


def test_main_accepts_the_ladder_the_study_actually_runs(tmp_path):
    """The control that keeps the one above honest: `1,2,4` must get PAST the
    ratio check. A guard that refused every ladder would satisfy every test in
    this file and make the study unrunnable — so require the refusal to end and
    the next stage to begin, which here is the missing openmc stack."""
    with pytest.raises((ModuleNotFoundError, ImportError, FileNotFoundError)):
        refinement.main(["--snapshot", str(tmp_path / "nope.h5"),
                         "--config", str(tmp_path / "nope.toml"),
                         "--outdir", str(tmp_path / "out"),
                         "--ratios", "1,2,4"])


@pytest.mark.parametrize("publish,to_canonical,complete,refuses", [
    (True,  False, False, True),    # --publish with the default scenario set
    (True,  False, True,  False),
    # THE ALIAS. A COMPLETED run with --outdir data/calibration and no
    # --publish: the stopped-early guard permits it, and the completeness check
    # used to skip it because args.publish was false. It would then overwrite
    # the canonical CSVs with the default, partial scenario set and drop the
    # rows that published claims cite.
    (False, True,  False, True),
    (False, True,  True,  False),
    (False, False, False, False),   # somewhere harmless: nothing to protect
])
def test_a_partial_scenario_set_cannot_replace_a_complete_one(
        publish, to_canonical, complete, refuses, tmp_path):
    """EVERY PUBLICATION PREREQUISITE KEYS ON THE DESTINATION.

    `--publish` is one route into data/calibration/; `--outdir` naming it is
    another. Guarding the flag protects one of them, and this repository has
    now made that mistake twice in the same function — once for the
    stopped-early check and once for this one.
    """
    from types import SimpleNamespace

    args = SimpleNamespace(
        publish=publish,
        outdir=(pilot.CANONICAL_TABLES if to_canonical else tmp_path),
        levers=complete, include_near_threshold=complete)

    if refuses:
        with pytest.raises(SystemExit, match="canonical tables"):
            pilot.refuse_incomplete_scenarios(args)
    else:
        assert pilot.refuse_incomplete_scenarios(args) is None


def test_main_refuses_an_incomplete_scenario_set_before_it_imports_openmc(
        tmp_path):
    """THE WIRING, AGAIN. Same shape as the ratio guard.

    `refuse_incomplete_scenarios` alone proves the predicate; delete its call
    from `main()` and the suite stays green while a run with the default
    scenario set proceeds to overwrite the canonical tables.

    openmc is not installed here, so reaching the import raises
    ModuleNotFoundError rather than SystemExit — which makes this a test of the
    ORDERING as much as of the call.
    """
    with pytest.raises(SystemExit, match="canonical tables"):
        pilot.main(["--snapshot", str(tmp_path / "nope.h5"),
                    "--config", str(tmp_path / "nope.toml"),
                    "--outdir", str(pilot.CANONICAL_TABLES)])

    # and the complete set gets PAST it, or the pilot could never publish
    with pytest.raises((ModuleNotFoundError, ImportError, FileNotFoundError)):
        pilot.main(["--snapshot", str(tmp_path / "nope.h5"),
                    "--config", str(tmp_path / "nope.toml"),
                    "--outdir", str(pilot.CANONICAL_TABLES),
                    "--levers", "--include-near-threshold"])
