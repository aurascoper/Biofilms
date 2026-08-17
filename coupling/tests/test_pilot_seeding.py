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
