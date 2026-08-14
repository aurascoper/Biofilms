"""Detectability and replicate count.

Two questions that get conflated, at the cost of a whole campaign: whether the
measurement is possible at all, and how many times to repeat it.
"""

from __future__ import annotations

import pytest

from biofilm_calibration.precision import (Detectability, PrecisionError,
                                           detectability, plan_problems,
                                           required_replicates)


# --- detectability ---------------------------------------------------------

def test_a_mass_buried_under_the_blank_is_not_a_replicate_problem():
    """THE detectability test. Replicates shrink the standard error of a mean;
    they do not recover a signal buried under a subtraction. The fix is the
    substrate or the biomass."""
    d = detectability(biofilm_mass_g=0.002, blank_uncertainty_g=0.001)
    assert not d.detectable
    assert d.signal_to_blank_noise == pytest.approx(2.0)
    assert "SUBSTRATE GEOMETRY" in d.reason
    assert "not a replicate-count problem" in d.reason


def test_an_adequate_mass_passes():
    d = detectability(biofilm_mass_g=0.05, blank_uncertainty_g=0.001)
    assert d.detectable
    assert d.signal_to_blank_noise == pytest.approx(50.0)


def test_the_ratio_target_is_declared_not_fixed():
    assert detectability(0.02, 0.001, minimum_ratio=10.0).detectable
    assert not detectability(0.02, 0.001, minimum_ratio=50.0).detectable


def test_a_zero_blank_uncertainty_is_flagged_not_trusted():
    d = detectability(1.0, 0.0)
    assert d.detectable
    assert "measured rather than assumed" in d.reason


def test_detectability_refuses_nonsense():
    with pytest.raises(PrecisionError, match="must be positive"):
        detectability(0.0, 0.001)
    with pytest.raises(PrecisionError, match="non-negative"):
        detectability(1.0, -1.0)


# --- replicate count -------------------------------------------------------

def test_more_between_batch_variance_needs_more_batches():
    """n is derived, not declared. n = 3 is a convention, not a calculation."""
    kw = dict(target_rel_uncertainty=0.05, mean_value=1.0,
              var_within_batch=0.0, n_within=1)
    low = required_replicates(var_between_batch=0.001, **kw)
    high = required_replicates(var_between_batch=0.01, **kw)
    assert high.required_batches > low.required_batches


def test_the_dominant_variance_component_is_named():
    """Which component dominates says WHERE to spend effort, and the three
    answers are different actions."""
    between = required_replicates(
        target_rel_uncertainty=0.05, mean_value=1.0,
        var_between_batch=0.01, var_within_batch=0.0001, n_within=4)
    assert between.dominant_component == "between_batch"
    assert "more independent cultures" in between.note

    within = required_replicates(
        target_rel_uncertainty=0.05, mean_value=1.0,
        var_between_batch=0.0001, var_within_batch=0.01, n_within=1)
    assert within.dominant_component == "within_batch_effective"
    assert "more fields per coupon are cheaper" in within.note

    blank = required_replicates(
        target_rel_uncertainty=0.05, mean_value=1.0,
        var_between_batch=0.0001, var_within_batch=0.0001, var_blank=0.02)
    assert blank.dominant_component == "blank"
    assert "lighter or more reproducible substrate" in blank.note


def test_fields_per_batch_only_help_within_batch_variance():
    """The asymmetry that makes the dominant component worth reporting."""
    kw = dict(target_rel_uncertainty=0.05, mean_value=1.0,
              var_between_batch=0.01, var_within_batch=0.01)
    one = required_replicates(n_within=1, **kw)
    many = required_replicates(n_within=100, **kw)
    assert many.required_batches < one.required_batches
    # ... but between-batch variance sets a floor that fields cannot cross
    only_between = dict(target_rel_uncertainty=0.05, mean_value=1.0,
                        var_between_batch=0.01, var_within_batch=0.0)
    assert (required_replicates(n_within=1, **only_between).required_batches
            == required_replicates(n_within=1000, **only_between).required_batches)


def test_a_tighter_target_needs_more_batches_quadratically():
    """Halving the target uncertainty quadruples n — away from the n >= 2
    floor, which a low-variance case hits immediately."""
    kw = dict(mean_value=1.0, var_between_batch=0.1, var_within_batch=0.0)
    loose = required_replicates(target_rel_uncertainty=0.10, **kw)
    tight = required_replicates(target_rel_uncertainty=0.05, **kw)
    assert loose.required_batches == 10
    assert tight.required_batches == 40


def test_at_least_two_batches_are_always_required():
    """A variance estimate needs two."""
    plan = required_replicates(target_rel_uncertainty=0.5, mean_value=1.0,
                               var_between_batch=1e-9, var_within_batch=0.0)
    assert plan.required_batches >= 2


def test_an_unreachable_target_is_capped_and_says_so():
    plan = required_replicates(target_rel_uncertainty=0.001, mean_value=1.0,
                               var_between_batch=1.0, var_within_batch=0.0,
                               max_batches=50)
    assert plan.required_batches == 50
    assert "CAPPED" in plan.note
    assert plan_problems(plan, 0.001), "an unmet target must be a gate problem"


def test_a_met_target_reports_no_problem():
    plan = required_replicates(target_rel_uncertainty=0.05, mean_value=1.0,
                               var_between_batch=0.001, var_within_batch=0.0)
    assert plan_problems(plan, 0.05) == []
    assert plan.achieved_rel_uncertainty <= 0.05


def test_required_replicates_refuses_nonsense():
    kw = dict(mean_value=1.0, var_between_batch=0.01, var_within_batch=0.0)
    with pytest.raises(PrecisionError, match=r"in \(0, 1\)"):
        required_replicates(target_rel_uncertainty=0.0, **kw)
    with pytest.raises(PrecisionError, match="must be positive"):
        required_replicates(target_rel_uncertainty=0.05, mean_value=0.0,
                            var_between_batch=0.01, var_within_batch=0.0)
    with pytest.raises(PrecisionError, match="non-negative"):
        required_replicates(target_rel_uncertainty=0.05, mean_value=1.0,
                            var_between_batch=-1.0, var_within_batch=0.0)


# --- the selected dynamic observable ---------------------------------------

def test_the_dynamic_observable_is_now_selected(data_dir):
    """Selected before culturing, so the campaign knows it needs time-resolved
    stacks rather than discovering it after a static one."""
    from biofilm_calibration.schema import read_table
    from biofilm_calibration.spatial.time_observable import (TIME_OBSERVABLE,
                                                             problems, selectable)
    rows = read_table(data_dir / "spatial" / "time_observable.csv",
                      TIME_OBSERVABLE)
    chosen = [r for r in rows if r["selected"] == "true"]
    assert len(chosen) == 1
    assert chosen[0]["observable_id"] == "biomass_autocorrelation_decay"
    assert selectable(chosen[0])
    assert problems(rows) == []


def test_the_selected_observable_is_field_level_on_both_sides(data_dir):
    """The property that disqualified single-cell MSD: no cell-to-parcel
    mapping is needed."""
    from biofilm_calibration.schema import read_table
    from biofilm_calibration.spatial.time_observable import TIME_OBSERVABLE
    rows = {r["observable_id"]: r for r in read_table(
        data_dir / "spatial" / "time_observable.csv", TIME_OBSERVABLE)}
    chosen = rows["biomass_autocorrelation_decay"]
    assert chosen["semantically_compatible_with_parcel"] == "true"
    assert "no cell-to-parcel mapping" in chosen["notes"]
    assert rows["single_cell_msd"]["selected"] == "false"
