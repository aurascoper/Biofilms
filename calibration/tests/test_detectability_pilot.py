"""The detectability pilot, and the four things it refuses.

The load-bearing assertion here is `test_replicate_sizing_is_not_attempted_...`:
if this script ever returns a replicate count for a coupon whose signal is under
its blank's noise, it will have designed a campaign that cannot succeed, and the
number will look like a plan.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "calibration" / "scripts"))

import detectability_pilot as pilot  # noqa: E402

_BLANK = {"blank_sample_id": "B1", "wet_mass_g": 0.100,
          "wet_mass_uncertainty_g": 0.0002}


def _row(**over):
    row = {
        "sample_id": "S1", "replicate_id": "1", "blank_sample_id": "B1",
        "culture_batch_id": "batch-a",
        "wet_mass_sample_plus_substrate_g": 0.140,
        "wet_mass_uncertainty_g": 0.0002,
        "hydrated_volume_cm3": 0.05,
        "volume_basis": "whole_biofilm_envelope",
        "volume_support": "full_coupon",
        "scaling_method": "", "scaling_uncertainty": None,
        "imaged_area_cm2": None, "weighed_area_cm2": None,
    }
    row.update(over)
    return row


# --- refusal 1: the volume must COVER what the balance weighed -----------

def test_unresolved_support_is_refused():
    problems = pilot.support_problems(_row(volume_support="unresolved"))
    assert any("volume_support is unresolved" in p for p in problems)


def test_scaled_volume_needs_its_method_and_uncertainty():
    problems = pilot.support_problems(_row(volume_support="stereological_scaling"))
    joined = " ".join(problems)
    assert "scaling_method" in joined
    assert "scaling_uncertainty" in joined
    # and it must say what was scaled to what
    assert "imaged_area_cm2" in joined


def test_full_coupon_support_passes():
    assert pilot.support_problems(_row()) == []


# --- refusal 2: the volume must ENCLOSE what the balance weighed ---------

def test_cells_only_volume_is_refused_under_a_whole_biofilm_mass():
    problems = pilot.basis_problems(_row(volume_basis="cells_only"))
    assert problems and "cells_only" in problems[0]


def test_cells_and_matrix_is_accepted_only_with_a_pore_fraction():
    assert pilot.basis_problems(_row(volume_basis="cells_and_matrix"))
    ok = pilot.basis_problems(_row(volume_basis="cells_and_matrix",
                                   pore_volume_fraction=0.3))
    assert ok == []


# --- refusals 3 and 4, and the sizing that must not happen ---------------

def _evaluate(row, blank=None):
    return pilot.evaluate_row(row, {"B1": blank or _BLANK},
                              minimum_ratio=10.0, seed=1, draws=4000)


def test_a_clean_coupon_clears_every_stage():
    out = _evaluate(_row())
    assert out["support_status"] == "PASS"
    assert out["basis_status"] == "PASS"
    assert out["detectability_status"] == "PASS"
    assert out["distribution_status"] == "VALID"
    assert out["refusals"] == ""
    assert out["biofilm_wet_mass_g"] == pytest.approx(0.040)


def test_signal_under_the_blank_noise_is_refused():
    # 2 mg of biofilm against a 1 mg blank sigma: a ratio of 2, target 10.
    out = _evaluate(_row(wet_mass_sample_plus_substrate_g=0.102),
                    blank={**_BLANK, "wet_mass_uncertainty_g": 0.001})
    assert out["detectability_status"] == "FAIL"
    assert "not a replicate-count problem" in out["refusals"]


def test_a_distribution_crossing_zero_is_not_identifiable():
    """The case the ratio test alone cannot catch.

    `detectability()` compares the point estimate against the BLANK's sigma, so
    a tiny blank uncertainty gives a reassuring 40x. But the difference also
    carries the SAMPLE's own weighing uncertainty, and here that is wide enough
    to put a material share of the distribution at or below zero. A 40x ratio
    and an unidentifiable mass at the same time is exactly why both checks run.
    """
    out = _evaluate(_row(wet_mass_sample_plus_substrate_g=0.140,
                         wet_mass_uncertainty_g=0.025),
                    blank={**_BLANK, "wet_mass_uncertainty_g": 0.001})
    assert out["detectability_status"] == "PASS"
    assert out["signal_to_blank_noise"] == pytest.approx(40.0, rel=0.01)
    assert out["distribution_status"] == "NOT_IDENTIFIABLE"
    assert out["negative_blank_corrected_mass_probability"] > \
        pilot.MAX_INVALID_DRAW_FRACTION


def test_an_unresolvable_blank_is_refused():
    out = pilot.evaluate_row(_row(blank_sample_id="MISSING"), {},
                             minimum_ratio=10.0, seed=1, draws=100)
    assert "does not resolve" in out["refusals"]


def test_replicate_sizing_is_not_attempted_on_undetectable_coupons():
    """The one that matters. A replicate count returned here would read as a
    campaign plan for a measurement the apparatus cannot make."""
    bad = [_evaluate(_row(sample_id=f"S{i}",
                          wet_mass_sample_plus_substrate_g=0.102),
                     blank={**_BLANK, "wet_mass_uncertainty_g": 0.001})
           for i in range(4)]
    assert all(r["distribution_status"] != "VALID" for r in bad)
    assert pilot.size_replicates(bad, target_rel_uncertainty=0.1,
                                 n_within=1) is None


def test_replicate_sizing_runs_once_coupons_clear():
    good = [_evaluate(_row(sample_id=f"S{i}", culture_batch_id=f"b{i % 2}",
                           wet_mass_sample_plus_substrate_g=0.140 + 0.002 * i))
            for i in range(4)]
    plan = pilot.size_replicates(good, target_rel_uncertainty=0.10, n_within=1)
    assert plan is not None
    assert plan["required_batches"] >= 2      # precision.py's floor
    assert plan["dominant_component"] in ("between_batch", "within_batch_effective",
                                          "blank", "none")


def test_the_demo_runs_and_exercises_all_four_refusals(capsys):
    assert pilot.demo_refusals() == 0
    out = capsys.readouterr().out
    assert "INCOMMENSURATE SUPPORT" in out
    assert "INCOMPATIBLE BASIS" in out
    assert "SIGNAL UNDER THE BLANK" in out
    assert "DISTRIBUTION CROSSES ZERO" in out


def test_blank_uncertainty_kind_has_no_default():
    """Wet, dry and blank masses share a balance calibration, so a sigma whose
    meaning is undeclared cannot be propagated honestly."""
    with pytest.raises(SystemExit):
        pilot.main(["--out", "/dev/null"])
