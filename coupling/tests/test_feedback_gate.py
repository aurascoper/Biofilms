"""The feedback gates, and the many ways they must refuse.

The load-bearing property of this file is that EVERY path out of the decision
engine that is not an explicit pass is a block. A gate that defaults to
permissive on a missing input is worse than no gate: it converts an absence of
evidence into an authorization.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_openmc.feedback_gate import (BELOW_THRESHOLD, BLOCKED_BIOLOGICAL,
                                          BLOCKED_NUMERICAL, BLOCKED_PHYSICAL,
                                          NOT_EVALUATED, OFFLINE_PASS,
                                          OFFLINE_UNCERTAIN, ONLINE_ENABLED,
                                          ONLINE_OUT_OF_DOMAIN,
                                          PRODUCTION_MASS_DENOMINATOR,
                                          UNSUPPORTED, EffectPolicy,
                                          OnlinePrerequisites, UncertaintyBudget,
                                          evaluate_offline, evaluate_online)
from biofilm_openmc.feedback_uq import (PairedEffect, absorbed_energy_rate,
                                        common_random_number_benefit,
                                        decompose_variance, effect_draws,
                                        mass_weighted_mean_dose,
                                        seed_sufficiency)

POLICY = EffectPolicy(qoi_id="dQ_E", effect_direction="increase",
                      effect_threshold=0.05, policy_id="test_v1",
                      declared_before_results=True)
BUDGET = UncertaintyBudget(transport_numerical_u99=0.005,
                           material_volume_u99=0.002, surrogate_error=0.0,
                           mass_denominator_method=PRODUCTION_MASS_DENOMINATOR)


def _draws(mean, sd, n=20000, seed=1):
    return np.random.default_rng(seed).normal(mean, sd, n)


# --- the threshold is the scientist's, not the gate's --------------------

def test_an_undeclared_threshold_is_not_evaluated():
    """No uncertainty calculation can choose what magnitude of feedback
    matters. A gate that supplies its own relevance threshold has decided the
    science it was built to test."""
    naked = EffectPolicy(qoi_id="q", effect_direction="increase",
                         policy_id="v1", declared_before_results=True)
    v = evaluate_offline(_draws(0.5, 0.01), naked, BUDGET,
                         physical_calibration_ready=True)
    assert v.verdict == NOT_EVALUATED
    assert any("effect_threshold is undeclared" in r for r in v.reasons)


def test_a_threshold_declared_after_results_is_refused():
    late = EffectPolicy(qoi_id="q", effect_direction="increase",
                        effect_threshold=0.05, policy_id="v1",
                        declared_before_results=False)
    v = evaluate_offline(_draws(0.5, 0.01), late, BUDGET,
                         physical_calibration_ready=True)
    assert v.verdict == NOT_EVALUATED
    assert any("not a threshold" in r for r in v.reasons)


def test_direction_is_required_because_a_two_sided_test_passes_on_wrong_sign():
    wrong = EffectPolicy(qoi_id="q", effect_direction="either",
                         effect_threshold=0.05, policy_id="v1",
                         declared_before_results=True)
    assert evaluate_offline(_draws(0.5, 0.01), wrong, BUDGET,
                            physical_calibration_ready=True).verdict == NOT_EVALUATED


def test_an_effect_of_the_wrong_sign_does_not_pass():
    """A decrease of the same magnitude is not the declared effect."""
    v = evaluate_offline(_draws(-0.5, 0.01), POLICY, BUDGET,
                         physical_calibration_ready=True)
    assert v.verdict != OFFLINE_PASS


# --- the three outcomes that are not a pass ------------------------------

def test_a_large_well_resolved_effect_passes():
    v = evaluate_offline(_draws(0.20, 0.02), POLICY, BUDGET,
                         physical_calibration_ready=True)
    assert v.verdict == OFFLINE_PASS
    assert v.effect_q001 > POLICY.effect_threshold
    assert v.probability_exceeds > 0.99


def test_too_small_and_too_uncertain_are_different_verdicts():
    """Conflating them hides which one more histories can fix. A small effect
    stays small however well it is resolved; a wide one may not."""
    small = evaluate_offline(_draws(0.01, 0.005), POLICY, BUDGET,
                             physical_calibration_ready=True)
    wide = evaluate_offline(_draws(0.12, 0.10), POLICY, BUDGET,
                            physical_calibration_ready=True)
    assert small.verdict == BELOW_THRESHOLD
    assert wide.verdict == OFFLINE_UNCERTAIN
    # The wide one's MEAN clears the threshold; its lower bound does not.
    assert wide.effect_mean > POLICY.effect_threshold > wide.effect_q001


# --- fail closed on every missing input ----------------------------------

def test_uncalibrated_material_blocks_before_any_statistics():
    v = evaluate_offline(_draws(0.5, 0.01), POLICY, BUDGET)
    assert v.verdict == BLOCKED_PHYSICAL


def test_unresolved_numerics_block_a_pass_that_would_otherwise_happen():
    """A pass driven by barely resolved numerics is a pass about the mesh."""
    sloppy = UncertaintyBudget(transport_numerical_u99=0.04,   # > 0.25 * 0.05
                               material_volume_u99=0.002, surrogate_error=0.0,
                               mass_denominator_method=PRODUCTION_MASS_DENOMINATOR)
    v = evaluate_offline(_draws(0.20, 0.02), POLICY, sloppy,
                         physical_calibration_ready=True)
    assert v.verdict == BLOCKED_NUMERICAL
    assert any("25%" in r for r in v.reasons)


def test_the_full_bin_denominator_is_a_negative_control_not_an_alternative():
    """It weighs the circumscribing cube — 4/pi too much on A0. A model
    ensemble that averages it with the exact method averages in a known error."""
    v = evaluate_offline(_draws(0.20, 0.02), POLICY,
                         UncertaintyBudget(0.005, 0.002, 0.0, "full_bin"),
                         physical_calibration_ready=True)
    assert v.verdict == BLOCKED_NUMERICAL
    assert any("negative control" in r for r in v.reasons)


def test_an_uncharacterised_budget_blocks_rather_than_assuming_zero():
    v = evaluate_offline(_draws(0.20, 0.02), POLICY,
                         UncertaintyBudget(mass_denominator_method=PRODUCTION_MASS_DENOMINATOR),
                         physical_calibration_ready=True)
    assert v.verdict == BLOCKED_NUMERICAL
    assert any("not characterised" in r for r in v.reasons)


def test_an_unresolved_model_branch_blocks():
    v = evaluate_offline(_draws(0.20, 0.02), POLICY,
                         UncertaintyBudget(0.005, 0.002, 0.0,
                                           PRODUCTION_MASS_DENOMINATOR,
                                           ("occupancy_mapping",)),
                         physical_calibration_ready=True)
    assert v.verdict == BLOCKED_NUMERICAL
    assert any("averaging over an undeclared choice" in r for r in v.reasons)


def test_empty_draws_are_not_evaluated_rather_than_passed():
    v = evaluate_offline([], POLICY, BUDGET, physical_calibration_ready=True)
    assert v.verdict == NOT_EVALUATED


# --- the online gate is an authorization, not a measurement --------------

def _passing_offline():
    return evaluate_offline(_draws(0.20, 0.02), POLICY, BUDGET,
                            physical_calibration_ready=True)


def test_online_defaults_are_all_blocking():
    """Every prerequisite defaults to the refusing value, so a caller that
    forgets one cannot accidentally authorize feedback."""
    v = evaluate_online(OnlinePrerequisites(), POLICY, BUDGET)
    assert v.verdict != ONLINE_ENABLED


def test_online_requires_the_offline_gate_to_have_passed():
    blocked = evaluate_offline(_draws(0.01, 0.005), POLICY, BUDGET,
                               physical_calibration_ready=True)
    v = evaluate_online(OnlinePrerequisites(offline=blocked), POLICY, BUDGET)
    assert v.verdict != ONLINE_ENABLED


def test_online_blocks_on_missing_biological_calibration():
    v = evaluate_online(OnlinePrerequisites(
        offline=_passing_offline(), physical_calibration_ready=True,
        numerical_resolution_ready=True, response_pathway="hamiltonian"),
        POLICY, BUDGET)
    assert v.verdict == BLOCKED_BIOLOGICAL


def test_online_fails_closed_outside_the_validated_envelope():
    """Offline evidence does not license extrapolation."""
    v = evaluate_online(OnlinePrerequisites(
        offline=_passing_offline(), physical_calibration_ready=True,
        numerical_resolution_ready=True, biological_posterior_ready=True,
        seconds_per_mcs_ready=True, state_within_validity_domain=False,
        response_pathway="hamiltonian"), POLICY, BUDGET)
    assert v.verdict == ONLINE_OUT_OF_DOMAIN


@pytest.mark.parametrize("pathway", ["growth_survival", "membrane_response"])
def test_unsupported_pathways_can_never_be_authorized(pathway):
    """The CPM has no birth and no death, and the membrane constitutive law is
    unresolved. No amount of evidence authorizes a response the model cannot
    represent."""
    v = evaluate_online(OnlinePrerequisites(
        offline=_passing_offline(), physical_calibration_ready=True,
        numerical_resolution_ready=True, biological_posterior_ready=True,
        seconds_per_mcs_ready=True, state_within_validity_domain=True,
        response_pathway=pathway), POLICY, BUDGET)
    assert v.verdict == UNSUPPORTED


def test_everything_ready_enables_exactly_one_verdict():
    v = evaluate_online(OnlinePrerequisites(
        offline=_passing_offline(), physical_calibration_ready=True,
        numerical_resolution_ready=True, biological_posterior_ready=True,
        seconds_per_mcs_ready=True, state_within_validity_domain=True,
        response_pathway="hamiltonian"), POLICY, BUDGET)
    assert v.verdict == ONLINE_ENABLED


def test_a_live_effect_that_stopped_clearing_reblocks():
    v = evaluate_online(OnlinePrerequisites(
        offline=_passing_offline(), physical_calibration_ready=True,
        numerical_resolution_ready=True, biological_posterior_ready=True,
        seconds_per_mcs_ready=True, state_within_validity_domain=True,
        response_pathway="hamiltonian"), POLICY, BUDGET,
        current_effect_draws=_draws(0.005, 0.002))
    assert v.verdict != ONLINE_ENABLED


# --- the paired-effect machinery ----------------------------------------

def test_the_source_rate_cancels_from_a_relative_effect():
    """Which is why the offline material screen can run per source particle,
    before any assay certificate exists."""
    mass = np.ones(8)
    base, feed = np.full(8, 2.0), np.full(8, 2.4)
    for rate in (1.0, 1e12):
        e = PairedEffect("q", np.array([absorbed_energy_rate(base * rate, mass)]),
                         np.array([absorbed_energy_rate(feed * rate, mass)]))
        assert e.relative_effect == pytest.approx(0.2)


def test_variance_decomposition_separates_what_seeds_can_fix():
    rng = np.random.default_rng(3)
    # Wide spread ACROSS outer points, tight within: parameter-dominated, so
    # more seeds would buy nothing and better measurements would buy everything.
    effects = [PairedEffect("q", np.zeros(5), rng.normal(mu, 0.001, 5))
               for mu in rng.normal(0.2, 0.05, 40)]
    d = decompose_variance(effects)
    assert d.parameter > d.transport
    assert d.transport_share < 0.1
    assert "only measurements help" in d.render()

    # And the reverse.
    noisy = [PairedEffect("q", np.zeros(5), rng.normal(0.2, 0.05, 5))
             for _ in range(40)]
    assert decompose_variance(noisy).transport_share > 0.5


def test_effect_draws_feed_the_gate_directly():
    effects = [PairedEffect("q", np.zeros(3), np.full(3, 0.2)) for _ in range(30)]
    v = evaluate_offline(effect_draws(effects), POLICY, BUDGET,
                         physical_calibration_ready=True)
    assert v.verdict == OFFLINE_PASS


def test_seed_sufficiency_says_when_a_variance_is_its_own_evidence():
    assert seed_sufficiency(1)["sd_relative_error"] == float("inf")
    assert seed_sufficiency(5)["sd_relative_error"] == pytest.approx(0.354, abs=0.01)
    assert not seed_sufficiency(5)["adequate_for_budget"]
    assert seed_sufficiency(20)["sd_relative_error"] == pytest.approx(0.162, abs=0.01)
    assert seed_sufficiency(20)["adequate_for_budget"]


def test_common_random_numbers_are_adopted_only_if_measured():
    """Matched seeds cancel common noise when the two states transport
    similarly. A large material change makes histories diverge, and assuming
    the benefit would understate exactly the effects that matter most."""
    helped = common_random_number_benefit(
        PairedEffect("q", np.zeros(5), np.full(5, 0.2) + 1e-6),
        PairedEffect("q", np.zeros(5), np.array([0.1, 0.3, 0.15, 0.25, 0.2])))
    assert helped["adopt_common_random_numbers"]

    did_not = common_random_number_benefit(
        PairedEffect("q", np.zeros(5), np.array([0.1, 0.3, 0.15, 0.25, 0.2])),
        PairedEffect("q", np.zeros(5), np.full(5, 0.2) + 1e-6))
    assert not did_not["adopt_common_random_numbers"]


def test_mass_weighted_mean_dose_is_mass_weighted():
    dose = np.array([1.0, 3.0])
    mass = np.array([1.0, 3.0])
    # (1*1 + 3*3)/4 = 2.5, not the unweighted 2.0
    assert mass_weighted_mean_dose(dose, mass) == pytest.approx(2.5)


# --- the uncertainty ledger ---------------------------------------------

def test_the_distribution_ledger_never_samples_an_engineering_choice():
    """The single most important rule in the ledger. Giving a mesh factor or a
    mass-denominator method a probability distribution mixes an engineering
    decision into the epistemic answer, and the resulting probability cannot be
    interpreted by anyone."""
    import csv
    from pathlib import Path

    from physical_contract import distribution_row_problems

    ledger = (Path(__file__).resolve().parents[2] / "data" / "uncertainty"
              / "feedback_parameter_distributions.csv")
    rows = [r for r in csv.DictReader(
        l for l in ledger.read_text(encoding="utf-8").splitlines(keepends=True)
        if not l.startswith("#"))]
    assert rows
    problems = [(r["parameter_id"], p) for r in rows
                for p in distribution_row_problems(r)]
    assert not problems, problems

    by_id = {r["parameter_id"]: r for r in rows}
    # Numerical controls are bounded by convergence, never drawn.
    for pid in ("mesh_coarsening_factor", "transport_histories",
                "material_volume_samples"):
        assert by_id[pid]["sampling_role"] == "convergence_axis"
    # The seed is a nested replicate, not a prior over a physical quantity.
    assert by_id["transport_seed"]["sampling_role"] == "inner_transport_replicate"
    # And nothing the model cannot represent is sampled at all.
    for pid in ("growth_survival_response", "membrane_response"):
        assert by_id[pid]["uncertainty_type"] == "unsupported"
        assert by_id[pid]["sampling_role"] == "excluded"


def test_correlated_inputs_are_grouped_rather_than_declared_independent():
    """propagate() samples independently and says so. Anything sharing a
    correlation_group must be combined before it reaches the sampler — wet and
    dry mass of one coupon share a balance, and a closed composition is
    negatively correlated by construction."""
    import csv
    from pathlib import Path

    ledger = (Path(__file__).resolve().parents[2] / "data" / "uncertainty"
              / "feedback_parameter_distributions.csv")
    rows = [r for r in csv.DictReader(
        l for l in ledger.read_text(encoding="utf-8").splitlines(keepends=True)
        if not l.startswith("#"))]
    groups = {}
    for r in rows:
        g = (r.get("correlation_group") or "").strip()
        if g:
            groups.setdefault(g, []).append(r["parameter_id"])
    assert "gravimetry_coupon" in groups and len(groups["gravimetry_coupon"]) > 1
    # The melanin identifiability problem is recorded as a shared group, since
    # only the PRODUCT of production and coupling reaches the dynamics.
    assert sorted(groups["melanin_product"]) == ["melanin_coupling",
                                                 "melanin_response"]


def test_no_gate_verdict_has_been_recorded():
    """The artifact ships empty on purpose: no measured config exists and no
    threshold has been declared, so any row here would be a decision nobody
    made."""
    import csv
    from pathlib import Path

    path = (Path(__file__).resolve().parents[2] / "data" / "uncertainty"
            / "feedback_gate_verdicts.csv")
    rows = [r for r in csv.DictReader(
        l for l in path.read_text(encoding="utf-8").splitlines(keepends=True)
        if not l.startswith("#"))]
    assert rows == []
