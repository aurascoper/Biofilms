"""Phase 2A: ten synthetic cases whose correct verdict is declared in advance.

The success criterion for this phase is not that the gate passes something. It
is that the gate REFUSES, DETECTS AND PASSES EXACTLY the cases it should — so
each fixture below states its expected verdict before the machinery runs, and a
gate that returned "pass" for all ten would fail this file completely.

No OpenMC here, on purpose. Gate logic must be shown correct independently of
transport, or a wrong verdict later is ambiguous between the decision rule and
the physics that fed it.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from biofilm_openmc.fingerprint import dose_state_hash, needs_rerun, transport_state_hash
from biofilm_openmc.lineage import aggregate_by_label
from biofilm_openmc.synthetic_gate import (
    EFFECT_BELOW_THRESHOLD, EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT,
    FORBIDDEN_IN_S0, INDETERMINATE_CALIBRATION, INDETERMINATE_MODEL_FORM,
    INDETERMINATE_NUMERICS, INDETERMINATE_TRANSPORT, PASS_SYNTHETIC_GATE,
    S0_VERDICTS, S0Verdict, ThresholdPolicy, VarianceBudget, decide,
    scenario_envelope)

POLICY = ThresholdPolicy()          # synthetic_gate_contract_v1, tier S0
RNG = np.random.default_rng(20260816)

ARTIFACT = Path(__file__).resolve().parents[2] / "artifacts" / "synthetic_gate_verdicts.json"


def _draws(mean, sd, n=40000):
    return RNG.normal(mean, sd, n)


# Each entry: (fixture id, effect draws, variance budget, expected verdict).
# The expected verdict is the CONTRACT; changing one to match an observed
# outcome would defeat the purpose of the file.
def _cases():
    return [
        # 1. baseline and feedback identical -> no effect at all
        ("identical_states", np.zeros(40000),
         VarianceBudget(transport=0.0), EFFECT_BELOW_THRESHOLD),

        # 2. a large, well-resolved monotone material effect
        ("large_monotone_effect", _draws(0.40, 0.01),
         VarianceBudget(transport=1e-4), PASS_SYNTHETIC_GATE),

        # 3. real effect drowned in REPLICATE noise -> more histories would help
        ("hidden_by_transport_noise", _draws(0.12, 0.09),
         VarianceBudget(transport=8.1e-3, numerics=1e-5, calibration=1e-5),
         INDETERMINATE_TRANSPORT),

        # 4. same effect smeared by a coarse mesh. A FLOOR: histories cannot
        #    remove resolution loss, which is the A0 lesson in one fixture.
        ("hidden_by_coarse_mesh", _draws(0.12, 0.09),
         VarianceBudget(transport=1e-5, numerics=8.1e-3, calibration=1e-5),
         INDETERMINATE_NUMERICS),

        # 5. same effect, but the spread comes from calibration draws -> no
        #    amount of compute helps; only better measurements do
        ("hidden_by_calibration_draws", _draws(0.12, 0.09),
         VarianceBudget(transport=1e-5, numerics=1e-5, calibration=8.1e-3),
         INDETERMINATE_CALIBRATION),

        # 6. statistically resolved, practically irrelevant
        ("resolved_but_small", _draws(0.05, 0.002),
         VarianceBudget(transport=4e-6),
         EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT),
    ]


@pytest.mark.parametrize("fixture_id,draws,budget,expected", _cases(),
                         ids=[c[0] for c in _cases()])
def test_each_fixture_returns_its_declared_verdict(fixture_id, draws, budget, expected):
    got = decide(draws, budget, POLICY)
    assert got.verdict == expected, f"{fixture_id}: {got.verdict} != {expected}\n{got.reason}"
    assert got.target_calibration is False
    assert got.tier == "S0"
    assert got.threshold_policy_id == "synthetic_gate_contract_v1"


def test_disagreeing_model_forms_are_not_averaged():
    """Fixture 7. Two admissible scenarios on opposite sides of the threshold.

    They are not two samples from a distribution over scenarios. Assigning them
    equal probability would manufacture a number that looks like evidence, so
    the gate reports the envelope and refuses."""
    scenarios = {"m_mech_law": _draws(0.04, 0.005),
                 "P_eff_law": _draws(0.30, 0.005)}
    pooled = np.concatenate(list(scenarios.values()))
    got = decide(pooled, VarianceBudget(model_form=1.7e-2, transport=1e-5),
                 POLICY, scenario_effects=scenarios)
    assert got.verdict == INDETERMINATE_MODEL_FORM
    assert got.dominant_variance == "model_form"

    env = scenario_envelope(scenarios)
    assert env["envelope_low"] < POLICY.effect_threshold <= env["envelope_high"]
    # The pooled mean sits between them and would have looked like an answer.
    assert 0.04 < float(np.mean(pooled)) < 0.30


def test_a_label_only_change_reuses_transport(snapshot, config):
    """Fixture 8. Labels never become materials, so relabelling cannot change
    the transport identity — only the attribution that reads it."""
    th = transport_state_hash(snapshot, config, "endfb-viii.0")
    relabelled = type(snapshot)(**{**snapshot.__dict__,
                                   "lineage_id": snapshot.lineage_id + 100,
                                   "label_state_hash": "different"})
    assert not needs_rerun(th, relabelled, config, "endfb-viii.0")
    assert transport_state_hash(relabelled, config, "endfb-viii.0") == th


def test_a_source_rate_only_change_reuses_transport(snapshot, config):
    """Fixture 9. Heating is per source particle, so the activity is not part
    of the transport identity — it rescales the result instead of redoing it."""
    from dataclasses import replace

    th = transport_state_hash(snapshot, config, "endfb-viii.0")
    faster = replace(config, photons_per_second=(config.photons_per_second or 1.0) * 7.0)
    assert transport_state_hash(snapshot, faster, "endfb-viii.0") == th
    # ...but the DOSE identity must move, or a cached Gy/s field would be reused
    # across activities.
    assert dose_state_hash(th, 1.0) != dose_state_hash(th, 7.0)


def test_founder_generation_zero_survives_attribution(snapshot):
    """Fixture 10. Generation 0 is a founder and collides with background."""
    generation = np.zeros_like(snapshot.cell_id)
    occupied = snapshot.cell_id > 0
    dose = np.full(snapshot.cell_id.shape, 3.0)
    mass = np.ones_like(dose)

    without = aggregate_by_label(dose, generation, mass)
    with_mask = aggregate_by_label(dose, generation, mass, occupied=occupied)
    assert without == {}                       # every founder dropped
    assert set(with_mask) == {0}
    assert with_mask[0]["n_voxels"] == int(occupied.sum())


# --- the boundary this tier exists to hold ------------------------------

def test_no_fixture_can_emit_an_offline_feedback_pass():
    """Structural, not conventional: the constant is not in the vocabulary, so
    no input and no policy can produce it."""
    assert not (S0_VERDICTS & FORBIDDEN_IN_S0)
    for _, draws, budget, _ in _cases():
        assert decide(draws, budget, POLICY).verdict not in FORBIDDEN_IN_S0
    # even an enormous, perfectly resolved effect
    huge = decide(_draws(1e6, 1.0), VarianceBudget(transport=1e-9), POLICY)
    assert huge.verdict == PASS_SYNTHETIC_GATE
    assert huge.verdict not in FORBIDDEN_IN_S0


def test_a_policy_claiming_target_calibration_is_refused():
    bad = ThresholdPolicy(target_calibration=True)
    got = decide(_draws(0.4, 0.01), VarianceBudget(transport=1e-4), bad)
    assert got.verdict == "NOT_EVALUATED"
    assert "target_calibration" in got.reason


def test_a_threshold_chosen_after_the_results_is_refused():
    late = ThresholdPolicy(declared_before_results=False)
    assert decide(_draws(0.4, 0.01), VarianceBudget(), late).verdict == "NOT_EVALUATED"


# --- the machine-readable artifact --------------------------------------

def test_verdict_json_is_nan_free_and_schema_valid():
    """A non-finite number is not JSON. A verdict only some parsers can read is
    not machine-readable, so the writer refuses rather than emitting bare NaN."""
    records = {}
    for fixture_id, draws, budget, expected in _cases():
        v = decide(draws, budget, POLICY)
        assert v.verdict == expected
        records[fixture_id] = {"expected": expected, **v.as_dict()}

    scenarios = {"a": _draws(0.04, 0.005), "b": _draws(0.30, 0.005)}
    mf = decide(np.concatenate(list(scenarios.values())),
                VarianceBudget(model_form=1.7e-2), POLICY,
                scenario_effects=scenarios)
    records["disagreeing_model_forms"] = {"expected": INDETERMINATE_MODEL_FORM,
                                          **mf.as_dict()}

    payload = {
        "schema_version": 1,
        "threshold_policy_id": POLICY.threshold_policy_id,
        "tier": POLICY.tier,
        "target_calibration": POLICY.target_calibration,
        "claim_class": "gate_logic_validation",
        "biological_calibration": "NOT_EVALUATED",
        "reference_d_verdict": "NOT_EVALUATED",
        "fixtures": records,
    }
    # allow_nan=False is the check: it raises rather than writing bare NaN.
    text = json.dumps(payload, indent=2, allow_nan=False)
    round_tripped = json.loads(text)

    required = {"verdict", "reason", "effect_median", "effect_low",
                "effect_high", "dominant_variance", "variance", "n_draws",
                "threshold_policy_id", "tier", "target_calibration"}
    for fixture_id, rec in round_tripped["fixtures"].items():
        assert required <= set(rec), f"{fixture_id} is missing {required - set(rec)}"
        assert rec["verdict"] == rec["expected"], fixture_id
        assert rec["target_calibration"] is False, fixture_id
        assert rec["tier"] == "S0"
        for key in ("effect_median", "effect_low", "effect_high"):
            assert rec[key] is None or math.isfinite(rec[key])

    assert round_tripped["reference_d_verdict"] == "NOT_EVALUATED"
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(text + "\n", encoding="utf-8")


def test_every_fixture_verdict_is_in_the_declared_vocabulary():
    for _, draws, budget, _ in _cases():
        assert decide(draws, budget, POLICY).verdict in S0_VERDICTS


def test_the_indeterminate_verdicts_are_distinguished_by_source_alone():
    """The same effect distribution, three different dominant variances, three
    different verdicts. If these collapsed to one, a compute budget could be
    spent on a calibration problem."""
    draws = _draws(0.12, 0.09)
    got = {
        decide(draws, VarianceBudget(transport=8.1e-3), POLICY).verdict,
        decide(draws, VarianceBudget(numerics=8.1e-3), POLICY).verdict,
        decide(draws, VarianceBudget(calibration=8.1e-3), POLICY).verdict,
    }
    assert got == {INDETERMINATE_TRANSPORT, INDETERMINATE_NUMERICS,
                   INDETERMINATE_CALIBRATION}


def test_the_numerics_verdict_says_histories_will_not_help():
    v = decide(_draws(0.12, 0.09), VarianceBudget(numerics=8.1e-3), POLICY)
    assert v.verdict == INDETERMINATE_NUMERICS
    assert "more histories cannot remove" in v.reason
