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
