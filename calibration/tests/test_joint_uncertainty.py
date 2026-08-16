"""The two Phase 2B gate conditions the fixture suite cannot check on its own:
correlated draws must stay correlated, and closed compositions must stay closed.

Both failures are silent. An uncorrelated wet/dry pair still produces a water
fraction, and an open composition still produces a number — they just describe a
measurement nobody made and a material that cannot exist.
"""

from __future__ import annotations

import numpy as np
import pytest

from biofilm_calibration.joint_uncertainty import (CorrelationGroup, JointSample,
                                                   closure_error,
                                                   observed_correlation,
                                                   sample_composition,
                                                   sample_joint)

WATER = {"H": 0.111894, "O": 0.888106}
DRAWS = 20000


def test_a_shared_balance_calibration_keeps_wet_and_dry_correlated():
    """The gate condition. Wet and dry mass of ONE coupon are weighed on one
    balance; sampling them independently would let the water fraction absorb an
    error that physically cancels."""
    rng = np.random.default_rng(1)
    shared = CorrelationGroup(
        name="gravimetry_coupon",
        values={"wet_mass_g": 1.400, "dry_mass_g": 0.180},
        shared_relative_sigma=0.02,             # the balance
        independent_relative_sigma={"wet_mass_g": 0.002, "dry_mass_g": 0.002})
    s = shared.sample(rng, DRAWS)
    r = observed_correlation(s["wet_mass_g"], s["dry_mass_g"])
    assert r > 0.95, f"shared calibration did not survive sampling (r={r:.3f})"

    # The physical consequence: the water fraction is far TIGHTER than it would
    # be under independent sampling, because the shared error cancels in a ratio.
    w_shared = 1.0 - s["dry_mass_g"] / s["wet_mass_g"]
    independent = CorrelationGroup(
        name="wrong", values={"wet_mass_g": 1.400, "dry_mass_g": 0.180},
        shared_relative_sigma=0.0,
        independent_relative_sigma={"wet_mass_g": 0.02, "dry_mass_g": 0.02})
    si = independent.sample(np.random.default_rng(1), DRAWS)
    w_independent = 1.0 - si["dry_mass_g"] / si["wet_mass_g"]
    assert w_shared.std() < 0.5 * w_independent.std(), (
        "treating a shared balance error as independent inflates the water "
        "fraction's uncertainty, which is the failure this test exists for")


def test_a_group_with_no_shared_error_is_not_spuriously_correlated():
    """The converse, so the test above is not passing for a trivial reason."""
    rng = np.random.default_rng(2)
    g = CorrelationGroup(name="unrelated", values={"a": 1.0, "b": 1.0},
                         shared_relative_sigma=0.0,
                         independent_relative_sigma={"a": 0.05, "b": 0.05})
    s = g.sample(rng, DRAWS)
    assert abs(observed_correlation(s["a"], s["b"])) < 0.05


def test_every_composition_draw_stays_on_the_simplex():
    """The gate condition. Not 'most draws' and not 'after renormalisation' —
    every draw, by construction, because repair is what destroys the induced
    correlations."""
    rng = np.random.default_rng(3)
    comp = sample_composition(WATER, rng, DRAWS, log_sigma=0.10)
    err = closure_error(comp)
    assert err.max() < 1e-12, f"largest closure error {err.max():.3e}"
    for element, values in comp.items():
        assert (values >= 0).all(), f"{element} went negative"
        assert (values <= 1).all(), f"{element} exceeded 1"


def test_composition_components_are_negatively_correlated_by_construction():
    """They sum to one, so raising one must lower the others. A sampler giving
    independent fractions would describe a different material each draw."""
    rng = np.random.default_rng(4)
    comp = sample_composition({"C": 0.5, "H": 0.3, "O": 0.2}, rng, DRAWS)
    r = observed_correlation(comp["C"], comp["O"])
    assert r < -0.2, f"expected negative correlation on the simplex, got {r:.3f}"


def test_a_declared_zero_component_stays_zero():
    """A component declared absent is a statement, not a small number to jitter
    around: perturbing it would invent material that was measured to be absent."""
    rng = np.random.default_rng(5)
    comp = sample_composition({"C": 0.6, "H": 0.4, "N": 0.0}, rng, 500)
    assert np.all(comp["N"] == 0.0)
    assert closure_error(comp).max() < 1e-12


def test_a_wide_perturbation_still_closes():
    rng = np.random.default_rng(6)
    for sigma in (0.01, 0.25, 1.0):
        comp = sample_composition(WATER, rng, 2000, log_sigma=sigma)
        assert closure_error(comp).max() < 1e-12, f"sigma={sigma}"


def test_an_invalid_central_composition_is_refused():
    rng = np.random.default_rng(7)
    with pytest.raises(ValueError, match="negative"):
        sample_composition({"C": -0.1, "H": 1.1}, rng, 10)
    with pytest.raises(ValueError, match="sums to zero"):
        sample_composition({"C": 0.0}, rng, 10)


def test_one_seed_reproduces_the_whole_design():
    """An interval that changes between runs cannot be reviewed or cited, and
    the whole outer design must be recoverable from one recorded integer."""
    groups = [CorrelationGroup("gravimetry_coupon",
                               {"wet_mass_g": 1.4, "dry_mass_g": 0.18},
                               shared_relative_sigma=0.02)]
    comps = {"baseline_biomass": WATER}
    a = sample_joint(groups, comps, seed=42, draws=8)
    b = sample_joint(groups, comps, seed=42, draws=8)
    c = sample_joint(groups, comps, seed=43, draws=8)
    assert isinstance(a, JointSample) and a.seed == 42 and a.draws == 8
    for key, values in a.flat().items():
        assert np.array_equal(values, b.flat()[key]), key
    assert not np.array_equal(a.flat()["wet_mass_g"], c.flat()["wet_mass_g"])


def test_the_flat_view_names_composition_elements_unambiguously():
    groups = [CorrelationGroup("g", {"rho": 1.05}, shared_relative_sigma=0.01)]
    s = sample_joint(groups, {"baseline_biomass": WATER}, seed=1, draws=4)
    flat = s.flat()
    assert "rho" in flat
    assert "baseline_biomass.H" in flat and "baseline_biomass.O" in flat
    stacked = np.vstack([flat["baseline_biomass.H"], flat["baseline_biomass.O"]])
    assert np.allclose(stacked.sum(axis=0), 1.0)
