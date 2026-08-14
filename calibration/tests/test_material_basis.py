"""Basis tagging, conversions, and the occupancy defect it exists to stop.

The load-bearing test is `test_site_occupancy_cannot_become_a_concentration`:
that substitution is live in biofilms_potts.jl right now, and no range check
can catch it because both quantities can read 0.05.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from biofilm_calibration.materials.basis_conversion import (
    BasisError, Quantity, active_from_taxonomic, ash_mass_fraction,
    dry_basis_to_wet_fraction, dry_biomass_concentration,
    metal_reducing_concentration, occupancy, water_mass_fraction,
    wet_bulk_density)


def test_site_occupancy_cannot_become_a_concentration():
    """THE test. compute_radial_biomass returns a site fraction and its mean is
    assigned to X_total, declared g/cm3 and multiplied by cm3/g/s. Both
    quantities can plausibly read 0.05, so only the label can catch it."""
    occ = occupancy(0.05)
    assert occ.unit == "dimensionless"
    with pytest.raises(BasisError, match="not a biomass concentration"):
        occ.require("g_cm3", "wet_bulk")

    # and it cannot be smuggled in as the active-reducer fraction either
    x_total = dry_biomass_concentration(dry_mass_g=0.5, hydrated_volume_cm3=10.0)
    with pytest.raises(BasisError, match="not a biomass concentration"):
        metal_reducing_concentration(x_total, occ)


def test_occupancy_is_bounded():
    for bad in (-0.01, 1.5):
        with pytest.raises(BasisError, match="outside"):
            occupancy(bad)


def test_wet_and_dry_definitions():
    rho = wet_bulk_density(wet_mass_g=10.0, hydrated_volume_cm3=10.0)
    assert rho.value == pytest.approx(1.0)
    assert rho.unit == "g_cm3"

    x = dry_biomass_concentration(dry_mass_g=0.98, hydrated_volume_cm3=10.0)
    assert x.value == pytest.approx(0.098)

    w = water_mass_fraction(wet_mass_g=10.0, dry_mass_g=0.98)
    assert w.value == pytest.approx(0.902)


def test_dry_mass_cannot_exceed_wet_mass():
    """Usually a surface-water or drying-protocol error, so the message says so."""
    with pytest.raises(BasisError, match="surface-water removal"):
        water_mass_fraction(wet_mass_g=1.0, dry_mass_g=2.0)


def test_wet_dry_basis_round_trip():
    """X_total = rho_wet * (1 - w_water), so the three definitions must agree."""
    wet, dry, vol = 12.5, 1.1, 9.0
    rho = wet_bulk_density(wet, vol).value
    w = water_mass_fraction(wet, dry).value
    x = dry_biomass_concentration(dry, vol).value
    assert x == pytest.approx(rho * (1.0 - w))


def test_x_red_is_the_product_of_active_fraction_and_x_total():
    x_total = dry_biomass_concentration(dry_mass_g=1.0, hydrated_volume_cm3=10.0)
    f = Quantity(0.25, "dimensionless", "dry_biomass")
    x_red = metal_reducing_concentration(x_total, f)
    assert x_red.value == pytest.approx(0.025)
    assert x_red.unit == "g_cm3"


def test_taxonomic_abundance_is_not_functional_activity():
    """Cells present but not respiring the metal do not reduce it, so the
    conversion refuses without a declared activity fraction."""
    with pytest.raises(BasisError, match="not functional activity"):
        active_from_taxonomic(0.4, None)
    got = active_from_taxonomic(0.4, 0.5)
    assert got.value == pytest.approx(0.2)
    assert got.basis == "dry_biomass"


def test_unmeasured_input_is_refused_not_defaulted():
    """None means not measured. Substituting zero would produce a confident
    wrong answer."""
    with pytest.raises(BasisError, match="not measured"):
        wet_bulk_density(None, 10.0)
    with pytest.raises(BasisError, match="not measured"):
        dry_biomass_concentration(1.0, None)


def test_dry_basis_to_wet_fraction_uses_the_water_fraction():
    w = water_mass_fraction(wet_mass_g=10.0, dry_mass_g=1.0)   # 0.9
    got = dry_basis_to_wet_fraction(0.02, w)                   # 2% of dry mass
    assert got.value == pytest.approx(0.002)
    assert got.basis == "wet_bulk"


def test_wrong_basis_is_named_in_the_error():
    q = Quantity(0.3, "dimensionless", "ash")
    with pytest.raises(BasisError, match=r"expected \(dimensionless, dry_biomass\)"):
        q.require("dimensionless", "dry_biomass")


def test_ash_fraction_is_on_a_dry_basis():
    a = ash_mass_fraction(ash_mass_g=0.1, dry_mass_g=1.0)
    assert a.value == pytest.approx(0.1)
    assert a.basis == "dry_biomass"
    with pytest.raises(BasisError, match="exceeds dry mass"):
        ash_mass_fraction(2.0, 1.0)


def test_unknown_basis_is_rejected_at_construction():
    with pytest.raises(BasisError, match="unknown basis"):
        Quantity(1.0, "g_cm3", "vibes")
