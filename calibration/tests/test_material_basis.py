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
    BasisError, HydratedVolume, Quantity, active_from_taxonomic,
    ash_mass_fraction, dry_basis_to_wet_fraction, dry_biomass_concentration,
    metal_reducing_concentration, occupancy, to_whole_biofilm,
    water_mass_fraction, wet_bulk_density)


def _vol(cm3):
    """A volume that encloses what the balance weighed."""
    return HydratedVolume(cm3, "whole_biofilm_envelope")


def test_site_occupancy_cannot_become_a_concentration():
    """THE test. compute_radial_biomass returns a site fraction and its mean is
    assigned to X_total, declared g/cm3 and multiplied by cm3/g/s. Both
    quantities can plausibly read 0.05, so only the label can catch it."""
    occ = occupancy(0.05)
    assert occ.unit == "dimensionless"
    with pytest.raises(BasisError, match="not a biomass concentration"):
        occ.require("g_cm3", "wet_bulk")

    # and it cannot be smuggled in as the active-reducer fraction either
    x_total = dry_biomass_concentration(dry_mass_g=0.5, hydrated_volume=_vol(10.0))
    with pytest.raises(BasisError, match="not a biomass concentration"):
        metal_reducing_concentration(x_total, occ)


def test_occupancy_is_bounded():
    for bad in (-0.01, 1.5):
        with pytest.raises(BasisError, match="outside"):
            occupancy(bad)


def test_wet_and_dry_definitions():
    rho = wet_bulk_density(wet_mass_g=10.0, hydrated_volume=_vol(10.0))
    assert rho.value == pytest.approx(1.0)
    assert rho.unit == "g_cm3"

    x = dry_biomass_concentration(dry_mass_g=0.98, hydrated_volume=_vol(10.0))
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
    rho = wet_bulk_density(wet, _vol(vol)).value
    w = water_mass_fraction(wet, dry).value
    x = dry_biomass_concentration(dry, _vol(vol)).value
    assert x == pytest.approx(rho * (1.0 - w))


def test_x_red_is_the_product_of_active_fraction_and_x_total():
    x_total = dry_biomass_concentration(dry_mass_g=1.0, hydrated_volume=_vol(10.0))
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
        wet_bulk_density(None, _vol(10.0))
    with pytest.raises(BasisError, match="must be a HydratedVolume"):
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


# --- volume basis: the mass and the volume must describe one material -------

def test_cell_only_volume_against_whole_biofilm_mass_is_refused():
    """THE volume-basis test. A balance weighs cells, matrix, interstitial
    water and retained solutes; a cell-only segmentation encloses none of the
    last three. Both numbers are positive floats of plausible size, so nothing
    but the label can catch it, and the error is systematic — rho_wet comes out
    high by whatever fraction of the biofilm is matrix and water."""
    cells_only = HydratedVolume(10.0, "cells_only")
    with pytest.raises(BasisError, match="excludes the extracellular"):
        wet_bulk_density(wet_mass_g=10.0, hydrated_volume=cells_only)
    with pytest.raises(BasisError, match="excludes the extracellular"):
        dry_biomass_concentration(dry_mass_g=1.0, hydrated_volume=cells_only)


def test_cells_and_matrix_still_needs_a_declared_pore_fraction():
    """A union of stained voxels excludes water-filled internal voids, so it is
    closer but still not the envelope."""
    v = HydratedVolume(10.0, "cells_and_matrix")
    with pytest.raises(BasisError, match="water-filled internal voids"):
        wet_bulk_density(10.0, v)

    with pytest.raises(BasisError, match="no declared pore volume fraction"):
        to_whole_biofilm(v, None)

    envelope = to_whole_biofilm(v, 0.2)
    assert envelope.segmentation_basis == "whole_biofilm_envelope"
    assert envelope.value_cm3 == pytest.approx(12.5)          # 10 / (1 - 0.2)
    assert wet_bulk_density(10.0, envelope).value == pytest.approx(0.8)


def test_the_correction_moves_rho_wet_in_the_right_direction():
    """Ignoring the pore volume inflates the density; correcting it must
    reduce it."""
    stained = HydratedVolume(10.0, "cells_and_matrix")
    naive = 10.0 / stained.value_cm3                     # what the bug produces
    corrected = wet_bulk_density(10.0, to_whole_biofilm(stained, 0.3)).value
    assert corrected < naive
    assert corrected == pytest.approx(naive * 0.7)


def test_an_unresolved_basis_is_refused():
    with pytest.raises(BasisError, match="must be declared"):
        wet_bulk_density(10.0, HydratedVolume(10.0, "unresolved"))


def test_a_bare_number_is_no_longer_accepted_as_a_volume():
    """The mass side has carried a basis since the unit contract; the volume
    was the untagged half."""
    with pytest.raises(BasisError, match="must be a HydratedVolume"):
        wet_bulk_density(10.0, 10.0)


def test_unknown_segmentation_basis_is_rejected_at_construction():
    with pytest.raises(BasisError, match="unknown segmentation basis"):
        HydratedVolume(1.0, "eyeballed")
    with pytest.raises(BasisError, match="must be positive"):
        HydratedVolume(0.0, "whole_biofilm_envelope")


def test_whole_biofilm_volume_passes_through_unchanged():
    v = HydratedVolume(4.0, "whole_biofilm_envelope")
    assert to_whole_biofilm(v, None) is v
    assert v.require_whole_biofilm() == 4.0
