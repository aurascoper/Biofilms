"""Basis-tagged quantities, and conversions that refuse the wrong basis.

The defect this exists to prevent is live in the model right now:
`compute_radial_biomass` returns dimensionless site occupancy fractions, and
their unweighted mean is assigned to `RadiolysisParams.X_total`, a field
declared `g cm^-3`, then multiplied by a rate in `cm^3 g^-1 s^-1`
(biofilms_potts.jl L1341-1342, L1445-1446, L1797, L1241). The product is `s^-1`
only if X is a mass concentration.

A RANGE CHECK CANNOT CATCH THIS. A biofilm dry-biomass concentration of about
0.05 g/cm3 and a site occupancy of 0.05 are both entirely plausible numbers.
Only a label distinguishes them, so every quantity here carries one and every
converter checks it.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..acquisition import SEGMENTATION_BASIS, WHOLE_BIOFILM_BASIS

# What a number is measured PER. Not a unit — a reference frame.
BASES = frozenset({
    "site_occupancy",       # dimensionless lattice occupancy; never a concentration
    "wet_bulk",             # per wet mass, or per hydrated volume
    "dry_biomass",          # per dry mass
    "ash",                  # per ash mass
    "medium",               # per medium mass
    "membrane_dry",
    "membrane_hydrated",
})

CONCENTRATION_UNIT = "g_cm3"
DIMENSIONLESS = "dimensionless"


class BasisError(ValueError):
    """Raised when a quantity is used against a basis it does not have."""


@dataclass(frozen=True)
class Quantity:
    """A number that knows what it is measured in and measured per."""
    value: float
    unit: str
    basis: str

    def __post_init__(self):
        if self.basis not in BASES:
            raise BasisError(
                f"unknown basis {self.basis!r}; expected one of {sorted(BASES)}")

    def require(self, unit: str, basis: str | tuple) -> float:
        """Return the value, or raise naming both what was wanted and what
        arrived."""
        wanted = (basis,) if isinstance(basis, str) else tuple(basis)
        if self.unit != unit or self.basis not in wanted:
            extra = ""
            if self.basis == "site_occupancy":
                extra = (" — site occupancy is a DIMENSIONLESS lattice fraction, "
                         "not a biomass concentration. This is the conversion "
                         "that is wrong in biofilms_potts.jl L1445/L1797; see "
                         "docs/calibration/material_basis_contract.md.")
            raise BasisError(
                f"expected ({unit}, {'|'.join(wanted)}), got "
                f"({self.unit}, {self.basis}).{extra}")
        return self.value


def occupancy(value: float) -> Quantity:
    """Wrap a CPM site-occupancy fraction so it can never be mistaken for a
    concentration downstream."""
    if not 0.0 <= value <= 1.0:
        raise BasisError(f"site occupancy {value} is outside [0, 1]")
    return Quantity(float(value), DIMENSIONLESS, "site_occupancy")


def _positive(name: str, x: float) -> float:
    if x is None:
        raise BasisError(f"{name} is not measured — refusing to derive from it")
    if x <= 0:
        raise BasisError(f"{name} must be positive, got {x}")
    return float(x)


@dataclass(frozen=True)
class HydratedVolume:
    """A volume that knows what it encloses.

    THE MASS AND THE VOLUME MUST DESCRIBE THE SAME MATERIAL. A balance weighs
    the whole biofilm — cells, matrix, interstitial water, retained solutes.
    A cell-only fluorescence segmentation encloses none of the last three, so
    dividing a whole-biofilm mass by a cell-only volume overstates rho_wet by
    whatever fraction of the biofilm is matrix and water. Nothing about the
    two numbers reveals the mismatch: both are positive floats of plausible
    size, and the error is systematic rather than noisy.
    """
    value_cm3: float
    segmentation_basis: str

    def __post_init__(self):
        if self.segmentation_basis not in SEGMENTATION_BASIS:
            raise BasisError(
                f"unknown segmentation basis {self.segmentation_basis!r}; "
                f"expected one of {sorted(SEGMENTATION_BASIS)}")
        if self.value_cm3 is None or self.value_cm3 <= 0:
            raise BasisError(f"hydrated volume must be positive, got {self.value_cm3}")

    def require_whole_biofilm(self) -> float:
        """The volume, or a refusal naming what it fails to enclose."""
        if self.segmentation_basis == WHOLE_BIOFILM_BASIS:
            return float(self.value_cm3)
        if self.segmentation_basis == "cells_only":
            raise BasisError(
                "hydrated volume was segmented as 'cells_only' but the mass is "
                "the whole biofilm: the volume excludes the extracellular "
                "matrix and the interstitial water that the balance weighed, so "
                "rho_wet would be overstated by that fraction. Segment the "
                "biofilm envelope, or supply a pore fraction via "
                "to_whole_biofilm().")
        if self.segmentation_basis == "cells_and_matrix":
            raise BasisError(
                "hydrated volume was segmented as 'cells_and_matrix': a union "
                "of stained voxels still excludes water-filled internal voids. "
                "Convert it with to_whole_biofilm(pore_volume_fraction=...) and "
                "declare the pore fraction, or segment the envelope directly.")
        raise BasisError(
            "hydrated volume has an unresolved segmentation basis — what the "
            "volume encloses must be declared before it can be divided into a "
            "mass")


def to_whole_biofilm(volume: HydratedVolume,
                     pore_volume_fraction: float) -> HydratedVolume:
    """Inflate a stained-voxel volume to the enclosing envelope.

        V_envelope = V_stained / (1 - pore_fraction)

    The pore fraction is a DECLARED measurement, not a default: it is exactly
    the quantity that distinguishes the two volumes, so assuming it would
    reintroduce the error this class exists to catch.
    """
    if volume.segmentation_basis == WHOLE_BIOFILM_BASIS:
        return volume
    if pore_volume_fraction is None:
        raise BasisError(
            "no declared pore volume fraction — it is precisely the quantity "
            "that separates a stained-voxel volume from the biofilm envelope")
    p = float(pore_volume_fraction)
    if not 0.0 <= p < 1.0:
        raise BasisError(f"pore volume fraction {p} is outside [0, 1)")
    return HydratedVolume(float(volume.value_cm3) / (1.0 - p),
                          WHOLE_BIOFILM_BASIS)


def wet_bulk_density(wet_mass_g: float, hydrated_volume: HydratedVolume) -> Quantity:
    """rho_wet = m_wet / V_hydrated, with both sides describing one material."""
    m = _positive("wet_mass_g", wet_mass_g)
    v = _require_volume(hydrated_volume)
    return Quantity(m / v, CONCENTRATION_UNIT, "wet_bulk")


def dry_biomass_concentration(dry_mass_g: float,
                              hydrated_volume: HydratedVolume) -> Quantity:
    """X_total = m_dry / V_hydrated.

    Dry mass over HYDRATED volume. The mixed basis is deliberate — it is what
    the kinetics need — and is exactly what gets silently 'simplified' to dry
    mass over dry volume, which is larger by roughly 1/(1 - w_water).
    """
    m = _positive("dry_mass_g", dry_mass_g)
    v = _require_volume(hydrated_volume)
    return Quantity(m / v, CONCENTRATION_UNIT, "wet_bulk")


def _require_volume(hydrated_volume) -> float:
    """A bare float is refused: an untagged volume is the untagged half of the
    density, and the mass side has been tagged since the basis contract."""
    if not isinstance(hydrated_volume, HydratedVolume):
        raise BasisError(
            f"hydrated volume must be a HydratedVolume carrying its "
            f"segmentation basis, got {type(hydrated_volume).__name__}. A bare "
            "number cannot say whether it encloses the same material the "
            "balance weighed.")
    return hydrated_volume.require_whole_biofilm()


def water_mass_fraction(wet_mass_g: float, dry_mass_g: float) -> Quantity:
    """w_water = (m_wet - m_dry) / m_wet."""
    wet = _positive("wet_mass_g", wet_mass_g)
    dry = _positive("dry_mass_g", dry_mass_g)
    if dry > wet:
        raise BasisError(
            f"dry mass {dry} exceeds wet mass {wet} — check the drying and "
            "surface-water removal protocols")
    return Quantity((wet - dry) / wet, DIMENSIONLESS, "wet_bulk")


def ash_mass_fraction(ash_mass_g: float, dry_mass_g: float) -> Quantity:
    """Ash as a fraction of DRY mass."""
    ash = _positive("ash_mass_g", ash_mass_g)
    dry = _positive("dry_mass_g", dry_mass_g)
    if ash > dry:
        raise BasisError(f"ash mass {ash} exceeds dry mass {dry}")
    return Quantity(ash / dry, DIMENSIONLESS, "dry_biomass")


def metal_reducing_concentration(x_total: Quantity,
                                 active_fraction: Quantity) -> Quantity:
    """X_red = f_red,dry * X_total.

    `active_fraction` must be on a dry-biomass basis and must be the ACTIVE
    REDUCER fraction. Taxonomic abundance of S. oneidensis is an upper bound on
    it at best: cells present but not respiring the metal do not reduce it.
    """
    x = x_total.require(CONCENTRATION_UNIT, "wet_bulk")
    f = active_fraction.require(DIMENSIONLESS, "dry_biomass")
    if not 0.0 <= f <= 1.0:
        raise BasisError(f"active reducer fraction {f} is outside [0, 1]")
    return Quantity(f * x, CONCENTRATION_UNIT, "wet_bulk")


def active_from_taxonomic(taxonomic_fraction: float,
                          activity_fraction: float | None) -> Quantity:
    """Convert a taxonomic dry-mass fraction into an ACTIVE reducer fraction.

    Refuses without a declared activity fraction rather than assuming that
    every cell of a metal-reducing species is reducing metal.
    """
    if activity_fraction is None:
        raise BasisError(
            "no declared activity fraction — taxonomic abundance is not "
            "functional activity, and assuming they are equal would silently "
            "set every present cell to fully active")
    t = float(taxonomic_fraction)
    a = float(activity_fraction)
    for name, v in (("taxonomic_fraction", t), ("activity_fraction", a)):
        if not 0.0 <= v <= 1.0:
            raise BasisError(f"{name} {v} is outside [0, 1]")
    return Quantity(t * a, DIMENSIONLESS, "dry_biomass")


def dry_basis_to_wet_fraction(value_per_g_dry: float,
                              water_fraction: Quantity) -> Quantity:
    """Convert a per-dry-mass fraction (e.g. mg metal per g dry) to a fraction
    of WET mass: w_wet = w_dry * (1 - w_water)."""
    w = water_fraction.require(DIMENSIONLESS, "wet_bulk")
    if not 0.0 <= w < 1.0:
        raise BasisError(f"water mass fraction {w} is outside [0, 1)")
    return Quantity(float(value_per_g_dry) * (1.0 - w), DIMENSIONLESS, "wet_bulk")


def blank_corrected_mass(sample_plus_substrate_g: float,
                         blank_g: float | None, *, label: str) -> float:
    """m_biofilm = m_sample+substrate - m_blank.

    Blank subtraction is the DEFINITION of a biofilm mass, not a refinement,
    so a missing blank is refused rather than treated as zero. A hydrated
    membrane retains substantial water; counting it as biofilm water inflates
    rho_wet and deflates w_water, and nothing downstream would notice.
    """
    m = _positive(f"{label} sample+substrate mass", sample_plus_substrate_g)
    if blank_g is None:
        raise BasisError(
            f"{label}: no blank mass — blank subtraction defines the biofilm "
            "mass, and omitting it silently attributes the substrate (and any "
            "water it retains) to the biofilm")
    b = float(blank_g)
    if b < 0:
        raise BasisError(f"{label}: negative blank mass {b}")
    if b >= m:
        raise BasisError(
            f"{label}: blank mass {b} is not less than the sample mass {m} — "
            "the correction would leave no biofilm")
    return m - b
