"""Where the two calibration branches meet.

The spatial branch says what one occupied CPM voxel *represents*. The material
branch says what material a voxel is *made of*. Neither is meaningful without
the other, and the combination can be incoherent in ways that neither branch
can detect alone — which is why this check lives here rather than in either of
them.

The contract is one sentence:

    One occupied CPM voxel of physical volume a^3 represents ______
    and is assigned OpenMC material ______ at density ______.

Coherent:

    part of a hydrated computational biomass parcel,
    assigned the hydrated-effective-biofilm material,
    at a measured wet bulk density.

Incoherent, and the reason this module exists:

    a literal cell fragment, whose density comes from bulk wet biofilm,
    whose lineage is claimed as a biological cell lineage.

Every clause there is individually defensible. Together they are not: a cell
interior is not bulk biofilm, and a lineage of parcels is not a lineage of
organisms. Each branch would pass its own gate.
"""

from __future__ import annotations

from dataclasses import dataclass

# What the occupied voxel is claimed to BE. From the spatial branch.
ENTITY_KINDS = frozenset({"biomass_parcel", "biological_cell",
                          "hyphal_compartment", "conidium", "unresolved"})
LINEAGE_SEMANTICS = frozenset({"computational", "compartmental", "biological",
                               "unresolved"})

# What the assigned material is claimed to REPRESENT. From the material branch.
MATERIAL_REPRESENTS = frozenset({
    "hydrated_bulk_biofilm",   # cells + water + EPS blended: a parcel's contents
    "cytoplasm",               # the inside of one cell
    "dry_solid",               # biomass with the water removed
    "unresolved",
})
MATERIAL_MODEL_KINDS = frozenset({"hydrated_effective_medium",
                                  "explicit_components"})

# Entity kinds that name an individual organism or organism part. These are the
# ones whose contents are NOT bulk biofilm.
_ORGANISM_SCALE = frozenset({"biological_cell", "hyphal_compartment", "conidium"})


@dataclass(frozen=True)
class VoxelBinding:
    """The completed integration sentence, as data."""
    entity_kind: str
    lineage_semantics: str
    material_class: str
    material_represents: str
    material_model_kind: str
    density_basis: str
    lattice_pitch_um: float | None = None
    density_g_cm3: float | None = None
    notes: str = ""

    def sentence(self) -> str:
        a = f"{self.lattice_pitch_um} um" if self.lattice_pitch_um else "a (UNSET)"
        rho = f"{self.density_g_cm3} g/cm3" if self.density_g_cm3 else "density (UNSET)"
        return (f"One occupied CPM voxel of edge {a} represents part of a "
                f"{self.entity_kind}, and is assigned OpenMC material "
                f"'{self.material_class}' ({self.material_represents}) at {rho}.")


def coherence_problems(b: VoxelBinding) -> list[str]:
    """Every way this binding contradicts itself. Empty means coherent."""
    problems: list[str] = []

    for field, allowed in (("entity_kind", ENTITY_KINDS),
                           ("lineage_semantics", LINEAGE_SEMANTICS),
                           ("material_represents", MATERIAL_REPRESENTS),
                           ("material_model_kind", MATERIAL_MODEL_KINDS)):
        value = getattr(b, field)
        if value not in allowed:
            problems.append(f"{field}={value!r} is not one of {sorted(allowed)}")

    # THE mixed-semantics rejection: an organism-scale entity whose material is
    # bulk biofilm. Both halves are defensible; together they are not.
    if b.entity_kind in _ORGANISM_SCALE and \
            b.material_represents == "hydrated_bulk_biofilm":
        problems.append(
            f"entity_kind={b.entity_kind!r} names an individual organism or "
            "organism part, but the material represents bulk biofilm (cells, "
            "water and EPS blended). A cell interior is not bulk biofilm — "
            "pick one scale")

    # ... and its mirror image.
    if b.entity_kind == "biomass_parcel" and b.material_represents == "cytoplasm":
        problems.append(
            "entity_kind='biomass_parcel' bundles cells with their water and "
            "extracellular material, so its contents cannot be cytoplasm")

    # Lineage must match the entity it is a lineage OF.
    if b.lineage_semantics == "biological" and b.entity_kind == "biomass_parcel":
        problems.append(
            "lineage_semantics='biological' with entity_kind='biomass_parcel': "
            "parcel ancestry is not organism ancestry. The CPM has no growth "
            "dynamics from which a biological generation could be claimed "
            "(biofilms_potts.jl:1578-1583)")
    if b.lineage_semantics == "computational" and b.entity_kind in _ORGANISM_SCALE:
        problems.append(
            f"lineage_semantics='computational' with entity_kind="
            f"{b.entity_kind!r}: if the entity is an organism, its lineage is "
            "not merely computational — say which one is intended")

    # A dry solid in a voxel that also has to hold the water.
    if b.material_represents == "dry_solid" and \
            b.material_model_kind == "hydrated_effective_medium":
        problems.append(
            "material_represents='dry_solid' contradicts "
            "material_model_kind='hydrated_effective_medium': the effective "
            "medium includes the water")

    # Mapper compatibility, restated here because it is a property of the
    # BINDING, not of the material alone.
    if b.material_model_kind == "explicit_components":
        problems.append(
            "material_model_kind='explicit_components' cannot be bound to a "
            "voxel: the mapper assigns one material class per occupied voxel "
            "(coupling/biofilm_openmc/materials.py), so separate solid and "
            "water phases would need sub-voxel homogenization that does not "
            "exist")

    if b.density_basis == "site_occupancy":
        problems.append(
            "density_basis='site_occupancy' is not a density at all — it is a "
            "dimensionless lattice fraction. See "
            "docs/calibration/material_basis_contract.md")
    elif b.density_basis not in {"wet_bulk", "unresolved"}:
        problems.append(
            f"density_basis={b.density_basis!r}: a voxel-averaged material "
            "density is a wet bulk density")

    return problems


def is_coherent(b: VoxelBinding) -> bool:
    return not coherence_problems(b)


@dataclass(frozen=True)
class IntegrationVerdict:
    can_emit_config: bool
    binding_coherent: bool
    blockers: tuple
    reasons: tuple

    def render(self) -> str:
        head = ("READY to emit a provisional biofilm transport config"
                if self.can_emit_config else
                "REFUSING to emit a biofilm transport config")
        lines = [head,
                 f"  binding coherent: {'yes' if self.binding_coherent else 'NO'}"]
        for r in self.reasons:
            lines.append(f"  + {r}")
        for b in self.blockers:
            lines.append(f"  - {b}")
        return "\n".join(lines)


def evaluate(binding: VoxelBinding, spatial_verdict: str,
             material_openmc_verdict: str) -> IntegrationVerdict:
    """May a provisional biofilm transport config be emitted?

    Three independent conditions, all required. The binding being coherent is
    NOT sufficient: a coherent sentence with unmeasured blanks in it is still
    a sentence with blanks.
    """
    problems = coherence_problems(binding)
    blockers = list(problems)
    reasons: list[str] = []

    if not problems:
        reasons.append("voxel binding is internally coherent")

    if spatial_verdict != "READY_FOR_TIME_CALIBRATION":
        blockers.append(
            f"spatial gate is {spatial_verdict}, not READY_FOR_TIME_CALIBRATION "
            "— without a selected lattice pitch there is no voxel volume, so "
            "no mass per voxel")
    else:
        reasons.append("spatial gate ready")

    if material_openmc_verdict != "READY_FOR_OPENMC_BIOFILM_TRANSPORT":
        blockers.append(
            f"material OpenMC gate is {material_openmc_verdict}, not "
            "READY_FOR_OPENMC_BIOFILM_TRANSPORT — no measured density or "
            "closed composition exists to assign")
    else:
        reasons.append("material OpenMC gate ready")

    if binding.lattice_pitch_um is None:
        blockers.append("lattice_pitch_um is unset")
    if binding.density_g_cm3 is None:
        blockers.append("density_g_cm3 is unset")

    return IntegrationVerdict(
        can_emit_config=not blockers,
        binding_coherent=not problems,
        blockers=tuple(blockers), reasons=tuple(reasons))


# `emit_transport_config` lives in `reference_system.py`, not here. It needs the
# apparatus as well as the voxel — the membrane thickness, the cylinder, the
# boundaries, the source and all three material classes — and this module's
# subject is deliberately only the voxel. `evaluate()` above still answers the
# binding-and-gates question on its own, which is what the gate reports print.
