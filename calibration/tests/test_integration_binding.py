"""The integration contract: coherence, and the refusal to emit.

The mixed-semantics rejections are the point. Each incoherent binding below is
built from clauses that pass their OWN branch's gate — that is exactly why the
check has to live at the join.
"""

from __future__ import annotations

import pytest

from biofilm_calibration.integration import (VoxelBinding, coherence_problems,
                                             evaluate, is_coherent)
from biofilm_calibration.integration_schema import VOXEL_BINDING
from biofilm_calibration.materials import report as material_report
from biofilm_calibration.schema import read_table
from biofilm_calibration.spatial import report as spatial_report

COHERENT = VoxelBinding(
    entity_kind="biomass_parcel",
    lineage_semantics="computational",
    material_class="baseline_biomass",
    material_represents="hydrated_bulk_biofilm",
    material_model_kind="hydrated_effective_medium",
    density_basis="wet_bulk",
)


def test_the_declared_binding_is_coherent():
    assert is_coherent(COHERENT), coherence_problems(COHERENT)


def test_organism_scale_entity_with_bulk_biofilm_material_is_rejected():
    """The reviewer's incoherent example. Both halves pass their own gate: the
    spatial branch would accept a literal cell, the material branch would
    accept a bulk wet composition. Only the join sees the contradiction."""
    b = VoxelBinding(**{**COHERENT.__dict__, "entity_kind": "biological_cell",
                        "lineage_semantics": "biological"})
    problems = coherence_problems(b)
    assert any("cell interior is not bulk biofilm" in p for p in problems)


def test_parcel_with_cytoplasm_material_is_rejected():
    """The mirror image: a parcel bundles cells with water and EPS, so its
    contents cannot be cytoplasm."""
    b = VoxelBinding(**{**COHERENT.__dict__,
                        "material_represents": "cytoplasm"})
    assert any("cannot be cytoplasm" in p for p in coherence_problems(b))


def test_biological_lineage_on_a_parcel_is_rejected():
    b = VoxelBinding(**{**COHERENT.__dict__, "lineage_semantics": "biological"})
    problems = coherence_problems(b)
    assert any("parcel ancestry is not organism ancestry" in p for p in problems)
    assert any("1578-1583" in p for p in problems)


def test_computational_lineage_on_an_organism_is_rejected():
    b = VoxelBinding(**{**COHERENT.__dict__,
                        "entity_kind": "hyphal_compartment",
                        "material_represents": "cytoplasm"})
    assert any("not merely computational" in p for p in coherence_problems(b))


def test_explicit_components_cannot_be_bound_to_a_voxel():
    b = VoxelBinding(**{**COHERENT.__dict__,
                        "material_model_kind": "explicit_components"})
    assert any("one material class per occupied voxel" in p
               for p in coherence_problems(b))


def test_dry_solid_contradicts_the_hydrated_effective_medium():
    b = VoxelBinding(**{**COHERENT.__dict__, "material_represents": "dry_solid"})
    assert any("includes the water" in p for p in coherence_problems(b))


def test_site_occupancy_is_not_a_density():
    """The live model defect, caught at the binding as well as in the
    converters."""
    b = VoxelBinding(**{**COHERENT.__dict__, "density_basis": "site_occupancy"})
    problems = coherence_problems(b)
    assert any("not a density at all" in p for p in problems)


def test_coherent_is_not_sufficient_to_emit():
    """A coherent sentence with unmeasured blanks is still a sentence with
    blanks."""
    v = evaluate(COHERENT, spatial_report.PROVISIONAL,
                 material_report.PROVISIONAL)
    assert v.binding_coherent
    assert not v.can_emit_config
    joined = " ".join(v.blockers)
    assert "spatial gate is PROVISIONAL" in joined
    assert "material OpenMC gate is PROVISIONAL" in joined


def test_unset_pitch_or_density_still_blocks_a_ready_pair_of_gates():
    v = evaluate(COHERENT, spatial_report.READY, material_report.READY_OPENMC)
    assert not v.can_emit_config
    assert any("lattice_pitch_um is unset" in b for b in v.blockers)
    assert any("density_g_cm3 is unset" in b for b in v.blockers)


def test_live_gates_still_refuse(data_dir, config_dir):
    """Run against the repository's actual state, not a fixture."""
    s = spatial_report.evaluate(
        data_dir / "spatial", config_dir / "cpm_spatial_acceptance_template.toml")
    m = material_report.evaluate(data_dir / "materials")
    v = evaluate(COHERENT, s.verdict, m.openmc)
    assert v.binding_coherent
    assert not v.can_emit_config


def test_shipped_binding_matches_the_declared_one(data_dir):
    """The contract in data/ must be the contract the code checks."""
    from biofilm_calibration.integration import VoxelBinding as VB
    from biofilm_calibration.spatial.schema import ENTITY_SEMANTICS

    rows = read_table(data_dir / "voxel_binding.csv", VOXEL_BINDING)
    assert len(rows) == 1
    row = rows[0]
    shipped = VB(entity_kind=row["entity_kind"],
                 lineage_semantics=row["lineage_semantics"],
                 material_class=row["material_class"],
                 material_represents=row["material_represents"],
                 material_model_kind=row["material_model_kind"],
                 density_basis=row["density_basis"])
    assert is_coherent(shipped), coherence_problems(shipped)

    # and it must agree with what the spatial branch declared per species
    species = read_table(data_dir / "spatial" / "entity_semantics.csv",
                         ENTITY_SEMANTICS)
    assert {r["entity_kind"] for r in species} == {shipped.entity_kind}
    assert {r["lineage_semantics"] for r in species} == {shipped.lineage_semantics}

