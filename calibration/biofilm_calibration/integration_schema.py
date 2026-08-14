"""Schema for the shipped voxel-binding contract."""

from __future__ import annotations

from .integration import (ENTITY_KINDS, LINEAGE_SEMANTICS, MATERIAL_MODEL_KINDS,
                          MATERIAL_REPRESENTS)
from .schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

DENSITY_BASIS = frozenset({"wet_bulk", "dry_biomass", "site_occupancy",
                           "unresolved"})

VOXEL_BINDING = TableSchema(
    name="voxel_binding",
    doc="the integration contract: what an occupied CPM voxel is, and is made of",
    comments=(
        "One row. It completes the sentence the two calibration branches meet",
        "at, and biofilm_calibration.integration.coherence_problems() rejects",
        "the combinations that contradict themselves — combinations whose",
        "individual clauses each pass their own branch's gate.",
        "lattice_pitch_um and density_g_cm3 are blank because neither has been",
        "measured; blank means NOT MEASURED and is never read as zero.",
    ),
    key=("binding_id",),
    columns=(
        Column("binding_id"),
        Column("entity_kind", vocabulary=ENTITY_KINDS),
        Column("lineage_semantics", vocabulary=LINEAGE_SEMANTICS),
        Column("material_class"),
        Column("material_represents", vocabulary=MATERIAL_REPRESENTS),
        Column("material_model_kind", vocabulary=MATERIAL_MODEL_KINDS),
        Column("density_basis", vocabulary=DENSITY_BASIS),
        Column("lattice_pitch_um", numeric=True, required=False),
        Column("density_g_cm3", numeric=True, required=False),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id"),
        Column("status", vocabulary=STATUS),
        Column("notes", required=False),
    ),
)
