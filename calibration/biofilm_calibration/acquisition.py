"""Acquisition contract: the baseline condition, sample pairing, and blanks.

Cross-domain on purpose. A lattice pitch comes from imaging and a density comes
from gravimetry, and unless those two measurements are linked to one culture
batch they describe two different biofilms — with nothing downstream able to
notice.

The six linkage columns below are added to BOTH branches' sample tables. They
are the whole mechanism:

    culture_batch_id      the unit of biological replication
    paired_sample_group   links an imaged coupon to its weighed sibling
    blank_sample_id       which blank corrects this sample
    measurement_order     detects drift within a weighing session
    growth_condition_id   joins to the frozen baseline condition
    medium_batch_id       medium salts end up in the ash and the dry mass
"""

from __future__ import annotations

from .schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

BLANK_KIND = frozenset({
    "dry_substrate",           # substrate tare and manufacturing variation
    "hydrated_substrate",      # water retained by the membrane or support
    "medium_exposed_abiotic",  # salts or precipitates deposited without biology
})
DOMAIN_SEMANTICS = frozenset({"microvolume", "representative_segment",
                              "full_apparatus", "unresolved"})
BOOLEAN = frozenset({"true", "false"})


def linkage_columns() -> tuple:
    """The six cross-domain columns, defined once so the two branches cannot
    drift apart on their names or their meaning.

    All OPTIONAL at the schema level and enforced at the GATES instead. The
    tables must be able to record reality — including a pilot sample with no
    sibling coupon — while the gate refuses to calibrate from it. Making them
    required here would prevent an unpaired sample from being written down at
    all, which loses information rather than protecting anything.
    """
    return (
        Column("culture_batch_id", required=False,
               doc="unit of biological replication; gate-enforced"),
        Column("paired_sample_group", required=False,
               doc="imaged coupon and its weighed sibling share this; gate-enforced"),
        Column("blank_sample_id", required=False,
               doc="which blank corrects this sample; gate-enforced"),
        Column("measurement_order", numeric=True, required=False,
               doc="position within a weighing session; detects drift"),
        Column("growth_condition_id", required=False,
               doc="joins to baseline_condition.csv; gate-enforced"),
        Column("medium_batch_id", required=False,
               doc="medium salts reach the ash and the dry mass"),
    )


BASELINE_CONDITION = TableSchema(
    name="baseline_condition",
    doc="the one frozen calibration specimen, declared before collection",
    comments=(
        "Ships UNSET. This is a protocol contract, not a claim about a",
        "specimen that exists.",
        "domain_semantics = microvolume for the first calibration: the CPM",
        "forces L = 2R (R = N/2 recomputed as a local at seven sites in",
        "biofilms_potts.jl, no length field anywhere), so no pitch can make",
        "the cubic lattice reproduce a full apparatus.",
        "A single-species or surrogate specimen exercises the harness but",
        "cannot clear the gate for the seven-species target system.",
    ),
    key=("growth_condition_id",),
    columns=(
        Column("growth_condition_id"),
        Column("consortium_or_surrogate"),
        Column("strain_identities"),
        Column("growth_medium"),
        Column("temperature_C", numeric=True),
        Column("pH", numeric=True),
        Column("oxygen_condition"),
        Column("substrate_or_membrane"),
        Column("biofilm_age_h", numeric=True),
        Column("flow_condition"),
        Column("irradiation"),
        Column("domain_semantics", vocabulary=DOMAIN_SEMANTICS),
        Column("is_target_system", vocabulary=BOOLEAN,
               doc="false for a surrogate: exercises the harness, cannot clear the gate"),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id", required=False),
        Column("status", vocabulary=STATUS),
        Column("notes", required=False),
    ),
)

BLANKS = TableSchema(
    name="blanks",
    doc="matched blanks, weighed through the identical procedure",
    comments=(
        "Blank correction is the DEFINITION of a biofilm mass, not a",
        "refinement:",
        "  m_wet,biofilm = m_wet,sample+substrate - m_wet,blank_substrate",
        "The hydrated_substrate blank is the one most easily skipped and the",
        "most costly to skip: a hydrated membrane retains substantial water,",
        "which would otherwise be counted as biofilm water — inflating rho_wet",
        "and deflating w_water.",
        "The medium_exposed_abiotic blank separates medium salts from biomass",
        "ash.",
    ),
    key=("blank_sample_id",),
    columns=(
        Column("blank_sample_id"),
        Column("blank_kind", vocabulary=BLANK_KIND),
        Column("culture_batch_id", required=False),
        Column("growth_condition_id"),
        Column("substrate_or_membrane"),
        Column("wet_mass_g", numeric=True, required=False),
        Column("dry_mass_g", numeric=True, required=False),
        Column("wet_mass_uncertainty_g", numeric=True, required=False),
        Column("dry_mass_uncertainty_g", numeric=True, required=False),
        Column("measurement_order", numeric=True, required=False),
        Column("source_id"),
    ),
)

ALL_TABLES = (BASELINE_CONDITION, BLANKS)
