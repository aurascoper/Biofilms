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
        "BIOSAFETY FOLLOWS STRAINS, NOT SPECIES. Several components are",
        "commonly BSL-1 while C. neoformans strains are commonly BSL-2, and",
        "collection material under one species name can carry different",
        "assigned levels. Current guidance is protocol-driven risk assessment",
        "rather than a level label, so the chain is: exact strain IDs ->",
        "institutional review -> containment decision -> approved procedures.",
        "biosafety_level_by_strain is KEYED BY THOSE EXACT IDs, not by an",
        "abbreviation of them: the level has to attach to the organism the",
        "committee reviewed, and only the identifier says which that is.",
        "AN APPROVAL IS EVIDENCE, NOT A DECLARATION. The identifier alone is",
        "not enough: the gate also requires an issuing authority, a resolvable",
        "approval document, dates bracketing the culturing, and a scope digest",
        "that matches the conditions on this row. Filler text is refused by",
        "value -- 'approved', 'pending', 'TBD' name nothing.",
    ),
    key=("growth_condition_id",),
    columns=(
        Column("growth_condition_id"),
        Column("consortium_or_surrogate"),
        # Exact strain identity comes FIRST, because biosafety follows strains
        # and not species: several components are commonly BSL-1 while
        # C. neoformans strains are commonly BSL-2, so assuming one facility
        # class from a species list is wrong.
        Column("strain_identities"),
        # THE KEY IS THE STRAIN IDENTIFIER, VERBATIM. Stated here because the
        # gate cannot invent it: an abbreviation like `DR:BSL1` declares no
        # relationship to `D. radiodurans R1`, so a consumer that read one from
        # initials would be inventing semantics this file never gave it. With
        # the identifier as its own key there is no convention to infer -- the
        # two fields name the same strings or the approval does not cover the
        # organism it claims to.
        Column("biosafety_level_by_strain",
               doc="per-strain, keyed by the strain_identities entry VERBATIM, "
                   "e.g. 'D. radiodurans R1:BSL1; C. neoformans H99:BSL2' — "
                   "never a single level for the consortium, and never an "
                   "abbreviation. Matched ignoring case and surrounding "
                   "whitespace only, so a strain identifier containing ':' "
                   "cannot be expressed and is refused at the gate"),
        # AN APPROVAL IS EVIDENCE, NOT A DECLARATION, so the identifier alone
        # is not enough: it must name an authority, resolve to a registered
        # document, carry dates that bracket the culturing, and bind the scope
        # it was issued for. All OPTIONAL here and enforced at the GATE, the
        # same split `linkage_columns()` uses and for the same reason — a
        # schema that refuses a blank forces every hand-built fixture row to
        # invent one, which is how placeholders get written in the first place.
        Column("institutional_approval_id",
               doc="biosafety committee approval; the gate refuses filler text "
                   "as well as a blank"),
        Column("institutional_approval_authority", required=False,
               doc="the committee or office that issued it — an identifier "
                   "with no issuer names nothing"),
        Column("approval_source_id", required=False,
               doc="resolves in data/calibration/spatial/sources.csv, and that "
                   "row must be document_type = institutional_biosafety_approval"),
        Column("approval_effective_date", required=False,
               doc="ISO date; must not be after culturing_start_date"),
        Column("approval_expiration_date", required=False,
               doc="ISO date; must be after culturing_start_date and after today"),
        Column("approved_protocol_version", required=False,
               doc="the protocol version the approval was issued against"),
        Column("approval_scope_hash", required=False,
               doc="digest of the approved scope columns. Detects a growth "
                   "condition edited after approval; does not detect forgery"),
        Column("culturing_start_date", required=False,
               doc="ISO date; the approval must already be effective on it"),
        Column("containment_facility"),
        Column("risk_assessment_reference", required=False,
               doc="protocol-driven assessment, not a level label alone"),
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


# What a segmented volume actually encloses. This is a DIFFERENT axis from the
# mass basis in materials/basis_conversion.py: that one says what a number is
# measured per, this one says what a volume contains.
SEGMENTATION_BASIS = frozenset({
    # Stained cells only. Excludes the matrix AND the pore water between
    # cells, so it is smaller than what a balance weighs — by whatever
    # fraction of the biofilm is EPS and interstitial liquid.
    "cells_only",
    # Cells plus stained extracellular matrix. Closer, but a union of stained
    # voxels still excludes water-filled internal voids, so it needs a declared
    # pore fraction before it can stand in for the envelope.
    "cells_and_matrix",
    # The enclosing biofilm envelope, internal voids included. The only basis
    # that matches the material a balance weighs.
    "whole_biofilm_envelope",
    "unresolved",
})

# The one basis that may be divided into a whole-biofilm wet mass directly.
WHOLE_BIOFILM_BASIS = "whole_biofilm_envelope"
