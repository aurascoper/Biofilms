"""Table schemas for the CPM spatial calibration.

Sample-level conditions are kept separate from object-level measurements on
purpose: a morphology number is meaningless without the medium, temperature,
growth phase and segmentation method that produced it, and repeating those on
every object row invites them to drift.
"""

from __future__ import annotations

from physical_contract import APPROVAL_DOCUMENT_TYPES

from ..acquisition import SEGMENTATION_BASIS, linkage_columns
from ..schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

ENTITY_KIND = frozenset({
    "biological_cell", "hyphal_compartment", "conidium", "biomass_parcel",
    "unresolved",
})
LINEAGE_SEMANTICS = frozenset({"biological", "compartmental", "computational",
                               "unresolved"})
SHAPE_TARGET = frozenset({"isotropic_blob", "rod", "filament", "unresolved"})
BOOLEAN = frozenset({"true", "false"})
QUALITY = frozenset({"ok", "suspect", "excluded"})

ENTITY_SEMANTICS = TableSchema(
    name="entity_semantics",
    doc="what one CPM cell ID represents, per species",
    comments=(
        "Declared in docs/calibration/cpm_entity_semantics.md.",
        "literal_generation_claim_allowed is false wherever the model cannot",
        "support the claim: compute_delta_H has no shape or connectivity term,",
        "V_target is one global scalar for all seven species, and divide_cell!",
        "has no trigger (biofilms_potts.jl:1578-1583).",
    ),
    key=("species_id",),
    columns=(
        Column("species_id"),
        Column("species_name"),
        Column("entity_kind", vocabulary=ENTITY_KIND),
        Column("division_unit"),
        Column("lineage_semantics", vocabulary=LINEAGE_SEMANTICS),
        Column("shape_target", vocabulary=SHAPE_TARGET),
        Column("literal_generation_claim_allowed", vocabulary=BOOLEAN),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id", required=False),
        Column("status", vocabulary=STATUS),
        Column("notes", required=False),
    ),
)

SAMPLE_METADATA = TableSchema(
    name="sample_metadata",
    doc="imaging and growth conditions for one morphology sample",
    comments=(
        "One row per imaged sample. Morphology measured under conditions other",
        "than those the model represents is a near-analog at best; the columns",
        "exist so that can be judged rather than assumed.",
        "held_out is a column, not an afterthought: a pitch fitted and",
        "validated on the same stacks has not been validated.",
        "segmentation_basis records WHAT the mask encloses. V_hydrated is",
        "the integral of this field, and it becomes the denominator of",
        "rho_wet — so a cells_only mask against a whole-biofilm wet mass is",
        "a silent, systematic error.",
    ),
    key=("sample_id",),
    columns=(
        Column("sample_id"),
        *linkage_columns(),
        Column("held_out", vocabulary=BOOLEAN),
        Column("species"),
        Column("strain", required=False),
        Column("medium"),
        Column("temperature_C", numeric=True),
        Column("pH", numeric=True),
        Column("oxygen_condition"),
        Column("growth_phase"),
        Column("biofilm_age_h", numeric=True, required=False),
        Column("imaging_method"),
        Column("voxel_size_x_um", numeric=True),
        Column("voxel_size_y_um", numeric=True),
        Column("voxel_size_z_um", numeric=True),
        Column("segmentation_method"),
        # WHAT was segmented, not how. The volume derived from this field
        # becomes the denominator of rho_wet, and a cell-only mask does not
        # enclose the material the balance weighs.
        Column("segmentation_basis", vocabulary=SEGMENTATION_BASIS),
        Column("segmentation_version", required=False),
        Column("replicate_id", required=False),
        Column("source_id"),
    ),
)

OBJECT_MORPHOLOGY = TableSchema(
    name="object_morphology",
    doc="per-segmented-object morphology, in micrometres",
    comments=(
        "One row per segmented object. Distributions are the unit of",
        "comparison, not means: a mean volume cannot distinguish a population",
        "of uniform cells from a mixture of two sizes.",
        "boundary_contact flags objects touching the image edge, whose volume",
        "and axis lengths are truncated and must not enter a size distribution.",
    ),
    key=("sample_id", "object_id"),
    columns=(
        Column("sample_id"),
        Column("object_id"),
        Column("volume_um3", numeric=True),
        Column("surface_area_um2", numeric=True, required=False),
        Column("equivalent_diameter_um", numeric=True, required=False),
        Column("major_axis_um", numeric=True, required=False),
        Column("intermediate_axis_um", numeric=True, required=False),
        Column("minor_axis_um", numeric=True, required=False),
        Column("aspect_ratio_major_minor", numeric=True, required=False),
        Column("connected_components", numeric=True, required=False),
        Column("boundary_contact", vocabulary=BOOLEAN),
        Column("measurement_uncertainty", numeric=True, required=False),
        Column("quality_flag", vocabulary=QUALITY),
    ),
)

BIOFILM_STRUCTURE = TableSchema(
    name="biofilm_structure",
    doc="sample-level biofilm structure: thickness, coverage, porosity",
    comments=(
        "Thickness and coverage are strongly heterogeneous across a biofilm,",
        "so a single sample row is a sample statistic and needs its own",
        "standard deviation, not just a mean.",
    ),
    key=("sample_id",),
    columns=(
        Column("sample_id"),
        Column("thickness_mean_um", numeric=True),
        Column("thickness_sd_um", numeric=True),
        Column("surface_coverage", numeric=True, required=False),
        Column("porosity", numeric=True, required=False),
        Column("biovolume_fraction", numeric=True, required=False),
        Column("aggregate_size_median_um", numeric=True, required=False),
        Column("radial_position_um", numeric=True, required=False),
        Column("axial_position_um", numeric=True, required=False),
        Column("measurement_uncertainty", numeric=True, required=False),
    ),
)

# Kinds of document this registry can pin. APPROVAL_DOCUMENT_TYPES comes from
# the contract so the vocabulary is shared rather than restated.
DOCUMENT_TYPES = frozenset(APPROVAL_DOCUMENT_TYPES | {
    "declaration", "publication", "dataset", "certificate", "evaluated_data"})


SOURCES = TableSchema(
    name="spatial_sources",
    doc="branch-local source registry for the spatial calibration",
    comments=(
        "Branch-local on purpose: merging this into data/sources.csv belongs",
        "to the integration branch, so the two calibration branches never",
        "touch the same registry file.",
    ),
    key=("source_id",),
    columns=(
        Column("source_id"),
        Column("title"),
        Column("publisher", required=False),
        Column("document_version", required=False),
        Column("url", required=False),
        Column("accessed_date", required=False),
        # WHAT KIND OF DOCUMENT THIS IS. Without it a source_id resolves happily
        # to any registered document, so an approval field could cite a
        # cross-section table and pass the check -- the "belongs to another
        # protocol" failure with nothing in the data able to see it.
        # `declaration` is the default for the in-repo declarations already
        # registered here.
        Column("document_type", required=False, vocabulary=DOCUMENT_TYPES),
        # Of the retrieved file, when one was retrieved. An approval that
        # cannot be located is not evidence, so the gate requires a url or a
        # digest on the row it resolves to.
        Column("sha256", required=False),
        Column("notes", required=False),
    ),
)

ALL_TABLES = (ENTITY_SEMANTICS, SAMPLE_METADATA, OBJECT_MORPHOLOGY,
              BIOFILM_STRUCTURE, SOURCES)
