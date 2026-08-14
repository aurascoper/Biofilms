"""Table schemas for the biofilm material calibration."""

from __future__ import annotations

from ..acquisition import linkage_columns
from ..schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

# What a number is measured PER. Mirrors basis_conversion.BASES, plus the
# analytical bases that only appear in tables.
ANALYSIS_BASIS = frozenset({
    "wet_bulk", "dry_biomass", "ash", "medium", "membrane_dry",
    "membrane_hydrated", "site_occupancy",
})
ANALYSIS_METHOD = frozenset({
    "measured", "by_difference", "manufacturer_stated", "stoichiometric",
    "unresolved",
})
QUALITY = frozenset({"ok", "suspect", "excluded"})
BOOLEAN = frozenset({"true", "false"})
MATERIAL_MODEL_KIND = frozenset({"hydrated_effective_medium", "explicit_components"})

SAMPLE_METADATA = TableSchema(
    name="material_sample_metadata",
    doc="growth and handling conditions for one material sample",
    comments=(
        "Composition varies with growth conditions and developmental stage, so",
        "a composition without its conditions cannot be judged applicable to",
        "anything.",
    ),
    key=("sample_id",),
    columns=(
        Column("sample_id"),
        *linkage_columns(),
        Column("description"),
        Column("species_composition", required=False),
        Column("medium"),
        Column("temperature_C", numeric=True),
        Column("pH", numeric=True),
        Column("oxygen_condition"),
        Column("growth_phase"),
        Column("biofilm_age_h", numeric=True, required=False),
        Column("replicate_id", required=False),
        Column("source_id"),
    ),
)

BULK_MEASUREMENTS = TableSchema(
    name="bulk_measurements",
    doc="wet/dry/ash masses and hydrated volume",
    comments=(
        "The surface-water procedure is LOAD-BEARING and is therefore six",
        "columns rather than one free-text note: a wet mass without a stated",
        "drain orientation, drain time, blot material, contact time, ambient",
        "temperature and time-to-weighing is not reproducible, and the surface",
        "water it includes goes straight into rho_wet and w_water.",
        "drying_endpoint likewise — 'dried overnight' and 'dried to constant",
        "mass' are different measurements.",
        "Masses are AS-WEIGHED, including the substrate. The biofilm mass is",
        "m_sample+substrate - m_blank, so the blank is part of the definition",
        "and blank_sample_id must resolve in blanks.csv.",
    ),
    key=("sample_id", "replicate_id"),
    columns=(
        Column("sample_id"),
        Column("replicate_id"),
        *linkage_columns(),
        # As-weighed, INCLUDING the substrate. The biofilm mass is derived by
        # blank subtraction, so the raw reading is what gets recorded.
        Column("wet_mass_sample_plus_substrate_g", numeric=True),
        Column("dry_mass_sample_plus_substrate_g", numeric=True),
        Column("ash_mass_g", numeric=True, required=False),
        Column("hydrated_volume_cm3", numeric=True),
        Column("wet_mass_uncertainty_g", numeric=True, required=False),
        Column("dry_mass_uncertainty_g", numeric=True, required=False),
        Column("volume_uncertainty_cm3", numeric=True, required=False),
        Column("volume_method"),
        # "Wet mass" is not a measurement until the surface water is defined.
        Column("drain_orientation"),
        Column("drain_time_s", numeric=True),
        Column("blot_material"),
        Column("blot_contact_time_s", numeric=True),
        Column("ambient_temperature_C", numeric=True),
        Column("time_to_weighing_s", numeric=True),
        Column("drying_protocol"),
        Column("drying_endpoint"),
        Column("ash_protocol", required=False),
        Column("quality_flag", vocabulary=QUALITY),
        Column("source_id"),
    ),
)

ELEMENTAL_ANALYSIS = TableSchema(
    name="elemental_analysis",
    doc="elemental mass fractions, on a declared basis",
    comments=(
        "Never mix dry-basis CHNS values with wet-basis water fractions in one",
        "sum without an explicit conversion — that is what the basis column is",
        "for.",
        "Oxygen is usually obtained BY DIFFERENCE, which is legitimate and must",
        "be recorded: by-difference oxygen absorbs the error of every other",
        "element and of the ash determination.",
        "These are ELEMENTS, not isotopes. An OpenMC material built by element",
        "uses natural abundance.",
    ),
    key=("sample_id", "basis", "element"),
    columns=(
        Column("sample_id"),
        Column("basis", vocabulary=ANALYSIS_BASIS),
        Column("element"),
        Column("mass_fraction", numeric=True),
        Column("uncertainty", numeric=True, required=False),
        Column("analysis_method", vocabulary=ANALYSIS_METHOD),
        Column("detection_limit", numeric=True, required=False),
        Column("below_detection_limit", vocabulary=BOOLEAN),
        Column("source_id"),
    ),
)

METAL_LOADING = TableSchema(
    name="metal_loading",
    doc="sorbed, precipitated and intracellular metal, per declared basis",
    comments=(
        "taxonomic_fraction and active_fraction are SEPARATE columns on",
        "purpose. The abundance of a metal-reducing species is an upper bound",
        "on the active reducer fraction, not equal to it: cells present but not",
        "respiring the metal do not reduce it.",
        "speciation_state distinguishes dissolved, sorbed, precipitated and",
        "intracellular metal, which are not interchangeable for transport.",
    ),
    key=("sample_id", "element", "speciation_state"),
    columns=(
        Column("sample_id"),
        Column("element"),
        Column("speciation_state"),
        Column("basis", vocabulary=ANALYSIS_BASIS),
        Column("loading_value", numeric=True),
        Column("loading_unit"),
        Column("uncertainty", numeric=True, required=False),
        Column("taxonomic_fraction", numeric=True, required=False),
        Column("active_fraction", numeric=True, required=False),
        Column("analysis_method", vocabulary=ANALYSIS_METHOD),
        Column("counted_in_ash", vocabulary=BOOLEAN),
        Column("source_id"),
    ),
)

COMPONENT_DEFINITIONS = TableSchema(
    name="component_definitions",
    doc="the phases the hydrated effective medium is blended from",
    comments=(
        "One row per component of the wet bulk. mass_fraction values must sum",
        "to 1 across the components of a scenario, or a forgotten component is",
        "absorbed as a silent renormalization.",
    ),
    key=("scenario_id", "component"),
    columns=(
        Column("scenario_id"),
        Column("component"),
        Column("mass_fraction", numeric=True),
        Column("density_g_cm3", numeric=True, required=False),
        Column("density_evidence", vocabulary=EVIDENCE_BASIS, required=False),
        Column("contributes_metals", required=False),
        Column("material_model_kind", vocabulary=MATERIAL_MODEL_KIND),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id", required=False),
        Column("status", vocabulary=STATUS),
        Column("notes", required=False),
    ),
)

SOURCES = TableSchema(
    name="material_sources",
    doc="branch-local source registry for the material calibration",
    comments=(
        "Branch-local: merging into data/sources.csv is the integration",
        "branch's job, so the two calibration branches never touch the same",
        "registry.",
    ),
    key=("source_id",),
    columns=(
        Column("source_id"),
        Column("title"),
        Column("publisher", required=False),
        Column("document_version", required=False),
        Column("url", required=False),
        Column("accessed_date", required=False),
        Column("notes", required=False),
    ),
)

ALL_TABLES = (SAMPLE_METADATA, BULK_MEASUREMENTS, ELEMENTAL_ANALYSIS,
              METAL_LOADING, COMPONENT_DEFINITIONS, SOURCES)
