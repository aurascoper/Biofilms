"""Two independent material gates.

They are separate because they need different things. OpenMC transport needs a
bulk density and a closed elemental composition. The radiodialysis model needs
mass CONCENTRATIONS — X_total and X_red in g/cm3 — and a rule for mapping them
onto space. A material can be good enough for transport and useless for
kinetics, so one verdict would have to be the worse of the two and would hide
which half was missing.

Ending `OPENMC: PROVISIONAL / RADIODIALYSIS: BLOCKED` is a valid, expected
result, not a failure of the branch.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..schema import SchemaError, read_table
from .schema import (BULK_MEASUREMENTS, COMPONENT_DEFINITIONS,
                     ELEMENTAL_ANALYSIS, METAL_LOADING, SAMPLE_METADATA)

READY_OPENMC = "READY_FOR_OPENMC_BIOFILM_TRANSPORT"
READY_RADIODIALYSIS = "READY_FOR_RADIODIALYSIS_MAPPING"
PROVISIONAL = "PROVISIONAL"
BLOCKED = "BLOCKED"
NOT_EVALUATED = "NOT_EVALUATED"


@dataclass(frozen=True)
class MaterialGates:
    openmc: str
    radiodialysis: str
    openmc_blockers: tuple
    radiodialysis_blockers: tuple
    reasons: tuple

    def render(self) -> str:
        lines = [f"OPENMC:        {self.openmc}",
                 f"RADIODIALYSIS: {self.radiodialysis}"]
        for r in self.reasons:
            lines.append(f"  + {r}")
        for b in self.openmc_blockers:
            lines.append(f"  - [openmc] {b}")
        for b in self.radiodialysis_blockers:
            lines.append(f"  - [radiodialysis] {b}")
        return "\n".join(lines)


def _load(path, schema):
    try:
        return read_table(path, schema), None
    except (OSError, SchemaError) as e:
        return None, str(e)


def evaluate(data_dir) -> MaterialGates:
    data_dir = Path(data_dir)
    reasons: list[str] = []
    om: list[str] = []
    rd: list[str] = []

    tables = {}
    for name, schema in (("sample_metadata", SAMPLE_METADATA),
                         ("bulk_measurements", BULK_MEASUREMENTS),
                         ("elemental_analysis", ELEMENTAL_ANALYSIS),
                         ("metal_loading", METAL_LOADING),
                         ("component_definitions", COMPONENT_DEFINITIONS)):
        rows, err = _load(data_dir / f"{name}.csv", schema)
        if err:
            return MaterialGates(NOT_EVALUATED, NOT_EVALUATED,
                                 (f"{name} unreadable: {err}",), (), ())
        tables[name] = rows

    if not tables["bulk_measurements"]:
        om.append("no bulk measurements — hydrated bulk density is undetermined")
        rd.append("no bulk measurements — X_total cannot be computed")
    else:
        mismatched = [r["sample_id"] for r in tables["bulk_measurements"]
                      if (r.get("volume_basis") or "").strip()
                      not in ("whole_biofilm_envelope",)
                      and not (r.get("pore_volume_fraction") or "").strip()]
        if mismatched:
            om.append(
                f"{len(mismatched)} bulk measurement(s) have a volume_basis "
                "that is not the whole-biofilm envelope and no declared "
                "pore_volume_fraction — the balance weighed cells, matrix and "
                "interstitial water, so dividing that mass by a cell-only or "
                "stained-voxel volume inflates rho_wet systematically")

        unblanked = [r["sample_id"] for r in tables["bulk_measurements"]
                     if not (r.get("blank_sample_id") or "").strip()]
        if unblanked:
            om.append(
                f"{len(unblanked)} bulk measurement(s) name no blank — blank "
                "subtraction defines the biofilm mass, so an unblanked sample "
                "attributes the substrate and its retained water to the biofilm")
    if not tables["elemental_analysis"]:
        om.append("no elemental analysis — composition is undetermined")
    if not tables["sample_metadata"]:
        om.append("no sample metadata — no composition can be judged applicable")
    else:
        unpaired = [r["sample_id"] for r in tables["sample_metadata"]
                    if not (r.get("paired_sample_group") or "").strip()]
        if unpaired:
            om.append(
                f"{len(unpaired)} material sample(s) have no paired_sample_group "
                "— an unpaired density describes a different biofilm from the "
                "one the pitch was fitted to")

    # Radiodialysis needs strictly more, and needs it in the right units.
    metal = tables["metal_loading"]
    if not metal:
        rd.append("no metal loading data — the active reducer fraction is "
                  "undefined, so X_red cannot be derived")
    else:
        untagged = [r for r in metal if r.get("active_fraction") is None]
        if untagged:
            rd.append(
                f"{len(untagged)} metal-loading row(s) have no active_fraction "
                "— taxonomic abundance is not functional activity and must not "
                "be substituted for it")

    rd.append("no spatial mapping rule declared: compute_radial_biomass returns "
              "site occupancy and its unweighted mean is assigned to X_total "
              "(biofilms_potts.jl L1341-1342, L1445-1446, L1797). Until that is "
              "replaced by a mass concentration, the radiodialysis mapping is "
              "blocked by a units error, not by missing data")

    if tables["component_definitions"]:
        kinds = {r["material_model_kind"] for r in tables["component_definitions"]}
        if kinds - {"hydrated_effective_medium"}:
            om.append(
                f"component definitions declare {sorted(kinds)}; only "
                "hydrated_effective_medium is compatible with the current "
                "one-material-class-per-occupied-voxel mapper")
        else:
            reasons.append("material_model_kind = hydrated_effective_medium")
    else:
        om.append("no component definitions — nothing to blend")

    reasons.append("basis tagging enforced: site occupancy cannot reach a "
                   "concentration slot")

    openmc = PROVISIONAL if om else READY_OPENMC
    radiodialysis = BLOCKED if rd else READY_RADIODIALYSIS
    return MaterialGates(openmc, radiodialysis, tuple(om), tuple(rd),
                         tuple(reasons))
