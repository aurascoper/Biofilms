"""The spatial calibration gate.

Returns one of four verdicts. The default is NOT_EVALUATED and the burden is on
the evidence to move it, not the other way round.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

from ..schema import SchemaError, read_table
from .schema import (BIOFILM_STRUCTURE, ENTITY_SEMANTICS, OBJECT_MORPHOLOGY,
                     SAMPLE_METADATA)

# Thresholds that must be declared before any pitch can be selected. Kept in
# step with config/cpm_spatial_acceptance_template.toml and with
# scale_candidates.select(), which refuses on the same set.
REQUIRED_THRESHOLDS = (
    "minimum_minor_axis_voxels",
    "maximum_volume_quantization_error",
    "maximum_axis_ratio_error",
    "maximum_memory_bytes",
    "minimum_objects_per_species",
)

READY = "READY_FOR_TIME_CALIBRATION"
PROVISIONAL = "PROVISIONAL"
MODEL_REVISION_REQUIRED = "MODEL_REVISION_REQUIRED"
NOT_EVALUATED = "NOT_EVALUATED"


@dataclass(frozen=True)
class SpatialGate:
    verdict: str
    reasons: tuple
    blockers: tuple

    def render(self) -> str:
        lines = [f"spatial gate: {self.verdict}"]
        for r in self.reasons:
            lines.append(f"  + {r}")
        for b in self.blockers:
            lines.append(f"  - {b}")
        return "\n".join(lines)


def _load(path, schema):
    try:
        return read_table(path, schema), None
    except (OSError, SchemaError) as e:
        return None, str(e)


def evaluate(data_dir, acceptance_path=None) -> SpatialGate:
    """Evaluate the gate against whatever evidence exists on disk."""
    data_dir = Path(data_dir)
    reasons: list[str] = []
    blockers: list[str] = []

    entities, err = _load(data_dir / "entity_semantics.csv", ENTITY_SEMANTICS)
    if err:
        return SpatialGate(NOT_EVALUATED, (), (f"entity semantics unreadable: {err}",))

    unresolved = [r["species_id"] for r in entities
                  if r["entity_kind"] == "unresolved"
                  or r["status"] not in {"ready", "provisional"}]
    if unresolved:
        blockers.append(
            f"entity semantics unresolved for: {', '.join(unresolved)}")
    else:
        kinds = {r["entity_kind"] for r in entities}
        reasons.append(f"entity semantics declared for all "
                       f"{len(entities)} species: {', '.join(sorted(kinds))}")

    # A literal-cell claim the model cannot support is a model problem, not a
    # calibration one, and it outranks every other verdict.
    literal = [r["species_id"] for r in entities
               if r["literal_generation_claim_allowed"] == "true"]
    if literal:
        return SpatialGate(
            MODEL_REVISION_REQUIRED,
            tuple(reasons),
            (f"literal generation claims declared for {', '.join(literal)}, "
             "but compute_delta_H has no shape or connectivity term, V_target "
             "is one global scalar, and divide_cell! has no trigger. Per-species "
             "target volumes, a shape constraint and a growth-linked division "
             "rule are model changes, not calibration.",))

    samples, err = _load(data_dir / "sample_metadata.csv", SAMPLE_METADATA)
    objects, err_o = _load(data_dir / "object_morphology.csv", OBJECT_MORPHOLOGY)
    structure, err_s = _load(data_dir / "biofilm_structure.csv", BIOFILM_STRUCTURE)
    for label, table, e in (("sample_metadata", samples, err),
                            ("object_morphology", objects, err_o),
                            ("biofilm_structure", structure, err_s)):
        if e:
            blockers.append(f"{label} unreadable: {e}")
        elif not table:
            blockers.append(f"{label} is empty — no measured morphology exists")

    acceptance = None
    if acceptance_path and Path(acceptance_path).exists():
        with open(acceptance_path, "rb") as fh:
            acceptance = tomllib.load(fh)
        rep = acceptance.get("representation", {})
        # An unset threshold is an ABSENT key, not a None value — the template
        # ships them commented out, exactly as config/coupling_template.toml
        # does. Checking for None would read the template as fully declared.
        unset = [k for k in REQUIRED_THRESHOLDS
                 if rep.get(k) is None]
        if unset:
            blockers.append(
                f"acceptance thresholds unset in {Path(acceptance_path).name}: "
                + ", ".join(unset))
        else:
            reasons.append("acceptance thresholds declared")
        dom = acceptance.get("domain", {}).get("semantics")
        if not dom or dom == "unresolved":
            blockers.append("domain semantics undeclared — the CPM forces "
                            "L = 2R, so a full apparatus cannot be claimed")
        else:
            reasons.append(f"domain semantics declared: {dom}")
    else:
        blockers.append("no acceptance configuration supplied")

    if blockers:
        return SpatialGate(PROVISIONAL, tuple(reasons), tuple(blockers))
    return SpatialGate(READY, tuple(reasons), ())
