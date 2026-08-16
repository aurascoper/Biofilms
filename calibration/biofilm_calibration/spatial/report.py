"""The spatial calibration gate.

Returns one of four verdicts. The default is NOT_EVALUATED and the burden is on
the evidence to move it, not the other way round.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

from ..schema import SchemaError, read_table
from .. import approval
from ..acquisition import BASELINE_CONDITION
from .schema import (BIOFILM_STRUCTURE, ENTITY_SEMANTICS, OBJECT_MORPHOLOGY,
                    SAMPLE_METADATA, SOURCES)
from .time_observable import TIME_OBSERVABLE
from .time_observable import problems as observable_problems

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


def evaluate(data_dir, acceptance_path=None, *, today=None) -> SpatialGate:
    """Evaluate the gate against whatever evidence exists on disk."""
    data_dir = Path(data_dir)
    reasons: list[str] = []
    blockers: list[str] = []

    entities, err = _load(data_dir / "entity_semantics.csv", ENTITY_SEMANTICS)
    if err:
        return SpatialGate(NOT_EVALUATED, (), (f"entity semantics unreadable: {err}",))

    # Explicit allowlist: any status outside these two blocks the gate,
    # including unsupported_by_current_model.
    usable = {"ready", "provisional"}
    unresolved = [r["species_id"] for r in entities
                  if r["entity_kind"] == "unresolved"
                  or r["status"] not in usable]
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

    if samples:
        unpaired = [r["sample_id"] for r in samples
                    if not (r.get("paired_sample_group") or "").strip()]
        if unpaired:
            blockers.append(
                f"{len(unpaired)} imaged sample(s) have no paired_sample_group "
                "— a pitch fitted to an unpaired stack describes a different "
                "biofilm from the one the density was measured on")
        unresolved_basis = [r["sample_id"] for r in samples
                            if (r.get("segmentation_basis") or "").strip()
                            in ("", "unresolved")]
        if unresolved_basis:
            blockers.append(
                f"{len(unresolved_basis)} imaged sample(s) have no declared "
                "segmentation_basis — V_hydrated is the integral of this field "
                "and becomes the denominator of rho_wet, so what the mask "
                "encloses decides whether the density is meaningful")
        if not any((r.get("held_out") or "").strip() == "true" for r in samples):
            blockers.append(
                "no sample is designated held_out — a pitch fitted and "
                "validated on the same stacks has not been validated")

    # A surrogate exercises the harness and cannot clear the target gate,
    # and no campaign may be designed without an institutional approval id.
    baseline, err_b = _load(data_dir.parent / "baseline_condition.csv",
                            BASELINE_CONDITION)
    if err_b:
        blockers.append(f"baseline_condition unreadable: {err_b}")
    elif not baseline:
        blockers.append(
            "no baseline condition declared — every calibration value inherits "
            "the conditions it was measured under")
    else:
        # AN APPROVAL IS EVIDENCE, NOT A DECLARATION. This used to be a presence
        # check, and since the schema already refuses a blank it was very nearly
        # a no-op: any non-empty string cleared it, including the literal words
        # "d approved". `approval.problems` checks what a reviewer would ask --
        # authority, a resolvable document of the right kind, dates bracketing
        # the culturing, and a scope digest that still matches this row.
        sources, _ = _load(data_dir / "sources.csv", SOURCES)
        blockers.extend(approval.problems(baseline, sources or [], today=today))

    # READY_FOR_TIME_CALIBRATION asserts the time calibration can BEGIN, and
    # dt_MCS = a^2 * S_sim / S_exp needs both a pitch and an observable the
    # model actually represents. A pitch alone is not enough.
    observables, err_t = _load(data_dir / "time_observable.csv", TIME_OBSERVABLE)
    if err_t:
        blockers.append(f"time_observable unreadable: {err_t}")
    else:
        blockers.extend(observable_problems(observables))

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
