#!/usr/bin/env python
"""What Reference D still needs, asked of the contract rather than of a document.

    python calibration/scripts/reference_d_status.py
    python calibration/scripts/reference_d_status.py --check     # CI mode

A measurement protocol written as prose drifts from the code the moment either
changes, and the drift is invisible: the document keeps listing a field that is
no longer required, or stops listing one that is, and nobody notices until a
campaign comes back missing something. So the protocol is not the authority
here. The contract is:

    - the two gates, evaluated against whatever is actually on disk
    - emit_transport_config's own refusal, under evidence_policy = measured_only

and `data/calibration/reference_d_requirements.csv` supplies only the half the
code cannot know — who can produce each input, by what artifact, and what
counts as acceptance. `--check` fails if the two disagree in either direction:
a required field with no row, or a row naming a field nothing requires.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from biofilm_calibration.integration import VoxelBinding, coherence_problems
from biofilm_calibration.integration_schema import VOXEL_BINDING
from biofilm_calibration.materials import report as material_report
from biofilm_calibration.reference_system import ReferenceSystemSpec, emit_problems
from biofilm_calibration.schema import read_table
from biofilm_calibration.spatial import report as spatial_report

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data" / "calibration"
REQUIREMENTS = DATA / "reference_d_requirements.csv"

SUPPLIED_BY = ("measured", "evaluated_data", "derived", "declared")

# What each requirement is waiting FOR, because "not started" hid the only
# distinction that matters for scheduling: an institutional approval and a
# numerical study are not measurements, and neither is a modelling decision.
STATUSES = (
    "awaiting_declaration",          # a modeller must decide it
    "awaiting_approval",             # an institution must issue it
    "awaiting_measurement",          # a laboratory must produce it
    "awaiting_derivation",           # computable once its inputs exist
    "awaiting_analysis",             # needs a study, not a measurement
    "awaiting_evaluated_data",       # published data, looked up
    "unsupported_by_current_model",  # needs a model change, not data
    "satisfied",
)

# Which states block which threshold. THE CAMPAIGN MUST NOT WAIT ON ITS OWN
# OUTPUTS: D-NUMERICS and D-AXIALFACE are `awaiting_analysis` precisely because
# they depend on data the campaign produces, so they cannot gate its start.
# An unmade modelling decision does gate it — measurements would answer a
# question nobody had framed. So does approval, which gates culturing itself.
_BLOCKS_CAMPAIGN = ("awaiting_declaration", "awaiting_approval")

CAMPAIGN_READY = "CAMPAIGN_READY"
CONFIG_READY = "REFERENCE_D_CONFIG_READY"
SWEEP_READY = "BIOFILM_SWEEP_READY"
NOT_READY = "NOT_READY"

# Frozen acceptance configs, read for their [enforcement] tables. A criterion
# nothing compares against may be declared, but it may never close a gate.
ACCEPTANCE_CONFIGS = ("reference_d_spatial_acceptance.toml",
                      "reference_d_material_acceptance.toml")

# Fields the emitter names that the ledger addresses under a different key,
# because one artifact closes several at once and splitting the row would make
# the protocol read as more work than it is.
_ALIASES = {
    "cylinder_radius_cm": "cylinder_radius_cm / cylinder_length_cm",
    "cylinder_length_cm": "cylinder_radius_cm / cylinder_length_cm",
    "spectrum_energies_eV": "spectrum_energies_eV / spectrum_probabilities",
    "spectrum_probabilities": "spectrum_energies_eV / spectrum_probabilities",
    "source_spatial": "source_spatial / source_angular",
    "source_angular": "source_spatial / source_angular",
    "axial_bc": "axial_bc / radial_outer_bc / origin_cm",
    "radial_outer_bc": "axial_bc / radial_outer_bc / origin_cm",
    "origin_cm": "axial_bc / radial_outer_bc / origin_cm",
    "batches": "batches / particles / seed / mesh_coarsening_factor",
    "particles": "batches / particles / seed / mesh_coarsening_factor",
    "seed": "batches / particles / seed / mesh_coarsening_factor",
    "mesh_coarsening_factor": "batches / particles / seed / mesh_coarsening_factor",
}


def load_requirements() -> list[dict]:
    with open(REQUIREMENTS, encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    problems = []
    seen = set()
    for r in rows:
        rid = r["requirement_id"]
        if rid in seen:
            problems.append(f"duplicate requirement_id {rid!r}")
        seen.add(rid)
        if r["supplied_by"] not in SUPPLIED_BY:
            problems.append(f"{rid}: supplied_by {r['supplied_by']!r} not in {SUPPLIED_BY}")
        if r["status"] not in STATUSES:
            problems.append(f"{rid}: status {r['status']!r} not in {STATUSES}")
        for col in ("artifact", "acceptance_criterion", "unblocks"):
            if not (r.get(col) or "").strip():
                problems.append(f"{rid}: {col} is empty")
    if problems:
        raise ValueError("reference_d_requirements.csv is malformed:\n  - "
                         + "\n  - ".join(problems))
    return rows


def enforcement_report() -> dict[str, list[str]]:
    """Which declared criteria a consumer actually compares against.

    A gate that closes on a criterion nobody evaluated is a gate that closed on
    nothing. The frozen configs classify every threshold, and this surfaces the
    ones that are decoration so a reader is not misled by their presence.
    """
    import tomllib

    out: dict[str, list[str]] = {"hard": [], "advisory": [], "derived": [],
                                 "not_implemented": []}
    for name in ACCEPTANCE_CONFIGS:
        path = REPO / "config" / name
        if not path.exists():
            continue
        with open(path, "rb") as fh:
            data = tomllib.load(fh)
        for key, level in (data.get("enforcement") or {}).items():
            out.setdefault(level, []).append(f"{name.split('_')[2]}:{key}")
    return out


def readiness(requirements, spatial_verdict, material_verdict,
              binding) -> dict[str, tuple[bool, list[str]]]:
    """The three thresholds, which gate different things and must not be one
    flag. A campaign can be ready while a config is not; a config can be ready
    while the sweep has not run."""
    by_status: dict[str, list[str]] = {}
    for r in requirements:
        by_status.setdefault(r["status"], []).append(r["requirement_id"])

    campaign_blockers = [f"{s}: {', '.join(sorted(by_status[s]))}"
                         for s in _BLOCKS_CAMPAIGN if by_status.get(s)]

    config_blockers = []
    unsatisfied = [r["requirement_id"] for r in requirements
                   if r["status"] not in ("satisfied", "unsupported_by_current_model")
                   and r["unblocks"] != "radiodialysis gate"]
    if unsatisfied:
        config_blockers.append(f"{len(unsatisfied)} requirement(s) unsatisfied")
    if spatial_verdict != spatial_report.READY:
        config_blockers.append(f"spatial gate is {spatial_verdict}")
    if material_verdict != material_report.READY_OPENMC:
        config_blockers.append(f"material OpenMC gate is {material_verdict}")
    if binding.lattice_pitch_um is None or binding.density_g_cm3 is None:
        config_blockers.append("voxel binding is not populated")

    sweep_blockers = list(config_blockers) or [
        "no measured Reference D config has been emitted or loaded"]

    return {
        CAMPAIGN_READY: (not campaign_blockers, campaign_blockers),
        CONFIG_READY: (not config_blockers, config_blockers),
        SWEEP_READY: (False, sweep_blockers),
    }


def declared_binding() -> VoxelBinding:
    """The binding as the repository actually declares it, blanks and all."""
    rows = read_table(DATA / "voxel_binding.csv", VOXEL_BINDING)
    r = rows[0]
    return VoxelBinding(
        entity_kind=r["entity_kind"], lineage_semantics=r["lineage_semantics"],
        material_class=r["material_class"],
        material_represents=r["material_represents"],
        material_model_kind=r["material_model_kind"],
        density_basis=r["density_basis"],
        lattice_pitch_um=r["lattice_pitch_um"],
        density_g_cm3=r["density_g_cm3"],
        notes=r.get("notes", ""))


def target_spec() -> ReferenceSystemSpec:
    """Reference D as it stands: identified, and otherwise empty.

    Deliberately constructed from blanks rather than from a stored file. There
    is no measured Reference D config to read, and inventing one to ask what is
    missing would answer the wrong question.
    """
    blank = {f.name: None for f in dataclasses.fields(ReferenceSystemSpec)}
    blank.update(reference_system_id="D_engineered_composite",
                 system_provenance="engineered_composite",
                 evidence_policy="measured_only",
                 execution_class="target_calibration",
                 target_calibration=True)
    return ReferenceSystemSpec(**blank)


def required_fields(binding, spec, gates) -> set[str]:
    """The emitter's own answer, reduced to field names."""
    spatial_verdict, material_verdict = gates
    problems = emit_problems(binding, spec, spatial_verdict=spatial_verdict,
                             material_openmc_verdict=material_verdict,
                             stage="dosimetry")
    fields = set()
    names = [f.name for f in dataclasses.fields(ReferenceSystemSpec)]
    names += ["lattice_pitch_um", "density_g_cm3"]
    for p in problems:
        for name in names:
            if p.startswith(f"{name} is unset") or p.startswith(f"{name} =") \
                    or p.startswith(f"{name} is not"):
                fields.add(_ALIASES.get(name, name))
    return fields


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the ledger and the contract disagree")
    args = ap.parse_args(argv)

    requirements = load_requirements()
    binding = declared_binding()

    # The FROZEN configs, not the templates: the templates ship unset by
    # design and reading them here would report the freeze as never happening.
    spatial = spatial_report.evaluate(
        DATA / "spatial", REPO / "config" / "reference_d_spatial_acceptance.toml")
    material = material_report.evaluate(
        DATA / "materials", REPO / "config" / "reference_d_material_acceptance.toml")

    print("REFERENCE D — engineered radiotropic composite\n")
    print(f"binding coherent : {'yes' if not coherence_problems(binding) else 'NO'}")
    print(f"spatial gate     : {spatial.verdict}")
    print(f"material OpenMC  : {material.openmc}")
    print(f"material radiodx : {material.radiodialysis}")

    fields = required_fields(binding, target_spec(),
                             (spatial.verdict, material.openmc))

    by_source: dict[str, list[dict]] = {k: [] for k in SUPPLIED_BY}
    for r in requirements:
        by_source[r["supplied_by"]].append(r)

    by_status: dict[str, list[str]] = {}
    for r in requirements:
        by_status.setdefault(r["status"], []).append(r["requirement_id"])
    print(f"\n{len(requirements)} requirements, "
          f"{len(by_status.get('satisfied', []))} satisfied")
    print("\nWAITING ON")
    for state in STATUSES:
        ids = by_status.get(state)
        if state != "satisfied" and ids:
            print(f"  {state:<30} {len(ids):>2}  {', '.join(sorted(ids))}")
    print()
    for source in SUPPLIED_BY:
        rows = by_source[source]
        if not rows:
            continue
        print(f"  {source.upper()} ({len(rows)})")
        for r in rows:
            mark = "x" if r["status"] == "satisfied" else " "
            print(f"    [{mark}] {r['requirement_id']:<14} {r['artifact'][:88]}")
        print()

    # THE FREEZE. Every requirement must be justified by something the code
    # actually refuses on, and every such refusal must have a requirement.
    #
    # Two kinds of justification, because the emitter is not the only authority:
    # a row may name a config field the emitter demands, OR it may name a GATE
    # that is not yet READY. `cpm.seconds_per_mcs` is the case that forces the
    # distinction — it is a cpm_feedback key the transport emitter never reads,
    # but the spatial gate will not reach READY_FOR_TIME_CALIBRATION without a
    # selectable dynamic observable, so the measurement is required all the same.
    open_gates = {"spatial gate": spatial.verdict != spatial_report.READY,
                  "material OpenMC gate": material.openmc != material_report.READY_OPENMC,
                  "radiodialysis gate": material.radiodialysis != material_report.READY_RADIODIALYSIS,
                  "both gates": spatial.verdict != spatial_report.READY
                  or material.openmc != material_report.READY_OPENMC}

    ledgered = {r["config_key"] for r in requirements}
    unledgered = sorted(fields - ledgered)
    unjustified = sorted(
        r["requirement_id"] for r in requirements
        if r["status"] != "satisfied"           # a met requirement needs no open gate
        and r["config_key"] not in fields
        and not open_gates.get(r["unblocks"], False))

    ok = True
    if unledgered:
        ok = False
        print("THE CONTRACT REQUIRES FIELDS THIS LEDGER DOES NOT NAME:")
        for f in unledgered:
            print(f"  - {f}")
    if unjustified:
        ok = False
        print("THESE REQUIREMENTS ARE JUSTIFIED BY NOTHING THE CODE REFUSES ON:")
        for rid in unjustified:
            print(f"  - {rid}: neither a required config field nor an open gate")
    if ok:
        print("every requirement is justified by an open gate or a required "
              "field, and every required field has a requirement.")

    blocked = [r for r in requirements
               if r["status"] == "unsupported_by_current_model"]
    if blocked:
        print(f"\n{len(blocked)} requirement(s) await a MODEL CHANGE, not a "
              "measurement — do not put them in a campaign:")
        for r in blocked:
            print(f"  - {r['requirement_id']}: {r['notes'][:110]}")

    enforcement = enforcement_report()
    decorative = enforcement.get("not_implemented", [])
    if decorative:
        print(f"\n{len(decorative)} declared criteria have NO CONSUMER and may "
              "not close a gate:")
        for key in sorted(decorative):
            print(f"  - {key}")

    print("\nREADINESS")
    # `reached` deliberately does NOT shadow `ok`: --check reports whether the
    # ledger and the contract agree, which is a different question from whether
    # the programme is ready. Conflating them made --check fail forever, since
    # BIOFILM_SWEEP_READY is false by construction until measurements exist.
    for verdict, (reached, blockers) in readiness(
            requirements, spatial.verdict, material.openmc, binding).items():
        print(f"  {verdict:<26} {'YES' if reached else 'no'}")
        for b in blockers:
            print(f"      - {b}")

    if args.check:
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
