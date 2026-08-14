"""The dynamic-observable contract, which gates seconds/MCS.

A selected lattice pitch is NECESSARY for the time calibration:

    dt_MCS = a^2 * S_sim / S_exp

but it is not sufficient. S_exp has to be measured on something the experiment
and the model BOTH represent, and under parcel semantics that is not automatic.
Ordinary single-cell tracking gives a mean-squared displacement for an
individual organism; a computational biomass parcel is not an observable
bacterium, so a per-cell MSD is not directly a parcel MSD.

So the dependency is:

    pitch selected -> dynamic observable declared -> seconds/MCS

not pitch -> time.

WHAT THE MODEL CAN PRODUCE. compute_delta_H has four terms — adhesion, volume,
radiation, melanin — so the CPM has motility and shape rearrangement driven by
surface minimisation, and nothing else. It has no nutrient-driven biomass
creation: `state.nutrient` is written but never read by the acceptance test,
divide_cell! has no trigger, and there is no growth rate anywhere. Any
observable whose experimental counterpart is dominated by GROWTH is therefore
not a candidate, however well it can be measured.
"""

from __future__ import annotations

from ..schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

BOOLEAN = frozenset({"true", "false"})

TIME_OBSERVABLE = TableSchema(
    name="time_observable",
    doc="candidate dynamic observables for the seconds/MCS calibration",
    comments=(
        "represented_by_model is judged against the FOUR terms actually in",
        "compute_delta_H (adhesion, volume, radiation, melanin). The CPM has",
        "motility and shape rearrangement but no nutrient-driven biomass",
        "creation, so growth-dominated observables are excluded regardless of",
        "how well they can be measured.",
        "semantically_compatible_with_parcel is the separate question the",
        "entity semantics force: a per-cell MSD is not a parcel MSD.",
        "An observable must satisfy ALL THREE to be selectable.",
    ),
    key=("observable_id",),
    columns=(
        Column("observable_id"),
        Column("description"),
        Column("model_quantity", doc="what the CPM would produce"),
        Column("experimental_quantity", doc="what would be measured"),
        Column("represented_by_model", vocabulary=BOOLEAN),
        Column("measurable", vocabulary=BOOLEAN),
        Column("semantically_compatible_with_parcel", vocabulary=BOOLEAN),
        Column("selected", vocabulary=BOOLEAN),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id", required=False),
        Column("status", vocabulary=STATUS),
        Column("notes", required=False),
    ),
)


def selectable(row: dict) -> bool:
    """All three conditions, not any of them."""
    return (row.get("represented_by_model") == "true"
            and row.get("measurable") == "true"
            and row.get("semantically_compatible_with_parcel") == "true")


def problems(rows) -> list[str]:
    """Why no observable is selected yet, or why a selection is invalid."""
    out: list[str] = []
    if not rows:
        return ["no candidate dynamic observables declared"]
    chosen = [r for r in rows if r.get("selected") == "true"]
    if not chosen:
        n_ok = sum(1 for r in rows if selectable(r))
        out.append(
            f"no dynamic observable selected ({len(rows)} candidates, "
            f"{n_ok} satisfy all three conditions) — seconds/MCS cannot be "
            "calibrated against an observable the model does not represent")
    if len(chosen) > 1:
        out.append(
            f"{len(chosen)} observables marked selected; exactly one defines "
            "S_exp")
    for r in chosen:
        if not selectable(r):
            missing = [k for k in ("represented_by_model", "measurable",
                                   "semantically_compatible_with_parcel")
                       if r.get(k) != "true"]
            out.append(
                f"observable {r['observable_id']!r} is selected but fails: "
                + ", ".join(missing))
    return out
