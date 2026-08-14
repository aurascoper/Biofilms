"""Screened public dataset candidates, including the rejections.

A search that only records what it accepted is not a provenance trail. The
rejected candidates are the more informative half: they say which modalities
and scales were considered and found incompatible, so the same ground is not
re-searched and a later reader can challenge the judgement.

SURROGATE STATUS IS ENFORCED, NOT CONVENTIONAL. A dataset whose organisms are
not the target consortium cannot clear the target gate however good its data
are, so `target_relationship` and `clears_target_gate` are separate columns and
the second may never be true unless the first is `exact_species` AND the
community matches.
"""

from __future__ import annotations

from ..schema import EVIDENCE_BASIS, STATUS, Column, TableSchema

TARGET_RELATIONSHIP = frozenset({
    "target_consortium",   # the seven-species system itself
    "exact_species",       # one named component species
    "surrogate",           # a different organism, used to exercise machinery
    "unrelated",
})
VERDICT = frozenset({
    "use_first",           # best ratio of size, quality and metadata
    "inspect_first",       # promising, but check the manifest before downloading
    "diagnostic_only",     # useful for a different question than pitch selection
    "rejected",
})
BOOLEAN = frozenset({"true", "false", "unverified"})

DATASET_CANDIDATES = TableSchema(
    name="dataset_candidates",
    doc="public datasets screened for the spatial calibration pilot",
    comments=(
        "Includes REJECTED candidates on purpose: a search that records only",
        "its accepts is not a provenance trail.",
        "size_bytes and file names were verified against the repository API",
        "where possible; anything unverified says so rather than being",
        "transcribed as fact.",
        "NOTE on Dryad: a dataset's storageSize is CUMULATIVE ACROSS VERSIONS.",
        "The current version's file total can be several times smaller, and",
        "budgeting a download from storageSize will badly overestimate it.",
        "clears_target_gate can never be true for a surrogate, however good",
        "the data are — the pitch and density must refer to the target",
        "consortium under its declared conditions.",
    ),
    key=("candidate_id",),
    columns=(
        Column("candidate_id"),
        Column("repository"),
        Column("identifier", doc="DOI or stable URL"),
        Column("title"),
        Column("organism"),
        Column("target_relationship", vocabulary=TARGET_RELATIONSHIP),
        Column("modality"),
        Column("size_bytes", numeric=True, required=False),
        Column("size_basis", required=False,
               doc="current_version | all_versions | unverified"),
        Column("n_files", numeric=True, required=False),
        Column("license", required=False),
        Column("has_3d_stacks", vocabulary=BOOLEAN),
        Column("has_voxel_calibration", vocabulary=BOOLEAN),
        Column("has_time_course", vocabulary=BOOLEAN),
        Column("verdict", vocabulary=VERDICT),
        Column("clears_target_gate", vocabulary=BOOLEAN),
        Column("rationale"),
        Column("verified_via", required=False),
        Column("accessed_date", required=False),
        Column("evidence_basis", vocabulary=EVIDENCE_BASIS),
        Column("source_id", required=False),
        Column("status", vocabulary=STATUS),
    ),
)


def problems(rows) -> list[str]:
    """Contradictions a candidate ledger can contain."""
    out: list[str] = []
    for r in rows:
        cid = r.get("candidate_id", "?")
        if r.get("clears_target_gate") == "true" and \
                r.get("target_relationship") != "target_consortium":
            out.append(
                f"{cid}: clears_target_gate=true but target_relationship="
                f"{r.get('target_relationship')!r} — only the target consortium "
                "under its declared conditions can clear the target gate")
        if r.get("verdict") == "rejected" and not r.get("rationale", "").strip():
            out.append(f"{cid}: rejected with no recorded rationale")
        if r.get("verdict") in ("use_first", "inspect_first") and \
                r.get("has_3d_stacks") == "false":
            out.append(
                f"{cid}: selected for the spatial pilot but has_3d_stacks=false")
    return out


def usable_for_pitch(rows) -> list[dict]:
    """Candidates that could contribute to a pitch, target or surrogate."""
    return [r for r in rows
            if r.get("verdict") in ("use_first", "inspect_first")
            and r.get("has_3d_stacks") in ("true", "unverified")]
