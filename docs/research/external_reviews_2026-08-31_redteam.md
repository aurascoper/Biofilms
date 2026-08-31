# Red-team: two commissioned external reviews vs. what this repository already holds

2026-08-31 · not a calibration, not a claims-ledger entry, not a dataset candidate.

## What this is

Two long external research reviews were commissioned from outside this repository — one on
symplectic integration and biofilm rheology, one on whether biofilm morphology is
mechanically determined. Both make specific claims about this codebase and issue
recommendations. This document checks each claim against the source and records where it
was already known, where it is sharper here, and where it is wrong.

**This document populates no calibration value and creates no ledger row.** Verdict
vocabulary borrowed from `data/claims_ledger.csv`'s closed set (`keep | restate |
requalify | delete | needs_calibration | needs_verification`, plus `supported` and
`must_not_be_claimed`).

**Two of the reviews' errors originated in the brief, not in the reviews.** The briefs were
written here, and they asserted that §3.10 is marked specification-not-method and that the
model is CPM-with-adhesion-only. Both came back reported as findings. **A premise returned
as a finding is not independent confirmation of anything, and it reads exactly like one** —
which is why it is recorded rather than quietly corrected. The two failures take different
repairs: an *asserted-X* defect is repaired in the review; a *told-X* defect is repaired in
**how the next review is commissioned** — state the premise as a premise to be checked, not
as background. Collapsing them loses the only lesson that transfers outward.

**Everything the reviews say about the outside literature is `needs_verification`.** Purcell,
Shaw's 18-minute relaxation time, Graner–Glazier, Yan 2019, the rheology-transfer analysis,
the modulus ranges: none of it was checked here, both documents self-flag items as resting
on aggregated reference lists, and none of it is adopted below as established.

---

## Claim 1 — §3.4 should be reclassified to state that no momentum state exists

**As stated (review):** "Reclassify the section. State explicitly that … no
independently-integrated momentum state exists or *should* exist."

**What the repo already holds** (`preprint/…tex:527`):

> **This subsection is specification, not method.** No program in this work performs
> symplectic integration. … no simulation reported here holds an independently-integrated
> momentum state … No identifier named `momentum` exists anywhere in this repository's
> simulation sources.

The manuscript further distinguishes `biofilms_3d.R`'s velocity array — recomputed from
position and forces each step, never advanced by its own equation of motion — a nuance the
review does not mention. `tests/manuscript_claims_tests.jl:85-99` enforces the identifier
claim with a planted-positive control.

**Verdict: `keep`.** Already done, in more detail than the recommendation asks for.

---

## Claim 2 — a momentum state would be physically inappropriate, not merely absent

**As stated (review):** Re ≈ 10⁻⁴–10⁻⁵, inertial relaxation ~10⁻⁷ s, ten orders below any
biological timescale, so adding momentum would be *physically wrong*.

**Checked here by computation**, not accepted: Re = 2.0e-5 (1 µm at 20 µm/s) to 2.0e-3
(100 µm feature); τ_p = 5.6e-8 s at a = 0.5 µm; T_bio/τ_p = 1.9e10 against the 18-minute
figure. `analysis/overdamped_regime.py` reproduces all of it.

**Verdict: `supported`, and it is the reviews' one genuine addition.** §3.4 said the absence
is factual; it now says the absence is *appropriate*, which is the stronger claim. Applied
to §3.4 on 2026-08-31 with the inputs stated alongside the numbers, because Re is a
property of the flow and not of the organism.

---

## Claim 3 — the Diele citation is a three-way bibliographic conflation

**As stated (review):** the citation conflates three papers; the title belongs to Diele &
Marangi, *Mathematics* 8(1):25.

**What the repo already held** (`data/claims_ledger.csv`, `PP-REF-02`, verdict `delete`,
prescribing the same replacement):

> EVERY FIELD WRONG. The DOI resolves to a neural-network quadratic-programming paper by
> Yang, Cao & Xu. Different authors, journal, year, volume and pages.

**Verdict: `requalify`.** The defect is real and the prescribed fix matches, but the
review's *diagnosis* is not what is wrong here, and it is the weaker one. The review
analysed title/journal/pages and **did not check the DOI**, which resolves to an unrelated
paper — a fourth error, and the worst of them. It also missed `PP-REF-01` entirely, whose
DOI resolves to an unrelated 3D MHD regularity paper.

**What the review did surface, by rediscovery: both fixes were prescribed and never
applied.** That is the finding, and it is ours rather than theirs.

---

## Claim 4 — §3.10 is marked specification-not-method

**As stated (review):** "Keep §3.10 marked specification-not-method."

**What the source showed:** it was not marked at all. The disclaimer appeared at §3.2,
§3.4 and §3.12 only — corroborated by `data/claims_ledger.csv:498` (`PP-V11-05`), and by
`PP-310-01`, an open row already recommending "Label as a proposed component, not an
implemented one."

**Verdict: `delete`.** **This premise came from the brief.** It was supplied as background,
returned as a finding, and would have read as confirmation that no action was needed. The
disclaimer was added on 2026-08-31.

---

## Claim 5 — the nutrient field is inert with respect to configuration

**As stated (review):** `compute_delta_H_terms` and `mcs_step!` never read `state.nutrient`,
so the field cannot influence configuration.

**Verified**, and already enforced: `tests/manuscript_claims_tests.jl:256-343` greps every
lattice-mutating function for nutrient reads, asserts none, and separately asserts
`update_nutrient!` *does* read it.

**Verdict: `keep`, with one scope correction.** The narrow claim holds. **"The field is
unused" would be false** — it is live, integrated every MCS, consumed, checkpointed, and
GPU-parity-tested. Inert with respect to acceptance is not inert.

---

## Claim 6 — the CPM is adhesion-and-volume only

**As stated (review):** a bare adhesion+volume CPM produces differential-adhesion sorting,
not biofilm morphology.

**What the source showed:** the Hamiltonian carries four terms — `adh`, `vol`, `rad`, `mel`
(`biofilms_potts.jl:590`). The radiation and melanin terms are exactly the field bias the
review's own analysis says would modify sorting.

**Verdict: `requalify`. This premise also came from the brief.** The conclusion — that the
output is directional redistribution rather than a biofilm morphology — is one the
manuscript already states; the route to it was supplied rather than found.

---

## Claim 7 — the §6.3 radial result needs an SE and direction caveat

**Verdict: `keep`.** Already the manuscript's own position, in the abstract (`…tex:112`),
in §6.3, in `README.md:719`, and in ledger rows `RM-KR-07` / `RM-KR-08`.

---

## Claim 8 — scope the CPM-for-biofilms novelty claim

**Verdict: `delete`.** No such claim exists. The novelty claim is narrowly about coupling
OpenMC to a CPM and already self-limits: "This work therefore extends an established
paradigm to a new transport code and a new biological target rather than introducing the
paradigm."

---

## Claim 9 — "no flow" should be a declared limitation

**Verdict: `needs_verification` → open.** The reviews are right. §7.5 declares that
"viscoelastic mechanics are absent from the executed simulations", but flow, advection,
shear and hydrodynamics appear nowhere as model or as declared limitation. SCOPE: searched
`preprint/…tex`, `README.md` and `data/claims_ledger.csv` on 2026-08-31 for those four terms;
not a claim about every file in the repository. **Not fixed in
this pass** and recorded here so it is not lost.

---

## The structural finding, which neither review made

Four ledger verdicts were recorded, prescribed, and never applied: `PP-REF-01`,
`PP-REF-02`, `PP-310-01`, `RM-G04-01`. One cause was assumed; there are **two**, and they
must not be filed as one.

**Cause 1 — the scope of enforcement.** A verdict is enforced only in the file its
`document` column names. `RM-G04-01` was restated in README while
`biofilms_potts.jl:16` went on listing a five-term Hamiltonian, because no guard reads
source comments. Closing this needs a `RETRACTED_IN_SOURCES` scan over `SIM_FILES`.
**Not built. Open.**

**Cause 2 — the phrase extractor.** `PP-REF-01/02` were invisible for an unrelated reason:
a citation splits on commas into runs under `MIN_WORDS = 5`, so `distinguishing_phrase`
returns `''`. **This is a property of the extractor, not of which files are scanned**, and
no widening of the scanned set reaches it. Closed on 2026-08-31 by `RETRACTED_CITATIONS`.

**`RETRACTED_CITATIONS` reaches `PP-REF-01` and `PP-REF-02` and nothing else. `RM-G04-01`
remains unenforced pending `RETRACTED_IN_SOURCES`.** Stated in those terms so a later
reader does not see "the systemic gap has a guard" and stop checking.

A third defect surfaced while applying the fixes: §3.4 referred implementing symplectic
integration to §7.5 as future work, and §7.5 does not list it. The cross-reference was
dangling and is now removed.

### A fourth, and it is a new member of the use-versus-mention family rather than a repeat

The first draft of the reference correction named both withdrawn DOIs **inside LaTeX
comments in the manuscript**. `RETRACTED_CITATIONS` would have caught it. **The design
would not have**, and that gap is the finding.

Use-versus-mention was resolved at the **file** level — which documents may name a
withdrawn DOI. A `.tex` comment defeats that resolution without violating it: it sits
inside a *permitted* file, at a location that is not part of the rendered document. It is
neither use nor mention in the sense the allow-list encodes. It is **invisible to a reader
of the PDF and fully visible to a reader of the source — and the source reader is the one
who copies the string**, which is the whole reason withdrawn identifiers are kept out of
the manuscript.

So a file-level allow-list is not sufficient for artifacts that have a rendered form and a
source form. The repair taken here was to remove the strings rather than to teach the guard
about comments, because the manuscript has no reason to carry a withdrawn DOI in any form;
the comments now point at the ledger row. **Recorded because the next guard of this shape
should decide its scope at the level of *rendered versus source*, not only of file.**

---

## Summary

| # | Claim | Verdict |
|---|---|---|
| 1 | §3.4 should say no momentum state exists | `keep` — already done, in more detail |
| 2 | Momentum would be physically inappropriate | `supported` — the one genuine addition; applied |
| 3 | The Diele citation is a three-way conflation | `requalify` — real defect, weaker diagnosis, missed the DOI and `PP-REF-01` |
| 4 | §3.10 is marked specification-not-method | `delete` — it was not; **premise from the brief** |
| 5 | The nutrient field is inert | `keep` — verified and already tested; "unused" would be false |
| 6 | The CPM is adhesion-and-volume only | `requalify` — four terms; **premise from the brief** |
| 7 | §6.3 needs an SE/direction caveat | `keep` — already the manuscript's position |
| 8 | Scope the CPM-for-biofilms novelty claim | `delete` — no such claim is made |
| 9 | "No flow" is undeclared | open — correct, and not fixed in this pass |

Kept so that a commissioned review's premises are not mistaken later for its findings, and
so the two causes behind four unapplied verdicts stay distinguishable.
