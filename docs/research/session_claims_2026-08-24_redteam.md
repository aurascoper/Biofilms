# Red-team: unsourced session claims vs. the existing evidence audit

2026-08-24 · not a calibration, not a claims-ledger entry, not a dataset candidate.

## What this is

A Claude Code session (outside this repository, in a general-knowledge conversation about
radiotrophic fungi) produced several confident explainers about radiotropism,
radiosensitivity, the "radiosynthesis" mechanism, an ISS shielding experiment, and
bioremediation potential — with no citation, no dosimeter check, no cross-reference to this
repository's own audit. This document runs that session output back through the audit this
repository already did on nearly the same claims
(`docs/research/radiotrophic_compatibility_audit.md`,
`docs/research/radiotrophic_pilot_recommendations.md`, both dated 2026-08-15) and records
where it was wrong, overstated, or simply unevaluated.

**This document populates no calibration value and creates no ledger row.** It is a record
of what an outside conversation got wrong or left unverified about this repository's subject
matter, kept here so the discrepancy doesn't quietly resurface as if it had been checked.

Verdict vocabulary borrowed from `data/claims_ledger.csv`'s closed set (`keep | restate |
requalify | delete | needs_calibration | needs_verification`, plus `supported` and
`must_not_be_claimed`), used here as an assessment label only — none of these rows exist in
that file, because `claims_ledger.csv` is scoped to claims already asserted in this
repository's own carrier documents (README, preprint), and a general conversation session
is not one.

---

## Claim 1 — ISS fungal layer "blocked about 2% of incoming cosmic radiation"

**As stated (session):** "a layer of the fungus just 1.7 millimeters thick successfully
blocked about 2% of incoming cosmic radiation."

**What the repo's audit found**
(`docs/research/radiotrophic_pilot_recommendations.md` §1.2, item 2 — "Not
`radiation_shielding`"):

> "Withdrawn by audit. 147 vs 151 CPM (2.6%), n = 1 dish, one sensor pair, phase-3 p = 0.069
> in SI-1 model 3.1 against a proper phase-1 negative control at p = 0.970. Counts dropped;
> the effect was not established."

**Verdict: `must_not_be_claimed`.** A p = 0.069 result at n = 1, already withdrawn by this
repository's own audit of the primary source, was presented in-session as an established,
successful shielding result. The 2.6% figure in the audit and the "about 2%" figure in the
session claim are the same underlying number; the audit's point is precisely that it never
cleared statistical significance and the counts were later dropped from the analysis.

## Claim 2 — "a layer of the fungus just 1.7 millimeters thick"

**As stated (session):** presented as a measured lawn thickness.

**What the repo's audit found** (`radiotrophic_pilot_recommendations.md` §1, preamble):

> "The `~1.67 mm` figure that circulates as a lawn thickness is `15 mm − 20 mL/75 cm²`, a
> geometric residual explicitly described in the supplement as 'for the fungal growth layer
> AND/OR GASEOUS HEADSPACE' — it may be entirely air."

**Verdict: `must_not_be_claimed`.** Not a measurement of fungal biomass thickness at all —
a subtraction with an explicit caveat, in the source's own supplement, that it may describe
headspace rather than fungus.

## Claim 3 — the melanin → NAD⁺ → NADH → ATP "radiosynthesis" mechanism

**As stated (session):** presented as a mapped, if not fully complete, metabolic pathway
("here is exactly how it works"), attributed to Albert Einstein College of Medicine
research.

**What the repo's audit found**
(`docs/research/radiotrophic_compatibility_audit.md`, headline verdict):

> "radiotrophy is not established for any of the seven species. The claim rests on one 2007
> paper whose exposure was characterised by analogy rather than by a stated dosimeter, whose
> growth effects are single-digit percentages, and whose integrated cell dose is not stated.
> Against it stand a same-species independent report of no substantial melanin protection
> under gamma, a same-species negative melanin-pathway result at 1 and 3 kGy, a same-genus
> null under a TLD-anchored ¹³⁷Cs field, an isogenic null in *A. niger* across three ionizing
> modalities, and a formal energy-budget critique. Peer review struck the word from the
> flagship ISS title. No author contact resolves this, because the discriminating experiment
> has not been run by anyone."

**Verdict: `requalify`.** The electron-transfer mechanism is real published work (the same
2007 paper this repo's audit already scoped), but presenting it as "exactly how it works" —
a settled pathway — omits that (a) the paper the mechanism rests on lacks a stated dosimeter
and integrated cell dose, and (b) multiple independent same-species and same-genus follow-up
studies found no melanin-pathway effect, to the point that peer review removed the word
"radiotrophic" from the flagship ISS paper's own title. This should be stated as one
contested mechanistic hypothesis, not as established biochemistry.

## Claim 4 — "grew significantly faster and accumulated more biomass"

**As stated (session):** a single, unqualified growth-rate claim.

**What the repo's audit found** (`radiotrophic_pilot_recommendations.md` §1.2, item 15):

> "Not a single headline growth ratio. The abstract says 1.21 ± 0.37-fold, SI-1 §D says
> 1.64-fold ± 0.279 SE on a logistic-slope basis, and the article's own exponential k values
> give 1.24 / 1.29 / 1.53. The basis of 1.21 ± 0.37 is undocumented. Record the disagreement;
> do not average, do not choose."

**Verdict: `needs_verification`.** Directionally consistent with the source, but the
session claim states a single clean effect where the source itself reports three
inconsistent figures depending on which section is read, with an undocumented headline
number. The repo's explicit instruction is to record that disagreement rather than round it
away.

## Claim 5 — Chernobyl-graphite bioremediation, and MURR facility specifics

**As stated (session):** fungi "found growing on and slowly breaking down highly radioactive
graphite" inside the Chernobyl reactor; MURR irradiation chambers, beamports, pneumatic
tubes, hot cells, and isotope-production capability as a research venue.

**Status: `not_evaluated_by_repo`.** A repo-wide search —
`grep -rin "MURR\|Missouri Research Reactor\|reactor pool\|hot cell\|beamport\|Lu-177\|
Lutetium\|isotope production\|graphite" data/ docs/` — returns no prior hits for the MURR
terms in this repository (see `docs/research/murr_facility_candidate.md` for the forward
note on that). The graphite-biodegradation claim is likewise not addressed anywhere in
`radiotrophic_compatibility_audit.md`. Neither claim is contradicted by anything in the
repo; neither is supported by anything in the repo either. Treat both as unverified,
outside-the-audit popular-science claims until they get their own sourced pass.

## Claim 6 — *D. radiodurans* as the reputational example of extreme radioresistance

**As stated (session):** "The bacterium *Deinococcus radiodurans* is famous for its extreme
radioresistance. It can survive radiation doses thousands of times higher than what would be
lethal to a human because its DNA repair mechanisms are incredibly fast and efficient" —
stated as established fact, no citation, no dose figures, no source.

**What the repo's audit found**
(`docs/research/radiotrophic_compatibility_audit.md`, "Deinococcus radiodurans" §, and §12
"Integrated"):

> "No record in this audit applied ionizing radiation to *D. radiodurans*. The species'
> radioresistance is established elsewhere in the literature; it is not established by
> anything catalogued here, and the label is not imported on reputation."

> "An exhaustive sweep for volumetric *D. radiodurans* biofilm data returned nothing across
> Dryad, Zenodo, BioImage Archive, BioStudies, EMPIAR, IDR, OSF, Figshare and PubMed. The
> *D. radiodurans* arm of any spatial campaign has no public fallback and must be acquired."

> "Two of the seven, *D. radiodurans* and *B. subtilis*, have no ionizing-radiation record in
> this audit at all."

**Verdict: `needs_verification`, not `must_not_be_claimed`.** This is a different failure
mode from Claims 1–3: it isn't contradicted, and it's very likely true in the broader
literature — *D. radiodurans*'s radioresistance is one of the best-known facts in radiation
microbiology. But the session stated it with reputational confidence and zero citation, and
this repository's own eleven-repository, per-species sweep — the most careful look this
specific claim has gotten *in this project* — came up with **nothing**: no dose curve, no
D10/D37 figure, no irradiated-vs-control contrast, nothing volumetric, for the exact species
the model carries as `DR` in its seven-species registry
(`data/calibration/spatial/entity_semantics.csv`). The audit's own discipline is the
instructive part: it refuses to import a true-sounding, well-known claim on reputation alone,
and neither should this document. If *D. radiodurans* radioresistance is ever asserted in
this repository as a source-backed fact rather than background color, it needs its own
citation — the audit found none to reuse.

---

## Summary

| # | Claim | Verdict |
|---|---|---|
| 1 | ISS shielding "successfully blocked ~2%" | `must_not_be_claimed` — audit withdrew it (p = 0.069, n = 1) |
| 2 | 1.7 mm measured lawn thickness | `must_not_be_claimed` — a geometric residual, may be headspace |
| 3 | Radiosynthesis as settled mechanism | `requalify` — one contested 2007 paper, several independent nulls |
| 4 | Single clean growth-rate figure | `needs_verification` — source reports three inconsistent ratios |
| 5 | Graphite bioremediation / MURR specifics | `not_evaluated_by_repo` — genuinely outside existing audit scope |
| 6 | *D. radiodurans* radioresistance by reputation | `needs_verification` — likely true, zero record in this repo's own sweep |

Nothing in this table is a calibration input, a dataset candidate, or a claims-ledger row.
It exists so the next time any of these six claims surfaces — in this repo, in the
career-ops interview prep, or in conversation with the MURR contact Caixia Wan is arranging
— it can be answered against sourced text instead of memory.
