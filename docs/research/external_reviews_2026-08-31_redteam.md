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

**Verdict: `needs_verification` → CLOSED 2026-08-31, and the framing changed on the way.**
The reviews are right that flow is undeclared, but adding "no flow" beside the existing
entry would have overstated its independence. Flow reaches morphology through mechanics, and
a model with no growth generates no growth-induced stress for mechanics to act on — so the
three are one structural absence under three names. §7.5 now says that, and says what it
costs: no route by which morphology could arise, so the spatial output is directional
redistribution in an imposed field rather than a morphology.

**The literature half stays unverified and out of the manuscript.** The reviews argue
growth-induced stress is the dominant still-condition morphogen (Yan 2019, Asally 2012).
That is not checked here, §7.5 says so explicitly, and nothing in the added sentence depends
on it. The claim carrying the weight is internal: this model has no growth, no mechanics and
no flow.

The original finding, for the record. §7.5 declared that
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

**Cause 1 — the scope of enforcement. CLOSED 2026-08-31, and building the guard split the
cause in two.** (a) The `document` column names *another* file: `RM-G04-01` was restated in
README while `biofilms_potts.jl:16` went on listing a five-term Hamiltonian. (b) Worse, and
not seen when this was first written: the column names *no* file. `document = repository` is
a `PSEUDO_DOCUMENT`, so `document_path` returns `None` and the row is never searched at all.
Measured: 37 rows sit in that class, 5 carrying an unresolved verdict that names a source.
`PP-DEFF-01` was one, and its prescribed comment fix sat unapplied at
`biofilms_radiodialysis.R:242` — **the fifth instance of the unapplied-verdict pattern,
found by building the guard for the first four.**

`RETRACTED_IN_SOURCES` is a declared vocabulary rather than a phrase sweep, and that was
measured before it was chosen: `distinguishing_phrase` over every `delete` row against all
144 source files returns **zero** hits, because claims living in source comments are short —
a formula, a range, a symbol — and split into runs under `MIN_WORDS`, exactly as figure
labels and citations do. Third vocabulary tier, same reason as the first two. `restate`
verdicts are not swept generally: AGENTS.md is explicit that a restate is a judgement about
meaning and stays with a reviewer; only the specific withdrawn string is mechanical.

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

# Third review — melanin radiotrophy

A third commissioned review arrived the same day, on whether melanin converts ionising
radiation into biologically usable energy. **The pattern holds for a third time:** almost
everything substantive in it was already in `docs/research/radiotrophic_compatibility_audit.md`
(2026-08-15) and carried into §2.6 before the review was commissioned.

## Claim 10 — the manuscript should not claim more than "not demonstrated"

**What the repo already holds** (`preprint/…tex:325`):

> "``Not demonstrated'' is not ``disproved'', and we state it as the former."

**Verdict: `keep`.** The manuscript is already more careful than the recommendation, and
draws a distinction the review does not.

## Claim 11 — no growth means radiotrophy cannot be tested

**Verdict: `keep`.** Stated in four places — abstract `:121`, §2.6 `:315`, Discussion
`:1244`, Conclusion `:1350` — and independently in `README.md:49-63`.

## Claim 12 — §2.6 draws the radiolytic-H₂ contrast

**Verdict: `keep`**, `:336-361`, and the manuscript already labels its own limits:
"This is a contrast and not a defence: it does not address the carbon-budget critique."

## Claim 13 — the energy budget, and the antioxidant alternative

**Verdict: `keep` on both.** The energy-budget critique is Walberg 2015, quoted verbatim in
the audit at `:669-677`. The antioxidant/repair account is §2.1 `:210-227`, with Cortesão's
NHEJ series.

**One arithmetic check of ours, stated rather than quietly dropped.** The review's budget
reproduces exactly on two of three figures: 5e-18 J/h deposited into a 1e-13 kg cell at
0.05 mGy/h, and ~60 ATP-equivalents/cell/h at perfect conversion. Its third — "~500 Gy/h to
supply 1 % of maintenance" — **recomputes here as ~96 Gy/h**, a factor of about 5. Working:
1 % of 1 mmol ATP/g/h at 2e-11 g dry mass is 1.2e8 ATP/cell/h, ×8e-20 J = 9.6e-12 J/h, ÷
1e-13 kg = 96 Gy/h. **Differs by 5×; not repeated.** The conclusion is unaffected — both
figures sit orders above any survivable chronic dose rate — which is exactly why the
disagreement is safe to state instead of omitting the number and leaving the next reader to
rediscover it.

## Claim 14 — "your own reference [36]"

**As stated:** "Your own reference [36] is the closest thing to evidence for this — Turick
et al. 2011."

**Verdict: `delete`.** `turick2011` is numerically the 36th `\bibitem`, but **it had no
`\cite` in the body at all.** It was one of fifteen defined-but-uncited entries whose count
`tests/manuscript_claims_tests.jl` already pinned and printed on every run. The review cites
it as something the manuscript relies on; the repository had already classified that exact
entry as unused.

**It is cited now, and for the review's own reason rather than to close the gap.** Turick
2011 is the one positive measurement in a section otherwise built on nulls, and it runs
against the mechanism: melanin "is continuously oxidized in the presence of gamma
radiation." SCOPE: verified from the abstract, cross-indexed on PubMed 21632287, Johns
Hopkins Pure and Semantic Scholar. The review's account of the apparatus — melanin ghosts,
carbon paste, no living cells, nanoampere currents — comes from a full text paywalled from
here and is **not** repeated in §2.6.

`casadevall2017` was cited in the same edit, verified by extracting the ASM PDF and reading
it rather than trusting the review's quote.

## Claim 15 — no fungus has a CO₂-fixation pathway

**Verdict: `needs_verification`, and it is the review's one real contribution.** Searched
the manuscript, `README.md`, `data/claims_ledger.csv` and all four
`docs/research/radiotrophic_*` files: no occurrence of Calvin-Benson-Bassham,
Wood-Ljungdahl, reductive TCA, or autotrophy. **Genuinely absent from this repository.**

It is also the sharper form of an argument §2.6 already gestures at. The manuscript names
the carbon-budget critique as unaddressed, but in Walberg's form — the medium supplied
enough carbon — which is a confounder. "No fungus has a fixation pathway" is structural: it
does not depend on any energetic estimate, and it would not be closed by a better
experiment.

**Not adopted, and not written into the manuscript.** It is a claim about fungal genomics,
outside anything this repository can check, and the review self-flags three of its citations
plus the hydrated-electron potential as unverified against primary sources. Given that two
entries in this bibliography have already turned out to be DOIs resolving to unrelated
papers, nothing from that report enters the manuscript without a primary read.

---

# What three reviews establish about commissioning, which is the durable finding

Three commissioned reviews, on three unrelated subjects — symplectic integration and
rheology, morphology and mechanics, melanin radiotrophy. **All three substantially
rediscovered work this repository already held**, and in the radiotrophy case rediscovered
a single document, `radiotrophic_compatibility_audit.md` of 2026-08-15, whose conclusions
were already carried into §2.6 and the ledger before the review was commissioned.

Across all three, the genuine additions were:

1. the overdamped physical argument (Re, τ_p) — applied to §3.4;
2. that no fungus has a CO₂-fixation pathway — recorded `needs_verification`, **not
   adopted**, because it is fungal genomics and nothing here can check it.

**That is the finding, and it is about commissioning rather than about any review.** A
brief written from a repository that already holds the answer produces a report that
returns the answer. Two of the false claims in the first two reviews came from premises
supplied in the briefs and came back reported as findings — a premise returned as a finding
reads exactly like independent confirmation, which is what makes it worth recording rather
than quietly correcting.

And the one thing all three added was **in a domain this repository has no purchase on**:
low-Reynolds hydrodynamics, and fungal carbon metabolism. Neither is checkable from here;
both are checkable by someone who knows the field.

**So the shape of what external review is good for here is now visible, and so is the shape
of what to ask for.** Not "audit our claims" — the audits exist, and a brief describing them
will get them back. Ask instead for the thing the repository structurally cannot do: name
the domain, state the premise as a premise to be checked rather than as background, and
expect the answer to arrive as `needs_verification` because verifying it is outside what
this repository can reach. A review whose findings this repo *could* have produced is a
review that was told what to find.

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
| 9 | "No flow" is undeclared | **closed 2026-08-31** — declared in §7.5, but as a *relation* rather than a list item: flow is downstream of mechanics, mechanics downstream of growth |
| 10 | Do not claim more than "not demonstrated" | `keep` — the manuscript is already more careful |
| 11 | No growth ⇒ radiotrophy untestable | `keep` — four places plus README |
| 12 | §2.6 draws the radiolytic-H₂ contrast | `keep` — and already labels its own limits |
| 13 | Energy budget; antioxidant alternative | `keep` — both predate the review; one figure recomputes 5× smaller |
| 14 | "Your own reference [36]" | `delete` — it was an uncited orphan the repo already pinned; now cited, for the review's reason |
| 15 | No fungus fixes carbon | `needs_verification` — genuinely absent here, and the review's one real contribution |

Kept so that a commissioned review's premises are not mistaken later for its findings, and
so the two causes behind four unapplied verdicts stay distinguishable.
