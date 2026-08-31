# The spatial dose–community program — a design, a gate, and three places it can fail

2026-08-30 · forward-looking · **no `source_id` claimed for anything below.** This document proposes
work. It is not evidence, nothing in it may be quoted elsewhere in this repository, and nothing here
may be promoted into `data/claims_ledger.csv` or any `data/research/*_candidates.csv` row until a
citable source exists and is pinned the way every source in `data/sources.csv` is. It follows the
posture of `docs/research/murr_facility_candidate.md`, which records a candidate lead honestly
rather than treating it as more solid than it is.

Items below are marked **verified** where a primary was read on the stated date, and **unverified**
otherwise. Per this repository's boolean convention, unverified is never read as false.

## Why this exists

`preprint/hoffman_memo.tex` asks questions and does not state a program. This is the program behind
them. `docs/research/radiotrophic_lab_gap.md` specifies a different experiment — the smallest
irradiation public data cannot replace, gated on an `institutional_approval_id`.
`docs/research/radiotrophic_pilot_recommendations.md` carries Pilots A and B. Section 7.3 of the
manuscript carries the facility-regime argument. None of them makes the spatial-field claim below.

## Scope: wet storage, essentially only

In-core is flowing borated coolant at temperature in a large flux and nothing lives there, so MIC is
not an operational concern during irradiation. Dry cask storage is helium-backfilled and hot.
That leaves the pool: the one stage where an assembly sits in water cool enough to host anything,
for years to decades. *(In-core and dry-cask conditions: **unverified**, stated from general
knowledge and not from a source read for this document.)*

Biofilms are documented there. The Sarró coupon study and the Karley spent-fuel-pool isolates are
**verified** and carry ledger rows (`HOFFMAN-01`, `HOFFMAN-02`).

## The precedent, and the debt it creates

Microbiology is already formally inside a nuclear materials safety case. SKB's Copper Sulfide Model
(TR-18-08) is a reactive transport model extended to account for microbial reduction of sulfate and
its effect on copper canister corrosion; sulfate-reducing bacteria are placed at the rock–clay
interface using dissolved organic matter and sulfate to generate sulfide, and physical, chemical,
microbial and mass-transport processes are coupled to the interfacial electrochemical reactions.
Best-estimate general corrosion is under 10 µm at one million years. **Verified 2026-08-30**
against SKB TR-18-08 and King et al., *Materials and Corrosion* (2021).

**The counter, stated here rather than left for a committee.** Microbial activity entered that
safety case because sulfide attacks copper through a specific, identified mechanism coupled to an
electrochemical model — not because biology was added as a general term. The precedent therefore
*demands* an equivalent mechanism for FeCrAl, and the memo's own position is that no such mechanism
is established. This is a debt the program takes on. It is not a warrant it inherits.

## What the model does today, and why it is not the test

The implemented phenotype is radiotropism: directed redistribution along a gradient. It predicts
colonisation *toward* the source, which inverts the usual assumption that proximity to fuel is
protective. That is a claim about **total biomass near the fuel**, which is exactly the monotonic
variable the falsification test below declares uninformative. It is recorded here as current
behaviour and is deliberately not part of the argument's spine.

Radiotropism is also not radiotrophy. Section 2.6 of the manuscript states that radiotrophy is not
established for any of the seven species modelled.

## What this model cannot do

Section 2.6 says it in the paper's own words: no dose-dependent survival process and no growth, so
**radiosensitivity cannot be expressed by this model at all.** The pool application is a motivation
for building that. It is not a capability on offer, and any sentence describing a comparison between
a dose-coupled model and a null is describing two models that do not exist.

## The novelty claim, scoped to the corpus actually read

Marciales et al., "Mechanistic microbiologically influenced corrosion modeling — a review,"
*Corrosion Science* (2018), reports that most mechanistic MIC models take SRB as the main player,
that few consider biofilm-specific environmental conditions, that none correlates sessile and
planktonic populations, and that most lack lab or field validation. Its classification is organised
around MIC mechanism types.

**Scope, and it is narrower than the claim wants to be.** The full text returns HTTP 403 from here,
so what was read is the review's abstract and reported findings, not its complete taxonomy. No
dose-coupled or radiation-coupled category appears in what was read. That is *source inaccessible,
not category absent* — the same state `HOFFMAN-11` records for NACE-2019-12944, and it changes the
day someone with access reads the category list. A second review covering radiation or repository
microbiology would be needed to extend the claim across that literature; none has been read, and one
review bounds one field.

## Phase 0 — the pre-check that gates the handover

**0-pre. Derive the noise floor first, as a curve. Nothing here is runnable until this is done.**
The stop condition compares a residual against a compositional noise floor. A floor that does not
exist makes the gate undecidable: Phase 0 could run in full, with pinned geometry and a passing
control, and return no verdict. The floor needs no pool — but it must be a **sensitivity curve and
not a number**, because the floor moves with the assumed effect size and the effect size is the
quantity Phase 1 exists to estimate. Assuming it large makes the design feasible by construction.
Take replicate count and read depth from a published 16S study of a low-biomass surface community,
report detectable residual **across a range of plausible effect sizes**, and state which range the
design survives. The conclusion is then "feasible if the compositional shift is at least X," which
someone with pool access can evaluate against what they already know.

**0a. Inputs pinned, not chosen.** Rack pitch, array size, burnup spread and cooling times all move
the residual, so a design failing at one setting can be rescued by widening the spread or shrinking
the array. Geometry and the burnup/cooling distribution come from a published rack design or a
stated source, are fixed before the run, and are recorded with the result. A second configuration
may be reported, never substituted.

**0b. The statistic is 0-pre's floor**, not a correlation. A threshold on `r` would pass 0.94 and
fail 0.96 on nothing.

**0c. Covariates enumerated now**, because choosing them after seeing the residual is the freedom 0a
exists to remove: distance to the nearest assembly; depth below the water surface; distance to the
nearest pool wall; distance to the nearest rack-array edge; height above the floor. Additions are
recorded as additions.

**0d. The control, because Phase 0 cannot check itself.** Dose at a rack surface is dominated by the
nearest assembly's gamma field, so the computed structure resembles inverse-square in distance
almost by construction and a small residual is expected physics rather than a finding. Run a
configuration where dissociation **must** occur — an extreme burnup contrast, or one assembly
position left empty:

| | control dissociates | control stays flat |
|---|---|---|
| **residual small** | the design does not work; say so | instrumentation is wrong; conclude nothing about the design |
| **residual large** | proceed, conditional on geometry class | **the residual reports discriminating power while the field fails to respond to a configuration that must change it. Its source is a tally, geometry or mesh artifact, and it reads as success. Nothing proceeds.** |

The bottom-right cell is the one a control that only guards nulls would miss.

**0e. Phase 0 gates the handover, not the science.** 0a pins a published rack design; Phase 1 would
sample a different pool with different burnups and cooling times, so Phase 0 establishes feasibility
**in a geometry nobody will sample**, and recomputing for the actual pool is Phase 1's first step
rather than a refinement. It is still worth running, for three reasons needing no access: it debugs
the transport pipeline, it exercises 0d's control, and it produces 0-pre's curve. **What Phase 0
decides is whether this document should be handed to anyone.** A failure means nothing is sent.

## Phase 1 — the empirical study, designed to be handed over

**Dose must vary within a substrate class.** Several positions on the same rack material, same
installation date, same cleaning history, at different distances. Sampling racks against a distant
liner cannot attribute anything to dose: that varies substrate, distance, hydrodynamics, surface
finish, cleaning history and installation date at once, and composition would differ under the null
too. That design is a correlational study reported as a falsification test.

**The discriminating sample is where computed dose and geometric distance disagree**, since within a
rack array they are otherwise near-collinear. This is where the coupling does work a distance
covariate could not: OpenMC makes dose computed rather than a surrogate. Whether any other role for
the coupling would also be load-bearing is not a question this document has surveyed.
Whether such positions exist at usable magnitude is exactly what Phase 0 asks, and there is reason
to expect them to be scarce.

**Temperature is a competing term, not a caveat.** It is a strong driver of community composition
and plausibly dominates any dose signal. If it can be measured at each position it enters the
comparison as a rival explanation; otherwise the null this tests against is weaker than the one a
reviewer would insist on, and the document says so. Flow is named the same way and is harder to
measure.

**The null is the primary accepted outcome.** Individual radiosensitivity is measured directly —
$D_{10}$ values for several of the seven species are real measurements, per Table 2 of the
manuscript — so the untested claim is not organism sensitivity but community composition across a
spatial gradient. If composition does not track computed dose where it diverges from distance, the
coupling earns nothing, and this program commits in advance to accepting that.

Executing needs an operating pool nobody here can reach.

## Phase 2 — build the survival process, which is the thesis

Contingent on Phase 1. Building dose-dependent survival is what would let the model predict **which
taxa shift and in which direction** — a prediction rather than a correlation. The comparison is
dose-dependent survival plus tropism against a null with no dose term, which is a model comparison
against a stated alternative rather than a methods label. **Both halves of that comparison are
unbuilt today.**

## Three places this fails before anyone spends money

1. **0-pre**: no range of plausible effect sizes leaves a detectable residual. Nothing is runnable.
2. **Phase 0**: the residual is below the floor with the control dissociating. Nothing is sent.
3. **Phase 1**: composition does not track computed dose. The coupling earns nothing.

Each has its disposition written before the run. That is the program's most defensible feature.
