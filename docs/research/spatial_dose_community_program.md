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

### 0-pre, first pass — four searches on 2026-08-30, a candidate answer and not a closed one

This is not a systematic read. **It does not close the gate; it moves the failure mode.**

**The plausible effect size spans about eleven orders of magnitude, and the program must pick an
end.** Two anchors, both verified:

- **Chronic, decade-scale selection.** A Th-232-contaminated site near Chengdu, irradiated over ten
  years, four sample groups spanning **193 ± 5 to 911 ± 41 nGy/h**. PCA clusters High and Medium
  fungi away from Low and Blank, putting a threshold near **480 nGy/h**; the effect on fungal
  diversity exceeds that on bacterial, and bacterial α-diversity does not differ significantly
  across groups. **Verified 2026-08-30.**
- **Within-experiment population response.** Shuryak et al., *PLoS One* (2017): effects reported at
  **36–126 Gy/h**, with 78 of 145 phylogenetically diverse fungal strains growing at 36 Gy/h.
  **Verified 2026-08-30.**

Those measure different things — selection integrated over generations against population density
inside an experiment — and a spent fuel pool sits on the first mechanism. **So the program must
declare which mechanism it predicts, because the expected residual differs by orders of magnitude
depending on the answer, and a prediction that accommodates both predicts nothing.**

**The pool sits between the two anchors, which is the awkward part.** A dose rate of **0.416 Gy/h**
is reported for areas holding irradiated fuel elements — roughly 10⁶ above the Chengdu threshold and
about 10² below Shuryak's tested rates. Neither anchor transfers directly, and no prior effect size
from the actual system was found — see the scoped absence under the novelty claim below, which names
the search and its date.

**One finding runs in favour of feasibility.** Chengdu detected a composition threshold across a
total dose range of about **4.7×**, not orders of magnitude. A detectable effect did not require an
enormous gradient in that system.

**And one runs against it, as a mechanism rather than a power problem.** Shuryak reports that
resistance is **cell-concentration dependent** — *D. radiodurans* grew at 67 Gy/h at high density and
growth was extinguished at a tenfold lower concentration, a roughly threefold shift in critical dose
rate — and that resistant cells **protect neighbours across phyla**: wild-type *D. radiodurans*
enhanced *E. coli* survival by more than 200-fold, apparently through catalase-mediated detoxification
of reactive oxygen species in the shared medium. If composition is buffered by density-dependent
cross-protection, it is **not decomposable into per-species radiosensitivity**, and "composition
tracks computed dose" has no monotonic form to test. That belongs here rather than being discovered
at Phase 1.

**It also happens to be the strongest available argument for a spatial model, and the document should
not hide that it cuts both ways.** A well-mixed or per-species model cannot represent protection that
depends on local cell density and a shared extracellular medium. A lattice model with local
neighbourhoods can. So the mechanism that breaks the simple test is the one that would justify
Phase 2 — which is a reason to state it plainly rather than a reason to lead with it.

### 0-pre, pre-registration — the inputs are fixed before the curve is computed

**A power analysis whose inputs the analyst chooses is a device for producing whatever answer is
wanted:** pick an effect size and get a required *n*; pick an *n* and get a detectable effect. 0-pre's
entire job is to be able to refuse, so the inputs and the refusal threshold are written down first,
the same discipline as pinning the observed acceptance band before running.

**The effect-size range comes from the anchors, not from what is achievable.** Chengdu is the low
anchor and the usable one — a threshold at 480 nGy/h inside a 4.7× total range, so composition
tracked a modest gradient. The pool's 0.416 Gy/h sits about 10⁶ above it, so there is a gradient
magnitude and no dose-response measured at that magnitude. The curve is therefore drawn across an
**assumed** range, and **the range and the end that triggers refusal are both written before
computing.**

### D2 — the horizontal contrast, and the refusal threshold written before anything is looked up

**The threshold and its basis are fixed first, because a look-up whose closing line is chosen
afterwards is the same instrument as a power analysis with chosen inputs, only faster.** The earlier
"about 2×" was gestured at by comparison to Chengdu and never derived, so it is replaced by a form
that can be audited:

> **Refusal threshold.** The design closes if the achievable horizontal dose contrast is below
> **one half** of the total dose range across which a composition threshold has actually been
> observed — Chengdu's 4.7× — giving a closing contrast of **2.35×**. The fraction is one half, and
> it is a declared choice rather than a derived one; what is derived is that the comparator is the
> only observed gradient, not an achievable one.

Written before the look-up so that the number cannot be moved by what the look-up returns.

**The comparator is sufficient-not-necessary, and the sharper one is unavailable — checked before
D2b, deliberately.** 4.7× is the *total span* across four groups with the threshold sitting inside
it, which establishes the span was **enough**, not that it was **needed**. The tighter comparator is
the adjacent-group bracket around 480 nGy/h, and that bracket is **not published**: the paper gives
Blank as 192.906 ± 5.05 and High as 910.964 ± 41.09 nGy/h, and reports Low and Medium only as
sampling distances (2 ± 0.5 cm and 7 ± 0.5 cm from source). The threshold sits between Low and
Medium, so the contrast that actually resolved it is **≤ 4.7× and not otherwise determined**.

**Direction of conservatism, which is now settled even though the magnitude is not:** since the
resolving bracket is adjacent-group and therefore no wider than the total span, 4.7× is an *upper
bound* on the necessary contrast, so **2.35× is likely conservative** — it may demand more contrast
than the effect required, erring toward refusing a workable design rather than accepting a broken
one. That is the direction to err in, and it is recorded rather than assumed. The threshold stands
at 2.35× until the two group means become available, and **it was looked up before D2b so that a
revision could not be selection dressed as sharpening.**

**The conservatism is CONDITIONAL ON A TRANSFER, and the transfer is declared rather than
established.** The first version of this sentence said the threshold is conservative *by
construction*, which overstates it. 4.7× is an upper bound on necessary contrast **at Chengdu's
absolute dose scale**, and the pool sits about 10⁶ above that scale. **Necessary contrast is not
scale-invariant if the response saturates**, and a saturating response requires *more* contrast at
high absolute dose, not less — which would flip the direction outright: 2.35× conservative against
Chengdu and **permissive against the pool.**

That does not resolve without a dose-response at pool magnitude, which is exactly what this
document already records as unavailable. So the honest statement is: **the threshold is conservative
if the Chengdu response transfers, and whether it transfers is assumed here rather than shown.** The
mechanism declaration in 0-pre — which mechanism the program predicts — is the same assumption
looked at from the other side, and if that declaration lands on saturation the threshold has to be
re-derived upward before D2b, not after.

Recorded rather than left implicit, because a later reader who finds the adjacent-group bracket
published elsewhere would otherwise see a declared fraction with no reason attached and revise it in
whichever direction the day required — and because a reason that is *stronger than the evidence*
invites exactly that revision when someone notices the gap.

**D2a — field steepness. Arithmetic, and it appears to be satisfied.** Water's linear attenuation
coefficient at Cs-137 and ~1 MeV is of order 0.07–0.09 cm⁻¹, so the bare half-value layer is roughly
8–10 cm; buildup from scattered photons stretches the effective falloff, plausibly to 15–30 cm per
halving at several mean free paths. Either way a 2.35× contrast needs **tens of centimetres** of
lateral separation, not metres.

*Provenance, in the register's own idiom: this is a **derivation** from tabulated attenuation
coefficients with named inputs. It is not a facility measurement and any real dose map supersedes it.
Filed as `coverage = derivation` so it cannot later be mistaken for something measured.*

**Scope: D2a establishes a NEGATIVE, not a dose.** What it shows is that *steepness is not the
binding constraint*, across a range of roughly 15–30 cm per halving with buildup and depending on
geometry. **It is not a dose estimate and must not be cited as one** for any other purpose — a
range that answers "is the gradient steep enough" is not a value that answers "what dose is at this
position."

**D2b — samplable positions. Facility access, and unresolved.** If the field is steep enough, the
binding constraint is not steepness but whether **two or more samplable surfaces exist at adequate
lateral separation at the same depth**. That is a question about a pool, not about physics.

**Depth tolerance and lateral separation are ONE parameter, and it is a ratio.** Dose varies
vertically at a comparable rate, so positions differing by ±10 cm in depth while separated by 20 cm
laterally would carry depth-induced variation comparable to the signal — the confound returning as
noise rather than as collinearity. **Depth tolerance ≤ 20% of lateral separation**, stated as a
ratio so that changing one constrains the other automatically.

**The look-up is asymmetric: it can close the gate, and it cannot open it.** Conditions (1) and (3) —
no detectable residual across plausible effect sizes, and buffering with no monotonic form to test —
are untouched by any contrast measurement. So the pass branch is **"proceed to D2b"**, never
"feasible", and a favourable steepness result advances nothing on its own. The gate's failure branch
is a disjunction of three conditions, so its pass condition is their conjunction.

**If no published horizontal dose map exists for an accessible pool, that is UNRESOLVED, not a
close.** An operational dose map is not the kind of thing that gets published — it is facility data,
not a finding — so absence of one is a **search-scope claim and not a physical one**, and reading it
as "the gradient must be inadequate" would close a gate on the wrong grounds. That branch is written
now precisely because it is the likely outcome and would otherwise be decided on the day. It
converts into a question for Hoffman under the rule that already governs the memo: state what the
design needs, ask whether the sampling exists, assert nothing about his facility from secondhand
reading.

**"Detectable residual" has no definition without a stated metric and test.** Bray–Curtis
dissimilarity with PERMANOVA and dose as a continuous covariate is the conventional shape, and its
power depends on within-group dispersion. That number is available: the amplicon survey verified on
2026-08-30 sampled 13 locations across a real basin. It carries no dose axis, but **dispersion across
a real basin's sites is the one number in it that transfers**, and the document says it was taken
from there.

**And that study is a near-miss rather than a bare absence, which is the better framing and the same
fact.** It is not evidence that nobody can sample a basin. It is evidence that **someone already
has** — 13 locations and depths, no dose measured, no field computed. The gap is that the group who
sampled a basin did not compute a field, and the missing half is exactly the half this repository
holds. That is a collaboration framing rather than an absence.

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

**The sharper risk is that dose is not separable from DEPTH, and the sampling geometry decides it.**
Water is the dominant attenuator — of order half a metre per tenth-value thickness for Co-60 — and
assemblies sit tens of feet below the surface, so the field spans many orders of magnitude
vertically. Depth also carries light, settling organic carbon, oxygen from the surface, and operator
proximity. **A regression of composition on computed dose over a vertical sample is a regression on
depth wearing a dose label**, and the residual-after-covariates framing in 0c is the right response
but may be unachievable in that geometry. Not hypothetical for the design: the nearest study found
in the same dated search sampled "13 different locations and depths."

**So the design that breaks the confound is horizontal.** Sample at **fixed depth** and varying
lateral distance from the racks, where dose falls off and every depth-linked covariate holds
approximately constant. **If horizontal sampling at fixed depth is not physically available in the
pool on offer, 0-pre closes on collinearity even though it does not close on effect size** — and that
is a sufficient reason to send nothing.

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

1. **0-pre**, whose failure branch is a **disjunction** of three conditions — so its pass condition
   is their **conjunction**, and satisfying one advances nothing:
   (a) no range of plausible effect sizes leaves a detectable residual;
   (b) the horizontal contrast is below 2.35× (D2a), **or** no two samplable surfaces exist at
       adequate lateral separation at the same depth (D2b);
   (c) cross-protection buffering leaves composition with no monotonic form to test.
   **A fourth outcome is not a close:** no published dose map for an accessible pool is UNRESOLVED,
   a search-scope claim rather than a physical one, and becomes a question for Hoffman.
2. **Phase 0**: the residual is below the floor with the control dissociating. Nothing is sent.
3. **Phase 1**: composition does not track computed dose. The coupling earns nothing.

Each has its disposition written before the run. That is the program's most defensible feature.
