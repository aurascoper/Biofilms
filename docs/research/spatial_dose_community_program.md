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

> **Refusal threshold — NO AVAILABLE ANCHOR BOUNDS THIS IN THE REQUIRED DIRECTION.** This is a
> finding with a decision attached, not a to-do.
>
> **The arithmetic is not repairable in the direction the original derivation assumed.** A
> conservative gate needs a **lower** bound on the necessary contrast, and Chengdu supplies only an
> upper one: a threshold at 480 nGy/h inside a 4.7× span, with Low and Medium given as sampling
> distances rather than dose rates, so the smallest contrast across which the effect was actually
> resolved is unpublished. No fraction of an upper bound is a floor — if the necessary contrast is
> 3×, a 2.35× gate admits an achievable 2.5×, a design that looks like it cleared a pre-registered
> bar.
>
> **Why unset is an improvement and not a regression:** an unset gate refuses everything, while a
> wrong gate admits exactly the cases it should not. Blocking D2b outright is the safe state.
>
> **Three options, and one has to be chosen before D2b runs:**
> 1. **Obtain the bracket** — author contact for the Low and Medium group means, or another study
>    reporting a resolved adjacent pair.
> 2. **Derive the threshold from dispersion rather than from an anchor's span** — the within-group
>    variance route, which needs no anchor at all and connects to the same 13-site dispersion 0-pre
>    already uses for the noise floor.
> 3. **Declare that 0-pre cannot set this threshold, and close on that.** A real outcome, not a
>    failure to decide.

Written before the look-up so that the number cannot be moved by what the look-up returns.

**The comparator is sufficient-not-necessary, and the sharper one is unavailable — checked before
D2b, deliberately.** 4.7× is the *total span* across four groups with the threshold sitting inside
it, which establishes the span was **enough**, not that it was **needed**. The tighter comparator is
the adjacent-group bracket around 480 nGy/h, and that bracket is **not published**: the paper gives
Blank as 192.906 ± 5.05 and High as 910.964 ± 41.09 nGy/h, and reports Low and Medium only as
sampling distances (2 ± 0.5 cm and 7 ± 0.5 cm from source). The threshold sits between Low and
Medium, so the contrast that actually resolved it is **≤ 4.7× and not otherwise determined**.

**What the bracket does and does not give:** it is adjacent-group and therefore no wider than the
total span, so 4.7× bounds the necessary contrast from **above**. An upper bound alone cannot set a
floor — see below — so no direction of conservatism follows from it. The comparator was looked up
before D2b so that a revision could not be selection dressed as sharpening; what that look-up
established is that the interval is (0, 4.7×] and nothing narrower.

**Two independent defects sat in one sentence, and fixing the second first left the first
untouched.** The withdrawal of "conservative by construction" weakened the claim along the
*transfer* axis — 4.7× bounds necessary contrast at Chengdu's absolute scale and the pool sits ~10⁶
above it, so a saturating response would need *more* contrast, not less. That argument stands. But
the sentence was **also** wrong arithmetically, for a reason independent of scale: **halving an
upper bound does not yield a lower bound.** Knowing only that the necessary contrast is at most 4.7×
tells you nothing about whether 2.35× clears it. Raised by Codex on pull request #23, and the
threshold is now unset rather than re-justified.

**What would pin it:** the adjacent-group bracket around 480 nGy/h. The paper gives Blank as
192.906 ± 5.05 and High as 910.964 ± 41.09 nGy/h and reports Low and Medium only as sampling
distances (2 ± 0.5 cm and 7 ± 0.5 cm from source), so the two means needed are unpublished. Until
they are obtained the necessary contrast is known only to lie in (0, 4.7×], and no gate can be set
on that interval.

**D2a — field steepness. Arithmetic, and it appears to be satisfied.** Water's linear attenuation
coefficient at Cs-137 and ~1 MeV is of order 0.07–0.09 cm⁻¹, so the bare half-value layer is roughly
8–10 cm; buildup from scattered photons stretches the effective falloff, plausibly to 15–30 cm per
halving at several mean free paths. Either way a contrast of a few-fold needs **tens of centimetres**
of lateral separation, not metres — which is the shape of the answer regardless of where the
threshold eventually lands.

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

### D3 — the neutron weight on a microbial target, and a gate that refuses without closing

**The four dose figures this program reasons from are all gamma, and a research reactor is a mixed
field.** They are: Chengdu's ambient rates at the Th-232 site, Shuryak's Cs-137 exposure rates, the
0.416 Gy/h reported for areas holding irradiated fuel elements, and the 248 Gy–2 kGy D₁₀ range for
the six spent-fuel-pool isolates cited in the Hoffman memo. That is the whole set — an enumeration
of what this document cites, not a claim about dosimetry at large. A neutron component does not enter any of those through absorbed
dose alone; it needs a biological weight, and **no weight for a microbial target has been
verified here. So the gate refuses: no neutron-bearing exposure is designed, costed or asked for
until one is sourced.**

> **SUPERSEDED 2026-09-01 as to existence, not as to the refusal.** A microbial value was
> verified on that date — see the UPDATE below. The gate still refuses; the reason is now
> transfer rather than absence. The sentence stands as the record of what was true on
> 2026-08-31, which is the point of dating it.

**The refusal is typed OPERATIONAL-pending, not STRUCTURAL, and the distinction is the reason for
writing it down.** STRUCTURAL is the one label in this register that says stop looking, and it
would be the most consequential verdict in the document resting on the least scoped claim — in a
document going to a facility whose own people are the ones who would know. Nothing found so far
earns it. Operationally nothing changes today: the gate refuses either way. What changes is that
the refusal stays reopenable by someone searching the way the field is actually organised.

**Four web searches, 2026-08-31, term sets verbatim.**

1. `microbial neutron RBE relative biological effectiveness bacteria fungi`
2. `bacterial spore inactivation D10 fast neutron mixed field gamma comparison reactor sterilization dosimetry`
3. `neutron RBE bacterial spore inactivation D10 mixed neutron gamma field`
4. `fungi neutron irradiation relative biological effectiveness spacecraft shielding microbial inactivation`

**The term-set diagnosis is confirmed, and it is only half the barrier.** A microbial neutron RBE
is largely not indexed under "RBE": it sits in sterilisation dosimetry, BNCT bystander work,
spacecraft-shielding and planetary-protection validation, and reactor inactivation studies, and it
is reported there as a D₁₀ under a stated mixed field with a gamma comparison alongside — an RBE in
substance with no RBE in the title. Searches 1 and 3, which asked for the conjunction
*microbial + neutron + RBE*, returned nothing usable. Searches 2 and 4, which named those
literatures instead, returned three on-point documents. **So an empty result under the first term
set discriminates nothing** — it separates "the term is not how the field names it" from "no value
exists" not at all. **But the second term set did not clear the gate either, and the reason is now
retrieval rather than naming:** each of the three is a document whose full text did not open from
here.

### UPDATE 2026-09-01 — a microbial neutron RBE exists, and the refusal changes grounds

**Verified from a primary index, not from a review's summary.** NCBI eutils `esummary` and
`efetch` on PMID 3533817 return: Hannan MA, Paul M, Phillips RL, *"Fast neutron r.b.e. for
lethality and genotoxicity in a wild-type and a repair-deficient strain of yeast"*,
*Int. J. Radiat. Biol.* **50**(5):811–24, 1986 Nov. Abstract, verbatim: the r.b.e. *"varied
from 2.7 to 4.1 for lethality, 2.8 to 7.1 for reverse mutation and 3.5 to 7.8 for mitotic
gene conversion"*, for 11 MeV cyclotron neutrons against ⁶⁰Co gamma in *S. cerevisiae* D7.

**This meets the clearing condition this subsection wrote for itself:** a value with a
stated endpoint (lethality by colony formation), a stated neutron energy, and a stated
photon reference. **The existence premise stated above is therefore superseded**, and is
marked as such where it stands rather than restated here.

**The gate still refuses, on better-scoped grounds.** Yeast is a eukaryote and a fungus;
three of the seven modelled species are fungi, so this is nearer than a mammalian import,
and it is not a bacterial value. Nor is 11 MeV cyclotron the fission spectrum a pool or a
reactor presents. **The open question is now transfer — eukaryote to prokaryote, and one
neutron energy to another — rather than existence**, which is a materially different and
much smaller gap.

**And one thing the same review claimed does NOT reproduce.** It reported DTIC ADA307995 as
Distribution A, freely readable, the earlier 403 being transient bot-detection. **A second
fetch on 2026-09-01 returned 403 again**, and the AFRRI boilerplate surfaced by search reads
*"available to qualified users from DTIC, or others may contact the National Technical
Information Service"* — which is in tension with unlimited public release. The review is
also internally inconsistent: it states authors and distribution were *"confirmed verbatim
from the report's SF298 page"* while stating that the year on that same page and the results
tables were not extractable. **A claim to have read a page that also reports that page's
contents as unreadable is not adopted here.** The document's status is unchanged: unread,
with retrieval unresolved.

**Three leads, named so that reopening is a task and not a re-search.**

1. **DTIC ADA307995**, *Neutron and γ-Ray Radiation Killing of Bacillus Species Spores: Dosimetry,
   …* — on point by title, and the closest thing found to the exact quantity. **HTTP 403 from here,
   unread.** That is *source inaccessible, not value absent* — the same state `HOFFMAN-11` records
   for NACE-2019-12944, and it changes the day someone with DTIC access opens it.
2. **NASBEE**, the neutron exposure accelerator for biological-effects experiments at NIRS
   (*Radiat. Phys. Chem.*, ScienceDirect, paywalled). An **RBE of 3.54 at D₁₀ for 2 MeV average
   neutrons** appears in the abstract-level material returned. **The target is not established as
   microbial from what was read**, and NIRS's facility is built for mammalian radiobiology, so this
   is a lead and explicitly not a number to use. Recorded because it fixes the *form* the answer
   takes — an RBE quoted at D₁₀ against a stated neutron energy — which is what a search should
   look for.
3. **arXiv 2408.10929**, a comparison of spallation, reactor and compact neutron sources for
   genetic mutation, reporting fast neutrons as an effective mutagen in plants and microorganisms.
   The PDF returned image streams and no extractable text from here.

SCOPE: the four web searches quoted above, run on 2026-08-31, plus one fetch attempt each at the
DTIC and arXiv documents, both of which failed to return readable text. No controlled-vocabulary
database was searched — not PubMed under MeSH, not INIS, not OSTI — no sterilisation-dosimetry
handbook was opened, no BNCT dosimetry review was read, and no author was contacted. **This is a
claim about a search, not about a literature**, and the three leads above are evidence that the
literature is there.

**What would clear this, in the order it should be tried.** OPERATIONAL-pending is only meaningful
if the pending part names a task, so it does:

1. **DTIC ADA307995, by whatever access exists.** It is the on-point document and the cheapest
   thing on this list. Public access is 403 from here; a `.mil`/`.edu` affiliation, an
   interlibrary-loan request, or an NTIS order all reach it. **Clearing condition: the report
   states a D₁₀ under a characterised neutron field with a gamma comparison, for a named
   organism.**
2. **Ask, before buying.** The Hoffman question costs nothing and a radiation-materials person at a
   reactor either knows the quantity or knows who does. Note this is now a *retrieval* question and
   not an existence one — "can you get at this" is a much easier thing to answer than "does this
   exist," and it should be asked in that form.
3. **NASBEE, through an institutional subscription**, to settle whether the 3.54 figure is
   microbial or mammalian. **Clearing condition either way: it resolves whether that number may be
   cited at all**, which is worth the look even though a mammalian answer clears nothing else.
4. **A controlled-vocabulary search**, which this program has not run: PubMed under MeSH, INIS,
   OSTI. This
   sits below the three above because the leads already in hand are more specific than anything a
   fresh search returns, but it is what would convert the SCOPE line's disclaimer into a bounded
   enumeration.
5. **Write to the authors.** Last because it spends an introduction, and the four above are
   cheaper.

**None of these is scheduled and none is anyone's assignment yet.** Listing them makes the gate a
task rather than a state; it does not make the task started.

**Why NASBEE in particular is deferred, and what would change that.** The obvious argument is
ordering — D3 fires ahead of D2's threshold and of D2b's sampling look-up, both unresolved,
so paying for NASBEE clears the third gate while the first two are open. **That argument is
contingent and should not be the one on record:** if D2b came back favourable tomorrow it
evaporates, and the purchase would still be the same purchase.

The argument that holds under any D2b outcome is about the paper. **NASBEE's RBE of 3.54 at
D₁₀ has no established microbial target** — NIRS builds that facility for mammalian
radiobiology — so buying it resolves this gate *only if* the target turns out to be
microbial. SCOPE ON THAT: the figure and its context were read at abstract level only, from
the search results of 2026-08-31; ScienceDirect returns the full text behind a paywall from
here. **"Target not established" is a statement about what was readable, not a finding that
the target is mammalian** — which is precisely why the purchase would resolve something, and
equally why the cheaper routes should be exhausted before paying to resolve it. If it does not, the purchase returns the gate exactly where it stands, having
established that one lead does not apply. That is a real outcome and a poor thing to pay for
when the two routes ahead of it in this list would establish the same thing at no cost.

> **Trigger for revisiting:** DTIC ADA307995 and the Hoffman question both returning nothing,
> **and** NASBEE's target independently confirmed as microbial *before* purchase rather than
> after.

The second clause is the load-bearing one. It converts "buy it and see" into a precondition,
which is the shape the rest of this register uses: a cost paid against a stated outcome
rather than against a hope.

**What enforces this gate, precisely, because "an unenforced paragraph" understates it.**
Nothing in the suite reads D3 — it is prose with a declared type. But the *precondition* it
rests on is guarded: `coupling/biofilm_openmc/model.py:79` sets `particle="photon"` and
`coupling/tests/test_model_build.py:25` asserts it, so switching the transport to a neutron
or mixed field fails a test loudly. **The gate cannot be walked past silently, only
deliberately** — and the deliberate case is a person reading a failing assertion and deciding
what to do, which lands them here. For a judgement that cannot be mechanised, a guarded
precondition whose failure routes to the paragraph is the strongest available form.

**It converts into a question rather than sitting here**, under the rule that already governs the
memo: state what the design needs, ask whether the quantity exists, assert nothing about the
facility from secondhand reading. The question is now in the memo's list — whether a usable neutron
RBE for bacteria or fungi exists at all, or whether that is simply not how dosimetry is done for
microbial targets. A radiation-materials person at a reactor either knows or knows who does, and
that is a cheaper route to the answer than any search run from here.

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
   (b) the horizontal contrast is below the threshold — **currently unset, see D2** — **or** no two
       samplable surfaces exist at adequate lateral separation at the same depth (D2b);
   (c) cross-protection buffering leaves composition with no monotonic form to test.
   **A fourth outcome is not a close:** no published dose map for an accessible pool is UNRESOLVED,
   a search-scope claim rather than a physical one, and becomes a question for Hoffman.
2. **Phase 0**: the residual is below the floor with the control dissociating. Nothing is sent.
3. **Phase 1**: composition does not track computed dose. The coupling earns nothing.

Each has its disposition written before the run. That is the program's most defensible feature.
