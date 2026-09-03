# Working in this repository

This repository couples a cellular-Potts biofilm model to OpenMC transport in
order to ask whether a dose gradient measurably changes biofilm morphology. It
is a **scientific instrument**, and its central risk is not a crash — it is a
number that looks measured and is not.

Read this before reviewing or changing anything. It is shared by every agent
that works here: Codex reads `AGENTS.md`, Claude reads `CLAUDE.md`, and
`CLAUDE.md` imports this file so the two cannot drift apart.

## The six rules

Each of these is here because it was violated, in this repository, and the
violation shipped. They are ordered by how often that has happened.

### 1. A check that cannot fail is not a check

Every guard needs a **negative control**: a known-bad input it must reject. A
test that passes when it finds nothing looks identical to a test that can find
nothing.

- The claims-ledger guard compared prose to raw LaTeX. It reported "0 deleted
  claims survive" and would have reported that against a manuscript containing
  every one of them. Fixed by running it against the pre-revision manuscript
  recovered from git and requiring ≥18 detections.
- `authorization_criteria` mapped gate refusals to criteria by substring and
  treated an unmatched refusal as satisfied.
- `plot_layer` referenced two names it never imported.
- The figure phrase guard was "confirmed" by appending the two withdrawn
  sentences to a `.txt` sidecar as clean lines and watching the guard fail.
  `pdftotext -layout` reads a three-column figure **across**, so the real
  extraction interleaves the columns and the sentence is contiguous nowhere.
  `FIG-09` and `FIG-10` were enforced by nothing for a full commit while a
  control reported otherwise. Fixed by taking the known-bad from the artifact
  path — `c2219a2:preprint/figures/phase2_diffusion_cell.txt` — instead of
  writing one.

When you add a guard, add the input that proves it bites. When you review one,
ask what would have to be true for it to fail, and check that something makes
it so.

**Two failure modes are invisible from reading the assertion, and they are a pair.**
Everything above asks whether a check *can* fire. These ask something else.

*The assertion that was never load-bearing.* **Weaken it to something trivially true and
re-run the controls that are supposed to fail.** If those still fail, something else was
catching them and the assertion is not what it appears to be; if they now pass, the
assertion was load-bearing.

State it that way and not as "weaken it and see if the suite stays green" — **a green suite
stays green under any weakening**, because the assertion was passing anyway. The weakening
only discriminates against an input the assertion is meant to reject. That correction was
itself made here after writing the looser form into this file and running it: weakening the
bibliography set-assertion to `@test true` left 40 passing, which established nothing.

The real instance: a same-run identity assertion in `test_figure_staleness.py` that cannot
fail while the function it guards is correct, because that function locates the element by
content and replaces only within the located slice. It is kept, and relabelled as
documentation of an invariant rather than presented as a control.

*The input that never changed.* A mutation check can be correct in every line and still
prove nothing, because the mutation did not apply. Reproducing a finding here, a `replace`
against normalised text missed a source storing XML entities (`&#215;`), changed nothing,
and the guard's truthful "nothing missing" over an unmodified file read exactly like
confirmation. **Assert that the mutation applied, and assert WHERE it applied** — a bare
`mutated != source` passes on a replacement that landed in a different element entirely,
which was the second defect in the same three lines.

Neither is visible by reading the assertion. One is an assertion that cannot fail; the
other is an input that never changed.

**And "cannot fail" is not the same as "tautological", which is a distinction worth keeping
straight before deleting anything** — dimensional analysis and conservation arguments are
tautological and are among the most productive checks here. See *Necessary truths:
constraining versus restating* below.

**The assertion must consume the finest-grained thing the function returns.** If a
function computes categories and the test compares their union, the categories are
documentation. Three guards here have failed this way and each fix was one line:
a count compared where the function returned a set; a subset counted where full
coverage was asserted; and category membership printed on mismatch while only the
union was checked, so a real citation gap could be relabelled "deliberate
context" with the suite green. **In all three the discriminating data was already
computed and discarded at the assertion line.**

**The victim of a mutation is SELECTED FROM THE ARTIFACT, never recalled.** Extract
the runs, the rows, the citations from the file as it is now; pick one from that
extraction; mutate it. Do not reach for a string you remember being there.

This has failed twice with different causes and the same outcome — an attack that
reports ESCAPES having tested nothing, which reads as confirmation of the very
finding being chased. Once the mutation could not apply because the source stores
XML entities and the replacement was written against normalised text. Once the
victim phrase had been withdrawn from that figure months earlier and was no longer
in the file at all. "Assert the mutation applied" catches both, and did not become
reflexive; selecting from the artifact removes the opportunity.

**And it is the same root as reading a paper's title as its finding.** Both are
reaching from memory where the artifact was available: a withdrawn phrase recalled
into an attack, an implication recalled out of a title. The repair is one habit in
two places — open the thing before writing the sentence about it.

**Restore a mutation from a scratchpad copy, never with `git checkout`.** A
mutation check edits a real file; undoing it with `git checkout -- <file>` reverts
*everything* uncommitted in that file, not the mutation. That destroyed four
applied reference fixes here in one command. `checkout` reads as an undo and is
not one when the file carries unstaged work, and the loss is silent — the suite
goes green because the mutation is gone, and so is the work. Copy to the
scratchpad before mutating, restore from the copy after.

**A paper's title is not its finding, and the audit already knows the difference.**
The sharpest version of the summary defect has no external symptom at all: the
citation is correct, the DOI resolves, nothing drifts, every guard passes — and
the sentence still misreports the source, because the title was read as the
result. §2.6 cited Khajo 2011 as reporting protection from lethal-dose gamma; its
endpoints are EPR, absorbance, TBARS and Bi binding, and
`docs/research/radiotrophic_compatibility_audit.md:230` says in as many words that
the survival comparison is not performed there. That row was in this repository,
unread, when the sentence was written.

**No guard reaches this, and proposing one is the wrong lesson.** The obligation is
procedural: a paper entering a section on evidence has its ENDPOINT read from the
audit before the sentence is written. The audit exists because titles and endpoints
diverge, and it already holds both.

**And draw the known-bad from the artifact path, never from your own hand.** A
control you construct tests your idea of the failure; one recovered from where
the real input comes from tests the pipeline. Where a guard reads a *derived*
artifact — an extraction, a render, a serialisation — **the derivation's shape
is part of the contract**, and a hand-written input silently opts out of it. The
sidecar case is the clean example and it is retroactive: `.txt` sidecars exist
because the FIG-01–07 lesson was that images hide claims, so the extraction is
what gets read — and nobody had written down that the *shape* of that extraction
was therefore load-bearing too. This failure mode is not domain-versus-range,
not a coverage gap, and not a vacuous assertion. It is a control that never met
the pipeline.

**A control enters by the production door.** If production discovers files by
walking a root, normalises them in a scanner, and classifies the scanner's
result, the control must place its known-bad under a root and call that scanner.
Calling the normaliser directly, or comparing the vocabulary to a string in the
test, proves a helper and leaves the production path free to bypass it. This is
not satisfied by sharing an assertion: the control must traverse the same input
boundary and return the same finest-grained result production consumes.

**An exemption is keyed on the property that earns it, at the location where it
is true.** A directory name, suffix, token value or other surface marker can
select a candidate; it cannot justify skipping one. Verify the relation itself:
a cache file is current bytecode derived from the sibling source, a clipped token
is the one measured run at the one figure location, and a recording document is
the exact declared path. Broad markers turn a local exception into an input
channel production never reads, which is the exclusion outliving its reason.

### 2. A skipped test is uncovered surface, not a neutral fact

`plot_layer` shipped a `NameError` on its first rendering line behind eight
skips, because pyvista is installed nowhere that runs the suite. The report said
`215 passed, 8 skipped` and was green over a dead function.

Where a missing dependency gates a test, add a bare-tier check that exercises
whatever does not need it — name resolution, schema shape, argument contracts.
`test_every_global_the_render_path_uses_actually_resolves` is the pattern.

Read a skip count as a question, never as a pass.

### 3. Never default to pass

A lookup, a substring map, or a dispatch table must **refuse** the case it does
not recognise. An unmapped refusal is not an absence of refusal.

This one is dangerous specifically where the expected state is "not ready":
an unreachable milestone and a correctly-withheld milestone print the same
thing, so nothing ever looks wrong. `authorization_criteria` was unreachable
for one release and nobody could have noticed by reading its output.

### 4. The producer declares semantics; the consumer must not assume them

`cell_id == 0` is empty space. `generation == 0` is a founder. The renderer
cannot tell these apart, so it must not try: `Layer.background` is declared by
the code that built the field.

**And the default must not be a statement.** `background=None` says "every cell
carries information, hide nothing" — true for `generation`, wrong for `omega_b`.
While that was also the *default*, a producer who simply forgot made the claim
silently, and two layers did: `omega_b` drew its false cells as part of the
metric's region, and an upsampled overlay grew a shell of empty space. The
default is `UNDECLARED` and is refused at write time. Saying "every cell is
informative" is still available; it has to be said.

**A consumer that reads the declaration must read all of it.** `occupied_mask`
tests `> 0`, which is one encoding — 0 empty, negatives out-of-domain — so a
source declaring `background = -1` is refused rather than silently mis-masked.
Validate every half of a contract, and give the control data that can fail each.

Applies broadly — units, sentinels, fill values, ordering, extensive versus
intensive. If a consumer needs to know, the producer states it in the data.

### 5. A ladder holds the object fixed

A convergence study varies **sampling** and nothing else. A parameter that
changes the system under study is a different experiment.

A pitch of 3.2 µm does not tile the spheres' 24 µm axis, and rounding the voxel
count silently resized the box to 25.6 µm while the analytic truth kept dividing
by 24 — a 6.7% denominator shift reported as rasterisation error. Refuse it, and
**name what you refused**: a row that is merely absent is indistinguishable from
a case nobody tried. (The first fix for this dropped the pitch from the default
list, which reintroduced exactly that ambiguity.)

### 6. Some claims are not ours to make

- **`D-APPROVAL` is institutional biosafety evidence.** No agent, and no human
  working through an agent, can declare it. It was once asserted as the words
  "d approved"; the refusal that stopped it was human and the code was a
  presence check. Biosafety follows **strains, not species**.
- **δ stays unset**, and must never be chosen from any effect measured in this
  repository. It is a prospective decision policy; selecting it from observed
  data is how a threshold comes to certify the noise that produced it.
- **Reference D stays `NOT_EVALUATED`** and `CAMPAIGN_READY` stays no until an
  institution says otherwise, in a document, with an identifier.

If a change would move any of these, stop and say so. Do not implement it.

## Necessary truths: constraining versus restating

Rule 1 can be misread as "a tautology is a defect." It is not, and flattening three
different things under that word has already caused trouble here. **The working distinction
is not tautological versus not. It is whether a necessary truth CONSTRAINS — forbids
something — or merely RESTATES.**

**1. The vacuous assertion. A defect, full stop.** An assertion that stays green when
weakened to something trivially true forbids nothing. It consumed a slot in the suite and
reported confidence it had not earned. Delete it, or — if it documents an invariant the
surrounding construction actually provides — say so in the file and stop calling it a
control. The mechanical test is in rule 1.

**2. Contingent circularity. Not a tautology, and not empty.** A model that recovers its
own assumption is not true by logical form; it is a correct implementation of a premise.
The accurate statement is that **the model is not independent evidence for its own
assumption** — non-informative about the thing it appears to be about, which is different
from vacuous. This repository has several: `README.md:413` records that the melanin
ordering among the three producers "is set by the input `α_M`", and
`docs/preprint_revision_plan.md:107` records a §6.3 scatter that "is a scatter of Table 2's
own entries." The `β_ion` column has the same shape, spending survival evidence on a
tropism coefficient. None of these is a bug in the code. Each is a claim that must not be
reported as a result.

**3. The load-bearing necessary truth. Tautological and indispensable.** Dimensional
analysis tells you nothing you did not put in — Buckingham π adds no physics — and it is
how the `D_eff` collision was caught and how the units error behind
`RADIODIALYSIS: BLOCKED` was found. Conservation arguments are tautological in the same
sense. So is the shell derivation in `docs/guides/calculus_in_this_code.md`: two face areas
and a volume, no new physics, and the most useful page in the document. **What makes these
productive is that the consequences are not obvious from the premises even though they
follow necessarily.** They forbid things — a film-scale `D_eff` above `D₀`, an occupancy
fraction in a `g cm⁻³` slot, a slab operator with an `r` in it.

**A check that forbids nothing and a model that predicts only its inputs fail the same
way**, which is why they belong in one place rather than two.

And the same test applies to the apparatus as a whole. A repository of guards that
collectively forbid nothing is the first failure at a larger scale. **Naming which claims
are pinned and which are not is what stops that** — the coverage output that lists
undetectable ledger rows by name, the guide's paragraph saying which of its prose claims no
table pins, the SCOPE line on every absence claim. Those exist so the apparatus cannot
quietly become a restatement of its own inputs.

## Correcting a published number

Any number this repository has published and later found wrong gets **all** of:

1. the code fix, with a negative control,
2. a row in `data/claims_ledger.csv` with verdict `delete` / `restate` /
   `requalify` / `keep`, stating what was wrong and how it was caught,
3. a marked correction in the document that carried it — not a silent edit,
   **which names the ledger row and never repeats the withdrawn content.** The
   row holds the wording; the artifact must not. Repeating a withdrawn DOI, range
   or formula inside a correction comment leaves it in the file someone copies
   from, and every guard written afterwards will fire on the correction itself.
   This happened three times in one session — a withdrawn DOI in a LaTeX comment,
   a second in the same file, a withdrawn range in an R comment — each time in a
   correction whose whole purpose was removing that string.

   **The reason it recurs is structural, not careless: writing a correction and
   writing a guard against the thing corrected are the same act, and the
   correction is authored in the window where the guard does not exist yet.** So
   the convention has to be a habit rather than something the suite catches,
   because at the moment of writing there is nothing to catch it.

4. **and, when what you closed was an ABSENCE, an audit of everything that
   asserted that absence.** Step 5 audits what was computed *from* a corrected
   number. This is its sibling pointed backwards: **a document whose job is
   recording what is unenforced becomes wrong the moment the enforcement ships.**
   Every "remains unenforced pending X", "not built", "named, not fixed" entry is
   a claim with an expiry date nobody sets.

   It recurs by construction, not by carelessness — the entry is written when the
   gap is real, and closing the gap is a separate act in a separate commit that
   has no reason to touch the document. On 2026-08-31 a red-team line read
   *"RM-G04-01 remains unenforced pending RETRACTED_IN_SOURCES"* for several hours
   after that guard shipped with RM-G04-01's string in its vocabulary.

   **Attribution-with-correction is permitted only while the string is in no
   retraction vocabulary.** Naming what an earlier report claimed and correcting
   it in the same breath is ordinary practice, and
   `coupling/biofilm_openmc/mesh.py:47` does it well. But once a string enters
   `RETRACTED_IN_SOURCES` or `RETRACTED_CITATIONS`, the artifact carries it in no
   form and the row holds it alone. Vocabulary membership is the trigger because
   it tracks the two harms separately: the copy risk is low while the string is
   visibly marked withdrawn in place, and the guard-collision risk is *created* by
   adding it. **A comment that was compliant becoming a failure the day its string
   is guarded is the guard working, not a regression.**

   **A superseded record is exempt, and it declares that with a fixed token so
   the scope of this audit is a grep rather than a memory.** The first line of a
   superseded document is:

   ```
   SUPERSEDED: <date> — <what it is the record of>
   ```

   Step 4's scope is then `grep -L '^SUPERSEDED:' docs/**/*.md`, which is the live
   set. A live record may not be stale; a superseded one may.

   **The token is deliberately not the English word, because that word is already
   taken here and means the opposite.** `docs/calibration/reference_d_measurement
   _protocol.md` opens `**Frozen:** 2026-08-15` — frozen as in *pinned and
   binding*, a live authoritative document. Matching on "frozen" or "superseded"
   as prose would exempt it, and would also match
   `docs/correspondence/wan_v11_note.md`, where "superseded" appears in the body
   describing the note's own content. Surveyed: 3 of 38 documents carry anything
   banner-like, and all three mean different things. A natural-language marker
   would have been actively wrong rather than merely unreliable.

5. **and, if the row prescribes a NUMBER, it states the derivation rather than
   the result.** "The argument must equal the bibitem count" survives drift;
   "correct the argument to 44" does not. `PP-REF-07` said 44, was right when
   written, and was wrong when applied months later against 61 bibitems —
   applying the row's own number would have introduced the defect it exists to
   remove.

   **This is the third variant of one family and it is not the same as rule 4.**
   Rule 4 is *closing a gap makes the records of that gap stale*. This is *a
   prescribed fix carries a derived value, and the derivation's inputs move while
   the prescription sits unapplied*. The fix was never wrong; its inputs were.

   **Nothing mechanical caught it, and nothing could have.** No guard fired — the
   number was noticed while applying the row by hand. A prescribed value with a
   stale derivation is invisible to every tier here, because **the row is
   internally consistent, the artifact is internally consistent, and only their
   relation is wrong.** There is no artifact to compare against: the prescription
   is not yet applied, so nothing has both halves. Hence a convention rather than
   a check. Any row prescribing a number carries the rule that produced it.

6. **and an audit of everything computed FROM it.** A withdrawn measurement
   takes its dependents with it. The 3.2 µm ladder row was withdrawn and its
   own headline corrected twice, while the interface-area range published a few
   lines below — `0.41–0.51 across 16×` — kept both of its halves from that row:
   0.41 was its error, and 16× was 3.2 → 0.2. The correction stood next to a
   number that contradicted it. Grep the value, and grep what it was divided by.
7. If an earlier correction claimed "every other conclusion is unchanged" and
   this one disproves it, **delete that sentence — do not reword it.** A
   paraphrase is the same claim: "the conclusions themselves are unchanged"
   replaced the original here and sat one line above the paragraph documenting
   the second number that moved. Say which numbers moved and stop. A correction
   is not the place for a reassurance whose scope nobody has checked.

A correction that only changes the number teaches nobody why it was wrong, and
this repository has now corrected its own corrections twice.

**A write to an external system is read back before it is reported as done.**
This is the identifier rule pointed at a different verb: there, an identifier is
read back against the registry that issues it; here, a write is read back from
the system that accepted it. Both exist because the failure is silent.

**The general form, because there are now three instances and one shape: an
operation's result is confirmed before anything is concluded from it. Writes get
read back; reads get checked non-empty.**

The read direction's instance: `gh run view --log --job=<id>` returned **0
lines**, and a grep over it for `undefined` found nothing. That silence was one
step from being reported as a clean LaTeX build. **An empty file and a clean file
grep identically**, and nothing in the pipeline flags that nothing was produced.
The same job's log fetched through `gh api
repos/<owner>/<repo>/actions/jobs/<id>/logs` returned **121,872 bytes**, and only
then did the greps mean anything -- the final `latexmk` pass had zero undefined
references, which the empty log could not have told anyone.

The rule is not "prefer `gh api`"; that is a workaround for one tool on one day.
It is: **assert the artifact is non-empty before searching it**, whatever
produced it. A search over nothing returns nothing, and nothing is
indistinguishable from a clean result at the point where the conclusion gets
drawn. This is the no-op mutation in a different costume -- the command
succeeded, produced nothing, and the emptiness was invisible exactly where it
mattered.

The instance: `gh pr edit` failed on a deprecated Projects-classic GraphQL call,
printed that deprecation as a *notice*, **exited zero, and left the body
unchanged**. Nothing in the output distinguished it from success. The correction
only landed because the live body was grepped afterwards, and it took a REST
`PATCH` to actually write. **So the class is not `gh pr edit`** — it is any
subcommand that can warn, exit zero, and not do the thing, which is why the
read-back attaches to the write rather than to the command.

A repo-wide grep found **no scripted external writes at all**: every `gh` call
in `scripts/` and `CLAUDE.md` is a read (`pr view`, `repo view`, eight
`query(` and zero `mutation`), and CI's one `curl` is already followed by
`sha256sum -c -`, so a truncated fetch fails loudly. **The exposure is entirely
in ad-hoc agent practice**, which is precisely why it belongs in this file
rather than in a lint.

**The audit-dependents step covers artifacts outside the repository, and it
needs an enumeration rather than a grep.** A withdrawn claim has now survived in
two dependents that a repository-wide search could not reach: the Wan
correspondence draft, and the PR description itself — which carried the retracted
zero, *and* the bound argument that was wrong by three orders of magnitude, for
three days after the manuscript was corrected. `PP-62-04` records both.

The step had been run as a grep over tracked files, so its scope was silently
*what I can grep* rather than *what asserts this claim* — the same substitution
that makes an absence claim scoped to one file look like a property of the
repository. **The PR description is the first thing a reader of this work sees
and it is not a file.** So the enumeration is written down when a claim is
published, not reconstructed when it is withdrawn: PR and issue bodies, review
comments, correspondence already sent, slide decks, anything quoting a number
that a later verdict can move.

**Record the sink type, not only the path, because sinks fail in different ways
— and not the ways you would guess.** The withdrawn `PP-62-04` zero was traced
into every sink it reached:

| sink | editable | reachable by grep from a work tree | corrected after |
|---|---|---|---|
| `docs/correspondence/wan_v11_note.md` | yes | yes | days, on a Codex P1 |
| PR #23 description | yes | **no** | **three days**, when a human read it |
| commit `d404438`'s message | **no** | no — needs `--all` and `%b` | **~10 hours**, by `5788db5` |
| bioRxiv posting | — | — | never posted; `README.md` declares it, dated, as a search result rather than proof |

**The append-only sink was corrected fastest and the editable one slowest**,
which is the opposite of the intuition that says an editable artifact is easier
to fix. Editability is not what decides whether a sink gets corrected; **whether
correcting it is already part of a workflow** is. Writing another commit is the
normal act, so the commit-message error was answered the same day by the ordinary
motion of the work. Editing a PR description is out-of-band, triggered by nothing,
so it sat until someone happened to read it.

Two consequences. **Prefer sinks whose correction rides an existing motion**, and
where a claim must go into an out-of-band sink, the enumeration is the only thing
that will bring anyone back to it. And **do not try to rewrite an append-only
sink**: `d404438` stands, `5788db5` corrects it, and the pair is the record. A
history rewrite would destroy the evidence that the correction happened.

**An identifier that routes a reader to an external retrieval is read back
against its source before the paragraph ships.** Accession numbers, DOIs, PMIDs,
PDB and NCBI ids — anything whose only job is to send someone somewhere. Read it
back from the registry that issues it, not from the document that quoted it to
you, and not from a search engine's gloss.

The failure mode is specific and it is worse than an ordinary typo: **a wrong
identifier produces a false negative that is attributed to you.** The reader
looks, finds nothing, and concludes the document does not exist — the opposite
of what the paragraph was arguing. This repository shipped `ADA303995` inside
the paragraph arguing that `ADA307995` should be obtained first, so the
paragraph's own recommendation, followed literally, would have failed and looked
like confirmation.

**Both directions of the check have now paid.** A review reported Hannan et al.
1986 as a microbial neutron RBE source. Reading `PMID 3533817` back through the
NCBI eutils API returned title, authors, journal, volume, issue, pages and the
abstract's RBE values — all matching, so the claim was adopted on the index's
authority rather than the review's. The *same* review reported a DTIC report as
freely readable; that claim had no identifier-backed check available, did not
reproduce on retrieval, and was not adopted. **The rule is cheap enough to run
on every identifier and it discriminates in both directions**, which is why it
is worth more than the ordinary caution it looks like.

### Address an entry by its key, never by its printed number

A `\bibitem` key, a ledger `claim_id`, a `\label` — these are stable. The number a renderer
prints for them is a function of everything above them, so one insertion silently re-points
every reference below it, and **nothing errors**: the row still reads plausibly and now names a
different paper.

Measured, not asserted. On 2026-09-01 all six `PP-REF-*` rows located themselves as
"References [11]", "[12]", "[15]", "[26]", "[4]", "[32]". Resolved against the `.tex` those
ordinals are `slade2011`, `daly2009`, `xavier2005`, `karley2018`, `meskauskas2004` and
`kauffman1989` — none of which is the work the row describes. Six for six, and already wrong at
`HEAD` before any entry was deleted, so the cause was not an edit; it was addressing a render
artifact. `docs/preprint_revision_plan.md` had the same defect, naming `audit2026` as "[41]"
when it is entry #61.

This is the identifier read-back rule pointed at an internal target: the `.tex` is the issuing
registry, the key is the identifier, and the printed number is a rendering of it. Enforced by
`test_no_ledger_row_locates_itself_by_printed_reference_number` in
`calibration/tests/test_claims_ledger.py`, with its own control.

**And when a number cannot be recovered, say so rather than guessing.** `PP-REF-06` also cited
"[20]" and "[35]"; those two keys were not recoverable, and the row now records them as
unrecovered. Guessing which entries they were would be the failure this rule exists to stop.

## The suites

Four, and they are separate because they answer different questions.

```sh
cd coupling     && .venv/bin/pytest -q tests        # transport coupling, viewer
cd calibration  && ../coupling/.venv/bin/pytest -q tests   # biology, ladders, gates
cd contract     && ../coupling/.venv/bin/pytest -q tests   # the shared vocabulary
julia --project=. tests/runtests.jl                 # the CPM kernels
```

`openmc-integration` is `workflow_dispatch`-gated: the hosted runner has no
photon library, and tests needing nuclear data must report SKIPPED rather than
be represented as passed.

## Before merging

```sh
scripts/preflight_merge.sh <pr>
```

It refuses while any review thread is unresolved and while the newest review
covers a commit that is no longer the head. Green CI is not sufficient — CI
cannot tell an addressed finding from an ignored one.

Resolve a thread when it is genuinely addressed. Resolving is the record that
someone read it.

**So commit freely while a review request is outstanding, and repost the request
when the committing stops.** The alternative — freezing the branch once a
request is posted — was considered and rejected, because the freshness check it
would enforce by hand is already enforced above: `preflight_merge.sh` refuses
while the newest review covers a commit that is no longer the head. A hand
convention duplicating a tool that already refuses is a second rule to drift
from, and drifting between two conventions is worse than either.

**What the request must carry is the head SHA it was posted against.** A review
scoped to a named commit stays interpretable after the branch moves — it is a
review of that commit, and the reader can diff forward. An unscoped review of an
unnamed head is the thing that silently misleads.
