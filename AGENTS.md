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

## Correcting a published number

Any number this repository has published and later found wrong gets **all** of:

1. the code fix, with a negative control,
2. a row in `data/claims_ledger.csv` with verdict `delete` / `restate` /
   `requalify` / `keep`, stating what was wrong and how it was caught,
3. a marked correction in the document that carried it — not a silent edit,
4. **and an audit of everything computed FROM it.** A withdrawn measurement
   takes its dependents with it. The 3.2 µm ladder row was withdrawn and its
   own headline corrected twice, while the interface-area range published a few
   lines below — `0.41–0.51 across 16×` — kept both of its halves from that row:
   0.41 was its error, and 16× was 3.2 → 0.2. The correction stood next to a
   number that contradicted it. Grep the value, and grep what it was divided by.
5. If an earlier correction claimed "every other conclusion is unchanged" and
   this one disproves it, **delete that sentence — do not reword it.** A
   paraphrase is the same claim: "the conclusions themselves are unchanged"
   replaced the original here and sat one line above the paragraph documenting
   the second number that moved. Say which numbers moved and stop. A correction
   is not the place for a reassurance whose scope nobody has checked.

A correction that only changes the number teaches nobody why it was wrong, and
this repository has now corrected its own corrections twice.

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
