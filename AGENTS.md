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

When you add a guard, add the input that proves it bites. When you review one,
ask what would have to be true for it to fail, and check that something makes
it so.

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

## Registering a question

The rule above is about a number this repository published and later found
wrong. This one is about a question it asked and never wrote down.

**An open question stated in prose gets a row in `data/claims_ledger.csv` when it
is asked, not when it is answered.**

It needs no new machinery, and inventing some would miss the point. `status =
needs_verification` is already one of the eight pinned verdicts and already means
exactly this; `required_to_fix` already holds what would settle it. RM-KR-07 is
the worked example, registering `mean r/R = 0.65` as open because "a multi-seed
run COULD resolve this one either way - which is why it stays
needs_verification" (quoted with its own hyphen). The rule is only that the row
is written at the moment of asking.

The case it comes from is MEL-N40-02. A write-up asked, publicly and in
falsifiable form, whether a melanin ordering that held 11 of 16 times at N = 20
would still hold at the preprint's own N = 40 six-parcel configuration. It asked
three times, in one HTML file, on one unmerged branch. A grep for its wording
over every local and remote ref returned that file and nothing else: no row, no
doc, no pull-request body.

**That is an absence claim with no scope wearing different clothes.** Its
durability was the durability of one file. Delete it, abandon the branch, or move
the post somewhere this repository does not track, and the question is gone with
no trace that it was ever asked, while the work it was asking about carries on.

It applies **regardless of the answer**. That question came back 16 of 16,
confirming the hypothesis it was asked about, and registering it was still right:
a question is worth recording while it is open, and a question nobody else can
find was never open to anyone else.

The asymmetry with the correction rule is the whole reason this one is separate.
A wrong number announces itself eventually — something recomputes it, or a guard
fires on it. An unregistered question announces nothing, ever. **Nothing
recomputes a question.**


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
