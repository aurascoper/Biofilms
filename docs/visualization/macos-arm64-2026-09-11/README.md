# ParaView views — macOS/arm64, 2026-09-11

First ParaView read of this exporter's output on Apple Silicon. Prior verification
(`Biofilms/docs/visualization/lattice_viewer.md:239`) was a Linux x86-64 6.2.0-RC1 binary.

- **ParaView 6.1.1**, native arm64, `brew install --cask paraview`, macOS 26.6.2, Apple M4.
- **Data:** `../evidence/biofilm-signal-evidence/inert_signal/paraview/signal_trajectory.pvd`
  (run `inert_signal_seed42_v2`, parent `manuscript42_seed42_8cbb8ee`).
- **Nothing was rerun.** No Julia, no repository edit, no figure regeneration. Both protected
  CSVs and the authoritative figures are untouched. `tools/render_label_trajectory.py` was
  **not** invoked — the parent's `status` is `manuscript_built`, the exact state its standing
  P2 overwrites.

## Verified on read

ParaView reported 101 timesteps, 0.0–100.0, and all ten cell arrays: `species`, `cell_id`,
`lineage_id`, `generation`, `interior_mask`, `radiation_cpm`, `melanin`,
`accumulated_dose_Gy`, `signal`, `occupied_above_threshold`. `Spacing 1.0 1.0 1.0` and
`WholeExtent 0 40 0 40 0 40` — the lattice-unit pitch refusal survives into the artifact.

All **209** artifact hashes in `derived_manifest.json` were recomputed: **209/209 match**,
zero mismatches, zero missing. `parent_manifest_sha256` and the `derived_manifest.sha256`
receipt both match.

## The slice plane is data-chosen, not geometric

**The z=20 mid-plane cannot show the threshold onset.** The 3D maximum sits near **z ≈ 28**
throughout.

**Coordinate convention, because two mid-plane numbers circulate and both are right.** Cells
span `[k, k+1]`, so an integer z is a cell *face*, not a cell *layer*:

| MCS 8 selection | max |
|---|---|
| VTK cutter plane at `z = 20.0` (cuts between layers 19 and 20) | 4.5761 |
| zero-based cell layer `k = 20`, centre `z = 20.5` | 4.4515 |

Both are below the 5.0 threshold, so the conclusion holds either way. This is also why the
chosen slice is written `SLICE_Z = 28.5` and not `28` — 28.5 is the centre of layer 28.

| MCS | 3D max | at (x,y,z) | sites > 5.0 | z=20 plane max | z=28.5 slice max |
|---|---|---|---|---|---|
| 7 | 4.8976 | (29,29,12) | 0 | 4.0691 | 4.8778 |
| 8 | 5.2971 | (26, 8,28) | 174 | 4.4515 | **5.2971** |
| 30 | 8.4933 | (24, 9,29) | 3052 | 7.0098 | 8.3851 |
| 100 | 8.9034 | (23,10,28) | 3384 | 6.9602 | **8.9034** |

**z = 28.5 is a plane selected after inspecting the trajectory**, because it shows the onset
and the final peak — at MCS 8 and MCS 100 its maximum equals the 3D maximum exactly. It is
**not** canonical and does **not** contain every timestep's global maximum: at MCS 30 the 3D
maximum is 8.4933 at z = 29, while this slice gives 8.3851. Any figure using it must say it
was chosen this way. `signal_z20_mcs*.png` are kept as the documented negative case.

## Files

| File | What |
|---|---|
| `build_views.py` | pvpython driver. Palette taken verbatim from `Biofilms/viewer/paraview_species.py`. |
| `verify_state.py` | Reopens the `.pvsm` in a **fresh** process and asserts 20 properties (`--receipt` pins 2 of them). Exits non-zero on failure. |
| `verify_artifacts.py` | Checks the data tier: pvd timestep set, referenced files, all 209 hashes, and the frozen pins. `--receipt` required for a trustworthy result. |
| `make_render_manifest.py` → `render_manifest.json` | Freezes the expected hashes and the full rendering recipe, beside the renderers. |
| `controls.py` → `control_verification.json` | The seven negative controls, runnable. |
| `biofilm_signal_views.pvsm` | Saved state, two views side by side. |
| `species_mcs{000,007,008,030,100}.png` | Parcel geometry: threshold `species` 1–7, Surface With Edges, categorical colours. |
| `signal_z28_mcs{...}.png` | Signal: threshold `interior_mask`=1, slice at z=28.5, `signal` fixed 0–9. |
| `signal_z20_mcs{...}.png` | Same at the geometric mid-plane — retained to show why it was rejected. |
| `build_report.json`, `state_verification.json` | Machine-readable receipts. |
| `html-viewer-dependencies.md` | External-dependency record for `biofilms-4d-viewer.html`. |

## Reopen check

`verify_state.py` loads the state in a fresh ParaView process: **20 checks, 0 failures**
pinned (18 unpinned) — reader path, file existence, the `.pvd` bound to content rather than
to a path suffix, 101 unique contiguous timesteps, species threshold 1–7 on `CELLS/species`,
`interior_mask` 1–1, slice origin z=28.5, categorical LUT with 7 annotated categories,
signal LUT 0–9 with `AutomaticRescaleRangeMode = Never`, and ray tracing off in both views.

An earlier revision of this file said "14 checks." That was true when written and stopped
being true when the script grew, with nothing tying the number to the code producing it.
Both counts here are now emitted by the run and copied from its output.

## The frozen receipt, and the clone that certified itself

`verify_artifacts.py` used to read `derived_manifest.sha256` from **beside the data**. So a
clone tampered *consistently* — bit flipped, manifest regenerated, receipt regenerated —
passed every check, because the manifest it was checked against was the one the tamperer
rewrote. Reproduced here: one bit flipped in `signal_mcs000050.vti` at byte 1506123, then
both manifest files regenerated, yields **0 of 209 artifact hashes altered** and "manifest
matches its own receipt: True."

`make_render_manifest.py` writes `render_manifest.json` **beside the renderers**, freezing
`derived_manifest_sha256`, `parent_manifest_sha256`, `pvd_sha256` and the artifact count,
alongside the rendering recipe (camera, slice as a full numeric coordinate, threshold
predicate, **adjacency rule**, LUT ranges, ray-tracing status, movie invocation) and hashes
of every exporter, renderer and output. `--receipt` pins both verifiers to it.

**The pin that catches the consistent tamper is not the obvious one.** That tamper leaves the
`.pvd` byte-identical — `f1a04725…` before and after — so pinning the `.pvd` alone misses it
entirely. Only the pinned `derived_manifest_sha256` bites. The `.pvd` reads like the index of
record; the manifest is the thing worth freezing.

Regenerating `render_manifest.json` is a deliberate act with its own timestamp. Regenerating
the manifest beside the data is what an attacker, or a careless rsync, gets for free.

**Protected CSVs.** `data/claims_ledger.csv` and `data/parameter_provenance.csv` are recorded
by `git hash-object` against `HEAD:`, not as a before/after pair. A before/after pair only
proves *this* pipeline was innocent; the blob comparison proves nothing has touched them since
the commit. Both `unmodified_vs_HEAD = true`.

## Staleness controls

`controls.py` — **runnable**, not narrated. An earlier revision of this file documented these
controls firing in prose; they had been run ad hoc in a shell and written up, with no script
anyone could re-run. In a repository whose own rule is *"a control that cannot fire is not a
control"*, a control that cannot be re-run is a claim, not a receipt.

```
python3 controls.py <inert_signal_dir> --receipt render_manifest.json \
        --with-state --state /path/to/biofilm_signal_views.pvsm
```

Every control runs on an APFS copy-on-write clone (`cp -c -R`); the evidence is opened
read-only and never written. Seven controls, 632 MB cloned per control, **12 s total**.
Results land in `control_verification.json` (`--out <dir>` to put them elsewhere).

**The two state controls need `biofilm_signal_views.pvsm`**, which is 550 KB of saved
ParaView state and is deliberately not committed alongside these scripts. Without
`--state` they report **BLOCKED**, and BLOCKED is counted as a non-pass: a control that
could not run has established nothing. The five data controls need only the evidence
bundle and run anywhere.

`verify_artifacts.py` checks the **data**; `verify_state.py` checks the **saved state**.
Pinned baselines against the real evidence: **16** and **20** checks, 0 failures.

Both now share one receipt contract: a supplied `--receipt` must load **and** carry the four
required pinned keys, or the run exits 2. `verify_state.py` previously loaded the receipt only
`if os.path.isfile(...)` with no `else`, and printed nothing about pinning — so a typo'd path
ran unpinned and produced output byte-identical to a deliberate unpinned run. Each run now
declares its mode: `pinned` (trusted-identity) or `self-consistency`, in the log and in
`state_verification.json`. Unpinned is 18 checks; pinned is 20.

Two adjacent repairs in the same pass: `verify_state.py` establishes every proxy's presence
**before** reading any property on it (it used to read `thr_sp.LowerThreshold` immediately
after the presence check and abort only later, turning a missing filter into an
`AttributeError` instead of a reported failure); and `animate_4d.py` now parses `config.toml`
as TOML and **refuses** a missing or undeclared threshold. Its regex form did
`if m and abs(...) > tol: raise` and then printed "threshold reconciled" unconditionally — so
a non-matching regex claimed a reconciliation it had not performed.

| # | Control | Result |
|---|---|---|
| — | baseline: unmutated clone | **0 failures** — the controls below are measured against this |
| 0 | `.pvsm` repointed at a decoy `.pvd` | **fires** — 4 failures (reader path, file existence, timestep count, frozen pvd hash) |
| 1 | timestep 50.0 removed from the `.pvd` | **fires** — 4 failures (100 entries, incomplete 0..100, self-consistency, frozen pvd hash) |
| 2 | timestep 50.0 duplicated | **fires** — 4 failures (102 entries, 101 unique of 102, self-consistency, frozen pvd hash) |
| 3 | one bit flipped in `signal_mcs000050.vti` | **fires** — 1 failure (208/209 hashes, names the altered file) |
| 4a | signal LUT narrowed 0–9 → 0–5 without re-render | **fires** — 1 failure (`signal LUT range 0..9 -> 0.0..5.0`) — see the correction below |
| 4b | ray tracing switched on without re-rendering | **not implementable against the state** — see below |
| 5 | **consistent tamper**: data + manifest + receipt all regenerated | **fires only when pinned** — unpinned **0 failures** (the hole), pinned **1 failure** (the fix) |

Control 5 is the one this harness exists for, and it is the only control whose *expected*
result includes a pass: unpinned it must pass, or the hole it documents is not real.

Controls 0, 1 and 2 report **4** failures where an earlier revision said 3. Not drift — the
pinned `.pvd` hash is a new failure mode catching the same tamper by a second, independent
route.

**A third, found later, and the most instructive: the LUT control fired without ever
performing its mutation.** It used a regex over the serialized state,
`(<Property name="RGBPoints".*?)(9)(\b)` with `re.S`. That is leftmost-first and lazy, so in a
1024-element Viridis table it matched the `9` in `<Element index="9" .../>` — an **index
attribute**, eight elements in — and rewrote it to `index="5"`, producing a duplicate index
and destroying index 9. The array was mangled, the verifier reported
`signal LUT range 0..9 -> 0.0..0.004874`, and the control counted one failure and looked
healthy. An earlier revision of the table above printed `0.0..5.0`, **a value no stored run
produced**; the committed receipt said `0.004874` the whole time.

The control now mutates through the parsed document, asserts the element it changed was at
`x = 9.0` before and `x = 5.0` after, and produces the `0.0..5.0` the label always claimed.
`sub_once` asserts a *count*; a count cannot tell you the edit landed on the right element.

**Control 5 also asserted only half its contract.** Its expectation is "unpinned PASSES, pinned
FAILS", but the verdict read `"FIRES" if nf_p else …` — the unpinned result fed a display
string and nothing else. A regression that broke the unpinned tier, or closed the hole the
control exists to document, would have left it green. It now requires `nf_u == 0 and nf_p`.

**Two harness bugs caught while building these, both worth recording.**

*Controls 1 and 2 initially reported "did not fire" and the verifier looked weak.* It was
not: the regex targeted `timestep="50"` while the file writes `timestep="50.0"`, so the
control silently no-opped and the verifier correctly passed an unmodified file. Every
control now carries `assert n == 1` on its own edit. **A control that cannot modify its
target is indistinguishable from a check that cannot fail.**

*The ray-tracing assertion was never load-bearing.* ParaView 6.1.1 **does not serialise
`EnableRayTracing` into a `.pvsm` at all** — verified by saving one state with it at the
default 0 and another with it explicitly set to 1, and grepping both: neither file mentions
the property. A reloaded state therefore always reports 0 no matter how the frames were
rendered, so the assertion would pass on a state saved with ray tracing on. It is kept,
**relabelled as documentation of an invariant rather than presented as a control**. The
render-time value is the actual evidence and lives in
`build_report.json["raytracing_props_zeroed"]`.

## macOS notes for the doc

The four ParaView-6.1.1-on-macOS issues and the absent `libopenvkl_module_cpu_device.dylib`,
each of which cost a run here, are stated **once**: in the `macOS / Apple Silicon, 2026-09-11`
section of `docs/visualization/lattice_viewer.md`. This heading has said "for the doc" since it
was written, and the account now lives in the doc rather than in both places.

They were duplicated across the two files until the copies had drifted — "Preset is" against
"The preset is", "both views" against "the views" — before either branch had merged. One
account and a pointer is the repair. Two accounts and an intention to keep them aligned is what
produced the drift.

The account sits there rather than here because `lattice_viewer.md` is already on the base both
branches share, so a pointer to it is at worst **early**: if this branch merges before #27, that
section is not in the file yet and the pointer does not resolve until #27 lands. A pointer the
other way would have been to a file that does not exist on #27's base at all, which is the
finding that started this.

The same environment fact also appears as a `note` field in `render_manifest.json`
(`make_render_manifest.py:126`). That is a receipt recording the conditions it was produced
under, not a third copy of the account, and it is left alone.


## 4D voxels, and a threshold definition that matters

`animate_4d.py` renders all 101 frames into two 3D voxel sequences with a **fixed camera and
fixed 0–9 colour scale**, so growth reads as growth and not as rescaling. `ffmpeg` encodes
each to MP4 at 12 fps.

| Output | What |
|---|---|
| `voxels_4d_parcels.mp4` | the 42-parcel consortium, categorical species colours, MCS 0–100 |
| `voxels_4d_quorum.mp4` | occupied voxels at or above the declared 5.0 threshold |
| `encode_movies.sh` | the exact ffmpeg invocation, so both MP4s are reproducible |
| `anim/connectivity_by_mcs.json` | connected-region count per MCS |
| `anim/anim_{parcels,quorum}_NNN.png` | the 202 source frames |
| `anim/anim_report.json` | per-frame voxel counts |

**Empty through MCS 7. At MCS 8, 174 voxels appear.** By MCS 100 there are 3,384.

**Cluster counts are computed, not eyeballed** — `anim/connectivity_by_mcs.json`, from a
ParaView `Connectivity` filter over the thresholded region, one value per MCS. An earlier
draft of this file said "five separate clusters" at MCS 8 and that the region "never merges."
The first half was **wrong** — counted off a single camera angle where one cluster occluded
another. The corrected trajectory is more interesting than the claim it replaces:

| MCS | connected regions |
|---|---|
| 7 | 0 |
| 8 | **6** |
| 9 | 13 |
| 11 | **18** (peak) |
| 15 → 100 | **10**, flat |

The region **fragments** first — 6 at onset, peaking at 18 by MCS 11 — then **consolidates**
to 10 and holds there for the remaining 85 frames. The "never merges into one body" half does
survive: the count never reaches 1 at any frame.

### The array to threshold on is `occupied_above_threshold`, not `signal`

First render thresholded `signal >= 5.0` directly. That disagrees with `metrics.json` on
**16 of 101 frames**, always by exactly +1.

Not a floating-point boundary convention — checked, and no voxel anywhere equals exactly
5.0. It is a **definitional** difference. `InertSignal.jl:82` states the endpoint as
*"`>= threshold` among currently occupied interior voxels"*, and `:86` implements
`count(occupied .& (A .>= p.threshold))`. Thresholding `signal` alone also catches
**unoccupied** interior voxels the field has diffused into, which is a different quantity.

The exported VTI already ships the right array. Summing `occupied_above_threshold`
reproduces `metrics.json` on **101/101 frames, zero disagreements**, so the quorum view now
thresholds that array and the rendered counts match the run's own endpoint exactly.

**The general lesson for any figure built off this bundle:** an exported field and a reported
endpoint can differ by a qualifier that lives only in the producer's docstring. Reproducing a
published number through a second, independent path is what surfaced it — the agreement check
was worth more than the render.
