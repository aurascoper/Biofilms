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

Four ParaView-6.1.1-on-macOS issues, each cost a run:

1. `GetLayout()` raises `No active view was found` in batch pvpython — no implicit layout exists. Use `CreateLayout()`.
2. `hasattr(view, "EnableOSPRay")` **raises** `NotSupportedException` rather than returning `False`; the back-compat helper intercepts the 5.7-era name. Probe only `EnableRayTracing`.
3. Preset is `Viridis`, not `Viridis (matplotlib)`.
4. Setting `Background` alone leaves scalar-bar and axis text on the previous palette — white on white. Use `LoadPalette(paletteName='WhiteBackground')`.

Also: `libopenvkl_module_cpu_device.dylib` is **absent from the macOS bundle**, so every run
prints `[openvkl] INITIALIZATION ERROR`. That is the OSPRay volume path, unused here — ray
tracing is off and both views are surface renders — but it is noise on every invocation and
means the macOS build cannot do OSPRay volume rendering as shipped.


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

| MCS | 26-conn (point-sharing) | 6-conn (face-sharing) |
|---|---|---|
| 7 | 0 | 0 |
| 8 | **6** | **7** |
| 9 | 13 | 14 |
| 11 | **18** (peak) | **19** (peak) |
| 15 | 10 | 13 |
| 16 | 10 | 12 |
| 17 → 100 | **10** | **10** |

The region **fragments** first — peaking at MCS 11 — then **consolidates** to 10 and holds
for the rest of the record. The "never merges into one body" half survives: the count never
reaches 1 at any frame.

Every number in that table carries its adjacency rule, because the two columns differ.
ParaView's `Connectivity` filter joins voxels touching at a corner; `scipy.ndimage.label`'s
default structure joins only faces. **The count stabilises permanently at MCS 15 under
point-sharing and MCS 17 under face-sharing** — an earlier revision of this file gave only
the point-sharing column and wrote "15 → 100: 10, flat", which is false under the other
convention.

### A flat count is not a stable set — so this was measured

`component_overlap.py` answers the question the counts cannot. Ten regions every frame is
equally consistent with ten *stable* regions and with regions being born and dying every
frame while the total happens to stay 10. It intersects each component's voxel set with
every component's voxel set in the next frame, inherits lineage by largest shared volume,
and records births, deaths, merges and splits.

| over the whole record | 26-connectivity | 6-connectivity |
|---|---|---|
| distinct lineages ever | 21 | 24 |
| births | 21 | 24 |
| **disappeared** (no overlapping successor) | **0** | **0** |
| **retired by merge** (overlapped, lost the claim) | **11** | **14** |
| merges | 10 | 13 |
| splits | 0 | 0 |
| alive at MCS 100 | 10 | 10 |
| of those, born before MCS 15 | **10** | **10** |
| churn within frames 16..100 | **0** | 3 |

**A correction, and it matters.** An earlier revision of this table reported "deaths after
MCS 15: 0" and the prose read "no region ever dies". `deaths` counted only predecessors with
**no overlapping successor**. In a merge both predecessors overlap, so neither was a death —
but only the largest claim is inherited, so the rest stopped existing and were counted
nowhere. **11 of 11 retirements (26-conn) and 14 of 14 (6-conn) were invisible.** The
conclusion survives and is now measured rather than inferred: nothing ever *disappears*, and
every lineage that stops does so by merge.

**The set is stable, not merely the count.** The same ten regions persist to the end under
both conventions, and consolidation is **monotone by merge** — nothing disappears, nothing
splits. The final state is adjacency-independent; the frame at which it is reached is not.
Under 26-connectivity the transition **into** frame 15 retires 3 lineages (13 → 10) and
nothing changes thereafter; under 6-connectivity the consolidation completes at frame 17,
which is the churn of 3 within frames 16..100.

**What this measures.** Persistence under **largest-overlap greedy inheritance** — not
material identity. Lineage integers from the two adjacency analyses are **not** the same
histories: with an identical final mask 6-connectivity refines 26-connectivity, so equal final
counts imply the same final partition, but say nothing about prior lineage assignments.

`component_overlap_v1_superseded.json` is retained as historical evidence of what was
published. Its `counts_by_mcs` are **identical** to the corrected run — the component counts
were never wrong; only the lineage accounting was.

`test_component_overlap.py` holds the accounting fixtures — single-frame lifespan, all-vanish
frame, two- and three-into-one merges, split, competing overlap, and a transition mixing
disappearance with merge. **26 checks.** Against the pre-repair implementation they report
`frames_seen == 1 -> 2`, `deaths == 1 -> 0`, and then a `KeyError` for a retirement counter
that did not exist.

`vti_read.py` is a read-only `.vti` reader in numpy alone — the exporter writes uncompressed
appended binary, so no VTK is needed. That is deliberate: an audit of the renderer's output
should not require the renderer's toolchain. It is validated against the run's own endpoint
before being trusted for anything new — summing `occupied_above_threshold` reproduces
`metrics.json` on **101/101 frames, zero mismatches**.

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
