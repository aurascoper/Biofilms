# Local prediction: synthetic survival, passive phase and growth

This prediction precedes production outcomes. The commit containing this file is the local freeze; there is no external registration receipt. `protocol.toml` fixes the configuration hash, seeds and horizons. Development used seeds 9001–9003 and analytical fixtures. Short unit tests also used seeds 42 and 43 (up to MCS 15); no primary endpoint was inspected and no parameter was tuned from those prefixes.

V1 executes the previously specified 256 seeds (42–297), N=40, six parcels per species, four arms and MCS 100 endpoint. No predecessor production result exists to replicate. V3 independently executes all four arms for 16 matched seeds (42–57), MCS 200. Growth OFF and ON are separate studies. The shared seed numbers pair initial conditions; they do not make these ensembles independent of one another.

The configuration is entirely synthetic. The seven hazard sensitivities are assumed, not drawn from species D10 evidence. The clock is 60 seconds/MCS and axis dose rate 0.01 Gy/s; the old mechanical model remains uncalibrated. The phase distribution has mean period 20 MCS, Normal frequency SD 0.05 rad/MCS and coupling 0.5/MCS. Frequencies and initial phases use a separate seeded stream; forcing and noise are off. V1 covers five nominal cycles and 50 coupling times; V3 covers ten and 100. The initial nutrient-limited doubling time is about 21 MCS, so V3 spans about nine initial doubling times. Nutrient depletion changes this timescale. The slowest imposed hazard has an axis half-life about 578 MCS: neither horizon establishes long-time survival. The horizons are finite observation windows, not asymptotic convergence claims.

Predictions: enabled hazard produces commitments in the ensemble; fixed-position selection can change survivor radius with exactly zero parcel displacement. Phase attachment changes no nonphase trajectory in either growth mode. Coupling may reduce phase-frequency spread inside connected components; global locking is not predicted on a fragmented moving graph. Selective removal and inherited daughter phases may increase coherence without dynamical locking. Growth transfers nutrient into biomass and permits division without violating the recorded source/sink budgets. Null or contrary scientific outcomes are retained; failed identity, budget or export checks invalidate that run.

Primary survival contrasts are survival_only minus baseline at the endpoint, per species and per seed. Report viable parcel count, founder survival, mean survivor radius, initial radius of those same survivors, and identified-parcel displacement. In V3, daughter displacement is relative to birth and is labeled so; founder survival cannot be read as population survival because division retires founders. Empty survivor radii are missing, never zero. Report evaluable and excluded seed counts. Seed bootstrap uses 5,000 samples, seed 9004, percentile 95% intervals; these are synthetic ensemble uncertainty, not biological intervals. No voxel or timepoint is an independent replicate. Shared nonphase factorial interactions must be exactly zero by the paired identities.

Phase inference uses the second half of each run, including both boundary observations (50–100 or 100–200). The cohort is the intersection of living IDs at every boundary in that window, minimum two. Every unordered cohort pair is tested, including nonedges. Finite-window locking requires absolute endpoint relative-phase drift at most 0.01 rad/MCS, relative-phase range at most 0.5 rad and zero full 2pi windings relative to the first window value. Evidence for an interaction contribution additionally requires at least 50% reduction in effective-frequency SD against the exact zero-K replay on the same recorded membership, initial phases, frequencies and inherited daughter phases. There is no common driver or phase noise to replay. This survivor-conditioned comparison is not a population effect or a test of biological synchronization. Record count, edges, connected components, isolates and per-component coherence with every global coherence. Report selection, birth and evolution increments separately; an extinction increment is missing. Between-arm phase contrasts use matched seeds, plus shared founder-ID diagnostics; ancestry does not count as independent synchronization. No four-arm phase interaction is defined when phase is off.

Controls use seeds 7000–7007 at MCS 100: fixed baseline/survival and moving survival with direct tropism on/off, all other paths held fixed. These are diagnostics, not extra primary arms. All eight are reported. Dose-zero, alpha-zero, disabled survival, phase free drift, two-node locking/drift, synchronous update, graph normalization, timestep refinement, selection-only, extinction, mass conservation, ID integrity and bad exports are covered by independent fixtures. Budgets must close within 1e-6 synthetic units; exact trajectory and restart checks allow no tolerance. Frequency and nutrient substeps stay fixed; analytical refinement is a separate development control.

Stop only after all declared runs complete or on a failed verification, reporting the failed case. Do not reduce seed count, horizon or change coupling after inspection. Diagnostic seed 42 exports every MCS, including t=0 and the terminal state; other seeds export 0, midpoint and terminal. Every MCS has compact metrics and canonical hashes. Every arm restarts at its midpoint and independently reproduces the remaining complete state. Parallel workers execute disjoint seeds with one Julia thread. The measured development growth benchmark took 49.97 seconds for four arms at 40 MCS including checkpoint replay and sparse exports; allow several hours total CPU time, eight workers, and tens of GB output. Rendering and VTK validation have separate receipts.

## Results (append after execution)

Pending. Preserve all text above when appending results.

Appended 2026-09-30 after complete execution. Nothing above this line changed; `git diff` on this file shows additions only.

### Receipts

- Execution source 841d052; configuration hash `eb090aa3c0a971776a6e262713eb0cf5c130b53cf2b8d7720ad93dec19f50b01`; analysis tool `analyse.py` at 732a57a (provenance note 1).
- Raw archive: `/home/aurascoper/Developer/Biofilms-sp-run/artifacts/fitness_phase/841d052/` (17 GB, this machine). `execution_jobs.json` lists 17 jobs, every return code 0. No `ABORTED.json`.
- Invariants: 272 of 272 `invariants.json` records are `exact` and `restart_exact`, one source, one config hash. `analyse.py` re-checked every seed's paired nonphase hashes on load.
- Exports: 1,088 of 1,088 `export_validation.json` passed (`verification/export_ensemble.json`).
- Views: 8 of 8 reopened from the relocated package path under ParaView 6.2.0-RC1; `reopen_receipt.json` sits beside each view in `views/growth_{0,1}_{arm}/`.
- Package: `/home/aurascoper/Developer/Biofilms-sp-run/artifacts/fitness_phase/841d052_review_package/` (2.3 GB), assembled by `package_review.py` at 732a57a; `SHA256SUMS` at its root.
- Analysis outputs: `v1_analysis.json` and `v3_analysis.json` (with `.png`) in both the archive and the package. Each analysis ran in 22 s.

### Provenance notes

1. `analyse.py` changed after the run, at 1158a0b. The change is confined to `plot_diagnostic`: the max relative-phase panel now takes the range over all cohort parcels instead of the difference from parcel 0. Every verdict below (`locking`, `coupling_supported`, cohort `count`, `frequency_spread_reduction`) is read from each arm's `phase_window.json`, which the Julia run wrote at run time. `analyse()` counts those values and recomputes none. No quantity named above passes through the changed lines. `analyse.py` was not edited during or after this analysis.
2. The 32 control runs (8 seeds x 4 cases) are metric-only: `controls.jl` writes one `controls.json` and no frames. They sit outside the 1,088 export count, which covers every run that exports.
3. Growth off (256 seeds, MCS 100) and growth on (16 seeds, MCS 200) are separate studies with separate horizons, as declared above. Nothing below compares one to the other.
4. Growth on has 16 seeds. A growth-on interval that includes zero reports the seed count, not the model.
5. The closeout plan named two relocated reopens. All eight were run, at 8 to 14 s each.

### Growth off: survival (survival_only minus baseline at MCS 100)

Radius is centre-of-mass distance from the axis over N/2 (0 axis, 1 wall). Displacement is the mean over living parcels of radius minus initial radius. Intervals are 95% seed-bootstrap (5,000 samples, seed 9004) over 256 paired seeds unless excluded. An excluded seed is one where survival_only had no living parcel of that species at MCS 100.

| Species (alpha, 1/Gy) | Alive | Survivor radius (evaluable / excluded) | Initial radius of those survivors | Displacement |
|---|---|---|---|---|
| 1 (0.01) | -1.74 [-1.88, -1.61] | +0.0086 [0.0019, 0.0153] (256 / 0) | +0.0075 [0.0003, 0.0145] | +0.0011 [-0.0025, 0.0047] |
| 2 (0.002) | -0.43 [-0.51, -0.35] | +0.0017 [-0.0017, 0.0052] (256 / 0) | +0.0004 [-0.0021, 0.0029] | +0.0013 [-0.0010, 0.0036] |
| 3 (0.015) | -2.62 [-2.76, -2.47] | +0.0192 [0.0103, 0.0282] (255 / 1) | +0.0210 [0.0111, 0.0305] | -0.0018 [-0.0061, 0.0026] |
| 4 (0.04) | -4.48 [-4.62, -4.34] | +0.0382 [0.0240, 0.0525] (207 / 49) | +0.0434 [0.0284, 0.0581] | -0.0052 [-0.0112, 0.0005] |
| 5 (0.02) | -3.18 [-3.32, -3.04] | +0.0191 [0.0095, 0.0290] (251 / 5) | +0.0212 [0.0102, 0.0320] | -0.0021 [-0.0068, 0.0026] |
| 6 (0.06) | -5.26 [-5.35, -5.16] | +0.0493 [0.0289, 0.0704] (145 / 111) | +0.0538 [0.0327, 0.0740] | -0.0045 [-0.0117, 0.0027] |
| 7 (0.03) | -3.82 [-3.97, -3.67] | +0.0356 [0.0221, 0.0489] (240 / 16) | +0.0359 [0.0219, 0.0500] | -0.0003 [-0.0056, 0.0049] |

Hazard commitments: 5,511 in survival_only and 5,511 in combined, out of 10,752 parcels; 0 in baseline and phase_only. Founder survival equals alive in this mode, because nothing divides. The survivor-radius shift is positive for every species, and its interval excludes zero for six of seven. The initial-radius shift of the same survivors matches it within 0.006. Every displacement interval includes zero. So the survivor population sits farther from the axis because parcels nearer the axis died, not because survivors moved. This is the prediction's first case: selection changes survivor radius with no measurable displacement.

### Growth off: phase

| Arm | Cohort evaluable | Cohort size (median, min) | Finite-window locking | Coupling supported | Spread reduction vs zero-K replay (median, range) | Endpoint coherence | Components | Isolates |
|---|---|---|---|---|---|---|---|---|
| phase_only | 256 / 256 | 42, 42 | 0 / 256 | 0 / 256 | 0.492 (0.172 to 0.818) | 0.411 [0.388, 0.435] | 11.3 | 5.2 |
| combined | 256 / 256 | 21, 12 | 0 / 256 | 0 / 256 | 0.350 (0.098 to 0.807) | 0.379 [0.358, 0.400] | 9.7 | 5.4 |

No seed meets the finite-window locking criterion in either arm. The spread-reduction bound alone (at least 50%) is met in 121 of 256 phase_only cohorts and 48 of 256 combined cohorts. The pre-registered evidence criterion requires locking as well, so no seed supports an interaction contribution. Combined minus phase_only endpoint coherence is -0.032 [-0.060, -0.004] over all current IDs and -0.050 [-0.076, -0.024] over the same founder IDs. The selection increment in combined sums to +0.024 [0.012, 0.037] per seed; it is 0 in phase_only by construction. Reading: on this fragmented moving graph (about 10 components, 5 isolates), coupling narrows the frequency spread inside components but produces no locking within 50 MCS, and mortality lowers coherence slightly through membership. Null result retained.

### Growth on: survival and growth (survival_only minus baseline at MCS 200, 16 seeds)

Budgets: max abs nutrient residual 7.2e-9, max abs biomass residual 7.2e-9, against tolerance 1e-6. Divisions: 5,227 in baseline and 5,227 in phase_only; 3,823 in survival_only and 3,823 in combined. Hazard commitments: 1,485 in survival_only and 1,485 in combined. Division retires founders in every arm, so founder survival is 0 or -0.06 for every species and is not a population survival fraction.

| Species (alpha, 1/Gy) | Alive | Biomass (synthetic mass units) | Survivor radius (evaluable / excluded) | Displacement |
|---|---|---|---|---|
| 1 (0.01) | +1.6 [-1.4, +4.8] | +6,142 [1,021, 11,415] | -0.005 [-0.099, 0.086] (8 / 8) | -0.044 [-0.100, 0.018] |
| 2 (0.002) | +9.4 [5.3, 13.8] | +34,137 [27,900, 40,681] | -0.052 [-0.115, 0.012] (14 / 2) | -0.057 [-0.119, -0.005] |
| 3 (0.015) | -4.3 [-7.0, -1.6] | -3,002 [-8,424, 2,098] | +0.009 [-0.073, 0.090] (11 / 5) | -0.009 [-0.065, 0.052] |
| 4 (0.04) | -8.4 [-10.4, -6.1] | -13,293 [-17,223, -9,449] | one evaluable seed (1 / 15) | one evaluable seed |
| 5 (0.02) | -5.2 [-7.8, -2.7] | -6,032 [-10,204, -1,805] | +0.017 [-0.059, 0.089] (9 / 7) | -0.036 [-0.108, 0.039] |
| 6 (0.06) | -6.3 [-8.5, -4.0] | -10,348 [-14,553, -6,117] | missing (0 / 16) | missing |
| 7 (0.03) | -6.5 [-8.8, -4.1] | -10,264 [-14,877, -5,774] | two evaluable seeds (2 / 14) | two evaluable seeds |

Radius and displacement here are over living parcels of any generation; daughter displacement is relative to birth. They are mostly missing because one arm often has no living parcel of that species at MCS 200, and they carry no verdict. Species 1 and 2, the two lowest assumed sensitivities, gain parcels and biomass under the hazard relative to baseline. The other five lose both. The design identifies no mechanism for the gain. Nutrient freed by the deaths of other species is one reading, and it is untested here. At 16 seeds the species-1 alive interval includes zero and reports the seed count.

### Growth on: phase

| Arm | Cohort evaluable | Cohort size (median, min) | Finite-window locking | Coupling supported | Spread reduction (median, range) | Endpoint coherence | Components | Isolates |
|---|---|---|---|---|---|---|---|---|
| phase_only | 16 / 16 | 46.5, 26 | 15 / 16 | 15 / 16 | 0.982 (0.940 to 0.992) | 0.994 | 1 | 0 |
| combined | 16 / 16 | 26, 17 | 15 / 16 | 15 / 16 | 0.982 (0.955 to 0.994) | 0.994 | 1 | 0 |

Seed 44 fails the range bound only, in both arms: relative-phase range 0.755 rad (phase_only) and 0.848 rad (combined) against the 0.5 rad bound, with drift within 0.0085 rad/MCS of zero and no windings. Combined minus phase_only endpoint coherence is -0.0003 [-0.0023, 0.0015]. The same-founder contrast is undefined because no founder survives division. The zero-K replay shares membership, frequencies, initial phases and daughter inheritance, so the 98% spread reduction is attributable to coupling on this graph and not to inheritance. Reading: with growth the lattice fills, the contact graph is one component with no isolates at MCS 200, and K = 0.5 locks it within the 100 to 200 window. This is a synthetic-parameter result. It says nothing about biological synchronisation and nothing about the fragmented growth-off graph.

### Controls (seeds 7000 to 7007, MCS 100, all eight reported)

| Case | Alive at MCS 100 per seed | Hazard commitments per seed | Max abs displacement | Mean radius change |
|---|---|---|---|---|
| fixed_baseline | 42 in every seed | 0 | 0 | 0 |
| fixed_survival | 19 to 24 | 18 to 23 | 0 | 0 |
| moving_tropism_on | 19 to 23 | 19 to 23 | 0.155 | -0.028 |
| moving_tropism_off | 19 to 23 | 19 to 23 | 0.159 | -0.029 |

Fixed-position survival removes 18 to 23 of 42 parcels per seed with exactly zero displacement. Turning the direct tropism term off changes the alive count in one seed by one parcel and leaves displacement and radius change within 0.004 of the on case.

### Limitations

Every parameter is synthetic; no biological calibration is claimed or supplied. Horizons are observation windows. Growth on has 16 seeds. The two growth modes are separate studies. Phase never feeds back into any nonphase quantity, so no phase-to-fitness claim is available from this design. The original movement result is untouched.

## Execution amendment: deterministic removal records

The first batch at `ae67192` was stopped when the growth arms failed restart equality at MCS 101 (seeds 44, 46 and 48). A field-by-field probe of seed 44 found identical lattice, fields, parcel values and RNGs; only the order of mechanical-removal events in the death log and archive differed. Julia serialization can rebuild Dict buckets in a different order after division. The default legacy cleanup remains unchanged. The enabled extension now requests cleanup in ascending parcel-ID order. This is an event-order correction, not a failure to remove entries while iterating. A regression replays the actual failed seed across MCS 100–103.

All outputs of the stopped batch remain under `artifacts/fitness_phase/ae67192`, with `ABORTED.json`; none enters outcome analysis. No primary outcome analysis was performed. The same configuration hash, seeds, horizons, observation windows, thresholds and export schedules will be rerun from the corrected source. This amendment preserves the original prediction rather than replacing it.
