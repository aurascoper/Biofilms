# V3: specify synthetic survival and defer phase execution

This specification separates loss of viable parcels from movement. It does not change the reported movement result or claim biological calibration.

Status: corrected design, 2026-09-29; **not a production preregistration**. The accompanying patch validates physical dose and explicit response inputs. It implements no hazard, growth or phase stepper. The production sensitivity vector, analysis implementation and run manifest must be committed in a later PREDICTION before those experiments start. The existing coefficient, horizon and clustering studies retain their own protocols.

Corrections to the V2 note are recorded in ledger row SYN-NOTE-01. The original note is retained as a historical input; its SHA-256 is `0f45128a3958a30f36b4a8248ce7d73d30cbab159b50a386dbcba17c1eb978fa`. The initial review and numerical checks are in Hunter's Downloads folder, `Exercised_fitness_phase_v2_review`. The model reviewed there was `66ffe9f`; the PRRT specification below is pinned to `a84e9cc`.

## 1. Three outcomes answer different questions

Founders have persistent parcel IDs; this extension has no growth or division. For species s, let N_s(0) be its initial parcel count, A_s(t) its viable IDs, and r_i(t) the parcel centre-of-mass radius. Use the existing cylinder radius R.

| Outcome | Definition | Interpretation |
|---|---|---|
| Viable fraction | number of viable IDs / N_s(0) | Loss of synthetic survival units |
| Mean radius of viable parcels | mean of r_i(t)/R over A_s(t), equal parcel weights | Composition and position of the remaining population |
| Displacement of identified viable parcels | mean of [r_i(t) − r_i(0)]/R over the same A_s(t) | Actual radial displacement among the remaining IDs |

Also report the initial radius of those same remaining IDs. The mean-radius change decomposes exactly into their mean displacement plus the difference between their initial mean radius and the initial mean radius of all founders. This exposes selection by location. It does not turn a survivor-conditioned displacement into an unbiased effect on all founders.

If A_s(t) is empty, viable fraction is zero and both radius statistics are missing, with `zero_viable` recorded. Never substitute zero radius or silently omit the seed. Report the number of evaluable seeds and extinct populations. A paired radius contrast uses the same declared set of evaluable seeds in both arms; missing arms prevent an all-seed radius verdict. Cross-arm matched-ID displacement uses only IDs viable in both arms and reports the intersection size and conditioning explicitly.

The prospective primary endpoint is MCS 100 at N = 40, six parcels per species, seeds 42–297. Report per-seed values before ensemble summaries. Exact outcome definitions, uncertainty intervals, missing-data rules and any decision threshold must appear in the committed PREDICTION and executable analysis before production. Ledger SYN-SURV-01 records the open question.

## 2. The physical clock and hazard are explicit assumptions

Each parcel is a **synthetic survival unit**. Applying a cell-level survival curve to whole-parcel removal assumes correlated loss within that parcel. A cells-per-parcel interpretation and appropriate experimental endpoints are required before making biological claims.

Choose the single-hit exponential family as a design assumption:

`p_i = -expm1(-alpha_s * delta_D_i)`.

Here alpha_s is finite and nonnegative in Gy⁻¹, and delta_D_i is the current step's occupied-site mean dose in Gy. This stable expression gives zero at zero exposure. Never reapply cumulative lifetime dose each step. The product of step survival probabilities is `exp(-alpha_s * sum(delta_D_i))` for a prescribed exposure history. Equal total exposure gives the same survival law, not identical individual Bernoulli draws under different step partitions.

Assuming this family, a D10 endpoint sets `alpha_s = ln(10)/D10_s`. A single endpoint does not uniquely select the family or justify a zero repair rate. Growth, repair and a quadratic dose term are excluded from this synthetic model by declaration.

The proposed primary clock is 60 seconds/MCS and the axis rate is 36 Gy/hour = 0.01 Gy/second. At 100 MCS the axis dose is 60 Gy. These are synthetic design choices. They do not calibrate the original mechanical clock; a growth-free model cannot fit a doubling time it does not compute. The baseline `seconds_per_mcs = NaN` stays unset and active accrual still refuses it.

Install the physical field as `axis_rate * I(r)/I0`, with I0 positive and the original dimensionless field fixed. Do not alter the radiation acceptance signal or melanin drive by importing that field. The primary uses the mean dose map; its SD is validated and recorded, not sampled as an extra biological noise process.

`diagnostics/synthetic_response/config.jl` requires a seven-entry alpha vector in CN, DR, CS, BS, AN, SO, OI order, an explicit positive clock, and a finite nonnegative axis rate. It stores an immutable copy and supplies no biological defaults. A later manifest must give each coefficient's provenance or label it assumed. Zero is permitted for controls; it is not a substitute for an unknown species coefficient. No seven-species calibrated vector is established here.

### 2.1 A protected phenotype changes damage yield, not absorbed dose

The first mechanism study adds a dimensionless factor `q_i` to the hazard,

`p_i = -expm1(-alpha_s * q_i * delta_D_i)`.

The unprotected control has `q_i = 1`; a protected synthetic parcel has a declared `0 < q_i < 1`. The physical dose map, clock, direct radiation acceptance term and melanin field stay fixed. The factor has no biological default. A Dsup-inspired value remains assumed until it is measured in a matching host, exposure and endpoint. Whole-tardigrade survival, cultured-cell DNA damage and yeast oxidative-stress survival cannot calibrate a microbial parcel directly.

This contrast tests damage prevention only. TDR1 and TRID1 motivate a later explicit damage-and-repair state, while CAHS proteins motivate a separate osmotic-stress dormancy and recovery model. Biomolecular phase separation in those studies is not oscillator phase locking. The evidence boundary, controls and order of work are in `docs/research/tardigrade_protection_repair_spec.md`.

## 3. Controls isolate the declared paths

The experiment design must distinguish the physical dose map, the direct acceptance term and the radiation-derived melanin path. Use paired initial states and declared seeds. Keep adhesion, target volumes and all unvaried parameters fixed.

| Contrast | What changes | What it tests |
|---|---|---|
| Fixed-position exposure | Disable all lattice motion; apply the hazard to unchanged occupied sites | Radius can change by selective loss while identified-parcel displacement is exactly zero |
| Direct tropism on/off | Change only beta_ion; keep physical dose, alpha, melanin drive and production settings fixed | Effect of the direct acceptance term, including later exposure changes caused by motion |
| Physical geometry | Gradient versus a declared uniform dose map; keep acceptance signals fixed | Effect of dose heterogeneity, without demanding equal total loss or disappearance of mechanical sorting |
| Equal sensitivity | Set all alpha entries equal; keep the dose gradient | Spatial selection can persist without species-specific sensitivity |
| Disabled extension | Original driver versus the extension disabled | Identical lattice, fields, events, observables and mechanical RNG bytes |
| Zero exposure or sensitivity | Set the physical dose or alpha to zero | No hazard commitments, with the original mechanical trajectory preserved |
| Protected phenotype | Change only the declared damage-yield factor `q`; keep physical dose and alpha fixed | Damage prevention without changing absorbed dose or movement terms |

The uniform-map value and any common alpha must be fixed before results. A melanin on/off contrast, if added, is a separate arm; equal alpha does not isolate that path. Do not require exact stochastic extinction order across species. The fixed-position control measures viable counts without clearing or shrinking sites, so the geometry remains fixed by construction.

## 4. PRRT supplies invariants, not a transferable acceptance ladder

The ten gates below come from `spec/response_growth_survival.md` at PRRT `a84e9cc`, including amendments. A gate's status here concerns this survival-only analogue, not the correctness of the PRRT implementation.

| Gate | Status | Requirement here |
|---|---|---|
| G-N: zero dose | Applicable | Zero dose gives zero hazard commitments; also test zero alpha |
| G-O: obliteration | Adapted | Test a declared numerical large-alpha-times-dose limit and forced draws; no universal 1000-Gy threshold |
| G-P: protraction | Inapplicable | No repair or quadratic term; instead test equal-total-dose survival-law invariance |
| G-C: repair consistency | Inapplicable | No repair half-time or Lea–Catcheside factor |
| G-S: trait-selection null | Inapplicable | No inherited expression trait or refill; equal-alpha controls test a different question |
| G-D: trait drift | Inapplicable | No activity-allocation feedback or trait ratchet |
| G-Q: arrested cells exposed | Adapted | Every viable parcel receives the hazard independently of a division attempt |
| G-M: dying parent | Applicable | A parcel committed to loss cannot divide, even through the manual division API |
| G-B: cycle dose banking | Adapted | Use only the current increment, exactly once; no cycle/reset interpretation |
| G-H: activity hoarding | Inapplicable | No activity allocation among parcels |

Additional required tests, each with a defective input or implementation it rejects:

- Unequal successive increments and a moved parcel establish incremental, occupied-site dose attribution. A cumulative-dose reuse mutant must fail.
- A continuous run, split windows and save/restore continuation must agree exactly, including both RNG states and event order. A missing or altered hazard RNG in the checkpoint must fail.
- A fixed-exposure sampling control must recover its imposed alpha statistically. Proposed test: 100,000 independent synthetic units, alpha = ln(2) Gy⁻¹, total dose 1 Gy, hazard seed 9001. The known survival probability 0.5 must lie inside the two-sided 99.9% Wilson interval; transform that interval through `-ln(S)/D` for alpha. Require nondegenerate observed counts. A doubled-alpha mutant must be rejected. This is an implementation check, not biological validation; pin it before execution and report a failure without changing the seed.
- Registry/lattice consistency, commitment/removal distinction, the below-two-site rule and refusal of a dying division parent must each have a targeted defective case.
- For zero exposure and disabled mode, compare the actual mechanical RNG state as well as outputs. A hidden extra draw must be detected.

H2 in the Odin port assigns draws to its shared `s.rng`, protecting the separate replay stream `jr`. Here use a dedicated hazard generator with its own explicit seed and checkpoint state. That is a local design choice. Draw in ascending viable ID order, and skip draws for zero alpha or zero increment. Separate generators prevent direct stream consumption; paths and conditional draw counts may diverge after a commitment changes mechanics.

## 5. Commitment and removal are separate events

Use the resorption semantics described by PRRT Semantics B. A failed hazard draw marks a parcel `dying` and sets its individual target volume to zero. It stays occupied until resorption or the below-two-site cleanup. No further hazard draw or division is permitted for it. Do not replenish lost parcels.

As an explicit minimal mechanical assumption, dying parcels retain their species contact energies and direct/melanin acceptance coefficients. Only their target volume changes. From the following step onward they produce no melanin and consume no nutrient. Existing extracellular melanin continues to diffuse. This is a synthetic rule, not a claim about dead biomass chemistry.

Declare this future per-MCS order:

1. Validate configuration and state before stepping.
2. Apply mechanical proposals using the prior viable/dying states and targets; refresh parcel geometry and verify lattice counts.
3. Advance fields, with uptake and production restricted to viable parcels.
4. Accrue physical dose on the resulting occupied sites. Record the current increment per parcel. Occupied dying material may continue to accrue dose, but no hazard reads it.
5. Advance the physical clock to the interval end. Draw once for each eligible viable parcel in ascending ID order. Record commitments, including ID, species, MCS, time, increment and cumulative dose at commitment.
6. Remove dying parcels below two occupied sites; clear their lattice sites, archive their identity and record physical removal separately. Distinguish any ordinary zero-volume mechanical removal from hazard commitment. Refresh geometry before output.

Freeze the commitment-dose record. Keep occupied-material exposure until removal separately from viable survival accounting. Report viable counts, dying counts and occupied volume independently. A parcel still resorbing at the horizon is already nonviable; its removal time is censored, not invented.

The existing deletion-during-iteration allegation is **closed: not reproduced**. The supplied Julia test removed and logged all 1,000 zero-volume entries. This patch does not change that cleanup. Reopen only on new reproducing evidence; do not present precautionary refactoring as a demonstrated fix.

## 6. Checkpoints must reproduce the complete extension

A future extension checkpoint must include the model revision, Julia/RNG version, configuration and species order, seed declarations, both current RNG states, lattice and fields, dose maps and accumulators, per-parcel status/target/geometry, cumulative and commitment doses, clocks, next ID, archives, and separate event logs. Store a schema version and an explicit enabled/disabled mode.

Restore must validate before returning a runnable simulation. Missing extension state cannot silently mean “disabled.” Legacy checkpoints can be restored only as explicit baseline mode; they cannot resume an enabled extension. Do not silently change the existing transport-exchange schema. A normal baseline restart remains a separate test from an enabled-extension restart.

The present patch stages imports, both optional transforms, and computed accrual before writing. It checks exact shapes, finite nonnegative rates/SDs/sensitivities, and a finite positive active clock. Overflow in computed dose is refused. A window checks physical inputs before its first mechanical step. Windows are not transactions across multiple MCS; an error later in a window does not roll back earlier valid steps. Transform callbacks must not mutate captured state.

## 7. Correct the arithmetic without claiming new calibration

Under the assumed exponential family, `alpha = -ln(S)/D` and `D10 = D ln(10)/[-ln(S)]`.

| Quoted CN endpoint | Derived alpha, Gy⁻¹ | Derived D10 | 4300 Gy / derived D10 |
|---|---:|---:|---:|
| S = 0.5 at 120 Gy | 0.00577623 | 398.63 Gy | 10.79 |
| S = 0.2 at 360 Gy | 0.00447066 | 515.04 Gy | 8.35 |

These are conditional conversions, not fitted or arm-specific biological coefficients. [Khajo et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3178601/) attributes the approximate killing percentages to earlier work. The strain, pigmentation and endpoint provenance remain unresolved for a production coefficient. Ghosal's D10 table also credits earlier experiments; cite those experiments before calling a value primary.

As an arithmetic check, D10 = 12,000 Gy implies survival 0.825404 at 1000 Gy under this chosen family. At D10 = 70 Gy and axis dose 60 Gy, fixed axis/wall populations have expected survival 0.138950/0.483810. Equal initial counts give a radius weighted by expected survivors of 0.776881 R, from an initial 0.5 R, with zero displacement. This ratio of expectations is not the expected ratio for a finite population with possible extinction.

## 8. Phase execution remains deferred

The manuscript's phase term is unexercised and unimplementable as specified. Mutable parcel state can hold phases; the missing evolution law and coupling definitions are specification defects. No unrestricted claim about absence of precedents is required.

A possible **new** phase model, not a discretisation of the manuscript, would use synchronous updates from the old phase vector:

`theta_i(t+dt) = wrap(theta_i(t) + dt * [omega_i + sum_j A_ij K_s(i),s(j) sin(theta_j(t)-theta_i(t)) / d_i])`.

Here `d_i = sum_j A_ij`; isolated parcels have zero coupling. The graph, weights, units and step size must be declared before a phase experiment. Noise is absent in this candidate equation. A future noise process needs a separate stochastic definition.

For scalar positive coupling on a fixed graph, the normalized operator is `D_degree^-1 A`, not raw A. Its largest eigenvalue is 1 on a connected graph without isolated nodes. The [Restrepo–Ott–Hunt threshold](https://arxiv.org/abs/cond-mat/0411202) cannot be treated as an exact test for a small moving contact graph without checking its assumptions. Do not require strict monotonic synchrony in every finite run, or apply an infinite fixed-lattice conclusion as a finite-CPM verdict.

For independent uniform phases, `E[r²] = 1/N`: `1/sqrt(N)` is a root-mean-square value, not a floor. At N = 42 it is 0.154303; the large-N mean approximation is 0.136748. The familiar closed expression `sqrt(1-Kc/K)` is not a general order-parameter curve for arbitrary frequency distributions. A proposed Normal-frequency control must use its own reference calculation.

[Una and Glimm 2024](https://peerj.com/articles/16974.pdf) and PMC10909357 identify the same paper. Its phase increment is `omega * [1 + K mean(sin(delta_theta))]`; translate that normalization explicitly before comparing K values. Its adhesion depends on phase. [Hester et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3188485/) models chick somitogenesis, so the cited clock precedents are not all mammalian.

The bounded Appendix A wording is applied to the manuscript source in this branch. It identifies the absent fitness-field equation, missing phase closures and unimplemented Hamiltonian field source without claiming that a CPM cannot represent phase or survival. The stable review PDF remains on its separate branch.
