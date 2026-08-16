# OpenMC nested pilot — variance, cost, and which lever actually moves dose

**Date:** 2026-08-16, revised twice — instrument, then metric · **Tier:** S0 ·
**`target_calibration = false`** · **Policy:** `synthetic_gate_contract_v1` ·
**Reference D verdict:** `NOT_EVALUATED` ·
**Driver:** `coupling/scripts/openmc_nested_pilot.py` ·
**Estimator:** `coupling/biofilm_openmc/feedback_uq.py` ·
**Gate input:** `data/calibration/openmc_debiased_effects.csv` ·
**Diagnostic:** `data/calibration/openmc_effect_samples.csv` ·
**Budget/verdict:** `artifacts/pilot/openmc_nested_pilot_{budget,verdict}.json`
(gitignored — the committed CSVs carry their own provenance headers)

> **Revision 2 — the metric itself was replaced.** Fixing the instrument (below)
> established that the primary metric had a positive bias under the null. It is
> now superseded as the gate statistic by a cross-replicate **debiased squared
> effect**, which is unbiased by construction. The raw norm survives as a
> declared diagnostic. Consequences: the weak levers shrink by 48–70% while the
> strong ones move by 1–4%, the noise-floor control is now *indistinguishable
> from zero*, and **the ~40× history recommendation from revision 1 is itself
> superseded** — see "Correcting the correction".
>
> **Revision 1 — the instrument.** An adversarial re-reading of the pilot's own
> code found four defects in the *instrument*, not in the transport. Two headline
> numbers were wrong in kind. Set out in "What was wrong" at the end.

## The measurement that was missing: the metric's own noise floor

The metric that was primary through revision 1, and is now a declared diagnostic,
is an **unsigned** norm:

$$E=\frac{\sqrt{\sum_b m_b (D_1-D_0)^2}}{\sqrt{\sum_b m_b D_0^2}},$$

so it cannot return zero even when the true effect is zero: residual Monte Carlo
noise enters as positive bias that never cancels. **How large that bias is was
never measured**, because the old `zero_effect` control ran both material states
on *identical material with identical seeds*, making the two models bit-identical
and forcing the metric to ~10⁻¹⁶ by construction — at any energy, for any
physics, however broken the transport. That is a determinism check, not a
false-positive control.

`noise_floor` replaces it: identical material, **decorrelated seeds**. What it
reports *is* the floor.

| histories / run | f = 1 | f = 2 | f = 4 | f = 5 |
|---|---|---|---|---|
| 1 × 10⁵ (pilot default) | **0.0545** | 0.0316 | 0.0151 | 0.0167 |
| 4 × 10⁵ | 0.0279 | 0.0151 | 0.0066 | 0.0085 |
| 1.6 × 10⁶ | 0.0132 | 0.0074 | 0.0031 | 0.0040 |

Successive ratios are 1.96 and 2.10 against the 2.0 expected from 1/√N, at every
mesh factor. **The floor is pure Monte Carlo noise and scales as predicted**, so
it can be projected rather than guessed.

At the pilot's own settings the floor at the finest mesh is **0.0545** — larger
than four of the six material levers the first version of this report tabulated
as measured effects.

### Why the budget's `transport` term did not catch this

The budget's `transport` entry is the variance of the effect *across replicates*,
5.2 × 10⁻⁶ for `near_threshold` — a standard deviation of 0.0023. The noise floor
is 0.0545, **24× larger**. They disagree because they measure different things:
replicate scatter is how much the estimate *moves*, and the floor is how far it
sits from zero *on average*. An unsigned metric has a noise **bias**, and the
pilot was reporting only its **scatter**, which is small and stable precisely
because the bias is systematic.

*(Those figures describe the budget as it stood in revision 1. The budget is now
denominated in E² throughout, and its `transport` term is the mean
delete-one-replicate jackknife variance of the debiased statistic — 7.7 × 10⁻⁷
against a discretisation bound of 9.2 × 10⁻⁴. The section below on the estimator
is what made the old comparison unnecessary.)*

## The estimator that removes the bias instead of measuring it

Knowing the floor is 0.0545 lets you *discount* an effect. It does not let you
*measure* one, because the raw norm's bias depends on the noise level, which
depends on the history count. The fix is an estimator with no bias to discount.

With `d_r = D₁,ᵣ − D₀,ᵣ` the paired difference at replicate `r`, and replicates
independent of one another:

$$\hat S=\frac{1}{R(R-1)}\sum_{r\neq s} d_r^{\mathsf T} W d_s,
\qquad E[\hat S]=\lVert E[d]\rVert_W^2 .$$

The cross terms carry no noise contribution, because `E[dᵣᵀ W d_s] = E[dᵣ]ᵀ W E[d_s]`
for `r ≠ s`. The diagonal `r = s` is exactly where the noise power sits, so it is
**excluded rather than subtracted afterwards** — and excluded by zeroing, not by
`G.sum() − tr G`, which would form and then cancel the large diagonal total in
precisely the near-null regime the estimator exists for.

The denominator gets the same treatment, because `E[‖D₀,ᵣ‖²_W]` carries the
*baseline's* noise power and a naive squared norm is biased high:

$$E^2=\hat S \Big/ \frac{1}{R(R-1)}\sum_{r\neq s} D_{0,r}^{\mathsf T} W D_{0,s}$$

Everything is compared against `δ²`. Taking a square root near zero would
reintroduce exactly the bias being removed.

**Ŝ may be negative, and must be.** An unbiased estimate of a non-negative
quantity has to fall below zero about half the time when the truth is zero. The
`noise_floor` scenario returned 9 negative and 7 positive values across its 16
draws — a statistic that never went negative would still be biased.

**It is already unbiased at R = 3.** Simulation at this geometry recovers a true
effect of 0.0200 at every R from 3 to 20, where the raw norm on the same draws
reads 0.0223–0.0294. More replicates buy *precision*, not correctness — which is
what supersedes the history recommendation below.

### What it does to the measured effects

| scenario (f = 1) | debiased E | raw E | change |
|---|---|---|---|
| `noise_floor` | **0.0103** *(not distinguishable from 0)* | 0.05449 | **−81%** |
| 5% Fe | 0.00253 | 0.00488 | −48% |
| density ×1.35 | 0.00743 | 0.02443 | **−70%** |
| 5% Gd | 0.03652 | 0.03786 | −4% |
| dehydration, H 0.112 → 0.06 | 0.04607 | 0.04707 | −2% |
| 20% Gd (`near_threshold`) | 0.12240 | 0.12372 | −1% |
| 40% Gd (`clearly_detectable`) | 0.22477 | 0.22595 | −1% |

**The correction is monotone in the effect size** — largest where the effect is
smallest — which is the signature of removing an additive noise floor rather than
of a physics change. The ranking is untouched, so every physical conclusion in
this report stands; what changes is how much of the weak levers was ever real.

## Common random numbers work, and that is now demonstrated rather than assumed

Every lever runs **both ways** — matched seeds and decorrelated seeds — so what
pairing buys is measured rather than assumed. Under the debiased estimator the
answer is unusually clean.

| Lever (f = 1) | paired E | unpaired E | jackknife sd (E², paired → unpaired) |
|---|---|---|---|
| 5% Fe | 0.00253 | 0.0117 *(consistent with 0)* | 2.35e-06 → 3.29e-04 (**140×**) |
| density ×1.35 | 0.00743 | 0.0116 *(consistent with 0)* | — |
| 5% Gd | 0.03652 | 0.03344 | 6.37e-05 → 5.33e-04 (**8.4×**) |
| dehydration | 0.04607 | 0.05524 | — |

**Common random numbers now move the precision and not the estimate**, which is
what a variance-reduction technique is supposed to do. Where the unpaired
measurement has the precision to say anything — 5% Gd — the two agree to 9%.
Where it does not, the unpaired value is statistically consistent with zero and
its apparent magnitude is noise.

Under the raw norm the same comparison showed paired and unpaired differing by
13–83%, and revision 1 read that as CRN "removing noise" from the estimate. It
was really the bias contaminating both by different amounts. The corrected
reading is stronger: pairing tightens the standard error by roughly one to two
orders of magnitude and leaves the answer alone.

## The finding: which lever moves dose, and by how much

Debiased, at **f = 1**, the finest mesh, paired, with rows behind every number in
`data/calibration/openmc_debiased_effects.csv`. `z` is the pooled estimate over
its own jackknife standard error — how far it stands from *no effect at all*:

| Material lever | Debiased E | z vs zero |
|---|---|---|
| `noise_floor` (identical material) | 0.0103 | **0.1 — not distinguishable from zero** |
| 5% Fe (replacing O) | 0.00253 | 7.9 ← the weakest real lever |
| density ×1.35 | 0.00743 | 4.2 |
| 5% Gd | 0.03652 | 46.1 |
| dehydration, H 0.112 → 0.06 | 0.04607 | 26.4 |
| 20% Gd (`near_threshold`) | 0.12240 | 60.6 |
| 40% Gd (`clearly_detectable`) | 0.22477 | 79.0 |

The ranking is **unchanged from every previous version of this table**. Every
lever is now resolved against zero on its own evidence, without reference to any
other scenario — including 5% Fe, whose 0.00253 is about **one twentieth** of the
raw norm's 0.0545 noise floor and was unmeasurable under the old metric.

`z` is deliberately not monotone in `E`: dehydration is the larger effect but
carries the larger standard error, so 5% Gd is the more firmly resolved one.
Magnitude and certainty are different questions, which is why the gate needs both
δ and a spread rather than a single test statistic.

**The physical conclusions are unchanged.** Dose is energy per unit mass, so a
uniform density change raises numerator and denominator together — heating
×1.306 against mass ×1.350, dose ×0.967 — and only second-order self-shielding
survives. At 661.7 keV Compton dominates and μ/ρ tracks Z/A ≈ 0.5 for nearly
every element, so swapping oxygen for iron changes almost nothing. What moves
the dose is a **departure from Z/A ≈ 0.5**: hydrogen at 1.0, or high-Z material
where the photoelectric term still contributes.

Two consequences carry forward. The strongest *physically realistic* lever is
the **hydrogen fraction**, which is exactly what the paired wet/dry gravimetry in
the Reference D protocol measures — an independent argument for the campaign's
priority ordering. And the metal-loading story is weaker than the project's
framing implies: 5% Fe is the least effective perturbation tested, and
`S. oneidensis` reducing iron does not by itself make a transport-visible
material change at these energies.

`clearly_detectable` is an extreme synthetic Gd 40% material with no biological
pretension; its only job is to prove the gate can pass. **The 0.10 threshold was
not adjusted to reach it**, which would have been a stop condition.

## The gate returned what was declared

Thresholds are declared in **E²** now, so δ = 0.10 becomes 0.01 and the practical
floor 0.02 becomes 4 × 10⁻⁴. `ThresholdPolicy` carries a `metric_id` and
`decide()` refuses draws from a different metric, because a threshold declared
for E and compared against E² is a silent factor-of-ten error that no inspection
of the number alone would catch.

| Scenario | Declared | Returned | E² | z |
|---|---|---|---|---|
| `noise_floor` | *(control)* | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 4.81e-05 | **0.1** |
| `lever_Fe_5pct` | *(measured)* | `EFFECT_BELOW_THRESHOLD` | 1.86e-05 | 7.9 |
| `near_threshold` | `DECIDED_BY_PILOT` | `PASS_SYNTHETIC_GATE` | 2.48e-02 | 60.6 |
| `clearly_detectable` | `PASS_SYNTHETIC_GATE` | ✅ same | 6.94e-02 | 79.0 |

1 056 OpenMC runs, 1.056 × 10⁸ histories, 520 s wall.

**A limitation this exposes, reported rather than patched.** `noise_floor`'s
verdict says *effect detected* on a scenario whose true effect is exactly zero
— and z = 0.1 says the opposite. Both are right about different things.
`decide()` takes the 99th percentile of the draws, and there are 16 draws (4
outer × 4 mesh factors), so that quantile is effectively the **maximum**. One
draw landed at +5.3 × 10⁻⁴, about 1.6σ, which is unremarkable noise but sits
just above the 4 × 10⁻⁴ practical floor.

So the verdict vocabulary is being driven by the upper tail of a very small
sample. This is the same class of problem as everything else in this report: a
quantile needs enough draws to mean anything, and `seed_sufficiency()` already
sets 20 as the floor for estimating a *variance* — a 99% quantile needs far more.
**The `z` column, not the verdict, is the trustworthy significance statement at
this sample size.** Raising the outer-draw count is the production fix.

## Numerics is a bound, not a variance

The old budget computed `np.var(by_factor, ddof=1)` over exactly **two** mesh
factors. A two-sample variance is algebraically $(d_1-d_0)^2/2$, so the published
3.8186 × 10⁻⁴ was a squared difference wearing a variance's name — it reproduces
$(0.15125673-0.12362119)^2/2$ to every digit. Dividing it by a genuine
multi-replicate transport variance produced the reported "70–115×", which
compared two different kinds of quantity.

This is also the project's own declared rule, not a matter of taste:
`mesh_coarsening_factor` carries `sampling_role = convergence_axis` in
`data/uncertainty/feedback_parameter_distributions.csv`, which forbids drawing
it and instructs carrying a **residual discretisation bound** instead.

The budget now runs **four** factors and reports the worst squared deviation from
the finest mesh, with `numerics_basis` and `numerics_n_factors` recorded beside
it. The effect rises monotonically under coarsening — it is systematic error, not
scatter about a centre, which is why the deviation is referenced to f = 1:

Debiased E by mesh factor — the effect rises monotonically with coarsening in
every scenario:

| Scenario | f = 1 | f = 2 | f = 4 | f = 5 |
|---|---|---|---|---|
| `near_threshold` | 0.1224 | 0.1506 | 0.1685 | 0.1821 |
| `clearly_detectable` | 0.2248 | 0.2596 | 0.2806 | 0.2844 |

Debiasing barely moves these — the noise bias mattered for the weak levers, not
for the strong ones — so the discretisation finding is independent of the metric
change, which is the best evidence that it is real.

Comparing sources on a common scale rather than in variance units, the
discretisation bound is 0.0304 against a transport standard deviation of
8.8 × 10⁻⁴ (both in E²) — a factor of **35**, and not the 115 that compared two
different kinds of quantity. Discretisation dominates scatter by more than an
order of magnitude; the A0 lesson holds. What does not hold is the inference
drawn from it, below.

## Cost, and a reversed recommendation

| Quantity | Measured |
|---|---|
| seconds per run | 0.50 |
| histories per second | 200 300 |
| statepoint per run | 58 952 B |
| raytrace (5 × 10⁵ rays) | 0.20 s, once per geometry |

**Hardware provenance:** AMD Ryzen AI 9 HX 370, OpenMC 0.15.3 CPU-only
(MPI/OpenMP), ENDF/B-VIII.0. Wall times are **machine-specific and not
portable**; the history counts are. No GPU throughput is reported because the
pinned OpenMC has no GPU path.

Projected full Reference D study — 128 outer draws × 2 material states × 20 seeds
× 2 mesh factors = **10 240 runs**:

| Histories per run | Serial wall | On 32 cores | Noise floor at f = 1 |
|---|---|---|---|
| 4 × 10⁶ (A0's setting) | ~57 h | ~1.8 h | ~0.009 |
| 1 × 10⁶ | ~14 h | ~27 min | ~0.017 |
| 1 × 10⁵ (this pilot) | ~1.4 h | ~3 min | 0.055 |

### Correcting the correction

This recommendation has now moved twice, so it is worth being explicit about what
each version got wrong.

**Revision 0** recommended the cheapest row: "with numerics at 97.7% of the
variance, quadrupling the histories buys a 0.7% component." That compared a
discretisation *bound* against a replicate *scatter* and omitted the noise
**bias** entirely.

**Revision 1** reversed it to ~40× histories, on the ground that the bias is what
decides whether a 0.02–0.05 effect is measurable. That was right *about the raw
norm*.

**Revision 2 supersedes it.** The debiased estimator has no bias to buy down. It
is unbiased at R = 3 and at the pilot's existing 10⁵ histories per run, and the
run confirms it: 5% Gd resolves at z = 46 and dehydration at z = 26 without
buying anything. **Histories now buy precision, not correctness**, and the
40× figure sized a problem that no longer exists.

What remains worth buying is different and smaller: histories tighten the
standard error as 1/√N, and the **discretisation bound still dominates the
budget** — 9.2 × 10⁻⁴ against a transport term of 7.7 × 10⁻⁷ in E², a factor of
1 200. The honest recommendation is therefore the **middle row for precision**,
with the real money in mesh resolution rather than in histories, and a
larger outer-draw count so the gate's quantile means something.

The general lesson is worth more than the number: **a metric with a bias makes
compute look necessary that a better estimator makes unnecessary.** Two of the
three recommendations in this document's history were sized against an artifact.

## What was wrong

Four defects, all in the instrument, all found by re-reading shipped code.

**1. The false-positive control was vacuous.** Seeds depended only on
`(draw, replicate)`, never on the material state, so `zero_effect` compared a
model with itself. Fixed by decorrelating the states in `noise_floor`; the floor
turned out to be 0.0545, not 10⁻¹⁶.

**2. "Numerics variance" was a two-point difference.** See above. Now a bound
over four factors, labelled as one.

**3. The material lever table was unprovenanced.** Six constants were hardcoded
into `_report()` and written into the budget artifact unconditionally, whatever
the run had actually done — no sample row, seed, replicate count or mesh factor
existed for any of them. They are re-measured under `--levers` and now have rows.

*A correction to the correction:* the apparent conflict between the dict's
`Gd_20pct = 0.1241` and the registry's `0.1384` was **not** fabrication. 0.1241
is the f = 1 value and 0.1384 the mean pooled over f = 1 and f = 2; reproduction
gives 0.1237 at f = 1. Both were right and neither said which it was, which is
the actual hazard of a number with no provenance. Reproduction also showed the
old `density ×1.35 = 0.0137` to be low by 78% against a measured 0.0244.

**4. Ω_b required `rel_err < 0.25`, a statistical cut on a region definition.**
It bit the mesh factors unequally — at f = 1 it removed roughly a third of the
bins (983/949/968 of 8 000, median rel_err 0.1837 against the 0.25 cut) while at
f = 2 it never bound (a constant 215/1 000, median 0.0837). Membership also
*moved with the history count*, so the metric's domain changed with its
precision. Removed; rel_err is now a reported diagnostic. Bin counts at f = 1
rose from ~965 to ~1 500.

The first report's claim that the 22% factor-1-to-2 spread was "genuine
discretisation rather than a masking artifact" was therefore **half right**: the
occupancy criterion had been equalised, the uncertainty criterion had not. With
both fixed the f = 1 → f = 2 movement is 0.1237 → 0.1509, still 22%, so the
conclusion survives — but it did not follow from the evidence given for it. It
survives the debiasing too, at 0.1224 → 0.1506.

**5. And then the metric itself.** Fixing the four above made the fifth problem
legible rather than solving it: the primary metric was still a norm, so it was
still biased under the null, and the best revision 1 could do was *measure* that
bias and discount against it. Measuring a bias is strictly worse than not having
one — the discount depends on the history count, so every comparison across
resolutions or budgets inherits it.

The pattern across all five is the same: **each defect was invisible to the
statistic that was supposed to catch it.** Replicate scatter could not see the
noise bias because the bias is systematic. A two-point spread could not see that
it was not a variance. A hardcoded table could not disagree with a run that never
happened. An estimator that could not fail its own control could not report that
it had. Nothing here was caught by the numbers looking wrong; everything was
caught by reading what produced them.

## Known residuals

- **Ω_b is volume-consistent across factors, not identical by construction.**
  Bin counts are 1 500 / 215 / 23 / 12 against the 187 / 23 / 12 that exact
  volume conservation would give — spot on at f = 4 and f = 5, 15% high at f = 2.
  Defining Ω_b once at the finest mesh and mapping it down would make this exact.
- **The dose floor still depends on the estimator**, through
  `nanmax(dose0)`. Bin counts at f = 1 vary 1 492–1 516 across draws, a residual
  1.6% against the ~33% the `rel_err` cut was causing.
- **The floor is measured unpaired.** It bounds the raw diagnostic, not the
  debiased statistic, which needs no floor.
- **The gate's 99% quantile is taken over 16 draws**, where it is effectively the
  maximum — see the gate section. `z` is the reliable significance statement
  until the outer-draw count rises.
- **`sobolev_smooth` exists and is barred from every gate.** It solves
  `(I − α∇²)u = f` exactly under zero-flux boundaries, numpy-only via an
  even-extension FFT, for display and for inspecting where a difference field
  concentrates. It must not enter a statistic: smoothing a difference field can
  both hide real localized structure and manufacture apparent structure out of
  noise, and α would become a second tunable with none of δ's discipline. A test
  asserts neither gate module references it.

## What this does not establish

No target parameter was inferred. Reference D's verdict is `NOT_EVALUATED`, the
offline feedback gate is untouched, online feedback stays disabled, and every
artifact carries `target_calibration = false`. The Gd 40% material is invented
and is not a claim about any organism.

The sensitivity ranking is measured **on the synthetic geometry at 661.7 keV**
and is not transferable to another source spectrum: the Z/A argument that makes
Fe uninteresting here would not survive a lower-energy source, where
photoelectric absorption scales steeply with Z. Note that the noise floor also
rises as energy falls, so a spectrum study needs the floor measured at each
energy before its effects can be read.

## Next

**Discretisation is now the whole of the problem.** It dominates the budget by a
factor of 1 200 over transport, and it is systematic rather than random: the
effect rises monotonically under coarsening at every scenario, 0.1224 → 0.1821
across factors 1 → 5 for `near_threshold`.

**"f = 1 is the finest tally available" is an implementation restriction, not an
OpenMC one** — an earlier version of this report stated it as though it were
physical. An `openmc.RegularMesh` sets its extent and its bin count
independently, and `biofilm_mesh_extent_cm()` already ignores the mesh dimension,
so a tally *finer* than the CPM material lattice is possible. What blocks it is
narrow: `resolve_mesh_dimension()` integer-divides and guards `factor >= 1`, the
biofilm base dimension is a literal `(n, n, n)` in `model.py`, and two packages
independently validate divisibility. The exact-CSG mass denominator
(`mesh_material_masses_kg`) is already resolution-agnostic and needs no change.

A subvoxel refinement study at ratios 1×, 2×, 4× on a fixed snapshot is therefore
the next measurement, and it answers a question the coarsening sweep cannot: how
much of the movement is tally discretisation inside a fixed piecewise-constant
geometry, as opposed to the geometry's own resolution.

**Snapshot resolution is a separate study** and must not be confused with it. A
naive N-sweep changes what a parcel physically *is*: `V_target = 120` in
`biofilms_potts.jl` is in **sites**, so raising N at fixed `V_target` shrinks
every parcel. The first version should rasterize one fixed analytic morphology at
several resolutions, with no CPM dynamics involved.
