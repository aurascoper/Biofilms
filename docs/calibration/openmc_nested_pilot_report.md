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
> effect**, whose numerator and denominator are each unbiased by construction
> (their ratio carries a bounded O(1/R) term — see below). The raw norm survives as a
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
`noise_floor` scenario returned 8 negative and 8 positive values across its 16
draws — a statistic that never went negative would still be biased.

**It is usable at R = 3.** Simulation at this geometry recovers a true effect of
0.0200 at every R from 3 to 20, where the raw norm on the same draws reads
0.0223–0.0294. More replicates buy *precision*, not a correction for the raw
norm's bias — which is what supersedes the history recommendation below.

**But the ratio is not itself unbiased, and an earlier version of this report
said it was.** The numerator and denominator are each unbiased; their ratio is
not, because `E[X/Y] ≠ E[X]/E[Y]` and the two are correlated through the
baseline rows they share. The right thing to record is the scale, not a point
estimate that depends on an unmeasured correlation:

    |bias| / SD  ~  relerr · √( 2(1−ρ) / (R · B_eff) )    < 3 × 10⁻²

across the whole envelope this pilot operates in. `B_eff` is the participation
ratio (Σwμ²)²/(Σw²μ⁴) ≈ 400 — **not** the bin count ≈ 1500, since a dose field
concentrated near a line source has far fewer effective bins than bins. No
plug-in correction is applied: both the delta-method and jackknife corrections
were measured and both are harmful at small R, turning −35% into +170% in the
stress regime. A correction whose own variance exceeds the bias it removes makes
the number worse.

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

## The gate, after an external review found three faults in it

Thresholds are declared in **E²**, so δ = 0.10 becomes 0.01 and the practical
floor 0.02 becomes 4 × 10⁻⁴. `ThresholdPolicy` carries a `metric_id` and
`decide()` refuses draws from a different metric.

An automated review of the merged branches then found three further faults, two
of which reached the verdict itself. All three are fixed and the run below is
regenerated.

**The gate was pooling draws across tally resolutions.** That put a systematic,
monotone discretisation trend inside a distribution the gate reads as sampling
spread — and `mesh_coarsening_factor` is declared `sampling_role =
convergence_axis` precisely so it is never drawn from. The trend is not small:
`clearly_detectable` runs 5.05 × 10⁻² at f = 1 to 8.09 × 10⁻² at f = 5, a **60%
drift and ~71× the correct standard error**.

Reading only the finest factor would have been a *relaxation*, not a fix:
`lever_density_x1.35`'s 99% upper bound falls 7.59 × 10⁻⁴ → 7.12 × 10⁻⁵, under
the practical floor, converting "real but not worth acting on" into "below
threshold" **by looking at less**. So `decide()` now takes the residual
discretisation bound and widens the interval symmetrically. Under finest+bound
no verdict became easier and three became stricter.

**The significance flag was wrong three ways in one line.** It pooled over draws
and factors; `sqrt(mean(jackknife_var))` is the RMS standard error of *one row*,
carrying no 1/√M and no between-draw term, so it was not the standard error of
the mean it divided into; and "3 σ" on M−1 = 3 degrees of freedom is a ~94%
one-sided test wearing a 99.9% label. It is now a cluster-on-outer-draw SE at
the finest mesh with the correct t critical value of 10.21.

The direction here is worth recording, because the first analysis got it
backwards: measured against the cluster SE the old statistic was **conservative
in 10 of 11 scenarios** (2.06× to 10.15×) and anti-conservative in exactly one
(0.87×). It produced no false positives — it simply was not measuring what its
name claimed.

**Ω_b required a dose floor taken from replicate 0.** The same error this module
had already learned once with `rel_err`, in a subtler form: the region became a
random set drawn from one replicate, so 2(R−1) of the R(R−1) ordered pairs the
U-statistic averages involve a row that helped choose the weights — four of six
at R = 3. Ω_b is now two geometric criteria and takes no field at all.

### The regenerated run

| Scenario | Returned | E² | E | t | vs zero |
|---|---|---|---|---|---|
| `noise_floor` | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 1.05e-04 | 0.0103 | 0.4 | **no** |
| `lever_Fe_5pct` | `EFFECT_BELOW_THRESHOLD` | 6.49e-06 | 0.00255 | 10.3 | yes |
| `lever_density_x1.35` | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 5.52e-05 | 0.00743 | 1.8 | **no** |
| `lever_Gd_5pct` | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 1.34e-03 | 0.0366 | 39.9 | yes |
| `lever_dehydration_H_0.06` | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 2.12e-03 | 0.0461 | 15.7 | yes |
| `near_threshold` | **`INDETERMINATE_NUMERICS`** | 1.50e-02 | 0.1224 | 81.0 | yes |
| `clearly_detectable` | `PASS_SYNTHETIC_GATE` | 5.05e-02 | 0.2248 | 84.0 | yes |

1 056 OpenMC runs, 1.056 × 10⁸ histories, 610 s wall. Ω_b is 1 631 bins, the
geometric count.

**`near_threshold`'s PASS is withdrawn.** Its interval is −0.0035 to 0.0336
against a threshold of 0.01, with a resolution bound of 0.0182 — it straddles,
and numerics dominates. The earlier PASS was carried by coarse meshes whose
spread exceeds the threshold it passed, so it was never resolution-converged.
The honest response is more factors near f = 1, not averaging the axis away.

**The lever magnitudes are unchanged by all of this** — every value moved by
≤ 0.7% when the dose floor was dropped — so the physical ranking and every
conclusion drawn from it stand. What changed is which of them the gate is
entitled to call resolved.

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
carries no history-dependent bias to buy down at R = 3 and at the pilot's
existing 10⁵ histories per run, and the
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

## The gate and the significance test can disagree, and a run that does cannot publish

They are computed independently — `decide()` reads quantile bounds,
`distinguishable_from_zero` runs a t-test — and nothing reconciled them. Two draws
of 0.05 with se = 0.01 return `PASS_SYNTHETIC_GATE` **and**
`distinguishable_from_zero = false`: the run asserts the effect cleared the
threshold and, one field later, that it is not separable from no effect.

`publication_block()` refuses that combination at publication time.
**Refusing is the point: it does not pick a winner.** Downgrading the PASS would
settle the disagreement for the t-test; keeping it settles it for the bounds.
Neither is the pipeline's call. The row is kept in
`openmc_nested_pilot_verdict.json` with the contradiction **named** in a
`publication_block` field — flagged, not retired — and the canonical tables are
refused. It fires in one direction only: a non-PASS verdict that is *not*
distinguishable is coherent, and `noise_floor` (t = 0.369) and
`lever_density_x1.35` (t = 1.840) are exactly that. A symmetric check would
refuse most of the table above. Checked against all 11 scenarios of the run
above: it blocks none of them.

**The minimum outer-draw count was raised from 2 to 4, and the old number was
derived from the wrong property.** It came from asking where `t_critical_999`
stops returning `None` — the function's *domain*. A floor needs its usable
*range*: at df = 1 it returns 318.3, a value and not an attainable one, so two
draws can PASS and can essentially never be significant. The replacement reads
the shape of the curve instead:

| draws | df | t_crit(0.999) | ratio to next |
|---|---|---|---|
| 2 | 1 | 318.3 | 14.3× |
| 3 | 2 | 22.33 | 2.19× |
| 4 | 3 | 10.21 | 1.42× |
| 5 | 4 | 7.173 | |

m = 4 is the first count at which one more draw no longer more than halves the
bar. This is a knee-of-the-curve heuristic about design stability, **not a power
calculation**, and it deliberately does not ask which of the findings above
survive — a floor chosen to protect existing results is chosen after seeing them,
and its margin here would have been 1.2%.

**α = 0.999 one-tailed is half of why the floor is needed.** At 0.95,
t_crit(df = 1) is 6.31 rather than 318.3. Raising the floor and lowering α would
address the same symptom; α is a declared decision policy and is recorded here
unchanged rather than moved as a side effect.

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
