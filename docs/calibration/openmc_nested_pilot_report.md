# OpenMC nested pilot — variance, cost, and which lever actually moves dose

**Date:** 2026-08-16, revised after an instrumentation repair · **Tier:** S0 ·
**`target_calibration = false`** · **Policy:** `synthetic_gate_contract_v1` ·
**Reference D verdict:** `NOT_EVALUATED` ·
**Driver:** `coupling/scripts/openmc_nested_pilot.py` ·
**Artifacts:** `artifacts/pilot/openmc_nested_pilot_{budget,verdict}.json`,
`data/calibration/openmc_effect_samples.csv`

> **This revision supersedes the first version of this report.** An adversarial
> re-reading of the pilot's own code found four defects in the *instrument*, not
> in the transport. The physical conclusions largely survive — the lever ranking
> is unchanged and now has provenance and an independent cross-check — but two
> headline numbers were wrong in kind, and **the compute recommendation is
> reversed.** The corrections are set out in "What was wrong" at the end.

## The measurement that was missing: the metric's own noise floor

The primary metric is an **unsigned** norm,

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

## Common random numbers work, and that is now demonstrated rather than assumed

A paired effect below the unpaired floor is not automatically noise, because
matched seeds cancel common histories. Every lever therefore runs **both ways**,
and the two give independent estimates that can be checked against each other:
if the unpaired result is effect and noise in quadrature, then
$\sqrt{u^2-n^2}$ should reproduce the paired value.

| Lever (f = 1) | paired | unpaired | √(u² − n²) | agree? |
|---|---|---|---|---|
| 5% Fe | 0.00488 | 0.05490 | 0.00672 | ✅ |
| density ×1.35 | 0.02443 | 0.05061 | — (u < n) | n/a |
| 5% Gd | 0.03786 | 0.06486 | 0.03518 | ✅ |
| dehydration, H 0.112 → 0.06 | 0.04707 | 0.07565 | 0.05247 | ✅ |

Three of four agree within 35%. The fourth is uninformative rather than
contradictory: the density lever's unpaired value sits *below* the floor, which
is ordinary scatter in the floor estimate itself and leaves the quadrature
subtraction with nothing to subtract.

**So the lever values are real, and pairing is what makes them measurable.** The
first version of this report asserted them with no such support.

## The finding: which lever moves dose, and by how much

Measured at **f = 1**, the finest mesh, paired, with sample rows behind every
number in `data/calibration/openmc_effect_samples.csv`:

| Material lever | Effect (rel. L2) | × noise floor |
|---|---|---|
| 5% Fe (replacing O) | 0.00488 | 0.09 ← the weakest tested |
| density ×1.35 | 0.02443 | 0.45 |
| 5% Gd | 0.03786 | 0.69 |
| dehydration, H 0.112 → 0.06 | 0.04707 | 0.86 |
| 20% Gd (`near_threshold`) | 0.12372 | 2.27 |
| 40% Gd (`clearly_detectable`) | 0.22595 | 4.15 |

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

| Scenario | Declared | Returned | Effect (pooled median) | vs floor |
|---|---|---|---|---|
| `noise_floor` | *(new)* | `EFFECT_DETECTED_BUT_NOT_PRACTICALLY_IMPORTANT` | 0.02546 | — |
| `clearly_detectable` | `PASS_SYNTHETIC_GATE` | ✅ same | 0.26894 | 4.1× |
| `near_threshold` | `DECIDED_BY_PILOT` | `PASS_SYNTHETIC_GATE` | 0.15994 | 2.5× |

1 056 OpenMC runs, 1.056 × 10⁸ histories, 527 s wall.

**`noise_floor`'s verdict is the most informative line in the table.** On a
scenario whose true effect is exactly zero, the gate reports *effect detected,
not practically important* — a real false positive, above the 0.02 practical
floor. That is the demonstration the old control could never produce, and it is
why `above_noise_floor` is now emitted beside every verdict.

`near_threshold` still passes, and at 2.5× the floor that pass is defensible;
but it is no longer the comfortable margin the first report implied.

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

| Scenario | f = 1 | f = 2 | f = 4 | f = 5 |
|---|---|---|---|---|
| `near_threshold` | 0.1237 | 0.1509 | 0.1685 | 0.1821 |
| `clearly_detectable` | 0.2260 | 0.2599 | 0.2806 | 0.2845 |

Comparing sources in **effect units** rather than variance units, the
discretisation bound is 0.058 against a transport standard deviation of 0.0023
for `near_threshold` — a factor of **25**, not 115. Discretisation still
dominates scatter by more than an order of magnitude; the A0 lesson holds. What
does not hold is the inference drawn from it, below.

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

**The recommendation is now the upper row.** The first version of this report
recommended the cheapest, reasoning that "with numerics at 97.7% of the variance,
quadrupling the histories buys a 0.7% component". That reasoning was wrong,
because the 97.7% compared a discretisation *bound* with a replicate *scatter*
and omitted the noise **bias** entirely. The bias is 0.055 at the pilot's
settings, it scales as 1/√N, and it is the quantity that decides whether an
effect of 0.02–0.05 — which is where every physically realistic lever sits — can
be measured at all.

Histories are the thing to buy after all, and roughly 40× more of them than the
pilot used are needed to push the floor below 0.01. That is still a
laptop-overnight or few-core job; it does not need a cluster.

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
conclusion survives — but it did not follow from the evidence given for it.

## Known residuals

- **Ω_b is volume-consistent across factors, not identical by construction.**
  Bin counts are 1 500 / 215 / 23 / 12 against the 187 / 23 / 12 that exact
  volume conservation would give — spot on at f = 4 and f = 5, 15% high at f = 2.
  Defining Ω_b once at the finest mesh and mapping it down would make this exact.
- **The dose floor still depends on the estimator**, through
  `nanmax(dose0)`. Bin counts at f = 1 vary 1 492–1 516 across draws, a residual
  1.6% against the ~33% the `rel_err` cut was causing.
- **The floor is measured unpaired.** For the large-effect scenarios, pairing
  retains little, so the unpaired floor is the right comparison; for the small
  levers it is conservative by design.

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

Two open items, in order.

**The magnitude is not converged, and f = 1 is the finest tally available** — the
mesh cannot be refined below the CPM lattice pitch. Converging it requires a
finer *snapshot*, which is a Reference D question about the selected pitch rather
than a transport-numerics one. Note that a naive snapshot refinement changes what
a parcel physically *is*: `V_target = 120` in `biofilms_potts.jl` is in **sites**,
not physical volume, so raising N at fixed `V_target` shrinks every parcel.

**The floor must be driven below the effects of interest** before any lever
below ~0.05 is quoted as resolved. That is ~40× the pilot's histories, and it is
now a projection with a measured scaling law behind it rather than an assumption.
