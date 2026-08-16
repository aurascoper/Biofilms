# OpenMC nested pilot — variance, cost, and which lever actually moves dose

**Date:** 2026-08-16 · **Tier:** S0 · **`target_calibration = false`** ·
**Policy:** `synthetic_gate_contract_v1` · **Reference D verdict:** `NOT_EVALUATED`
**Driver:** `coupling/scripts/openmc_nested_pilot.py` ·
**Artifacts:** `artifacts/pilot/openmc_nested_pilot_{budget,verdict}.json`,
`data/calibration/openmc_effect_samples.csv`

The pilot's job was to buy a budget before spending one, and to check that the
gate returns the verdicts declared for it in advance. It did both, and it
overturned the assumption the scenario registry was originally built on.

## The gate returned exactly what was declared

| Scenario | Declared verdict | Returned | Effect (rel. L2) |
|---|---|---|---|
| `zero_effect` | `EFFECT_BELOW_THRESHOLD` | ✅ same | 0.0000 |
| `clearly_detectable` | `PASS_SYNTHETIC_GATE` | ✅ same | 0.2442 |
| `near_threshold` | `DECIDED_BY_PILOT` | `PASS_SYNTHETIC_GATE` | 0.1384 |

108 OpenMC runs, 1.08 × 10⁷ histories, 60 s wall. The zero-effect case returned
**exactly** zero, which is the false-positive control that licenses every other
row. `near_threshold` had no declared verdict by design — see below.

## The finding: density is a weak lever, and it is weak for a structural reason

The registry originally defined "clearly detectable" as a **35% biomass density
increase**. Measured, that moves the mass-weighted relative L2 dose effect by
**1.4%** — nowhere near the declared 0.10 threshold.

That is not a bug. **Dose is energy per unit mass**, so raising the density
raises the numerator and the denominator together:

| Quantity | Baseline | ×1.35 density | Ratio |
|---|---|---|---|
| heating in biomass (eV/src) | 170 682 | 222 876 | **×1.306** |
| mass in biomass (kg) | 2.95922 | 3.99477 | **×1.350** |
| dose = heating/mass | 57 678 | 55 792 | **×0.967** |

The two nearly cancel and only second-order self-shielding survives. **A uniform
density change is close to invisible in an absorbed-dose field.**

### Composition is not automatically stronger either

At Cs-137's 661.7 keV, Compton scattering dominates, and the Compton mass
attenuation coefficient tracks **Z/A**, which is ≈ 0.5 for nearly every element.
So swapping oxygen for iron changes almost nothing. What moves the dose is a
**departure from Z/A ≈ 0.5** — hydrogen at 1.0, or high-Z material where the
photoelectric term still contributes.

| Material lever | Relative L2 effect |
|---|---|
| 5% Fe (replacing O) | **0.0033** ← the weakest tested |
| density ×1.35 | 0.0137 |
| 5% Gd | 0.0376 |
| dehydration, H 0.112 → 0.06 | **0.0459** |
| 20% Gd | 0.1241 |
| 40% Gd | 0.2242 |

**Two consequences worth carrying forward.**

The strongest *physically realistic* lever here is the **hydrogen fraction** —
which is exactly what the paired wet/dry gravimetry in the Reference D protocol
measures. That is an independent argument for the campaign's priority ordering:
the water fraction is not merely a bookkeeping input, it is the material
property the transport is most sensitive to.

And the metal-loading story is weaker than the project's framing implies. 5% Fe
is the *least* effective perturbation tested. Metal loading moves absorbed dose
only when the metal is high-Z; `S. oneidensis` reducing iron does not, by
itself, make a transport-visible material change at these energies.

`clearly_detectable` was therefore redefined as an extreme synthetic Gd 40%
material with no biological pretension — its only job is to prove the gate can
pass. **The 0.10 threshold was not adjusted to reach it**, which would have been
a stop condition.

## Variance: numerics dominates, by two orders of magnitude

Per scenario, because a global maximum charges every case the largest case's
discretisation error:

| Scenario | transport | numerics | calibration | dominant |
|---|---|---|---|---|
| `clearly_detectable` | 8.5 × 10⁻⁶ | **6.0 × 10⁻⁴** | 9.4 × 10⁻⁶ | numerics |
| `near_threshold` | 3.3 × 10⁻⁶ | **3.8 × 10⁻⁴** | 1.6 × 10⁻⁶ | numerics |
| `zero_effect` | 1.6 × 10⁻³² | 6.6 × 10⁻³⁴ | 3.2 × 10⁻³³ | — (no effect) |

This is the A0 lesson reproduced on the biofilm geometry and quantified: **the
resolution floor is 70–115× the Monte Carlo noise.** Buying more histories would
buy almost nothing, and a single pooled "uncertainty" number would have hidden
that completely — which is the reason the gate reports the dominant source in
its verdict rather than a total.

Note `zero_effect`'s budget: ~10⁻³³ across the board, because an identical
material state produces an identical field at every resolution. Under the old
global-maximum budget it was charged 6.0 × 10⁻⁴ of somebody else's numerics.

## Measured cost, and the projection

| Quantity | Measured |
|---|---|
| seconds per run | 0.55 |
| histories per second | 181 230 |
| statepoint per run | 94 440 B |
| raytrace (5 × 10⁵ rays) | 0.22 s, once per geometry |

**Hardware provenance:** AMD Ryzen AI 9 HX 370, OpenMC 0.15.3 CPU-only
(MPI/OpenMP), ENDF/B-VIII.0. Wall times are **machine-specific and not
portable**; the history counts are. No GPU throughput is reported because the
pinned OpenMC has no GPU path and a borrowed number would describe a different
program.

Projected full Reference D study — 128 outer draws × 2 material states × 20
seeds × 2 mesh factors = **10 240 runs**:

| Histories per run | Serial wall | On 32 cores | Storage |
|---|---|---|---|
| 4 × 10⁶ (A0's setting) | ~62 h | ~2 h | ~1 GB |
| 1 × 10⁶ | ~15 h | ~30 min | ~1 GB |

**The recommendation is the lower row.** With numerics at 97.7% of the variance,
quadrupling the histories buys a 0.7% component and leaves the 97.7% one
untouched. The budget belongs in mesh resolution and raytrace samples.

This is a laptop-overnight or a few-core job. It does not need a cluster, and
sizing it against one would be answering a question the measurement does not
ask.

## What this does not establish

No target parameter was inferred. Reference D's verdict is `NOT_EVALUATED`, the
offline feedback gate is untouched, online feedback stays disabled, and every
artifact carries `target_calibration = false`. The Gd 40% material is invented
and is not a claim about any organism.

The sensitivity ranking above is measured **on the synthetic geometry at
661.7 keV** and is not transferable to another source spectrum: the Z/A argument
that makes Fe uninteresting here would not survive a lower-energy source, where
photoelectric absorption scales steeply with Z.

## `near_threshold`, run 2026-08-16

The case whose verdict was declared `DECIDED_BY_PILOT` came back
**`PASS_SYNTHETIC_GATE`**, and the way it passed is the interesting part.

| | effect | 1% lower bound | verdict |
|---|---|---|---|
| `near_threshold` (20% Gd) | median 0.1384 | **0.1202** | `PASS_SYNTHETIC_GATE` |

Its variance is numerics-dominated — 3.8 × 10⁻⁴ against 3.3 × 10⁻⁶ transport,
a factor of 115 — and the effect still differs by **22% between mesh factors**
(0.1236 at f = 1, 0.1513 at f = 2). Yet the verdict is a pass, because the
distribution does not *straddle* the threshold: even the coarse end sits above
0.10.

**A dominant variance source does not by itself make a verdict indeterminate.**
Indeterminacy requires the spread to reach across the decision boundary, and
here it does not. The gate distinguished "noisy" from "noisy enough to matter",
which is the distinction it exists for.

The honest summary is therefore split: **the decision is resolution-robust; the
magnitude is not converged.** Anyone quoting 0.1384 as *the* effect would be
over-reading it by 22%.

## Two defects this scenario exposed

Running a third scenario made both visible, and both had been silently wrong.

**The variance budget took a global maximum across scenarios.** Every scenario
was judged against the largest scenario's discretisation error, so
`near_threshold` inherited `clearly_detectable`'s spread and `zero_effect` was
credited with a numerics variance it does not have. Budgets are now computed
per scenario: `near_threshold`'s own numerics variance is 3.8 × 10⁻⁴, not the
6.0 × 10⁻⁴ it was being charged, and `zero_effect`'s is ~10⁻³³.

**`omega_b` subsampled occupancy under coarsening.** Taking `occ[::factor]`
picks one representative voxel per coarse bin, so a bin that is 12% biomass
counted as fully biological and the two mesh factors measured different regions
— discretisation error that was really a definitional artifact. Occupancy now
comes from the exact CSG volumes already computed for the mass denominator: a
bin is biological when at least half its volume is biomass. Factor-2 bins went
from 199 to 215, and the per-factor spread survived, which is how we know the
22% is genuine discretisation rather than a masking artifact.

## Next

The mesh disagreement is the open item, and it has a floor: **factor 1 is the
finest mesh available**, since the tally cannot be refined below the CPM lattice
pitch. Converging the magnitude therefore requires a finer *snapshot*, not a
finer tally — which is a Reference D question about the selected pitch, not a
transport-numerics one.
