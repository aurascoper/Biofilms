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
| `clearly_detectable` | `PASS_SYNTHETIC_GATE` | ✅ same | 0.2441 |

72 OpenMC runs, 7.2 × 10⁶ histories, 39 s wall. The zero-effect case returned
**exactly** zero, which is the false-positive control that licenses every other
row.

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

## Variance: numerics is 97.7% of it

| Source | Variance | Share | What reduces it |
|---|---|---|---|
| **numerics** | 5.92 × 10⁻⁴ | **97.7%** | finer mesh, more raytrace samples. **A floor.** |
| calibration | 9.35 × 10⁻⁶ | 1.5% | better measurements only |
| transport | 4.28 × 10⁻⁶ | 0.7% | more histories and replicates |
| model form | 0 | 0% | a decision |

This is the A0 lesson reproduced on the biofilm geometry and quantified: **the
resolution floor is ~140× the Monte Carlo noise.** Buying more histories would
buy almost nothing, and a single pooled "uncertainty" number would have hidden
that completely — which is the reason the gate reports the dominant source in
its verdict rather than a total.

## Measured cost, and the projection

| Quantity | Measured |
|---|---|
| seconds per run | 0.542 |
| histories per second | 184 612 |
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

## Next

`near_threshold` (20% Gd, effect 0.1241) is the honest candidate for the
uncertainty-dominated case: it sits just above the threshold, so its verdict is
a finding rather than a declaration. It should be run once the mesh-resolution
question is addressed, since numerics currently dominates the spread it would be
judged against.
