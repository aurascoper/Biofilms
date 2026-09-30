# Draft registration: does the lost melanin crossing return when volume stiffness scales with production?

**Status: draft for review, 2026-09-30. Not registered.** Registration takes effect at the first clean commit that removes this banner and adds the producer and its tests together. Until then no production run may start, and nothing here is a result.

The experiment below uses a different model, not a control arm of the sweep. It changes the volume stiffness λ_V (the weight of the volume penalty) for every species. It asks whether the melanin crossing depends on the ratio of the melanin term to the volume term, rather than on melanin alone. That is ledger question MEL-COLLAPSE-01, and this registration replaces the larger draft in PR #67 if both are approved.

## 1. What is already known

These facts were read before this draft, so no readout below may treat them as blind predictions.

- **Sweep result (66ffe9f, 256 seeds, MCS 100).** MCS is a Monte Carlo step, the model's time unit. The paired full-minus-null mean radius for *C. neoformans* is −0.0899 at a production multiplier of 10⁵ and −0.0031 at 10⁶. The inward criterion is −0.05, in units of the cylinder radius.
- **Volume.** *C. neoformans* parcels average 1,466 sites at 10⁵ and 3,665 at 10⁶, against a target of 120. The two radiotropic species fill 35.1% of the interior at 10⁵ and 87.6% at 10⁶. Radiotropic here means *C. neoformans* or *C. sphaerospermum*.
- **Within-arm check.** At 10⁶, the slope of *C. neoformans* radius on joint occupancy across seeds is +0.079, with a 95% interval of −0.214 to +0.374. So occupancy alone does not predict radius within that arm.
- **Energy terms.** A radiotropic parcel that takes a site gains −0.5·M, where M is the melanin at that site. Adding one site to a parcel of volume V costs λ_V·(2(V − 120) + 1), with λ_V = 10 by default.
- **Linearity.** For a fixed occupancy history, the melanin field is proportional to the production multiplier. This is derived from the update rule and is not yet tested.

## 2. Why a volume cap is not used

A hard volume cap would freeze the radiotropic parcels at 10⁶. A parcel at its cap cannot gain a site. Any move that takes a site from it costs about 0.5·M, which is near 900,000 at MCS 100, so such moves are almost never accepted. The crossing would then disappear from trapping, which is a different mechanism.

## 3. The intervention

At production multiplier 10⁶, set λ_V = 100 instead of 10. For a fixed occupancy history, the melanin and volume terms then both have ten times their 10⁵ values. Adhesion and the direct radiation term keep their default values.

The scaled model is not an exact copy of the 10⁵ model. Adhesion becomes smaller relative to the scaled terms, and every species, not only the radiotropic ones, has a stiffer volume. The registration tests the ratio, not an equivalence.

## 4. Fixed design

- **Seeds:** every integer 42 to 297, paired within seed.
- **Configuration:** N = 40, six parcels per species, MCS 100, the sweep's snapshot schedule.
- **Arms:** full, with production ×10⁶ and λ_V = 100; and null, with both radiation vectors zero and λ_V = 100. That makes 512 runs.
- **Producer:** a new file in this folder that reuses the sweep's parameter builder and adds one keyword for λ_V. The sweep producer is not edited, because the archived shards pin its checksum.
- **Identity gate:** before production, the new producer with λ_V = 10 must reproduce the archived 10⁶ full row and null row exactly for every seed. A seed that fails blocks the whole run.
- **Output:** one atomic CSV per seed and arm, outside git, with the sweep's columns and provenance header, plus λ_V in the header.

## 5. Readouts and prediction

- **Primary:** the paired full-minus-null mean radius for *C. neoformans* at MCS 100, with a normal 95% interval over seeds. The criterion is −0.05, the same value the sweep used.
- **Prediction:** the crossing returns. The paired mean is −0.05 or lower.
- **Descriptive:** joint radiotropic occupancy, volume per parcel, and the same contrast for *C. sphaerospermum*.

Possible outcomes, stated before any run:

- **Crossing returns and occupancy falls near the 10⁵ value.** A result like that supports the ratio of the two terms as the reason the crossing was lost.
- **Crossing does not return although occupancy falls.** Space filling is then not enough to explain the loss.
- **Occupancy stays high.** λ_V = 100 did not restrain expansion, and the result is inconclusive.

A null result is retained. No threshold, arm or seed set changes after inspection.

## 6. Controls

- **Scaling applied.** Each shard header must state λ_V = 100 and production ×10⁶. A shard produced with the default λ_V must fail validation.
- **Identity gate.** A mutated archived row must make the reproduction check fail.
- **Smoke test.** N ≤ 16, one parcel per species, at most four MCS, and seeds outside 42 to 297.

## 7. What even a clear answer would not show

A returned crossing would show that this model's crossing depends on the balance of two of its terms. It would not show that space filling is the only route, and it says nothing about real melanised fungi.

## 8. Cost

The sweep averaged about 4.4 s per run as a whole process. At that rate, 512 runs take about 38 minutes on one worker, or about 10 minutes on four.
