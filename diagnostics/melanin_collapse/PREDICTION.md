# Draft protocol: why the melanin crossing is lost at a million-fold production

**Status: draft for review, 2026-09-30. Not registered.** Registration takes effect at the first clean commit that removes this banner and adds the producer, the analyser and their tests together. Until then no production run may start, and nothing here is a result.

The question is ledger row MEL-COLLAPSE-01. In the melanin-production sweep, the *C. neoformans* inward crossing was present at multipliers of ten thousand and one hundred thousand, then disappeared at one million. This protocol records what the model does at each multiplier, so the three candidate explanations can be told apart. It changes no model term and makes no biological claim.

## 1. What is already known

These facts were read before this draft, so no readout below may treat them as blind predictions.

- **Sweep result (66ffe9f, 256 seeds, MCS 100).** MCS is a Monte Carlo step, the model's time unit. The paired full-minus-null mean radius for *C. neoformans* is −0.0313 at 10³, −0.0815 at 10⁴, −0.0899 at 10⁵ and −0.0031 at 10⁶, in units of the cylinder radius.
- **Melanin keeps rising.** Mean melanin on *C. neoformans* sites is about 16,200, 176,000 and 1,774,000 at 10⁴, 10⁵ and 10⁶. The update has no upper cap.
- **Volume takeover.** Mean volume at MCS 100 is about 700 sites per species in the null arm. For *C. neoformans* it is 2,015, 8,798 and 21,992 at the three highest multipliers. For *C. sphaerospermum* it is 2,114, 8,831 and 21,992. The interior has 50,200 sites.
- **The melanin term reads one site.** In `compute_delta_H_terms`, a radiotropic source gains −0.5·M at the target site, and a radiotropic target loses +0.5·M. The term uses no gradient. Radiotropic here means *C. neoformans* or *C. sphaerospermum*.
- **Acceptance.** A proposal with ΔH ≤ 0 is accepted without a draw. Otherwise it is accepted with probability exp(−ΔH/T), where T = 5.
- **Scale.** Half the mean melanin, divided by T, is about 1,600 at 10⁴ and 177,000 at 10⁶. So at 10⁴ the melanin term is already far larger than the other terms for each radiotropic proposal, and the crossing is still present there.
- **Linearity (derived, not yet tested).** For a fixed occupancy history, the melanin field is linear in the production vector. Scaling production by k scales the field by k and leaves its normalised shape unchanged. Any change of shape between multipliers must therefore come through a different occupancy history.

## 2. Candidate explanations

1. **Field flattening:** the normalised radial melanin profile loses its slope at 10⁶, so no site is preferred by radius.
2. **Acceptance saturation:** favoured radiotropic moves are accepted with certainty at every radius, so acceptance stops depending on radius.
3. **Space filling:** the two radiotropic species occupy most of the interior, so few moves remain that the melanin term can reward.

These can occur together. The volume facts in section 1 were seen before this draft, so the third explanation gets descriptive readouts only.

## 3. Fixed design

- **Seeds:** every integer 42 to 297, the same 256 seeds as the sweep.
- **Configuration:** N = 40, six parcels per species, MCS 100, full arm, production multipliers 1, 10³, 10⁴, 10⁵ and 10⁶, plus the null arm. That makes 1,536 runs.
- **Only archived levels.** Every level already has an archived shard, so every run can be checked against it.
- **Stepper.** `run_simulation` has no hooks, so the producer drives the same loop by hand. It seeds one `MersenneTwister(seed)`, calls `init_state`, and then runs each MCS in the archived order. That order is `mcs_step!`, `update_melanin!`, `update_nutrient!` and `update_centers_of_mass!`.
- **Observation hook.** A wrapper passed as `energy_terms` calls `compute_delta_H_terms` once and returns its result unchanged. The wrapper records site coordinates and computes the acceptance probability by arithmetic. It draws no random number.
- **Identity gate.** For every run, the ten printed columns at MCS 0 and MCS 100 must equal the archived 66ffe9f shard exactly. A run that fails is void and gets rerun, never analysed. This gate is how a new revision proves that the default path is unchanged.
- **Storage.** No per-proposal record is kept. Counters are accumulated for each MCS and written as one atomic file per run, outside git.

## 4. Readouts

Each readout is recorded for every MCS from 1 to 100. The interior is split into three equal-area radial bands: inner, middle and outer.

- **R1, field shape.** The mean melanin in each band, divided by the interior mean. The slope is the inner value minus the outer value.
- **R2, radial selectivity.** Take every proposal with a *C. neoformans* source and a non-radiotropic or empty target. Selectivity is the mean acceptance probability for inner-band targets minus that for outer-band targets.
- **R3, saturation share.** Among the proposals in R2, the share with ΔH ≤ 0. Also, for proposals in which a non-radiotropic source takes a *C. neoformans* site, the share with acceptance probability below 10⁻⁶.
- **R4, time to saturation.** The first MCS at which the saturation share in R3 is 0.99 or more.
- **R5, occupancy (descriptive only).** Volume of each radiotropic species, their joint share of the interior, and the count of *C. neoformans* faces that touch a non-radiotropic or empty site.

## 5. Declared readings

Each reading uses the seed mean at MCS 100 with a 95% percentile bootstrap: 5,000 samples, analysis seed 9005. The cut-offs are declared choices, not derived values. None was chosen from R1 to R4, which no one has computed.

- **Flattening is supported** if the R1 slope at 10⁶ is less than half the slope at 10⁴, and the upper bound of the ratio's interval is below 1.
- **Saturation is supported** if the R2 selectivity at 10⁶, averaged over MCS 1 to 100, is less than a quarter of its value at 10⁴. R4 must also be earlier at 10⁶ than at 10⁴.
- **Space filling gets no verdict.** R5 is reported beside the other readouts, but its values were partly seen before this draft.

Each reading is evaluated separately, and more than one can be supported. A null result on both is retained.

## 6. Controls

- **Inert hook.** A run with the wrapper must equal a run without it: lattice, fields and the next 20 draws of the random stream. This follows `tests/rad_proposals_tests.jl`.
- **Hook that draws.** A wrapper that consumes one random number must change the trajectory. If it does not, the inert-hook test cannot fail and is void.
- **Linearity.** With the lattice fixed, scaling production by k must scale the melanin field by k within floating-point rounding. A second case, with occupancy changed between the two fields, must fail the same check.
- **Selectivity sign.** A hand-built lattice with more melanin on the inner band must give a positive R2 at a small multiplier.

## 7. What even a clear answer would not show

These readouts report model quantities that change together. They do not show which one causes the lost crossing. A causal claim needs a separately registered intervention, such as a cap on the melanin term's magnitude or a fixed radiotropic volume. Nothing here applies to real melanised fungi.

## 8. Cost

The unhooked sweep took 3 h 9 min for 4,352 runs on two workers, about 4.4 s per process. The hook adds work to each proposal, so the time per run must be measured in a smoke test before production. The smoke test uses N ≤ 16, one parcel per species, at most four MCS and seeds outside 42 to 297.
