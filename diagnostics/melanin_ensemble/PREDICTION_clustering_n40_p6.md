# Registered prediction: melanin and contact between distinct parcels

Registration date: 2026-09-29. This plan takes effect at the first clean commit
containing this file, its producer, analyser and tests. Production runs must
record that commit and this file's SHA-256. Only synthetic controls and small
off-plan smoke runs may precede that commit.

## What was already known

The committed radial ensemble already reports seeds 42–297, N = 40, six parcels
per species and MCS 100 under four arms. Its full-minus-no_mel mean-radius
contrasts for C. neoformans (CN) and C. sphaerospermum (CS) are small. These
radius data informed this hypothesis and are not a new independent confirmation.
No contact endpoint for these runs has been inspected in planning this study.
The melanin path has no positive control. The model favors melanin-rich sites
for CN and CS, but its radial result alone establishes neither retention nor
clustering. This study tests contact between distinct parcels; it does not
measure temporal retention of a parcel's own sites or an effective adhesion
coefficient.

## Fixed experiment and observable

Run the existing serial model at N = 40, six parcels per species, to MCS 100,
with all 256 seeds 42–297. Compare full with no_mel, paired within seed.
The full arm uses shipped CPMParams. The no_mel arm sets every
alpha_M_species entry to zero; beta_ion and all other model parameters stay
fixed. The model's original run_simulation initializes and consumes its own
MersenneTwister stream, exactly as radial_sweep.jl does. No alternative stepper
or simulation loop is introduced. Read MCS 0 and 100.

For each occupied parcel i, count all six face directions from each of its
sites. Ignore neighbors with the same parcel ID. Each remaining face is one
boundary face B_i, including faces toward medium, a wall, outside the array,
or any distinct parcel. Among these faces, H_i counts only neighbors belonging
to a DISTINCT parcel of the SAME species. Define h_i = H_i / B_i. The species
endpoint C_s is the arithmetic mean of h_i across its occupied parcels, giving
each parcel equal weight. An interface between two same-species parcels enters
each parcel's own count once. Internal faces of one parcel never count as
contact between parcels. A species with no parcel has an undefined endpoint;
it is never assigned zero and never silently omitted.

The denominator measures the parcel's boundary, rather than its volume or the
entire species' pooled surface. Thus a larger parcel does not receive a larger
weight. The ratio is unchanged under uniform voxel refinement of the same
geometry; it is not invariant to a real change in parcel shape or exposure.
The six face directions measure lattice surface contact. The acceptance rule
continues to use its existing 26 neighbors; this endpoint is not its energy.
The CSV retains each parcel's ID, volume, H_i and B_i so C_s can be recomputed.

The companion radius endpoint is the mean of parcel centre-of-mass radii,
divided by R = N/2, as in the original radial ensemble. Record the original
snapshot statistics as well as contact counts. Every original radius,
volume, parcel count and mean-melanin field must reproduce the committed
radial CSV at its six-decimal precision, for both arms and both times.

## Primary prediction and decision rule

The compound prediction is: under full minus no_mel, CN AND CS each have
increased same-species contact while their mean-radius contrasts meet the
existing equivalence bound. For EACH species require ALL of:

- Mean paired contact difference at MCS 100 is at least 0.05 (five percentage
  points). This is an uncalibrated design choice, selected before seeing any
  contact results, not an effect derived from this model's output.
- Exact one-sided paired sign test for positive contact difference gives
  p < 0.05. Exact ties carry no sign and leave both numerator and denominator;
  all ties give p = 1. The upper binomial tail uses p = 1/2.
- The paired radius mean's normal 90% interval, mean plus or minus
  1.6448536269514722 times paired SD / sqrt(256), lies strictly inside
  [-0.05, +0.05] r/R. This matches the earlier radial analysis's interval
  convention and existing equivalence bound. It permits smaller nonzero effects.

The overall prediction is supported only if all conditions hold for both CN
and CS. This is a conjunction, not a search for any successful species. No
separate family-wide claims for individual species will be made. Otherwise
report that the compound prediction did not meet its registered criterion,
and identify each condition that failed; do not call that proof of no contact
effect. Report paired contact mean, SD, sign counts/ties and exact p; radius
mean, SD and interval. All other species and MCS 0 contacts are descriptive.
No post-hoc threshold, endpoint, horizon or species substitution is allowed.

The analyser's constants are bound to these declared values by the tests:

| Constant | Value |
|---|---|
| CONTACT_ALPHA | 0.05 |
| CONTACT_MEOI | 0.05 |
| CONTACT_RADIUS_BOUND | 0.05 |
| CONTACT_Z90 | 1.6448536269514722 |

## Completeness, controls and limits

Each seed must contain full and no_mel, both times and all seven species.
Duplicate or missing rows, incompatible metadata, nonfinite endpoints,
orphan parcel IDs, invalid counts, unequal MCS 0 lattice hashes, and baseline
disagreement stop analysis. Partial runs produce no primary verdict. Shards
may be validated for operational restart without issuing a hypothesis result;
the primary analysis requires all 256 seeds. There is no interim stopping rule.

Controls include single-parcel shapes, face-touching distinct parcels,
diagonal-only proximity, different species, unequal parcel sizes, and uniform
voxel refinement. Known-bad inputs test orphan detection, pairing and radius
reference checks. Smoke tests use N <= 16, one parcel per species, at most four
MCS, and seeds outside 42–297; their results cannot support this prediction.

Even a successful result would show increased inter-parcel contact in this
model under these declared coefficients. It would not establish adhesion in
organisms, temporal retention, a biological mechanism, or radiotrophy.
