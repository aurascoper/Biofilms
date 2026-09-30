RECORD: line numbers refer to the preprint tex at d7b7c28, the commit that added this file.

The same commit holds the preprint tex revision it cites. At `56b2380` the tex differed by 213 added and 200 removed lines and eight of the cited passages did not exist. Each citation also names the enclosing label and a phrase `grep -F` finds, so it can be re-resolved if the tex changes again.

# An outside generator cited this preprint for four real passages, distorted five and invented one

2026-09-30 · not a calibration, not a claims-ledger entry, not a dataset candidate.

## What this is

An AI writing tool was given the preprint (`preprint/modeling_radioresistance_and_radiotropic_fitness.tex`, "Does a Dose Gradient Move Biomass Parcels?", version 2.0) together with an overview of Israeli space technology. It produced three things: statements about what "the paper" says, a rewrite of the overview citing the paper, and two Python "telemetry engines" plus two equations said to follow the paper's "strict modeling constraints". A second message from the same tool then retracted most of it and proposed an audit file for this directory.

This document checks each attribution against the tex and the executed model (`diagnostics/fitness_phase/`), and checks the retraction too. It populates no calibration value, creates no ledger row and changes no code. It is kept so the attributions do not resurface as if someone had checked them.

The revision the tool read is the one committed with this record, not `56b2380`. Four of its quotations (the seven assumed sensitivities, the phase parameters $2\pi/20$ and $K=0.5$, the word "stencil") do not exist at `56b2380`, where the survival and phase subsections are still marked "specification, not a reported method".

Verdict vocabulary is the closed set from `data/claims_ledger.csv` (`keep | restate | requalify | delete | needs_calibration | needs_verification`, plus `supported` and `must_not_be_claimed`), used as an assessment label only. No row is added to that file.

---

## Claim 1 — "A nine-order-of-magnitude gap between laboratory exposures (Gy min⁻¹) and chronic deep-space settings (mGy yr⁻¹)"

**As stated (tool):** parameters or sensor calibrations from short, intense exposures are "non-transferable to deep-space chronicity".

**What the tex says** (`sec:planned_feedback`, lines 1689–1697, phrase "nearly nine orders of magnitude apart"; also `sec:discussion`, lines 959–962): a reactor irradiation delivers of order Gy min⁻¹, the environmental settings the paper cites as motivation are closer to mGy yr⁻¹, and "a sensitivity constant calibrated in one regime is not transferable to the other".

**Verdict:** `restate`. The number and the non-transfer statement are the paper's. The two poles are reactor irradiation and contaminated field sites. Deep space and sensors are the tool's.

---

## Claim 2 — "Damage is indirect, via radiolysis products, in a race between generation and the body's cellular repair and recombination"

**What the tex says** (`sec:planned_feedback`, lines 1690–1692, phrase "competing generation and recombination kinetics"): in aqueous systems the dominant damage pathway is indirect, through radiolysis products whose steady-state concentration is set by competing generation and recombination kinetics, both rate-dependent.

**Verdict:** `restate`. Recombination here is radical chemistry. "Repair" occurs eight times in the tex, all in the species-evidence appendix, never as a model term. "Cellular repair" is the tool's.

---

## Claim 3 — "Compton scattering depends on electron density, so light-element matter is compositionally transparent; use this to test whether regolith or ice changes are transport-visible"

**What the tex says** (`sec:planned_feedback`, lines 1645–1656, phrase "dominated by Compton scattering"): at the photon energies of interest, attenuation through C, H, N and O is Compton-dominated, scales with electron density per unit volume and is insensitive to atomic number; the one-way simplification holds "only while the sorbate is dilute and low-Z".

**Verdict:** `requalify`. A photon-transport statement. Solar-particle-event protons and galactic-cosmic-ray nuclei are charged particles, as the tool's own table says, and nothing in the paper covers their transport.

---

## Claim 4 — "The paper's offline counterfactual: if material changes are transport-invisible, an expensive coupling is unjustified"

**What the tex says** (`sec:discussion`, lines 948–951 and 974–980, phrases "counterfactual question is logically prior" and "not as a proposed deployment"): the same sentence, about a biofilm in a spent-fuel pool perturbing its own local field, recorded "as an observation about regime and available inputs, not as a proposed deployment".

**Verdict:** `restate`. Sensors, CubeSats and "multi-million-dollar sensor adjustments" are the tool's.

---

## Claim 5 — "The paper adds highly specific physics: lead shatters into secondary radiation, hydrogen-rich polyethylene, a tissue-shielding paradigm, launch-weight constraints, stem-cell niches"

**What the tex says:** nothing. Zero lines match lead, polyethylene, HDPE, spallation, secondary particle, launch or stem cell. "Shield" occurs four times, all about radiotrophic-fungus experiments. "Niche" occurs once, as "radiotrophic niche" (`sec:limitations`, line 1151).

**Verdict:** `must_not_be_claimed` as a statement about this paper. The physics is textbook and can be cited elsewhere.

---

## Claim 6 — "Fitness is modelled as a discrete survival/damage hazard, not a continuous field"

**What the tex says:** `sec:psde`, lines 1281–1300, phrase "specification, not method": no program integrates the $F_s$ field equation. `app:synthetic_survival`, `eq:survival_hazard`, line 1329: what runs is $p_i=-\operatorname{expm1}(-\alpha_s\,\Delta D_i)$ per parcel. Line 1342, phrase "nonviable immediately": a committed parcel is nonviable at once.

**Verdict:** `supported` as a reading. The tool's continuous 0–1 index, `fitness *= exp(-hazard*dt)`, is its own construct. The paper's parcels are binary.

---

## Claim 7 — "Eq. 15 is the survival/damage hazard"

**What the tex says:** the survival hazard is `eq:survival_hazard`, the seventh numbered display in both revisions. The fifteenth is `eq:fene` (line 1591) in the committed revision and `eq:kelvinvoigt` at `56b2380`. Neither is a hazard. Both mechanics equations are unexecuted: `sec:specified`, line 449, phrase "FENE and Kelvin--Voigt mechanics, for which", says no stress, strain or modulus exists in any source, and `biofilms_potts.jl` has zero hits for FENE, Kelvin or Voigt.

**Verdict:** `delete`.

---

## Claim 8 — "α_s damage coefficient = 0.015 Gy⁻¹" (second engine: 1.2e-4)

**What the tex says** (`app:synthetic_survival`, line 1335, phrase "(0.01,0.002,0.015,0.04,0.02,0.06,0.03)"): seven assumed sensitivities in `tab:params` species order, and "These choices exercise a response law; they are not fitted species values or a calibration of mechanical time".

**Verdict:** `restate`. 0.015 is the third species' assumed value. 1.2e-4 has no source. The executed values are read from `diagnostics/fitness_phase/config.toml` line 12 (`alpha_Gy_inv`).

---

## Claim 9 — "Kuramoto intrinsic frequencies ~ N(2π/20, 0.05) rad/s, K = 0.5"

**What the tex says** (`app:phase_passenger`, lines 1375–1377, phrase "frequencies with mean"): mean $2\pi/20$ and SD 0.05 rad per MCS, scalar $K=0.5$ per MCS, $\Delta t=0.05$ MCS, "All parameters are synthetic". The clock is an assumed 60 s per MCS (line 1336, phrase "The assumed clock is 60"). The code reads the values from `diagnostics/fitness_phase/config.toml` lines 18–20 (`phase_omega_mean`, `phase_omega_sd`, `phase_K`).

**Verdict:** `delete` the unit conversion. Rad per MCS became rad per second, so a 20-minute period became 20 seconds. The values were never physical.

---

## Claim 10 — "Radiation is a forcing on phase, scaled by (1 − F); the γ gain follows radical load; fitness gates the disruption"

**What the tex says** (`app:phase_passenger`, lines 1377–1380, phrases "radiation forcing are off" and "no temporal phase for entrainment"): "Noise and radiation forcing are off. Phase never feeds back into survival, dose, growth, movement or acceptance." "A static dose gradient supplies no temporal phase for entrainment." Line 247, present in both revisions: "a read-only phase diagnostic with radiation coupling disabled". The terms $\gamma_s\Delta_s$, $\alpha_{s,\mathrm{nir}}N$ and $\beta_{s,\mathrm{ion}}I_\gamma$ belong to the unexecuted field equation `eq:psde`.

**Verdict:** `must_not_be_claimed`. The tool's central coupling is the one the paper says it did not run.

---

## Claim 11 — "The order parameter R = |Σ eⁱᶿ| / N measures coherence"

**What the code and tex say:** diagnostics/fitness_phase/phase.jl@`coherence(ids,theta) = isempty(ids) ? NaN : abs(sum(cis(theta[id])` computes exactly that. `app:phase_passenger`, lines 1384–1392, phrase "Coherence is reported with active count": coherence with active count, edges, isolates and components, and a locking verdict that requires endpoint drift at most 0.01 rad/MCS, relative-phase range at most 0.5 rad, no full winding, and a fall of at least 50 % in effective-frequency SD against a zero-coupling replay.

**Verdict:** `supported` for the formula. The paper's verdict is stricter than "R drops below baseline", and it never uses the words "order parameter".

---

## Claim 12 — "Spatial network stencil mapping", "edge-couplings and spatial coherence metrics"

**What the tex says** (`app:phase_passenger`, lines 1369–1370, phrase "26-neighbour adhesion stencil"): the Boolean adjacency records contacts under the 26-neighbour adhesion stencil; edges, isolates and components are reported.

**Verdict:** `restate`. Traceable, garbled.

---

## Claim 13 — "Circadian, cardiac, HRV, RMSSD and sleep-state modelling adapted from the paper"

**What the tex says:** one clause in the introduction's Kuramoto literature sentence, line 260, "from neural networks to circadian rhythms". No other line matches circadian, cardiac, heart or sleep.

**Verdict:** `must_not_be_claimed`.

---

## Claim 14 — "Deep-space background 3.17e-11 Gy/s (about 1 mGy/yr); vest factor 0.418 applied to it"

**Not a paper claim.** Deep-space galactic-cosmic-ray absorbed dose is about 0.5 mGy per day (MSL RAD cruise measurement, Zeitlin et al. 2013), about 180 mGy per year. 0.418 is the solar-particle-event-scaled AstroRad figure, and the tool's own table separates SPE from GCR.

**Verdict:** `delete`. Wrong by about two orders of magnitude.

---

## Claim 15 — AstroRad: 58.2 %, 193 days, August 1972 event, phantoms Helga and Zohar, Science Advances 10.1126/sciadv.adz1892

**Not a paper claim.** The tool cited Science Advances and ESA for these figures, not the preprint. ESA's Orion blog of 28 August 2026 (fetched 2026-09-30) gives all four. The DOI returned HTTP 403 on fetch.

**Verdict:** `supported` by ESA; `needs_verification` for the DOI.

---

## The retraction restates the findings, adds one correct claim and three wrong ones

The tool's second message withdrew claims 1, 3, 5, 7, 8, 9, 10, 13, 14 and 15 and proposed an audit file titled "Post-Payload Verification Gate" for this directory.

- Correct addition: Eq. 15 is unexecuted (claim 7 above).
- Wrong: the phase values are "hard-coded" and time conversion is "forbidden by contract". They are configuration (`config.toml` lines 18–20). The paper calls them synthetic and the clock assumed; it forbids nothing.
- Wrong: "AstroRad plagiarism ... falsely attributed as a derivative insight from the preprint". The first message cited Science Advances and ESA for those figures and attributed only the shielding physics (claim 5) to the paper. This retracts a claim never made.
- Wrong by omission: it drops claims 2, 6, 11 and 12 and the unsourced 1.2e-4, and it drops the finding that claims 1 to 4 can still cite the paper once restated.
- Its draft file does not follow this directory's form (a claim as the title, the date in the filename, no author line, numbered claim headings) and uses words the paper does not ("priors", "industrial reactor", "production streams"). It was not adopted. This record replaces it.

## What can still cite the paper

- The dose-rate gap, the radiolysis kinetics and the offline counterfactual (claims 1, 2, 4), restated as above.
- The Compton argument, for photon transport only (claim 3).
- Nothing in claims 5, 10, 13 or 14, and neither Python engine.

## Reproduce

Every line number above was read with `sed -n`, not inferred from a match count. To re-resolve the numbers if the tex changes again:

```sh
cd preprint && for ph in 'nearly nine orders of magnitude apart' 'competing generation and recombination kinetics' \
  'dominated by Compton scattering' 'counterfactual question is logically prior' 'radiotrophic niche' \
  'specification, not method' 'operatorname{expm1}' '(0.01,0.002,0.015,0.04,0.02,0.06,0.03)' \
  'frequencies with mean' 'radiation forcing are off' 'no temporal phase for entrainment' \
  'Coherence is reported with active count' '26-neighbour adhesion stencil' \
  'from neural networks to circadian rhythms' 'FENE and Kelvin--Voigt mechanics, for which' 'label{eq:fene}'; do
  printf '%-50s %s\n' "$ph" "$(grep -n -F -- "$ph" modeling_radioresistance_and_radiotropic_fitness.tex | cut -d: -f1 | tr '\n' ',')"; done
```

Equation numbers: count `\begin{equation}` and `\begin{align}` in document order; the `align` at `eq:psde` has one `\notag` and so one number.
