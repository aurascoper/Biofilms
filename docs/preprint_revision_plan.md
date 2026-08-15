# Preprint revision plan

**Subject:** `preprint/modeling_radiotrophic_fitness.tex` (source of record) and
`preprint/modeling_radiotrophic_fitness.md` (derivative, shorter).
**New title:** *Modeling Radioresistance and Radiotropic Fitness*.
**Governing records:** `docs/calibration/integration_contract.md` §"Which radiation phenotype this
model actually represents"; `docs/research/radiotrophic_compatibility_audit.md`;
`data/parameter_provenance.csv`.

This plan exists because the preprint sources are being revised in stages. This file records what
has landed, what is deliberately deferred, and what is currently a trap for a reader.

---

## 1. Live hazards, in priority order

### 1.1 The PDF is stale and must be rebuilt or withdrawn

`preprint/modeling_radiotrophic_fitness.pdf` is a 21-page build of the **pre-revision** `.tex`. It
carries the old title, the withdrawn radiotrophy claims, the fabricated Computational Methods
section, and the incorrect attribution of the radial stratification to `β_ion`. No LaTeX build was
attempted in the revision pass (no `pdflatex` on the machine), so the PDF was left in place rather
than half-updated.

A PDF that no longer matches its source is worse than a missing one: anyone who downloads it gets
the retracted claims with none of the corrections, and nothing in the file says so. **Before this
repository is shared or linked in any form, the PDF must either be rebuilt from the current `.tex`
or deleted.** Deleting it is the cheaper correct action; rebuilding it requires a LaTeX toolchain
with `booktabs`, `tabularx`, `longtable`, `rotating`, `pdflscape`, `microtype` and `lmodern`, plus
the four figures in `preprint/figures/`.

### 1.2 The file names still say `radiotrophic`

Both sources are still named `modeling_radiotrophic_fitness.{tex,md,pdf}`. The rename was
deliberately **not** done in the same commit as the rewrite, because a rename plus a full-file
rewrite produces a diff no reviewer can read.

Intended new basename: **`modeling_radioresistance_and_radiotropic_fitness`**
(`.tex`, `.md`; the `.pdf` follows whatever 1.1 decides). Do the rename as its own commit with no
content change, and update `graphicspath` only if the figures directory moves (it should not).

### 1.3 Posting status is unverified and must not be asserted either way

The header of both sources reads *"Preprint — submitted to bioRxiv (Systems Biology), Version 1.0 —
April 2026."* Two web searches on 2026-08-15 found **no record** of this preprint being posted, no
DOI appears anywhere in the repository, and the document's own final section is titled *"Notation
Fixes Required Before Submission"*, which is not what a submitted manuscript carries.

"No record found" is a search result, not a proof of non-posting. Do not write "not posted" as a
fact in any committed text. The header line was left untouched for that reason. Resolve it by
asking the authors, then state the true status once — with a DOI if one exists.

---

## 2. What landed in the phenotype-fold pass

Mechanical and non-deferrable items only. The section-by-section revision is §3.

1. **Retitle** in both files: `.tex` header comment, `\title` block, `pdftitle`; `.md` H1.
2. **Phenotype fold.** Every use of *radiotrophic*/*radiotrophy* that asserted radiation-enhanced
   growth was replaced with the word for the endpoint actually meant — *melanized*,
   *radioresistant*, or *radiotropic*. Occurrences that legitimately discuss radiotrophy as a
   contested literature claim were kept and marked contested.
3. **New subsection** (`§2.5` / `sec:radiotrophy_status`) recording that radiotrophy is not
   established for any of the seven modelled species, quoting the audit's four verdicts.
4. **Results correction.** The radial stratification is **melanin-mediated**, not `β_ion`-mediated:
   `ΔH_rad` biases Metropolis acceptance by 1.000010 for the radiotropic species against 1.155 for
   `ΔH_mel` at the observed `M = 1.44`. New table `tab:term_magnitudes`. Radiation still reaches
   the dynamics, indirectly, via `melanin_drive`.
5. **Computational Methods rewritten** to describe the code that exists. The old section named
   `DifferentialEquations.jl`, `JuMP.jl`, `Ipopt`, `Agents.jl`, `PlotlyJS.jl`, Sobol indices and
   Morris screening; a case-insensitive grep for each returns zero files outside the preprint
   itself. No optimisation and no sensitivity analysis was ever run.

Three results claims were **withdrawn in full** because no code in the repository computes them,
and each withdrawal is stated in-text rather than silently deleted:

| Withdrawn claim | Where it was | Why |
|---|---|---|
| `γ_s ≈ 0.07 / 0.06` are "positive gamma sensitivity coefficients" letting the fungi "exploit ionizing radiation as a metabolic stimulus" | §6.1 | In `biofilms.R` every `rad_sensitivity` enters as `-β·I₀`, a constant scalar with no spatial gradient, same sign for all seven species. 0.06 belongs to "New Extremophile", not *C. sphaerospermum* (0.05) |
| The kNN decision tree "correctly predicts species dominance transitions" above 1 kGy | §6.4 | The only k-NN in the repo (`reactor_decision_tree.R`) classifies four fuel isotopes by cost and energy density. No species, no dose axis, no phase variable. And no model here has birth or death |
| The motility/diffusion/sensitivity scatter "reveals" an inverse relationship | §6.3 | It is a scatter of Table 2's own entries. Reported as a finding it is circular |

---

## 2b. What landed in the adversarial verification pass (2026-08-15)

The phenotype-fold pass corrected the vocabulary and the headline attribution. This pass read the
revised sources against `README.md`, `data/claims_ledger.csv`,
`docs/calibration/integration_contract.md` and the compatibility audit, and found claims that
survived the fold because nobody had checked them against the code. All of the following are now
fixed in **both** sources unless noted.

**The dose gradient was reported backwards.** The CPM radiation field is maximal on the cylinder
*axis* and decays monotonically outward — pinned by `tests/deterministic_radiation.jl` — so the core
is the high-dose region and the wall is the low-dose one. The abstract ("concentrating in
high-radiation zones"), §6.5's melanin paragraph, both CPM figure captions and the Conclusion
("displaced toward the high-dose outer zone and *B. subtilis* retreating to the low-dose core") all
had it inverted. The fungi finish the run *further out*, which is the direction opposite to the one
their negative `β_ion` implements — a fact the corrected text now states rather than hides.

**The stratification was reported as established.** §6.5 said "it is not imposed as an initial
condition. That much stands", and the Conclusion led with it. Parcel centres are rejection-sampled
uniformly over `r < R − 4` = 26 = 0.867R at `N = 60`, `R = 30`, whose expected mean radius is 17.3
lattice units with a per-parcel SD of 6.1; at `n = 6` the SE on a species mean is ≈ 2.5, so the
4.5-unit gap is ≈ 1.3 SE of the difference and *both* species means sit within one SE of the
seeding mean (0.578R). The null is approximate anyway — parcels occupy volume and grow, so it must
be simulated, not derived. No multi-seed run against that null
exists. The radial figure is now a **diagnostic**, and replicate seeds are named as the cheapest
outstanding experiment. The melanin-dominance finding is unaffected — it is arithmetic on shipped
constants, not an outcome of the run — and the Conclusion now leads with that instead.

**§6.6 was quoting dimensionally dishonest numbers** (`.tex` only; the `.md` has no §6.6). "50 Gy
cumulative dose" against a `seconds_per_mcs` that ships `NaN` and a hard-coded placeholder
`Ḋ = 1.0`; "`P_eff = 0.027 cm s⁻¹`", which is the Nafion-117 literature prior `P₀` rescaled by a
model output. Both withdrawn in-text. `P_eff/P₀ = e ≈ 2.72` is retained and stated as *exact by
construction* from four hand-set constants. The "thin annular shell" depletion reading is withdrawn
(the biomass profile is averaged to a scalar, so the sink is spatially uniform), as is the
"% depleted" reading (`c(t=0) ≡ 0`, so it measures non-penetration). The transport solution itself is
kept and defended: the PDE is correctly derived, solved twice independently, and the two
implementations agree. This closes deferred item 7 below.

**The "self-regulating remediation loop" was reported as a result.** §6.6 opened with it and §7.3
elaborated it as a positive feedback loop. `m` and `P_eff` are independent open-loop functions of
time and neither reads the other; the retained STOP CONDITION is now stated in-text. It is
presented as the design hypothesis the framework exists to make testable, which is the honest and
more interesting version.

**§7.3 predicted an attenuation the model cannot compute.** "*O. intermedium* AM7 establishes a
protective biofilm matrix … that attenuates local dose rates sufficiently to permit *S. oneidensis*
colonisation" — the gamma field is initialised once and never updated, biomass is not a term in it,
there is no death term for dose to act through, and *O. intermedium* is absent from both Langevin
models, so the "emergent from the Langevin dynamics" attribution named a simulation that never
instantiated the organism. Restated as a hypothesis with the three gaps named.

**§7.4 asserted mechanical behaviour from unimplemented equations.** No FENE potential and no
Kelvin-Voigt element exists in any source file; no stress, strain or modulus is a state variable.
The protected-microniche mechanism is restated as motivation for the formalism.

**§7.1 and §3.4 asserted symplectic integration.** Nothing here performs it and no simulation
carries a momentum state. Demoted to stated design intent in both places, closing deferred item 5.

**§6.2 still claimed a comparison against observation.** "The Th-232 decay model compares actual
versus predicted radioactive decay" — there are no measured decay data in any `data/` file and no
decay code in the repository. Withdrawn; the constant-background premise, which is what the rest of
the paper actually uses, is kept.

**Provenance fixes.** The `audit2026` / [41] entry cited the audit at revision `f4e8dbf`, which
predates the audit's own commits; now `904b9a4` / `a03bc32` at repository revision `5d1b777`. §5's
"six cells per species" now reads "six biomass parcels", with the parcel semantics stated.

Corrected in `README.md` in the same pass: stale HEAD `f4e8dbf` (twice), a citation subtitle that
did not match the manuscript's, the `data/research/` ledger sizes (given as 180/155/133 total lines
where the data rows are 64/64/26), `data/a0_transport_convergence.csv` (117 lines, 100 data rows),
the calibration test-module count, the skip inventory (all four skips are coupling modules skipping
on `import openmc`; no calibration test skips), `\bibitem` count 44 → 47, and the claim that the
`.md` still lacks the CPM results. Corrected in `data/claims_ledger.csv` and
`docs/calibration/integration_contract.md`: the statement that `cpm.melanin_coupling` is "the only
row in that ledger not ranked `unknown`" — it is one of five, and the only one outside the A0
transport sweep, whose four ranked rows carry `sensitivity_domain = transport_numerical`.

---

## 3. Deferred to the section-by-section revision

Ordered by how much a reader is misled if it is left.

1. **Figure inventory does not match the text.** §6.1–6.4 reference "Figure 1" through "Figure 4"
   which do not exist as floats in the `.tex`; the four `\includegraphics` figures are the CPM and
   radiodialysis outputs used by §6.5–6.6. Either produce the Langevin figures or renumber.
2. **Table 2 rows that the audit refuted.** `α_M = 0.03–0.10 µg cell⁻¹ Gy⁻¹` for *A. niger* rests
   on no located measurement, and a per-cell basis is not one the CPM entity semantics admit. Both
   are now flagged in-text under Section 4; the rows themselves are unchanged and should be
   removed or re-based.
3. **The sign of `α_M`.** Eq. `eq:melanin` has melanin production increasing with dose. The only
   same-species same-modality data (LAC1/LAC2 down at 1 and 3 kGy) and the only
   pigmentation-versus-dose measurement (Bland 2022, pigmentation down under gamma) both run the
   other way. Flagged in-text; the equation is unchanged. Changing it is a model decision.
4. **The PSDE is specified but never integrated.** Eq. `eq:psde` governs a fitness field `F_s`. No
   simulation in the repository integrates it — the R codes advance positions only, and the CPM is
   a lattice Metropolis model. Sections 3.2, 3.4 and 3.7 should say which parts of the framework
   are implemented and which are specification.
5. ~~**Symplectic integration.**~~ **CLOSED** in the verification pass — demoted to stated design
   intent in §3.4 and §7.1. Implementing it remains future work, listed in §7.5.
6. **`H_pair`.** The old Methods listed a five-term CPM Hamiltonian; `compute_delta_H` has four
   (adhesion, volume, radiation, melanin). Methods now says four. Section 3.3's mapping should be
   reconciled with the code.
7. ~~**The RADIODIALYSIS results (§6.6)...**~~ **CLOSED** in the verification pass. The physical
   units are withdrawn in-text and the dimensionless quantities are retained with their status; the
   underlying units defect at `biofilms_potts.jl:1341-1342` and `:1445-1446`, and the hard-coded
   `X_red = 0.3` labelled "(Shewanella proxy)" that the audit's active-reducer finding argues
   against, are **code** defects and remain open. `RADIODIALYSIS: BLOCKED` stands.
8. **`k = 4` versus `k = 3`.** `biofilms.R` clusters with `centers = 4`; the CPM's own
   `kmeans_cells` uses `k = 3`. Fixed in Methods and §6.1; check every remaining mention.
9. **Lattice size.** Struct default is `N = 60`; the reported runs use `N = 40`. Methods now says
   both. Verify against whatever run actually produced the shipped figures.
10. **"Notation Fixes Required Before Submission"** (the `.md`'s final section) is a to-do list
    living inside the manuscript. Work it or move it here.
11. **The `.md` and `.tex` have diverged in content, not only in format.** The `.md` lacks the
    radiodialysis mathematics (§3.11 in the `.tex`) and, before this pass, lacked the CPM results
    entirely. Decide whether the `.md` is a maintained artefact or a generated one; if generated,
    generate it.

---

## 4. Citations added in this pass

- `cortesao2020` — *A. niger* N402 vs ΔfwnA LD₉₀ across three ionizing modalities. Front Microbiol
  11:560. Cited for the withdrawal of the *A. niger* radioprotection claim.
- `walberg2015` — UW-La Crosse thesis, `handle:1793/73406`. Closed-energy-budget critique of
  radiosynthesis. Licence not specified in the repository record; contains no new measurement.
- `audit2026` — the in-repo compatibility audit, cited by path and revision.

`\begin{thebibliography}{38}` now under-counts the entries. Harmless (it only sets label width),
but bump it when the bibliography is next touched.
