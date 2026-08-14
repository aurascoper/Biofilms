# Audit: coupling-readiness of the Biofilms repository

**Branch:** `feat/openmc-dose-coupling` · **Baseline:** `ec8929d` (clean master) ·
**Date:** 2026-08-13 · **Scope:** `biofilms_potts.jl` (serial CPM), `biofilms_potts_jacc.jl`
(GPU port), the three R models, preprint prose, and repository infrastructure.

Every claim below was verified against source at the baseline SHA; line numbers refer to
that revision. This document is the contract for the OpenMC coupling work: it records what
the code *actually* computes (as opposed to what headers, figures, or the preprint imply),
the contradictions that remain open, and the stop conditions derived from them.

---

## 1. The Hamiltonian: advertised five terms, implemented four

The file header (L7) and the coupled-section header (L1284) advertise

```
H_CPM = H_adhesion + H_volume + H_radiation + H_pairwise + H_melanin
```

but the Metropolis acceptance path `compute_delta_H` (L393–454) returns

```julia
return ΔH_adh + ΔH_vol + ΔH_rad + ΔH_mel     # L452
```

**The pairwise term is not part of the dynamics.** `total_pairwise_energy` (L687,
`E += -p.γ_mutual * exp(-r2 / p.σ_mutual^2)` at L698) is called exactly once, inside
`take_snapshot` (L767), and stored as a reporting statistic (`Snapshot.total_energy_pairwise`,
printed at L777). `γ_mutual`/`σ_mutual` appear nowhere in `compute_delta_H`. Cell–cell
attraction enters the dynamics only indirectly through the static adhesion matrix `J`.

Additionally, the melanin acceptance term hard-codes its coupling constant
(`0.5 * M_local`, L447/L450) rather than using the per-species `α_M_species` vector the
params struct declares — matching the JACC port's `MEL_COEF`, but another place where the
declared parameterization and the implemented one differ.

## 2. There is no biology to derive generations from

- **No division.** `next_id` appears only inside `init_state` (L288–323); after
  initialization it is dead state. `CellInfo(` is constructed at exactly one site (L304).
  Cells can only shrink and die (`volume <= 0 → delete!(state.cells, id)`, L525–529) —
  and death leaves **no lifecycle record**: the Dict entry is deleted, nothing is archived.
- **No growth.** `CellInfo.volume` changes only by ±1 on accepted copies (L516–521).
  `V_target` is a compile-time constant (L58) used only in the volume penalty (L420/426);
  under `λ_V = 10` a cell essentially never reaches `2·V_target`. Any volume-threshold
  division rule is therefore a *synthetic demonstration policy*, not an implementable
  biological cell cycle.
- **Nutrient is decorative.** `state.nutrient` is written by its own init/update functions
  (L239, L616/630, L1243/1259) and **read by nothing else**: not by `compute_delta_H`,
  not by any growth/target-volume rule (none exists), and it is never plotted.
- **`CellInfo` (L172–176)** carries only `species, volume, com` — no parent, lineage,
  generation, or birth-time fields.

Consequence for the coupling work: a CPM copy event must never be treated as mitosis, and
"generation" tracking requires new, clearly-labeled infrastructure (typed division/death
events, an archival cell table) with the division primitive manually invoked only.

## 3. Radiation semantics: the code is self-consistent, the figure is not

**Field** (`init_radiation!`, L220–226): `I(r) = I0 · exp(-κ·r/N)` with `r = radial_dist`
from the lattice axis — **maximum on the central axis**, decaying outward. Static after
init; dimensionless (`I0 = 1.0`, no Gy anywhere in the CPM half).

**Response** (`compute_delta_H`, L432–441): `ΔH_rad += β_ion[sp] · I_local` when a cell
gains a site. `β_ion` (L74–81) is **negative** for the radiotrophs C. neoformans and
C. sphaerospermum (energy *gain* in high flux) and strongly positive for S. oneidensis
(`7.5e-2`, "extremely radiosensitive"). So radiotrophs are driven **toward the axis**,
radiosensitive species outward. The section-13 header agrees (L1465–66: "radiotrophic
species drifting toward the axis").

**Figure 1 contradicts both** (L1568–1587): the hard-coded shaded zones label
"radiosensitive core" at r/R ∈ [0, 0.45] and "radiotrophic niche" at r/R ∈ [0.65, 1.0] —
the reverse of the implemented physics and of the section's own header. The committed
preprint figure `fig1_radial_stratification.*` embeds these reversed labels.

Also hard-coded: the Fig-3 annotation "(50 Gy cumulative)" (L1743) is literal text, not
computed from `Ddot_R · t` — and *cannot* be computed inside `export_figures`, whose
signature receives only `(trajectory, coupled_traj, params)`; neither `rp` nor `rd` is in
scope, and `CoupledSnapshot` carries no physical time or cumulative dose.

**Resolution adopted (user-approved):** the biological zone labels are replaced with
neutral physical labels ("inner high-radiation region" / "outer lower-radiation region")
and the ungrounded 50 Gy annotation is removed (commit 2). Committed preprint figure files
are left untouched; regenerating them after commit 2 will therefore differ from the
preprint versions. Deterministic sign tests (`tests/deterministic_radiation.jl`) pin the
implemented response direction.

## 4. The R 3D model's radiation force is broken (dead-code sign)

`biofilms_3d.R:159–174`:

```r
if (parameters$radiotrophic) {
  F_rad <- parameters$rad_sensitivity * I_r * (-radial_unit)   # :169
} else {
  F_rad <- -parameters$rad_sensitivity * I_r * radial_unit     # :172
}
```

The two branches are **algebraically identical** (scalar negation commutes), so every
species receives an inward radial force; the `radiotrophic` flag has no dynamical effect
and affects only plot legends/markers (:439, :450) — the figures label a stratification
the physics does not produce. This contradicts the file's own header (:10–12) and
README.md:101. The correct else branch is `+parameters$rad_sensitivity * I_r * radial_unit`
(fixed in commit 2). `assets/preview_bioreactor_3d.png` and preprint Figure 1 derive from
this code path; committed figures are not regenerated on this branch.

## 5. No physical unit contract exists

`RadiolysisParams` (L1050–1080) declares cm/s/g/Gy units, but every physically-labeled
slot receives a dimensionless or lattice-unit value:

| Declared | Actually supplied |
|---|---|
| radial domain R [cm] | `R = Float64(params.N)/2.0` = **20 lattice units** (L1295), straight into `r_grid` and the cm²/s diffusion operator |
| `X_total` [g cm⁻³], `X_red` [fraction of dry mass] | **occupancy fractions**: `compute_radial_biomass` (L1199) normalizes to fraction-of-sites-per-bin (L1223–25); `X_red` counts only `species == SO` sites. The radial profiles are then **collapsed by `mean()`** (L1330–31) before entering the solver — the radial resolution is computed and discarded |
| `dt_rd` [s] | `0.5` with the comment "s per CPM MCS **(tune to match CPM time scale)**" — an unset knob; the repo has **two independent clocks** (MCS count and radiodialysis seconds) with no declared conversion |
| `Ddot_R` [Gy s⁻¹] | constant `1.0`, never derived from `state.radiation` despite its comment |

Nothing in the repository defines: voxel pitch in cm, seconds per MCS, material mass
densities, elemental compositions, photon spectrum, physical source rate, coordinate
origin/axis conventions, or a dose→CPM-energy conversion. **All of these ship as
REQUIRED-but-unset config fields that fail loudly** (user decision); none may be invented.

The dose contract must also be *per-consumer*: the current dimensionless field has three
distinct consumers — the Hamiltonian (L432), melanin production (`update_melanin!`, L600),
and the membrane model's separate scalar `Ddot_R` — so a single conversion coefficient is
insufficient. Physical `dose_rate_Gy_s` stays a separate field from whatever normalized
signals feed each consumer.

## 6. Membrane: integrity never changes geometry, and m vs P_eff contradict

Membrane integrity `m` evolves as `dm/dt = -k_dam · Ddot_R · m` (L1126) and couples into
the CPM in exactly one place: `C_wall_eff = p.C_wall * m` (L1247–48). The cylinder radius,
`interior` mask, and wall sites (`lattice == -1`) are set once in `init_state` (L280–286)
and never modified. **Geometry is static**; the CPM `-1` region is "outside the biological
domain", *not* a finite-thickness membrane material, and must not be translated literally
into an OpenMC wall material (it would create solid corner wedges around the cylinder).

**Open constitutive contradiction:** the section header (L1038–39) claims "P_eff(t) scales
the wall boundary in `update_nutrient_coupled!`", but that function takes `m`.
`P_eff = P0 · exp(+α_P · D_cum)` (L1123) **grows** with dose while `m` **decays** — the
two are not interchangeable, and the intended membrane response is undeclared.

> **STOP CONDITION:** OpenMC dose must not be fed into the radiodialysis/membrane model
> until the intended m-versus-P_eff coupling is explicitly selected and the membrane dose
> statistic (mass-weighted volume average vs area-weighted surface measure) is declared.

**Narrowed by declaration, 2026-08-14.** The membrane is declared **fixed** for the
fixed-membrane stage: `m_mech = 1`, `P_j = P_j0` per solute, no dose-dependent membrane
update (`docs/physical_reference_system.md` §3). The two quantities are separated
conceptually — `m_mech(D, Ḋ)` is mechanical integrity, `P_j(D, Ḋ, T, hydration, pH)` is
per-solute permeability entering the flux as `J_j = P_j (c_ext,j − K_j c_wall,j)` — so the
eventual choice is well-posed rather than a choice between two spellings of one variable.

This **removes the feedback path; it does not resolve the contradiction.** Both laws above
still stand in the code exactly as described, and reactivate the moment feedback is wired,
so the stop condition is retained verbatim for any dose-responsive membrane model. The
validity envelope is material-specific and does not transfer: the SRNL Nafion-117 Co-60
coupon study applies to Nafion, and a silicone membrane has **no** validated dose envelope,
making its fixity a declared simplifying assumption rather than a literature-backed one.

## 7. Species registries are inconsistent across the repository

| Source | Roster | Encoding |
|---|---|---|
| `biofilms_potts.jl` L21–29 + `biofilms_potts_jacc.jl` :29–30 | CN, DR, CS, BS, AN, SO, OI — **identical, authoritative** | `RADIOTROPHIC = {CN, CS}` |
| `biofilms.R` :6–14 | includes placeholder **"New Extremophile"**; no O. intermedium | no radiotrophic flag; only `rad_sensitivity` (+ a one-off `rad_aversion`); DR given sensitivity 0.0 |
| `biofilms_3d.R` :14–72 | includes **Pseudoalteromonas** (nowhere else); omits O. intermedium | `radiotrophic` flag (CN, CS) — dead code, §4 |
| `biofilms_radiodialysis.R` | no registry (aggregate biomass, SO proxy comment :148) | — |
| preprint prose (`modeling_radiotrophic_fitness.md:200`, `.tex:580`) | names **Polaribacter** and **Flavobacterium** (a "marine taxa" k-means cluster) that **no simulation instantiates**, attributed to the flat-domain run that contains neither | — |

Orderings also differ between all three code registries, so any positional cross-file
mapping (colors, legends, parameter tables) is silently mismatched.

**Contract for the coupling work:** the Julia seven-species registry is authoritative;
the R models are legacy/exploratory. The preprint-prose discrepancy is recorded here, not
edited (preprint results are out of scope for this branch).

## 8. JACC port: statistically, not pathwise, equivalent — serial first

The GPU port's own header states its 8-color checkerboard dynamics are validated
*statistically* against the serial random-sequential model (CSV comparison across seeds),
not trajectory-identical. Serial-first coupling is therefore correct; the JACC coupling is
a design document only on this branch (`docs/jacc_coupling_port.md`, commit 10). Facts
that document must respect: device arrays are fixed-size at init (no growth); the RNG is
counter-based (splitmix64 keyed by site/MCS/color); `vols` reads inside a color pass are
deliberately stale (marked `ponytail:` at :167) — acceptable for a Metropolis energy,
**not** for a division trigger; event ordering within a color pass is nondeterministic, so
genealogy events must be keyed by `(site, mcs, color)` and sorted host-side.

## 9. Infrastructure baseline

At `ec8929d` the repository has: no `tests/`, no CI, no Python, no `docs/`, no declared R
dependencies, and one stale `.gitmodules` entry (`Lumped-Uncertainty-SLS-MPC`, not checked
out, referenced by nothing). The only assertion-bearing code is the JACC `--selftest`.
The literal comment line `#  13. Figure export` (two spaces) is a **load-bearing split
marker**: both `validate_serial.jl` and the JACC selftest load the serial monolith by
splitting the source at that string and evaluating the first half in a sandbox module that
imports only `LinearAlgebra, Statistics, Random, Printf`. Hence two hard rules for every
edit on this branch: never touch that comment line, and add nothing above it that needs a
non-stdlib package.

Environment at audit time: Julia 1.12 (JACC, AMDGPU); system Python 3.14 with **no**
OpenMC/h5py/numpy. OpenMC installation is planned (pinned conda-forge stack, commit 5);
until it succeeds, every OpenMC-dependent result is reported as blocked/skipped, never
passed.

## 10. Contradiction register (open items, not silently fixed)

1. Header/preprint five-term Hamiltonian vs four-term implementation (§1) — recorded.
2. Fig-1 zone labels vs implemented radiation response (§3) — labels neutralized in
   commit 2; committed preprint figures left as-is and now differ from regenerated output.
3. "50 Gy cumulative" annotation vs absence of any time/dose contract (§3) — removed in
   commit 2.
4. `biofilms_3d.R` radiation branches identical (§4) — fixed in commit 2; committed
   figures/animations derived from the old dynamics left as-is.
5. `m` vs `P_eff` membrane semantics (§6) — **narrowed by declaration (2026-08-14):**
   membrane fixed for the fixed-membrane stage, feedback deferred. Still **open in code**,
   and the formal stop condition is retained for any dose-dependent membrane model. See
   `docs/online_coupling_design.md` gate 1.
6. Occupancy fractions in g cm⁻³ slots; lattice R in cm slot; two clocks (§5) — open;
   closed only by the physical config contract (REQUIRED fields, fail loudly).
   Partially addressed 2026-08-14: requirements are now staged, so transport no longer
   demands the CPM clock or the biological response scales it never reads, and each key's
   evidence is tracked per stage in `data/parameter_provenance.csv`.
7. Preprint prose species vs implemented registries (§7) — recorded; out of scope.

## 11. Tests added with this audit (commit 1)

- `tests/contract_csv.jl` — byte-level regression: `validate_serial.jl 42` CSV output vs
  golden fixture generated at the baseline (pins the serial RNG stream; fixture header
  notes Julia-1.12 dependence).
- `tests/deterministic_radiation.jl` —
  (a) field shape: maximum on the axis, monotone radial decay;
  (b) exact local ΔH sign: with `J = 0`, `λ_V = 0`, melanin ≡ 0, a single-site cell
  gaining a high-I site has `sign(ΔH) = sign(β_ion[sp])` for every species, with magnitude
  decaying with radius;
  (c) multi-seed drift: with a deliberately **amplified synthetic β** (production values
  are far too small for a short-run single-seed test to be stable), the mean radial
  position of a radiotroph-signed cell decreases and a radiosensitive-signed cell
  increases, averaged across seeds.

Run: `julia --project=. tests/runtests.jl`
