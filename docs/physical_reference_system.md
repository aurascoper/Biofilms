# Physical reference system

**Date:** 2026-08-14 · **Status:** A0 partially ready, C evidence-only, B/D declared ·
**Ledger:** `data/parameter_provenance.csv` · **Source registry:** `data/sources.csv`

Nothing in this repository defines a physical system. The legacy R visualization uses
normalized `R = 1`, `L = 2`, so no value in it can be converted to centimetres by assigning a
unit, and the coupling config ships every physical key unset on purpose (audit §5). Before any
value can be sourced, the apparatus has to be frozen; before any value can be trusted, its
evidence basis, conditions, and *model compatibility* have to travel with it.

This document freezes what the reference systems are. The per-parameter values live in the
ledger, one row per `(reference_system_id, required_for_stage, config_key)`.

---

## 1. Provenance taxonomy

Two vocabularies that are often conflated, kept separate here because they answer different
questions.

**`system_provenance`** — what kind of thing the *reference system* is:

| Value | Meaning |
|---|---|
| `published_replica` | Every essential dimension and operating value comes from one experiment. You may say you reproduced that apparatus. |
| `certified_component` | One component is reproduced exactly (a certified source, a named membrane grade); the surrounding geometry is a standardized test geometry. |
| `engineered_composite` | Geometry, membrane, source, and biology come from different references. Scientifically legitimate, but it is a **new design**, not a reproduction. |
| `declared` | A modeling choice with no external referent. |

**`evidence_basis`** — what backs *one parameter*: `direct_measurement`, `assay_certificate`,
`manufacturer_datasheet`, `primary_literature`, `evaluated_nuclear_data`, `derived`, `declared`,
`synthetic`.

These compose. A photon emission rate can be `derived` from an `assay_certificate` inside an
`engineered_composite` — three independent facts that one column cannot carry.

**`model_compatibility`** — whether the implemented model can represent the value at all:
`direct`, `requires_transform`, `unsupported`. This is the column that stops a true value from
being used in the wrong model, and it is orthogonal to whether the value is known.

**`status` vs `mapping_status`** are also orthogonal. `status = ready, mapping_status = unmapped`
is a normal, expected combination: a completed seconds-per-MCS calibration or a measured
permeability law is scientifically ready while the software has no key to put it in.

---

## 2. Staged reference program

### A0 — idealized Cs-137 benchmark · `model_compatibility: direct`

> A sealed-source-free, idealized ¹³⁷Cs photon source (single 661.655 keV line, isotropic,
> on the cylinder axis) irradiating liquid water at 20 °C inside the modeled cylinder, with
> vacuum outer boundaries, run to validate geometry → tally → normalization only.

The only reference system the current code can build. It validates the transport chain — CSG,
lattice orientation, mesh congruence, spectrum handling, and the eV/source → Gy/s conversion —
with no biology and no membrane feedback.

Two declared simplifications, both recorded in the ledger rather than hidden:

- **Single line.** The evaluated data also give X-rays totalling ≈7.92 photons per 100
  disintegrations (Ba K-shell at 31.8–37.4 keV, L-shell at 3.95–5.81 keV) and two negligible
  gamma lines. Retaining only the 661.655 keV line is a declared simplification, not the full
  emission spectrum. The omitted photons are low-energy and strongly attenuated in water, so
  they matter near the source and can be added as extra spectrum lines when that is the point of
  the run.
- **No capsule.** Manufacturer capsule geometry is deliberately *not* part of A0 — see A1.

Two traps in the evaluated data that the ledger's `derivation` and `source_locator` columns
exist to prevent, both live in this one nuclide:

- The **transition** energy is 661.657(3) keV (§2.2 of the source) but the **emission** energy
  is 661.655(3) keV (§5.2) — they differ by the nuclear recoil. A photon source uses §5.2.
- The transition probability is 94.57(26) % (§2.2); the **photon** yield is 85.01(20) per 100
  disintegrations (§5.2). The difference is internal conversion,
  `P_γ = P_(γ+ce) / (1 + α_T)` with `α_T = 0.1124(16)`. Using 94.57 % as an emission
  probability overestimates the source by ~11 %.

Cs-137 is chosen over Co-60 for the first anchor because a single line is the smallest thing
that validates the chain. A **Co-60 anchor is worth adding second**, precisely because its two
lines and more-than-one-photon-per-decay yield exercise the `S_γ` vs `p_j` split that Cs-137
cannot.

### A1 — sealed-capsule benchmark · `model_compatibility: unsupported`

A certified sealed capsule cannot currently be modeled. `build_model()` offers only
`line_z_axis` (a discrete x,y with uniform z) and `point_origin` (a point at the cylinder
centre). There is no active core, no encapsulation, no capsule dimensions, no self-attenuation,
and no arbitrary placement or orientation.

Knowing a capsule's outer dimensions therefore does **not** make it model-ready. Implementing A1
needs a source-geometry contract roughly like:

```toml
[source]
kind = "sealed_cylinder"
position_cm = [...]
axis = [...]

  [source.active_core]
  radius_cm = ...
  length_cm = ...
  material = "source_core"

  [source.capsule]
  outer_radius_cm = ...
  outer_length_cm = ...
  material = "capsule_steel"
```

plus the corresponding CSG. Note that the **active core is a separate dimension from the outer
capsule**, and a product catalogue routinely gives the latter without the former. A catalogue
also gives an activity *range*, never a particular source's activity: that requires an assay
certificate with a reference date.

This is the next transport work item. It is not built here.

### B — Nafion N115 Donnan module · declared, unfilled

The fixed-membrane transport anchor: a named commercial grade with measured wet thickness,
known exposed area, explicit feed and draw compartments, and a described conditioning protocol.
Planar geometry, so it validates the membrane boundary condition against concentration-vs-time
data before returning to a cylindrical solver.

Its value as a lesson is already usable: a nominal datasheet thickness and a measured wet
thickness are *both correct under different conditions*. The ledger must therefore never carry a
single `membrane_thickness`; it carries the nominal value, the measured wet value, the ionic
form, the conditioning protocol, the temperature, and the external ionic strength as separate
rows.

### C — Syron–Casey MABR geometry anchor · `model_compatibility: unsupported`

Recorded as a **geometry evidence anchor only**. Two independent reasons it cannot be built:

1. **Inverted membrane topology.** The implemented model puts the lattice-filled biological
   domain inside `-inner_cyl` and the membrane in the annulus `+inner_cyl & -outer_cyl` — biomass
   inside, membrane outside. A MABR is the opposite: gas lumen, then the silicone tube wall, then
   biofilm on the tube's *exterior*, then bulk liquid, then the glass vessel. The code has no
   central gas lumen, no surrounding bulk-liquid annulus, and no outer vessel.
2. **Unrepresentable aspect ratio.** The lattice is cubic `(n,n,n)` with a single isotropic
   pitch, so the modeled domain is a cube of side `n·pitch`, while the cylinder radius and length
   are configured independently. A millimetre-bore, metre-long apparatus satisfies containment
   only at an infeasible voxel count. `check_lattice_congruence()` now rejects this rather than
   silently padding the cylinder with medium and tallying nothing there.

Three ways forward, none chosen here: keep it as an evidence anchor; declare a **representative
axial segment** with a stated segment length and axial boundary interpretation; or add
anisotropic pitch and an `inner_membrane_outer_biofilm` topology.

Its dimensions are also **not verified**. A search for the specific combination (3.0 mm OD
silicone tube, 1 mm wall, 1100 mm active length, 18 mm ID glass outer tube) did not locate a
primary source, and nearby lab MABRs report materially different geometries (e.g. 350 mm length,
4 mm OD / 3 mm ID). The stated 3.0 mm OD with a 1 mm wall would also imply a 1 mm bore, which is
not obviously consistent with those. Every C row is therefore `provisional` with a
`source_locator` of `NOT LOCATED`, pending the original Methods section.

### D — engineered radiotrophic composite · declared, unfilled

Where a certified source and a cylindrical biofilm reactor could legitimately meet — and the
**only** place they may. Combining A's source with C's reactor anywhere else produces exactly
the silent composite the ledger's uniqueness key
`(reference_system_id, required_for_stage, config_key)` exists to prevent: a ledger keyed on
`config_key` alone would either duplicate rows or quietly merge two incompatible systems into one
configuration and call it reproduced.

---

## 3. Fixed-membrane declaration

For the fixed-membrane stage the membrane is held constant:

```
m_mech = 1
P_j     = P_j0        (per solute j)
no dose-dependent membrane update
```

This resolves the modeling ambiguity by **removing the feedback path**, not by choosing between
the contradictory laws. The two variables are separated conceptually so the eventual choice is
well-posed:

- `m_mech(D, Ḋ)` — mechanical integrity. Calibrated against tensile modulus or strength vs dose.
- `P_j(D, Ḋ, T, hydration, pH)` — permeability of solute *j*, entering the boundary flux as
  `J_j = P_j (c_ext,j − K_j c_wall,j)`.

These are not interchangeable and need not move together: radiation damage can weaken a material
while opening transport pathways, and a membrane can become *less* selective, *more* permeable to
one species and *less* conductive to another. The original contradiction — `P_eff = P₀·exp(+α_P·D_cum)`
growing while `m` decays, with the nutrient coupling multiplying by `m` — is left described as-is
in audit §6, per that document's rule that contradictions are recorded rather than silently fixed.

### Validity envelope, split by material

The envelope is **material-specific and does not transfer**:

- **Nafion (references B, and any Nafion-based D).** The SRNL Nafion-117 Co-60 coupon study is
  the applicable benchmark. Its own design shows why total dose alone is an insufficient
  descriptor: it varied atmosphere (air vs sealed with deionized water) and dose rate
  (1 vs 460 krad/h) independently, and the degradation indicators depended on both.
- **Silicone (reference C).** There is **no applicable dose envelope**. Nafion irradiation data
  do not establish anything about a silicone tube. Holding the silicone membrane fixed is a
  *declared simplifying assumption with no experimentally validated dose envelope*. That is
  legitimate; calling it literature-validated would not be.

The STOP CONDITION is unchanged for any dose-responsive membrane model: OpenMC dose must not
feed the radiodialysis model until the `m`-vs-`P_j` constitutive choice is explicitly selected
and the membrane dose statistic is declared. See `docs/online_coupling_design.md`.

---

## 4. The transport mesh is not the biological pitch

The feasibility run reported 75 % of occupied bins unscored at the 10 µm demo tally pitch, with
median relative error ≈ 1 on the rest. That is a statement about tally statistics, not about
biology.

`geometry.voxel_pitch_cm` is therefore a **numerical parameter** to be chosen by a
resolution/history sweep. The CPM lattice pitch is a **separate, still-blocked** quantity: it
follows from what a CPM ID represents (a literal cell, a fungal compartment, or a coarse biomass
parcel) and from `a = (V_physical,s / V_target,s)^(1/3)`, which needs microscopy under the
intended growth conditions. The 10 µm feasibility mesh must never silently become the biological
pitch — at 120 sites it would imply a 120 000 µm³ target volume, a sphere about 61 µm across.

Because every metric in the sensitivity comparison is a ratio, the sweep can now run at the
`transport` stage in Gy per source particle, before any source activity or biological
calibration exists.

---

## 5. Why there is still no populated configuration TOML

Even A0 — the one `direct` reference system — cannot be built end to end, for a reason
independent of the staged loader: `required_classes()` hard-requires `medium`,
`baseline_biomass`, and `membrane` for every model, so a pure water phantom would need
fabricated biomass and membrane entries just to satisfy the biofilm builder.

**Do not add fake biomass to make it build.** Two honest routes, neither taken here:

1. a separate water-phantom model builder; or
2. a geometry-kind discriminator with stage- and geometry-specific required materials.

Until one exists, A0's `materials.baseline_biomass.*` and `materials.membrane.*` rows stay
`blocked`, and `source.photons_per_second` stays `blocked` pending an assay certificate with a
reference date. A per-source transport TOML becomes legitimate for A0 the moment the
water-phantom question is resolved; a full-coupling TOML remains out of scope.
