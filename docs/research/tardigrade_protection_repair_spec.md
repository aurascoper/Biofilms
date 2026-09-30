# Tardigrade-inspired protection, repair and dormancy

Tardigrade biology supports three distinct research mechanisms: damage prevention, damage repair and reversible metabolic suppression. The first can be tested as a synthetic survival contrast. The other two require new state variables. None changes physical absorbed dose, supplies an oscillator phase law, or justifies importing whole-animal radiation tolerance into a microbial parcel.

Status: research design, 2026-09-29; not a production preregistration. No tardigrade-inspired term has been implemented or exercised.

## 1. The evidence supports three different model changes

| Candidate | Evidence boundary | Model translation |
|---|---|---|
| Dsup | Dsup reduced X-ray damage in cultured human cells and protected chromatin from hydroxyl-radical damage. In *E. coli*, the *R. cf. varieornatus* ortholog increased survival after 2,180 Gy gamma irradiation, whereas the *H. exemplaris* ortholog did not. Heterologous expression also improved survival under oxidative damage in yeast. Effects depend on protein, host and stress; expression of the tested constructs in cultured primary rat neurons increased DNA damage and toxicity. | A fixed protected/unprotected survival contrast, calibrated only from a matched ortholog, host, exposure and endpoint. |
| TDR1 | Tardigrade DNA damage Response 1 is induced after ionising radiation, binds DNA, forms aggregates at high concentration, and reduced bleomycin-induced phospho-H2AX foci when expressed in human cells. A chromosome-organisation role during repair remains a proposal. | A candidate for a later damage-response model. It does not yet supply a repair rate. |
| TRID1 | A distinct tardigrade-specific radiation-induced disordered protein in *Hypsibius henanensis* was linked to phase separation and increased double-strand-break repair efficiency. | A stronger lead for an explicit damage-and-repair state with measured recovery kinetics. |
| CAHS proteins | CAHS proteins assemble under osmotic or dehydration-like stress. CAHS D condensation, metabolic slowing and improved survival were associated during osmotic stress, and the condensation and metabolic effect reversed after recovery. | A later dormancy branch with separate viability, metabolic activity and recovery outcomes. |

Primary studies: Hashimoto et al. 2016, <https://doi.org/10.1038/ncomms12808>; Clark-Hachtel et al. 2024, <https://doi.org/10.1016/j.cub.2024.03.019>; Anoud et al. 2024, <https://doi.org/10.7554/eLife.92621>; Li et al. 2024, <https://doi.org/10.1126/science.adl0799>; Sanchez-Martinez et al. 2024, <https://doi.org/10.1002/pro.4941>. Host dependence is shown by Diaz Escarcega et al. 2023, <https://doi.org/10.1016/j.mcn.2023.103826>. The yeast protection result is Aguilar et al. 2025, <https://doi.org/10.1038/s41467-025-63652-3>.

### 1.1 Source-to-claim matrix

| Source and construct | Host | Stressor and endpoint | Control | Demonstrated | Proposed or unresolved |
|---|---|---|---|---|---|
| Clark-Hachtel et al. 2024; *R. cf. varieornatus* Dsup and *H. exemplaris* Dsup | *E. coli* BL21-AI | One 2,180 Gy gamma exposure; CFU survival relative to the unirradiated culture of the same construct | Empty and GFP expression vectors; six replicates per Dsup construct | The *R. cf. varieornatus* ortholog increased survival; the *H. exemplaris* ortholog did not | No microbial dose response, D10, LD90, protection factor or transfer to another host |
| Hashimoto et al. 2016; *R. cf. varieornatus* Dsup | HEK293 cells | X rays; DNA damage and viability | Non-transfected cells | Reduced damage and improved viability in the tested cells | Transfer to microbial hosts or a survival coefficient |
| Aguilar et al. 2025; *R. varieornatus* Dsup | *S. cerevisiae* | Hydrogen peroxide; survival after acute oxidative stress | Uninduced or vector controls | Improved survival under the tested oxidative stress | Ionizing-radiation protection or a coefficient for the fungal panel |
| Diaz Escarcega et al. 2023; tested Dsup constructs | Primary rat cortical neurons | Expression without an imposed radiation survival assay; DNA damage and toxicity | Matched expression controls | Increased damage and toxicity in the tested neuronal system | A universal adverse effect of Dsup |
| Anoud et al. 2024; TDR1 | Human cells and *H. exemplaris* assays | Bleomycin or ionizing radiation; expression, binding and phospho-H2AX | Study-specific untreated and expression controls | Radiation induction, DNA binding and reduced bleomycin-induced damage signal | Chromosome preservation during repair and a repair rate |
| Li et al. 2024; TRID1 | *H. henanensis* and study assays | Ionizing radiation; phase separation and double-strand-break repair | Study-specific controls | Association with increased repair efficiency | Transferable recovery kinetics for a parcel model |
| Sanchez-Martinez et al. 2024; CAHS D | Cultured cells | Osmotic stress; condensation, metabolism, survival and recovery | Study-specific expression controls | Reversible assembly and biostasis in the tested system | Radiation protection, intact-organ preservation or an oscillator law |

Biomolecular phase separation in the TRID1 and CAHS studies is a material-state change. It is unrelated to oscillator phase locking and supplies no coupling law for the phase passenger.

## 2. First experiment: a Dsup-inspired protection contrast

Keep the exponential synthetic survival family and introduce a dimensionless biological damage-yield factor `q_i`:

`p_i = -expm1(-alpha_s * q_i * delta_D_i)`.

The unprotected control has `q_i = 1`. A protected parcel has a declared `0 < q_i < 1`. The factor changes biological damage attributed to a fixed physical dose increment; it never rescales the imported dose map, the physical clock, the radiation acceptance signal or the melanin field.

The first exercise is a paired protected/unprotected contrast with the same initial state, physical field, `alpha_s`, mechanical parameters and hazard seed. `q` has no default. A synthetic value may test the implementation, but it must be labelled assumed. A biological value requires a matched ortholog, host, radiation modality, dose-rate regime and endpoint. Clark-Hachtel et al. establish that a Dsup ortholog can protect a bacterial host, but the single-dose experiment does not yield a microbial D10, LD90 or `q`. Its 2,180 Gy acute exposure is also far above the 60 Gy synthetic axis exposure and exceeds several panel endpoints. The biological parameter remains a no-go until a matched dose response exists. Whole-tardigrade survival and cultured human-cell DNA damage cannot be inserted as a microbial parcel coefficient. The *S. cerevisiae* result is fungal-host evidence and is taxonomically closest to the basidiomycetous yeast *C. neoformans*, but it uses oxidative stress rather than ionizing radiation and remains farther from the filamentous fungi *A. niger* and *C. sphaerospermum*.

Report viable fraction, survivor-weighted radius and identified-parcel displacement under the existing missing-data rules. Also report damage commitments by dose bin. A fixed-position control must show that protection changes survival without creating movement. With `q = 1`, the extension must reproduce the unprotected hazard events and checkpoint bytes under the same hazard seed.

Required implementation controls reject a default-filled or nonfinite `q`, `q <= 0`, `q > 1` in the protection-only model, a code path that edits dose, and a path that changes mechanics before a hazard commitment. A later sensitisation study may allow `q > 1`, but it is a separate declared interpretation.

## 3. Repair needs an explicit damage state

Repair cannot be represented by silently lowering `alpha_s`. A repair model needs, at minimum, a per-parcel damage burden, an increment rule from current dose, a recovery law and timescale, a hazard or function threshold, and checkpoint state. Acute and protracted exposures with equal total dose must then be separate controls because repair can make their outcomes differ.

TDR1 supplies a radiation-responsive DNA-binding candidate, but its proposed chromosome-preservation role does not yield those parameters. TRID1 supplies a closer mechanistic lead through phase separation and double-strand-break repair, but its kinetics and host context still must be read from the experiment before a rate law is chosen. Protection and repair must be varied independently so that reduced initial damage is not misreported as faster recovery.

## 4. Dormancy needs activity and recovery states

A CAHS-inspired model needs an external dehydration or osmotic-stress variable, a transition into dormancy, a metabolic-activity state, survival during dormancy and a recovery transition. Viable, metabolically active and physically present are separate outcomes. The current CPM has none of these inputs or states, so CAHS does not enter the first survival extension.

CAHS gelation should not be called a universal protective bioglass. Hydrated condensates, fibrillar networks, gels and dry material states have different measurements, and protection can occur in non-gelling variants. Evidence for cultured-cell or molecular preservation does not establish room-temperature preservation of intact organs.

## 5. Whole animals and adjacent applications remain outside this model

Adding tardigrades as organisms would require an animal-ecology model with feeding, reproduction, development, locomotion and life-stage-specific stress responses. A microbial biomass parcel cannot stand in for a whole tardigrade.

Tardigrades have been cultured in NASA's Bioculture System on the International Space Station, which establishes controlled-hardware feasibility. It does not by itself supply radiation dosimetry or a response coefficient. *Thulinius ruffoi* in activated sludge is a proposed nutrient-load bioindicator that still needs experimental validation. The cryo-electron-microscopy precedent concerns the nematode AavLEA1 and truncated tardigrade RvLEAM protecting fragile complexes from air-water-interface damage during plunge freezing; it is a sample-preparation result, not an organ-preservation result.

## 6. Order of work

1. Add and validate the fixed `q` protection contrast without changing the dose field.
2. Obtain matched protection measurements before assigning a biological value.
3. Specify damage and repair only after a recovery time course supports its state and rate law.
4. Open a dormancy branch only when dehydration or osmotic stress is an explicit model input.
5. Keep whole-animal tardigrade ecology in a separate model.
