# Tardigrade mechanisms remain future work

No protection coefficient is justified for this parcel model by the sources checked here. The executed synthetic survival design uses assumed sensitivities and no tardigrade-derived term. A proposed protection factor would modify biological response at fixed physical dose; it would not modify the dose field. A later damage-and-repair model need not reduce to a constant sensitivity multiplier.

Status: corrected research scope, 2026-09-30; not a production prediction. Correction record: `TARDIGRADE-AUDIT-02` in `data/claims_ledger.csv`. This replaces the earlier treatment of uncompleted audit items as established mechanisms.

## 1. Source checks constrain the claim

| Primary source | Endpoint and result | Limit for this model |
|---|---|---|
| [Clark-Hachtel et al. 2024, Figure 5B and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11078613/) | E. coli BL21-AI CFU survival after 2,180 Gy gamma irradiation improved with R. cf. varieornatus Dsup, but not the H. exemplaris ortholog. | One dose and one bacterial host do not establish a transferable protection coefficient. A negative audit that missed this panel cannot establish bacterial impossibility. |
| [Aguilar et al. 2025, Figures 1b and 3b](https://www.nature.com/articles/s41467-025-63652-3) | Yeast gained peroxide protection but became more sensitive to MMS, Zeocin and UV. Figure 3b also measures acute peroxide survival by colony counts. | The paper is not limited to spot plates. These chemical-stress assays do not provide an ionizing-radiation survival curve or Gy-based protection ratio. |
| [Diaz Escarcega et al. 2023](https://pubmed.ncbi.nlm.nih.gov/36858083/) | Expression of the tested constructs caused DNA-damage signals, chromatin condensation and toxicity in cultured primary rat cortical neurons without an imposed irradiation challenge. | Construct toxicity is host-specific evidence, not a universal effect of Dsup. |

The Aguilar paper names Zeocin in its results, Figure 1b and methods, but bleomycin in the discussion. This note follows the assay description and records the inconsistency. The acute assay uses 0, 4, 6 and 8 mM peroxide for 90 minutes, followed by colony counting; the approximately 25-fold statement refers to serial-dilution plates. Chromatin association tracked protection in the tested alleles; lowering the radical signal alone did not suffice. None of these results calibrates a microbial radiation hazard.

This is a bounded source check, not evidence that no relevant study exists. Whole-animal tolerance, human-cell damage readouts and chemical-oxidant survival cannot be substituted for a matched microbial ionizing-radiation response.

## 2. A protection contrast needs its own prediction

A possible future synthetic contrast is `p_i = -expm1(-alpha_s * q_i * delta_D_i)`, with `q=1` as control and a declared `0<q<1` as protection. There is no default or biological value of `q`. No such contrast is part of the fitness/phase ensembles.

Before implementation, freeze the contrast, values, seeds and outcomes. Require fixed-position selection controls, identified-parcel displacement, missing radii after extinction, unchanged physical dose and movement parameters, and exact recovery of the unprotected trajectory when `q=1`. A sensitising `q>1` would require an explicitly named interpretation. A biological value requires a matched ortholog, host, radiation regime and survival endpoint.

## 3. Uncompleted audit items stay unverified

TDR1, TRID1, CAHS, the activated-sludge bioindicator, NASA Bioculture cultivation, the cryo-EM application and whole-animal LD50 values remain an unverified research queue in this review. Their earlier mechanism claims are held from the manuscript. This status does not mean the studies are false or absent.

Repair would need a damage state, recovery law and timescale. Dormancy would need activity, viability and recovery states plus a declared stress input. Biomolecular phase separation would not, by itself, supply an oscillator phase-coupling law. These are requirements for future models, not claims that those proteins have supplied the missing parameters.
