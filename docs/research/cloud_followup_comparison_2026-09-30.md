# The cloud follow-up outputs match at printed precision

The 12 supplied Downloads files were read, hashed and checked against the existing local study. All 68,096 comparable CSV rows match in all ten printed fields after mapping arm names. This is agreement at six-decimal output precision, not proof of identical unrounded simulation states or independent authentication of the remote machine.

## 1. The numerical designs match; the predictions differ

Both producers scale the entire signed `beta_ion` vector or the entire `alpha_M_species` vector. Neither changes the guarded melanin coupling 0.5. They use N=40, six parcels per species and seeds 42–297, with the same parcel-centroid radius observable. The supplied files report Julia 1.12.7; local shards record 1.12.6. Both are x86_64 Linux. The reported model revisions are 24cd3ec and 66ffe9f; the local git model blobs were also compared.

| Study | Matched rows | Compared MCS | Other supplied rows |
|---|---:|---|---:|
| Beta | 28,672 | 0, 100 | 0 |
| Melanin production | 28,672 | 0, 100 | 0 |
| Horizon | 10,752 | 0, 100, 400 | 7,168 at MCS 200 and 300, absent from local archive |

The repository prediction defines an inward paired-mean crossing at −0.05 r/R, without a significance gate. The cloud prediction adds a sign test, different threshold predictions, and a CS primary. Its claim that failure at a million-fold beta multiplier proves a broken ladder is invalid: the positive control has a different configuration and observable. The cloud readout acknowledges that failure. Do not replace the repository prediction or import the cloud's “two failed predictions” scorecard into this paper.

The baseline agrees row by row. CN at alpha ×10,000 has mean −0.0814887890625 r/R with 11/256 positive differences. At ×1,000,000 it has −0.0031252265625, with 121/256 positive. Those values reproduce the local result.

## 2. Mean melanin does not resolve the lost crossing

At alpha factors 10,000, 100,000 and 1,000,000, the CN mean-melanin values are 16,199.79, 176,148.14 and 1,774,081.60. They continue to grow. The melanin update clips below zero but has no upper saturation cap. These species means contain no spatial derivative, so they cannot establish or reject flattening of the local field. Nor do they distinguish acceptance-probability saturation from geometric trapping. Those remain hypotheses requiring new spatial and acceptance diagnostics, not an explanation demonstrated by this CSV.

## 3. The proposed table correction needs narrower wording

All five bibliographic records were checked against Europe PMC/PubMed. The reconstructed titles match after typographic normalization; Dadachova 2004's missing title is now supplied from its record. The primary texts, not the citation strings, determine the following limits.

- **B. subtilis:** [Moeller 2007, Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC1855867/) gives dry-spore X-ray D10=838.1±98.0 Gy (mean±SD, n=3). This is not a vegetative value. The same paper points to Weller 2002 for stationary-phase cell X-ray experiments, so “no vegetative primary” is too broad. The replacement says no vegetative D10 is established by this experiment. A slight shoulder in the wild-type survival curve also prevents treating ln(10)/D10 as a demonstrated exponential law.
- **C. neoformans:** [Dadachova 2004](https://jnm.snmjournals.org/content/45/2/313) reports about 10% survival near 4 kGy. [Khajo 2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0025092) cites earlier work for lethality at 120/360 Gy and does not itself measure survival. [Dadachova 2008](https://doi.org/10.1111/j.1755-148X.2007.00430.x) reports low-dose CFU curves, a low-dose melanin protection factor and different dose-rate regimes. The chains are unresolved as a common calibration; the assays are not interchangeable. The 81.28% calculation at 360 Gy assumes an exponential through D10=4 kGy. It is not a model-free biological contradiction or grounds to assign a unique sensitivity.
- **A. niger:** [Cortesão 2020, Table 2](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2020.00560/full) reports conidial LD90=366 Gy in saline and 187 Gy when air-dried from saline. This does not calibrate hyphal biomass.

The remaining sensitivity notes remain unverified or explicitly proxy-based. An audit's failure to locate a primary is not proof that none exists. None of these changes alters the CPM tropism coefficients or the new ensemble's assumed hazard sensitivities.

## 4. Receipts and reproducible comparison

Input inventory: `receipts/downloads_manifest_20260930.json`. Comparison: `receipts/cloud_followup_comparison_20260930.json`. Bibliographic records: `receipts/species_sources_20260930.json`. The comparator is `diagnostics/radial_followups/compare_external.py`; its test rejects a changed value, a missing row and a duplicated row selected from the committed production CSV artifact. No supplied Julia program was executed as part of this comparison.

Local raw archive: `/home/aurascoper/Developer/Biofilms-runs/2026-09-29_66ffe9f/radial`. Supplied files remain in Downloads. Ledger records: `CROSS-MACHINE-01`, `MEL-COLLAPSE-01`, `SPECIES-NOTES-01`.
