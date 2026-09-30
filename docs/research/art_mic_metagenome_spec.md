# ART in radioactive environments: public-data study specification

This study asks whether array-associated reverse transcriptases (ARTs) occur in public metagenomes from wet nuclear storage or other radioactive environments. Its operational setting is microbially influenced corrosion (MIC). It does not treat an ART hit as radioresistance, host protection, phase locking or corrosion control.

Status: scoped study, 2026-09-29; not a production preregistration. Dataset discovery, accession freezing, profile recovery and power checks must be completed before the primary comparison is registered. No sequence outcome has been inspected for this project.

## 1. Why the first outcome is presence

The ART preprint describes a jumbo-phage RT family with an upstream repeat array and a downstream partner gene. It reports no radiation, stress, biofilm, growth or survival phenotype. Its tested hosts are *Staphylococcus lentus* during phage SA1 infection and *Escherichia coli* as an expression chassis; other family members include *Listeria* phages and environmental viral contigs. None is one of the seven organisms in the CPM.

The first useful question is therefore whether an ART-like locus is detectable in the environments the next model would represent. Public sequence data can answer that question without sampling a reactor pool. MIC supplies the applied reason to ask it: a phage system could alter phage success, infected-cell loss or competition among phages, which could change the composition of a surface biofilm. The beneficiary and sign cannot be inferred from the locus.

The study also audits the current seven-species panel. That panel was assembled from organisms with usable dose-response literature. It is not a measured spent-fuel-pool community. Taxa detected in the selected pool datasets, with the limits of their sequencing methods, define the candidate list for a later pool-specific model.

## 2. Candidate data and matched controls

Freeze every accession, sample description, library type, read layout, sequencing platform, nominal radiation or contamination measure, sampling surface or matrix, and available negative or background sample before assembly. The first confirmed public pool set is MG-RAST 4698718.3, 4698717.3, 4698715.3 and 4698714.3: spent-fuel pool, fuel-transfer channel, dry drain and wet drain samples from one nuclear power plant. These are bulk shotgun metagenomes, not a dedicated virome, and the paper reports low radiation at the sampled wall.

Other radioactive-environment candidates may include uranium-contaminated sediment or groundwater projects such as PRJNA566790. They are a separate stratum from wet fuel storage and cannot be pooled into one treatment label.

For each radioactive sample, select controls in this order:

1. a paired background or neighboring sample from the same study;
2. the same matrix, extraction method, platform, read length and comparable depth from a different study;
3. a matched subset selected without looking at ART outcomes.

Record failed matches. A dry drain is a useful within-study comparator but is not automatically non-radioactive or chemically matched. The primary analysis must be stratified by study and environment; it cannot treat all radioactive sites as replicates of one condition.

## 3. Detection pipeline

Use the ART paper's family search and architecture test as the reference procedure.

1. Recover the authors' ART profile and reference sequences if released. Otherwise rebuild a profile from the 12 named GenBank phage RTs, recording every accession and alignment command. Search predicted proteins with HMMER and retain both raw scores and the pre-registered ART threshold. The paper used `hmmsearch` at E <= 1e-5 for family discovery and later required 207 bits for one partner-centered ART call; choose one use before production and validate it against ART, other phage RT, retron and diversity-generating retroelement references.
2. Assemble each sample with a fixed short-read pipeline. Run a declared alternative assembly only as a sensitivity analysis. Call viral contigs with one pinned classifier and retain unclassified contigs for the initial RT search.
3. For every candidate RT, recover up to 3 kb upstream and the downstream gene neighborhood. Run the paper's repeat scan: recurring short words, at least three copies, regular 100--450 nt start-to-start spacing and a shuffle-based excess test. Confirm final arrays with the length-free delimitation procedure or a documented reproduction of it.
4. Search the first downstream co-oriented genes against the three ART partner-family profiles. Record the RT, long N-terminal region, catalytic `[YF].DD` motif, array, partner and viral-classification evidence separately. No single motif is sufficient for an ART call.
5. Re-search raw reads around every candidate and every contig-end case. Targeted assembly may extend a locus; it cannot turn a failure to extend into a biological absence.

Run the full pipeline first on withheld known ART loci and near-neighbor phage RT, retron and diversity-generating retroelement loci. Downsample and fragment the positive controls to each study's read-length, depth and assembly-contiguity distribution. Freeze a minimum recovery criterion before reading environmental group labels.

## 4. Every candidate has one of three architecture outcomes

The primary unit is an ART-like RT locus, not a sample-level yes/no flag. Give each locus one state:

| State | Rule |
|---|---|
| Architecture present | ART-family RT plus a qualifying upstream array; report downstream partner evidence separately |
| Architecture absent and assessable | no qualifying array with at least 1.5 kb assembled upstream of the RT |
| Unassessable | less than 1.5 kb upstream, an RT or array at a contig boundary, or failed recovery below the frozen sensitivity criterion |

The 1.5 kb rule follows the ART paper. It searched up to 3 kb upstream and labelled no-array loci with less than 1.5 kb upstream as not assessed. In its broader RT-neighborhood census, 85.9% of recovered 20 kb windows were truncated by a contig end. A negative result is interpretable only after the assessable fraction and spike-in recovery are reported.

Primary descriptive outputs are ART-family RTs per assembled gigabase, architecture-confirmed loci per assessable ART-family RT, the unassessable fraction, and exact intervals by study. Environmental contrasts use matched sets and assembly-quality covariates fixed in the preregistration. A result with no confirmed ART loci but inadequate recovery is inconclusive, not a negative finding.

## 5. Claims the study cannot make

An ART locus does not show that its RT is active, that array RNAs are substrates, that the RT binds its partner, or that the system benefits a bacterial host. The preprint establishes a conserved RT catalytic motif, an unusual long N-terminal region, repeat arrays and partner-family association. It reports abundant array-derived candidate RNAs, but no RT assay or biological function.

Before ART enters a biofilm equation, four items are required: presence in the relevant host--phage community; a causal phenotype; the beneficiary and sign; and measurements that support a rate law and a model observable. Host protection, infected-cell growth arrest and competition between phages require different state transitions.

DRT2 is a useful precedent, not evidence about ART. DRT2 reverse transcription can assemble a toxic repetitive gene and cause abortive infection. That example shows why infected-cell survival and population protection must be separate outcomes.

## 6. Work products before production

Commit an accession manifest, source and control selection rules, recovered profiles and reference sequences, container or environment lockfile, assembly-quality report, positive/near-neighbor recovery benchmark, PREDICTION file, pipeline tests and an analysis that refuses unassessable loci. The first production run starts only after those artifacts share one revision.

Pool sampling is outside this study. Any later reactor or spent-fuel-pool sampling requires the facility's approvals, radiological controls and sampling authority.
