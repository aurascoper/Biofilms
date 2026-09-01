# Radiotrophic compatibility audit

**Date:** 2026-08-15 · **Revision audited:** `f4e8dbf` · **Scope:** evidence audit, no calibration ·
**Populated:** nothing · **Verdicts:** four, at the end ·
**Vocabularies:** `docs/calibration/acquisition_contract.md`, `calibration/biofilm_calibration/schema.py`

This document records what public evidence exists for the radiation phenotypes the model already
names, what modality that evidence was measured under, whether any of it is retrievable as data
rather than as a figure, and where it contradicts the assumptions the repository currently carries.
It populates no parameter. Every quantity it reports is attributed to the record that measured it,
and every quantity it could not verify is recorded as unverified rather than estimated.

The audit was run as seven parallel retrieval tracks, each independently re-verified by a second
pass against the primary sources. Where the two passes disagreed, the verification pass governs and
the disagreement is recorded rather than silently resolved.

## 1. What this document is not

It is not a calibration. `lattice_pitch_um`, `density_g_cm3`, `seconds_per_mcs`,
`normalization.hamiltonian_scale` and `normalization.melanin_scale` are untouched, and no material
composition is proposed. Three separate findings below would have supplied a plausible-looking
number, and each is recorded as a refusal with its reason.

It is also not a proof of absence. Every "no candidate found" statement below is a search result on
one date across a named set of repositories. Institutional repositories, corresponding-author
servers, laboratory GitHub pages and dissertation appendices were sampled, not exhausted.

## 2. Repositories and databases searched

Queries are given verbatim where the exact string matters. Hit counts are as returned.

### Data repositories

| Repository | Query | Result |
|---|---|---|
| Dryad API v2 | `/search?q=Cladosporium` | 12 datasets, all plant-pathology or mycobiome, zero radiation |
| Dryad API v2 | `/search?q=Cryptococcus` | 8 datasets, all genomics/clinical/ecology |
| Dryad API v2 | `/search?q=Cryptococcus neoformans biofilm` | `total 0` |
| Dryad API v2 | `/search?q=Aspergillus niger` | 10 datasets, none radiation, imaging or composition |
| Dryad API v2 | `/search?q=Shewanella` | 7 datasets in the whole corpus |
| Dryad API v2 | `/search?q=Ochrobactrum intermedium` | 0 |
| Dryad API v2 | `/search?q=melanin radiation` | 14 datasets, all animal pigmentation |
| Dryad API v2 | `/datasets/doi:10.5061/dryad.8sc52` + `/versions/2048/files` | full 7-file manifest, digests, sizes |
| Dryad API v2 | `/datasets/doi:10.5061/dryad.f4qrfj71n` + `/versions/242440/files` | full 9-file manifest, SHA-256 for all nine |
| Zenodo REST | `q="Cladosporium sphaerospermum"` | 45 records, all taxonomy/phylogeny |
| Zenodo REST | `q="Cryptococcus neoformans"` | 1455 records, overwhelmingly MycoKeys taxonomic treatments |
| Zenodo REST | `q="Aspergillus niger" AND (biofilm OR confocal OR irradiat*)&type=dataset` | 0 |
| Zenodo REST | `q=Deinococcus` | 84 records, the datasets are protein-crystal diffraction images |
| Zenodo REST | `q=melanized AND fungi AND radiation` | 7 records, triaged; one usable |
| Zenodo REST | `/records/21921798`, `/records/14994137`, `/records/3667494`, `/records/4944335` | full manifests, checksums, licences |
| Figshare API | POST `/v2/articles/search {"search_for":"Cladosporium sphaerospermum radiation"}` | 2 hits, both Frontiers supplements |
| Figshare API | `{"search_for":"Deinococcus radiodurans biofilm"}` | 2 hits, a DOCX table and a PDF data sheet |
| Figshare API | `{"search_for":"Bacillus subtilis biofilm confocal z-stack"}` | 0 |
| EBI BioImage Archive | `/api/v1/BioImages/search?query=Cryptococcus` | `totalHits: 0`, `isTotalHitsExact: true` |
| EBI BioImage Archive | `?query=Shewanella`, `?query=Deinococcus` | 0; 1 (a BamA cryo-EM single-particle set) |
| EBI BioImage Archive | `?query=biofilm` | 11 collections, enumerated in full |
| EBI BioStudies | `?query=Shewanella oneidensis biofilm` | 16,574 hits, every one an `S-EPMC` literature mirror |
| Image Data Resource | `/searchengine/api/v1/.../searchvalues/?value=Cryptococcus` | 0 images, projects, screens, plates, wells |
| Image Data Resource | `/api/v0/m/projects/?limit=500` enumerated, regex `biofilm\|subtilis\|Deinococc\|bacteri` | 0 of 147 projects |
| EMPIAR | `/empiar/api/entry/search/biofilm/` | `{"detail":"Not found."}` |
| NASA OSDR | `/osdr/data/search?term=Cladosporium` | 6 studies, none the ISS radiation experiment |
| NASA OSDR | `?term=Aspergillus niger` | 3 studies (OSD-132, OSD-260, OSD-350), all sequencing |
| NASA OSDR | `?term=Shewanella` | 91 records, all omics |
| heiDATA (Dataverse) | `/api/datasets/:persistentId/?persistentId=doi:10.11588/data/KH6NDD` | 31 TIFF tomograms, 13,799,884,261 B |
| depositonce (DSpace 7) | `/api/pid/find?id=hdl:11303/19659` → bundles → bitstreams | 10 files, 46,491 B, per-file MD5 |
| MINDS@UW | `handle:1793/73406` | thesis record, 1.49 MB, licence not specified |

### Sequence and functional-genomics archives

| Archive | Query | Result |
|---|---|---|
| NCBI GEO (eutils `db=gds`) | `Cryptococcus neoformans AND radiation` | 7 Series; only GSE80230 and GSE117227 are ionizing-radiation experiments |
| NCBI GEO | `Exophiala dermatitidis radiation` | exactly 2 Series: GSE142318, GSE152116 |
| NCBI GEO | direct fetch GSE80230, GSE117227, GSE142318, GSE152116, GSE294405 | full landing pages, sample tables, file sizes |
| NCBI BioProject | `Cryptococcus neoformans AND radiation` | PRJNA481510, PRJNA481427, PRJNA318361 |
| NCBI BioProject | `Exophiala dermatitidis AND radiation` | PRJNA596544, 635404, 638366, 928362, 928692 |
| NCBI BioProject | `melanin AND gamma radiation` | PRJNA596544 only — an *Exophiala* project |
| NCBI nuccore / biosample / sra / gene / assembly | `AM7 AND (Ochrobactrum intermedium[Organism] OR Brucella intermedia[Organism])` | 1 / 0 / 0 / 0 / 0 |
| NCBI assembly | `Ochrobactrum intermedium AM7`, `Brucella intermedia AM7` | 0 each; 96 assemblies exist for taxid 94625, none is AM7 |
| NCBI SRA | `Cryptococcus neoformans[Organism] AND melanin` | SRP199743 / PRJNA545271 (melanin biosynthesis, no radiation) |

### Literature and metadata services

Europe PMC REST (`resultType=core`, raw `abstractText` requested) was the primary metadata route,
used specifically because PubMed's cookie interstitial and several publisher sites refuse automated
clients. Queries included `"Cladosporium sphaerospermum" AND (gamma OR "ionizing radiation" OR
radionuclide OR Chernobyl OR radiotrop*)` (114 hits), `"Exophiala dermatitidis" AND (radiation OR
irradiation)`, and DOI- and title-scoped lookups for every candidate. Crossref, OpenAlex and
DataCite were used to pin identifiers, licences and open-access status. Full text was read from PMC
where available: PMC9294542, PMC9287308, PMC5137497, PMC174092, PMC3178601, PMC4278410, PMC4692441,
PMC7146846, PMC7182175, PMC7772362, PMC7835796, PMC9486643, PMC10469600, PMC10581076, PMC10926633,
PMC11567843, PMC11788442, PMC12941786, PMC1866175, PMC1347324.

### Controls encountered and not defeated

| Control | Where | Handling |
|---|---|---|
| Dryad Anubis proof-of-work | `datadryad.org` download endpoints | not probed; `browser_manual` recorded. For `10.5061/dryad.8sc52` a legitimate second host (Zenodo 4944335) serves byte-identical files, verified by MD5 against the Dryad API |
| Dryad API `401` | `/api/v2/datasets/doi:10.5061/dryad.8sc52/download` | recorded, not worked around |
| Anubis v1.26.2, user-agent dependent | heiDATA | a plain `curl` UA returned HTTP 200 with `restricted: false`; a markdown-fetching client was challenged. No challenge defeated. Access recorded as verified-but-UA-dependent |
| HTTP 403 | bioRxiv landing pages and PDFs (three separate preprints) | not retried; abstracts read from independent mirrors, values marked unverified at source |
| HTTP 403 | ScienceDirect, Elsevier, ACS, OUP | abstracts read via Europe PMC; full texts recorded as unread |
| HTTP 402 | Wiley (`10.1111/1462-2920.14550`) | recorded as `metadata_only` |
| HTML bot-wall on a PDF endpoint | PMC174092 and the ASM mirror for Wang 1996 | the article is a scanned page image; `access_status` corrected to `browser_manual` and the numbers behind it marked unverified |
| Cookie interstitial | PubMed | routed to Europe PMC REST, the intended machine route |

## 3. Exact-species candidates, per species

The seven modelled species are `biofilms_potts.jl:21-31`. Coverage is uneven to the point of being
the finding.

### Cryptococcus neoformans

| Record | What it measures | Modality | Access |
|---|---|---|---|
| GSE80230 (Jung/Kwon 2016) | transcriptome, 3 kGy, 30/60/120 min recovery; spot-dilution survival; LAC1/LAC2 at 1 and 3 kGy | 60Co gamma | `direct_download`, 25,835,520 B verified against `filelist.txt` |
| GSE117227 (Rad53/Chk1) | RNA-seq 30 min after 0.5 kGy, H99 and *rad53Δ* | gamma, isotope unstated | `direct_download`, 850.0 kB |
| Khajo 2011 (PLOS ONE) | EPR, 263 nm breakdown products, TBARS, 213Bi binding after 120 and 360 Gy at 11.94 Gy/min | 137Cs | `direct_download`, no deposit |
| Malo 2018 (Fungal Biol 122:449-456) | ultrastructural damage, melanized vs non-melanized | gamma, deuterons, alphas | `metadata_only`, Elsevier |
| Schultzhaus 2019 (Environ Microbiol) | transcriptome, melanized vs non-melanized | gamma | `metadata_only`, `isOpenAccess: N`, no accession located |
| Dadachova 2007 (PLOS ONE) | CFU, dry weight, 14C-acetate under a 188Re field | 188Re beta at 0.05 mGy/h (cells); 137Cs 14 Gy/min (isolated melanin) | `direct_download`, no deposit |
| Wang 1996 (Infect Immun) | melanin 15.4% of dry cell mass, ATCC 24067, 10 d, 1.0 mM L-dopa | none | `browser_manual`, scanned pages |
| Ravi 2009 (Mycopathologia) | biofilm thickness ~17.0 µm at 37 °C, ~16 µm at 25-35 °C, 48 h | none | `direct_download`, rendered figures only |
| Zenodo 21921798 | solid-state NMR of *Cryptococcus* cell walls, 143,676,256 B, CC BY 4.0, deposited 2026-08-13 | none | `direct_download`, contents uninspected |

Repository sweep for volumetric *C. neoformans* biofilm imaging returned zero across IDR, BioImage
Archive and Dryad, each independently reproduced.

### Cladosporium sphaerospermum

| Record | What it measures | Modality | Access |
|---|---|---|---|
| Averesch 2022 (Front Microbiol 13:877625) | 147 vs 151 CPM under a Petri dish over 622.5 h; relative optical density from camera brightness | mixed low-Earth-orbit | `direct_download`, CC BY 4.0, 14,821,993 B of supplement |
| Vember & Zhdanova 2015 (Mikrobiol Zh 77(5):70-80) | survival of conidia and mycelial fragments at 10,080 ber, strains 5-1, 60, 852, 3176 | 90Sr/90Y beta under 0.05 mm Al | `direct_download` — a free 469,611 B publisher PDF, in Ukrainian |
| Zhdanova 2004 (Mycol Res 108:1089-1096) | directed hyphal growth quantified as mean return angle; *C. sphaerospermum* 3176 responded | 109Cd, directional | `inaccessible`, Elsevier paywall |
| Dadachova 2007 | colony diameter and radial growth rate, tricyclazole contrast | 188Re beta | `direct_download` |

### Aspergillus niger

| Record | What it measures | Modality | Access |
|---|---|---|---|
| Cortesão 2020 (Front Microbiol 11:560) | LD90 by CFU: N402 366 Gy X-ray, 506 Gy He, 112 Gy Fe; *ΔfwnA* 353 / 567 / 112 | 200 kV unfiltered X-ray; 150 MeV/n He (LET 2.2); 500 MeV/n Fe (LET 200); UV-C | `direct_download`, CC BY 4.0, 190,728 B PDF |
| J Radiat Res 65(1):28-35 (2024) | germination 69% against 5.5% colony formation at 0.4 kGy; OD600 every ~10 min; hourly time-lapse | 60Co, 1.87 kGy/h | `browser_manual`, supplement subscription-gated |
| Front Microbiol 14:1233740 (2023) | 81.13 ± 0.27 J/m² transmitted at 1038 J/m² incident through pyomelanin in KOH | UV-C 254 nm | `direct_download` |
| depositonce-18458 (Barthel thesis) | cell-wall chitin 214.1 ± 3.4 µg glucosamine per mg extracted wall, 11 strains, n=3 flasks | none | `direct_download`, CC BY 4.0, 46,491 B |
| Front Microbiol 13:975763 (2022) | paired wet and dry masses of colony-carrying filters; mycelium thickness 59 ± 1 µm | none (clinostat) | `direct_download` article; raw data author-request |

### Shewanella oneidensis

| Record | What it measures | Modality | Access |
|---|---|---|---|
| Qiu 2006 (J Bacteriol 188:1199) / GSE3876 | 170 induced and 87 repressed genes over 1 h recovery; D20 = 40 Gy in Davis, D10 ~70 Gy in TGY | 60Co, 109 Gy/min, on ice | `direct_download` |
| Brown 2015 (PLOS ONE 10:e0131249) | D10 ~84 Gy, D20 ~57 Gy; lag extended to ~7.5 h at 95 Gy; FT-IR and MALDI fingerprints; Fe(III) reduction more than doubled at 50 Gy | Faxitron CP-160, 160 kV, 6 mA, tungsten; 0.79 Gy/min by Fricke | `direct_download`, CC BY, no deposit |
| Dryad `10.5061/dryad.8sc52` | SERS chemical maps of an Ag-precipitating biofilm at 6, 9 and 35 d, 0.5 µm/px, 1600 × 63,000 per scan | none (532 nm probe) | `browser_manual` at Dryad; `direct_download` at the Zenodo mirror |
| Dryad `10.5061/dryad.p8cz8w9rw` | ZFC/FC magnetometry, XRD, FTIR of magnetite before and after MR-1 reduction | none | `browser_manual`, 768,305 B current version |
| Zenodo 14994137 (PEDOT:PSS) | wettability, roughness, FTIR, electrochemistry, LIVE/DEAD; Olympus `.oir` at 0.2072 µm/px | none | `direct_download`, 64,402,398 B, CC BY 4.0 |

### Bacillus subtilis

No ionizing-radiation record was located for *B. subtilis* in this audit. What exists is morphology:

| Record | What it measures | Access |
|---|---|---|
| S-BIAD474 | colony biofilm volumetric imaging, wild type plus six matrix mutants, 48 h and 72 h; z-step 3.87 µm, 1024² over 890 µm | `direct_download`, 50 files, ~1 TB; Fig 7 volumetric component is six files at 1.8-4.3 GB |
| heiDATA `10.11588/data/KH6NDD` | 31 soft-X-ray tomograms, label-free hydrated, 530 eV, 92 projections | `direct_download`, 13,799,884,261 B, CC BY 4.0 |
| S-BSST261 | biofilm edges at 24 h and 5 d, Leica `.lif` | `direct_download`, 18,529,647,313 B |
| Dryad `10.5061/dryad.f4qrfj71n` | confocal of 20 µm cryosections, dual transcriptional reporters, 48 h | `browser_manual`, 49,717,304,106 B current version, CC0-1.0 |
| Zenodo 3544513 | σB pulsing confocal, "data extracted from confocal microscopy" | `direct_download`, 2,454,298,484 B |

### Deinococcus radiodurans

No record in this audit applied ionizing radiation to *D. radiodurans*. The species' radioresistance
is established elsewhere in the literature; it is not established by anything catalogued here, and
the label is not imported on reputation. What exists:

- Shukla et al. 2022 (*J Appl Microbiol* 133:796-807) — CLSM of a **recombinant** biofilm-forming
  derivative carrying `gfp` and `kanR`, enhanced by 25 mM Ca²⁺. The abstract states verbatim:
  *"Deinococcus radiodurans R1 is known as a nonbiofilm former bacterium."* No image deposit exists.
- Manobala et al. 2020 (*Chemosphere* 269:128722) — uranium biosorption kinetics against biofilm
  dry mass; XRD indicates uranyl phosphate deposits. Uranium is a chemical sorbate here; nothing
  was irradiated. The strain is not named in the abstract.

An exhaustive sweep for volumetric *D. radiodurans* biofilm data returned nothing across Dryad,
Zenodo, BioImage Archive, BioStudies, EMPIAR, IDR, OSF, Figshare and PubMed. The *D. radiodurans*
arm of any spatial campaign has no public fallback and must be acquired.

### Ochrobactrum intermedium AM7

The entire public footprint of this strain is **1086 base pairs and one paywalled abstract.**

- GenBank `MK268754.2`, 1086 bp partial 16S, `isolation_source` = "soil around Kakrapar Atomic
  Power Station", `country` = "India: Surat".
- Shukla et al. 2020, `10.1016/j.jhazmat.2020.122047` — thorium tolerance at 1000 mg L⁻¹, EPS yield
  unchanged under Th, FTIR functional groups. OpenAlex: `is_oa: false`,
  `any_repository_has_fulltext: false`, `best_oa_location: null`. Cited 36 times, with no follow-up
  work on AM7 by any group in six years.

Zero genomes, BioSamples, SRA runs, proteomes, transcriptomes, images or deposited datasets. No
culture-collection accession was located, so whether the strain is obtainable at all is unverified.
AM7 has never been exposed to ionizing radiation in any public experiment located here; its only
radiation connection is its isolation source, and an isolation source is not an exposure. Thorium
tolerance is chemical-toxicity tolerance to a metal salt: no dose, no dose rate, no dosimetry was
reported. Whether the self-dose from Th-232 at that concentration over a culture period is
negligible was **not computed** in this audit and the conclusion does not require it.

**Taxonomic mapping, recorded and not collapsed.** NCBI taxid 94625 is now served as
`Brucella intermedia`, lineage `... Brucellaceae; Brucella/Ochrobactrum group; Brucella`, with
*Ochrobactrum intermedium* retained as a synonym. The reclassification is actively contested in
the literature (`10.1128/jcm.00438-23`, *"If You're Not Confused, You're Not Paying Attention:
Ochrobactrum Is Not Brucella"*). This audit assigns no risk group — that determination belongs to
the institutional biosafety committee and must precede culturing — but it records that the
containment driver named in `data/calibration/baseline_condition.csv` may not be *C. neoformans*,
and that both names must travel together, as the *Wangiella*/*Exophiala* mapping does.

## 4. Radiotrophic versus radioresistant: the measured basis for each label

This is the spine of the document. Labels below are assigned on the endpoint actually measured,
never on the title, never on the organism's reputation, never on an adjective in a strain
description.

| Record | Species | Endpoint measured | Label earned | Label refused, and why |
|---|---|---|---|---|
| Dadachova 2007 | *C. neoformans* H99, *W. dermatitidis*, *C. sphaerospermum* | CFU ~2.5× at 18/23/30 h (P = 0.006); 14C-acetate ~3×; dry weight +6.5% (CN, P = 0.02) and +8-10% (WD, P = 0.02) | **radiotrophic** — growth and metabolism measurably increased under exposure | `melanized_radioprotective`: the contrasts are growth contrasts, not damage or survival contrasts |
| Zhdanova 2004 | soil micromycetes incl. *C. sphaerospermum* 3176 | mean return angle < 90° toward a 109Cd source | **radiotropic** — directional growth toward a source | everything else; no survival, no growth rate, no melanin |
| Vember 2015 | *C. sphaerospermum* 5-1, 60, 852, 3176 | survival as % of unirradiated control; D37 curves at 14 d | **radioresistant** — elevated survival, no growth-enhancement claim | `radiotrophic`: no growth endpoint of any kind |
| Cortesão 2020 | *A. niger* N402 | LD90 366 / 506 / 112 Gy | **radioresistant** | `radiotrophic`: no growth, biomass or metabolic measurement was made *under irradiation*. This is an absence of measurement, not a measured null |
| Qiu 2006 | *S. oneidensis* MR-1 | 170/87 genes over 1 h recovery; D20 = 40 Gy | **radiation_responsive** | `radioresistant` is **affirmatively disconfirmed**: *"about 20 times less than that of Escherichia coli and about 200 times less than that of Deinococcus radiodurans"* |
| Brown 2015 | *S. oneidensis* MR-1 | FT-IR and MALDI fingerprints; D10 ~84 Gy; Fe(III) reduction more than doubled at 50 Gy | **radiation_responsive** | `radiotrophic`: total biomass yield unchanged, lag *extended* to ~7.5 h, the effect required exogenous riboflavin as a shuttle, and viability was a small fraction of control |
| JRR 2024 | *A. niger* NBRC 6342 | germination 69% against colony formation 5.5% at 0.4 kGy; OD600 and elongation kinetics | **radiation_responsive** | `radioresistant`: survival fell monotonically; the framing is inactivation kinetics |
| Khajo 2011 | *C. neoformans* cap67 | EPR, 263 nm absorbance, TBARS, 213Bi binding at 120 and 360 Gy | **radiation_responsive** | `melanized_radioprotective`: the survival comparison is **not performed here** — verbatim, *"These radiation doses are in a range known from prior work to demonstrate significant protection"*, citing Dadachova 2008 |
| Malo 2018 | *C. neoformans* | cell-wall integrity retained in melanized cells at mid and max dose; non-melanized cells showed two failure morphologies | **melanized_radioprotective**, **radiation_responsive** | `radioresistant`: the survival clause refers to the group's prior work, not to a curve in this paper |
| Averesch 2022 | *C. sphaerospermum* ATCC 11289 | 147 vs 151 CPM, n = 1 dish, p = 0.069; flight-versus-ground growth | **none_demonstrated** | see below |
| Vasileiou 2020 | melanin in solution | Sr-90 beta transmission versus a cellulose control of matched elemental composition | **none_demonstrated** (a measured null) | `radiation_shielding`: the term requires transmitted radiation reduced *by biomass*, and there is no biomass and no melanin-specific reduction |
| Robertson 2012 | *W. dermatitidis* (= *E. dermatitidis*) wild type and non-melanized *wdpks1* | RNA-seq, >3000 genes differentially expressed; cell division and cell size; survivability after long-term low-dose exposure; ROS; carotenoid | **pending full text** — see note | no label assigned: the abstract's growth clause reads *"we confirmed that ionizing radiation enhanced cell growth"*, and **whether that growth was measured here or restated from prior work is not decidable from the abstract** — the Khajo 2011 distinction exactly |

**Robertson 2012, and why it is not the source for the isogenic nulls.** Read from PubMed on
2026-09-01 (PMID 23139812, `10.1371/journal.pone.0048674`). It reports that low-dose exposure
*"significantly increased survivability of both the wild-type and the wdpks1 mutant"* and that
ribosomal biogenesis genes were up-regulated in the irradiated wild type *"but not in the
irradiated wdpks1 mutant"* — a melanin-**dependent** difference. That is a positive result, not a
null, and the manuscript's isogenic-null sentence rests on the three studies in §9.2 instead.
**SCOPE: the abstract alone, via NCBI eutils, with the full text unretrieved — so the growth
clause's provenance stays open and this record earns no label.**


### The ISS experiment, in detail

`10.3389/fmicb.2022.877625` is the flagship record behind the word "radiotrophic" as it appears in
this repository, and it earns no phenotype at all.

- **No dose axis exists.** Verbatim: *"Due to the nature of the employed radiation sensors (PIN
  photodiode), dosimetric data was not obtained."* Two mutually inconsistent count-to-dose
  conversions appear in the supplement (≈2.2 nGy per count; ≈0.5 nGy per CPM), differing by 4.4×,
  and both divide someone else's whole-body dose rate by this detector's count rate.
- **The shielding result is not significant.** Phase 1 gave p = 0.970 (a proper negative control);
  phase 3 gave p = 0.069 with lower counts beneath the biomass. n = 1 dish, one sensor pair, no
  replication of the flight arm, no attenuation coefficient, no areal density. The attenuation path
  is dominated by 13.33 mm of potato dextrose agar.
- **The growth comparison is flight versus ground**, with microgravity, transit history, temperature
  regime and dish-to-dish variance confounded. The authors say so: *"the potential role of
  microgravity and its impact on fungal growth ... cannot be assessed and/or gauged with the
  employed experimental setup."*
- **Melanin was never measured.** The paper says "presumably constant" and "putatively".
- **Peer review removed the claim.** Every page header of the published supplement still reads
  *"Supplementary to: A Self-Replicating Radiation-Shield for Human Deep-Space Exploration:
  Radiotrophic Fungi can Attenuate Ionizing Radiation aboard the International Space Station"*,
  while the published article is titled *"Cultivation of the Dematiaceous Fungus..."*. The words
  "radiotrophic" and "radiation-shield" were struck between preprint and publication. This project's
  own rule — phenotype follows measurement, not titles — is visibly what peer review applied.
- **The widely-quoted attenuation figure exists in no version.** bioRxiv v1 says 2.17 ± 0.35%;
  bioRxiv v7 says approximately 0.84%; the published article says neither. The commonly cited
  "2.17 ± 0.25%" is a mis-transcription of a superseded preprint.
- **The headline growth number does not reproduce internally.** Abstract: 1.21 ± 0.37-fold.
  Supplement §D: 1.64-fold ± 0.279 SE on a logistic-slope basis. The article's own exponential rate
  constants give 1.24 / 1.29 / 1.53, mean 1.35 ± 0.16 (inference, arithmetic mine). The basis of
  1.21 ± 0.37 is undocumented. Recorded, not repaired, not averaged, not chosen between. Note also
  that 1.21 ± 0.37 does not exclude unity.

### Two unsourced repository claims

`biofilms_potts.jl:22` and `:24` describe *C. neoformans* and *C. sphaerospermum* as "radiotrophic"
in inline comments with no citation; `data/calibration/spatial/entity_semantics.csv` repeats
"Radiotrophic melanized yeast" for CN with `evidence_basis=declared`. This audit locates the primary
source those comments trace to — Dadachova 2007 — and simultaneously locates three same-species or
same-genus results pointing the other way. The comments should become cited-and-contested or be
downgraded; that is a documentation finding, reported and not fixed.

Two further in-repo claims are refuted rather than merely unsourced:

- `preprint/modeling_radiotrophic_fitness.md:31` asserts *A. niger* melanin "confers substantial
  radioprotection". Cortesão 2020 states the opposite for every ionizing modality it tested:
  *"deficiency in pigmentation did not alter survivability after exposure to ionizing radiation
  (X-rays or cosmic radiation)."*
- `preprint/modeling_radiotrophic_fitness.md:178` tabulates a melanin production rate
  α_M = 0.03-0.10 µg/cell/Gy for *A. niger*. No located measurement supports it — Cortesão 2020
  quantifies no melanin at all — and a per-cell basis is not one the CPM entity semantics admit,
  since a cell ID is a computational biomass parcel.

## 5. Raw-data availability

The distinction that matters is between a deposit with a manifest, a supplementary document, and a
rendered figure. An open-access article is not an available dataset.

### Volumetric imaging deposits (raw, manifested, licensed)

| Deposit | Species | Content | Size | Licence |
|---|---|---|---|---|
| S-BIAD474 | *B. subtilis* | Leica `.lif`, colony biofilm volumetric imaging; Fig 7 has six files, 100% GFP condition is one acquisition per timepoint | ~1 TB total; Fig 7 files 1.8-4.3 GB | none stated in the record |
| heiDATA KH6NDD | *B. subtilis* | 31 ImageJ TIFF tomograms, label-free hydrated SXT | 13,799,884,261 B | CC BY 4.0 |
| S-BSST261 | *B. subtilis* | 5 Leica `.lif`, edges at 24 h and 5 d | 18,529,647,313 B | none stated |
| Dryad f4qrfj71n | *B. subtilis* | 4 × `.7z` confocal of 20 µm cryosections + 5 zips | 49,717,304,106 B current version | CC0-1.0 |

**Zero volumetric deposits exist for *C. neoformans*, *C. sphaerospermum*, *A. niger*,
*D. radiodurans*, *S. oneidensis* or *O. intermedium* AM7.** For *C. neoformans* the emptiness was
reproduced independently on three repository APIs.

### Non-volumetric raw deposits

| Deposit | Content | Size | Note |
|---|---|---|---|
| Dryad 8sc52 / Zenodo 4944335 | SERS hyperspectral maps, three biofilm ages | 5,141,559,882 B over 7 files | the third array axis is **spectral**, not spatial: each ASCII export is exactly 1600 × 63,001, verified by line-length arithmetic on a 2 MB range read. Filenames carry the token `z0` |
| Zenodo 14994137 | Olympus `.oir` plus `.ims` and renders | 64,402,398 B | every `.oir` reports `<BaseDimension>1024 1024 1 2 1</BaseDimension>` and `<ImplExtendMax>212.132 212.132 0.000</ImplExtendMax>` — Z = 1, z extent 0.000 µm. Single plane |
| depositonce-18458 | 10 CSV/text files, chitin and glucan assays | 46,491 B | see §8 for the basis defect |
| Zenodo 3667494 | Sr-90 beta transmission spectra, Geant4 10.5 applications, notebooks | 13,235,258 B, MD5 `641b6fb1…` | the only measured transmitted-radiation data in the whole audit |
| Zenodo 21921798 | `NMR_DepositionV2.zip`, *Cryptococcus* cell walls | 143,676,256 B | contents uninspected; must not be described as raw until unpacked |
| GEO series | GSE80230 (25,835,520 B), GSE117227 (850.0 kB), GSE142318, GSE152116, GSE294405 | small | array and RNA-seq only |

### Supplementary documents that are not datasets

Frontiers `Data_Sheet_1` PDFs and DOCX files (190,728 B; 2,204,940 B; 13,296,418 B) are figure
supplements. `Table_1.XLSX` for the ISS study (13,826,234 B) is a two-dimensional optical-density
and count table. A supplementary XLSX is not volumetric data.

### Figures only

*C. neoformans* biofilm thickness (Ravi 2009), *D. radiodurans* biofilm CLSM (Shukla 2022),
*A. niger* biofilm COMSTAT biovolume, and every melanin-versus-modality survival curve exist only as
rendered images inside articles. A rendered confocal figure in a paper is not a raw 3D dataset.

### Two size traps, both live

Dryad `storageSize` is cumulative across versions. For `p8cz8w9rw` it reports 835,554 B while the
current version (4) is 768,305 B over 5 CSVs. For `f4qrfj71n` it reports 215,767,147,677 B against a
current-version total of 49,717,304,106 B — a factor of 4.34. For `8sc52`, `storageSize` reports
5,141,603,486 B against a manifest sum of 5,141,559,882 B, a 43,604 B discrepancy that is
unreconciled. The manifest sum is the download budget; `storageSize` is not.

`data/calibration/spatial/dataset_candidates.csv` row `shewanella_dryad_8sc52` currently carries
`size_bytes=5142000000` with `size_basis=current_version`. That figure is a rounded `storageSize`.
The same row's `rationale` cell terminates mid-sentence at *"Do not use it to select a."*

## 6. Radiation-modality compatibility and the OpenMC photon benchmark question

Modalities located in this audit, held distinct:

| Modality | Records | Photon-benchmark compatible? |
|---|---|---|
| 60Co gamma (1.173, 1.332 MeV) | GSE80230 (3 kGy), Qiu 2006 (40 Gy at 109 Gy/min), JRR 2024 (0.2-0.8 kGy at 1.87 kGy/h), Yuzon 2023 (1000 and 8000 Gy at 250 Gy/min) | in principle yes; in practice no dosimetry chain is stated in any of them |
| 137Cs (661.657 keV) | Khajo 2011 (11.94 Gy/min), Bland 2022 (~350 µCi, TLD-anchored) | the only photon-dominated fields with a stated source and a spatial variation |
| X-ray, polychromatic | Cortesão 2020 (200 kV, 15 mA, **no filter**), Brown 2015 (Faxitron CP-160, 160 kV, 6 mA, tungsten) | bremsstrahlung spectra with no stated effective energy, filtration or HVL |
| Heavy ions | Cortesão 2020 (150 MeV/n He, LET 2.2; 500 MeV/n Fe, LET 200) | no — charged particles |
| Protons, deuterons, alphas | GSE152116, Malo 2018, Szarka 2026 | no |
| Beta | Vember 2015 (90Sr/90Y), Dadachova 2007 cell arm (188Re), Vasileiou 2020 (90Sr) | no |
| Mixed alpha/beta/gamma, internal | Malo 2021 (225Ac at 1 mCi/mL in the medium) | no; the source is distributed in the growth medium |
| Mixed low-Earth-orbit | Averesch 2022, MARSBOx | no |
| UV-A, UV-B, UV-C, pulsed light, visible | many | **not ionizing radiation for this project**, under any circumstances |

### What a reconciliation would require, case by case

**The ISS field cannot be reconciled, and the attempt is a stop condition.** Four inputs are needed
and all four are absent: an incident particle spectrum by species and energy at the sensor
location; a detector response function for the X100-7 SMD in that field; a mass model of everything
between space and the diode; and the attenuator's own areal density and elemental composition. The
supplement names OLTARIS, SPENVIS, GEANT4 and CREME96 and then states *"such analyses are out of the
scope of this study."* No Monte Carlo of any kind was run. A photon-only model additionally cannot
represent a field the authors attribute to trapped-belt electrons of hundreds of keV and protons
above 100 MeV. Supplying any of the four inputs would mean inventing it.

**Bland 2022 is the closest thing to a usable photon geometry, and it is still not a benchmark.**
It has a 137Cs source of stated activity, TLD measurements at two positions, and a real spatial dose
variation: for *C. cladosporioides* the on-target dose ran from 85.2 rad at the central mycelial plug
to 6.1 rad at the growing edge — a 14-fold measured gradient across a single colony. What blocks it:
the gradient is **MicroShield point-kernel modelling spot-validated by two TLD readings**, one of
which disagreed with the model by a factor of two (~100 rad measured against ~50 rad predicted); the
total dose is 0.5-1 Gy over seven days, three to four orders of magnitude below every other ionizing
record here; and the biological observable is colony diameter, an areal growth measurement.

**Cortesão 2020's X-ray arm** would need the tube spectrum characterised. It was run with no
filter at 200 kV, so it is a broad bremsstrahlung distribution whose effective energy is not stated;
pooling it with the He and Fe arms — which the paper reports alongside it — would be a modality
collapse.

**Vasileiou 2020 is the only measured transmission dataset in the audit**, and it measures **Sr-90
beta (electron) attenuation**. The 40 kVp X-ray case in the same archive appears in Geant4
simulation only and was never measured. The comparison engine is Geant4, not OpenMC. Agreement
against it would license the statement "the transport-comparison workflow was exercised" and nothing
stronger.

**Its result is also the substantive one.** Melanin in solution and suspension provided no improved
shielding compared to a cellulose control of similar elemental composition. That is the physically
expected outcome if attenuation is governed by areal density and Z, and it strips melanin of any
special shielding status absent a spatial-arrangement argument. The paper's own positive finding is
that spatial arrangement — thin films, hollow microspheres — could matter in Monte Carlo, which is a
geometry claim, not a chemistry claim.

**No candidate in this audit can serve as an OpenMC photon-transport validation.** The two records
that come closest fail on opposite axes: Bland 2022 has a photon field and no transport measurement;
Vasileiou 2020 has a transport measurement and no photons.

## 7. Does any exact target species have a calibrated 3D radiation-exposed biomass field?

**No.**

Not one. The three requirements — an exact modelled species, a physically calibrated three-
dimensional biomass field, and a characterised radiation exposure — are never satisfied together,
and are satisfied pairwise only twice:

- **3D and exact species, no radiation:** *B. subtilis* (S-BIAD474, heiDATA, S-BSST261). Both
  volumetric deposits have caveats that matter independently of radiation. S-BIAD474's only channels
  are constitutive cytoplasmic GFP and mKate2 — no matrix stain, no membrane dye, no general biomass
  stain anywhere in the study — so its `segmentation_basis` is `cells_only`, which
  `HydratedVolume.require_whole_biofilm()` refuses outright. Its 100% GFP volumetric condition is a
  single acquisition at 48 h and a single acquisition at 72 h, so it does not satisfy the held-out
  rule on its own. heiDATA's deposited TIFF headers carry a self-inconsistent calibration
  (`unit=micron` with `spacing=2.105e-07`, which is 0.21 picometre), and the article's field of view
  is a per-acquisition range of 15-17 µm, so the voxel size cannot be recovered even by arithmetic.
  Eight of its 31 files are cells in suspension at t = 0, not biofilm.
- **Exact species and radiation, no 3D:** *C. neoformans*, *C. sphaerospermum*, *A. niger*,
  *S. oneidensis*. Every one of these is transcriptomes, survival curves, plate photography or
  spectroscopy.

`raw_3d_images` is `false` on every radiation-exposed record in this entire audit. There is no
confocal z-stack, no calibrated volumetric segmentation and no parcel-compatible dynamic observable
anywhere in the radiation literature for the seven species. A radiation-exposed 3D biomass field for
a target species must be **acquired**, not retrieved.

The consequence for the pitch branch is unchanged from the repository's existing position and now
rests on stronger evidence: `pitch_selection.select_pitch()` refuses without declared tolerances and
without a held-out stack, `scale_candidates.select()` refuses without declared acceptance
thresholds, and `spatial/report.evaluate()` refuses a surrogate. None of those refusals is relaxed
by anything found here, and the ~16-17 µm *C. neoformans* biofilm thickness from Ravi 2009 is a
constraint observation, not a pitch — reading it as one would be the N=60 error wearing a biological
hat.

## 8. Does any source pair morphology with wet/dry material measurements?

**No source in this audit pairs a calibrated hydrated volume with a wet mass, a dry mass and matched
blanks, for any organism, under one growth condition.**

Four near-misses, each instructive about a different failure mode:

**Frontiers Microbiol 13:975763 (*A. niger*, 2022) — the closest approach, and it fails on the
blank.** Verbatim: *"each colony-carrying filter, of each biological replicate (n = 3), was
recovered and placed inside a pre-weighted aluminum paper. The weight of the aluminum with a colony
was measured before (wet weight = ww) and after desiccation (dry weight = dw) for 24 h at 60 °C."*
The weighed object is aluminium foil plus filter plus colony. Deriving w_water needs the aluminium
tare **and** the filter's own retained water, not merely a blank — a stronger obstruction than a
single missing blank. No hydrated volume was measured; thickness alone is not a volume, and
`_require_volume()` refuses a bare number precisely because it cannot say what it encloses. No
surface-water removal protocol is documented. The raw data are author-request.

**Frontiers Microbiol 12 (*S. oneidensis* electrospun anodes) — a derived quantity wearing a
gravimetric label.** Verbatim: *"The dry weight equivalent of the biofilms attached to the anodes was
quantified by analysis of the protein content of cell lysate."* That is a colorimetric protein
reading multiplied by an undeclared protein-to-dry-mass factor. The mass is per whole electrode with
no area normalisation, on a porous nanofibre mat whose biofilm reportedly fills ~8% of pore space.
No CLSM, no thickness, no water content, no wet mass.

**AEM 2017 (*S. oneidensis* flow cells) — the best method precedent and zero calibratable numbers.**
It has CLSM z-stacks with x-z reslicing, thickness measured (~10 µm at +0.4 V), cells and EPS stained
**separately** (AFP for cells, ConA-Alexa 594 for EPS), and an areal cell density of 40 ± 2.2 cells
per 100 µm² at ~10 µm from the electrode. But protein and EPS are relative comparisons in figures,
not absolute masses with matched blanks; thickness at five random points is not a
`whole_biofilm_envelope` volume; there is no water content; and the raw stacks are not deposited.

**Manobala 2020 (*D. radiodurans*) — the only record using biofilm dry mass as an independent
variable in a metal-loading experiment**, and it fails on basis: a sorption capacity is metal mass
per **dry** biomass, with no water fraction to convert it, no closed elemental set, and no
metal-free control biofilm characterised the same way.

### The composition record, and its basis defect

The one machine-readable compositional deposit for a modelled species is `10.14279/depositonce-18458`
(*A. niger*, CC BY 4.0, 46,491 B, 11 chitin-synthase strains, n = 3 flasks each with three
sub-samples in photometric duplicate, a reagent blank present and subtracted). Control MJK17.25
reads 214.1198505 ± 3.389406288 µg glucosamine per mg, flasks 215.388 / 209.483 / 217.489; the
across-strain range runs 156.64 to 352.13 µg/mg.

**The two assays in that deposit are on two different mass bases under one misleading label.** The
companion methods state that chitin was assayed from *"4 mg (± 0.5 mg)"* of **cell wall extracted
from ~50 mg freeze-dried biomass** by two 1 M NaCl washes and a water wash, while glucan was assayed
from *"5 mg (± 0.5 mg)"* of **freeze-dried, ground whole biomass**. The deposit's own
`Sample weight [mg]` rows confirm it: 3.5-4.5 mg for chitin, 4.6-5.5 mg for glucan. So chitin is
g/g cell wall and glucan is g/g dry mycelium. The plate files label both `[mg/mg BM]`; the two
summary files name **no denominator at all**, reading only `Glucosamin [µg/mg]` and
`Glucan [µg/mg]`. A reader taking a summary CSV alone has no denominator to check, which is worse
than an ambiguous one. A third file, `Chitin_CFW_assay.csv`, is calcofluor-white fluorescence in
arbitrary units normalised to the control — a third basis again. Three files, three bases, one
deposit, and the README discloses none of it.

Converting cell-wall mass fractions to a dry-cell basis needs the wall-to-cell dry mass ratio, which
is not measured in the deposit, the thesis or the peer-reviewed companion. Converting a dry-cell
basis to hydrated bulk needs w_water, also absent. Neither conversion was performed, and no
literature-typical value was substituted.

### Refusals recorded in this section

- **The 1.09 g cm⁻³ buoyant density** (Loferer-Krossbacher/Fagerbakke lineage) with a 30% dry-matter
  figure was located and **refused**. It is a cell buoyant density from density-gradient
  centrifugation of individual cells. `density_basis` must be `wet_bulk`, and the binding sentence
  assigns a hydrated bulk biofilm material to a voxel — a parcel bundles cells *and* their
  extracellular material *and* the interstitial water. A cell buoyant density excludes the matrix
  and the pore water and overstates the bulk by roughly the void fraction. It is dimensionally
  plausible, widely cited, easy to justify in a sentence, and wrong. `density_g_cm3` stays blank.
- **The Frontiers Env Sci 2017 biofilm concentrations** (cells 224 mM, bound EPS 168 mM, loosely
  associated EPS 56 mM, at 113 g biomass per mole from the empirical formula C₅H₇O₂N) were located
  and **refused**. They are cited to Nielsen 1997 and Laspidou & Rittmann 2002 — imported from prior
  literature, not measured on that system. Under the repository's vocabulary they are
  `evidence_basis=declared`, which sits in `BLOCKED_EVIDENCE`, and the export gate would refuse them
  with *"a proxy or declared composition is a placeholder, and shipping it would launder it into a
  measurement."* Importing them would be exactly that laundering. The paper also measures within a
  2 × 2 × 2 mm³ NMR voxel — a millimetre-scale bulk-averaging volume, three orders of magnitude
  coarser than any plausible lattice pitch.
- **Genome-scale metabolic model biomass objective functions** were the only route located to an
  *S. oneidensis* elemental composition, and were **refused**: a biomass reaction is a declared
  stoichiometry assembled from macromolecular assumptions and frequently borrowed wholesale from
  *E. coli*, on a dry basis, failing both the evidence rule and the basis rule.

### The active-reducer fraction

`f_red,dry` — the fraction of biofilm dry mass that is *active metal-reducing* biomass — was not
found for *S. oneidensis*, for AM7, or for any biofilm system inspected. Nothing found even attempts
the measurement. Combined with Brown 2015 — Fe(III) reduction more than doubled in irradiated
biomass at 50 Gy, in a population where viability was a small fraction of the non-irradiated control
— the record supports a stronger statement than "not measured": **reducing activity is decoupled
from viable cell count**, and therefore a fortiori from taxonomic abundance and from CPM site
occupancy. That is direct empirical support for the existing refusal in `active_from_taxonomic()`,
and against the model's hard-coded default `X_red = 0.3` labelled "(Shewanella proxy)". Note also
that even a perfect `f_red,dry` measurement would not unblock the RADIODIALYSIS gate, which is
blocked by a units error at `biofilms_potts.jl:1341-1342` and `:1445-1446`, not by missing data.

## 9. Evidence contradicting simplified melanin assumptions

This section exists because the melanin assumption is the one this project is most likely to import
unexamined, and because the evidence against it is stronger, more recent and more replicated than
the evidence for it.

### 9.1 The word "melanin" denotes at least three different polymers

Across the records audited: **induced DOPA-melanin** (*C. neoformans*, grown on an L-dopa
precursor), **constitutive DHN-melanin** (*C. antarcticus*, *K. petricola*, and the pks1-dependent
pigment of *A. niger* and *E. dermatitidis*), and **allomelanin** (*E. viscosa*). Pacelli 2017 states
the distinction explicitly for its own two arms. These polymers differ in elemental composition. No
record in this audit supplies an elemental composition for any of them, so no melanized-biofilm
material composition is derivable from the literature at all. Any future composition must name which
polymer.

### 9.2 Isogenic melanin-null contrasts find no survival effect under ionizing radiation

Three independent studies on *Exophiala dermatitidis* wild type against its non-melanized *pks1*
mutant, across four modalities:

| Study | Modality | Result, verbatim |
|---|---|---|
| GSE142318 / Environ Microbiol 2020 | gamma | *"We find that melanin does not protect E. dermatitidis from gamma-radiation."* |
| GSE152116 / Front Microbiol 2020 | alpha, proton, deuteron | *"No significant difference in survival was observed between these strains under any condition, suggesting that melanin does not impart protection to acute irradiation to these particles."* Fitted α: WT 6.38 kGy⁻¹ (95% CI 0.26, 16.51) against *pks* 2.10 kGy⁻¹ (95% CI 0.01, 10.04) — the intervals overlap across almost their entire range |
| Microbiol Spectr 2023 | 60Co at 250 Gy/min | *"melanin did not directly affect survival after gamma radiation exposure"* |

The same *E. dermatitidis* gamma study identifies what **does** dominate: *"environmental factors
such as nutrient availability, culture age and culture density are much greater determinants of cell
survival after exposure"* than melanin. That is a standing confounder warning: a melanin effect
measured without controlling culture age and density is not interpretable.

### 9.3 The same finding holds inside the modelled species

Cortesão 2020, on *A. niger* — one of the seven — with an isogenic pigmentation mutant:

| Modality | N402 LD90 | MA93.1 (*ΔfwnA*) LD90 |
|---|---|---|
| X-ray | 366 Gy | 353 Gy |
| Helium ion | 506 Gy | **567 Gy** |
| Iron ion | 112 Gy | 112 Gy |

Indistinguishable, and the He-ion pair runs the wrong way. The authors state it: *"Interestingly,
deficiency in pigmentation did not alter survivability after exposure to ionizing radiation (X-rays
or cosmic radiation)."* What does carry resistance in that dataset is DNA repair — the NHEJ mutant
MA78.6 collapses to LD90 57 / 55 / 50 Gy across the three modalities, a roughly sixfold drop, while
the pigment mutant shows none (inference, labelled).

The UV-C arm of the same study is the one place a ~2× melanin effect appears (N402 1038 J/m²
against MA93.1 512 J/m²), and it is **confounded**: the NHEJ mutant, with no pigment defect, drops
to 580 J/m². UV-C sensitivity tracks DNA repair at least as much as pigment — and UV-C is not
ionizing radiation for this project regardless. The same paper's crystal-violet biofilm arm runs the
other way entirely: at 1000 J/m² UV-C the reduction was 81% for wild type and **97%** for the
pigment mutant.

### 9.4 Pigmentation moves in the wrong direction under gamma

Bland 2022 (*Sci Rep* 12) is the only pigmentation-versus-dose measurement located, and its sign is
the opposite of the model prior. Under UV, pigmentation **increased**; under gamma, pigmentation
**decreased**; both at P < 0.001, in both species. Growth rate did not change under gamma at all
(P = 0.255 for *P. variotii*, P = 0.615 for *C. cladosporioides*). The authors note the result
*"contrasts prior reports that suggest melanin provides enhanced protection against gamma
irradiation."*

Independently, and inside a modelled species: the *C. neoformans* gamma study behind GSE80230
reports that *"the expression patterns of melanin-producing genes LAC1 and LAC2 were gradually
decreased after high (3 kGy) or low (1 kGy) doses of gamma radiation exposure ... suggesting that
gamma radiation itself did not trigger melanin formation."* That is a same-species, same-modality,
two-point dose axis with a **negative** melanin-pathway result.

If the model carries any term of the form `melanin_production = f(ionizing dose)` with `f`
increasing, both records are evidence against it. They also show that any pigmentation-versus-dose
result imported from a UV study would carry the **wrong sign** applied to gamma.

### 9.5 Melanin's effect is species-specific even between closely related black fungi

Catanzaro 2024 (*Environ Microbiol Rep*, CC BY, n ≥ 9, the best-replicated record in the audit)
knocked out the **same pigment** with the **same gene class** in two rock-inhabiting black fungi:
*"DHN melanin provided UV-B protection in C. antarcticus, whereas the same pigment or even
carotenoids proved ineffective in K. petricola"*, concluding that *"although the adaptive trait of
DHN melanization commonly occurs across black fungi, it is not equally functional."* Carotenoids
protected neither species when melanin was absent.

The modality is UV-B, so nothing transfers to the ionizing case. The **structural** finding does: a
per-species melanin constant derived from a different species is not defensible even between two
black fungi. This is independent support for the standing prohibition in
`biofilm_material_calibration.md` §1 — *"Do not create seven materials because there are seven
species"* — and it shows that even the biology does not transfer, let alone the transport-relevant
composition. The same record also carries its own confounder: *"C. antarcticus could tolerate higher
UV-B doses but was sensitive to white light, whereas K. petricola showed the opposite trend."*

### 9.6 Melanin changes which mechanism dominates, not how much protection there is

The *E. dermatitidis* experimental-evolution study found two separable resistance mechanisms
evolving in the same organism depending only on melanin status: constitutive DNA-repair expression
in melanized evolved lines, redox-homeostasis upregulation in non-melanized ones. The one place
melanin measurably mattered was **mutation accumulation, not survival**: *"melanized lines that
evolved to exhibit higher ionizing radiation resistance experienced fewer mutations, whereas
similarly evolved, non-melanized lines accumulated more mutations."*

That is an effect on damage with no effect on survival. It must not be generalised into a survival
term, and it is direct evidence that a single scalar "melanin protects" coefficient is the wrong
functional form.

### 9.7 The counterweight, and why it does not settle the question

Pacelli 2017 is retained deliberately because it points the other way. Verbatim: *"Melanin afforded
protection against high-dose (1.5 kGy) deuterons for both CN and CA (p-values < 10⁻⁴). For X-rays
(0.3 kGy), melanin protected CA (p-values < 10⁻⁴) and probably CN."* It is the one record in this
audit where melanin measurably alters survival under an ionizing modality.

Read closely it carries three internal contradictions worth more than its headline:

1. **Modality dependence.** Protection is firm for deuterons in both species but for X-rays it is
   firm only in *C. antarcticus* and merely *"probably"* in *C. neoformans* — so melanin protection
   is **not established for sparsely-ionizing photons** in the one species here that is also one of
   the seven. A hedge in an abstract is not a measurement.
2. **Opposite-sign metabolic readouts.** *"Deuterons increased XTT activity in melanized strains of
   both species, while the activity in non-melanized cells remained stable or decreased. For ATP
   levels the reverse occurred: it decreased in melanized strains, but not in non-melanized ones."*
   Two metabolic proxies moving in opposite directions in the same cells means "melanin protects
   metabolism" is not a coherent single claim.
3. **Chemistry dependence.** Induced DOPA-melanin against constitutive DHN-melanin — "melanin" is
   not one variable across the two arms.

### 9.8 Melanin is not a superior shield on the one occasion transmission was measured

Vasileiou 2020, measured against a cellulose control of matched elemental composition: *"melanin does
not provide improved shielding in comparison to cellulose from beta-radiation."* Melanin's shielding
reputation rests, on the only measurement located, on nothing intrinsic to the chemistry.

### 9.9 Melanin carries a baseline cost in the founding radiotrophy experiment

In Dadachova 2007 itself, non-irradiated non-melanized (Lac⁻) cells grew better and incorporated more
14C-acetate than non-irradiated melanized cells. Melanin carried a cost in the same experiment, so
the melanized-versus-albino comparison under radiation is confounded by a difference already present
without radiation. That baseline asymmetry is not usually quoted alongside the headline.

### 9.10 The energy budget is formally contested

Walberg 2015 (thesis, UW-La Crosse, advisor Volk; `handle:1793/73406`, licence **not specified**)
concludes that radiosynthesis has not been demonstrated, on two grounds: *"the energy deposited in
their irradiated systems was too low to show the amount of growth seen, and the carbon sources
provided in the media were sufficient to support the observed growth."* That is a closed-budget
argument — the same discipline this repository applies to mass fractions, applied to energy — plus a
confounder argument. It contains no new measurement, so it is `evidence_basis=derived`, and it is
the negative side of an open dispute, not a verdict. "Not demonstrated" is not "disproved."

### 9.11 There is no isogenic melanin-null contrast for any of the seven

Every clean melanin knockout located is in *E. dermatitidis*, *K. petricola*, *C. antarcticus* or
*C. cladosporioides*. The one exception in a modelled species is *A. niger* MA93.1 (*ΔfwnA*), and
its result is the **null** in §9.3. The *C. neoformans* melanized/non-melanized comparison in Pacelli
2017 is paywalled and its X-ray result for that species is hedged.

### 9.12 The two halves of `melanin_production` are in different papers

`melanin_production` requires quantitative melanin or pigmentation against dose, dose rate, modality,
nutrient state and time, with a melanized/non-melanized contrast. Scored against that:

| Record | Melanin quantified? | Dose axis? | Verdict |
|---|---|---|---|
| Wang 1996 | yes — 15.4% of dry cell mass, ESR kinetics over 14 d, three precursors | **no radiation at all** | fails |
| Khajo 2011 | structure only (EPR, 263 nm, TBARS), not amount | two doses, one rate, one modality | fails |
| GSE80230 | LAC1/LAC2 **transcript**, not pigment mass | 1 and 3 kGy | fails, and points negative |
| Dadachova 2007 | not quantified against dose | 188Re beta | fails |
| Bland 2022 | mean grey value, an optical proxy on an image | 0.5-1 Gy, 3-4 orders below every other record | fails |
| Malo 2018 | not quantified | doses not stated numerically | fails |

**No study supplies both axes.** The halves sit in different papers, on different strains (ATCC 24067
serotype D against cap67/B3501 against H99 serotype A), under different melanization regimes (10 d at
an unverified temperature against 20 d at 23 °C in darkness with 150 rpm shaking). Stitching them
would be a fabrication, and it was not done.

### 9.13 The 15.4% figure is not a species constant

The one melanin mass fraction located for a modelled species is Wang 1996's *"melanin comprised
15.4% of the dry mass of the cell after 10 days of growth in medium containing 1.0 mM L-dopa"*. Four
reasons it is not a constant, all measured **inside the same paper**:

1. **Strain** — ATCC 24067, serotype D, against an up-to-eightfold spread across seven strains of
   the same species. The repository's *C. neoformans* is not strain-pinned, and H99, the strain in
   both GEO deposits, was not among those profiled.
2. **Precursor** — relative melanin content 0.05 (catechol), 0.45 (dopamine), 1.0 (L-dopa). A
   twentyfold collapse on another precursor.
3. **Duration and temperature** — melanization was slower at 37 °C throughout, converging only by
   day 14, so a day-10 value sits on a still-rising curve; the temperature of the quantified culture
   is itself unverified.
4. **Medium** — defined minimal medium with a supplied catechol precursor is a maximal-melanization
   laboratory condition, not the nutrient state of any biofilm.

And a fifth, on basis: 15.4% is lyophilized-ghost mass over lyophilized-whole-cell mass, a
**dry-cell** quantity. The contract requires `hydrated_bulk_biofilm` on a `wet_bulk` basis, and the
conversion needs a water mass fraction from the same system under a documented surface-water removal
protocol, which does not exist. Access to the underlying mass pair is `browser_manual`: PMC174092 is
a scanned page-image article and both PDF endpoints returned HTML bot-walls, so the 153 mg / 994.7 mg
pair and the CHN percentages are recorded as unverified.

### 9.14 Summary of the melanin position

The public record does not support a scalar melanin-radioprotection term under ionizing radiation.
It supports the opposite for survival, it is silent or negative on pigment production against gamma
dose, it shows the effect is species-specific between closely related fungi, it shows melanin
changes which repair mechanism dominates rather than adding protection, and it shows melanin has no
special attenuation properties against a composition-matched control. `normalization.melanin_scale`
cannot be sourced from the literature as it stands, and this is a **missing measurement**, which
under this project's discipline stays missing.

## 10. Rejected candidates and reasons

A search that records only its accepts is not a provenance trail. These rejections are the part of
the audit most likely to be re-derived by a later pass if it is not written down.

### Rejected on modality

- **MARSBOx stratospheric balloon flight** (`10.3389/fmicb.2021.601713`, *A. niger* N402). The total
  mission ionizing dose was **20.9 µGy** at 75.5 ± 13 µGy/day measured with an M-42 silicon diode —
  the authors equate it to roughly ten days of natural background. The UV component was 1148 kJ/m²
  over 280-400 nm. The observed 2-log survival reduction is a UV effect with an ionizing
  contribution roughly seven orders of magnitude below the doses in the same group's X-ray and ion
  study. Calling it a mixed-space-radiation validation would attribute a UV result to an ionizing
  field.
- **PNAS 2025 melanin-PLA biocomposites in LEO** (`10.1073/pnas.2427118122`). The measured outcome is
  polymer mass loss and surface wrinkling — degradation of a plastic — not transmitted radiation
  reduced by biomass. The material is a dry engineered composite, not a hydrated biofilm. And the LEO
  exposure is confounded with atomic-oxygen erosion and thermal cycling, which the abstract names as
  co-acting mechanisms.
- **UV records, en bloc.** *C. sphaerospermum* wavelength-dependent UV survival (PMID 41183924),
  accumulated-melanin UV tolerance (PMID 39287919), *A. niger* pyomelanin UV-C shielding, the
  Evolution Canyon UV-A adaptive melanin response (`10.1371/journal.pone.0002993`), *C. neoformans*
  biofilm resistance to heat, cold and UV light (`10.1128/AEM.02506-06`), *D. radiodurans* UV-C
  stress panel (`10.3390/ijms25010421`), pulsed-light pigmentation, and the Catanzaro UV-B work.
  **UV is not ionizing radiation for this project.** Retaining any of them as radiation evidence
  would let a UV protection result be read as radioresistance. The Evolution Canyon paper fails a
  second time on basis: melanin was quantified as absorbance at 420 nm in Soluene-350 *per 10⁶
  cells*, an optical readout on a cell-count basis with no mass denominator, and the authors state
  their isolates were identified on morphology alone and *"could be heterogeneous and perhaps consist
  of molecularly different species"*.
- **GEO GSE21484, light-regulated genes in *Cryptococcus***. Visible light, and the taxon is
  *C. deneoformans* JEC21, serotype D — a separate species from *C. neoformans* var. *grubii*.
- **Zenodo *Deinococcus* crystallography deposits** (`10.5281/zenodo.8302395` and siblings). X-ray
  **diffraction** images from purified protein crystals. The X-rays are a structural probe. Recorded
  because a keyword search for "Deinococcus" plus "X-ray" returns them.

### Rejected on species or strain

- **Every *O. intermedium* strain that is not AM7** — SDCr-5 and BCR400 (Cr(VI) reduction), BB12 (Cd
  tolerance), CN3, MZV101, D-2, and the 96 genome assemblies for taxid 94625 whose reference is
  NCTC12171. These are real, quantified metal phenotypes in the same species, from different
  environments with different EPS. AM7's distinguishing claim is exceptional EPS production and
  thorium tolerance; nothing licenses transferring a chromate-reductase phenotype from a tannery
  isolate to it. *"Do not collapse related species"* applies a fortiori within a species.
- **Providencia thoriotolerans AM3** (`10.1038/s41598-021-82863-4`). Same laboratory, same isolation
  campaign, wrong genus. Read only for what it reveals about the AM7 methods: the group's panel is
  EPS yield in g% w/v, viscosity, flocculating activity, carbohydrate and protein in g% w/v, FTIR,
  and ICP-OES **on soil**. No dry weight, no XPS, no TEM-EDX, no thorium binding capacity in mg/g,
  and — like AM7 — a single 16S accession as its entire deposit.
- **Cladosporium cladosporioides** (Bland 2022) is not *C. sphaerospermum*. Different species in the
  same genus. Its strains are *"environmental isolates from facilities at Sandia National
  Laboratories"* with no culture-collection accession, so they are neither obtainable nor
  re-identifiable, and the paper does not state its biological replicate count. Its result is
  directionally load-bearing (§9.4) and quantitatively unusable.
- **Aspergillus nidulans** proton irradiation (GSE294405). Same genus as *A. niger*, not the same
  species, not melanized. Retained in the audit only for its structural finding: MSB oxidative
  pretreatment **before** irradiation increased tolerance while irradiation **before** MSB decreased
  it. The same two stressors in opposite order give opposite outcomes, which no dose-only scalar term
  can represent. It reports fluence and accumulated charge (1/2/5 × 10¹⁰ p cm⁻², 3200/6400/16000 ±
  20% pC) and **no absorbed dose in Gy anywhere**, so it cannot be placed on the axis the other
  records use.
- **Cryomyces antarcticus, Knufia petricola, Exophiala dermatitidis** — outside the seven. Retained
  as mechanistic surrogates for §9 only. Surrogate status is enforced, not conventional.

### Rejected on data availability

- **Bland 2022's underlying data.** *"The datasets used and/or analyzed during the current study
  available from the corresponding author on reasonable request."* No repository, no accession, no
  data DOI. An open-access article is not an available dataset.
- **Every *A. niger* biofilm CLSM/COMSTAT paper.** The right modality — confocal z-stacks with
  biovolume quantification — and no raw stacks deposited anywhere findable. What is published is
  rendered figures and derived biovolume-per-area scalars, which cannot be re-segmented, carry no
  voxel calibration, and cannot select a pitch.
- **ResearchGate figure pages** for *C. neoformans* B3501 biofilm confocal images. Extracted figure
  images on an aggregator with no manifest, no voxel calibration, no licence and no persistent
  identifier — the single thing most likely to be mistaken for the missing volumetric data.
- **Schultzhaus 2019's sequence data.** Three targeted NCBI queries located nothing. The same group's
  *Exophiala* data **are** deposited across five BioProjects, which makes the absence conspicuous
  rather than expected.
- **Cao 2011 EPS/U(VI) immobilisation** (`10.1021/es200095j`). Wrong species (*Shewanella* sp.
  HRCR-1, an environmental isolate), and the ACS figshare collection holds exactly one item — a
  supporting-information PDF. The PNNL Datahub node and OSTI record are bibliographic only.

### Rejected on representational grounds

- **NASA OSDR OSD-132, OSD-260, OSD-350** and the 91 *Shewanella* OSDR records. Genomes, proteomes
  and transcriptomes. A transcriptome alone does not define a CPM scale parameter and a proteome
  alone does not define an OpenMC material. The ISS *A. niger* isolate's exposure history is
  uncontrolled and uncharacterised, so even the qualitative radiation framing is not a dose-response.
- **S-BSST676** (*B. subtilis*, ~48 GB across three declared biological replicates at 24/48/72 h).
  The record's own file descriptions read *"Colony biofilm images obtained using stereo
  microscope"* — macro-scale photographs, not z-stacks. Painful, because it is the only exact-species
  deposit in the audit with explicitly labelled first/second/third biological replicates. Retain the
  citation as a precedent for how replicate structure should be declared.
- **S-BSST577**, four TIFF movies at exactly 377,886,841 B each — single-plane 2× macro time-lapse.
- **D. radiodurans metal-induced biofilm EM** (`10.1007/s00792-025-01403-4`) and *C. neoformans*
  ultrastructural sectioning. Two-dimensional ultrastructure sits orders of magnitude below the scale
  at which a parcel pitch is chosen.
- **EMPIAR-12837**, cryo-EM movies of *D. radiodurans* BamA. Exact species in name only; a purified
  membrane protein.
- **Quantitative modelling of microbial population responses to chronic irradiation** (PMID 26808049).
  A modelling synthesis over other people's data. Under the evidence hierarchy it is at best
  `derived`, which sits in `BLOCKED_EVIDENCE`.

### Rejected as superseded or non-citable

- **bioRxiv `10.1101/2020.07.16.205534`** (the ISS preprint). Superseded by the peer-reviewed article
  whose title deliberately removed both "radiotrophic" and the shielding claim, and whose reported
  result is 147 against 151 CPM rather than any percentage. Quoting the preprint's number would cite
  a withdrawn framing over the published one.
- **`github.com/chkern/space-radiation`.** The analysis code for the ISS study, and genuinely useful:
  `src/01_preprocess.R` reveals the supplement's sheet and column structure without downloading it,
  including the `_n` / `_nt` **noise channels** that must not be silently merged with the signal
  channels. But the GitHub API licence field is `null` and there is no LICENSE file — default all
  rights reserved — and its README points at bioRxiv v5 rather than the published supplement, so it
  does not pin the version of the data it was run against. Not a citable persistent artefact under
  the source-registry conventions. If ever cited, pin a commit SHA.
- **Tertiary sources** — Wikipedia "Radiotrophic fungus" and "Radiosynthesis (metabolism)",
  Grokipedia, microbewiki, a Stanford course page, several blogs, and a US patent on oral melanin
  administration. Read only to locate primary sources, which they did: they are how the Walberg
  thesis surfaced. Recorded so that a search-summary paraphrase from one of them is never mistaken
  for a measurement.
- **Zenodo `10.5281/zenodo.21329147`**, a 3.6 kB self-deposited HTML document proposing a melanin
  quantum-spin mechanism, and `10.5281/zenodo.15275837`, a self-deposited essay licensed GPL-3.0 on a
  PDF. No data, no method, no peer review. A resolving DOI is not evidence.

## 11. Author-contact opportunities

Named by the specific measurement wanted, not by paper. Ordered by what each would unblock.

**1. Zhdanova, Tugay, Dighton, Zheltonozhsky, McDermott — *Mycological Research* 108:1089-1096
(2004), via library or corresponding author.**
Wanted: whether the dose gradient **at the growth front** was measured or calculated, and by what
method; the source geometry and activity of the 109Cd emitter; and the return-angle data for
*C. sphaerospermum* 3176. This is the only architecture located anywhere that could in principle
identify a directional term: a directional source, a spatial gradient, a quantified directional
response, and a confirmed test on an exact modelled species. Three of four earlier blockers have
been resolved — the species membership of *C. sphaerospermum* is confirmed via secondary summaries
and corroborated by strain 3176's provenance in the Vember 2015 primary text. **One blocker remains,
and it is the dose-gradient characterisation.** Do not scrape; the DOI resolves to Elsevier
ScienceDirect (not Cambridge Core, as is sometimes stated).

**2. Cortesão / Moeller, DLR — Frontiers in Microbiology 13:975763 (2022).**
Wanted: the per-replicate wet and dry masses of the colony-carrying filters, **the aluminium tare
weight**, and a matched abiotic filter blank wet and dry under the same 24 h / 60 °C protocol.
This is the closest any public source comes to a water mass fraction for a modelled species, and
w_water needs no volume — only masses and a blank. It would still not yield a density, because no
hydrated volume exists. Data availability is *"will be made available by the authors, without undue
reservation."*

**3. Bachand / Bland, Sandia National Laboratories — Scientific Reports 12 (2022).**
Wanted: the biological replicate count per treatment group; the per-ROI grey-value data behind the
pigmentation result; and the full TLD map behind the 85.2 → 6.1 rad intra-colony gradient. The
replicate count is a hard blocker — a P-value without an n cannot be propagated — and the TLD map is
the only measured photon dose gradient located anywhere. Also worth asking whether the environmental
isolates were deposited anywhere, since without an accession the strains are unre-identifiable.

**4. Schultzhaus / Wang, U.S. Naval Research Laboratory — Environmental Microbiology 21:2613-2628
(2019).**
Wanted: the accession for the *C. neoformans* melanin-plus-gamma RNA-seq data, and the dose, dose
rate, dosimetry method and strain from the Methods. Three targeted NCBI queries found nothing while
the same group's *Exophiala* data are deposited across five BioProjects. This is the single most
important record for model revision in the audit and only its abstract is readable.

**5. Malo / Dadachova — *Fungal Biology* 124:368-375 (2020) and *Comput Struct Biotechnol J* 19
(2021).**
For the 2020 paper, wanted: the dose, dose rate, source geometry, replicate count and the growth
metric behind *"Gamma radiation caused increased colony growth irrespective of exposure history in
melanized fungus. Beta particles produced growth inhibition."* That sentence is the sign-flip
record — same organism, same melanin, three modalities, three different signs — and every number
behind it is paywalled. For the 2021 paper, wanted: resolution of the internal inconsistency between
*"183 Gy/5 weeks"* and *"a low dose rate of 0.02 mGy/min"*, which disagree by a factor of ~180
(183 Gy / 50,400 min = 3.63 mGy/min). Both numbers are recorded as printed; neither is usable until
the authors resolve it. Also wanted: where the study's sequence reads are, since PRJNA224192 is
correctly cited as the *E. dermatitidis* NIH/8656 **reference transcriptome** and the study deposits
nothing of its own.

**6. Saraf, Gujarat University — *J Hazard Mater* 388:122047 (2020).**
Wanted: a culture-collection accession or a strain transfer for *O. intermedium* AM7; the EPS yield
with its basis stated; the RSM run table with replicate structure; and any Th binding capacity in
mg per g. Whether the strain is obtainable at all is unverified, and an unobtainable strain is a
campaign-design problem, not a data problem. Note the biosafety item in §3 before requesting a
transfer.

**7. Ankur / Wang — Zenodo `10.5281/zenodo.21921798` (deposited 2026-08-13).**
Wanted: the strain designations, the culture conditions, and whether any melanized sample is in
`NMR_DepositionV2.zip`. The Zenodo description field is empty and the record carries no related
preprint identifier. Cheaper than an author request: unpack the 143.7 MB archive.

**8. Manobala / Rao — *Chemosphere* 269:128722 (2020).**
Wanted: the strain (the abstract names none; DR1-bf⁺ is an inference from their 2019 companion), the
biofilm dry-mass values with the blank protocol, and whether a uranium-free control biofilm was
characterised with the same measurements.

**9. Averesch — Frontiers in Microbiology 13:877625 (2022).**
Wanted: the basis of the 1.21 ± 0.37-fold headline, which reproduces from neither the article's own
rate constants nor the supplement's logistic-slope figure. Low priority: nothing in this repository
turns on it, and the record's value is as a corrective prior.

**10. A retrieval, not a contact.** `docs/calibration/dataset_screening.md` is cited by `source_id`
`DATASET_SCREEN_2026` in `data/calibration/spatial/sources.csv` and is absent from disk;
`dataset_candidates.csv` has no test asserting its `source_id` resolves against
`spatial/sources.csv`. Any ledger built from this audit inherits both gaps.

## 12. Verdicts

### Exact radiotrophic data

```
TARGETED_LAB_EXPERIMENT_REQUIRED
```

Radiotrophy — ionizing radiation measurably enhancing growth or metabolism — is claimed for two of
the seven modelled species, and the claim rests entirely on one 2007 paper whose cell-arm exposure
was a 188Re beta field characterised by analogy to Chernobyl rather than by a stated dosimeter, whose
growth effects are single-digit percentages (6.5% and 8-10% dry weight), whose integrated cell dose
is not stated (0.05 mGy/h over 15-30 h implies roughly 0.75-1.5 µGy — inference, labelled — some
eight orders below its own isolated-melanin arm), and in which melanin already carried a baseline
cost before any radiation was applied. Against it stand a same-species later independent report of no
substantial melanin protection under gamma, a same-species negative melanin-pathway result at 1 and
3 kGy, a same-genus null growth result under a TLD-anchored 137Cs field, an isogenic null in
*A. niger* across three ionizing modalities, and a formal energy-budget critique holding that the
deposited energy is insufficient and the media carbon sufficient. Peer review struck the word from
the flagship ISS title. No author contact resolves this: the discriminating measurement — dose-resolved
growth and metabolic rate under a characterised photon field, with an isogenic melanin-null contrast,
on the target strains, with the media carbon budget closed — does not exist in any archive and cannot
be requested from anyone, because nobody has made it. The one retrieval that would materially change
the picture, Zhdanova 2004, bears on *radiotropism*, not radiotrophy, and its own remaining blocker
is whether its dose gradient was measured or calculated. Until such an experiment is run, the
repository should carry `radiotrophic` for **no** species, and the unsourced inline comments at
`biofilms_potts.jl:22` and `:24` should become cited-and-contested or be removed.

### CPM calibration

```
PUBLIC_DATA_CAN_VALIDATE_PIPELINE_ONLY
```

Public exact-species volumetric imaging exists — for exactly one of the seven. Four *B. subtilis*
deposits are open, licensed and manifested, and the two volumetric ones can exercise the field
harness, the occupancy mapping and the pitch-evaluation machinery on a named modelled species rather
than on the *V. cholerae* surrogate the repository currently holds. That is a real gain and it is
the ceiling. S-BIAD474's only channels are cytoplasmic fluorescent proteins, so its
`segmentation_basis` is `cells_only`, which `HydratedVolume.require_whole_biofilm()` refuses; its
100% GFP volumetric condition is one acquisition per timepoint, so it does not satisfy the held-out
rule on its own; and its z-step of 3.87 µm against a lateral 0.869 µm is a 4.45:1 anisotropy that
floors any pitch derived from it. heiDATA's SXT tomograms are label-free and hydrated — the one
dataset in the audit that measures a biomass volume fraction the way the occupancy mapping needs —
but their deposited voxel calibration is self-inconsistent and the article's field of view is a
per-acquisition range, so no length derived from them is usable until the authors confirm it, and
eight of the 31 files are suspension cells rather than biofilm. Nothing here calibrates a radiation
term: `hamiltonian_radiation` requires spatial redistribution in a characterised gradient, and every
irradiation located is a uniform field, an areal growth assay or a survival curve. The one measured
photon dose gradient in the audit (Bland 2022, 85.2 → 6.1 rad across a colony) is applied to colony
diameter, an areal observable, and an areal experiment cannot identify a directional term. Two
model-revision findings are recorded separately and would change the shape of this verdict rather
than its token: hyphal-elongation tropism is measured on organisms whose CPM representation is an
isotropic biomass parcel with no elongation or connectivity term, and `growth_survival` remains
`unsupported_by_current_model` against genuinely excellent survival data.

### Materialization

```
TARGET_WET_BULK_MEASUREMENT_REQUIRED
```

No public source located in this audit pairs a wet mass, a dry mass and matched blanks with a
calibrated hydrated volume, for any organism, under one growth condition — and without a hydrated
volume there is no ρ_wet, no X_total, and no voxel binding, whatever the composition. The near-misses
each fail differently and instructively: the *A. niger* clinostat study weighed aluminium plus filter
plus colony with no tare and no filter blank, and measured a thickness rather than a volume; the
electrospun-anode study reported a "dry weight equivalent" derived from a colorimetric protein
reading per whole electrode; the *S. oneidensis* flow-cell study has the right staining design and
reports only relative comparisons; the *D. radiodurans* uranium work uses biofilm dry mass as an
independent variable on a dry basis with no water fraction. Composition is also required and is
separately unmet: the one machine-readable compositional deposit for a modelled species reports
chitin per gram of **extracted cell wall** and glucan per gram of **whole dry mycelium** under one
shared label, with the summary files naming no denominator at all, and neither the wall-to-cell dry
mass ratio nor w_water is measured anywhere in that deposit, its thesis or its companion paper. Cell
wall is not dry cell is not hydrated bulk. The ceiling of public data is component priors on
inconsistent bases, from which `mixture.blend()` would refuse a non-closing set and
`export.problems_for()` would refuse a declared or derived evidence basis. Three attractive numbers
were located and refused: a 1.09 g cm⁻³ cell buoyant density on the wrong basis, a set of biofilm
biomass concentrations imported from prior literature into a modelling paper, and genome-scale
biomass objective functions. `density_g_cm3` stays blank, and blank means not measured.

### Integrated

```
PUBLIC_COMPONENTS_FOUND_BUT_NOT_PAIRED
```

The components exist, separately, and none of them is paired with another. Exact-species radiation
response exists for four of the seven (*C. neoformans*, *C. sphaerospermum*, *A. niger*,
*S. oneidensis*), with real dosimetry in three of those cases. Exact-species 3D morphology exists for
one (*B. subtilis*) and for none of the four with radiation data. Exact-species composition exists
for one (*A. niger*, cell-wall basis, no radiation, no morphology, no biofilm — the culture was
submerged shake-flask mycelium). Two of the seven, *D. radiodurans* and *B. subtilis*, have no
ionizing-radiation record in this audit at all, and one, *O. intermedium* AM7, has a total public
footprint of 1086 base pairs and one paywalled abstract, having never been irradiated in any public
experiment. Nothing anywhere describes the seven-species consortium under its declared conditions,
and `clears_target_gate` can never be true for a single-species dataset however good its data are.
This reproduces and extends the screening verdict already recorded in
`data/calibration/spatial/dataset_candidates.csv` — the earlier finding is confirmed, now on
API-level verification across eleven repositories rather than on web search, and confirmed at the
level of individual named species rather than only at the level of the consortium. The paired
experiment must be designed: on one system under one growth condition, replicate-level wet mass and
dry mass with matched hydrated and dry substrate blanks, a hydrated volume on a
`whole_biofilm_envelope` basis or a declared `pore_volume_fraction`, ash content, CHNS plus
ICP-OES/ICP-MS closing to 1 on a wet-bulk basis, a calibrated 3D stack of the same specimen, and a
characterised dose field with a measured gradient. Until then the honest state of the material table
is header-only, and header-only is not a placeholder for invented rows.

---

**Nothing in this document populates a calibration value. Every number is attributed to the record
that measured it. Every quantity that could not be verified is recorded as unverified. Nine access
controls were encountered and none was circumvented.**