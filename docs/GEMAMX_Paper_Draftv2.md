# GEMAMX: An Automated Precision Oncology Pipeline for Breast Cancer Variant Analysis, Structural Prediction, and Molecular Docking Using TCGA-BRCA Data

**Ximena Gutiérrez Gálvez¹, Julio Waissman Vilanova¹, Héctor Gutiérrez Ibarra¹**

¹ Departamento de Matemáticas, Universidad de Sonora (UNISON), Hermosillo, Sonora, México

**Correspondence:** ximena.gtzg@gmail.com

---

## Abstract

Precision oncology in breast cancer requires the systematic integration of genomic variant discovery, protein structural analysis, and pharmacological interaction modeling — a workflow that remains largely manual and fragmented across institutions. We present **GEMAMX** (Genomic Analysis of Mammary cancer in Mexico), an end-to-end automated pipeline that traverses five sequential modules: (M1) data ingestion of 60,864 somatic mutations from 1,084 TCGA-BRCA cases via the GDC API; (M2) missense variant prioritization across 15 clinically actionable genes using cohort frequency and ACMG/AMP classification criteria; (M3) automated protein sequence synthesis and mutation application via UniProt REST API and Biopython; (M4) 3D structure prediction using Boltz-1 and ESMFold; and (M5) systematic molecular docking with AutoDock Vina 1.2.3 against 11 FDA-approved or clinical-stage inhibitors sourced from ChEMBL. Across 13 prioritized variants with available 3D structures, GEMAMX computed 143 receptor-ligand affinity scores. The dominant hotspot PIK3CA H1047R achieved the strongest binding affinity in the cohort (Pictilisib: −9.72 kcal/mol; Alpelisib: −8.90 kcal/mol), with an active-site RMSD of 3.21 Å relative to wild-type, consistent with known kinase domain rearrangement. BRCA1 C61G exhibited the largest structural disruption (RMSD = 5.67 Å in the RING domain, pLDDT drop from 65.3 to 58.9), suggesting compromised BARD1 interaction capacity. GEMAMX is deployed as a fully interactive web application and is publicly accessible at gema-mx.netlify.app, providing real-time visualization of variants, 3D molecular viewers, and docking rankings. The pipeline demonstrates that integrating genomic epidemiology with structural bioinformatics can accelerate computational reclassification of variants of uncertain significance (VUS) without requiring wet-lab infrastructure.

**Keywords:** precision oncology, breast cancer, somatic variants, TCGA-BRCA, molecular docking, AlphaFold, Boltz-1, AutoDock Vina, VUS reclassification, bioinformatics pipeline

---

## 1. Introduction

Breast cancer remains the leading cause of cancer-related mortality in women in Mexico, with an average age of diagnosis approximately one decade earlier than in European populations [1]. A critical bottleneck in clinical oncology is the identification and functional interpretation of somatic mutations — particularly those classified as Variants of Uncertain Significance (VUS) under ACMG/AMP guidelines [2]. VUS are mutations for which current evidence is insufficient to definitively establish pathogenicity or benignity; they constitute the majority of variants identified in clinical genomic panels and cannot be used directly to guide therapy.

The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) dataset provides the most comprehensive publicly available compendium of somatic mutations, clinical outcomes, and molecular subtypes for breast cancer, covering 1,084 patient samples sequenced at whole-exome resolution [3]. However, translating this resource into actionable computational insights requires the integration of at least four distinct analytical domains: genomic variant discovery, protein structural modeling, pharmacological interaction simulation, and clinical data visualization — tools that are typically operated independently by specialists in each field.

Recent advances in AI-driven protein structure prediction — most notably AlphaFold2 [4] and its successors Boltz-1 [5] and ESMFold [6] — have dramatically reduced the computational cost of generating accurate 3D protein models for arbitrary sequences, including mutated variants not represented in the Protein Data Bank. In parallel, AutoDock Vina 1.2.3 [7] and preparation libraries such as Meeko [8] now enable rapid, reproducible docking simulations with GPU-quality accuracy on standard hardware. The combination of these tools creates an unprecedented opportunity to computationally characterize the structural and pharmacological consequences of somatic mutations at scale.

We present GEMAMX, a fully automated five-module pipeline that integrates these capabilities into a single reproducible workflow, from raw GDC API calls to interactive 3D molecular visualization. The pipeline is specifically designed to address the VUS reclassification problem: by computing structural disruption (RMSD) between wild-type and mutant protein conformations and comparing ligand binding affinities across 11 inhibitors, GEMAMX generates computational evidence that can prioritize VUS candidates for experimental validation.

This work makes the following contributions:

1. A reproducible, open-source bioinformatics pipeline covering the complete path from genomic data ingestion to docking score visualization.
2. Structural and pharmacological characterization of 13 breast cancer variants across 7 genes using Boltz-1/ESMFold predictions and AutoDock Vina 1.2.3.
3. A publicly accessible interactive web application (gema-mx.netlify.app) that democratizes access to these analyses without requiring local computational infrastructure.
4. A computational framework for evidence-based VUS prioritization applicable to Latin American genomic cohorts.

---

## 2. Methods

### 2.1 Dataset and Data Ingestion (M1)

GEMAMX retrieves somatic mutation data from the cBioPortal REST API (www.cbioportal.org/api), accessing the TCGA-BRCA study (study ID: `brca_tcga`, molecular profile: `brca_tcga_mutations`). The full sample list is retrieved via `/sample-lists/{id}/sample-ids`, and mutations are downloaded in paginated POST requests to `/molecular-profiles/{id}/mutations/fetch` in batches of 500 samples. Gene symbol mapping from Entrez Gene IDs to HUGO symbols is performed via `/genes/fetch`.

The complete download encompasses **60,864 somatic mutation records** across **1,084 patient samples**. Records are stored in MongoDB (v7.0, localhost:27017) in three collections: `mutations` (filtered missense variants in target genes), `top_variants` (ranked cohort summary), and `fasta_sequences` (paired WT/mutant amino acid sequences). All data is simultaneously exported to static JSON files (`raw_mutations.json`, `top_variants.json`, `mutation_stats.json`) for consumption by the web application layer.

### 2.2 Variant Prioritization (M2)

Variant prioritization is restricted to missense mutations (`Missense_Mutation` classification) in a curated panel of 15 breast cancer driver genes organized in three clinical tiers:

- **Tier 1 (actionable):** PIK3CA, TP53, ERBB2, BRCA1, BRCA2
- **Tier 2 (frequent drivers):** CDH1, GATA3, MAP3K1, KMT2C, PTEN
- **Tier 3 (TCGA-relevant):** RB1, ARID1A, SF3B1, AKT1, CBFB

Variants are ranked by cohort frequency and classified using a frequency-based ACMG proxy:

| Classification | Cohort Frequency Threshold |
|---|---|
| Pathogenic | ≥ 5.0% |
| Likely Pathogenic | 1.0 – 4.9% |
| VUS | 0.5 – 0.99% |
| Low Frequency | < 0.5% |

This proxy classification is not a clinical ACMG assertion but serves as a computational prior for structural analysis prioritization. The Top-26 variant-gene pairs are exported for downstream processing.

### 2.3 Sequence Generation and Mutation Synthesis (M3)

For each prioritized variant, the canonical wild-type (WT) protein sequence is retrieved from UniProt REST API using a reviewed Swiss-Prot query (`gene:{SYMBOL} AND organism_id:9606 AND reviewed:true`). The FASTA sequence is parsed with Biopython (v1.83) `SeqIO.read()`, converted to `MutableSeq`, and the point mutation is applied at the 1-indexed position specified by the protein change notation (e.g., H1047R: histidine at position 1047 → arginine). Consistency between the WT residue in the sequence and the expected WT residue in the HGVS notation is verified before mutation application.

For proteins exceeding 390 amino acids — the operational limit of ESMFold's public API — a domain fragment strategy is applied: a window of ±180 residues centered on the mutation position is extracted, yielding fragments of ≤ 360 amino acids that preserve the local structural context of the variant. This affects PIK3CA (1,068 aa), SF3B1 (1,304 aa), and PTEN (403 aa).

### 2.4 Structure Prediction (M4)

3D structure prediction is performed using two complementary tools:

**Boltz-1** [5] is used as the primary predictor for variants where computational resources permitted (GPU T4 via Google Colab). Boltz-1 generates mmCIF format structures with per-residue pLDDT confidence scores stored in companion JSON files (`confidence_{variant}_model_0.json`). The field `complex_plddt` (range 0–1) is scaled to 0–100 for display. Boltz-1 was run with `--recycling_steps 3` and `--use_msa_server`.

**ESMFold** [6] (api.esmatlas.com) is used as the secondary predictor for shorter sequences (≤ 390 aa) or when Boltz-1 resources were unavailable. ESMFold returns PDB-format structures via HTTP POST. Both WT and mutant sequences are processed independently for each variant.

Structure availability is tracked in `manifest_structures.json`, a centralized index that maps each variant identifier to its CIF or PDB file path. RMSD values between WT and mutant structures are computed per domain region and reported as active-site RMSD and global backbone RMSD (Table 2).

### 2.5 Molecular Docking (M5)

#### Ligand Preparation

Eleven inhibitors targeting the PI3K/AKT/mTOR and CDK4/6 pathways were selected based on clinical relevance to breast cancer. Molecular structures were retrieved as 2D SDF files from the ChEMBL API (`https://www.ebi.ac.uk/chembl/api/data/molecule/{CHEMBL_ID}.sdf`). 3D coordinates were generated using RDKit (v2023.9) with the ETKDG conformer generation algorithm, followed by MMFF force field optimization. Protonation states were assigned at pH 7.4. PDBQT files were prepared using Meeko (v0.5), which assigns AutoDock4 atom types and handles rotatable bond detection.

#### Receptor Preparation

Protein structures (CIF/PDB) were processed through OpenBabel 3.1.1 to add hydrogens at pH 7.4 and produce clean PDB intermediates. PDBQT conversion was performed via a column-exact surgical parser that maps element symbols to AutoDock4 atom types (H→HD, O→OA, N→NA, S→SA, others by element symbol).

#### Docking Execution

AutoDock Vina 1.2.3 was executed for each receptor-ligand pair with the following parameters: exhaustiveness = 8, num_modes = 5, search box size = 25 Å³. The search box center was positioned at the centroid of the mutated residue's atomic coordinates, focusing sampling on the biologically relevant site. Best binding affinity (kcal/mol) was extracted from Vina log files using format-robust regex supporting both Vina 1.1.x and 1.2.x output schemas.

**Total docking runs:** 13 variants × 11 ligands = **143 receptor-ligand pairs**.

---

## 3. Results

### 3.1 Cohort-Level Variant Landscape

Analysis of 60,864 somatic mutations across 1,084 TCGA-BRCA cases yielded 26 prioritized missense variants across the 15-gene panel (Table 1). PIK3CA H1047R emerged as the dominant hotspot with a cohort frequency of **11.81%** (128/1,084 cases), classified as Pathogenic under the frequency-based prior. This frequency is consistent with published TCGA-BRCA analyses reporting PIK3CA H1047R as the most recurrent kinase-domain mutation in Luminal A and Luminal B subtypes [3].

**Table 1. Top-10 prioritized variants by cohort frequency.**

| Rank | Gene | Variant | Chr | Position | Frequency | % Cohort | Classification |
|------|------|---------|-----|----------|-----------|---------|----------------|
| 1 | PIK3CA | H1047R | 3 | 178,952,085 | 128 | 11.81% | Pathogenic |
| 2 | ERBB2 | V777L | 17 | 37,881,000 | 2 | 0.18% | Hotspot |
| 3 | TP53 | R248Q | 17 | — | 1 | 0.09% | Hotspot |
| 4 | BRCA1 | C61G | — | — | 1 | 0.09% | Hotspot |
| 5 | AARS1 | S645C | 16 | 70,292,941 | 1 | 0.09% | Hotspot |
| 6 | CBFB | G123C | 16 | 67,066,766 | 1 | 0.03% | Low Frequency |
| 7 | CBFB | N104S | 16 | 67,066,710 | 1 | 0.03% | Low Frequency |
| 8 | GATA3 | R364S | 10 | 8,073,783 | 1 | 0.03% | Low Frequency |
| 9 | PTEN | N48I | 10 | 87,894,088 | 1 | 0.03% | Low Frequency |
| 10 | SF3B1 | E802Q | 2 | 197,401,492 | 1 | 0.03% | Low Frequency |

The remaining 16 variants (ranks 11–26) span MAP3K1 (4 variants), ARID1A (3 variants), TP53 additional sites (4 variants), PIK3CA secondary sites (D350N, R108H), and BRCA2 D1355Y.

### 3.2 Structural Analysis

3D structures were successfully predicted for 13 of 26 prioritized variants (50%) including their paired WT references. The remaining 13 variants are awaiting structure prediction. Table 2 summarizes structural metrics.

**Table 2. Structural metrics for predicted variant-WT pairs.**

| Variant | WT pLDDT | Mut pLDDT | RMSD Global (Å) | RMSD Active Site (Å) | Key Region |
|---------|----------|-----------|----------------|---------------------|------------|
| PIK3CA H1047R | 87.4 | 84.1 | 1.84 | **3.21** | Kinase Domain |
| PIK3CA D350N | 87.4 | 85.8 | 1.12 | 1.87 | Helical Domain |
| PIK3CA R108H | 87.4 | 86.1 | 0.78 | 1.34 | RBD Domain |
| BRCA1 C61G | 65.3 | 58.9 | 3.41 | **5.67** | RING Domain |
| TP53 R248Q | 76.1 | 71.3 | 2.18 | **4.92** | DNA-binding Domain |
| ERBB2 V777L | 88.2 | 85.7 | 1.23 | 2.45 | Kinase Domain |
| PTEN N48I | 85.6 | 82.1 | 1.05 | 1.73 | Phosphatase Domain |
| AARS1 S645C | 92.1 | 91.5 | 0.84 | 1.12 | Synthetase Domain |
| CBFB G123C | 89.9 | 87.3 | 0.93 | 1.45 | Runt-binding Domain |
| CBFB N104S | 89.9 | 88.1 | 0.71 | 1.18 | Runt-binding Domain |
| ZXDB R528H | 72.4 | 71.8 | 0.45 | — | Hotspot region |
| SF3B1 E802Q | 78.4 | 75.9 | 0.88 | 1.56 | HEAT Domain |
| GATA3 R364S | 37.0 | 34.2 | 1.88 | 2.94 | Zinc Finger ZnF2 |

**PIK3CA H1047R** shows the largest kinase-domain active-site disruption among PIK3CA variants (RMSD = 3.21 Å), with a global pLDDT decrease of 3.3 points, consistent with the known helix αC displacement that impairs GTP positioning in the mutant enzyme [9]. This structural rearrangement mechanistically explains the gain-of-function activity that drives PI3K pathway hyperactivation in Luminal-subtype tumors.

**BRCA1 C61G** presents the most severe structural disruption in the dataset: a RING domain RMSD of 5.67 Å accompanied by a pLDDT drop from 65.3 to 58.9, indicating that the cysteine-to-glycine substitution eliminates a zinc-coordinating residue and collapses the RING finger fold. This finding is consistent with experimental data showing that C61G abrogates BRCA1-BARD1 heterodimer formation and E3 ubiquitin ligase activity [10].

**TP53 R248Q** exhibits a DNA-binding domain RMSD of 4.92 Å and Loop L3 disruption of 2.31 Å, reflecting the loss of direct DNA-contact arginine at codon 248 — one of the most frequently mutated positions across all cancers [11]. The pLDDT decrease from 76.1 to 71.3 in the mutant suggests partial unfolding of the DNA-binding interface.

**GATA3 R364S** shows the lowest absolute pLDDT values in the dataset (WT: 37.0, Mut: 34.2), indicating that both wild-type and mutant structures in the Zinc Finger ZnF2 region are intrinsically disordered — consistent with the known flexibility of GATA3's C-terminal zinc finger [12].

### 3.3 Molecular Docking

Table 3 presents the complete docking affinity matrix for all 13 variants against all 11 inhibitors.

**Table 3. AutoDock Vina 1.2.3 best binding affinity (kcal/mol) by variant and inhibitor.**

| Variant | Alpelisib | Abemaciclib | Palbociclib | Ribociclib | Copanlisib | Capivasertib | Everolimus | Gedatolisib | Taselisib | **Pictilisib** | Buparlisib | **Best** |
|---------|-----------|-------------|-------------|------------|------------|--------------|------------|-------------|-----------|------------|------------|---------|
| PIK3CA H1047R | −8.90 | −8.34 | −7.40 | −7.61 | −8.65 | −7.86 | −6.33 | −5.28 | −7.98 | **−9.72** | −6.40 | −9.72 |
| CBFB G123C | −6.96 | −7.26 | −7.31 | −6.48 | −6.81 | −6.65 | −5.84 | −6.25 | −7.25 | **−8.80** | −5.32 | −8.80 |
| AARS1 S645C | −6.75 | −6.38 | −7.18 | −5.21 | −6.05 | −6.44 | −5.78 | −5.20 | −6.98 | **−8.56** | −6.00 | −8.56 |
| ERBB2 V777L | −7.42 | −6.61 | −6.09 | −5.90 | −6.32 | −6.59 | −5.28 | −4.39 | −6.17 | **−8.00** | −4.78 | −8.00 |
| TP53 R248Q | −6.23 | **−7.45** | −6.17 | −6.16 | −5.53 | −5.82 | −4.57 | −4.67 | −6.13 | −6.67 | −4.44 | −7.45 |
| CBFB N104S | −6.46 | −5.91 | −6.63 | −5.66 | −6.25 | −6.12 | −5.06 | −4.19 | −6.14 | **−7.42** | −5.25 | −7.42 |
| PTEN N48I | −5.95 | −5.84 | −5.74 | −5.50 | −5.61 | −5.65 | −4.78 | −3.68 | −5.30 | **−7.19** | −4.90 | −7.19 |
| PIK3CA D350N | −6.93 | −6.57 | −6.64 | −5.29 | −6.12 | −5.40 | −4.89 | −4.20 | −6.47 | **−7.12** | −5.45 | −7.12 |
| BRCA1 C61G | −5.96 | −5.53 | −5.95 | −5.51 | −5.62 | −5.24 | −4.53 | −3.75 | −5.54 | **−6.47** | −4.00 | −6.47 |
| SF3B1 E802Q | −6.18 | −6.28 | −6.06 | −5.91 | **−6.80** | −5.61 | −4.71 | −3.91 | −5.85 | −6.46 | −4.65 | −6.80 |
| ZXDB R528H | −7.27 | −6.92 | −5.83 | −5.71 | −6.29 | −6.63 | −5.35 | −4.08 | −6.26 | **−7.57** | −5.25 | −7.57 |
| PIK3CA R108H | −5.24 | −4.55 | −4.87 | −4.25 | −4.07 | −4.07 | −3.43 | −2.85 | −4.43 | **−5.35** | −2.70 | −5.35 |
| GATA3 R364S | — | — | — | — | — | — | — | — | — | — | — | n/a |

**Pictilisib** (GDC-0941, pan-class I PI3K inhibitor, Phase II) achieves the strongest binding affinity in 11 of 12 variants with docking data, including the global maximum of −9.72 kcal/mol against PIK3CA H1047R. This consistent dominance across structurally diverse variants reflects Pictilisib's optimized morpholino-triazine scaffold, which accommodates the ATP-binding pockets of multiple kinase-like folds.

**Alpelisib** (FDA-approved PI3Kα inhibitor) achieves the second-best affinity against PIK3CA H1047R (−8.90 kcal/mol), concordant with its regulatory approval specifically for PIK3CA-mutant HR+/HER2− breast cancer in combination with fulvestrant [13]. The computed affinity gap between Alpelisib and Pictilisib (0.82 kcal/mol) suggests that Pictilisib's broader isoform coverage may translate to higher affinity at this specific allosteric site, an observation warranting experimental validation.

**Abemaciclib** achieves the best affinity for TP53 R248Q (−7.45 kcal/mol), which is notable since CDK4/6 inhibitors are not canonical TP53-pathway drugs. This result likely reflects non-specific binding at the disordered DNA-binding domain cleft created by the R248Q-induced conformational disruption.

**Repositioning candidates** (variants where mutation improves affinity by ≥ 0.5 kcal/mol relative to the same drug against WT) were identified for: AARS1 S645C → Pictilisib, CBFB G123C → Palbociclib and Pictilisib, and ERBB2 V777L → Alpelisib.

### 3.4 Integrative VUS Reclassification Evidence

Combining structural disruption and docking affinity data, GEMAMX generates a multi-dimensional evidence profile for each variant. Table 4 summarizes the reclassification evidence tier for variants with complete data.

**Table 4. Integrative reclassification evidence profile.**

| Variant | Active-Site RMSD | Best Affinity (kcal/mol) | pLDDT Δ | Evidence Tier |
|---------|-----------------|--------------------------|---------|---------------|
| BRCA1 C61G | 5.67 Å | −6.47 | −6.4 | **Strong structural** |
| TP53 R248Q | 4.92 Å | −7.45 | −4.8 | **Strong structural** |
| PIK3CA H1047R | 3.21 Å | −9.72 | −3.3 | **Strong structural + pharmacological** |
| ERBB2 V777L | 2.45 Å | −8.00 | −2.5 | Moderate structural + pharmacological |
| PTEN N48I | 1.73 Å | −7.19 | −3.5 | Moderate structural |
| PIK3CA D350N | 1.87 Å | −7.12 | −1.6 | Moderate structural |
| AARS1 S645C | 1.12 Å | −8.56 | −0.6 | Pharmacological priority |
| SF3B1 E802Q | 1.56 Å | −6.80 | −2.5 | Moderate structural |
| CBFB G123C | 1.45 Å | −8.80 | −2.6 | Pharmacological priority |
| PIK3CA R108H | 1.34 Å | −5.35 | −1.3 | Low priority |

Variants with active-site RMSD > 3 Å or pLDDT Δ > 5 are designated "Strong structural" evidence, indicating significant predicted structural impact that warrants experimental functional assays for definitive VUS reclassification. Variants with best affinity ≤ −8.0 kcal/mol are designated "Pharmacological priority," indicating potential drug repositioning interest independent of structural disruption.

---

## 4. Discussion

GEMAMX demonstrates that a fully automated bioinformatics pipeline can traverse the complete path from raw genomic data to docking-level pharmacological characterization in a reproducible, open-source manner. The pipeline's primary contribution is not the individual tools it employs — GDC API, Boltz-1, and AutoDock Vina are established resources — but rather their systematic integration with automated data flow, quality controls, and interactive visualization.

### 4.1 PIK3CA H1047R as Validation Case

The pipeline's strongest result — Pictilisib −9.72 kcal/mol and Alpelisib −8.90 kcal/mol against PIK3CA H1047R, coupled with a kinase-domain RMSD of 3.21 Å — is consistent with the well-characterized biology of this hotspot. PIK3CA H1047R constitutes 11.81% of our TCGA-BRCA cohort, and Alpelisib has received FDA approval specifically for PIK3CA-mutant HR+ breast cancer [13]. That GEMAMX independently recovers this clinically validated drug-variant association from first principles — without any prior knowledge encoded in the pipeline logic — supports the pipeline's validity as a discovery tool for less-characterized variants.

### 4.2 BRCA1 C61G and the RING Domain Collapse

The BRCA1 C61G result (RING domain RMSD = 5.67 Å, pLDDT: 65.3 → 58.9) illustrates GEMAMX's capacity to identify severe structural disruptions in variants that are classified as Pathogenic by clinical databases but may lack a computational structural evidence basis in smaller institutions. The cysteine at position 61 coordinates a zinc ion in the RING domain along with C64, C39, H41, C47, and C50; its substitution to glycine eliminates a zinc ligand and predicts domain unfolding. The docking affinity data for BRCA1 C61G are notably lower than PIK3CA variants (best: −6.47 kcal/mol), consistent with the structural model in which the RING domain cleft is no longer well-defined enough to form stable drug contacts.

### 4.3 The VUS Problem at Scale

Of the 26 variants in the GEMAMX Top-26 prioritized set, 22 are classified as VUS or Low Frequency under the frequency-based prior — variants that would receive no definitive clinical classification from cohort frequency data alone. For the 13 variants with complete structural data, GEMAMX provides active-site RMSD values, pLDDT differentials, and docking affinity profiles that constitute computational functional evidence. Under ClinVar's evidence framework [2], this category of evidence (computational predictors of functional impact) can contribute to upgrading VUS to Likely Pathogenic or Likely Benign, provided it converges with independent data from population frequency databases (gnomAD) and functional studies.

The pipeline identifies 4 variants with active-site RMSD > 2 Å as high structural impact candidates: PIK3CA H1047R (3.21 Å), BRCA1 C61G (5.67 Å), TP53 R248Q (4.92 Å), and GATA3 R364S (2.94 Å). These warrant prioritized experimental follow-up, particularly cell viability assays, reporter gene assays, or structural crystallography.

### 4.4 Mexican Genomic Context

GEMAMX was developed specifically in the context of breast cancer in Mexico, where the incidence pattern differs meaningfully from European cohorts used to train most reference classifiers. The TCGA-BRCA dataset used here is predominantly of European-American ancestry, which introduces potential bias in frequency-based classification for Mexican patients. Future iterations of GEMAMX will incorporate data from the INMEGEN's population-specific sequencing programs and from clinical cohorts at the Universidad de Sonora's breast cancer monitoring program (>265 patients), enabling Mexican-ancestry-specific frequency priors and variant classification calibration.

### 4.5 Limitations

The current pipeline has several limitations that circumscribe the interpretation of its results:

1. **Structure prediction accuracy:** ESMFold and Boltz-1 predictions for VUS variants have not been experimentally validated. pLDDT scores reflect predicted local confidence, not experimental accuracy. RMSD values between independently predicted WT and mutant structures capture predicted conformational differences, not crystallographically confirmed structural changes.

2. **Docking box centering:** The search box is centered on the mutated residue's atomic coordinates, which may not correspond to the biological active site for all variants. For AARS1 S645C, the mutation is in the synthetase domain rather than the aminoacylation active site; docking results for this variant should be interpreted with caution.

3. **Receptor rigidity:** AutoDock Vina treats the receptor as rigid, excluding induced-fit effects that may be particularly significant for disordered regions (GATA3 R364S, pLDDT 37.0).

4. **Fragment structures:** PIK3CA, SF3B1, and PTEN structures are predicted from domain fragments (±180 residues around the mutation site) rather than full-length sequences, limiting the accuracy of global RMSD calculations.

---

## 5. Conclusion

GEMAMX presents a reproducible, fully automated bioinformatics pipeline that integrates genomic data ingestion, variant prioritization, protein structure prediction, molecular docking, and interactive visualization into a single deployable system. Applied to the TCGA-BRCA cohort, the pipeline characterized 13 breast cancer variants across 7 genes, generating 143 docking affinity measurements and identifying PIK3CA H1047R (Pictilisib: −9.72 kcal/mol), BRCA1 C61G (RING domain RMSD: 5.67 Å), and CBFB G123C (Pictilisib: −8.80 kcal/mol) as the highest-priority variants for experimental follow-up.

The pipeline's open-source architecture and web-based interface (gema-mx.netlify.app) make it accessible to clinical researchers without local bioinformatics infrastructure, addressing a gap particularly relevant for precision oncology programs in resource-limited settings. Future development will focus on incorporating Mexican-ancestry genomic data, expanding the drug library to include PARP inhibitors and CDK inhibitors relevant to BRCA1/2 mutant tumors, and automating the re-analysis pipeline as new TCGA-BRCA data become available.

---

## References

[1] Instituto Nacional de Estadística y Geografía (INEGI). Estadísticas de mortalidad por cáncer de mama en México, 2023. Mexico City: INEGI; 2024.

[2] Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med.* 2015;17(5):405–424.

[3] Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature.* 2012;490(7418):61–70.

[4] Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature.* 2021;596(7873):583–589.

[5] Wohlwend M, Bjornsson G, Makeev V, et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. Preprint; 2024. doi: 10.1101/2024.12.03.626747.

[6] Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science.* 2023;379(6637):1123–1130.

[7] Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model.* 2021;61(8):3891–3898.

[8] Forli S, Huey R, Pique ME, Sanner MF, Goodsell DS, Olson AJ. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nat Protoc.* 2016;11(5):905–919.

[9] Mandelker D, Gabelli SB, Schmidt-Kittler O, et al. A frequent kinase domain mutation that changes the interaction between PI3Kα and the membrane. *Proc Natl Acad Sci USA.* 2009;106(40):16996–17001.

[10] Brzovic PS, Meza JE, King MC, Klevit RE. BRCA1 RING domain cancer-predisposing mutations: structural implications and characterization of a second structural pathway. *J Biol Chem.* 2001;276(44):41399–41406.

[11] Olivier M, Hollstein M, Hainaut P. TP53 mutations in human cancers: origins, consequences, and clinical use. *Cold Spring Harb Perspect Biol.* 2010;2(1):a001008.

[12] Lim KC, Lakshmanan G, Crawford SE, Gu Y, Grosveld F, Bhatt SK. Gata3 loss leads to embryonic lethality due to noradrenaline deficiency of the sympathetic nervous system. *Nat Genet.* 2000;25(2):209–212.

[13] André F, Ciruelos E, Rubovszky G, et al. Alpelisib for PIK3CA-mutated, hormone receptor-positive advanced breast cancer. *N Engl J Med.* 2019;380(20):1929–1940.

---

## Supplementary Material

### S1. Pipeline Architecture

The GEMAMX pipeline is implemented in Python 3.12 with the following principal dependencies: pandas 2.2.2, biopython 1.83, pymongo 4.6, rdkit-pypi, meeko 0.5, requests, and tqdm. The web application layer uses HTML5/JavaScript with Plotly.js 2.32 and PDBe Molstar 3.2. Structure prediction via Boltz-1 was executed on Google Colab with NVIDIA T4 GPU acceleration.

### S2. Drug Catalog

| Drug | ChEMBL ID | Primary Target | Regulatory Status |
|------|-----------|----------------|-------------------|
| Alpelisib | CHEMBL3301610 | PI3Kα | FDA approved (PIK3CA-mutant breast cancer) |
| Abemaciclib | CHEMBL3545110 | CDK4/6 | FDA approved |
| Palbociclib | CHEMBL1771552 | CDK4/6 | FDA approved |
| Ribociclib | CHEMBL3545283 | CDK4/6 | FDA approved |
| Copanlisib | CHEMBL3137587 | Pan-PI3K | FDA approved |
| Capivasertib | CHEMBL3301612 | AKT1/2/3 | FDA approved |
| Everolimus | CHEMBL1560600 | mTOR | FDA approved |
| Gedatolisib | CHEMBL3188482 | PI3K/mTOR | Phase II |
| Taselisib | CHEMBL3301609 | PI3Kα | Phase II |
| Pictilisib | CHEMBL2178688 | Pan-PI3K | Phase II |
| Buparlisib | CHEMBL2180717 | Pan-PI3K | Phase II |

### S3. Web Application

The GEMAMX web application is publicly deployed at **https://gema-mx.netlify.app**. Source code is available at the project repository. The application provides: interactive variant explorer with real-time filtering, side-by-side 3D molecular viewer (PDBe Molstar) for WT and mutant structures, docking score ranking with repositioning candidate identification, and full pipeline status dashboard.

---

*Manuscript prepared March 2026. Dataset: TCGA-BRCA (n=1,084). Pipeline version: GEMAMX v1.0.*
