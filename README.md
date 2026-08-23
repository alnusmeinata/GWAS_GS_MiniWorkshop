# Mini-Workshop: Accelerating Breeding with GWAS & Genomic Selection (GS)

Welcome! This 105-minute hands-on workshop is designed for breeders and plant/animal scientists with a basic background in breeding theory.

The goal of this session is to demonstrate how genome-wide association studies (GWAS) and genomic selection (GS) accelerate genetic gain ($\Delta G$) by shortening breeding cycles and improving selection accuracy—without getting bogged down in complex programming.

---

## Workshop Structure & Agenda

| Part | Module | Time | Focus / Key Deliverable |
| :--- | :--- | :--- | :--- |
| **Part 1** | **Briefing & Theory Review** | 15 min | Overview of core GWAS (MAS/QTL) vs. GS (GEBVs) principles and workshop workflow. |
| **Part 2** | **Software Installation & Setup Check** | 45 min | Verification and troubleshooting of R, RStudio, and required packages across participant environments. |
| **Part 3** | **Hands-On Pipeline & Interpretation** | 45 min | Running the 12-script workflow line-by-line: genotype QC, kinship, GWAS, GS models, and interpreting output plots/results. |

---

## Repository Structure

```text
gwas-gs-workshop/
├── README.md
├── Computer requirement and setup.pdf
├── 01_Theory review/
│   └── theoretical_bridge_and_data.pdf
├── 02_dataset/
│   ├── SampleGenotypeData.hmp.xlsx
│   └── Phenotype data.xlsx
└── 03_script/
    ├── EWINDO_0_GenoConversion.R
    ├── EWINDO_1_Imputation.txt
    ├── EWINDO_1a_VCF_Decompression.R
    ├── EWINDO_2_SNPs_thinning.txt
    ├── EWINDO_3_VCF2CSV.R
    ├── EWINDO_4_Population_Kinship.R
    ├── EWINDO_5_Phenotype_treatment.R
    ├── EWINDO_6_Genomic Heritability.R
    ├── EWINDO_7_GWAS.R
    ├── EWINDO_7a_ManhattanPlot.R
    ├── EWINDO_8_GS_Gblup.R
    └── EWINDO_8a_GS_plotting.R
