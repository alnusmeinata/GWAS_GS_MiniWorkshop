# Mini-Workshop: Accelerating Breeding with GWAS & Genomic Selection (GS)

Welcome! This 90-minute hands-on workshop is designed for breeders and plant/animal scientists with basic breeding theory background

The goal of this session is to demonstrate how genome-wide association studies (GWAS) and genomic selection (GS) accelerate genetic gain ($\Delta G$) by shortening breeding cycles and improving selection accuracy—without getting bogged down in complex programming.

---

## Workshop Structure & Agenda

| Part | Module | Time | Focus / Key Deliverable |
| :--- | :--- | :--- | :--- |
| **Part 1** | **Theory Review / Small Lecture** | 20 min | Core principles of GWAS (MAS/QTL) vs. GS (GEBVs) and data structure overview. |
| **Part 2** | **Software Installation Check** | 10 min | Verification of R, RStudio, and required package installations on participant laptops. |
| **Part 3** | **Initial Familiarity with R** | 15 min | Learning RStudio interface basics: script window, console, environment, and `Ctrl + Enter`. |
| **Part 4** | **Mini-Workshop: Follow Along** | 45 min | Running `run_gwas_gs.R` line-by-line: allele coding, Manhattan plots, kinship, and GEBVs. |

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
    └── GWAS_GS_demo.R
