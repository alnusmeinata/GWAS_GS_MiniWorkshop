# ==============================================================================
# MODULE 1: Genotype Import, Deduplication, Quality Control, & VCF Export
# ==============================================================================
# Inputs : SampleGenotypeData.hmp.xlsx
# Outputs: 1. Genotype_QC_Diagnostics.png  (QC diagnostic histograms)
#          2. geno_qc.vcf                   (QC-filtered VCF ready for Beagle)
# ==============================================================================

# --- 0. Package Setup ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_packages <- c("readxl", "tidyverse", "data.table", "gridExtra")
new_pkgs <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(data.table)
  library(gridExtra)
})

# Parameters
file_hapmap     <- "SampleGenotypeData.hmp.xlsx"
max_snp_missing <- 0.10  # Max 10% missing per SNP
min_maf         <- 0.05  # Min 5% MAF
max_ind_missing <- 0.20  # Max 20% missing per sample


# --- 1. Load Data & Deduplicate Sample Columns ---
cat("\n[1/4] Importing HapMap dataset and removing duplicated samples...\n")

hapmap_raw <- read_excel(file_hapmap, .name_repair = "minimal")
meta_cols  <- hapmap_raw[, 1:11]
geno_cols  <- hapmap_raw[, 12:ncol(hapmap_raw)]

deduplicate_samples <- function(geno_df) {
  sample_names <- colnames(geno_df)
  total_samples <- length(sample_names)
  
  dup_mask <- duplicated(sample_names)
  num_duplicates <- sum(dup_mask)
  unique_dup_names <- unique(sample_names[dup_mask])
  
  cat("==========================================\n")
  cat(" SAMPLE DEDUPLICATION SUMMARY REPORT\n")
  cat("==========================================\n")
  cat(" Total Sample Columns Input  :", total_samples, "\n")
  cat(" Duplicated Sample Columns   :", num_duplicates, "\n")
  cat(" Unique Duplicated Sample IDs:", length(unique_dup_names), "\n")
  if (num_duplicates > 0) {
    cat(" List of Duplicate Samples   :", paste(unique_dup_names, collapse = ", "), "\n")
  }
  cat(" Remaining Unique Samples    :", total_samples - num_duplicates, "\n")
  cat("==========================================\n\n")
  
  return(geno_df[, !dup_mask])
}

geno_clean   <- deduplicate_samples(geno_cols)
hapmap_dedup <- cbind(meta_cols, geno_clean)


# --- 2. Calculate QC Metrics ---
cat("[2/4] Converting alleles to dosage matrix for QC metric computation...\n")

geno_char <- as.matrix(hapmap_dedup[, 12:ncol(hapmap_dedup)])
rownames(geno_char) <- hapmap_dedup$`rs#`

convert_hapmap_to_numeric <- function(geno_mat) {
  num_mat <- matrix(NA, nrow = nrow(geno_mat), ncol = ncol(geno_mat))
  rownames(num_mat) <- rownames(geno_mat)
  colnames(num_mat) <- colnames(geno_mat)
  
  for (i in 1:nrow(geno_mat)) {
    vec <- geno_mat[i, ]
    valid <- vec[!vec %in% c("NN", "N", "XX", "-", "00", "0", "./.")]
    if (length(valid) == 0) next
    
    alleles <- unique(unlist(strsplit(valid, "")))
    if (length(alleles) == 2) {
      ref <- alleles[1]; alt <- alleles[2]
      num_mat[i, vec == paste0(ref, ref)] <- 0
      num_mat[i, vec %in% c(paste0(ref, alt), paste0(alt, ref))] <- 1
      num_mat[i, vec == paste0(alt, alt)] <- 2
    } else if (length(alleles) == 1) {
      ref <- alleles[1]
      num_mat[i, vec == paste0(ref, ref)] <- 0
    }
  }
  return(num_mat)
}

geno_num <- convert_hapmap_to_numeric(geno_char)

snp_missing <- rowMeans(is.na(geno_num))
p_freq      <- rowMeans(geno_num, na.rm = TRUE) / 2
maf         <- pmin(p_freq, 1 - p_freq)
snp_het     <- rowSums(geno_num == 1, na.rm = TRUE) / rowSums(!is.na(geno_num))

ind_missing <- colMeans(is.na(geno_num))
ind_het     <- colSums(geno_num == 1, na.rm = TRUE) / colSums(!is.na(geno_num))

df_snp_qc <- data.frame(SNP = rownames(geno_num), Missingness = snp_missing, MAF = maf, Heterozygosity = snp_het)
df_ind_qc <- data.frame(Sample = colnames(geno_num), Missingness = ind_missing, Heterozygosity = ind_het)


# --- 3. Save Diagnostic QC Plots ---
cat("[3/4] Generating and saving diagnostic QC plots...\n")

p1 <- ggplot(df_snp_qc, aes(x = Missingness)) +
  geom_histogram(bins = 30, fill = "#4E79A7", color = "black", alpha = 0.8) +
  geom_vline(xintercept = max_snp_missing, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(title = "A. SNP Missingness", x = "Missing Rate per SNP", y = "Count") + theme_bw(base_size = 11)

p2 <- ggplot(df_ind_qc, aes(x = Missingness)) +
  geom_histogram(bins = 30, fill = "#F28E2B", color = "black", alpha = 0.8) +
  geom_vline(xintercept = max_ind_missing, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(title = "B. Sample Missingness", x = "Missing Rate per Individual", y = "Count") + theme_bw(base_size = 11)

p3 <- ggplot(df_snp_qc, aes(x = MAF)) +
  geom_histogram(bins = 30, fill = "#E15759", color = "black", alpha = 0.8) +
  geom_vline(xintercept = min_maf, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(title = "C. Minor Allele Frequency (MAF)", x = "MAF", y = "Count") + theme_bw(base_size = 11)

p4 <- ggplot(df_ind_qc, aes(x = Heterozygosity)) +
  geom_histogram(bins = 30, fill = "#76B7B2", color = "black", alpha = 0.8) +
  labs(title = "D. Sample Observed Heterozygosity", x = "Heterozygosity Rate", y = "Count") + theme_bw(base_size = 11)

qc_panel <- grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("Genotype_QC_Diagnostics.png", qc_panel, width = 11, height = 9, dpi = 300)


# --- 4. Apply Filters & Export Clean VCF ---
cat("[4/4] Filtering dataset and exporting directly to VCF format...\n")

keep_snps <- (snp_missing <= max_snp_missing) & (maf >= min_maf)
keep_inds <- ind_missing <= max_ind_missing

hapmap_qc <- hapmap_dedup[keep_snps, c(TRUE, rep(TRUE, 10), keep_inds)]

cat("\n==========================================\n")
cat("       GENOTYPE QC FILTERING SUMMARY      \n")
cat("==========================================\n")
cat(" Initial SNPs            :", nrow(geno_num), "\n")
cat(" SNPs Filtered Out       :", nrow(geno_num) - sum(keep_snps), "\n")
cat(" Remaining Clean SNPs    :", nrow(hapmap_qc), "\n")
cat(" ---------------------------------------- \n")
cat(" Initial Samples         :", ncol(geno_num), "\n")
cat(" Samples Filtered Out    :", ncol(geno_num) - sum(keep_inds), "\n")
cat(" Remaining Clean Samples :", ncol(hapmap_qc) - 11, "\n")
cat("==========================================\n\n")

write_hapmap_to_vcf <- function(hmp_df, vcf_path) {
  meta <- hmp_df[, 1:11]
  geno <- hmp_df[, 12:ncol(hmp_df)]
  vcf_header <- c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
  vcf_rows <- character(nrow(hmp_df))
  
  for (i in 1:nrow(hmp_df)) {
    chrom <- as.character(meta[i, 3]); pos <- as.character(meta[i, 4]); id <- as.character(meta[i, 1])
    vec <- unlist(geno[i, ])
    valid <- vec[!vec %in% c("NN", "N", "XX", "-", "00", "0", "./.")]
    alleles <- unique(unlist(strsplit(valid, "")))
    
    ref <- if (length(alleles) >= 1) alleles[1] else "N"
    alt <- if (length(alleles) >= 2) alleles[2] else "N"
    
    gt_vec <- rep("./.", length(vec))
    gt_vec[vec == paste0(ref, ref)] <- "0/0"
    gt_vec[vec %in% c(paste0(ref, alt), paste0(alt, ref))] <- "0/1"
    gt_vec[vec == paste0(alt, alt)] <- "1/1"
    
    vcf_rows[i] <- paste(c(chrom, pos, id, ref, alt, ".", "PASS", ".", "GT", gt_vec), collapse = "\t")
  }
  
  col_names <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(geno)), collapse = "\t")
  writeLines(c(vcf_header, col_names, vcf_rows), vcf_path)
}

write_hapmap_to_vcf(hapmap_qc, "geno_qc.vcf")

cat("MODULE 1 COMPLETE: 'geno_qc.vcf' and 'Genotype_QC_Diagnostics.png' successfully created!\n")
