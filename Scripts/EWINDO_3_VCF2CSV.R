# ==============================================================================
# MODULE 1C: Convert Full & Thinned VCFs to Numeric Matrices (0, 1, 2)
# ==============================================================================
# Inputs : 1. geno_imputed.vcf (or geno_imputed.vcf.gz)
#          2. geno_thin.vcf (or geno_thin.recode.vcf)
# Outputs: 1. geno_full_numeric.csv  (Full matrix for Kinship & GWAS)
#          2. geno_thin_numeric.csv  (Thinned matrix for PCA Population Structure)
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# Helper Function: Parse VCF to Numeric Matrix (0, 1, 2)
parse_vcf_to_numeric <- function(vcf_path) {
  cat(" Reading VCF file:", vcf_path, "...\n")
  
  # Read VCF skipping header meta-lines
  vcf_dt <- fread(vcf_path, skip = "#CHROM", header = TRUE)
  
  # Extract SNP IDs (fallback to Chrom_Pos if ID is missing or '.')
  snp_ids <- vcf_dt$ID
  missing_id <- snp_ids == "." | is.na(snp_ids)
  if (any(missing_id)) {
    snp_ids[missing_id] <- paste(vcf_dt$`#CHROM`[missing_id], vcf_dt$POS[missing_id], sep = "_")
  }
  
  sample_cols <- colnames(vcf_dt)[10:ncol(vcf_dt)]
  gt_mat      <- as.matrix(vcf_dt[, 10:ncol(vcf_dt), with = FALSE])
  
  # Clean GT strings (removes likelihoods/depth if attached, e.g. "0|1:0.99" -> "0|1")
  gt_clean <- apply(gt_mat, 2, function(x) sub(":.*", "", x))
  
  # Initialize numeric output matrix
  num_mat <- matrix(NA_integer_, nrow = nrow(gt_clean), ncol = ncol(gt_clean))
  rownames(num_mat) <- snp_ids
  colnames(num_mat) <- sample_cols
  
  # Convert GT strings to numeric dosage
  num_mat[gt_clean %in% c("0/0", "0|0")] <- 0
  num_mat[gt_clean %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  num_mat[gt_clean %in% c("1/1", "1|1")] <- 2
  
  return(num_mat)
}

# --- 1. Process Full Imputed VCF ---
cat("\n[1/2] Converting Full Imputed VCF to Numeric Matrix...\n")

full_vcf <- if (file.exists("geno_imputed.vcf")) {
  "geno_imputed.vcf"
} else if (file.exists("geno_imputed.vcf.gz")) {
  "geno_imputed.vcf.gz"
} else {
  stop("Error: Could not find 'geno_imputed.vcf' or 'geno_imputed.vcf.gz'!")
}

geno_full_num <- parse_vcf_to_numeric(full_vcf)
write.csv(geno_full_num, "geno.csv", row.names = TRUE)
cat(" -> Successfully saved 'geno_full_numeric.csv' [", nrow(geno_full_num), "SNPs x", ncol(geno_full_num), "Samples]\n")


# --- 2. Process Thinned VCF ---
cat("\n[2/2] Converting Thinned VCF to Numeric Matrix...\n")

thin_vcf <- if (file.exists("geno_thin.vcf")) {
  "geno_thin.vcf"
} else if (file.exists("geno_thin.recode.vcf")) {
  "geno_thin.recode.vcf"
} else {
  stop("Error: Could not find 'geno_thin.vcf' or 'geno_thin.recode.vcf'!")
}

geno_thin_num <- parse_vcf_to_numeric(thin_vcf)
write.csv(geno_thin_num, "geno_thin.csv", row.names = TRUE)
cat(" -> Successfully saved 'geno_thin_numeric.csv' [", nrow(geno_thin_num), "SNPs x", ncol(geno_thin_num), "Samples]\n")

cat("\n==============================================================================\n")
cat(" CONVERSION COMPLETE: 'geno_full_numeric.csv' and 'geno_thin_numeric.csv' created!\n")
cat(" Ready for Kinship, PCA, and GWAS modules.\n")
cat("==============================================================================\n")

