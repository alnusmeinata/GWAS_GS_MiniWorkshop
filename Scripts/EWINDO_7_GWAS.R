# 0. Automatic Package Installation & Loading
required_packages <- c("mrMLM", "tidyverse", "data.table", "reshape")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(mrMLM)
library(tidyverse)
library(data.table)
library(reshape)

# 1. Settings & File Paths
method <- "FASTmrEMMA"
nPCs <- 2
filename_geno <- "hapmap_diploid_geno_imputed.hmp.txt"
filename_geno_thin <- "geno_thin.csv"
filename_pheno <- "pheno_filtered.csv"
output_prefix <- "gwas"
dir_qk <- "GWAS_QK"

dir.create(dir_qk, showWarnings = FALSE)

# 2. Read Phenotype Data
pheno_raw <- fread(filename_pheno, header = TRUE, data.table = FALSE)

if ("accession" %in% colnames(pheno_raw)) {
  pheno_raw <- pheno_raw %>% dplyr::rename(ID = accession)
}
pheno_raw$ID <- trimws(as.character(pheno_raw$ID))
traits <- setdiff(colnames(pheno_raw), "ID")

# 3. Read HapMap Genotype Data
hmp <- fread(filename_geno, header = FALSE, stringsAsFactors = TRUE)

# 4. Helper Functions
rename_outputs <- function(phe_set, output_prefix, dir_path) {
  output_files <- c(
    "_intermediate result.csv",
    "_Final result.csv",
    "_Manhattan plot.png",
    "_qq plot.png"
  )
  for (f in output_files) {
    result_file <- file.path(dir_path, paste0("1", f))
    if (file.exists(result_file)) {
      new_name <- file.path(dir_path, paste0(output_prefix, "_", phe_set, f))
      file.rename(result_file, new_name)
    }
  }
}

write_PCA <- function(samples, path, filename_geno_thin, nPCs = 2) {
  if (nPCs <= 0) return(NULL)
  
  # Read with check.names = FALSE to prevent R from modifying sample column headers
  snp <- read.csv(filename_geno_thin, header = TRUE, row.names = 1, check.names = FALSE)
  snp[] <- lapply(snp, as.numeric)
  
  geno_mat <- t(as.matrix(snp))
  storage.mode(geno_mat) <- "numeric"
  
  # Clean sample names
  rownames(geno_mat) <- trimws(rownames(geno_mat))
  samples <- trimws(samples)
  
  # Handle potential 'X' prefix added by R or missing in matching
  if (!any(samples %in% rownames(geno_mat)) && any(paste0("X", samples) %in% rownames(geno_mat))) {
    samples <- paste0("X", samples)
  } else if (!any(samples %in% rownames(geno_mat)) && any(gsub("^X", "", rownames(geno_mat)) %in% samples)) {
    rownames(geno_mat) <- gsub("^X", "", rownames(geno_mat))
  }
  
  avail_samples <- intersect(samples, rownames(geno_mat))
  
  if (length(avail_samples) == 0) {
    stop("write_PCA Error: No sample IDs matched between phenotype/HapMap and geno_thin.csv")
  }
  
  geno_mat <- geno_mat[avail_samples, , drop = FALSE]
  
  # Remove monomorphic (zero variance) SNPs
  col_vars <- apply(geno_mat, 2, var, na.rm = TRUE)
  geno_mat <- geno_mat[, !is.na(col_vars) & col_vars > 0, drop = FALSE]
  
  if (ncol(geno_mat) == 0) {
    stop("write_PCA Error: All SNPs have zero variance for the selected sample subset.")
  }
  
  pca <- prcomp(geno_mat, scale. = TRUE)
  pcs <- as.data.frame(pca$x)
  pcs <- pcs[, paste0("PC", 1:min(nPCs, ncol(pcs))), drop = FALSE]
  pcs <- pcs %>%
    rownames_to_column("<ID>") %>%
    as_tibble()
  
  pcs_out <- rbind(
    c("<PCA>", rep("", ncol(pcs) - 1)),
    colnames(pcs),
    pcs
  )
  
  write.table(pcs_out, file = path, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(path)
}

write_kinship <- function(samples, path, filename_geno_thin) {
  samples <- trimws(samples)
  
  if (file.exists("kinship.txt")) {
    kin_mat <- read.delim("kinship.txt", row.names = 1, check.names = FALSE)
  } else {
    # Compute VanRaden Kinship from geno_thin if kinship.txt is absent
    snp <- read.csv(filename_geno_thin, header = TRUE, row.names = 1, check.names = FALSE)
    geno_mat <- t(as.matrix(snp))
    storage.mode(geno_mat) <- "numeric"
    
    rownames(geno_mat) <- trimws(rownames(geno_mat))
    
    if (!any(samples %in% rownames(geno_mat)) && any(paste0("X", samples) %in% rownames(geno_mat))) {
      samples_match <- paste0("X", samples)
    } else if (!any(samples %in% rownames(geno_mat)) && any(gsub("^X", "", rownames(geno_mat)) %in% samples)) {
      rownames(geno_mat) <- gsub("^X", "", rownames(geno_mat))
      samples_match <- samples
    } else {
      samples_match <- samples
    }
    
    avail <- intersect(samples_match, rownames(geno_mat))
    geno_sub <- geno_mat[avail, , drop = FALSE]
    
    Z <- scale(geno_sub, center = TRUE, scale = FALSE)
    P <- colMeans(geno_sub, na.rm = TRUE) / 2
    K <- (Z %*% t(Z)) / (2 * sum(P * (1 - P), na.rm = TRUE))
    kin_mat <- as.data.frame(K)
  }
  
  rownames(kin_mat) <- trimws(rownames(kin_mat))
  colnames(kin_mat) <- trimws(colnames(kin_mat))
  
  if (!any(samples %in% rownames(kin_mat)) && any(paste0("X", samples) %in% rownames(kin_mat))) {
    samples <- paste0("X", samples)
  }
  
  avail_samples <- intersect(samples, rownames(kin_mat))
  kin_sub <- kin_mat[avail_samples, avail_samples, drop = FALSE] %>%
    rownames_to_column("ID") %>%
    as_tibble()
  
  top_row <- c(as.character(length(avail_samples)), rep("", ncol(kin_sub) - 1))
  kinship_out <- rbind(top_row, as.matrix(kin_sub))
  
  write.table(kinship_out, file = path, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(path)
}

# 5. Execution Loop (Q + K Only)
for (trait_name in traits) {
  cat("\n========================================\n")
  cat("Running Q+K GWAS (FASTmrEMMA) for trait:", trait_name, "\n")
  cat("========================================\n")
  
  # Filter non-NA samples for target trait
  pheno_trait <- pheno_raw %>%
    select(ID, all_of(trait_name)) %>%
    filter(!is.na(.[[trait_name]]))
  
  # Match sample IDs against HapMap header
  hmp_samples <- trimws(unlist(hmp[1, -c(1:11)]))
  valid_samples <- intersect(hmp_samples, pheno_trait$ID)
  
  if (length(valid_samples) == 0) {
    cat("Warning: No matching sample IDs found for trait:", trait_name, "- Skipping.\n")
    next
  }
  
  # Filter HapMap columns
  keep_cols <- c(rep(TRUE, 11), hmp_samples %in% valid_samples)
  hmp_filtered <- hmp[, keep_cols, with = FALSE]
  
  # Format Phenotype data structure
  pheno_matched <- pheno_trait %>% filter(ID %in% valid_samples)
  pheno_formatted <- data.frame(
    V1 = c("<Phenotypes>", pheno_matched$ID),
    V2 = c(trait_name, as.character(pheno_matched[[trait_name]]))
  )
  
  # Generate PCA (Q) and Kinship (K) matrices
  filePS_qk <- write_PCA(valid_samples, file.path(dir_qk, paste0(output_prefix, "_", trait_name, "_PopStr.csv")), filename_geno_thin, nPCs)
  fileKin_qk <- write_kinship(valid_samples, file.path(dir_qk, paste0(output_prefix, "_", trait_name, "_kin.csv")), filename_geno_thin)
  
  # Run mrMLM Q+K model
  mrMLM(
    fileGen = hmp_filtered,
    filePhe = pheno_formatted,
    filePS = filePS_qk,
    fileKin = fileKin_qk,
    Genformat = "Hmp",
    PopStrType = "PCA",
    method = method,
    trait = 1,
    SearchRadius = 20,
    CriLOD = 3,
    DrawPlot = FALSE,
    dir = dir_qk
  )
  
  # Standardize output filenames
  rename_outputs(trait_name, output_prefix, dir_qk)
}
