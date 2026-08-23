# 0. Automatic Package Installation & Loading
required_packages <- c("sommer", "tidyverse", "data.table", "parallel", "caret")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(sommer)
library(tidyverse)
library(data.table)
library(parallel)
library(caret)

# 1. Configuration & Settings
mc.cores <- min(8, detectCores() - 1)
filename_geno_thin <- "geno.csv"
filename_pheno <- "pheno_filtered.csv"
gwas_dir <- "GWAS_QK"
output_file <- "gp_gblup_10fold_cv_results.csv"

# 2. Read Phenotype Data
pheno_raw <- fread(filename_pheno, header = TRUE, data.table = FALSE)
if ("accession" %in% colnames(pheno_raw)) {
  pheno_raw <- pheno_raw %>% dplyr::rename(ID = accession)
}
pheno_raw$ID <- trimws(as.character(pheno_raw$ID))
traits <- setdiff(colnames(pheno_raw), "ID")

# 3. Read & Format Genotype Data
snp <- read.csv(filename_geno_thin, header = TRUE, row.names = 1, check.names = FALSE)
geno_full <- t(as.matrix(snp))
storage.mode(geno_full) <- "numeric"
rownames(geno_full) <- trimws(rownames(geno_full))

# Fix optional 'X' prefix naming mismatch between pheno and geno
if (!any(pheno_raw$ID %in% rownames(geno_full)) && any(paste0("X", pheno_raw$ID) %in% rownames(geno_full))) {
  pheno_raw$ID <- paste0("X", pheno_raw$ID)
}

total_snps <- ncol(geno_full)
n_loci_seq <- unique(c(seq(1000, total_snps, by = 1000), total_snps))

# 4. Prepare Fold Splits & Ranked Loci per Trait
set.seed(123) # For reproducible 10-fold CV partitions
trait_data_list <- list()
iteration_list <- list()

for (trait in traits) {
  # Extract non-NA phenotyped samples matching genotype matrix
  pheno_trait <- pheno_raw %>%
    dplyr::select(ID, all_of(trait)) %>%
    filter(!is.na(.[[trait]])) %>%
    filter(ID %in% rownames(geno_full))
  
  if (nrow(pheno_trait) < 10) {
    message(sprintf("Skipping trait %s due to insufficient samples (<10).", trait))
    next
  }
  
  # Read GWAS intermediate results to rank SNPs by significance
  gwas_path <- file.path(gwas_dir, sprintf("gwas_%s_intermediate result.csv", trait))
  if (!file.exists(gwas_path)) {
    message(sprintf("GWAS result missing for trait %s at %s - Skipping.", trait, gwas_path))
    next
  }
  
  gwas_res <- read.csv(gwas_path, check.names = FALSE)
  p_col <- grep("-log10\\(P\\)", colnames(gwas_res), value = TRUE)[1]
  
  ranked_loci <- gwas_res %>%
    arrange(desc(.data[[p_col]])) %>%
    pull(Marker) %>%
    intersect(colnames(geno_full))
  
  # Create 10-fold CV partition
  folds <- caret::createFolds(pheno_trait[[trait]], k = 10, list = TRUE)
  
  trait_data_list[[trait]] <- list(
    pheno = pheno_trait,
    ranked_loci = ranked_loci,
    folds = folds
  )
  
  # Build grid of parameter iterations for this trait
  for (fold_idx in seq_along(folds)) {
    for (nl in n_loci_seq) {
      iteration_list[[length(iteration_list) + 1]] <- list(
        trait = trait,
        fold = fold_idx,
        n_loci = nl
      )
    }
  }
}

# 5. Define GBLUP Worker Function
run_gblup_cv <- function(item, trait_data_list, geno_full) {
  tryCatch({
    trait <- item$trait
    fold_idx <- item$fold
    n_loci <- item$n_loci
    
    t_info <- trait_data_list[[trait]]
    pheno_df <- t_info$pheno
    test_indices <- t_info$folds[[fold_idx]]
    
    # Subset genotype to top n_loci
    selected_loci <- t_info$ranked_loci[1:min(n_loci, length(t_info$ranked_loci))]
    geno_sub <- geno_full[pheno_df$ID, selected_loci, drop = FALSE]
    
    # Prepare phenotype vector with NA for validation fold
    pheno_model <- pheno_df
    pheno_model$y <- pheno_model[[trait]]
    pheno_model$y[test_indices] <- NA
    pheno_model$id <- factor(pheno_model$ID)
    
    # Compute Kinship Matrix
    K <- sommer::A.mat(geno_sub)
    colnames(K) <- rownames(K) <- pheno_model$id
    
    # Fit GBLUP Model
    ans <- mmer(
      y ~ 1,
      random = ~ vsr(id, Gu = K),
      rcov = ~ units,
      data = pheno_model,
      verbose = FALSE
    )
    
    # Predict and calculate correlation on test set
    predictions <- ans$U$`u:id`$y[test_indices]
    actuals <- pheno_df[[trait]][test_indices]
    corr <- cor(predictions, actuals, use = "complete.obs")
    
    return(data.frame(
      trait = trait,
      fold = fold_idx,
      n_loci = n_loci,
      corr = corr
    ))
  }, error = function(e) {
    return(data.frame(
      trait = item$trait,
      fold = item$fold,
      n_loci = item$n_loci,
      corr = NA
    ))
  })
}

# 6. Execute Parallel 10-Fold Cross-Validation
cl <- makeCluster(mc.cores)
clusterExport(cl, varlist = c("run_gblup_cv", "trait_data_list", "geno_full"), envir = environment())
clusterEvalQ(cl, {
  library(sommer)
  library(dplyr)
})

cat("Starting GBLUP 10-fold Cross-Validation across all traits...\n")
results_list <- parLapply(cl, iteration_list, run_gblup_cv, trait_data_list = trait_data_list, geno_full = geno_full)
stopCluster(cl)

# 7. Format and Save Results
final_results <- bind_rows(results_list) %>% filter(!is.na(corr))
write.csv(final_results, output_file, row.names = FALSE)

cat("\nCross-validation complete. Results saved to:", output_file, "\n")
