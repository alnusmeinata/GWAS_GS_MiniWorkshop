# 0. Install and load required packages
required_packages <- c("BGLR", "tidyverse", "parallel", "data.table", "dplyr", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(BGLR)
library(tidyverse)
library(parallel)
library(data.table)
library(dplyr)
library(ggplot2)

# Settings
mc.cores <- 4
method <- "BRR"
output_prefix <- "gh_"
filename_geno <- "geno.csv"
filename_pheno <- "pheno_filtered.csv"

# Reduced iterations for fast execution testing
nIter <- 100
burnIn <- 20
thin <- 2

# Trait sets matching pheno_filtered.csv
phe_sets <- list(
  c("yield", "plant_height"),
  c("leaf_area", "seed_size")
)
traits <- unlist(phe_sets)

# 1. Read genotype data (No SNP duplicates)
geno_temp <- read.table(filename_geno, sep = ',', row.names = 1, header = TRUE, check.names = FALSE)
geno <- t(geno_temp)
rownames(geno) <- as.character(rownames(geno))

# 2. Read phenotype data and standardize ID column
pheno <- read.csv(filename_pheno, stringsAsFactors = FALSE)
if ("accession" %in% colnames(pheno)) {
  pheno <- pheno %>% rename(ID = accession)
}
pheno$ID <- as.character(pheno$ID)

# 3. Genomic Heritability Estimation via BGLR
results <- list()

for (phe_set in phe_sets) {
  
  # Subset phenotype and remove missing rows across target trait set
  pheno_subset <- pheno %>%
    select(ID, all_of(phe_set)) %>%
    filter(if_all(all_of(phe_set), ~ !is.na(.)))
  
  # Match IDs between genotype and phenotype
  matched_ids <- intersect(rownames(geno), pheno_subset$ID)
  
  # Subset and align genotype matrix
  geno_matched <- geno[matched_ids, , drop = FALSE]
  geno_scaled <- scale(geno_matched) / sqrt(ncol(geno_matched))
  
  # Subset and align phenotype matching genotype order
  pheno_matched <- pheno_subset %>%
    filter(ID %in% matched_ids) %>%
    arrange(match(ID, matched_ids))
  
  stopifnot(all(pheno_matched$ID == rownames(geno_scaled)))
  
  phe_results <- data.frame()
  
  for (phe in phe_set) {
    y <- pheno_matched[[phe]] %>% unlist()
    
    valid_ids <- pheno_matched$ID[!is.na(y)]
    
    X_clean <- geno_scaled[match(valid_ids, rownames(geno_scaled)), , drop = FALSE]
    X_clean <- X_clean[, colSums(is.na(X_clean)) == 0, drop = FALSE]
    
    y_clean <- y[match(valid_ids, pheno_matched$ID)]
    
    if (anyNA(X_clean)) {
      stop("NAs found in X_clean before BGLR for trait: ", phe)
    }
    
    for (j in 1:10) {
      set.seed(j)
      
      BGLR(
        y = y_clean,
        ETA = list(list(X = X_clean, model = method, saveEffects = TRUE)),
        saveAt = sprintf("%s%s_%d_", output_prefix, phe, j),
        nIter = nIter,
        burnIn = burnIn,
        thin = thin,
        verbose = FALSE
      )
      
      varU <- scan(sprintf("%s%s_%d_ETA_1_varB.dat", output_prefix, phe, j), quiet = TRUE)
      varE <- scan(sprintf("%s%s_%d_varE.dat", output_prefix, phe, j), quiet = TRUE)
      
      h2 <- varU / (varU + varE)
      discard_idx <- floor(burnIn / thin)
      if (discard_idx > 0 && discard_idx < length(h2)) {
        h2 <- h2[-(1:discard_idx)]
      }
      
      phe_results <- rbind(phe_results, data.frame(phe = phe, chain = j, h2 = mean(h2)))
    }
  }
  
  results[[paste(phe_set, collapse = "_")]] <- phe_results
}

# 4. Combine and export results
final_results <- bind_rows(results)
final_results$h2 <- as.numeric(final_results$h2)

print(final_results)
write.csv(final_results, sprintf("%sgenomic_heritability.csv", output_prefix), row.names = FALSE)

# 5. Summarize and Plot Max Heritability
result_summary <- final_results %>%
  group_by(phe) %>%
  summarise(max_h2 = max(h2), .groups = 'drop') %>%
  mutate(phe = factor(phe, levels = traits))

p <- ggplot(result_summary, aes(x = phe, y = max_h2)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = round(max_h2, 3)), vjust = -0.3, size = 4) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = "Trait",
    y = "Genomic Narrow-Sense Heritability (h²)",
    title = "Max Genomic Heritability across Traits"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

print(p)
ggsave('max_heritability_plot.png', plot = p, width = 8, height = 5, dpi = 300)
