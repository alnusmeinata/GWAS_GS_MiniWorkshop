# ==============================================================================
# MODULE 2: Population Structure (PCA with K-Means Elbow) & Kinship Matrix
# ==============================================================================
# Inputs : 1. geno_thin.csv          (Thinned dosage matrix)
#          2. geno.csv               (Full imputed dosage matrix)
# Outputs: 1. PCA.csv                (Taxa + Top 5 PCs only, for GWAS input)
#          2. Kinship.csv            (Genomic Relationship Matrix for GWAS)
#          3. PCA_Elbow_Diagnostics.png (2-panel plot: Elbow Curve + PCA Scatter)
#          4. Kinship_Heatmap.png    (Heatmap of sample relatedness)
# ==============================================================================

# --- 0. Package Setup ---
required_packages <- c("data.table", "ggplot2", "gridExtra")
new_pkgs <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
})


# --- 1. Compute PCA on Thinned Genotype Matrix ---
cat("\n[1/3] Computing Principal Component Analysis (PCA) on thinned SNPs...\n")

# Load thinned numeric matrix (SNPs in rows, Samples in columns)
geno_thin_file <- if (file.exists("geno_thin.csv")) "geno_thin.csv" else "geno_thin_numeric.csv"
geno_thin_raw  <- fread(geno_thin_file, header = TRUE)

snp_ids_thin <- geno_thin_raw[[1]]
sample_ids   <- colnames(geno_thin_raw)[-1]

# Transpose matrix to Samples x SNPs format for PCA
X_thin <- t(as.matrix(geno_thin_raw[, -1, with = FALSE]))
rownames(X_thin) <- sample_ids
colnames(X_thin) <- snp_ids_thin

# Run PCA
pca_res <- prcomp(X_thin, center = TRUE, scale. = FALSE)
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Format PCA dataframe for GWAS Export (Taxa + Top 5 PCs ONLY)
pca_df <- data.frame(
  Taxa = rownames(X_thin),
  PC1  = pca_res$x[, 1],
  PC2  = pca_res$x[, 2],
  PC3  = pca_res$x[, 3],
  PC4  = pca_res$x[, 4],
  PC5  = pca_res$x[, 5]
)

# Export PCA.csv for GWAS covariate input (Clean format)
write.csv(pca_df, "PCA.csv", row.names = FALSE)
cat(" -> Saved 'PCA.csv' (Taxa + Top 5 PCs only for", nrow(pca_df), "samples)\n")


# --- 2. K-Means Clustering & Elbow Method Diagnostic Plots ---
cat("[2/3] Performing K-Means Elbow analysis (K = 1 to 10)...\n")

set.seed(123)  # For reproducible K-means results
max_k  <- min(10, nrow(X_thin) - 1)
scores <- pca_res$x[, 1:min(5, ncol(pca_res$x))]  # Top 5 PCs used for clustering

# Calculate Within-Cluster Sum of Squares (WSS) across K values
wss <- numeric(max_k)
for (k in 1:max_k) {
  km <- kmeans(scores, centers = k, nstart = 25)
  wss[k] <- km$tot.withinss
}

df_elbow <- data.frame(K = 1:max_k, WSS = wss)

# Temporary K-means assignment (K = 3) purely for diagnostic plot coloring
opt_k       <- min(3, max_k)
km_opt      <- kmeans(scores, centers = opt_k, nstart = 25)
pca_plot_df <- cbind(pca_df, Cluster = as.factor(km_opt$cluster))

# Plot 1: K-Means Elbow Curve
p_elbow <- ggplot(df_elbow, aes(x = K, y = WSS)) +
  geom_line(color = "#E15759", linewidth = 1) +
  geom_point(color = "#E15759", size = 3) +
  scale_x_continuous(breaks = 1:max_k) +
  labs(
    title = "A. K-Means Elbow Plot",
    x = "Number of Clusters (K)",
    y = "Total Within-Cluster SS"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Plot 2: PCA Scatter Plot colored by Subpopulation Cluster
p_pca <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    title = sprintf("B. Population Structure (PCA, K = %d)", opt_k),
    x = sprintf("PC1 (%.2f%% Var)", var_explained[1]),
    y = sprintf("PC2 (%.2f%% Var)", var_explained[2]),
    color = "Cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Save combined 2-panel diagnostic figure
pca_panel <- grid.arrange(p_elbow, p_pca, ncol = 2)
ggsave("PCA_Elbow_Diagnostics.png", pca_panel, width = 11, height = 5, dpi = 300)
cat(" -> Saved 'PCA_Elbow_Diagnostics.png'\n")


# --- 3. Compute Kinship Matrix (VanRaden 2008) on Full Genotype Matrix ---
cat("\n[3/3] Calculating VanRaden Kinship Matrix on full imputed SNPs...\n")

# Load full imputed numeric matrix
geno_full_file <- if (file.exists("geno.csv")) "geno.csv" else "geno_full_numeric.csv"
geno_full_raw  <- fread(geno_full_file, header = TRUE)

snp_ids_full <- geno_full_raw[[1]]

# Transpose matrix to Samples x SNPs format
X_full <- t(as.matrix(geno_full_raw[, -1, with = FALSE]))
rownames(X_full) <- sample_ids
colnames(X_full) <- snp_ids_full

# Calculate allele frequencies (p) and centered dosage matrix (Z)
p <- colMeans(X_full, na.rm = TRUE) / 2
Z <- sweep(X_full, 2, 2 * p, "-")

# VanRaden Genomic Relationship Matrix (K = Z * Z' / (2 * sum(p * (1 - p))))
denom <- 2 * sum(p * (1 - p))
K <- (Z %*% t(Z)) / denom

# Export Kinship.csv for GWAS matrix input
write.csv(K, "Kinship.csv", row.names = TRUE)
cat(" -> Saved 'Kinship.csv' (", nrow(K), "x", ncol(K), "square matrix)\n")

# Generate and save Kinship Heatmap
png("Kinship_Heatmap.png", width = 900, height = 800, res = 130)
pheatmap(
  K,
  main = "Genomic Kinship Matrix (VanRaden Method)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100)
)
dev.off()
cat(" -> Saved 'Kinship_Heatmap.png'\n")

cat("\n==============================================================================\n")
cat(" MODULE 2 COMPLETE: PCA & Kinship files generated successfully!\n")
cat(" 'PCA.csv' now contains only 'Taxa' and 'PC1'-'PC5'.\n")
cat("==============================================================================\n")