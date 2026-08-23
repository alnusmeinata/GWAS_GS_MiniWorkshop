# 0. Automatic Package Installation & Loading
required_packages <- c("tidyverse", "ggplot2", "data.table", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(tidyverse)
library(ggplot2)
library(data.table)
library(scales)

# 1. Load Cross-Validation Results
results_file <- "gp_gblup_10fold_cv_results.csv"

if (!file.exists(results_file)) {
  stop(paste("File missing:", results_file, "- Please run the GBLUP cross-validation script first."))
}

df_results <- fread(results_file, header = TRUE, data.table = FALSE)

# 2. Summarize Metrics: Median and SD across 10 folds per Trait and SNP Count
df_summary <- df_results %>%
  filter(!is.na(corr)) %>%
  group_by(trait, n_loci) %>%
  summarise(
    median_corr = mean(corr, na.rm = TRUE),
    sd_corr = sd(corr, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Generate Line Plot with Median Points and Standard Deviation Error Bars
p <- ggplot(df_summary, aes(x = n_loci, y = median_corr, color = trait)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = median_corr - sd_corr, ymax = median_corr + sd_corr),
    width = max(df_summary$n_loci) * 0.015, # Dynamic whisker width
    alpha = 0.6,
    linewidth = 0.5
  ) +
  facet_wrap(~trait, scales = "free_y") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "GBLUP Genomic Prediction Accuracy across SNP Subsets",
    subtitle = "10-Fold Cross-Validation (Points = Median, Whiskers = ±1 SD)",
    x = "Number of Top SNPs (n_loci)",
    y = "Prediction Accuracy (Correlation, r)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.background = element_rect(fill = "gray95", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 11)
  )

# 4. Save High-Resolution Output Plot
output_plot <- "gblup_prediction_accuracy_plot.png"
ggsave(filename = output_plot, plot = p, width = 10, height = 7, units = "in", dpi = 300)

cat("Plot successfully generated and saved to:", output_plot, "\n")
