# 0. Install missing packages automatically
required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "GGally")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)

# 1. Load the data
df <- read_excel("Phenotype data.xlsx")

# 2. Rename columns to standardized, clean names
df_clean <- df %>%
  rename(
    accession    = Accession,
    yield        = Yield,
    plant_height = `Plant height`,
    leaf_area    = `Leaf area`,
    seed_size    = `Seed size`
  )

# 3. View individual distributions via histograms & save plot
p_hist <- df_clean %>%
  pivot_longer(cols = yield:seed_size, names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white") +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw() +
  labs(title = "Distribution of Phenotypic Traits", x = "Value", y = "Count")

print(p_hist)
ggsave("phenotype_histograms.png", plot = p_hist, width = 8, height = 6, dpi = 300)

# 4. Generate pairwise correlation matrix & save plot
p_pairs <- ggpairs(df_clean, columns = 2:5, title = "Phenotypic Trait Relationships")

print(p_pairs)
ggsave("phenotype_ggpairs.png", plot = p_pairs, width = 8, height = 8, dpi = 300)

# 5. Outlier detection using standard 1.5 * IQR method
remove_iqr_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  
  x[x < lower_bound | x > upper_bound] <- NA
  return(x)
}

# Apply outlier replacement across trait columns and drop NA rows
df_filtered <- df_clean %>%
  mutate(across(c(yield, plant_height, leaf_area, seed_size), remove_iqr_outliers)) %>%
  drop_na()

# 6. Export filtered dataset to CSV
write.csv(df_filtered, "pheno_filtered.csv", row.names = FALSE)
