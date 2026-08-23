# 0. Automatic Package Installation & Loading
required_packages <- c("tidyverse", "patchwork", "ggplot2", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)

# 1. Trait Initialization (Auto-detect from pheno_filtered.csv if not loaded)
if (!exists("traits")) {
  filename_pheno <- "pheno_filtered.csv"
  pheno_raw <- fread(filename_pheno, header = TRUE, data.table = FALSE)
  if ("accession" %in% colnames(pheno_raw)) {
    pheno_raw <- pheno_raw %>% dplyr::rename(ID = accession)
  }
  traits <- setdiff(colnames(pheno_raw), "ID")
}

folder <- "GWAS_QK"

# 2. Plotting Loop strictly for Q+K Model Results
for (trait in traits) {
  cat("Generating plots for trait:", trait, "\n")
  
  filename_inter <- file.path(folder, paste0("gwas_", trait, "_intermediate result.csv"))
  filename_final <- file.path(folder, paste0("gwas_", trait, "_Final result.csv"))
  
  if (!file.exists(filename_inter)) {
    message(paste("Intermediate file missing:", filename_inter, "- Skipping plot."))
    next
  }
  
  # Read intermediate results
  a <- read.csv(filename_inter, check.names = FALSE) %>% as_tibble()
  
  # Dynamic matching for -log10(P) column name
  p_col <- grep("-log10\\(P\\)", colnames(a), value = TRUE)
  if (length(p_col) > 0) {
    a <- a %>% dplyr::rename(P = !!p_col[1])
  } else if (!"P" %in% colnames(a)) {
    message(paste("Could not locate -log10(P) column in", filename_inter))
    next
  }
  
  # Ensure standard chromosome column naming
  if ("Chromosome" %in% colnames(a)) {
    a <- a %>% dplyr::rename(CHR = Chromosome)
  }
  
  a <- a %>%
    mutate(X = row_number()) %>%
    mutate(CHR_group = as.integer(as.factor(CHR)) %% 2)
  
  # Process final significant markers (LOD Scores)
  has_final <- FALSE
  if (file.exists(filename_final) && file.info(filename_final)$size > 0) {
    a_fin <- tryCatch(
      read.csv(filename_final, check.names = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(a_fin) && nrow(a_fin) > 0) {
      rs_col <- grep("RS|Marker|SNP", colnames(a_fin), value = TRUE, ignore.case = TRUE)[1]
      lod_col <- grep("LOD", colnames(a_fin), value = TRUE, ignore.case = TRUE)[1]
      
      if (!is.na(rs_col) && !is.na(lod_col)) {
        a_fin <- a_fin %>%
          dplyr::select(!!rs_col, !!lod_col) %>%
          dplyr::rename(Marker = !!rs_col, LOD.score = !!lod_col)
        
        a <- left_join(a, a_fin, by = "Marker")
        has_final <- TRUE
      }
    }
  }
  
  # Dynamic axis scaling
  log10p_upper <- max(c(7, a$P), na.rm = TRUE) * 1.05
  lod_upper <- if (has_final && any(!is.na(a$LOD.score))) max(c(6, a$LOD.score), na.rm = TRUE) * 1.05 else 6
  ratio <- lod_upper / log10p_upper
  
  # --- QQ Plot ---
  upper_qq <- max(-log10(ppoints(nrow(a))))
  ylim_qq <- max(upper_qq, a$P, na.rm = TRUE) * 1.05
  
  qq <- tibble(P = sort(a$P, decreasing = TRUE)) %>%
    filter(P > 0) %>%
    mutate(exp = -log10(ppoints(n()))) %>%
    ggplot(aes(exp, P)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_point(size = 1.2, color = "navyblue") +
    scale_y_continuous(limits = c(0, ylim_qq)) +
    scale_x_continuous(limits = c(0, ylim_qq)) +
    coord_fixed() +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # --- Manhattan Plot ---
  axis_set <- a %>%
    group_by(CHR) %>%
    summarise(center = mean(X), .groups = 'drop') %>%
    mutate(CHR = as.character(CHR))
  
  p3 <- ggplot(a, aes(X, P, color = factor(CHR_group))) +
    geom_point(size = 1) +
    {
      if (has_final && any(!is.na(a$LOD.score))) {
        list(
          geom_point(data = filter(a, !is.na(LOD.score)), aes(X, LOD.score / ratio), color = "magenta", size = 2),
          geom_segment(
            data = filter(a, !is.na(LOD.score)),
            aes(x = X, xend = X, y = P, yend = LOD.score / ratio),
            color = "magenta",
            linewidth = 0.4,
            linetype = "dashed"
          )
        )
      }
    } +
    geom_hline(yintercept = 3, color = "magenta", linetype = 'dotdash', linewidth = 0.8) +
    scale_x_continuous(
      labels = axis_set$CHR,
      breaks = axis_set$center,
      name = "Scaffold / Chromosome"
    ) +
    scale_y_continuous(
      sec.axis = sec_axis(~. * ratio, name = "LOD Score"),
      limits = c(0, log10p_upper)
    ) +
    scale_color_manual(values = c('darkorange','darkblue')) +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.text.y.right = element_text(color = "magenta"),
      axis.title.y.right = element_text(color = "magenta"),
      axis.line.y.right = element_line(color = "magenta"),
      axis.ticks.y.right = element_line(color = "magenta"),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold")
    ) +
    ggtitle(paste("Q+K Model (FASTmrEMMA):", trait))
  
  # --- Combine Layout with Patchwork ---
  ylab_left <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Observed~-log[10](P)", angle = 90, parse = TRUE) +
    coord_cartesian(clip = "off") + theme_void()
  
  ylab_center <- ggplot() +
    annotate("text", x = 1, y = 1, label = "-log[10](P)", angle = 90, parse = TRUE) +
    coord_cartesian(clip = "off") + theme_void()
  
  manhattan_combined <- ylab_left + qq + ylab_center + p3 +
    patchwork::plot_layout(widths = c(1, 6, 1, 18))
  
  # Save plot inside the GWAS_QK folder
  out_file <- file.path(folder, paste0("manhattan_", trait, "_GWAS_QK.png"))
  ggsave(filename = out_file, plot = manhattan_combined, width = 11, height = 3.5, units = "in", dpi = 300)
}
