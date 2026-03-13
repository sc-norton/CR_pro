# === PROJECT: USC - CR Muscle Proteomics =============================
# === Script: PCA Analysis

#### PURPOSE: To assess CREvPLA and SURVvCTL via PCA

#### OUTPUTS: 


#### WORKFLOW:
#     0. Setup - load required packages, establish directories
#     1. Load Data

# --- 0: Setup --------------------------------------------
cat("\n>> 0 - Setup\n")
cat("\n Loading required packages...")

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, limma, writexl, ggVennDiagram,
               clusterProfiler, org.Hs.eg.db, enrichplot)

cat("\n Packages successfully loaded.\n")

# Identify directories
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
report_dir <- file.path(base_dir, "05_analysis", "b_results")

dir_tree <- list(
  CREvPLA = list(
    root         = file.path(report_dir, "CREvPLA"),
    limma        = file.path(report_dir, "CREvPLA", "limma_results"),
    plots        = file.path(report_dir, "CREvPLA", "plots"),
    volcano      = file.path(report_dir, "CREvPLA", "plots", "volcanos"),
    venn         = file.path(report_dir, "CREvPLA", "plots", "venns"),
    da_class     = file.path(report_dir, "CREvPLA", "DA_classification"),
    gsea         = file.path(report_dir, "CREvPLA", "GSEA")
  ),
  SURVvCTL = list(
    root         = file.path(report_dir, "SURVvCTL"),
    limma        = file.path(report_dir, "SURVvCTL", "limma_results"),
    plots        = file.path(report_dir, "SURVvCTL", "plots"),
    volcano      = file.path(report_dir, "SURVvCTL", "plots", "volcanos"),
    venn         = file.path(report_dir, "SURVvCTL", "plots", "venns"),
    da_class     = file.path(report_dir, "SURVvCTL", "DA_classification"),
    gsea         = file.path(report_dir, "SURVvCTL", "GSEA")
  ),
  FvM = list(
    root         = file.path(report_dir, "FvM"),
    limma        = file.path(report_dir, "FvM", "limma_results"),
    plots        = file.path(report_dir, "FvM", "plots"),
    volcano      = file.path(report_dir, "FvM", "plots", "volcanos"),
    venn         = file.path(report_dir, "FvM", "plots", "venns"),
    da_class     = file.path(report_dir, "FvM", "DA_classification"),
    gsea         = file.path(report_dir, "FvM", "GSEA")
  )
)

# Create all directories at once
invisible(lapply(dir_tree, function(ds_dirs) {
  lapply(ds_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}))

cat("\n Base directory: ", base_dir, "\n")

# --- 1: Load Data ------------------
cat("\n>> 1 - Load Imputed Data\n")

SURV <- readRDS(file.path(base_dir, 
                          "03_imputation", 
                          "c_data", 
                          "01_imputed_DAList_SURV.RDS"))
CRE <- readRDS(file.path(base_dir, 
                          "03_imputation", 
                          "c_data", 
                          "01_imputed_DAList_CRE.RDS"))


# --- 2: PCA function -----------

plot_pca <- function(mat, metadata, color_by = "group", title = NULL) {
  # mat: proteins x samples matrix (log-transformed, imputed)
  pca <- prcomp(t(mat), scale. = TRUE)
  
  pca_df <- as.data.frame(pca$x[, 1:2]) |>
    tibble::rownames_to_column("sample_id") |>
    left_join(metadata, by = "sample_id")
  
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    stat_ellipse(
      geom  = "polygon",
      alpha = 0.1,
      level = 0.95,
      type  = "t"
    ) +
    geom_point(size = 3) +
    labs(
      title = title,
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)")) +
    theme_minimal()
}

sum(is.na(SURV$data))

plot_pca(SURV$data, SURV$metadata, color_by = "cancer_time", title = "SURV v CTL")
plot_pca(SURV$data, SURV$metadata, color_by = "cancer", title = "SURV v CTL")
plot_pca(CRE$data, CRE$metadata, color_by = "supp_time", title = "CRE v PLA")
