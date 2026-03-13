# === PROJECT: USC - CR Muscle Proteomics =============================
# === Script: DA Analysis

#### PURPOSE: To analyze the DA of proteins in the samples and save the results.

#### OUTPUTS: (a) analyzed datasets (.csv)


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
cat("\n>> 1 - Load Normalized Data\n")

SURV <- readRDS(file.path(base_dir, 
                          "02_normalization", 
                          "c_data", 
                          "01_normalized_DAList_SURV.RDS"))
CRE <- readRDS(file.path(base_dir, 
                          "02_normalization", 
                          "c_data", 
                          "01_normalized_DAList_CRE.RDS"))


# --- 2: Core analysis function -----------

run_limma_analysis <- function(dal,
                               contrast_specs,
                               group_col    = "supp_time",
                               block_col    = "pid",
                               label        = "analysis",     # used in output filenames
                               lfc_thresh   = 0.6,
                               pval_thresh  = 0.05,
                               dirs,
                               table_cols   = c("uniprot_id", "gene", "protein"),
                               title_col    = "gene") {
  
  cat(sprintf("\n>> Running limma: %s\n", label))
  
  # 1. Design & contrasts
  cat("  Step 1: Design & contrasts\n")
  dal <- add_design(dal, paste0("~ 0 + ", group_col, " + (1 | ", block_col, ")"))
  colnames(dal$design$design_matrix) <-
    gsub(paste0("^", group_col), "", colnames(dal$design$design_matrix))
  dal <- add_contrasts(dal, contrasts_vector = contrast_specs)
  
  design_mat   <- dal$design$design_matrix
  contrast_mat <- dal$design$contrast_matrix
  block_var    <- dal$metadata[[block_col]]
  mat          <- as.matrix(dal$data)
  
  # 2. Array weights
  cat("  Step 2: Array weights\n")
  aw       <- arrayWeights(mat, design_mat)
  names(aw) <- colnames(mat)
  cat(sprintf("    Weights: %.2f - %.2f (median %.2f)\n",
              min(aw), max(aw), median(aw)))
  write_csv(tibble(sample_id = names(aw), array_weight = round(aw, 4)),
            file.path(dirs$limma, paste0("array_weights_", label, ".csv")))
  
  # 3. Fit model
  cat("  Step 3: Fit model\n")
  dupcor <- duplicateCorrelation(mat, design_mat, block = block_var, weights = aw)
  cat(sprintf("    Within-subject correlation: %.3f\n", dupcor$consensus.correlation))
  
  fit <- lmFit(mat, design_mat,
               block       = block_var,
               correlation = dupcor$consensus.correlation,
               weights     = aw) |>
    contrasts.fit(contrast_mat) |>
    eBayes(robust = TRUE, trend = TRUE)
  
  dal$eBayes_fit             <- fit
  dal$eBayes_fit$correlation <- dupcor$consensus.correlation
  
  # 4. Extract results
  cat("  Step 4: Extract & write results\n")
  dal <- extract_DA_results(dal, pval_thresh = pval_thresh,
                            lfc_thresh = 0, adj_method = "BH")
  write_limma_tables(dal, output_dir = dirs$limma, overwrite = TRUE)
  file.rename(
    file.path(dirs$limma, "results.xlsx"),
    file.path(dirs$limma, paste0("results_", label, ".xlsx"))
  )
  write_limma_plots(dal, grouping_column = group_col, output_dir = dirs$plots,
                    table_columns = table_cols, title_column = title_col,
                    overwrite = TRUE)
  
  # 5. Read results back & compute pi-values
  per_contrast_dir <- file.path(dirs$limma, "per_contrast_results")
  results_list <- lapply(setNames(names(dal$results), names(dal$results)), function(cname) {
    read_csv(file.path(per_contrast_dir, paste0(cname, ".csv")),
             show_col_types = FALSE)
  })
  
  results_list <- lapply(results_list, function(df) {
    df$uniprot  <- df$uniprot_id %||% rownames(df)
    df$Pi.Val   <- abs(df$logFC) * (-log10(df$adj.P.Val))
    df$adj.Pi.Val <- 10^-df$Pi.Val
    df[, c("uniprot", setdiff(names(df), "uniprot"))]
  })
  
  # 6. DA classification (adj p, nominal p, pi-value)
  results_list <- classify_DA(results_list, lfc_thresh = lfc_thresh,
                              pval_thresh = pval_thresh)
  da_log       <- summarise_DA(results_list, lfc_thresh = lfc_thresh,
                               pval_thresh = pval_thresh)
  
  write.csv(da_log, 
            file.path(dirs$da_class, paste0("DA_criteria_", label, ".csv")),
            row.names = FALSE)
  
  # 7. Volcano plots
  cat("  Step 5: Volcano plots\n")
  volcanos <- make_volcano_plots(results_list, pval_thresh = pval_thresh,
                                 lfc_thresh = lfc_thresh)
  pdf(file.path(dirs$volcano, paste0("volcano_plots_", label, ".pdf")),
      width = 10, height = 8)
  invisible(lapply(volcanos, print))
  dev.off()
  
  list(dal          = dal,
       results_list = results_list,
       da_log       = da_log,
       volcanos     = volcanos)
  }

# ---- 3: DA helpers -----------

classify_DA <- function(results_list, lfc_thresh = 0.6, pval_thresh = 0.05) {
  lapply(results_list, function(df) {
    classify <- function(lfc, pval)
      case_when(lfc >  lfc_thresh & pval < pval_thresh ~ "UP",
                lfc < -lfc_thresh & pval < pval_thresh ~ "DOWN",
                TRUE                                    ~ "NDE")
    df$DE.adj <- classify(df$logFC, df$adj.P.Val)
    df$DE.nom <- classify(df$logFC, df$P.Value)
    df$DE.pi  <- classify(df$logFC, df$adj.Pi.Val)
    df
  })
}

summarise_DA <- function(results_list, lfc_thresh = 0.6, pval_thresh = 0.05) {
  cols    <- c(adj = "DE.adj", nom = "DE.nom", pi = "DE.pi")
  labels  <- c(adj = "BH-adj P-val", nom = "Nominal P-val", pi = "Pi-value")
  do.call(rbind, lapply(names(cols), function(key) {
    do.call(rbind, lapply(names(results_list), function(cname) {
      df <- results_list[[cname]]
      data.frame(comparison = cname, method = labels[[key]],
                 n_UP   = sum(df[[cols[[key]]]] == "UP"),
                 n_DOWN = sum(df[[cols[[key]]]] == "DOWN"),
                 stringsAsFactors = FALSE)
    }))
  }))
}

# --- Volcano Plots ----------------
make_volcano_plots <- function(results_list, pval_thresh = 0.05, lfc_thresh = 0.6) {
  lapply(names(results_list), function(name) {
    df <- results_list[[name]]
    ggplot(df, aes(logFC, -log10(adj.P.Val), color = DE.adj)) +
      geom_point() +
      geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed") +
      geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
      scale_color_manual(values = c(UP = "red1", DOWN = "turquoise3", NDE = "grey")) +
      coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) +
      annotate("text", x = -6, y = 7, hjust = 0, vjust = 1, size = 4,
               label = sprintf("UP: %d\nDOWN: %d",
                               sum(df$DE.adj == "UP"), sum(df$DE.adj == "DOWN"))) +
      labs(title    = name,
           subtitle = sprintf("BH-adj P-val < %.2f  |  |log2FC| > %.1f",
                              pval_thresh, lfc_thresh),
           x = expression("log"[2]*"FC"),
           y = expression("-log"[10]*"p-value"),
           color = "Direction") +
      theme_classic() +
      theme(text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
  }) |> setNames(names(results_list))
}

# --- Venn diagram helpers ----------

make_venn <- function(results_list, comparisons, direction_col = "DE.adj",
                      report_dir, label = "venn") {
  # comparisons: named list, e.g. list("CRE UP" = "Training_CRE", ...)
  # Each element: contrast name -> direction ("UP"/"DOWN")
  venn_sets <- lapply(names(comparisons), function(set_name) {
    spec     <- comparisons[[set_name]]   # list(contrast = ..., direction = ...)
    df       <- results_list[[spec$contrast]]
    df$uniprot[df[[direction_col]] == spec$direction]
  }) |> setNames(names(comparisons))
  
  p <- ggVennDiagram(venn_sets) +
    scale_fill_gradient(low = "white", high = "darkblue") +
    scale_x_continuous(expand = expansion(mult = 0.3)) +
    scale_y_continuous(expand = expansion(mult = 0.3)) +
    theme(plot.margin = margin(20, 20, 20, 20, "mm"))
  
  pdf(file.path(report_dir, paste0(label, ".pdf")), width = 10, height = 8)
  print(p + theme(plot.margin = margin(40,40,40,40, "mm")))
  dev.off()
  
  invisible(venn_sets)
}

# ---- Run Analysis --------

# Build contrasts for each dataset
contrast_specs_cre <- c(
  "Baseline_CREvPLA = CRE_T1 - PLA_T1",
  "Training_CREvPLA = CRE_T2 - PLA_T2",
  "Training_CRE     = CRE_T2 - CRE_T1",
  "Training_PLA     = PLA_T2 - PLA_T1",
  "Interaction_supp      = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

contrast_specs_cre_sex <- c(
  "Baseline_FvM = Female_T1 - Male_T1",
  "Training_FvM = Female_T2 - Male_T2",
  "Training_F   = Female_T2 - Female_T1",
  "Training_M   = Male_T2 - Male_T1",
  "Interaction_sex  = (Female_T2 - Female_T1) - (Male_T2 - Male_T1)"
)

contrast_specs_surv <- c(
  "Baseline_SURVvCTL = SURV_T1 - CTL_T1",
  "Training_SURVvCTL = SURV_T2 - CTL_T1",
  "Training_SURV     = SURV_T2 - SURV_T1"
)

# Run analysis on both datasets
out_cre  <- run_limma_analysis(CRE, contrast_specs_cre,
                               group_col = "supp_time",
                               label     = "CRE",
                               dirs      = dir_tree$CREvPLA)

out_cre_sex <- run_limma_analysis(CRE, contrast_specs_cre_sex,
                                  group_col = "sex_time",
                                  label     = "SEX",
                                  dirs      = dir_tree$FvM)

out_surv <- run_limma_analysis(SURV,  contrast_specs_surv,
                               group_col = "cancer_time",
                               label     = "SURV",
                               dirs      = dir_tree$SURVvCTL)


# Combine results for joint export / DA bar chart
results_combined <- c(out_cre$results_list, 
                      out_surv$results_list, 
                      out_cre_sex$results_list)
results_combined <- lapply(results_combined, function(x) {
  df <- as.data.frame(x)
  df[ , !names(df) %in% c("gene_symbol"), drop = FALSE]
})
da_combined      <- rbind(out_cre$da_log, out_surv$da_log, out_cre_sex$da_log)

write_xlsx(results_combined,
           file.path(report_dir, "limma_results_all_muscle.xlsx"))
write.csv(da_combined,
          file.path(report_dir, "DA_methods_combined.csv"), row.names = FALSE)

# --- DA bar chart across both datasets----------
plot_DA <- da_combined %>%
  pivot_longer(c(n_UP, n_DOWN), names_to = "direction", values_to = "count") %>%
  mutate(direction = gsub("n_", "", direction))

plot_DA2 <- da_combined %>%
  filter(method == "BH-adj P-val") %>%
  pivot_longer(c(n_UP, n_DOWN), names_to = "direction", values_to = "count") %>%
  mutate(direction = gsub("n_", "", direction))

plot_DA3 <- da_combined %>%
  filter(method == "Nominal P-val") %>%
  pivot_longer(c(n_UP, n_DOWN), names_to = "direction", values_to = "count") %>%
  mutate(direction = gsub("n_", "", direction))

p_DA <- ggplot(plot_DA, aes(x = comparison, y = count, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~method, scales = "free_x") +
  scale_fill_manual(values = c(UP = "#E96A2C", DOWN = "#0D234A")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Comparison", y = "Number of Proteins", fill = "Direction",
       title    = "Differentially Abundant Proteins by Criteria and Comparison",
       subtitle = sprintf("|log2FC| > 0.6, sig level 0.05"))

p_DA2 <- ggplot(plot_DA2, aes(x = comparison, y = count, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(UP = "#E96A2C", DOWN = "#0D234A")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Comparison", y = "Number of Proteins", fill = "Direction",
       title    = "Differentially Abundant Proteins by Comparison",
       subtitle = sprintf("|log2FC| > 0.6, BH-adj P < 0.05"))

p_DA3 <- ggplot(plot_DA3, aes(x = comparison, y = count, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(UP = "#E96A2C", DOWN = "#0D234A")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Comparison", y = "Number of Proteins", fill = "Direction",
       title    = "Differentially Abundant Proteins by Comparison",
       subtitle = sprintf("|log2FC| > 0.6, Nominal P < 0.05"))

ggsave(file.path(report_dir, "DA_methods_combined.png"), 
       p_DA,
       width = 14, height = 6, dpi = 300)

ggsave(file.path(report_dir, "DA_count_all_comparisons_BH.png"),
       p_DA2,
       width = 14, height = 6, dpi = 300)

ggsave(file.path(report_dir, "DA_count_all_comparisons_nom.png"),
       p_DA3,
       width = 14, height = 6, dpi = 300)

# Venn diagrams
make_venn(results_combined,
          comparisons = list(
            "CRE UP"   = list(contrast = "Training_CRE",  direction = "UP"),
            "CRE DOWN" = list(contrast = "Training_CRE",  direction = "DOWN"),
            "PLA UP"   = list(contrast = "Training_PLA",  direction = "UP"),
            "PLA DOWN" = list(contrast = "Training_PLA",  direction = "DOWN")
          ),
          report_dir = dir_tree$CREvPLA$venn, label = "venn_CREvPLA_proteins")

make_venn(results_combined,
          comparisons = list(
            "SURV_T1 UP"   = list(contrast = "Baseline_SURVvCTL", direction = "UP"),
            "SURV_T1 DOWN" = list(contrast = "Baseline_SURVvCTL", direction = "DOWN"),
            "SURV_T2 UP"   = list(contrast = "Training_SURVvCTL", direction = "UP"),
            "SURV_T2 DOWN" = list(contrast = "Training_SURVvCTL", direction = "DOWN")
          ),
          report_dir = dir_tree$SURVvCTL$venn, label = "venn_SURVvCTL_proteins")

make_venn(results_combined, 
          comparisons = list(
            "Female UP"   = list(contrast = "Training_F", direction = "UP"),
            "Female DOWN" = list(contrast = "Training_F", direction = "DOWN"),
            "Male UP"     = list(contrast = "Training_M", direction = "UP"),
            "Male DOWN"   = list(contrast = "Training_M", direction = "DOWN")
          ),
          report_dir = dir_tree$FvM$venn, label = "venn_FvM_proteins")


# --- Export UP/DOWN protein lists per comparison ----------
cat("\n>> Exporting UP/DOWN protein lists per comparison\n")

export_dir <- file.path(report_dir, "DA_protein_lists")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

invisible(lapply(names(results_combined), function(cname) {
  df <- results_combined[[cname]]
  
  up_df   <- df[df$DE.adj == "UP",   ]
  down_df <- df[df$DE.adj == "DOWN", ]
  
  write.csv(up_df,
            file.path(export_dir, paste0(cname, "_UP.csv")),
            row.names = FALSE)
  write.csv(down_df,
            file.path(export_dir, paste0(cname, "_DOWN.csv")),
            row.names = FALSE)
  
  cat(sprintf("  %s -> UP: %d, DOWN: %d\n", cname, nrow(up_df), nrow(down_df)))
}))

cat("\n  Protein lists saved to:", export_dir, "\n")

# --- Compare Concordant/Discordant Proteins --------

# Select F and M RT data
F_RT <- results_combined$Training_F %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_F      = logFC,
                P.Value_F    = P.Value,
                adj.P.Val_F  = adj.P.Val,
                adj.Pi.Val_F = adj.Pi.Val) %>%
  mutate(signif_F = adj.Pi.Val_F < 0.05)
M_RT <- results_combined$Training_M %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_M      = logFC,
                P.Value_M    = P.Value,
                adj.P.Val_M  = adj.P.Val,
                adj.Pi.Val_M = adj.Pi.Val) %>%
  mutate(signif_M = adj.Pi.Val_M < 0.05)

# Join data sets
cord <- inner_join(F_RT, M_RT, by = "gene")

# Compare agreement
cord <- cord %>%
  mutate(
    direction_F = ifelse(logFC_F > 0, "up", "down"),
    direction_M = ifelse(logFC_M > 0, "up", "down"),
    agreement   = ifelse(direction_F == direction_M, "Concordant", "Discordant"),
  ) %>%
  mutate(
    sig_group = case_when(
      signif_M & signif_F  ~ "Both significant",
      signif_M & !signif_F ~ "Male only",
      !signif_M & signif_F ~ "Female only",
      TRUE                 ~ "Neither"
    )
  )

# Scatter plot
ggplot(cord, aes(x = logFC_M, y = logFC_F, color = sig_group)) +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf,  
           fill = "green3", alpha = 0.15) +  # top-right:  both up
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0,    
           fill = "green3", alpha = 0.15) +  # bottom-left: both down
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = -Inf, ymax = 0,    
           fill = "steelblue3", alpha = 0.15) +  # bottom-right: disagree
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax = Inf,  
           fill = "hotpink2", alpha = 0.15) +  # top-left:    disagree
  annotate("text", x = Inf,  y = Inf,  label = "UP Both",   
           hjust = 1.1, vjust = 1.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = -Inf, y = -Inf, label = "DOWN Both",  
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = Inf,  y = -Inf, label = "UP Male \n DOWN Female",  
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "steelblue3", fontface = "italic") +
  annotate("text", x = -Inf, y = Inf,  label = "UP Female \n DOWN Male",  
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "hotpink", fontface = "italic") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Both significant" = "purple", 
                                "Male only" = "blue",
                                "Female only" = "hotpink")) +
  labs(
    x     = "Male RT Response",
    y     = "Female RT Response",
    title = "Concordance: Sex Comparison",
    color = "Signifant DE"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid       = element_blank()
  )

# Select CRE and PLA RT data
CRE_RT <- results_combined$Training_CRE %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_CRE      = logFC,
                P.Value_CRE    = P.Value,
                adj.P.Val_CRE  = adj.P.Val,
                adj.Pi.Val_CRE = adj.Pi.Val) %>%
  mutate(signif_CRE = adj.Pi.Val_CRE < 0.05)
PLA_RT <- results_combined$Training_PLA %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_PLA      = logFC,
                P.Value_PLA    = P.Value,
                adj.P.Val_PLA  = adj.P.Val,
                adj.Pi.Val_PLA = adj.Pi.Val) %>%
  mutate(signif_PLA = adj.Pi.Val_PLA < 0.05)

# Join data sets
cord2 <- inner_join(CRE_RT, PLA_RT, by = "gene")

# Compare agreement
cord2 <- cord2 %>%
  mutate(
    direction_CRE = ifelse(logFC_CRE > 0, "up", "down"),
    direction_PLA = ifelse(logFC_PLA > 0, "up", "down"),
    agreement   = ifelse(direction_CRE == direction_PLA, "Concordant", "Discordant"),
  ) %>%
  mutate(
    sig_group = case_when(
      signif_PLA & signif_CRE  ~ "Both significant",
      signif_PLA & !signif_CRE ~ "PLA only",
      !signif_PLA & signif_CRE ~ "CRE only",
      TRUE                     ~ "Neither"
    )
  )

# Scatter plot
ggplot(cord2, aes(x = logFC_PLA, y = logFC_CRE, color = sig_group)) +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf,  
           fill = "green3", alpha = 0.15) +  # top-right:  both up
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0,    
           fill = "green3", alpha = 0.15) +  # bottom-left: both down
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = -Inf, ymax = 0,    
           fill = "steelblue3", alpha = 0.15) +  # bottom-right: disagree
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax = Inf,  
           fill = "red", alpha = 0.15) +  # top-left:    disagree
  annotate("text", x = Inf,  y = Inf,  label = "Both UP",   
           hjust = 1.1, vjust = 1.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = -Inf, y = -Inf, label = "Both DOWN",  
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = Inf,  y = -Inf, label = "UP PLA \n DOWN CRE",  
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "steelblue3", fontface = "italic") +
  annotate("text", x = -Inf, y = Inf,  label = "UP CRE \n DOWN PLA",  
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "red", fontface = "italic") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Both significant" = "purple", 
                                "PLA only" = "blue",
                                "CRE only" = "red")) +
  labs(
    x     = "PLA RT Response",
    y     = "CRE RT Response",
    title = "Concordance: Supplement Comparison",
    color = "Signifant DE (∏)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid       = element_blank()
  )

# Select SURV v CTL data
SURV_T1 <- results_combined$Baseline_SURVvCTL %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_T1      = logFC,
                P.Value_T1    = P.Value,
                adj.P.Val_T1  = adj.P.Val,
                adj.Pi.Val_T1 = adj.Pi.Val) %>%
  mutate(signif_T1 = adj.Pi.Val_T1 < 0.05)
SURV_T2 <- results_combined$Training_SURVvCTL %>%
  dplyr::select(uniprot_id, gene, logFC, P.Value, adj.P.Val, adj.Pi.Val) %>%
  dplyr::rename(logFC_T2      = logFC,
                P.Value_T2    = P.Value,
                adj.P.Val_T2  = adj.P.Val,
                adj.Pi.Val_T2 = adj.Pi.Val) %>%
  mutate(signif_T2 = adj.Pi.Val_T2 < 0.05)

# Join data sets
cord3 <- inner_join(SURV_T1, SURV_T2, by = "gene")

# Compare agreement
cord3 <- cord3 %>%
  mutate(
    direction_T1 = ifelse(logFC_T1 > 0, "up", "down"),
    direction_T2 = ifelse(logFC_T2 > 0, "up", "down"),
    agreement   = ifelse(direction_T1 == direction_T2, "Concordant", "Discordant"),
  ) %>%
  mutate(
    sig_group = case_when(
      signif_T1 & signif_T2  ~ "Both significant",
      signif_T1 & !signif_T2 ~ "BL only",
      !signif_T1 & signif_T2 ~ "RT only",
      TRUE                   ~ "Neither"
    )
  )

# Scatter plot
ggplot(cord3, aes(x = logFC_T1, y = logFC_T2, color = sig_group)) +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf,  
           fill = "green3", alpha = 0.15) +  # top-right:  both up
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0,    
           fill = "green3", alpha = 0.15) +  # bottom-left: both down
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = -Inf, ymax = 0,    
           fill = "steelblue3", alpha = 0.15) +  # bottom-right: disagree
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax = Inf,  
           fill = "red", alpha = 0.15) +    # top-left:    disagree
  annotate("text", x = Inf,  y = Inf,  label = "UP T1 & T2",   
           hjust = 1.1, vjust = 1.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = -Inf, y = -Inf, label = "DOWN T1 & T2",  
           hjust = -0.1, vjust = -0.5, size = 3.5, color = "green4", fontface = "italic") +
  annotate("text", x = Inf,  y = -Inf, label = "UP T1 \n DOWN T2",  
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "steelblue3", fontface = "italic") +
  annotate("text", x = -Inf, y = Inf,  label = "UP T2 \n DOWN T1",  
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "red", fontface = "italic") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Both significant" = "purple", 
                                "BL only" = "blue",
                                "RT only" = "red")) +
  labs(
    x     = "SURV_T1 v CTL (logFC)",
    y     = "SURV_T2 v CTL (logFC)",
    title = "SURV v CTL",
    color = "Signifant DE (∏)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid       = element_blank()
  )
