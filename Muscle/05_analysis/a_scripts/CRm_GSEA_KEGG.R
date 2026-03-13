# =============================================================================
# GSEA & ORA PIPELINE — Modular Version
# =============================================================================
# WORKFLOW:
#   1. Run run_gsea_pipeline() once per dataset  → saves .rds to disk
#   2. Use plot_ora() / plot_gsea() anytime after → loads from .rds, no re-running
# =============================================================================

# --- Load packages ------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
  ggrepel, limma, writexl, ggVennDiagram,
  clusterProfiler, org.Hs.eg.db, enrichplot, msigdbr, enrichR
)
cat("\nPackages successfully loaded.\n")

# --- Directories --------------------------------------------------------------
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
report_dir <- file.path(base_dir, "05_analysis", "b_results")
organism   <- org.Hs.eg.db

# --- Load limma results -------------------------------------------------------
results_list <- list(
  Baseline_CREvPLA  = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 1),
  Training_CREvPLA  = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 2),
  Training_CRE      = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 3),
  Training_PLA      = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 4),
  Interaction_supp  = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 5),
  Baseline_SURVvCTL = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 6),
  Training_SURVvCTL = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 7),
  Training_SURV     = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 8),
  Baseline_FvM      = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 9),
  Training_FvM      = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 10),
  Training_F        = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 11),
  Training_M        = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 12),
  Interaction_sex   = read_excel(file.path(report_dir, "limma_results_all_muscle.xlsx"), sheet = 13)
)

diff_list <- list(
  Baseline_SURVvCTL_UP = read.csv(file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_UP.csv")),
  Baseline_SURVvCTL_DOWN = read.csv(file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_DOWN.csv")),
  Training_SURVvCTL_UP = read.csv(file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_UP.csv")),
  Training_SURVvCTL_DOWN = read.csv(file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_DOWN.csv"))
)

folder_map <- c(
  Baseline_SURVvCTL = "SURVvCTL",
  Training_SURVvCTL = "SURVvCTL",
  UP_BL_only = "SURVvCTL",
  UP_RT_only = "SURVvCTL",
  DOWN_BL_only = "SURVvCTL",
  DOWN_RT_only = "SURVvCTL",
  Training_SURV     = "SURVvCTL",
  Baseline_CREvPLA  = "CREvPLA",
  Training_CREvPLA  = "CREvPLA",
  Training_CRE      = "CREvPLA",
  Training_PLA      = "CREvPLA",
  Interaction_supp  = "CREvPLA",
  Baseline_FvM      = "FvM",
  Training_FvM      = "FvM",
  Training_F        = "FvM",
  Training_M        = "FvM",
  Interaction_sex   = "FvM"
)

# =============================================================================
# STEP 0 — Identify Unique DE Pros SURV_T1 v SURV_T2
# =============================================================================

# UP proteins
df1 <- diff_list$Baseline_SURVvCTL_UP
df2 <- diff_list$Training_SURVvCTL_UP

df1$source <- "BL"
df2$source <- "RT"

combined <- bind_rows(df1, df2) %>%
  group_by(uniprot) %>%
  mutate(
    presence = case_when(
      "BL" %in% source & "RT" %in% source ~ "Both",
      "BL" %in% source                    ~ "BL only",
      TRUE                                ~ "RT only"
    )
  ) %>%
  ungroup()

UP_BL_only <- combined %>%
  filter(presence == "BL only")

UP_RT_only <- combined %>%
  filter(presence == "RT only")

# DOWN proteins
df3 <- diff_list$Baseline_SURVvCTL_DOWN
df4 <- diff_list$Training_SURVvCTL_DOWN

df3$source <- "BL"
df4$source <- "RT"

combined <- bind_rows(df3, df4) %>%
  group_by(uniprot) %>%
  mutate(
    presence = case_when(
      "BL" %in% source & "RT" %in% source ~ "Both",
      "BL" %in% source                    ~ "BL only",
      TRUE                                ~ "RT only"
    )
  ) %>%
  ungroup()

DOWN_BL_only <- combined %>%
  filter(presence == "BL only")

DOWN_RT_only <- combined %>%
  filter(presence == "RT only")

# Combine as list
diff_only_list <- list(
  UP_BL_only = UP_BL_only,
  UP_RT_only = UP_RT_only,
  DOWN_BL_only = DOWN_BL_only,
  DOWN_RT_only = DOWN_RT_only
)

# =============================================================================
# STEP 1 — PIPELINE FUNCTION (run once per dataset, saves .rds)
# =============================================================================

run_gsea_pipeline <- function(df, label, folder) {
  
  cat(sprintf("\n>> Running GSEA pipeline for: %s\n", label))
  
  # --- Rank metric ---
  df$rank                   <- df$logFC
  original_gene_list        <- df$rank
  names(original_gene_list) <- df$Ensembl
  gene_list                 <- na.omit(original_gene_list)
  gene_list                 <- sort(gene_list, decreasing = TRUE)
  
  # --- Hallmark GSEA ---
  hallmark <- msigdbr(collection = "H") %>% dplyr::select(gs_name, ensembl_gene)
  h2       <- msigdbr(collection = "H") %>% dplyr::select(gs_name, gs_description)
  
  gsea <- GSEA(
    geneList  = gene_list,
    TERM2GENE = hallmark,
    TERM2NAME = h2
  )
  
  # --- GO GSEA (BP) ---
  gse_bp <- gseGO(
    geneList      = gene_list,
    ont           = "BP",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gse_bp_clean        <- gse_bp
  gse_bp_clean@result <- gse_bp@result %>% filter(!is.na(NES), !is.na(pvalue))
  
  # --- GO GSEA (CC) ---
  gse_cc <- gseGO(
    geneList      = gene_list,
    ont           = "CC",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gse_cc_clean        <- gse_cc
  gse_cc_clean@result <- gse_cc@result %>% filter(!is.na(NES), !is.na(pvalue))
  
  # --- GO GSEA (MF) ---
  gse_mf <- gseGO(
    geneList      = gene_list,
    ont           = "MF",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gse_mf_clean        <- gse_mf
  gse_mf_clean@result <- gse_mf@result %>% filter(!is.na(NES), !is.na(pvalue))
  
  # --- KEGG GSEA ---
  ids       <- bitr(names(original_gene_list), fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = organism)
  dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]
  df2       <- df[df$Ensembl %in% dedup_ids$ENSEMBL, ]
  df2$Y     <- dedup_ids$ENTREZID[match(df2$Ensembl, dedup_ids$ENSEMBL)]
  
  kegg_gene_list        <- df2$logFC
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list        <- na.omit(kegg_gene_list)
  kegg_gene_list        <- sort(kegg_gene_list, decreasing = TRUE)
  
  kk2 <- gseKEGG(
    geneList      = kegg_gene_list,
    organism      = "hsa",
    eps           = 0,
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    seed          = 26,
    pAdjustMethod = "BH",
    keyType       = "ncbi-geneid"
  )
  
  # --- GO ORA ---
  sig_genes  <- df %>% filter(adj.P.Val < 0.05)
  sig_entrez <- bitr(sig_genes$Ensembl, fromType = "ENSEMBL",
                     toType = "ENTREZID", OrgDb = organism)
  # bg_entrez  <- bitr(df$Ensembl, fromType = "ENSEMBL",
  #                    toType = "ENTREZID", OrgDb = organism)
  
  run_go_ora <- function(ont_type) {
    tryCatch(
      enrichGO(
        gene          = sig_entrez$ENTREZID,
        # universe      = bg_entrez$ENTREZID, # removed to have full GO background
        OrgDb         = organism,
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat(sprintf("   GO ORA (%s) failed for %s: %s\n", ont_type, label, e$message))
        NULL
      }
    )
  }
  
  ora_bp <- run_go_ora("BP")
  ora_cc <- run_go_ora("CC")
  ora_mf <- run_go_ora("MF")
  
  # --- Save .rds ---
  results <- list(
    label        = label,
    gsea         = gsea,
    gse_bp_clean = gse_bp_clean,
    gse_cc_clean = gse_cc_clean,
    gse_mf_clean = gse_mf_clean,
    kk2          = kk2,
    ora_bp       = ora_bp,
    ora_cc       = ora_cc,
    ora_mf       = ora_mf
  )
  
  rds_dir <- file.path(report_dir, folder, "rds")
  if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
  saveRDS(results, file.path(rds_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}


# =============================================================================
# STEP 2 — LOAD HELPER (used internally by plot functions)
# =============================================================================

load_results <- function(label, folder) {
  path <- file.path(report_dir, folder, "rds", sprintf("gsea_results_%s.rds", label))
  if (!file.exists(path)) {
    stop(sprintf(
      "No saved results found for '%s'. Run run_gsea_pipeline() first.", label
    ))
  }
  readRDS(path)
}


# =============================================================================
# STEP 3 — PLOT FUNCTIONS (call anytime after pipeline has been run)
# =============================================================================

# -----------------------------------------------------------------------------
# plot_ora()
#
# Arguments:
#   label     — dataset name, e.g. "Training_CRE"
#   folder    — subfolder, e.g. "CREvPLA"  (use folder_map[label] if unsure)
#   ont       — ontology: "BP", "CC", or "MF"
#   type      — plot type: "dot" or "bar"
#   show_n    — number of terms to show (default 20)
#   save      — if TRUE, saves PDF to ORA subfolder
# -----------------------------------------------------------------------------

plot_ora <- function(label, folder, ont = "BP", type = "dot", show_n = 20, save = FALSE) {
  
  res <- load_results(label, folder)
  
  obj <- switch(ont,
                BP = res$ora_bp,
                CC = res$ora_cc,
                MF = res$ora_mf,
                stop("ont must be 'BP', 'CC', or 'MF'")
  )
  
  plot_title <- sprintf("GO %s ORA - %s", ont, label)
  
  if (is.null(obj) || nrow(obj@result %>% filter(p.adjust < 0.05)) == 0) {
    
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("No significant results\n", label)) +
      theme_void() +
      labs(title = plot_title)
    
  } else if (type == "dot") {
    
    p <- dotplot(obj, showCategory = show_n) +
      labs(title = plot_title) +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 8))
    
  } else if (type == "bar") {
    
    p <- barplot(obj, showCategory = show_n) +
      labs(title = plot_title) +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 8))
    
  } else {
    stop("type must be 'dot' or 'bar'")
  }
  
  if (save) {
    ora_dir <- file.path(report_dir, folder, "ORA")
    if (!dir.exists(ora_dir)) dir.create(ora_dir, recursive = TRUE)
    out_path <- file.path(ora_dir,
                          sprintf("ORA_GO_%s_%s_%s.pdf", ont, type, label))
    ggsave(out_path, p, width = 10, height = 8)
    cat(sprintf("Saved: %s\n", out_path))
  }
  
  p
}


# -----------------------------------------------------------------------------
# plot_gsea()
#
# Arguments:
#   label     — dataset name, e.g. "Training_CRE"
#   folder    — subfolder, e.g. "CREvPLA"
#   type      — "hallmark", "BP", "CC", "MF", "kegg", or "kegg_ridge"
#   show_n    — number of terms to show (default 10)
#   save      — if TRUE, saves PDF to GSEA subfolder
# -----------------------------------------------------------------------------

plot_gsea <- function(label, folder, type = "hallmark", show_n = 10, save = FALSE) {
  
  res <- load_results(label, folder)
  
  obj <- switch(type,
                hallmark   = res$gsea,
                BP         = res$gse_bp_clean,
                CC         = res$gse_cc_clean,
                MF         = res$gse_mf_clean,
                kegg       = res$kk2,
                kegg_ridge = res$kk2,
                stop("type must be 'hallmark', 'BP', 'CC', 'MF', 'kegg', or 'kegg_ridge'")
  )
  
  plot_title <- sprintf("%s GSEA - %s", type, label)
  n_sig      <- nrow(obj@result %>% filter(p.adjust < 0.05))
  
  if (n_sig == 0) {
    
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("No significant results\n", label)) +
      theme_void() +
      labs(title = plot_title)
    
  } else if (type == "kegg_ridge") {
    
    p <- tryCatch(
      ridgeplot(obj) +
        labs(x = "Enrichment Distribution", title = plot_title),
      error = function(e) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste("Ridgeplot error\n", label)) +
          theme_void()
      }
    )
    
  } else {
    
    p <- dotplot(obj, showCategory = show_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = plot_title) +
      theme(axis.text.y = element_text(size = 8))
    
  }
  
  if (save) {
    gsea_dir <- file.path(report_dir, folder, "GSEA")
    if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
    out_path <- file.path(gsea_dir, sprintf("GSEA_%s_%s.pdf", type, label))
    ggsave(out_path, p, width = 12, height = 8)
    cat(sprintf("Saved: %s\n", out_path))
  }
  
  p
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Run pipeline for one dataset (EXAMPLE) ---
# run_gsea_pipeline(results_list$Training_CRE, "Training_CRE", "CREvPLA")

# --- Run pipeline for ALL datasets ---
lapply(names(results_list), function(nm) {
  run_gsea_pipeline(df = results_list[[nm]], label = nm, folder = folder_map[[nm]])
})

lapply(names(diff_only_list), function(nm) {
  run_gsea_pipeline(df = diff_only_list[[nm]], label = nm, folder = folder_map[[nm]])
})

# --- Plot ORA ---
p <- plot_ora("UP_BL_only", "SURVvCTL",  ont = "BP", type = "bar")
p

plot_ora("Training_CRE",      "CREvPLA",  ont = "BP", type = "bar")
plot_ora("Training_CRE",      "CREvPLA",  ont = "MF", type = "bar", show_n = 10)
plot_ora("Baseline_FvM",      "FvM",      ont = "CC", type = "dot", save = TRUE)

# --- Plot GSEA ---
plot_gsea("Training_CRE",     "CREvPLA",  type = "hallmark")
plot_gsea("Training_CRE",     "CREvPLA",  type = "BP")
plot_gsea("Training_CRE",     "CREvPLA",  type = "kegg")
plot_gsea("Training_CRE",     "CREvPLA",  type = "kegg_ridge", save = TRUE)
plot_gsea("Baseline_SURVvCTL","SURVvCTL", type = "MF",         save = TRUE)

p <- plot_ora("Baseline_SURVvCTL", "SURVvCTL", ont = "MF", type = "bar", show_n = 10)
p

# =============================================================================
# BATCH PLOT — All datasets × all types, saved to file
# =============================================================================

# --- Define what to plot per function ---
ora_types  <- c("BP", "CC", "MF")
gsea_types <- c("hallmark", "BP", "CC", "MF", "kegg", "kegg_ridge")

# --- Loop over all datasets ---
lapply(names(diff_only_list), function(nm) {
  folder <- folder_map[[nm]]
  cat(sprintf("\n>> Plotting: %s\n", nm))
  
  # ORA plots (dot + bar for each ontology)
  for (ont in ora_types) {
    for (type in c("dot", "bar")) {
      tryCatch(
        plot_ora(nm, folder, ont = ont, type = type, save = TRUE),
        error = function(e) cat(sprintf("   ORA %s %s failed: %s\n", ont, type, e$message))
      )
    }
  }
  
  # GSEA plots
  for (type in gsea_types) {
    tryCatch(
      plot_gsea(nm, folder, type = type, save = TRUE),
      error = function(e) cat(sprintf("   GSEA %s failed: %s\n", type, e$message))
    )
  }
})

cat("\nAll plots saved.\n")