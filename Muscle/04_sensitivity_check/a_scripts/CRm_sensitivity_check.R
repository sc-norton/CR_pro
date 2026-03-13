# === PROJECT: USC - CR Muscle Proteomics =============================
# === Script: Sensitivity Check

# === PURPOSE: 

# === INPUTS: 

# === OUTPUTS: 
# === WORKFLOW:
#     0. Setup 
#     1. Load & Restructure Data 
#     2. Contaminant Identification (with HPA)
#     3. Assemble DAList 
#     4. Filter Data 
#     5. Outlier Removal 
#     6. Normalize Data 
#     7. Export Normalized Data 


#!/usr/bin/env Rscript
# sensitivity_analysis.R — SURV + CRE datasets
# Top-3 norms x (top-3 imps + none) = 12 combos per dataset

# --- 0: Setup --------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, MsCoreUtils, msImpute, pcaMethods, imputeLCMD, missForest,
               limma, readr, dplyr, tidyr, tibble, stringr,
               ggplot2, patchwork, scales)

# --- Paths ---
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
report_dir <- file.path(base_dir, "04_sensitivity_check", "b_reports")
data_dir   <- file.path(base_dir, "04_sensitivity_check", "c_data")
norm_dir   <- file.path(base_dir, "02_normalization")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

# --- Dataset-specific config ---
# Contrasts differ: SURV has CTL samples, CRE does not
DATASET_CONFIG <- list(
  SURV = list(
    pre_norm_file    = file.path(norm_dir, "c_data", "00_DAList_pre_norm_SURV.rds"),
    norm_scores_file = file.path(norm_dir, "b_reports", "10_norm_quality_scores_SURV.csv"),
    group_col        = "cancer_time",
    block_col        = "pid",
    contrasts = c(
      Baseline_SURVvCTL      = "SURV_T1 - CTL_T1",
      Training_SURVvCTL      = "SURV_T2 - CTL_T1",
      Training_SURV          = "SURV_T2 - SURV_T1"
    )
  ),
  
  CRE = list(
    pre_norm_file    = file.path(norm_dir, "c_data", "00_DAList_pre_norm_CRE.rds"),
    norm_scores_file = file.path(norm_dir, "b_reports", "10_norm_quality_scores_CRE.csv"),
    group_col        = "supp_time",
    block_col        = "pid",
    contrasts = c(
      Baseline_CREvPLA = "CRE_T1 - PLA_T1",
      Training_CREvPLA = "CRE_T2 - PLA_T2",
      Training_CRE     = "CRE_T2 - CRE_T1",
      Training_PLA     = "PLA_T2 - PLA_T1",
      Interaction      = "(CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
    )
  )
)


# --- Helper: build impute specs from benchmark file ---
build_imp_specs <- function(bench_file, top_n = 3) {
  if (!file.exists(bench_file)) {
    warning(sprintf("Benchmark file not found: %s — using 'none' only", bench_file))
    return(list(none = NULL))
  }
  bench    <- read_csv(bench_file, show_col_types = FALSE) |> arrange(mean)
  top_imps <- bench$method[seq_len(min(top_n, nrow(bench)))]
  cat(sprintf("  Top %d imputations: %s\n", top_n, paste(top_imps, collapse = ", ")))
  
  specs <- list(none = NULL)
  for (imp_name in top_imps) {
    if (grepl("^mix_", imp_name)) {
      parts <- strsplit(sub("^mix_", "", imp_name), "_")[[1]]
      specs[[imp_name]] <- list(method = "mixed", mar = parts[1], mnar = parts[2])
    } else {
      specs[[imp_name]] <- list(method = imp_name)
    }
  }
  specs
}

# --- Helper: run one norm+impute+limma pipeline ---
run_pipeline <- function(dal, norm_method, imp_spec,
                         group_col, block_col, contrast_vec) {
  
  # Normalise
  dal_n <- tryCatch(
    normalize_data(dal, norm_method = norm_method),
    error = function(e) { warning(sprintf("Norm failed [%s]: %s", norm_method, e$message)); NULL }
  )
  if (is.null(dal_n)) return(NULL)
  
  mat  <- as.matrix(dal_n$data)
  meta <- dal_n$metadata
  
  # Impute
  if (!is.null(imp_spec)) {
    mat <- tryCatch({
      if (imp_spec$method == "mixed") {
        rn <- setNames(rep(FALSE, nrow(mat)), rownames(mat))  # fallback if no MAR/MNAR file
        impute_matrix(mat, method = "mixed", randna = rn,
                      mar = imp_spec$mar, mnar = imp_spec$mnar)
      } else {
        impute_matrix(mat, method = imp_spec$method)
      }
    }, error = function(e) {
      warning(sprintf("Impute failed [%s]: %s", imp_spec$method, e$message)); NULL
    })
    if (is.null(mat)) return(NULL)
  }
  
  # Design — use the dataset-specific grouping column
  groups <- factor(meta[[group_col]])
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  # Build contrast matrix from named character vector
  contrast_mat <- tryCatch(
    makeContrasts(contrasts = contrast_vec, levels = design),
    error = function(e) {
      warning(sprintf("makeContrasts failed: %s", e$message)); NULL
    }
  )
  if (is.null(contrast_mat)) return(NULL)
  
  # Weights + duplicate correlation
  aw <- tryCatch(arrayWeights(mat, design),
                 error = function(e) { rep(1, ncol(mat)) })
  
  dupcor <- tryCatch(
    duplicateCorrelation(mat, design, block = meta[[block_col]], weights = aw),
    error = function(e) list(consensus.correlation = 0)
  )
  
  fit <- lmFit(mat, design,
               block       = meta[[block_col]],
               correlation = dupcor$consensus.correlation,
               weights     = aw) |>
    contrasts.fit(contrast_mat) |>
    eBayes(robust = TRUE, trend = TRUE)
  
  # Extract per-contrast topTable
  setNames(lapply(colnames(contrast_mat), function(cn) {
    tt       <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
    tt$gene  <- rownames(mat)
    tt$Pi.Val     <- abs(tt$logFC) * (-log10(pmax(tt$adj.P.Val, .Machine$double.eps)))
    tt$adj.Pi.Val <- 10^-tt$Pi.Val
    tt
  }), colnames(contrast_mat))
}

# --- Helper: pairwise Spearman rho across methods ---
compute_spearman <- function(all_results, contrast_names) {
  all_genes    <- unique(unlist(lapply(all_results, function(r) r[[1]]$gene)))
  method_names <- names(all_results)
  
  spearman_long <- list()
  for (cname in contrast_names) {
    t_mat <- sapply(method_names, function(lab) {
      tt     <- all_results[[lab]][[cname]]
      t_vals <- setNames(rep(NA_real_, length(all_genes)), all_genes)
      t_vals[tt$gene] <- tt$t
      t_vals
    })
    sm    <- cor(t_mat, method = "spearman", use = "pairwise.complete.obs")
    pairs <- which(lower.tri(sm), arr.ind = TRUE)
    spearman_long[[cname]] <- tibble(
      contrast     = cname,
      method_a     = rownames(sm)[pairs[, 1]],
      method_b     = colnames(sm)[pairs[, 2]],
      spearman_rho = sm[pairs]
    )
    cat(sprintf("  %s: mean rho = %.4f, min rho = %.4f\n",
                cname,
                mean(sm[lower.tri(sm)], na.rm = TRUE),
                min(sm[lower.tri(sm)],  na.rm = TRUE)))
  }
  list(long = bind_rows(spearman_long), matrices = setNames(
    lapply(contrast_names, function(cn) {
      t_mat <- sapply(method_names, function(lab) {
        tt <- all_results[[lab]][[cn]]
        tv <- setNames(rep(NA_real_, length(all_genes)), all_genes)
        tv[tt$gene] <- tt$t; tv
      })
      cor(t_mat, method = "spearman", use = "pairwise.complete.obs")
    }), contrast_names)
  )
}

# --- Helper: DEP counts across criteria ---
compute_dep_counts <- function(all_results, contrast_names) {
  bind_rows(lapply(names(all_results), function(lab) {
    bind_rows(lapply(contrast_names, function(cname) {
      tt <- all_results[[lab]][[cname]]
      bind_rows(lapply(
        list(
          nominal_p05 = tt$P.Value   < 0.05,
          FDR_05      = tt$adj.P.Val < 0.05,
          FDR_10      = tt$adj.P.Val < 0.10,
          pi_05       = tt$adj.Pi.Val < 0.05
        ),
        function(mask) tibble(n_up = sum(tt$logFC[mask] > 0, na.rm = TRUE),
                              n_down = sum(tt$logFC[mask] < 0, na.rm = TRUE),
                              n_total = sum(mask, na.rm = TRUE))
      ), .id = "criterion") |>
        mutate(method = lab, contrast = cname)
    }))
  }))
}

# --- Helper: CAT curves ---
compute_cat_curves <- function(all_results, contrast_names, max_k = 300) {
  all_genes    <- unique(unlist(lapply(all_results, function(r) r[[1]]$gene)))
  method_names <- names(all_results)
  
  bind_rows(lapply(contrast_names, function(cname) {
    # Consensus ranking = mean rank across all methods
    consensus_ranks <- rowMeans(sapply(all_results, function(res) {
      tt <- res[[cname]]
      r  <- setNames(rep(NA_real_, length(all_genes)), all_genes)
      r[tt$gene] <- rank(-abs(tt$t))
      r
    }), na.rm = TRUE)
    consensus_ranks <- consensus_ranks[!is.na(consensus_ranks)]
    
    bind_rows(lapply(method_names, function(lab) {
      tt           <- all_results[[lab]][[cname]]
      method_ranks <- setNames(rank(-abs(tt$t)), tt$gene)
      shared       <- intersect(names(method_ranks), names(consensus_ranks))
      ord_m <- names(sort(method_ranks[shared]))
      ord_c <- names(sort(consensus_ranks[shared]))
      ks    <- seq_len(min(max_k, length(shared)))
      tibble(
        k           = ks,
        concordance = sapply(ks, function(k) length(intersect(ord_m[1:k], ord_c[1:k])) / k),
        method      = lab,
        contrast    = cname
      )
    }))
  }))
}

# =============================================================================
# MAIN: Loop over datasets
# =============================================================================
all_dataset_results <- list()   # stores full results per dataset for cross-dataset use

for (ds_name in names(DATASET_CONFIG)) {
  cfg <- DATASET_CONFIG[[ds_name]]
  cat(sprintf("\n\n==============================================\n"))
  cat(sprintf(" DATASET: %s\n", ds_name))
  cat(sprintf("==============================================\n"))
  
  # Output dirs per dataset
  ds_report_dir <- file.path(report_dir, ds_name)
  ds_data_dir   <- file.path(data_dir,   ds_name)
  dir.create(ds_report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_data_dir,   recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. Load data & rankings ---
  cat("Step 1: Load data & rankings\n")
  stopifnot(file.exists(cfg$pre_norm_file))
  
  dal_raw <- readRDS(cfg$pre_norm_file)
  cat(sprintf("  Loaded: %d proteins x %d samples\n",
              nrow(dal_raw$data), ncol(dal_raw$data)))
  
  norm_scores <- read_csv(cfg$norm_scores_file, show_col_types = FALSE)
  top_norms   <- norm_scores$norm[1:3]
  cat(sprintf("  Top 3 norms: %s\n", paste(top_norms, collapse = ", ")))
  
  # Imputation benchmark — optional per dataset (may share one file)
  bench_file <- file.path(base_dir, "03_imputation", "c_data",
                          sprintf("benchmark_summary_%s.csv", ds_name))
  if (!file.exists(bench_file)) {
    # Fall back to shared benchmark if dataset-specific one doesn't exist
    bench_file <- file.path(base_dir, "03_imputation", "c_data", "benchmark_summary.csv")
  }
  IMP_SPECS <- build_imp_specs(bench_file, top_n = 3)
  
  # --- 2. Build combo grid ---
  combos <- expand.grid(norm = top_norms, imp = names(IMP_SPECS),
                        stringsAsFactors = FALSE)
  combos$label <- paste0(combos$norm, "_", combos$imp)
  cat(sprintf("  Grid: %d combos\n", nrow(combos)))
  
  contrast_names <- names(cfg$contrasts)
  
  # --- 3. Run pipeline ---
  cat("Step 3: Run combo grid\n")
  set.seed(42)
  all_results <- list()
  
  for (i in seq_len(nrow(combos))) {
    lab <- combos$label[i]
    cat(sprintf("  [%d/%d] %s\n", i, nrow(combos), lab))
    
    res <- tryCatch(
      run_pipeline(
        dal           = dal_raw,
        norm_method   = combos$norm[i],
        imp_spec      = IMP_SPECS[[combos$imp[i]]],
        group_col     = cfg$group_col,
        block_col     = cfg$block_col,
        contrast_vec  = cfg$contrasts
      ),
      error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL }
    )
    if (!is.null(res)) all_results[[lab]] <- res
  }
  
  cat(sprintf("  Completed: %d / %d\n", length(all_results), nrow(combos)))
  if (length(all_results) == 0) { cat("  No results — skipping dataset.\n"); next }
  
  method_names <- names(all_results)
  
  # --- 4. Spearman rho ---
  cat("Step 4: Pairwise Spearman rho\n")
  spearman_out <- compute_spearman(all_results, contrast_names)
  write_csv(spearman_out$long, file.path(ds_data_dir, "pairwise_spearman.csv"))
  
  # --- 5. CAT curves ---
  cat("Step 5: CAT curves\n")
  cat_df <- compute_cat_curves(all_results, contrast_names)
  write_csv(cat_df, file.path(ds_data_dir, "cat_curves.csv"))
  
  # --- 6. DEP counts ---
  cat("Step 6: DEP counts\n")
  dep_df <- compute_dep_counts(all_results, contrast_names)
  write_csv(dep_df, file.path(ds_data_dir, "dep_counts.csv"))
  
  # --- 7. Summary table ---
  cat("Step 7: Summary table\n")
  sm_list <- spearman_out$matrices
  
  summary_df <- bind_rows(lapply(method_names, function(lab) {
    mean_rhos <- sapply(contrast_names, function(cn) {
      sm <- sm_list[[cn]]
      mean(sm[lab, setdiff(method_names, lab)], na.rm = TRUE)
    })
    
    dep_fdr05 <- dep_df |>
      filter(method == lab, criterion == "FDR_05") |>
      select(contrast, n_total) |>
      pivot_wider(names_from = contrast, values_from = n_total,
                  names_prefix = "n_DEP_")
    
    parts <- strsplit(lab, "_", fixed = TRUE)[[1]]
    bind_cols(
      tibble(
        method        = lab,
        norm          = parts[1],
        imputation    = paste(parts[-1], collapse = "_"),
        mean_spearman = round(mean(mean_rhos), 4)
      ),
      as_tibble(t(round(mean_rhos, 4))) |> setNames(paste0("rho_", contrast_names)),
      dep_fdr05 |> select(-any_of("method"))
    )
  })) |> arrange(desc(mean_spearman))
  
  print(summary_df)
  write_csv(summary_df, file.path(ds_data_dir, "sensitivity_comparison.csv"))
  
  # --- 8. Plots ---
  cat("Step 8: Plots\n")
  imp_order <- combos |>
    arrange(factor(imp, levels = names(IMP_SPECS)), norm) |>
    pull(label)
  
  pdf(file.path(ds_report_dir, "sensitivity_comparison.pdf"), width = 14, height = 10)
  
  # Spearman heatmap
  heat_df <- bind_rows(lapply(contrast_names, function(cn) {
    sm <- sm_list[[cn]]
    expand.grid(method_a = factor(rownames(sm), levels = imp_order),
                method_b = factor(colnames(sm), levels = imp_order),
                stringsAsFactors = FALSE) |>
      mutate(rho = as.vector(sm), contrast = cn)
  }))
  
  p_heat <- ggplot(heat_df, aes(x = method_a, y = method_b, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", rho)), size = 2) +
    facet_wrap(~ contrast, ncol = 2) +
    scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                         midpoint = 0.95, name = "Spearman\nrho") +
    labs(title = sprintf("[%s] Pairwise Spearman rho of t-statistics", ds_name),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  print(p_heat)
  
  # CAT curves
  cat_plots <- lapply(contrast_names, function(cn) {
    ggplot(cat_df |> filter(contrast == cn),
           aes(x = k, y = concordance,
               color = factor(method, levels = imp_order))) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0.7, linetype = "dashed", alpha = 0.4) +
      labs(x = "Top k proteins", y = "Concordance", title = cn, color = "Pipeline") +
      theme_minimal(base_size = 9) +
      scale_color_brewer(palette = "Paired") +
      ylim(0, 1)
  })
  print(wrap_plots(cat_plots, ncol = 1) +
          plot_annotation(title = sprintf("[%s] CAT curves vs consensus", ds_name)))
  
  # DEP bar chart
  dep_long <- dep_df |>
    filter(criterion == "FDR_05") |>
    mutate(method = factor(method, levels = imp_order)) |>
    pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "count") |>
    mutate(direction = ifelse(direction == "n_up", "Up", "Down"),
           count = ifelse(direction == "Down", -count, count))
  
  print(
    ggplot(dep_long, aes(x = method, y = count, fill = direction)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, color = "grey30") +
      facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
      labs(title = sprintf("[%s] DEP counts per pipeline (FDR < 0.05)", ds_name),
           x = NULL, y = "DEP count") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold"))
  )
  
  dev.off()
  
  # Store for any cross-dataset comparisons downstream
  all_dataset_results[[ds_name]] <- list(
    results    = all_results,
    spearman   = spearman_out,
    dep        = dep_df,
    cat        = cat_df,
    summary    = summary_df
  )
  
  cat(sprintf("\n  [%s] Done. Outputs in: %s\n", ds_name, ds_data_dir))
}

cat("\n\n==============================================\n")
cat(" ALL DATASETS COMPLETE\n")
cat("==============================================\n")
