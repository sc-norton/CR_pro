# === PROJECT: USC - CR Plasma Proteomics =============================
# === Script: Normalization 

#### PURPOSE: To load the raw data and metadata, assemble a DAList for use in
# proteoDA, filter the data for contaminants and missingness, and normalize the
# for downstream analyses.

#### OUTPUTS: (a) pre-normalized datasets (for normalization and imputation checks)
#             (b) normalized datasets
#             (c) normalized DALists
#             (d) reports

#### WORKFLOW:
#     0. Setup - load required packages, establish directories
#     1. Load & Restructure Data - load raw data, metadata; restructure and
#        align
####  2. Contaminant Identification -
####  3. Assemble DAList -
####  4. Filter Data -
####  6. Outlier Removal - 
#### 10. Normalize Data -
#### 11. Normalize Imputed Data -
#### 12. Assess the Normalized Data (QC) -
#### 13. Export Normalized Data - 

# --- 0: SETUP --------------------------------------------
cat("\n>> 0 - Setup\n")
cat("\n Loading required packages...")

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork)

cat("\n Packages successfully loaded.\n")

# Identify directories
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
input_dir  <- file.path(base_dir, "01_input")
report_dir <- file.path(base_dir, "02_normalization", "b_reports")
data_dir   <- file.path(base_dir, "02_normalization", "c_data")

cat("\n Base directory: ", base_dir, "\n")

# --- 1: LOAD & RESTRUCTURE DATA ---------------------------------------
cat("\n>> 1 - Load & Restructure Data\n")

# Load Data
raw <- read.csv(file.path(input_dir, "CRp_raw.csv"))
meta <- read.csv(file.path(input_dir, "CRp_meta.csv"))

# Rename annotation column names
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq", "sum_spectra")
colnames(raw)[1:6] <- annot_cols

# Make Uniprot_ID the row name
rownames(raw) <- raw$uniprot_id

# Split input data into protein intensity data and annotation data
intensity <- raw[ , setdiff(names(raw), annot_cols)] 
annotation <- raw[ ,annot_cols] 

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))
cat(sprintf("\n Meta data: %d samples\n", nrow(meta)))

# Align row names of metadata with column names of data
cat(" Aligning data and metadata...\n")
rownames(meta) <- meta[ ,1]

meta <- meta[order(meta$sample_id), ]

cat(" Data aligned to metadata.\n")

# Validate alignment
cat(" Validating alignment...\n")

data_samples <- colnames(intensity)
meta_samples <- meta$sample_id

stopifnot(
  "Sample mismatch between data and metadata" =
    setequal(data_samples, meta_samples)
)

# Reorder intensity columns to match metadata row order
intensity <- intensity[, meta_samples]

cat("   Validation passed: all", length(meta_samples), "samples aligned.\n")

# Sample summary
cat("\n   Sample distribution:\n")
print(meta %>% count(cancer, supp, timepoint, cancer_time, supp_time))

sample_summary <- as.data.frame(meta %>% 
                                  count(cancer, 
                                        supp, 
                                        timepoint, 
                                        cancer_time, 
                                        supp_time))

write.csv(sample_summary, file.path(report_dir, "00_raw_sample_summary.csv"))

# --- 2: CONTAMINANT IDENTIFICATION with HPA----------------------------------
cat("\n>> 2 - Contaminant Identification\n")

# Identify keratins
keratins <- annotation %>%
  filter(grepl("keratin", description, ignore.case = T)) %>%
  pull(uniprot_id)

# Join to annotation data
annotation <- annotation %>%
  mutate(is_keratin = uniprot_id %in% keratins)

n_flagged <- sum(annotation$is_keratin)

cat(sprintf("\n Proteins flagged as contaminants: %d\n", n_flagged))

# Print contaminants
cat("\n Flagged contaminants: \n")
print(annotation %>%
        filter(is_keratin) %>%
        dplyr::select(gene, uniprot_id, description)) 
  
# Export contaminants
flagged_contaminants <- annotation[annotation$is_keratin == T, ]
write.csv(flagged_contaminants, 
          file.path(report_dir, "01_flagged_contaminants.csv"))
cat("\n Saved flagged_contaminants.csv : \n", report_dir, "\n")


# --- 3: ASSEMBLE DALISTS -----------------------------------------
cat("\n>> 3 - Assemble DAList \n")

# Check for duplicate uniprot_ids
cat(" Checking for duplicate uniprot ids...\n")
dup_ids <- annotation$uniprot_id[duplicated(annotation$uniprot_id)]

if (length(dup_ids) > 0) {
  cat(sprintf("   Found %d duplicate uniprot_ids — deduplicating by highest mean intensity.\n",
              length(dup_ids)))
  
  intensity_num <- as.data.frame(lapply(intensity, as.numeric))
  annotation$row_mean <- rowMeans(intensity_num, na.rm = TRUE)
  
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  
  cat(sprintf("   After deduplication: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids found.\n")
}

# Assemble DAList (all samples)
dal <- DAList(
  data       = intensity,
  metadata   = meta,
  annotation = annotation
)

# DAList (w/out CTL)
dal2 <- filter_samples(dal, dal$metadata$cancer != "CTL")

# Combine
dals <- list(SURV = dal,
             CRE = dal2)

cat(sprintf(" DAList assembled: %d proteins x %d samples \n", 
            nrow(dal$data),
            ncol(dal$data)))

cat(sprintf(" DAList2 assembled: %d proteins x %d samples \n", 
            nrow(dal2$data),
            ncol(dal2$data)))

# --- 4: FILTER DATA -------------------------------------
cat("\n>> 4 - Filter Data \n")

grouping_cols <- c(SURV = "cancer_time", CRE = "supp_time")

# Filtering function
log_step <- function(log, step, n_before, n_after) {
  bind_rows(log, tibble(
    step = step, 
    n_before = n_before, 
    n_after = n_after, 
    n_removed = n_before - n_after))
}

# Build filter plot function
make_filter_plot <- function(log, label) {
  ggplot(mutate(log, step = factor(step, levels = step)), aes(x = step, y = n_after)) +
    geom_col(fill = "steelblue", width = 0.5) +
    geom_text(aes(label = n_after), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "Number of Proteins",
         title = paste("Proteins Retained Through Filtering -", label)) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

f.dals <- lapply(names(dals), function(nm) {
  d <- dals[[nm]]
  grp <- grouping_cols[[nm]]
  log <- tibble(step = character(),
                n_before = integer(),
                n_after = integer(),
                n_removed = integer())
  n_start <- nrow(d$data)
  
  # Step 1: Convert 0s to NA
  cat(sprintf(" [%s] Filter Step 1: Convert 0s to NA\n", nm))
  d <- zero_to_missing(d)
  cat(sprintf("   Converted zeros to NA. Proteins: %d\n", nrow(d$data)))
  
  # Step 2: Remove contaminants
  cat(sprintf(" [%s] Filter Step 2: Remove Contaminants\n", nm))
  n_before <- nrow(d$data)
  d <- filter_proteins_by_annotation(d, !is_keratin)
  n_after <- nrow(d$data)
  log <- log_step(log, "Contaminant removal", n_before, n_after)
  cat(sprintf("   Contaminant removal: %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  # Step 3: Missingness filter
  cat(sprintf(" [%s] Filter Step 3: Filter by Missingness\n", nm))
  n_before <- nrow(d$data)
  d <- filter_proteins_by_proportion(d, min_prop = 0.5, grouping_column = grp)
  n_after <- nrow(d$data)
  log <- log_step(log, "Missingness filter", n_before, n_after)
  cat(sprintf("   Missingness filter: %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  # Finalize filter log
  log <- bind_rows(
    tibble(step = "Raw input", n_before = NA_integer_,
           n_after = n_start, n_removed = NA_integer_),
    log
  ) %>%
    mutate(pct_retained = round(n_after / n_start * 100, 1))
  
  list(dal = d, log = log)
})

names(f.dals)<- names(dals)

# Unpack results
dals        <- lapply(f.dals, `[[`, "dal")
filter_logs <- lapply(f.dals, `[[`, "log")

# Print filter summaries
cat("\n   Filtering summary:\n")
lapply(names(filter_logs), function(nm) {
  cat(sprintf("\n  -- %s --\n", nm))
  print(filter_logs[[nm]])
})

# Save filter logs and filtered protein lists
walk(names(filter_logs), function(nm) {
  log <- filter_logs[[nm]]
  d   <- dals[[nm]]
  
  write_csv(log, file.path(report_dir, sprintf("02_filtering_effects_%s.csv", nm)))
  
  kept_uniprots    <- rownames(d$data)
  removed_uniprots <- setdiff(annotation$uniprot_id, kept_uniprots)
  
  filtered_proteins <- annotation %>%
    filter(uniprot_id %in% removed_uniprots) %>%
    dplyr::select(uniprot_id, gene, description, is_keratin) %>%
    mutate(reason = if_else(is_keratin, "Contaminant", "Missingness"))
  
  write_csv(filtered_proteins, file.path(report_dir, sprintf("03_filtered_proteins_%s.csv", nm)))
  cat(sprintf(" Saved: 03_filtering_effects_%s.csv\n", nm))
  cat(sprintf(" Saved: 03_filtered_proteins_%s.csv\n", nm))
})

# Save filter plots
walk(names(filter_logs), function(nm) {
  p <- make_filter_plot(filter_logs[[nm]], nm)
  ggsave(file.path(report_dir, sprintf("04_filtering_plot_%s.pdf", nm)),
         p, width = 8, height = 5)
  cat(sprintf(" Saved: 04_filtering_plot_%s.pdf\n", nm))
})


# --- 5: OUTLIER REMOVAL -----------------------------------------
cat("\n>> 6 - Outlier Removal\n")

run_outlier_removal <- function(d, nm) {
  cat(sprintf("\n [%s] Running outlier detection...\n", nm))
  
  # Per-sample percent missing
  pct_missing <- colMeans(is.na(d$data)) * 100
  
  # Paired missingness (delta)
  paired_subjects <- d$metadata %>%
    filter(cancer == "SURV") %>%
    count(pid) %>%
    filter(n == 2) %>%
    pull(pid)
  
  delta_missing <- tibble(sample_id = character(), pid = character(),
                          pct_missing = numeric(), delta_missing = numeric())
  
  for (subj in paired_subjects) {
    rows  <- d$metadata %>% filter(pid == subj)
    t1_id <- rows$sample_id[rows$timepoint == "T1"]
    t2_id <- rows$sample_id[rows$timepoint == "T2"]
    if (length(t1_id) == 1 && length(t2_id) == 1) {
      delta_missing <- bind_rows(delta_missing, tibble(
        sample_id     = c(t1_id, t2_id),
        pid           = subj,
        pct_missing   = c(pct_missing[t1_id], pct_missing[t2_id]),
        delta_missing = pct_missing[t2_id] - pct_missing[t1_id]
      ))
    }
  }
  
  # Add unpaired samples
  unpaired_ids <- setdiff(d$metadata$sample_id, delta_missing$sample_id)
  if (length(unpaired_ids) > 0) {
    delta_missing <- bind_rows(delta_missing, tibble(
      sample_id     = unpaired_ids,
      pid           = d$metadata$pid[match(unpaired_ids, d$metadata$sample_id)],
      pct_missing   = pct_missing[unpaired_ids],
      delta_missing = NA_real_
    ))
  }
  
  # IQR-based thresholds
  miss_threshold  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
  delta_vals      <- delta_missing$delta_missing[!is.na(delta_missing$delta_missing)]
  
  if (length(delta_vals) > 2) {
    delta_threshold <- quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals)
    delta_lower     <- quantile(delta_vals, 0.25) - 1.5 * IQR(delta_vals)
  } else {
    delta_threshold <- Inf
    delta_lower     <- -Inf
  }
  
  delta_missing <- delta_missing %>%
    mutate(miss_flag = pct_missing > miss_threshold |
             (!is.na(delta_missing) & (delta_missing > delta_threshold | delta_missing < delta_lower)))
  
  # PCA outlier detection
  data_for_pca <- d$data
  for (j in seq_len(ncol(data_for_pca))) {
    nas <- is.na(data_for_pca[, j])
    if (any(nas)) data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
  
  data_log2   <- log2(data_for_pca + 1)
  pca_res     <- prcomp(t(data_log2), center = TRUE, scale. = TRUE)
  pc_scores   <- pca_res$x[, 1:3]
  mahal_dist  <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
  
  pca_flags <- tibble(
    sample_id  = colnames(d$data),
    mahal_dist = mahal_dist,
    pca_flag   = mahal_dist > qchisq(0.99, df = 3)
  )
  
  # MAD-based intensity outlier detection
  sample_medians <- apply(log2(d$data + 1), 2, median, na.rm = TRUE)
  global_median  <- median(sample_medians)
  mad_val        <- mad(sample_medians)
  
  mad_flags <- tibble(
    sample_id     = names(sample_medians),
    sample_median = sample_medians,
    mad_deviation = abs(sample_medians - global_median),
    mad_flag      = abs(sample_medians - global_median) > 3 * mad_val
  )
  
  # Consensus
  outlier_diag <- delta_missing %>%
    dplyr::select(sample_id, pct_missing, delta_missing, miss_flag) %>%
    left_join(pca_flags, by = "sample_id") %>%
    left_join(mad_flags,  by = "sample_id") %>%
    mutate(n_flags = miss_flag + pca_flag + mad_flag,
           consensus_outlier = n_flags >= 2)
  
  write_csv(outlier_diag, file.path(report_dir, sprintf("05_outlier_diagnostics_%s.csv", nm)))
  
  n_outliers <- sum(outlier_diag$consensus_outlier)
  cat(sprintf("   Outlier consensus: %d sample(s) flagged by >= 2 methods\n", n_outliers))
  
  if (n_outliers > 0) {
    cat("\n   Consensus outliers:\n")
    print(outlier_diag %>%
            filter(consensus_outlier) %>%
            dplyr::select(sample_id, pct_missing, delta_missing, mahal_dist, miss_flag, pca_flag, mad_flag))
  }
  
  # Diagnostic plots
  p1 <- ggplot(outlier_diag %>% filter(!is.na(delta_missing)),
               aes(x = pct_missing, y = delta_missing)) +
    geom_point(aes(color = consensus_outlier), size = 3) +
    geom_hline(yintercept = c(delta_lower, delta_threshold),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = miss_threshold,
               linetype = "dashed", color = "red", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = "% Missing", y = "Delta Missing (T2 - T1)", title = "A: Paired Missingness") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  
  pc_df <- as.data.frame(pca_res$x[, 1:2])
  pc_df$sample_id <- rownames(pc_df)
  pc_df <- left_join(pc_df, outlier_diag %>% dplyr::select(
    sample_id, consensus_outlier), by = "sample_id")
  
  p2 <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = consensus_outlier), size = 3) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = "B: PCA Outliers (Mahalanobis)") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  
  p3 <- ggplot(outlier_diag, aes(x = reorder(sample_id, sample_median), y = sample_median)) +
    geom_point(aes(color = consensus_outlier), size = 2.5) +
    geom_hline(yintercept = global_median, color = "black") +
    geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
               linetype = "dashed", color = "red", alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
    labs(x = "Sample", y = "Median log2 intensity", title = "C: MAD Intensity Outliers") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          legend.position = "none")
  
  ggsave(file.path(report_dir, sprintf("06_outlier_diagnostic_plots_%s.pdf", nm)),
         p1 + p2 + p3 + plot_layout(ncol = 1), width = 10, height = 12)
  cat(sprintf("   Saved: 06_outlier_diagnostic_plots_%s.pdf\n", nm))
  
  # Remove outliers and return cleaned DAList
  if (n_outliers > 0) {
    outlier_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(sample_id)
    cat(sprintf("   Removing %d outlier sample(s): %s\n",
                length(outlier_ids), paste(outlier_ids, collapse = ", ")))
    keep_idx <- !colnames(d$data) %in% outlier_ids
    d$data     <- d$data[, keep_idx, drop = FALSE]
    d$metadata <- d$metadata[d$metadata$sample_id %in% colnames(d$data), ]
    cat(sprintf("   After outlier removal: %d samples remain\n", ncol(d$data)))
  } else {
    cat("   No outlier samples removed.\n")
  }
  
  d
}

# Run on both DALists
dals <- lapply(names(dals), 
               function(nm) run_outlier_removal(dals[[nm]], nm)) |> 
  setNames(names(dals))

# --- 6: NORMALIZE DATA ----------------------------
cat("\n>> 6 - Normalize Data\n")

# Config
grouping_cols <- c(SURV = "cancer_time", CRE = "supp_time")

# Save pre-normalized DALists
lapply(names(dals), function(nm) {
  saveRDS(dals[[nm]], file.path(data_dir, sprintf("00_DAList_pre_norm_%s.rds", nm)))
})

# Pre-normalized reports
lapply(names(dals), function(nm) {
  write_norm_report(dals[[nm]], grouping_column = grouping_cols[nm],
                    output_dir = report_dir,
                    filename   = sprintf("07_normalization_report_%s.pdf", nm),
                    overwrite  = TRUE)
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample_id",
                  output_dir = report_dir,
                  filename   = sprintf("08_qc_report_pre_norm_%s.pdf", nm),
                  overwrite  = TRUE)
})

# Normalize
dals <- lapply(dals, normalize_data, norm_method = "cycloess")

# Post-norm reports
lapply(names(dals), function(nm) {
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample_id",
                  output_dir = report_dir,
                  filename   = sprintf("09_qc_report_post_norm_%s.pdf", nm),
                  overwrite  = TRUE)
})

lapply(names(dals), function(nm) {
  cat(sprintf("  [%s] Final: %d proteins x %d samples\n", nm,
              nrow(dals[[nm]]$data), ncol(dals[[nm]]$data)))
})

# --- 6b: Normalization Quality Scoring ---------------
cat("Step 6b: Norm quality scoring\n")

compute_pcv <- function(mat, groups) {
  cvs <- unlist(lapply(unique(groups), function(grp) {
    sub <- mat[, groups == grp, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2 || mean(x) == 0) return(NA_real_)
      sd(x) / abs(mean(x))
    })
  }))
  median(cvs, na.rm = TRUE)
}

compute_pmad <- function(mat, groups) {
  mads <- unlist(lapply(unique(groups), function(grp) {
    sub <- mat[, groups == grp, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      mad(x, constant = 1)
    })
  }))
  median(mads, na.rm = TRUE)
}

compute_cor <- function(mat, groups) {
  cors <- unlist(lapply(unique(groups), function(grp) {
    idx <- which(groups == grp)
    if (length(idx) < 2) return(NULL)
    cm <- cor(mat[, idx], use = "pairwise.complete.obs")
    cm[lower.tri(cm)]
  }))
  mean(cors, na.rm = TRUE)
}

NORMS         <- c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi")
grouping_cols <- c(SURV = "cancer_time", CRE = "supp_time")   # matches rest of pipeline

all_norm_scores <- list()

for (ds_name in names(dals)) {
  cat(sprintf("\n  [%s] Scoring normalisation methods...\n", ds_name))
  
  pre_dal <- readRDS(file.path(data_dir, sprintf("00_DAList_pre_norm_%s.rds", ds_name)))
  grp_col <- grouping_cols[[ds_name]]
  grps    <- factor(pre_dal$metadata[[grp_col]])
  
  norm_scores <- purrr::map_dfr(NORMS, function(nm) {
    dal_n <- tryCatch(
      normalize_data(pre_dal, norm_method = nm),
      error = function(e) {
        cat(sprintf("    %s: FAILED (%s)\n", nm, e$message))
        NULL
      }
    )
    if (is.null(dal_n)) return(NULL)
    
    mat_n <- as.matrix(dal_n$data)
    
    pcv     <- compute_pcv(mat_n, grps)
    pmad    <- compute_pmad(mat_n, grps)
    cor_val <- compute_cor(mat_n, grps)
    
    cat(sprintf("    %s: PCV=%.4f  PMAD=%.4f  COR=%.4f\n", nm, pcv, pmad, cor_val))
    tibble(dataset = ds_name, norm = nm,
           PCV = round(pcv, 4), PMAD = round(pmad, 4), COR = round(cor_val, 4))
  })
  
  norm_scores <- norm_scores |>
    mutate(
      PCV_rank        = rank(PCV),
      PMAD_rank       = rank(PMAD),
      COR_rank        = rank(-COR),
      norm_composite  = (PCV_rank + PMAD_rank + COR_rank) / 3
    ) |>
    arrange(norm_composite)
  
  cat(sprintf("\n  [%s] Normalisation rankings:\n", ds_name))
  print(norm_scores |> dplyr::select(norm, PCV, PMAD, COR, norm_composite))
  
  write_csv(norm_scores,
            file.path(report_dir, sprintf("10_norm_quality_scores_%s.csv", ds_name)))
  
  all_norm_scores[[ds_name]] <- norm_scores
}
# --- Export Normalized Data ----------------------
cat("\n>> 13 - Export Normalized Data\n")

invisible(lapply(names(dals), function(nm) {
  d <- dals[[nm]]
  
  export <- bind_cols(
    d$annotation %>% dplyr::select(uniprot_id, protein, gene, description),
    d$data
  )
  
  write.csv(export, file.path(data_dir, sprintf("01_normalized_data_%s.csv", nm)))
  saveRDS(d,        file.path(data_dir, sprintf("01_normalized_DAList_%s.RDS", nm)))
  
  cat(sprintf("  [%s] Exported: %d proteins x %d samples\n", nm, nrow(d$data), ncol(d$data)))
}))

# --- Protein Assessments ------------

invisible(lapply(names(dals), function(nm) {
  d         <- dals[[nm]]
  miss_pct  <- rowMeans(is.na(d$data))
  
  report <- tibble(
    Total      = nrow(d$data),
    Complete   = sum(miss_pct == 0),
    Incomplete = sum(miss_pct > 0)
  )
  
  write.csv(report, file.path(report_dir, sprintf("10_protein_report_%s.csv", nm)))
  cat(sprintf("  [%s] Total: %d | Complete: %d | Incomplete: %d\n", 
              nm, report$Total, report$Complete, report$Incomplete))
}))

# --- Sample Assessments -----------------

invisible(lapply(names(dals), function(nm) {
  meta <- dals[[nm]]$metadata
  
  # Timepoint counts
  t1_count <- sum(meta$timepoint == "T1", na.rm = T)
  t2_count <- sum(meta$timepoint == "T2", na.rm = T)
  
  # Supplement counts
  CRE_count <- sum(meta$supp == "CRE", na.rm = T)
  PLA_count <- sum(meta$supp == "PLA", na.rm = T)
  
  # Supp_timee counts
  CRE_T1_count <- sum(meta$supp_time == "CRE_T1", na.rm = T)
  CRE_T2_count <- sum(meta$supp_time == "CRE_T2", na.rm = T)
  PLA_T1_count <- sum(meta$supp_time == "PLA_T1", na.rm = T)
  PLA_T2_count <- sum(meta$supp_time == "PLA_T2", na.rm = T)
  
  # Paired vs unpaired subjects
  paired_subjects <- meta %>% 
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_cre <- meta %>%
    filter(supp == "CRE") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_pla <- meta %>%
    filter(supp == "PLA") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  unpaired_subjects <- meta %>%
    count(pid) %>%
    filter(n == 1) %>%
    nrow()
  
  # Controls (only relevant for SURV v CTL)
  n_ctl <- sum(meta$cancer == "CTL", na.rm = T)
  n_surv <- sum(meta$cancer == "SURV", na.rm = T)
  
  # T1 and T2 only
  t1_only <- t1_count - paired_subjects
  t2_only <- t2_count - paired_subjects
  
  report <- tibble(
    total_samples       = nrow(meta),
    n_CRE               = CRE_count,
    n_PLA               = PLA_count,
    n_CRE_T1            = CRE_T1_count,
    n_CRE_T2            = CRE_T2_count,
    n_PLA_T1            = PLA_T1_count,
    n_PLA_T2            = PLA_T2_count,
    n_T1                = t1_count,
    n_T2                = t2_count,
    n_paired            = paired_subjects,
    n_paired_CRE        = paired_cre,
    n_paired_PLA        = paired_pla,
    n_unpaired          = unpaired_subjects,
    n_T1_only           = t1_only,
    n_T2_only           = t2_only,
    n_SURV              = n_surv,
    n_CTL               = n_ctl
  )
  
  write.csv(report, file.path(report_dir, sprintf("11_sample_report_%s.csv", nm)))
}))

