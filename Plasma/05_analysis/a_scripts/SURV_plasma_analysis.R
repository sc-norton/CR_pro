#### PROTEOMIC ANALYSIS SCRIPT 

#### PURPOSE: 

# Workflow: 
#  0. Setup -
#  1. 

# --- 0: SETUP ---------------------------------------------------------------
cat(">> 0 — Setup\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA,
  readxl, readr, dplyr, tidyr, stringr,
  httr, xml2,
  ggplot2, ggVennDiagram, patchwork, knitr, ggridges,
  limma, enrichplot, clusterProfiler, org.Hs.eg.db,
  writexl, openxlsx
)

# Establish directories
base_dir <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = T)
report_dir <- file.path(base_dir, "03_analysis", "b_reports")

dir.create(report_dir, recursive = T, showWarnings = F)

cat(" Base directory:", base_dir, "\n")

#--- 1: LOAD & SUBSET DATA ----------------------------------------------------
cat("\n>> 1 - Load Normalized Data\n")

dalist <- readRDS(file.path(base_dir, 
                          "02_normalization", 
                          "c_data", 
                          "01_normalized_DAList.RDS"))

data <- as.data.frame(dalist$data)
meta <- as.data.frame(dalist$metadata)
annot <- as.data.frame(dalist$annotation)

# Subset data for CREvPLA analyses
cre_data <- data %>%
  dplyr::select(-contains("PPS"))

cre_meta <- meta %>%
  filter(cancer == "SURV")

# Confirm cre_data matches cre_meta
cre_data_samples <- colnames(cre_data)
cre_meta_samples <- cre_meta$sample_id

stopifnot(
  "Sample mismatch between data and metadata" =
    setequal(cre_data_samples, cre_meta_samples)
)
cat("\n CREvPLA datasets match, filtering successful\n")

# Create separate DAList for CREvPLA
cre_dal <- DAList(
  data       = cre_data,
  metadata   = cre_meta,
  annotation = annot
)

#--- ANALYSIS WITH LIMMA --------------------------------------------------
# (1) DEFINE EXPERIMENTAL FACTORS
# (2) CREATE DESIGN MATRICES
# (3) CREATE CONTRASTS
# (4) 


#--- DEFINE EXPERIMENTAL FACTORS ---------------------------------------------

# for CRE v PLA:
supp <- factor(cre_meta$supp)
time <- factor(cre_meta$timepoint)
pid_cre <- factor(cre_meta$pid)
supp_time <- factor(cre_meta$supp_time)

#for SURV v CTL:
cancer <- factor(meta$cancer)
cancer_time <- factor(meta$cancer_time)
pid <- factor(meta$pid)

# sex <- meta$Sex
# age <- as.numeric(meta$Age)
# BF <- meta$perBFpre
# BMI <- meta$BMIpre
# fCSA <- meta$fCSApre
# nuclei <- meta$nuclei_pre
# SC <- meta$SC_pre
# VLmCSA <- meta$VLmCSA_pre

#--- BUILD DESIGN MATRICES & CONTRASTS ----------------------------------------

# CRE v PLA comparison
design_sup <- model.matrix(~ 0 + supp_time, data = cre_data)
colnames(design_sup) <- gsub("^supp_time", "", colnames(design_sup))

cat(sprintf(" Supplement Design Matrix: %d samples x %d groups\n",
            nrow(design_sup), ncol(design_sup)))
design_sup

contrast_sup <- makeContrasts(
  Baseline_CRE_PLA = CRE_T1 - PLA_T1,
  Training_CRE_PLA = CRE_T2 - PLA_T2,
  Training_CRE = CRE_T2 - CRE_T1,
  Training_PLA = PLA_T2 - PLA_T1,
  Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1),
  levels = design_sup
)

dupcor_sup <- duplicateCorrelation(cre_data, 
                                   design = design_sup, 
                                   block = pid_cre)


fit_sup <- lmFit(cre_data,
                 design = design_sup,
                 block = pid_cre,
                 correlation = dupcor_sup$consensus.correlation)

fit_sup2 <- contrasts.fit(fit_sup, contrast_sup)
fit_sup2 <- eBayes(fit_sup2)

results_list <- list(
  Baseline_CRE_PLA = topTable(fit_sup2, coef = "Baseline_CRE_PLA",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Training_CRE_PLA = topTable(fit_sup2, coef = "Training_CRE_PLA",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Training_CRE     = topTable(fit_sup2, coef = "Training_CRE",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Training_PLA     = topTable(fit_sup2, coef = "Training_PLA",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Interaction      = topTable(fit_sup2, coef = "Interaction",
                              number = Inf, adjust.method = "BH", sort.by = "P")
)


# SURV v CTL comparison (2,4,5):
design_can <- model.matrix(~0 + cancer_time, data = data)
colnames(design_can) <- gsub("^cancer_time", "", colnames(design_can))

cat(sprintf(" Cancer Design Matrix: %d samples x %d groups\n",
            nrow(design_can), ncol(design_can)))
design_can

contrast_can <- makeContrasts(
  Baseline_SURV_CTL = SURV_T1 - CTL_T1,
  Training_SURV_CTL = SURV_T2 - CTL_T1,
  Training_SURV = SURV_T2 - SURV_T1,
  levels = design_can
)

# Duplicate Correlation
dupcor_can <- duplicateCorrelation(data, design = design_can, block = pid)

# Fitting data to linear model
fit_can <- lmFit(data,
                 design = design_can,
                 block = pid,
                 correlation = dupcor_can$consensus.correlation)

fit_can2 <- contrasts.fit(fit_can, contrast_can)
fit_can2 <- eBayes(fit_can2)

# Compile results for SURV v CTL
results_can <- list(
  Baseline_SURV_CTL = topTable(fit_can2, coef = "Baseline_SURV_CTL",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Training_SURV_CTL = topTable(fit_can2, coef = "Training_SURV_CTL",
                              number = Inf, adjust.method = "BH", sort.by = "P"),
  Training_SURV     = topTable(fit_can2, coef = "Training_SURV",
                              number = Inf, adjust.method = "BH", sort.by = "P")
)

# Combine results from CREvPLA and SURVvCTL
results_list <- c(results_list, results_can)

# Add uniprot column to results
results_list <- lapply(results_list, function(df) {
  df$uniprot <- row.names(df)
  # Reorder columns - uniprot first
  df <- df[, c("uniprot", setdiff(colnames(df), "uniprot")), drop = FALSE]
  return(df)
})

# Calculate π-Value
results_list <- lapply(results_list, function(df) {
  df$Pi.Val <- abs(df$logFC) * (-(log10(df$adj.P.Val)))
  return(df)
})

# Calculate ∏-Value
results_list <- lapply(results_list, function(df) {
  df$adj.Pi.Val <- 10^-df$Pi.Val
  return(df)
})

# Export results to excel file
write_xlsx(results_list, file.path(report_dir, "limma_results.xlsx"))
cat("\n Saved: limma_results.xlsx\n", report_dir)

#--- ASSESS DA & ADD TO RESULTS --------------------------------------------

# DA Method Log
DA_log <- tibble(
  dataset = character(),
  criteria = character(),
  n_UP = integer(),
  n_DOWN = integer(),
  n_NDE = integer()
)
 
# Assess DA & Add to Results
results_list <- lapply(results_list, function(df) {
    df$DE.adj <- "NDE"
    df$DE.adj[df$logFC > 0.6 & df$adj.P.Val < 0.05]  <- "UP"
    df$DE.adj[df$logFC < -0.6 & df$adj.P.Val < 0.05] <- "DOWN"
    df
})

# Calculate counts adjusted p-value DA
DA_log1 <- do.call(rbind, lapply(names(results_list), function(name) {
  df <- results_list[[name]]
  data.frame(
    comparison = name,
    method = "adj P-val",  # or extract method from the name if needed
    n_UP = sum(df$DE.adj == "UP"),
    n_DOWN = sum(df$DE.adj == "DOWN"),
    stringsAsFactors = FALSE
  )
}))

results_list <- lapply(results_list, function(df) {
  df$DE.nom <- "NDE"
  df$DE.nom[df$logFC > 0.6 & df$P.Val < 0.05]  <- "UP"
  df$DE.nom[df$logFC < -0.6 & df$P.Val < 0.05] <- "DOWN"
  df
})

# Calculate counts nominal p-value DA
DA_log2 <- do.call(rbind, lapply(names(results_list), function(name) {
  df <- results_list[[name]]
  data.frame(
    comparison = name,
    method = "nom P-val",  # or extract method from the name if needed
    n_UP = sum(df$DE.nom == "UP"),
    n_DOWN = sum(df$DE.nom == "DOWN"),
    stringsAsFactors = FALSE
  )
}))

results_list <- lapply(results_list, function(df) {
  df$DE.pi <- "NDE"
  df$DE.pi[df$logFC > 0.6 & df$adj.Pi.Val < 0.05]  <- "UP"
  df$DE.pi[df$logFC < -0.6 & df$adj.Pi.Val < 0.05] <- "DOWN"
  df
})

# Calculate counts using Pi-value DA
DA_log3 <- do.call(rbind, lapply(names(results_list), function(name) {
  df <- results_list[[name]]
  data.frame(
    comparison = name,
    method = "Pi-val",  # or extract method from the name if needed
    n_UP = sum(df$DE.pi == "UP"),
    n_DOWN = sum(df$DE.pi == "DOWN"),
    stringsAsFactors = FALSE
  )
}))

DA_log_combined <- rbind(DA_log2, DA_log1, DA_log3)
write.csv(DA_log_combined, file.path(report_dir, "DA_methods.csv"))

# Reshape DA_log to long format
plot_DA <- pivot_longer(DA_log_combined, 
                          cols = c(n_UP, n_DOWN),
                          names_to = "direction",
                          values_to = "count")

# Clean up labels
plot_DA$direction <- gsub("n_", "", plot_DA$direction)

# Create grouped bar plot
p_DA <- ggplot(plot_DA, aes(x = comparison, y = count, fill = direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ method, scales = "free_x") +
      scale_fill_manual(values = c("UP" = "#E96A2C", "DOWN" = "#0D234A")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Comparison", y = "Number of Proteins", fill = "Direction",
         title = "Differentially Abundant Proteins by Criteria and Comparison",
        subtitle = "|log2FC| = 0.6, sig level 0.05")

ggsave(file.path(report_dir, "DA_methods.png"), p_DA)


### CREATE VOLCANO PLOTS ### ==================================================

## Create empty list to store results
volcano <- list()

## Run loop to create volcano plot for each dataframe in list
for (name in names(results_list)) {
  df <- results_list[[name]]
  volcano[[name]] <- ggplot(data = df, aes(x=logFC, y=-log10(adj.P.Val), col = DE.adj)) +
    scale_color_manual(values = c("UP" = "turquoise3", "DOWN" = "red1", "NDE" = "grey")) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed") +
    coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) + # set coord plane limits
    labs(color = 'Expression', #legend_title, 
         x = expression("log"[2]*"FC"), 
         y = expression("-log"[10]*"p-value"),
         title = paste(name),
         subtitle = "BH-adj P-val < 0.05") + 
    annotate("text", x = -6, y = 7,
             label = paste0("UP: ", sum(df$DE.adj == "UP"),
                            "\nDOWN: ", sum(df$DE.adj == "DOWN")),
             hjust = 0, vjust = 1, size = 4) +
    theme_classic() +
    theme(
      text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")) # All text bolded
  
}

# Store plots to build outputs
v1 <- volcano$Baseline_CRE_PLA
v2 <- volcano$Training_CRE_PLA
v3 <- volcano$Training_CRE
v4 <- volcano$Training_PLA
v5 <- volcano$Interaction
v6 <- volcano$Baseline_SURV_CTL
v7 <- volcano$Training_SURV_CTL
v8 <- volcano$Training_SURV

# Build outputs
v_CRE <- (v1 | v2 | v3) / (v4 | v5)
v_CRE

v_SURV <- v6 / v7 / v8
v_SURV

# Export plots
pdf(file.path(report_dir, "volcano_plots.pdf"), width = 10, height = 8)
invisible(lapply(volcano, print))
dev.off()

ggsave(file.path(report_dir, "volcano_CRE_PLA.png"), v_CRE,
       width = 12, height = 8, dpi = 300)

ggsave(file.path(report_dir, "volcano_SURV_CTL.png"), v_SURV, 
       width = 8, height = 12, dpi = 300)


## CREATE VENN DIAGRAMS ## ===================================================

# Assuming you have lists of significantly altered proteins
SURV_T1_UP_proteins <- rownames(results_list$Baseline_SURV_CTL[
  results_list$Baseline_SURV_CTL$DE.adj == "UP", ])
SURV_T1_DOWN_proteins <- rownames(results_list$Baseline_SURV_CTL[
  results_list$Baseline_SURV_CTL$DE.adj == "DOWN", ])
SURV_T2_UP_proteins <- rownames(results_list$Training_SURV_CTL[
  results_list$Training_SURV_CTL$DE.adj == "UP", ])
SURV_T2_DOWN_proteins <- rownames(results_list$Training_SURV_CTL[
  results_list$Training_SURV_CTL$DE.adj == "DOWN", ])

CRE_UP_proteins <- rownames(results_list$Training_CRE[
  results_list$Training_CRE$DE.adj == "UP", ])
CRE_DOWN_proteins <- rownames(results_list$Training_CRE[
  results_list$Training_CRE$DE.adj == "DOWN", ])
PLA_UP_proteins <- rownames(results_list$Training_PLA[
  results_list$Training_PLA$DE.adj == "UP", ])
PLA_DOWN_proteins <- rownames(results_list$Training_PLA[
  results_list$Training_PLA$DE.adj == "DOWN", ])

# Create Venn diagram
venn_cre <- list(
  "CRE UP"   = CRE_UP_proteins,
  "CRE DOWN" = CRE_DOWN_proteins,
  "PLA UP"   = PLA_UP_proteins,
  "PLA DOWN" = PLA_DOWN_proteins
)

p.v1 <- ggVennDiagram(venn_cre) +
  scale_fill_gradient(low = "white", high = "steelblue")

ggsave(file.path(report_dir, "venn_CREvPLA_proteins.png"),
       plot = p.v1,
       width = 10,
       height = 8,
       dpi = 300)

venn_surv <- list(
  "SURV_T1 UP"   = as.character(SURV_T1_UP_proteins),
  "SURV_T1 DOWN" = as.character(SURV_T1_DOWN_proteins),
  "SURV_T2 UP"   = as.character(SURV_T2_UP_proteins),
  "SURV_T2 DOWN" = as.character(SURV_T2_DOWN_proteins)
)

p.v2 <- ggVennDiagram(venn_surv) +
  scale_fill_gradient(low = "white", high = "steelblue")

ggsave(file.path(report_dir, "venn_SURVvCTL_proteins.png"),
       plot = p.v2,
       width = 10,
       height = 8,
       dpi = 300)

# --- Venn Breakdown Table ---------------------------------
# SURV v CTL
only_UP_T1 <- setdiff(SURV_T1_UP_proteins, SURV_T2_UP_proteins)
only_UP_T2 <- setdiff(SURV_T2_UP_proteins, SURV_T1_UP_proteins)
shared_UP <- intersect(SURV_T1_UP_proteins, SURV_T2_UP_proteins)

only_DOWN_T1 <- setdiff(SURV_T1_DOWN_proteins, SURV_T2_DOWN_proteins)
only_DOWN_T2 <- setdiff(SURV_T2_DOWN_proteins, SURV_T1_DOWN_proteins)
shared_DOWN <- intersect(SURV_T1_DOWN_proteins, SURV_T2_DOWN_proteins)

venn_table <- data.frame(
  Category = c("only_UP_T1", "only_UP_T2", "shared_UP",
               "only_DOWN_T1", "only_DOWN_T2", "shared_DOWN"),
  Count = c(length(only_UP_T1), length(only_UP_T2), length(shared_UP),
            length(only_DOWN_T1), length(only_DOWN_T2), length(shared_DOWN))
)

# SURV v CTL - Differential Protein Lists

SURV_T1_specific_UP <- data.frame(
  uniprot_id = c(only_UP_T1)
)

# Left join annotation data to unique proteins
SURV_T1_specific_UP <- left_join(SURV_T1_specific_UP,
                                 annot, 
                                 by = "uniprot_id")

# --- GSEA -----------------------------------
## SURV_T1 v CTL
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"

# Step 1: Create a ranked gene list
gene_list <- as.data.frame(results_list$Baseline_SURV_CTL)
# Left join HPAm ENSEMBL ids 
gene_list <- left_join(gene_list,
                       HPAm, 
                       by = c("uniprot" = "Uniprot"))

original_gene_list <- gene_list$logFC
names(original_gene_list) <- gene_list$Ensembl

# Left join HPAm ENSEMBL ids 
gene_list <- left_join(gene_list,
                       HPAm, 
                       by = c("uniprot" = "Uniprot"))

# Rank by log2FC * -log10(p-value) or just log2FC
ranked_genes_SURV_T1_vs_CTL <- results_list$Baseline_SURV_CTL %>%
  mutate(rank_metric = logFC * -log10(P.Value)) %>%
  arrange(desc(rank_metric))

# Add gene symbol column to ranked genes
ranked_genes_SURV_T1_vs_CTL$gene_symbol <- annot$gene[
  match(rownames(ranked_genes_SURV_T1_vs_CTL), annot$uniprot_id)]

# Create named vector (gene symbols as names, metric as values)
gene_list <- ranked_genes_SURV_T1_vs_CTL$rank_metric
names(gene_list) <- ranked_genes_SURV_T1_vs_CTL$gene_symbol  # or protein_id

# Remove duplicates and NAs
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- gene_list[!is.na(gene_list)]

# Step 2: Run GSEA
gsea_SURV_T1_vs_CTL <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",  # Biological Process (or "MF", "CC", "ALL")
      keyType = "SYMBOL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = TRUE
)

# View results
head(gsea_SURV_T1_vs_CTL@result)

write.csv(gsea_SURV_T1_vs_CTL@result, file.path(report_dir,"gsea_SURV_T1_vs_CTL.csv"))

pdf(file.path(report_dir,"gsea_dotplot_SURV_T1_vs_CTL.pdf"), width = 10, height = 8)
dotplot(gsea_SURV_T1_vs_CTL, showCategory = 20)
dev.off()


dotplot(gsea_SURV_T1_vs_CTL, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
gseaplot2(gsea_SURV_T1_vs_CTL, geneSetID = 1, pvalue_table = T)
ridgeplot(gsea_SURV_T1_vs_CTL, showCategory = 10) + labs(x= "enrichment distribution")

## SURV_T2 v CTL
# Step 1: Create a ranked gene list
# Rank by log2FC * -log10(p-value) or just log2FC
ranked_genes_SURV_T2_vs_CTL <- results_list$Training_SURV_CTL %>%
  mutate(rank_metric = logFC * -log10(P.Value)) %>%
  arrange(desc(rank_metric))

# Add gene symbol column to ranked genes
ranked_genes_SURV_T2_vs_CTL$gene_symbol <- annot$gene[
  match(rownames(ranked_genes_SURV_T2_vs_CTL), annot$uniprot_id)]

# Create named vector (gene symbols as names, metric as values)
gene_list2 <- ranked_genes_SURV_T2_vs_CTL$rank_metric
names(gene_list2) <- ranked_genes_SURV_T2_vs_CTL$gene_symbol  # or protein_id

# Remove duplicates and NAs
gene_list2 <- gene_list2[!duplicated(names(gene_list))]
gene_list2 <- gene_list2[!is.na(gene_list)]

# Step 2: Run GSEA
gsea_SURV_T2_vs_CTL <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological Process (or "MF", "CC", "ALL")
  keyType = "SYMBOL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

# View results
head(gsea_SURV_T2_vs_CTL@result)

write.csv(gsea_SURV_T2_vs_CTL@result, file.path(report_dir,"gsea_SURV_T2_vs_CTL.csv"))

pdf(file.path(report_dir,"gsea_dotplot_SURV_T2_vs_CTL.pdf"), width = 10, height = 8)
dotplot(gsea_SURV_T2_vs_CTL, showCategory = 20)
dev.off()


dotplot(gsea_SURV_T2_vs_CTL, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
gseaplot2(gsea_SURV_T2_vs_CTL, geneSetID = 1, pvalue_table = T)
ridgeplot(gsea_SURV_T2_vs_CTL, showCategory = 10) + labs(x= "enrichment distribution")

# --- KEGG ---------------------

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = gene_list[gene_list$Ensembl %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, 
        title = "KEGG Enriched Pathways", 
        split=".sign") + facet_grid(.~.sign)


