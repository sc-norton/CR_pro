#--- Variable Regressions ------------------

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
input_dir <- file.path(base_dir, "02_normalization", "c_data")
report_dir <- file.path(base_dir, "05_analysis", "b_results" , "regressions")

dir.create(report_dir, showWarnings = F, recursive = T)

# Load data
SURV <- readRDS(file.path(input_dir, "01_normalized_DAList_SURV.RDS"))

# Subset data
BL <- filter_samples(SURV, SURV$metadata$cancer_time == "SURV_T1")
RT <- filter_samples(SURV, SURV$metadata$cancer_time == "SURV_T2")

# --- Fat Pct (BL SURV) ------
# Ensure numeric
class(BL$metadata$DXA_percFat)

# Create design
design_fat <- model.matrix(~ DXA_percFat + age + sex + BMI, 
                           data = BL$metadata)

# Fit the model
mat <- as.matrix(BL$data)

fit_fat <- lmFit(mat, design_fat) |>
  eBayes(robust = TRUE, trend = TRUE)

results_fat <- topTable(fit_fat,
                       coef    = "DXA_percFat",
                       number  = Inf,
                       sort.by = "P") |>
  tibble::rownames_to_column("uniprot_id")

head(results_fat)

sig_proteins_fat <- results_fat %>%
  filter(P.Value < 0.05, abs(logFC) > 0)
head(sig_proteins_fat)

top10_fat <- sig_proteins_fat %>%
  arrange(P.Value) %>%
  slice_head(n = 10)

head(top10)
# Long format with metadata
mat_long <- as.data.frame(t(as.matrix(BL$data))) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(BL$metadata, by = "sample_id") %>%
  pivot_longer(cols      = all_of(top10_fat$uniprot_id),
               names_to  = "uniprot_id",
               values_to = "intensity") %>%
  left_join(top10_fat %>% dplyr::select(uniprot_id, logFC, P.Value), by = "uniprot_id") %>%
  left_join(BL$annotation %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
  mutate(label = paste0(gene, "  |  β=", round(logFC, 2), 
                        "  p=", signif(P.Value, 2)),
         label = factor(label, levels = unique(label[order(P.Value)])))  

ggplot(mat_long, aes(x = DXA_percFat, y = intensity)) +
  geom_point(alpha = 0.7, color = "#E96A2C") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~ label, scales = "free_y") +          
  labs(x     = "DXA Percent Fat (%)",
       y     = "Normalized Intensity",
       title = "Top 10 Proteins Associated with Body Fat %") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 7))


# --- Appendicular Lean Mass Index (BL SURV) -------
# Ensure numeric
class(BL$metadata$DXA_ALM_per_m2)

# Create design
design_almi <- model.matrix(~ DXA_ALM_per_m2 + age + sex + BMI, 
                           data = BL$metadata)

# Fit the model
mat <- as.matrix(BL$data)

fit_almi <- lmFit(mat, design_almi) |>
  eBayes(robust = TRUE, trend = TRUE)

results_almi <- topTable(fit_almi,
                        coef    = "DXA_ALM_per_m2",
                        number  = Inf,
                        sort.by = "P") |>
  tibble::rownames_to_column("uniprot_id")

head(results_almi)

sig_proteins_almi <- results_almi %>%
  filter(P.Value < 0.05, abs(logFC) > 0)
head(sig_proteins_almi)

top10_almi <- sig_proteins_almi %>%
  arrange(P.Value) %>%
  slice_head(n = 10)

head(top10_almi)

# Long format with metadata
mat_long_almi <- as.data.frame(t(as.matrix(BL$data))) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(BL$metadata, by = "sample_id") %>%
  pivot_longer(cols      = all_of(top10_almi$uniprot_id),
               names_to  = "uniprot_id",
               values_to = "intensity") %>%
  left_join(top10_almi %>% dplyr::select(uniprot_id, logFC, P.Value), by = "uniprot_id") %>%
  left_join(BL$annotation %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
  mutate(label = paste0(gene, "  |  β=", round(logFC, 2), 
                        "  p=", signif(P.Value, 2)),
         label = factor(label, levels = unique(label[order(P.Value)])))  

ggplot(mat_long_almi, aes(x = DXA_ALM_per_m2, y = intensity)) +
  geom_point(alpha = 0.7, color = "#E96A2C") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~ label, scales = "free_y") +          
  labs(x     = "DXA ALM Index (kg/m2)",
       y     = "Normalized Intensity",
       title = "Top 10 Proteins Associated with ALMI") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 7))
