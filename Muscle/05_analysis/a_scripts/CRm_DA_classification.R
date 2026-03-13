# --- GSEA --------------
# --- 0: Setup ----------------
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


cat("\n Base directory: ", base_dir, "\n")
# ============================================================
#  GO Biological Process Pie Chart — Upregulated Proteins
#  Organism: Human | IDs: UniProt
# ============================================================
#
#  BEFORE RUNNING:
#  1. Go to https://www.uniprot.org/id-mapping
#  2. Paste UniProt IDs → Map to UniProtKB
#  3. Download as CSV, rename and set the path below
#
# ============================================================


# ---- 1. DEFINE ALL DATASETS --------------------------------

datasets <- list(
  list(
    input_file  = file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_UP_mapped.csv"),
    output_file = file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_UP_chart.png"),
    title       = "Protein Categories - Higher Abundance in SURV at BL",
    subtitle    = "232 proteins UP at BL vs CTL"
  ),
  list(
    input_file  = file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_DOWN_mapped.csv"),
    output_file = file.path(report_dir, "DA_protein_lists", "Baseline_SURVvCTL_DOWN_chart.png"),
    title       = "Protein Categories - Lower Abundance in SURV at BL",
    subtitle    = "64 proteins DOWN at BL vs CTL"
  ),
  list(
    input_file  = file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_UP_mapped.csv"),
    output_file = file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_UP_chart.png"),
    title       = "Protein Categories - Higher Abundance in SURV after RT",
    subtitle    = "303 proteins UP after RT vs CTL"
  ),
  list(
    input_file  = file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_DOWN_mapped.csv"),
    output_file = file.path(report_dir, "DA_protein_lists", "Training_SURVvCTL_DOWN_chart.png"),
    title       = "Protein Categories - Lower Abundance in SURV after RT",
    subtitle    = "118 proteins DOWN after RT vs CTL"
  )
)

# Column names from the UniProt CSV export (check your file header if unsure)
uniprot_col <- "Entry"                             # UniProt ID column
go_bp       <- "Gene Ontology (biological process)"  # GO BP column

# ---- 2. BROAD CATEGORY MAPPING -----------------------------
# Each category is matched by keywords against GO BP term text.
# Add, remove, or rename categories here to suit your biology.

category_map <- list(
  "Metabolism"               = c("metabol", "biosynthes", "catabo", "lipid", "fatty acid",
                                 "glucos", "glycol", "oxidat", "tca cycle", "citric acid",
                                 "amino acid", "nucleotide synthesis", "energy homeostasis"),
  "Cell Cycle & Division"    = c("cell cycle", "mitosis", "meiosis", "cytokinesis",
                                 "chromosome segreg", "cell division"),
  "DNA Repair & Replication" = c("dna repair", "dna damage", "dna replicat",
                                 "nucleotide excision", "mismatch repair", "homologous recomb"),
  "Gene Expression"          = c("transcription", "mrna", "rna process", "splicing",
                                 "translation", "ribos", "chromatin remodel", "epigenet",
                                 "histone modif", "gene express", "rna polymerase"),
  "Signaling"                = c("signal transduct", "phosphorylat", "kinase cascade",
                                 "second messenger", "mapk", "pi3k", "wnt", "notch",
                                 "hedgehog", "jak-stat", "nf-kb", "camp", "calcium signal"),
  "Immune Response"          = c("immune", "inflammat", "cytokine", "interferon",
                                 "interleukin", "antigen", "innate immun", "adaptive immun",
                                 "lymphocyte", "nk cell", "toll-like receptor"),
  "Apoptosis & Cell Death"   = c("apoptos", "cell death", "necrosis", "autophagy",
                                 "caspase", "programmed cell"),
  "Transport"                = c("transport", "vesicle", "exocytosis", "endocytosis",
                                 "protein secretion", "ion transport", "membrane trafficking",
                                 "nuclear import", "nuclear export"),
  "Cell Adhesion & Migration" = c("cell adhesion", "cell migration", "cell motility",
                                  "cytoskeleton", "actin", "integrin", "extracellular matrix"),
  "Development"              = c("development", "differentiation", "morphogenesis",
                                 "embryo", "organogenesis", "stem cell", "cell fate"),
  "Stress Response"          = c("stress response", "heat shock", "unfolded protein",
                                 "oxidative stress", "hypoxia", "protein folding", "chaperone"),
  "Protein Degradation"      = c("proteasom", "ubiquitin", "protein degradat",
                                 "protein cataboli", "protein turnover")
)

assign_category <- function(go_string) {
  if (is.na(go_string) || str_trim(go_string) == "") return("Unclassified")
  
  go_lower <- tolower(go_string)
  
  for (cat in names(category_map)) {
    if (any(sapply(category_map[[cat]], function(kw) grepl(kw, go_lower, fixed = TRUE)))) {
      return(cat)
    }
  }
  return("Unclassified")
}


category_colours <- c(
  "Metabolism"               = "#4E9A6F",
  "Cell Cycle & Division"    = "#E07B54",
  "DNA Repair & Replication" = "#5B8DB8",
  "Gene Expression"          = "#A66DB8",
  "Signaling"                = "#E0B83A",
  "Immune Response"          = "#D94F4F",
  "Apoptosis & Cell Death"   = "#7A7A7A",
  "Transport"                = "#4ABFBF",
  "Cell Adhesion & Migration"= "#B8A96D",
  "Development"              = "#6DAF6D",
  "Stress Response"          = "#CF7A3A",
  "Protein Degradation"      = "#8C6BAE",
  "Unclassified"             = "#CCCCCC"
)

# ---- 3. Create pie chart for each dataset ---------------------
for (ds in datasets) {
  cat("\n Processing:", ds$input_file, "\n")
  
  df <- read_csv(ds$input_file)
  cat("Proteins loaded:", nrow(df), "\n")
  
  df <- df %>%
    mutate(Category = sapply(.data[[go_bp]], assign_category))
  
  summary_df <- df %>%
    count(Category, name = "Count") %>%
    arrange(desc(Count)) %>%
    mutate(
      Percentage = round(Count / sum(Count) * 100, 1),
      Label = paste0(Category, "\n", Count, " (", Percentage, "%)")
    )
  
  cat("Category breakdown:\n")
  print(summary_df %>% dplyr::select(Category, Count, Percentage))
  
  used_colours <- category_colours[names(category_colours) %in% summary_df$Category]
  
  pie_chart <- ggplot(summary_df, aes(x = factor(1), y = Count, fill = Category)) +
    geom_col(width = 1, colour = "white", linewidth = 0.1) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = used_colours,
      breaks = summary_df$Category,
      labels = setNames(
        paste0(summary_df$Category, " | n=", summary_df$Count, " (", summary_df$Percentage, "%)"),
        summary_df$Category
      )
    ) +
    labs(title = ds$title, subtitle = ds$subtitle) +
    theme_void()
  
  ggsave(
    filename = ds$output_file,
    plot     = pie_chart,
    width    = 9,
    height   = 7,
    dpi      = 300,
    bg       = "white"
  )
  
  cat("Saved:", ds$output_file, "\n")
}

cat("\nAll datasets processed.\n")

# --- UP : Identify proteins only in one list ----

df1 <- datasets[[1]]$input_file
df2 <- datasets[[3]]$input_file

df1 <- read.csv(df1)
df2 <- read.csv(df2)

df1$source <- "BL"
df2$source <- "RT"

combined <- bind_rows(df1, df2) %>%
  group_by(Entry) %>%
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

UP_BL_only <- UP_BL_only %>%
  mutate(Category = sapply(UP_BL_only$Gene.Ontology..biological.process., 
                           assign_category))

UP_RT_only <- UP_RT_only %>%
  mutate(Category = sapply(UP_RT_only$Gene.Ontology..biological.process., 
                           assign_category))

summary_df <- UP_BL_only %>%
  count(Category, name = "Count") %>%
  arrange(desc(Count)) %>%
  mutate(
    Percentage = round(Count / sum(Count) * 100, 1),
    Label = paste0(Category, "\n", Count, " (", Percentage, "%)")
  )

cat("Category breakdown:\n")
print(summary_df %>% dplyr::select(Category, Count, Percentage))

