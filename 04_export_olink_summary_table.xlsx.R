# ======================================================
# 04_export_olink_summary_table.xlsx.R
# Export protein-level Olink summary table (CP vs TD)
# ======================================================

# ---- Packages ----
library(dplyr)
library(tidyr)
library(writexl)

# ---- Assumes these objects already exist ----
# olink_wide   : cleaned wide matrix (SampleID, group, proteins)
# volcano_df  : results with protein, log2fc, padj
# source("00_load_and_clean_olink.R")
# source("02_volcano_CP_vs_TD_olink.R")

# ---- Ensure correct group order ----
olink_wide <- olink_wide %>%
  mutate(group = factor(group, levels = c("TD", "CP")))

# ---- Define sample order: TD first, then CP ----
sample_order <- olink_wide %>%
  distinct(SampleID, group) %>%
  arrange(group, SampleID) %>%
  pull(SampleID)

# ---- Identify protein columns ----
protein_cols <- setdiff(colnames(olink_wide), c("SampleID", "group"))

# ---- Build NPX table: one row per protein, columns = individuals ----
npx_by_protein <- olink_wide %>%
  select(SampleID, all_of(protein_cols)) %>%
  pivot_longer(
    cols      = all_of(protein_cols),
    names_to  = "protein",
    values_to = "NPX"
  ) %>%
  mutate(SampleID = factor(SampleID, levels = sample_order)) %>%
  arrange(protein, SampleID) %>%
  distinct() %>%
  pivot_wider(
    names_from  = SampleID,
    values_from = NPX
  )

# ---- Keep unique volcano stats ----
summary_stats <- volcano_df %>%
  select(protein, log2fc, padj) %>%
  distinct(protein, .keep_all = TRUE)

# ---- Combine stats + NPX values ----
final_table <- summary_stats %>%
  left_join(npx_by_protein, by = "protein") %>%
  arrange(protein)

# ---- Final sanity check (optional) ----
# Ensure one row per protein
stopifnot(
  nrow(final_table) == length(unique(final_table$protein))
)

# ---- Export to Excel ----
out_dir <- "/Users/sebastian.edman/Desktop/CP Olink/results"
dir.create(out_dir, showWarnings = FALSE)

write_xlsx(
  final_table,
  path = file.path(out_dir, "Olink_CP_vs_TD_summary.xlsx")
)