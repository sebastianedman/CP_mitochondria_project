# ======================================================
# 00_load_and_clean_olink.R
# Load and clean Olink Target 96 Inflammation data
# CP vs Typically Developing controls
# ======================================================

# ---- Packages ----
library(readr)
library(dplyr)
library(tidyr)

# ---- File paths ----
data_dir <- "/Users/sebastian.edman/Desktop/CP Olink"
file_name <- "2025058_CP-OLINK_EP_Standard_2025-12-02.csv"

# ---- Load raw Olink data (semicolon separated) ----
olink_raw <- readr::read_delim(
  file.path(data_dir, file_name),
  delim = ";"
)

# ---- Fix sample naming error ----
# Former EP215  -> EP214
# Former EP215_2 -> EP215 (rerun retained)
olink <- olink_raw %>%
  mutate(
    SampleID = case_when(
      SampleID == "EP215"   ~ "EP214",
      SampleID == "EP215_2" ~ "EP215",
      TRUE                  ~ SampleID
    )
  )

# ---- Remove Olink internal control samples ----
olink_control_ids <- c("Sample_12", "Sample_24")

olink <- olink %>%
  filter(!SampleID %in% olink_control_ids)

# ---- Define Typically Developing (TD) samples ----
td_ids <- c(
  "EP206", "EP207", "EP210", "EP212",
  "EP218", "EP219",
  "EP223", "EP224", "EP225", "EP226",
  "EP228", "EP229", "EP230", "EP231", "EP232", "EP233"
)

# ---- Assign group ----
olink <- olink %>%
  mutate(
    group = case_when(
      SampleID %in% td_ids ~ "TD",
      TRUE                 ~ "CP"
    )
  )

# ---- Sanity checks ----
# Number of unique samples
n_samples <- olink %>% distinct(SampleID) %>% nrow()

# Samples per group
group_counts <- olink %>%
  distinct(SampleID, group) %>%
  count(group)

print(n_samples)
print(group_counts)

# ---- Create wide NPX table (samples x proteins) ----
olink_wide <- olink %>%
  select(SampleID, group, Assay, NPX) %>%
  pivot_wider(
    names_from  = Assay,
    values_from = NPX
  )

# ---- Final check ----
dim(olink_wide)