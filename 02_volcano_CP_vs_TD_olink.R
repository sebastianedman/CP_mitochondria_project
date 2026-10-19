# ======================================================
# 02_volcano_CP_vs_TD_olink.R
# Volcano plot – Olink Target 96 Inflammation
# CP vs Typically Developing controls
# Up-regulated = higher in CP
# Significance uses BOTH padj and |log2FC| cutoffs
# ======================================================

# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(purrr)

# ---- Load cleaned data ----
# This assumes you have already run:
# source("00_load_and_clean_olink.R")

# Make sure group factor is correctly ordered
olink_wide <- olink_wide %>%
  mutate(group = factor(group, levels = c("TD", "CP")))

# ---- Metadata and protein matrix ----
meta <- olink_wide %>%
  select(SampleID, group)

olink_mat <- olink_wide %>%
  select(-SampleID, -group)

# ---- Differential analysis ----
padj_cutoff   <- 0.05
log2fc_cutoff <- 0.58   # ~1.5x change. Use 1.0 for ~2x if you prefer.

volcano_df <- purrr::map_dfr(colnames(olink_mat), function(prot) {
  
  df <- tibble(
    NPX   = olink_mat[[prot]],
    group = meta$group
  ) %>%
    filter(!is.na(NPX))
  
  if (nlevels(df$group) != 2) {
    return(tibble(
      protein = prot,
      log2fc  = NA_real_,
      p_value = NA_real_
    ))
  }
  
  tt <- t.test(NPX ~ group, data = df)
  
  means <- tapply(df$NPX, df$group, mean, na.rm = TRUE)
  log2fc <- means["CP"] - means["TD"]
  
  tibble(
    protein = prot,
    log2fc  = as.numeric(log2fc),
    p_value = tt$p.value
  )
}) %>%
  mutate(
    padj           = p.adjust(p_value, method = "BH"),
    neg_log10_padj = -log10(padj),
    # ---- UPDATED: direction uses BOTH padj and |log2FC| cutoffs ----
    direction = case_when(
      padj < padj_cutoff & log2fc >=  log2fc_cutoff ~ "Up",    # higher in CP
      padj < padj_cutoff & log2fc <= -log2fc_cutoff ~ "Down",  # lower in CP
      TRUE                                           ~ "Not significant"
    )
  )

# ---- Volcano plot ----
p_volcano <- ggplot(volcano_df, aes(x = log2fc, y = neg_log10_padj)) +
  geom_point(aes(color = direction), alpha = 0.65, size = 4.5, na.rm = TRUE) +
  
  # ---- padj cutoff line ----
geom_hline(
  yintercept = -log10(padj_cutoff),
  linetype  = "dotted",
  linewidth = 1.2,
  color     = "black",
  lineend   = "round"
) +
  
  # ---- NEW: log2FC cutoff lines ----
geom_vline(
  xintercept = c(-log2fc_cutoff, log2fc_cutoff),
  linetype   = "dotted",
  linewidth  = 1.2,
  color      = "black",
  lineend    = "round"
) +
  
  scale_color_manual(
    values = c(
      "Up"              = "#870052",
      "Down"            = "#4DB5BC",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    title = "Volcano plot: Olink Inflammation (CP vs TD)",
    x     = "log2 Fold Change (CP − TD)",
    y     = "-log10(adjusted p)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major  = element_line(color = "grey90", linewidth = 1.0),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.background  = element_rect(fill = NA, color = NA),
    plot.background   = element_rect(fill = NA, color = NA),
    axis.ticks        = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.30, "cm"),
    axis.text.x       = element_text(color = "black", size = 24),
    axis.text.y       = element_text(color = "black", size = 24),
    axis.title.x      = element_text(color = "black", size = 22),
    axis.title.y      = element_text(color = "black", size = 22)
  )

# ---- Optional: label top significant proteins ----
label_top     <- TRUE
top_n_labels  <- 15

if (label_top) {
  
  # ---- UPDATED: labels match the new definition (Up/Down only) ----
  top_to_label <- volcano_df %>%
    filter(direction != "Not significant") %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = top_n_labels)
  
  p_volcano <- p_volcano +
    ggrepel::geom_label_repel(
      data = top_to_label,
      aes(label = protein),
      
      size = 4.8,
      fontface = "plain",
      
      box.padding   = 0.75,
      point.padding = 0.95,
      max.overlaps  = 10,
      
      label.size = 1.0,
      label.r = unit(0.45, "lines"),
      fill = "white",
      color = "black",
      label.colour = "black",
      
      na.rm = TRUE
    )
}

print(p_volcano)