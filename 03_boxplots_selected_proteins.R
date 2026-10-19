# ======================================================
# 03_boxplots_selected_proteins.R
# Boxplots – Selected Olink inflammation proteins
# CP vs Typically Developing controls
# ======================================================

# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# ---- Assumes cleaned data already loaded ----
# source("00_load_and_clean_olink.R")

# ---- Proteins to plot (exact Olink names) ----
boxplot_proteins <- c(
  "TWEAK",
  "IL-10RB",
  "uPA",
  "IL-17C"
)

# ---- Prepare data (wide → long) ----
boxplot_data <- olink_wide %>%
  select(SampleID, group, all_of(boxplot_proteins)) %>%
  pivot_longer(
    cols      = all_of(boxplot_proteins),
    names_to  = "protein",
    values_to = "NPX"
  ) %>%
  mutate(
    group   = factor(group, levels = c("TD", "CP")),
    protein = factor(
      protein,
      levels = c("TWEAK", "IL-10RB", "uPA", "IL-17C")
    )
  )

# ---- Boxplot ----
p_box <- ggplot(boxplot_data, aes(x = group, y = NPX, fill = group)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = group),
    width  = 0.12,
    size   = 2.5,
    alpha  = 0.75
  ) +
  facet_wrap(~ protein, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("TD" = "#4DB5BC", "CP" = "#870052")
  ) +
  scale_color_manual(
    values = c("TD" = "#4DB5BC", "CP" = "#870052")
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.10)
  ) +
  labs(
    title = "Selected Olink inflammation proteins",
    x     = NULL,
    y     = "NPX"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    axis.title.y     = element_text(size = 22),
    strip.text       = element_text(size = 18),
    panel.spacing.x  = unit(1.5, "lines"),
    panel.spacing.y  = unit(1.2, "lines"),
    legend.position  = "none"
  )

print(p_box)