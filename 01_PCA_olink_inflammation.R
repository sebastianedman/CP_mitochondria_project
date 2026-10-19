# ======================================================
# PCA – Olink Target 96 Inflammation
# CP vs TD
# ======================================================

# ---- Load packages ----
library(dplyr)
library(ggplot2)
library(grid)

# ---- Assume olink_wide already exists ----
# (loaded + cleaned in earlier script)

# ---- Create scaled matrix ----
olink_mat <- olink_wide %>%
  select(-SampleID, -group) %>%
  as.matrix()

olink_mat_scaled <- scale(olink_mat)
rownames(olink_mat_scaled) <- olink_wide$SampleID

# ---- PCA ----
pca_res <- prcomp(olink_mat_scaled, center = FALSE, scale. = FALSE)

var_expl <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_var <- round(var_expl[1] * 100, 1)
pc2_var <- round(var_expl[2] * 100, 1)

pca_df <- data.frame(
  SampleID = olink_wide$SampleID,
  group    = olink_wide$group,
  PC1      = pca_res$x[, 1],
  PC2      = pca_res$x[, 2]
)

# ---- Plot ----
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 4.5, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("CP" = "#4DB5BC", "TD" = "#870052")) +
  labs(
    title = "PCA: Olink Target 96 Inflammation",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.text        = element_text(size = 24),
    axis.title       = element_text(size = 22)
  )

print(p_pca)