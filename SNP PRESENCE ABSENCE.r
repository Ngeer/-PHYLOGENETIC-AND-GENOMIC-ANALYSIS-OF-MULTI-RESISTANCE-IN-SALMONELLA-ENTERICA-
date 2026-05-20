

# ================================================================
#  SNP-Type Pattern Heatmap  — Publication Ready (v2)
#  Salmonella isolates × Efflux / Resistance genes
#  Requires: ggplot2, dplyr, tidyr, readxl, ggtext
# ================================================================

# ── 0. Install / load packages ───────────────────────────────
needed <- c("ggplot2", "dplyr", "tidyr", "readxl", "ggtext")
for (p in needed) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(ggtext)

# ── 1. Read & clean ──────────────────────────────────────────
setwd("C:\\Users\\Stuart\\Downloads")          # ← change if needed

df <- read_excel("SNP-1.xlsx") |>
  rename(Gene = Gene, SNP_raw = `SNP-Type`, Strain = Strain) |>
  filter(
    !is.na(Gene), Gene != "Gene",
    !is.na(Strain), Strain != "SNP-Type",
    !is.na(SNP_raw)
  ) |>
  mutate(
    Gene   = trimws(as.character(Gene)),
    Strain = trimws(as.character(Strain)),
    # Map raw text → clean 4-level category
    SNP_cat = case_when(
      SNP_raw == "No SNP detected"             ~ "None",
      SNP_raw == "Synonymous"                  ~ "Synonymous",
      SNP_raw == "Non-synonymous"              ~ "Non-synonymous",
      SNP_raw == "Synonymous / Non-synonymous" ~ "Both",
      TRUE                                     ~ "None"          # safety catch-all
    )
  )

# ── 2. Complete grid (every Gene × Strain combination) ───────
#    This is what prevents white / NA cells in the heatmap
all_combinations <- expand.grid(
  Gene   = unique(df$Gene),
  Strain = unique(df$Strain),
  stringsAsFactors = FALSE
)

df_complete <- all_combinations |>
  left_join(df, by = c("Gene", "Strain")) |>
  mutate(SNP_cat = if_else(is.na(SNP_cat), "None", SNP_cat))

# ── 3. Numeric matrix for clustering ─────────────────────────
snp_levels  <- c("None", "Synonymous", "Non-synonymous", "Both")
snp_numeric <- c("None" = 0L, "Synonymous" = 1L,
                 "Non-synonymous" = 2L, "Both" = 3L)

wide <- df_complete |>
  mutate(SNP_num = snp_numeric[SNP_cat]) |>
  select(Gene, Strain, SNP_num) |>
  pivot_wider(names_from = Strain, values_from = SNP_num,
              values_fn  = max, values_fill = 0L)

mat           <- as.matrix(wide[, -1])
rownames(mat) <- wide$Gene

# ── 4. Hierarchical clustering (Ward's D2) ───────────────────
row_ord      <- hclust(dist(mat),    method = "ward.D2")$order
col_ord      <- hclust(dist(t(mat)), method = "ward.D2")$order
gene_order   <- rownames(mat)[row_ord]
strain_order <- colnames(mat)[col_ord]

# ── 5. Final long data with ordered factor levels ─────────────
plot_df <- df_complete |>
  mutate(
    Gene    = factor(Gene,    levels = gene_order),
    Strain  = factor(Strain,  levels = strain_order),
    SNP_cat = factor(SNP_cat, levels = snp_levels)
  )

# ── 6. Colour palette ────────────────────────────────────────
snp_colours <- c(
  "None"            = "#0d1b2a",   # deep navy — absence
  "Synonymous"      = "#1565C0",   # strong blue
  "Non-synonymous"  = "#E65100",   # deep amber-orange
  "Both"            = "#B71C1C"    # dark red
)

# ── 7. Build plot ────────────────────────────────────────────
p <- ggplot(plot_df, aes(x = Strain, y = Gene, fill = SNP_cat)) +

  geom_tile(colour = "#2a2a2a", linewidth = 0.3) +

  scale_fill_manual(
    values = snp_colours,
    name   = "SNP Type",
    labels = c("None", "Synonymous",
               "Non-synonymous", "Both"),
    guide  = guide_legend(
      title.position = "top",
      title.hjust    = 0.5,
      nrow           = 1,
      keywidth       = unit(1.3, "cm"),
      keyheight      = unit(0.55, "cm"),
      label.theme    = element_text(size = 11, colour = "#1a1a1a"),
      override.aes   = list(colour = "#2a2a2a", linewidth = 0.5)
    )
  ) +

  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand   = c(0, 0)) +

  labs(
    title    = "SNP-Type Distribution Across Salmonella Isolates",
    subtitle = "Pattern of synonymous and non-synonymous SNPs in resistance genes",
    x        = "Isolate",
    y        = "Gene"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    # ---- Titles
    plot.title    = element_markdown(
      size   = 18, face = "bold", hjust = 0.5,
      colour = "#0d1b2a", margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 16, hjust = 0.5, colour = "#555555",face = "bold",
      margin = margin(b = 16)
    ),
    plot.margin   = margin(20, 24, 14, 14),

    # ---- Axis text
    axis.text.x = element_text(
      angle = 45, vjust = 1, hjust = 1,
      size  = 10, colour = "#1a1a1a"
    ),
    axis.text.y = element_text(
      face = "italic", size = 10, colour = "#1a1a1a"
    ),
    axis.title.x = element_text(
      size = 12, face = "bold", colour = "#1a1a1a",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 12, face = "bold", colour = "#1a1a1a",
      margin = margin(r = 10)
    ),
    axis.ticks = element_blank(),

    # ---- Panel
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),

    # ---- Legend
    legend.position    = "top",
    legend.title       = element_text(size = 11, face = "bold", colour = "#1a1a1a"),
    legend.margin      = margin(b = 8),
    legend.box.spacing = unit(0.15, "cm")
  )

# ── 8. Save ───────────────────────────────────────────────────
ggsave(
  filename = "snp_heatmap_pub_v2.png",
  plot     = p,
  width    = 18,
  height   = 13,
  dpi      = 600,
  units    = "in",
  bg       = "white"
)

cat("Done! Saved: snp_heatmap_pub_v2.png\n")











































