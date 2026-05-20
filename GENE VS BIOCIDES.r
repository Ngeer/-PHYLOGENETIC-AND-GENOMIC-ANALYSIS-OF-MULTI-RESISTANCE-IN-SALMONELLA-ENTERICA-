# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

setwd("C:\\Users\\Stuart\\Downloads")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
df_raw <- read.csv("Sheet 2.csv",
                   fileEncoding = "UTF-16",
                   check.names  = FALSE,
                   sep          = "\t")

# ── 2. RENAME FIRST COLUMN ────────────────────────────────────────────────────
colnames(df_raw)[1] <- "Biocide"

# ── 3. CLEAN BIOCIDE NAMES — fix spelling variants and trailing commas ─────────
#
#  Duplicates found:
#  - "Acriflavin" and "Acriflavine"            → same biocide, merge
#  - "Methylene Blue, Crystal Violet"           → appears twice, merge
#  - "Sodium Deoxycholate (SDC)"               → appears 3 times (one trailing comma, one with HCl row kept separate)
#  - "Cyclohexane, Diphenyl Ether,"            → trailing comma, merge with clean version
#  - "Sodium Deoxycholate (SDC),"              → trailing comma duplicate
#
df_raw$Biocide <- trimws(df_raw$Biocide)

df_raw$Biocide <- dplyr::recode(df_raw$Biocide,
  # Spelling fix
  "Acriflavin"                    = "Acriflavine",

  # Trailing comma fixes
  "Cyclohexane, Diphenyl Ether,"  = "Cyclohexane, Diphenyl Ether",
  "Sodium Deoxycholate (SDC),"    = "Sodium Deoxycholate (SDC)",
  "Carbonyl cyanide 3-chlorophenylhydrazone (CCCP) ,"
                                  = "CCCP"
)

# Also shorten the very long names for readability
df_raw$Biocide <- dplyr::recode(df_raw$Biocide,
  "Carbonyl cyanide 3-chlorophenylhydrazone (CCCP)"
    = "CCCP",
  "Tetraphenylphosphonium (TPP), Sodium Deoxycholate (SDC), Ethidium Bromide, Benzylkonium Chloride (BAC), Acriflavine"
    = "TPP; SDC; Ethidium Bromide;\nBAC; Acriflavine",
  "Acriflavine, Phenol, Triclosan, p-xylene, Cyclohexane, Pentane"
    = "Acriflavine; Phenol;\nTriclosan; p-xylene; Cyclohexane",
  "Acriflavine, Sodium Dodecyl Sulfate (SDS), Sodium Deoxycholate (SDC), Tetraphenylphosphonium (TPP), Benzylkonium Chloride (BAC), Methyl Viologen, Ethidium Bromide"
    = "Acriflavine; SDS; SDC;\nTPP; BAC; Methyl Viologen",
  "Benzylkonium Chloride (BAC), Sodium Dodecyl Sulfate (SDS)"
    = "BAC; SDS",
  "Cyclohexane, Diphenyl Ether"
    = "Cyclohexane; Diphenyl Ether",
  "Cyclohexane, Diphenyl Ether, n-hexane"
    = "Cyclohexane; Diphenyl Ether;\nn-hexane",
  "Cyclohexane, Pentane, n-hexane, Diphenyl Ether"
    = "Cyclohexane; Pentane;\nn-hexane; Diphenyl Ether",
  "Hydrogen Peroxide (H2O2), Benzylkonium Chloride (BAC), Chlorhexidine"
    = "H2O2; BAC; Chlorhexidine",
  "Methylene Blue, Crystal Violet"
    = "Methylene Blue;\nCrystal Violet",
  "Phenylmercury Acetate, 2-Chlorophenylhydrazine, Carbonylcyanide m-chlorophenyl hydrazone (CCCP), Tetrachlorosalicylanilide (TCS)"
    = "Phenylmercury Acetate;\n2-Chlorophenylhydrazine; CCCP; TCS",
  "Sodium Deoxycholate (SDC)"
    = "Sodium Deoxycholate\n(SDC)",
  "Sodium Deoxycholate (SDC), Hydrochloric acid (HCl)"
    = "SDC; HCl",
  "Sodium Dodecyl Sulfate (SDS), Sodium Deoxycholate (SDC)"
    = "SDS; SDC"
)

# ── 4. MERGE DUPLICATES BY SUMMING COUNTS ─────────────────────────────────────
#
#  After renaming, rows with the same Biocide label are now true duplicates.
#  We sum their gene counts so no information is lost.
#
df_merged <- df_raw %>%
  group_by(Biocide) %>%
  summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE)),
            .groups = "drop")

# ── 5. RESHAPE TO LONG FORMAT ─────────────────────────────────────────────────
df_long <- df_merged %>%
  pivot_longer(
    cols      = -Biocide,
    names_to  = "Gene",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count), Count > 0)

# ── 6. CHECK UNIQUE VALUES ────────────────────────────────────────────────────
cat("Unique count values:", sort(unique(df_long$Count)), "\n")
cat("Unique biocides after merging:\n")
print(unique(df_long$Biocide))

# ── 7. AXIS ORDER ─────────────────────────────────────────────────────────────
# X-axis: biocide labels alphabetical
biocide_order <- sort(unique(df_long$Biocide))

# Y-axis: genes alphabetical reversed (soxS top → acrA bottom)
gene_order <- rev(sort(unique(df_long$Gene)))

df_plot <- df_long %>%
  mutate(
    Biocide = factor(Biocide, levels = biocide_order),
    Gene    = factor(Gene,    levels = gene_order)
  )

# ── 8. BUBBLE PLOT ────────────────────────────────────────────────────────────
p <- ggplot(df_plot, aes(x = Biocide, y = Gene)) +

  geom_point(aes(size = Count, colour = Count),
             alpha = 0.92) +

  scale_colour_gradientn(
    colours = c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"),
    name    = "No. of\ndetections",
    limits  = c(1, max(df_long$Count)),
    guide   = guide_colorbar(
      barwidth = 2, barheight = 14,
      title.position = "top", title.hjust = 0.5,
      ticks.colour   = "grey60"
    )
  ) +

  scale_size_continuous(
    range  = c(5, 16),
    name   = "No. of\ndetections",
    breaks = sort(unique(df_long$Count)),
    labels = as.character(sort(unique(df_long$Count)))
  ) +

  guides(
    colour = guide_colorbar(
      barwidth = 2, barheight = 14,
      title.position = "top", title.hjust = 0.5
    ),
    size = guide_legend(
      override.aes   = list(colour = "#fd8d3c"),
      title.position = "top", title.hjust = 0.5
    )
  ) +

  labs(
    title    = "AMR Gene vs Biocide",
    subtitle = "The plot shows the number of strains (indicated by bubble size and colour) that have a given gene-biocide association",
    x        = "Biocide / Disinfectant Combination",
    y        = "Resistance Gene "
  ) +

  theme_minimal(base_size = 15) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = "grey80",
                                     linewidth = 0.3),
    legend.key        = element_rect(fill = "white", colour = NA),

    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_blank(),

    axis.text.x  = element_text(colour = "black", size = 14,
                                 angle = 45, hjust = 1, vjust = 1,
                                 lineheight = 0.85),
    axis.title.x = element_text(colour = "black", size = 15, face = "bold",
                                 margin = margin(t = 14)),

    axis.text.y  = element_text(colour = "black", size = 14, face = "italic"),
    axis.title.y = element_text(colour = "black", size = 15, face = "bold",
                                 angle = 90, margin = margin(r = 12)),

    axis.ticks = element_line(colour = "grey70", linewidth = 0.3),

    plot.title    = element_text(colour = "black", size = 20, face = "bold",
                                  hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = "grey30", size = 16,
                                  hjust = 0.5, margin = margin(b = 14)),
    plot.margin   = margin(20, 20, 20, 20),

    legend.text  = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12.5, face = "bold"),
    legend.position = "right",
    legend.box   = "vertical"
  )

# ── 9. SAVE ───────────────────────────────────────────────────────────────────
ggsave(
  filename  = "bubble_gene_vs_biocide.png",
  plot      = p,
  width     = 28,
  height    = 18,
  dpi       = 600,
  bg        = "white",
  limitsize = FALSE
)

message("Done — saved as bubble_gene_vs_biocide.png")