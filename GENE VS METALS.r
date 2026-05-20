# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

setwd("C:\\Users\\Stuart\\Downloads")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
# The file has a fake header row 0 ("Gene","Gene",...) and the real
# gene names are in row 1. We skip row 0 by reading with skip=1.
df_raw <- read.csv("GENE VS METAL.csv",
                   fileEncoding = "UTF-16",
                   check.names  = FALSE,
                   sep          = "\t",
                   skip         = 1)       # skip the fake "Gene Gene Gene..." row

# ── 2. RENAME FIRST COLUMN AND DROP SUMMARY ROW ───────────────────────────────
colnames(df_raw)[1] <- "Metal"

# Row 1 (after skip) is the "Metals" summary row with blank Metal name — drop it
df_raw <- df_raw %>%
  filter(trimws(Metal) != "" & trimws(Metal) != "Metals")

# ── 3. CLEAN METAL NAMES — fix trailing commas and duplicates ─────────────────
df_raw$Metal <- trimws(df_raw$Metal)

df_raw$Metal <- dplyr::recode(df_raw$Metal,
  # Fix trailing commas
  "Gold (Au),"         = "Gold (Au)",
  "Zinc (Zn),"         = "Zinc (Zn)",
  # Fix leading space
  " Zinc (Zn)"         = "Zinc (Zn)"
)

# ── 4. MERGE DUPLICATES BY SUMMING ────────────────────────────────────────────
df_merged <- df_raw %>%
  group_by(Metal) %>%
  summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE)),
            .groups = "drop")

# ── 5. RESHAPE TO LONG FORMAT ─────────────────────────────────────────────────
df_long <- df_merged %>%
  pivot_longer(
    cols      = -Metal,
    names_to  = "Gene",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count), Count > 0)

# Check
cat("Unique metals after merging:\n")
print(sort(unique(df_long$Metal)))
cat("\nUnique count values:", sort(unique(df_long$Count)), "\n")

# ── 6. AXIS ORDER ─────────────────────────────────────────────────────────────
metal_order <- sort(unique(df_long$Metal))
gene_order  <- rev(sort(unique(df_long$Gene)))

df_plot <- df_long %>%
  mutate(
    Metal = factor(Metal, levels = metal_order),
    Gene  = factor(Gene,  levels = gene_order)
  )

# ── 7. BUBBLE PLOT ────────────────────────────────────────────────────────────
p <- ggplot(df_plot, aes(x = Metal, y = Gene)) +

  geom_point(aes(size = Count, colour = Count),
             alpha = 0.92) +

  scale_colour_gradientn(
    colours = c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"),
    name    = "No. of\ndetections",
    limits  = c(min(df_long$Count), max(df_long$Count)),
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
    title    = "AMR Gene vs Metal",
    subtitle = "The plot shows the number of strains (indicated by bubble size and colour) that have a given gene-metal association",
    x        = "Metal / Metal Combination",
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

    axis.text.x  = element_text(colour = "black", size = 12,
                                 angle = 40, hjust = 1, vjust = 1),
    axis.title.x = element_text(colour = "black", size = 13, face = "bold",
                                 margin = margin(t = 14)),

    axis.text.y  = element_text(colour = "black", size = 13, face = "italic"),
    axis.title.y = element_text(colour = "black", size = 13, face = "bold",
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

# ── 8. SAVE ───────────────────────────────────────────────────────────────────
ggsave(
  filename  = "bubble_gene_vs_metal.png",
  plot      = p,
  width     = 18,
  height    = 16,
  dpi       = 600,
  bg        = "white",
  limitsize = FALSE
)

message("Done — saved as bubble_gene_vs_metal.png")