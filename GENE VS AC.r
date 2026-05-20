# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "tidyr", "ggplot2")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

setwd("C:\\Users\\Stuart\\Downloads")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
# New CSV: rows = antibiotics, columns = genes
df_raw <- read.csv("Sheet 1 (1).csv",
                   fileEncoding = "UTF-16",
                   check.names  = FALSE,
                   sep          = "\t")

# ── 2. RENAME FIRST COLUMN → "Antibiotic" ─────────────────────────────────────
colnames(df_raw)[1] <- "Antibiotic"

# ── 3. CLEAN THE ANTIBIOTIC ROW NAMES (fix garbage byte in row 6) ─────────────
df_raw$Antibiotic <- trimws(df_raw$Antibiotic)

# Fix the one row with a trailing garbage byte
df_raw$Antibiotic <- ifelse(
  startsWith(df_raw$Antibiotic,
             "ciprofloxacin;tigecycline;chloramphenicol;rifampin;tetracycline;ampicillin;cefalothin"),
  "ciprofloxacin;tigecycline;chloramphenicol;rifampin;tetracycline;ampicillin;cefalothin",
  df_raw$Antibiotic
)

# ── 4. MAP ANTIBIOTIC ROW NAMES → SHORT DISPLAY LABELS ────────────────────────
antibiotic_map <- c(
  "amikacin;gentamicin C;tobramycin"
    = "Amikacin; Gentamicin C;\nTobramycin",

  "aminocoumarin antibiotic, aminoglycoside antibiotic"
    = "Aminocoumarin +\nAminoglycoside",

  "bicyclomycin-like antibiotic"
    = "Bicyclomycin-like",

  "cefalothin;chloramphenicol;tigecycline;ampicillin;tetracycline;rifampin"
    = "Cefalothin; Chloramphenicol;\nTetracycline; Rifampin",

  "ciprofloxacin"
    = "Ciprofloxacin",

  "ciprofloxacin;tigecycline;chloramphenicol;rifampin;tetracycline;ampicillin;cefalothin"
    = "Ciprofloxacin + Multi",

  "fluoroquinolone antibiotic"
    = "Fluoroquinolone",

  "fluoroquinolone;macrolide;penam"
    = "Fluoroquinolone;\nMacrolide; Penam",

  "nalidixic acid"
    = "Nalidixic acid",

  "novobiocin"
    = "Novobiocin",

  "penicillin beta-lactam, cephalosporin, fluoroquinolone antibiotic"
    = "Penicillin \u03b2-lactam;\nCephalosporin; FQ",

  "peptide antibiotic"
    = "Peptide antibiotic",

  "phenicol antibiotic"
    = "Phenicol",

  "phenicol antibiotic, penicillin beta-lactam, cephalosporin, carbapenem, monobactam"
    = "Phenicol; \u03b2-lactam;\nCarbapenem; Monobactam",

  "tetracycline antibiotic, penicillin beta-lactam, cephalosporin, disinfecting agents and antiseptics, phenicol antibiotic, rifamycin antibiotic, glycylcycline, fluoroquinolone antibiotic,carbapenem,monobactam"
    = "TET; \u03b2-lactam; Phenicol;\nDisinfectants; FQ; Carbapenem",

  "tigecycline;rifampin;tetracycline;chloramphenicol;cefalotin;triclosan;ampicillin;acriflavine"
    = "Tigecycline; Rifampin;\nTriclosan; Acriflavine"
)

df_raw$Antibiotic <- antibiotic_map[df_raw$Antibiotic]

# ── 5. RESHAPE TO LONG FORMAT ─────────────────────────────────────────────────
df_long <- df_raw %>%
  pivot_longer(
    cols      = -Antibiotic,
    names_to  = "Gene",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count), Count > 0)

# ── 6. AXIS ORDER ─────────────────────────────────────────────────────────────
# X-axis: antibiotic labels in order they appear in data
xaxis_order <- c(
  "Amikacin; Gentamicin C;\nTobramycin",
  "Aminocoumarin +\nAminoglycoside",
  "Bicyclomycin-like",
  "Cefalothin; Chloramphenicol;\nTetracycline; Rifampin",
  "Ciprofloxacin",
  "Ciprofloxacin + Multi",
  "Fluoroquinolone",
  "Fluoroquinolone;\nMacrolide; Penam",
  "Nalidixic acid",
  "Novobiocin",
  "Penicillin \u03b2-lactam;\nCephalosporin; FQ",
  "Peptide antibiotic",
  "Phenicol",
  "Phenicol; \u03b2-lactam;\nCarbapenem; Monobactam",
  "TET; \u03b2-lactam; Phenicol;\nDisinfectants; FQ; Carbapenem",
  "Tigecycline; Rifampin;\nTriclosan; Acriflavine"
)

# Y-axis: genes alphabetical reversed (soxS top → acrA bottom)
gene_order <- rev(sort(unique(df_long$Gene)))

df_plot <- df_long %>%
  mutate(
    Antibiotic = factor(Antibiotic, levels = xaxis_order),
    Gene       = factor(Gene,       levels = gene_order)
  ) %>%
  filter(!is.na(Antibiotic))

# ── 7. PLOT ───────────────────────────────────────────────────────────────────
p <- ggplot(df_plot, aes(x = Antibiotic, y = Gene)) +

  geom_point(aes(size = Count, colour = Count),
             alpha = 0.92) +
scale_colour_gradientn(
  colours = c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"),
  name    = "No. of\ndetections",
  limits  = c(18, 54),
  guide   = guide_colorbar(
    barwidth = 2, barheight = 14,
    title.position = "top", title.hjust = 0.5,
    ticks.colour   = "grey60"
  )
) +

scale_size_continuous(
  range  = c(8, 18),
  name   = "No. of\ndetections",
  breaks = c(18, 36, 54),
  labels = c("18", "36", "54")
) +

guides(
  colour = guide_colorbar(
    barwidth       = 2,
    barheight      = 14,
    title.position = "top",
    title.hjust    = 0.5
  ),
  size = guide_legend(
    override.aes   = list(colour = "#feb24c", size = c(4, 10, 16)),
    title.position = "top",
    title.hjust    = 0.5
  )
) + 
  labs(
    title    = "AMR Gene vs Antibiotic Class",
    subtitle = "The plot shows the number of strains (indicated by bubble size and colour) that have a given gene-antibiotic association",
    x        = "Antibiotic / Drug Combination",
    y        = "Resistance Gene / Efflux System"
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

    axis.text.x  = element_text(colour = "black", size = 13,
                                 angle = 45, hjust = 1, vjust = 1,
                                 lineheight = 0.85),
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
  filename  = "bubble_gene_vs_antibiotic_WHITE.png",
  plot      = p,
  width     = 28,
  height    = 18,
  dpi       = 600,
  bg        = "white",
  limitsize = FALSE
)

message("Done — saved as bubble_gene_vs_antibiotic_WHITE.png")