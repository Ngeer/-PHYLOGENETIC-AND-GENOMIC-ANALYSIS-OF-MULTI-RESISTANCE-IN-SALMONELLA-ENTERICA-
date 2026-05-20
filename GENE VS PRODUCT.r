# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
install.packages(c("dplyr", "tidyr", "ggplot2", "stringr"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)




pkgs <- c("dplyr", "tidyr", "ggplot2", "stringr")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

setwd("C:\\Users\\Stuart\\Downloads")   # ← adjust path

# ── 1. LOAD RAW FILE ──────────────────────────────────────────────────────────
df_raw <- read.csv("GENE VS PRODUCT.csv",
                   fileEncoding = "UTF-16",
                   check.names  = FALSE,
                   sep          = "\t",
                   header       = FALSE)

# ── 2. EXTRACT GENE NAMES & BUILD DATA FRAME ──────────────────────────────────
gene_names          <- as.character(df_raw[1, 2:ncol(df_raw)])
data_rows           <- df_raw[2:nrow(df_raw), ]
colnames(data_rows) <- c("Product", gene_names)
rownames(data_rows) <- NULL

# ── 3. CLEAN PRODUCT COLUMN ───────────────────────────────────────────────────
data_rows$Product <- trimws(as.character(data_rows$Product))
data_rows         <- data_rows[data_rows$Product != "" & data_rows$Product != "0", ]

# Fix Excel formula corruption: "efflu+C8x" → "efflux"
data_rows$Product <- gsub("efflu\\+C8x", "efflux", data_rows$Product)
data_rows$Product <- trimws(data_rows$Product)

# ── 4. DUPLICATE CHECK ────────────────────────────────────────────────────────
cat("=== DUPLICATE CHECK (before merge) ===\n")
dup_tbl <- sort(table(data_rows$Product))
cat("Products with duplicates:\n")
print(names(dup_tbl[dup_tbl > 1]))

# ── 5. MERGE DUPLICATES BY SUMMING ────────────────────────────────────────────
data_rows[gene_names] <- lapply(data_rows[gene_names], function(x) as.numeric(x))

df_merged <- data_rows %>%
  group_by(Product) %>%
  summarise(across(all_of(gene_names),
                   ~ sum(., na.rm = TRUE)),
            .groups = "drop")

cat("\n=== AFTER MERGE ===\n")
cat("Unique products:", nrow(df_merged), "\n")

# ── 6. RESHAPE TO LONG FORMAT ─────────────────────────────────────────────────
df_long <- df_merged %>%
  pivot_longer(
    cols      = all_of(gene_names),
    names_to  = "Gene",
    values_to = "Count"
  ) %>%
  filter(!is.na(Count), Count > 0)

cat("\nUnique Count values:", sort(unique(df_long$Count)), "\n")

# ── 7. SHORT LABEL LOOKUP TABLE ───────────────────────────────────────────────
label_map <- c(
  "Aminoglycosides efflux system AcrAD-TolC, inner-membrane proton/drug antiporter AcrD (RND type)" = "AcrAD-TolC antiporter AcrD (aminoglycosides)",
  "Copper-sensing two-component system response regulator CpxR"                                      = "Two-component regulator CpxR",
  "CRP is a global regulator that represses MdtEF multidrug efflux pump expression."                 = "CRP - MdtEF efflux repressor",
  "DNA-binding transcriptional dual regulator Rob"                                                   = "Transcriptional regulator Rob",
  "DNA-binding transcriptional dual regulator SoxS"                                                  = "Transcriptional regulator SoxS",
  "EmrD is a multidrug transporter from the Major Facilitator Superfamily (MFS"                      = "MFS transporter EmrD",
  "Lipopolysaccharide core heptose(II)-phosphate phosphatase @ Polymyxin resistance protein PmrG"    = "Polymyxin resistance protein PmrG",
  "Multidrug efflux system AcrEF-TolC, membrane fusion component AcrE"                              = "AcrEF-TolC fusion component AcrE",
  "Multidrug efflux system AcrAB-TolC, inner-membrane proton/drug antiporter AcrB (RND type)"       = "AcrAB-TolC antiporter AcrB (RND)",
  "Multidrug efflux system AcrAB-TolC, membrane fusion component AcrA"                              = "AcrAB-TolC fusion component AcrA",
  "Multidrug efflux system AcrEF-TolC, inner-membrane proton/drug antiporter AcrF (RND type)"       = "AcrEF-TolC antiporter AcrF (RND)",
  "Multidrug efflux system EmrAB-OMF, inner-membrane proton/drug antiporter EmrB (MFS type)"        = "EmrAB-OMF antiporter EmrB (MFS)",
  "Multidrug efflux system EmrAB-OMF, membrane fusion component EmrA"                               = "EmrAB-OMF fusion component EmrA",
  "Multidrug efflux system MdtABC-TolC, inner-membrane proton/drug antiporter MdtB (RND type)"      = "MdtABC-TolC antiporter MdtB (RND)",
  "Multidrug efflux system MdtABC-TolC, inner-membrane proton/drug antiporter MdtC (RND type)"      = "MdtABC-TolC antiporter MdtC (RND)",
  "Multidrug efflux system MdtABC-TolC, membrane fusion component MdtA"                             = "MdtABC-TolC fusion component MdtA",
  "Multidrug efflux system, inner membrane proton/drug antiporter (RND type) => MexQ of MexPQ-OpmE system" = "MexPQ-OpmE antiporter MexQ (RND)",
  "Multidrug efflux system, membrane fusion component => MexP of MexPQ-OpmE system"                 = "MexPQ-OpmE fusion component MexP",
  "Multidrug efflux system, outer membrane factor lipoprotein of OprM/OprM family"                  = "OprM outer membrane lipoprotein",
  "Multidrug efflux transporter MdtK/NorM (MATE family)"                                            = "MATE transporter MdtK/NorM",
  "Multidrug resistance regulator EmrR (MprA)"                                                      = "Multidrug resistance regulator EmrR",
  "Multidrug resistance transporter => Bicyclomycin resistance protein Bcr-1"                       = "Bicyclomycin resistance protein Bcr-1",
  "Multiple antibiotic resistance protein MarA"                                                      = "MAR protein MarA",
  "Multiple antibiotic resistance protein MarR"                                                      = "MAR protein MarR",
  "RamA (resistance antibiotic multiple) is a positive regulator of AcrAB-TolC"                     = "AcrAB-TolC positive regulator RamA",
  "Response regulator BaeR"                                                                          = "Response regulator BaeR",
  "Sensory histidine kinase BaeS"                                                                    = "Sensory histidine kinase BaeS",
  "Transcriptional regulator of acrAB operon, AcrR"                                                 = "acrAB operon regulator AcrR"
)

# Apply short labels
df_long$Product <- ifelse(df_long$Product %in% names(label_map),
                          label_map[df_long$Product],
                          df_long$Product)

# ── 8. AXIS ORDER ─────────────────────────────────────────────────────────────
product_order <- sort(unique(df_long$Product))
gene_order    <- rev(sort(unique(df_long$Gene)))

df_plot <- df_long %>%
  mutate(
    Product = factor(Product, levels = product_order),
    Gene    = factor(Gene,    levels = gene_order)
  )

# ── 9. BUBBLE PLOT ────────────────────────────────────────────────────────────
p <- ggplot(df_plot, aes(x = Gene, y = Product)) +

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
    range  = c(4, 14),
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
    title    = "AMR Gene against Product",
    subtitle = "The plot shows the number of strains (indicated by bubble size and colour) that have a given gene-product association",
    x        = "Resistance Gene Efflux System",
    y        = "Product / Functional Description"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = "grey80",
                                     linewidth = 0.3),
    legend.key        = element_rect(fill = "white", colour = NA),

    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_blank(),

    axis.text.x  = element_text(colour = "black", size = 9,
                                 angle = 45, hjust = 1, vjust = 1,
                                 face = "italic"),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold",
                                 margin = margin(t = 14)),

    axis.text.y  = element_text(colour = "black", size = 10),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold",
                                 angle = 90, margin = margin(r = 12)),

    axis.ticks = element_line(colour = "grey70", linewidth = 0.3),

    plot.title    = element_text(colour = "black", size = 18, face = "bold",
                                  hjust = 0.5, margin = margin(b = 2)),
    plot.subtitle = element_text(colour = "grey30", size = 16,
                                  hjust = 0.5, margin = margin(b = 14)),
    plot.margin   = margin(20, 20, 20, 20),

    legend.text  = element_text(colour = "black", size = 11),
    legend.title = element_text(colour = "black", size = 11.5, face = "bold"),
    legend.position = "right",
    legend.box   = "vertical"
  )

# ── 10. SAVE ──────────────────────────────────────────────────────────────────
ggsave(
  filename  = "bubble_gene_vs_product.png",
  plot      = p,
  width     = 16,
  dpi = 600,
  height    = 12,
  bg        = "white",
  limitsize = FALSE
)

message("Done — saved as bubble_gene_vs_product.png")