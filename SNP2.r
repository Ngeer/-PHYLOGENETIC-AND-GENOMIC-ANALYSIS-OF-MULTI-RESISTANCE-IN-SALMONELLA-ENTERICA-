library(ggplot2)
library(dplyr)
library(scales)

# ── 1. Load ───────────────────────────────────────────────────────
raw <- read.csv("C:\\Users\\Stuart\\Downloads\\snp_distance_matrix.csv",
                row.names = 1, check.names = FALSE)

rownames(raw) <- gsub("^snippy_", "", rownames(raw))
colnames(raw) <- gsub("^snippy_", "", colnames(raw))

mat <- as.matrix(raw)
storage.mode(mat) <- "numeric"

# ── 2. Cluster ────────────────────────────────────────────────────
ord       <- hclust(as.dist(mat), method = "ward.D2")$order
iso_names <- rownames(mat)[ord]
mat_ord   <- mat[iso_names, iso_names]

# ── 3. Long format ────────────────────────────────────────────────
df <- expand.grid(Row = iso_names, Col = iso_names,
                  stringsAsFactors = FALSE)
df$SNPs    <- mapply(function(r, c) mat_ord[r, c], df$Row, df$Col)
df$is_diag <- df$Row == df$Col
df$is_YA   <- df$Row == "YA00509485" | df$Col == "YA00509485"

# Cell type: "diag", "ya", "normal"
df$cell_type <- ifelse(df$is_diag, "diag",
                ifelse(df$is_YA,   "ya", "normal"))

# Fill value for normal cells only (capped at 60k)
df$fill_normal <- ifelse(df$cell_type == "normal",
                         pmin(df$SNPs, 60000), NA)

# Label on every cell
df$Label <- formatC(df$SNPs, format = "d", big.mark = ",")

# Text colour
df$txt_col <- ifelse(df$cell_type == "diag",   "black",
              ifelse(df$cell_type == "ya",      "white",
              ifelse(df$SNPs > 27000,           "white", "black")))

df$Row <- factor(df$Row, levels = rev(iso_names))
df$Col <- factor(df$Col, levels = iso_names)

# Split into separate data frames — fixes the length mismatch error
df_normal <- df[df$cell_type == "normal", ]
df_ya     <- df[df$cell_type == "ya",     ]
df_diag   <- df[df$cell_type == "diag",   ]

# ── 4. Plot ───────────────────────────────────────────────────────
p <- ggplot(df, aes(x = Col, y = Row)) +

  # Normal cells: white to black
  geom_tile(data = df_normal,
            aes(fill = fill_normal),
            colour = "white", linewidth = 0.4) +

  # YA cells: solid dark red
  geom_tile(data = df_ya,
            fill = "#8B0000",
            colour = "white", linewidth = 0.4) +

  # Diagonal cells: solid grey
  geom_tile(data = df_diag,
            fill = "#BBBBBB",
            colour = "white", linewidth = 0.4) +

  # Labels on every cell
  geom_text(aes(label = Label, colour = txt_col),
            size = 2.4, fontface = "bold") +

  scale_colour_identity() +

  scale_fill_gradient(
    low      = "#FFFFFF",
    high     = "#000000",
    limits   = c(0, 60000),
    na.value = "transparent",
    name     = "SNP Distance",
    labels   = label_comma(),
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(9, "cm"),
      barheight      = unit(0.5, "cm"),
      frame.colour   = "black",
      ticks.colour   = "black"
    )
  ) +

  scale_x_discrete(expand = c(0, 0), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +

  labs(
    title    = expression(paste("Pairwise SNP Distance Matrix of ",
                                italic("Salmonella"), " Isolates")),
    subtitle = "White-Black scale: 0-60,000 SNPs  |  Dark red = YA00509485 (~295,000-299,000 SNPs)  |  Ward's D2 clustering",
    x = NULL, y = NULL
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(size=15, face="bold",
                                    hjust=0.5, margin=margin(b=4)),
    plot.subtitle    = element_text(size=9, hjust=0.5,
                                    colour="grey40", margin=margin(b=10)),
    plot.margin      = margin(15, 15, 15, 15),
    plot.background  = element_rect(fill="white", colour=NA),
    panel.background = element_rect(fill="white", colour=NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle=45, hjust=1, vjust=1,
                                    size=10, face="bold", colour="black"),
    axis.text.y      = element_text(size=10, face="bold", colour="black"),
    axis.ticks       = element_blank(),
    legend.position  = "top",
    legend.title     = element_text(size=10, face="bold"),
    legend.text      = element_text(size=9)
  )

# ── 5. Save ───────────────────────────────────────────────────────
ggsave("C:\\Users\\Stuart\\Downloads\\snp_heatmap.png",
       plot = p, width = 16, height = 15, dpi = 300, bg = "white")

ggsave("C:\\Users\\Stuart\\Downloads\\snp_heatmap.pdf",
       plot = p, width = 16, height = 15, bg = "white")

cat("Done!\n")