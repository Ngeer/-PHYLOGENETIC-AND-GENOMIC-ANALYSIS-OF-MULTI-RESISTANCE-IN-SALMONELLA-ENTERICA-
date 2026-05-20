library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# ── 1. Load ───────────────────────────────────────────────────────
raw <- read.csv("C:\\Users\\Stuart\\Downloads\\snp_distance_matrix.csv",
                row.names = 1, check.names = FALSE)

rownames(raw) <- gsub("^snippy_", "", rownames(raw))
colnames(raw) <- gsub("^snippy_", "", colnames(raw))

mat <- as.matrix(raw)
storage.mode(mat) <- "numeric"

# ── 2. Remove outlier ─────────────────────────────────────────────
keep <- setdiff(rownames(mat), "YA00509485")
mat  <- mat[keep, keep]

# ── 3. Cluster ────────────────────────────────────────────────────
ord       <- hclust(as.dist(mat), method = "ward.D2")$order
iso_names <- rownames(mat)[ord]
mat_ord   <- mat[iso_names, iso_names]

# ── 4. Long format ────────────────────────────────────────────────
df <- expand.grid(Row = iso_names, Col = iso_names,
                  stringsAsFactors = FALSE)
df$SNPs    <- mapply(function(r, c) mat_ord[r, c], df$Row, df$Col)
df$is_diag <- df$Row == df$Col

# KEY FIX: diagonal = 0, show as "0"; all others show number
df$Label <- ifelse(df$is_diag, "0", formatC(df$SNPs, format="d", big.mark=","))

# KEY FIX: diagonal gets WHITE fill so it stands out as reference
df$fill_val <- ifelse(df$is_diag, NA, df$SNPs)

df$Row <- factor(df$Row, levels = rev(iso_names))
df$Col <- factor(df$Col, levels = iso_names)

# Text colour: white on dark cells, black on light cells
max_snp    <- max(df$SNPs[!df$is_diag])
df$txt_col <- ifelse(df$SNPs > max_snp * 0.45, "white", "black")

# ── 5. Plot ───────────────────────────────────────────────────────
p <- ggplot(df, aes(x = Col, y = Row)) +

  # Background tile using fill_val (NA diagonal = grey)
  geom_tile(aes(fill = fill_val), colour = "white", linewidth = 0.5) +

  # Numbers on every single cell including diagonal
  geom_text(aes(label = Label, colour = txt_col),
            size = 2.6, fontface = "bold") +

  scale_colour_identity() +

  scale_fill_gradient(
    low      = "#FFFFFF",   # white = low SNPs
    high     = "#000000",   # black = high SNPs
    na.value = "#E8E8E8",   # light grey = diagonal (self vs self)
    name     = "SNP Distance",
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(10, "cm"),
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
    subtitle = "Whole-genome SNP distances · Ward's D2 clustering · YA00509485 excluded (outlier ~298,000 SNPs)",
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

# ── 6. Save ───────────────────────────────────────────────────────
ggsave("C:\\Users\\Stuart\\Downloads\\snp_heatmap.png",
       plot = p, width = 15, height = 14, dpi = 300, bg = "white")

ggsave("C:\\Users\\Stuart\\Downloads\\snp_heatmap.pdf",
       plot = p, width = 15, height = 14, bg = "white")

cat("Done!\n")