
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!require("viridis", quietly = TRUE)) install.packages("viridis")

library(ggplot2)
library(reshape2)
library(viridis)

setwd("C:\\Users\\Stuart\\Downloads")

ani_raw <- read.csv("ANB.csv", row.names = 1, check.names = FALSE)

ani_mat <- as.matrix(ani_raw)
storage.mode(ani_mat) <- "numeric"

ani_mat <- ifelse(ani_mat >= 9700, ani_mat / 100, ani_mat)
diag(ani_mat) <- 100

diag(ani_mat)[is.na(diag(ani_mat))] <- 100
for (i in 1:nrow(ani_mat))
  ani_mat[i, is.na(ani_mat[i,])] <- mean(ani_mat[i,], na.rm = TRUE)

cat("NAs left:", sum(is.na(ani_mat)), "\n")

row_clust <- hclust(dist(ani_mat), method = "complete")
col_clust <- hclust(dist(t(ani_mat)), method = "complete")

ani_mat_clust <- ani_mat[row_clust$order, col_clust$order]

ani_long <- melt(ani_mat_clust)
colnames(ani_long) <- c("Row", "Col", "Value")

p <- ggplot(ani_long, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Value), 
                color = ifelse(Value >= 99, "white", "black")), 
            size = 2.5, fontface = "bold") +
  scale_fill_gradient(low = "#e0f7da", high = "#003d1a", 
                      limits = c(97, 100),
                      name = "ANIb (%)", guide = guide_colorbar(title.position = "top")) +
  scale_color_identity() +
  labs(title = "Genomic Relatedness of Salmonella Isolates",
       subtitle = "Hierarchical Clustering of Salmonella Isolates Based on Average Nucleotide Identity (ANIb) Values",
       x = "Isolates", y = "Isolates") +
  theme_minimal() +
  theme(
    aspect.ratio = 0.6,
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, 
                              margin = margin(t = 15, b = 5), color = "#1a1a1a"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, 
                                 margin = margin(b = 15), color = "#555555"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, color = "#333333"),
    axis.text.y = element_text(size = 10, color = "#333333"),
    axis.title = element_text(size = 12, face = "bold", color = "#333333"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#f8f8f8", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

ggsave("anib_heatmap.png", plot = p, width = 20, height = 15, dpi = 600, units = "in", bg = "white")

cat("Done! File saved in:", getwd(), "\n")