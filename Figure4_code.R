############################################################
## Figure 4. RNA-seq reproducibility analysis
## (A) Correlation density plot (Same vs Different)
## (B) Correlation heatmap
## Output: TIFF (600???800 dpi)
## Reproducible & journal-ready
############################################################

# ==========================================================
# 0. Packages
# ==========================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

# ==========================================================
# 1. Input files
# ==========================================================
## Figure 4A
tpm_file   <- "combine.All.2024.ECnt.txt"
sinfo_file <- "sinfo_all.txt"

## Figure 4B
expr_file   <- "2024_data_expression_vehicle.txt"
header_file <- "2024_data_expression_vehicle_header.txt"

stopifnot(
  file.exists(tpm_file),
  file.exists(sinfo_file),
  file.exists(expr_file),
  file.exists(header_file)
)

# ==========================================================
# 2. Figure 4A ??? Correlation density plot
# ==========================================================

## ---- 2-1. Load data ----
tpm   <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1)
sinfo <- read.table(sinfo_file, sep = "\t", header = TRUE)

## ---- 2-2. Define experimental group ----
sinfo$Group <- apply(
  sinfo[, c("Treatment_1", "Dos", "Solvent", "Cell")],
  1, paste, collapse = "_"
)

## ---- 2-3. Correlation matrix ----
cm <- cor(tpm, use = "pairwise.complete.obs")
diag(cm) <- NA

## ---- 2-4. Same vs Different correlations per cell ----
cordf_list <- list()

for (cell in unique(sinfo$Cell)) {
  
  cell_idx <- sinfo$Cell == cell
  groups   <- unique(sinfo$Group[cell_idx])
  
  same_cor <- c()
  diff_cor <- c()
  
  for (g in groups) {
    idx <- sinfo$Group == g
    
    ## Same condition
    if (sum(idx) > 1) {
      same_cor <- c(
        same_cor,
        cm[idx, idx][upper.tri(cm[idx, idx])]
      )
    }
    
    ## Different condition
    diff_cor <- c(
      diff_cor,
      as.vector(cm[idx, !idx])
    )
  }
  
  cordf_list[[cell]] <- rbind(
    data.frame(Cell = cell, Condition = "Same",      Correlation = same_cor),
    data.frame(Cell = cell, Condition = "Different", Correlation = diff_cor)
  )
}

cordf_all <- bind_rows(cordf_list)

## ---- 2-5. Mean correlation (for dashed line) ----
mu <- cordf_all %>%
  group_by(Condition) %>%
  summarise(grp.mean = mean(Correlation, na.rm = TRUE))

## ---- 2-6. Theme ----
theme_set(
  theme_classic() +
    theme(
      legend.position = "top",
      axis.title = element_text(size = 13),
      axis.text  = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    )
)

## ---- 2-7. Save Figure 4A ----
tiff(
  filename = "Figure4A_RNAseq_correlation_density.tiff",
  width = 8,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

ggplot(cordf_all, aes(x = Correlation)) +
  geom_density(aes(fill = Condition), alpha = 0.4) +
  geom_vline(
    data = mu,
    aes(xintercept = grp.mean, color = Condition),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c(Different = "#1f78b4", Same = "#EFC000FF")) +
  scale_color_manual(values = c(Different = "#1f78b4", Same = "#EFC000FF")) +
  labs(
    x = "Pearson correlation coefficient",
    y = "Density",
    fill  = "Condition",
    color = "Condition"
  )

dev.off()

# ==========================================================
# 3. Figure 4B ??? Correlation heatmap
# ==========================================================

## ---- 3-1. Load data ----
expr_data   <- read.table(expr_file, sep = "\t", header = TRUE, row.names = 1)
header_data <- read.table(header_file, sep = "\t", header = TRUE, row.names = 1)

## ---- 3-2. Correlation matrix ----
corr_matrix <- cor(expr_data, use = "pairwise.complete.obs")

## ---- 3-3. Save Figure 4B heatmap ----
tiff(
  filename = "Figure4B_correlation_heatmap.tiff",
  width = 5,
  height = 5,
  units = "in",
  res = 800,
  compression = "lzw"
)

pheatmap(
  corr_matrix,
  color = colorRampPalette(c("lightgray", "gray", "darkred"))(50),
  breaks = seq(0.75, 1, length.out = 50),
  border_color = "black",
  angle_col = 315,
  cellheight = 5,
  cellwidth  = 5,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row  = 5,
  legend = FALSE
)

dev.off()

# ==========================================================
# 4. (Optional) Figure 4B annotation heatmap
# ==========================================================
annotation <- header_data[, 3:4]

tiff(
  filename = "Figure4B_annotation_heatmap.tiff",
  width = 2.5,
  height = 5,
  units = "in",
  res = 800,
  compression = "lzw"
)

pheatmap(
  annotation,
  color = colorRampPalette(
    c("red", "orange", "darkgreen", "blue",
      "pink", "cyan", "green", "hotpink")
  )(70),
  border_color = "gray",
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  legend = FALSE
)

dev.off()

message("Figure 4A and 4B successfully generated.")