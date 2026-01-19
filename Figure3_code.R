############################################################
## Mapping reads stacked bar plots (4 cell lines)
## Unique vs Multiple mapped reads per sample
##
## Input files:
##   - HepG2.mappingReads.All.txt
##   - SW1783.mappingReads.All.txt
##   - A549.mappingReads.All.txt
##   - HT29.mappingReads.All.txt
##
## Output:
##   - *_mapping_reads_stacked_bar.tiff
##
## Description:
##   - Samples are ordered by total mapped reads (descending)
##   - X-axis sample labels are hidden
##   - Axis font sizes increased (+1)
##   - Legend font sizes increased (+1)
##
## Author: ---
## Date  : ---
############################################################

## ---- 1. Load libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)

## ---- 2. Define input files ----
input_files <- c(
  "HepG2.mappingReads.All.txt",
  "SW1783.mappingReads.All.txt",
  "A549.mappingReads.All.txt",
  "HT29.mappingReads.All.txt"
)

## ---- 3. Loop over cell line files ----
for (file in input_files) {
  
  ## ---- 3-1. Extract cell line name ----
  cell_line <- sub("\\.mappingReads\\.All\\.txt$", "", file)
  message("Processing: ", cell_line)
  
  ## ---- 3-2. Load data ----
  mapping_df <- read.table(
    file      = file,
    sep       = "\t",
    header    = TRUE,
    stringsAsFactors = FALSE
  )
  
  ## Expected columns:
  ## ID | Unique_Reads | Multiple_Reads
  
  ## ---- 3-3. Reshape data (wide → long) ----
  mapping_long <- mapping_df %>%
    pivot_longer(
      cols = c(Unique_Reads, Multiple_Reads),
      names_to  = "ReadType",
      values_to = "Reads"
    )
  
  ## ---- 3-4. Order samples by total mapped reads ----
  mapping_long <- mapping_long %>%
    mutate(
      ID = fct_reorder(ID, Reads, .fun = sum, .desc = TRUE)
    )
  
  ## ---- 3-5. Create stacked bar plot ----
  p <- ggplot(mapping_long, aes(x = ID, y = Reads, fill = ReadType)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(
      values = c(
        Unique_Reads   = "#1f78b4",  # blue
        Multiple_Reads = "grey70"    # grey
      )
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      x = "Sample",
      y = "Number of mapped reads",
      fill = "Read type"
    ) +
    theme_classic() +
    theme(
      ## X-axis
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      
      ## Axis font sizes (+1)
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text.y  = element_text(size = 13),
      
      ## Legend font sizes (+1)
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    )
  
  ## ---- 3-6. Save figure (publication-ready) ----
  output_fig <- paste0(cell_line, "_mapping_reads_stacked_bar.tiff")
  
  ggsave(
    filename = output_fig,
    plot     = p,
    width    = 8,
    height   = 5,
    dpi      = 600
  )
}

############################################################
## End of script
############################################################