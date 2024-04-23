#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: bam_coverage.R <depth file>", call.=FALSE)
}

depth <- args[1]
sampleName <- unlist(str_split(basename(depth), pattern = "\\."))[1]
major_minor <- unlist(str_split(basename(depth), pattern = "\\."))[3]

df <- read_tsv(depth, col_names = FALSE)

# If there are no mapped reads df is empty
if (nrow(df) > 0) {
  df <- df %>% 
    # Rename columns
    rename("Reference" = X1,
           "Position" = X2,
           "Coverage" = X3)

# Get the mapped reference
reference <- df %>% distinct(Reference) %>% pull(Reference)

plot <- df %>% 
    ggplot() +
    aes(x = Position, y = Coverage) + 
    geom_line() +
    geom_hline(yintercept = 10, color = "darkgreen", linetype = "dotted") +
    annotate("text", x=900, y=-150, label="Coverage cutoff = 10") +
    ggtitle(paste0(sampleName, ".", major_minor, ".", reference))

# Save the plot
ggsave(plot = plot, 
       file = paste0(sampleName, ".", major_minor, ".", reference, ".png"), 
         device = "png", 
         dpi = 300)
}