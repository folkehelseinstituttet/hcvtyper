#!/usr/bin/env Rscript

library(tidyverse)

# List all files from PARSE_PHYLOGENY
parse_phylo_files <- list.files(pattern = "parse_phylo",
                                full.names = TRUE)

# List all files from CALCULATE_PAIRWISE_ALIGNMENT_METRICS
alignment_metrics_files <- list.files(pattern = "alignment_metrics",
                                full.names = TRUE)

