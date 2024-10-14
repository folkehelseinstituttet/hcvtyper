#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: summarize_stats_rov.R <prefix> <stats file>", call.=FALSE)
}

prefix <- args[1]
gene   <- args[2]
depth  <- args[3]

# Get the sample name from the depth file. Compare with the prefix
prefix_2 <- unlist(str_split(basename(depth), pattern = "\\."))[1]
if (prefix != prefix_2) {
    stop("Prefixes do not match: ", prefix, " vs ", prefix_2, call.=FALSE)
}

# Empty df
df <- as.data.frame(matrix(nrow = 1, ncol = 6))
colnames(df) <- c("sampleName", "reference", "cov_breadth_min_1", "cov_breadth_min_5", "cov_breadth_min_10", "avg_depth")

# Add the sample name
df$sampleName <- prefix

# Read the depth per position
cov <- read_tsv(depth, col_names = FALSE)

# Reference length
ref_length <- nrow(cov)

# Reference name
reference <- cov %>% distinct(X1) %>% pull(X1)
df$reference <- reference

# Average depth
df$avg_depth <- mean(cov$X3)

# Nr. of positions with coverage >=1, >= 5 and > 9
# If ref_length is zero it means that no reads were mapped. Set coverage to zero.
# Coverage may also be zero if there are reads mapped, but never more than 5 per position
if (ref_length > 0) {
    pos_1 <- nrow(
        cov %>%
            filter(X3 >= 1)
    )
    pos_5 <- nrow(
        cov %>%
            filter(X3 >= 5)
    )
    pos_10 <- nrow(
        cov %>%
            filter(X3 > 9)
    )

    # Coverage breadth
    breadth_1 <- round(pos_1 / ref_length * 100, digits = 2)
    df$cov_breadth_min_1 <- breadth_1

    breadth_5 <- round(pos_5 / ref_length * 100, digits = 2)
    df$cov_breadth_min_5 <- breadth_5

    breadth_10 <- round(pos_10 / ref_length * 100, digits = 2)
    df$cov_breadth_min_10 <- breadth_10
} else if (ref_length == 0) {
    df$cov_breadth_min_1 <- 0
    df$cov_breadth_min_5 <- 0
    df$cov_breadth_min_10 <- 0
}

df <- tibble(df)

write_csv(df, paste0(prefix, ".", gene, ".", reference, ".markdup", ".csv"))
