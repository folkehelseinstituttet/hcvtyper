#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: summarize_stats_rov.R <prefix> <stats file>", call.=FALSE)
}

prefix <- args[1]
stats  <- args[2]

df <- as.data.frame(matrix(nrow = 1, ncol = 3))
colnames(df) <- c("sampleName", "contig", "mapped_reads")

# Read the stats file, extract the name of the sequence that was mapped to and the number of mapped reads
mapped_seq <- read_tsv(stats, col_names = FALSE, n_max = 7) %>%
    # Extract the command used
    filter(str_detect(X1, "The command line was")) %>%
    # Separate out the pattern after the string "reference"
    separate(X1, into = c(NA, "tmp"), sep = "reference ") %>%
    # Separate out the pattern before the string ".fasta"
    separate(tmp, into = c("seq", NA), sep = "\\.fasta")

mapped_reads <- read_tsv(stats, col_names = FALSE, comment = "#") %>%
    filter(str_detect(X1, "reads mapped:")) %>%
    separate(X1, into = c(NA, "reads"), sep = "reads mapped:   ") %>%
    pull(reads)

df$sampleName <- str_remove(prefix, "\\.withdup")
df$contig <- mapped_seq
df$mapped_reads <- mapped_reads

write_csv(paste0(prefix, ".", mapped_seq, ".csv"), df)
