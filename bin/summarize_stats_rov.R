#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: summarize_stats_rov.R <prefix> <stats file>", call.=FALSE)
}

prefix <- args[1]
stats  <- args[2]

df <- as.data.frame(matrix(nrow = 1, ncol = 4))
colnames(df) <- c("sampleName", "gene", "contig", "mapped_reads")

# Read the stats file, extract the name of the sequence that was mapped to and the number of mapped reads
tmp <- read_tsv(stats, col_names = FALSE, n_max = 7) %>%
    # Extract the command used
    filter(str_detect(X1, "The command line was")) %>%
    # Separate out the pattern after the string "reference"
    separate(X1, into = c(NA, "tmp"), sep = "reference ") %>%
    # Separate out info from the fasta file name
    separate(tmp, into = c("sample", "gene,", "seq", NA), sep = "\\.")

mapped_seq <- tmp %>% pull(seq)
prefix_2 <- tmp %>% pull(sample)
gene <- tmp %>% pull(`gene,`)

# Compare the input prefix with prefix_2. Error if they are not identical
# Allow for the possibility of the input prefix having a suffix of ".withdup" or ".markdup"
if (str_remove(prefix, "\\.withdup|\\.markdup") != prefix_2) {
    stop("Prefixes do not match: ", prefix, " vs ", prefix_2, call.=FALSE)
}

mapped_reads <- read_table(stats, col_names = FALSE, comment = "#", n_max = 10) %>%
    filter(X3 == "mapped:") %>%
    #separate(X1, into = c(NA, "reads"), sep = "reads mapped:   ") %>%
    pull(X4)

df$sampleName <- prefix_2
df$contig <- mapped_seq
df$gene <- gene
df$mapped_reads <- mapped_reads

df <- tibble(df)

write_csv(df, paste0(prefix, ".", mapped_seq, ".csv"))
