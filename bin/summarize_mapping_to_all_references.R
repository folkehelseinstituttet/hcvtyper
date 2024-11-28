#!/usr/bin/env Rscript

library(tidyverse)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: summarize_mapping_to_all_references.R <idxstats file> <depth file> <sample name> <references>", call.=FALSE)
}

idxstats   <- args[1]
depth      <- args[2]
sampleName <- args[3]
references <- args[4]

# First calculate coverage for all references
# Read the depth file from the first mapping.
# The file can be empty and the reading fails
cov <- read_tsv(depth, col_names = FALSE) %>% 
  group_by(X1) %>% # Group by name of the reference
  summarise(
    total_rows = n(), # Get the total number of positions for the reference (genome length)
    count_gt_4 = sum(X3 > 4), # Get the number of positions with coverage >= 5
    percent_gt_4 = (count_gt_4 / total_rows) * 100,
    percent_gt_4_int = round(percent_gt_4, digits = 0) # Round to nearest integer. Need integer for groovy/nextflow filtering later
  )

# Then read mapped reads from the first mapping
# Read the idsxtats output of the first mapping
df <- read_table(idxstats, col_names = FALSE) %>%
  # Separate subtype and reference from the sequence names
  separate(X1, into = c("Subtype", "Reference"), sep = "_", remove = FALSE) %>%
  # Discard the unmapped reads marked by an * (more precisely these are unmapped reads without coordinates)
  filter(X1 != "*") %>%
  # Separate the genotype from the subtype.
  # For 2k1b we use the whole name for genotype also
  mutate(Genotype = if_else(Subtype == "2k1b", Subtype, substr(Subtype, 1, 1))) #%>%
  # Rename Genotype 2k1b to 1 as a preparation for detecting minor genotypes
  #mutate(Genotype = str_replace(Genotype, "2k1b", "1"))

# Join the percent coverage to the mapping statistics
df <- left_join(df, cov, by = c("X1" = "X1"))

# Create empty final dataframe to populate
df_final <- as.data.frame(matrix(nrow = 1, ncol = 8))
colnames(df_final) <- c("sample", "total_mapped_reads", "major_ref", "major_reads", "major_cov", "minor_ref", "minor_reads", "minor_cov")

# Add sample name
df_final$sample[1] <- sampleName

# Sometimes the mappings stats are completely empty
if (nrow(df) > 0) {

# First get the total number of mapped reads to all references
df_final$total_mapped_reads[1] <- sum(df$X3, na.rm = TRUE)

# Count number of reads per subtype
summary <- df %>%
  # Group by Genotype also to retain that column
  group_by(Subtype, Genotype) %>%
  summarise(reads = sum(X3)) %>%
  arrange(desc(reads))

## Major
# Find major reference to use
major_tmp <- summary$Subtype[1]
major_ref <- df %>%
  filter(Subtype == major_tmp) %>%
  # Choose the reference with most mapped reads
  arrange(desc(X3)) %>%
  head(n = 1) %>%
  pull(X1)

df_final$major_ref[1] <- major_ref

# How many read pairs mapped to the major subtype
major_reads <- summary %>%
  filter(Subtype == major_tmp) %>%
  pull(reads)

df_final$major_reads[1] <- major_reads

# Add the major ref coverage
df_final$major_cov[1] <- df %>% filter(X1 == major_ref) %>% pull(percent_gt_4_int)

## Minor
# Only execute if two or more references have reads mapped
# And require that the reference with second most reads belong to a different genotype than the major
major_genotype <- head(summary$Genotype, n = 1)

if (major_genotype == "1" | major_genotype == "2k1b") { # If major genotype is 1 or 2k1b, discard gt 1 and 2k1b as potential minor
  tmp <- df %>%
    # Remove major Genotype
    filter(Genotype != major_genotype) %>%
    filter(Genotype != "2k1b") %>%
    # Choose the reference with the higest coverage
    arrange(desc(percent_gt_4)) %>%
    head(n = 1)
} else { # If the major gt is not 1 or 2k1b, then allow for these genotypes to be present among potential minor
  tmp <- df %>%
    # Remove major Genotype
    filter(Genotype != major_genotype) %>%
    # Choose the reference with the highest coverage
    arrange(desc(percent_gt_4)) %>%
    head(n = 1)
}
minor_tmp <- tmp %>% pull(Subtype)
minor_ref <- tmp %>% pull(X1)

if (length(minor_ref > 0)) {
  df_final$minor_ref[1] <- minor_ref
}

# How many reads mapped to the minor subtype
minor_reads <- summary %>%
  filter(Subtype == minor_tmp) %>%
  pull(reads)

df_final$minor_reads[1] <- minor_reads

# Add minor coverage
df_final$minor_cov[1] <- df %>% filter(X1 == minor_ref) %>% pull(percent_gt_4_int)
}

# Write results
write_csv(df_final, file = paste0(sampleName, ".parsefirstmapping.csv"))

# Write out the fasta sequences
# Read the reference fasta file
fasta <- read.fasta(file = references)
write.fasta(sequences = fasta[major_ref], names = major_ref, file.out = paste0(sampleName, ".", major_ref, "_major.fa"))
if (length(minor_ref > 0)) {
write.fasta(sequences = fasta[minor_ref], names = minor_ref, file.out = paste0(sampleName, ".", minor_ref, "_minor.fa"))
}
