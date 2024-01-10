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

# Create empty dataframe to populate
df_final <- as.data.frame(matrix(nrow = 1, ncol = 7))
colnames(df_final) <- c("sample", "major_ref", "major_read_pairs", "major_cov", "minor_ref", "minor_read_pairs", "minor_cov")

df_final$sample[1] <- sampleName

# Read the summary of the first mapping
df <- read_tsv(idxstats, col_names = FALSE) %>% 
  # Identify genotype and subtype
  separate(X1, into = c("Subtype", "Reference"), sep = "_", remove = FALSE) %>% 
  # Discard the unmapped reads marked by an *
  filter(X1 != "*")

# Sometimes the mappings stats are completely empty
if (nrow(df) > 0) {
# Count number of reads per subtype
summary <- df %>% 
  group_by(Subtype) %>% 
  summarise(read_pairs = sum(X3)) %>% 
  arrange(desc(read_pairs))

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

# How many read pairs mapped to the minor subtype
major_read_pairs <- summary %>% 
  filter(Subtype == major_tmp) %>% 
  pull(read_pairs)

df_final$major_read_pairs[1] <- major_read_pairs

# Read the depth file from the first mapping. 
# The file can be empty and the reading fails
  cov <- read_tsv(depth, col_names = FALSE) %>% 
  # Filter out the minority subtype
  filter(X1 == major_ref) 


# Reference length
ref_length <- nrow(cov)

# Calculate percentage of positions with a coverage of 5 or more
# Nr. of positions with coverage >= 5
pos <- nrow(
  cov %>% 
    filter(X3 > 4)
)

# Coverage breadth. Convert to integer to be used in a bash if statement later
breadth <- round(pos / ref_length * 100, digits = 2)
breadth_int <- as.integer(pos / ref_length * 100)

df_final$major_cov[1] <- breadth_int

## Minor
minor_tmp <- summary$Subtype[2] 
minor_ref <- df %>% 
  filter(Subtype == minor_tmp) %>% 
  # Choose the reference with most mapped reads
  arrange(desc(X3)) %>% 
  head(n = 1) %>% 
  pull(X1)

df_final$minor_ref[1] <- minor_ref

# How many reads mapped to the minor subtype
minor_read_pairs <- summary %>% 
  filter(Subtype == minor_tmp) %>% 
  pull(read_pairs)

df_final$minor_read_pairs[1] <- minor_read_pairs

# Read the depth file from the first mapping
cov <- read_tsv(depth, col_names = FALSE) %>% 
  # Filter out the minority subtype
  filter(X1 == minor_ref) 

# Reference length
ref_length <- nrow(cov)

# Calculate percentage of positions with a coverage of 5 or more
# Nr. of positions with coverage >= 5
pos <- nrow(
  cov %>% 
    filter(X3 > 4)
)

# Coverage breadth. Convert to integer to be used in a bash if statement later
breadth <- round(pos / ref_length * 100, digits = 2)
breadth_int <- as.integer(pos / ref_length * 100)

df_final$minor_cov[1] <- breadth_int


# Write results
write_csv(df_final, file = paste0(sampleName, ".parsefirstmapping.csv"))
#write_lines(major_ref                         , file = paste0(sampleName, ".major_ref.txt"))
#write_lines(c(minor_ref, minor_reads, breadth_int, breadth), file = paste0(sampleName, ".minor_ref.txt"))

# Write out the fasta sequences
# Read the reference fasta file
fasta <- read.fasta(file = references)
write.fasta(sequences = fasta[major_ref], names = major_ref, file.out = paste0(sampleName, ".", major_ref, "_major.fa"))
write.fasta(sequences = fasta[minor_ref], names = minor_ref, file.out = paste0(sampleName, ".", minor_ref, "_minor.fa"))

}
# Write out sessionInfo() to track versions
# session <- capture.output(sessionInfo())
# write_lines(session, file = "R_versions.txt")