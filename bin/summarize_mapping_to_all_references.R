#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: summarize_mapping_to_all_references.R <idxstats file> <depth file> <sample name>", call.=FALSE)
}

idxstats   <- args[1]
depth      <- args[2]
sampleName <- args[3]

# Read the summary of the first mapping
df <- read_tsv(idxstats, col_names = FALSE) %>% 
  # Identify genotype and subtype
  separate(X1, into = c("Subtype", "Reference"), sep = "_", remove = FALSE) 

# Count number of reads per subtype
summary <- df %>% 
  group_by(Subtype) %>% 
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

## Minor
minor_tmp <- summary$Subtype[2] 
minor_ref <- df %>% 
  filter(Subtype == minor_tmp) %>% 
  # Choose the reference with most mapped reads
  arrange(desc(X3)) %>% 
  head(n = 1) %>% 
  pull(X1)

# How many reads mapped to the minor subtype
minor_reads <- summary %>% 
  filter(Subtype == minor_tmp) %>% 
  pull(reads)

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

# Write files
write_lines(major_ref                         , file = paste0(sampleName, ".major_ref.txt"))
write_lines(c(minor_ref, minor_reads, breadth_int, breadth), file = paste0(sampleName, ".minor_ref.txt"))

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")