#!/usr/bin/env Rscript
  
library(tidyverse)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: blast_parse.R <prefix> <blast_out> <scaffolds> <references> <agens>", call.=FALSE)
}

prefix     <- args[1]
blast_out  <- args[2]
scaffolds  <- args[3]
references <- args[4]
agens      <- args[5]

# Read the reference fasta file
fasta <- read.fasta(file = references)

# NB! Dette er spades scaffolds som blir lest?? Blir ikke riktig å splitte slik?
# Split scaffolds per subtype
# Read scaffolds fasta
scaffolds_fa <- read.fasta(file = scaffolds)

scaf <- read_tsv(blast_out, col_names = FALSE) %>% 
  # Rename columns
  rename("qseqid" = "X1",
    "sseqid" = "X2",
    "pident" = "X3",
    "length" = "X4",
    "mismatch" = "X5",
    "gapopen" = "X6",
    "qstart" = "X7",
    "qend" = "X8",
    "sstart" = "X9",
    "send" = "X10",
    "evalue" = "X11",
    "bitscore" = "X12")
  
# Separate the genotype from the subject header
# Add a column for the genotype
scaf <- scaf %>% 
  separate(sseqid, into = c("genotype", NA), remove = FALSE) 

# Write out the reformatted blast result
write_csv(scaf, file = paste0(prefix, "_blast_out.csv"))

# Which genotypes are present?
genotypes <- scaf %>% 
  distinct(genotype)

write_tsv(genotypes, paste0(prefix, ".genotypes.tsv"))

# What is the most common reference sequence (sseqid) per genotype?
# Need to choose when two or more subtypes are equally frequent. Choose the one hit by the longest scaffold
ref_info <- scaf %>%
  # Get scaffold length info
  separate(qseqid, c(NA, NA, NA, "sc_length", NA, NA), sep = "_", remove = FALSE) %>% 
  mutate(sc_length = as.numeric(sc_length)) %>% 
  # Select the row with the longest scaffold lengths for each genotype/blast hit
  group_by(genotype) %>% 
  slice_max(order_by = sc_length, n = 1)

# Write out a fasta file for the reference file of each genotype
for (i in 1:nrow(ref_info)) {
  #write_tsv(subtypes[i, 1], file = paste0(prefix, ".", subtypes$sseqid[i], ".txt"), col_names = FALSE)
  write.fasta(sequences = fasta[[ref_info$sseqid[i]]], names = ref_info$sseqid[i], file.out = paste0(prefix, ".", ref_info$sseqid[i], "_ref.fa"))
}

# Split scaffolds per genotype
# First extract the scaffold names that matches the different genotypes from the blast output
split <- scaf %>% 
  group_split(genotype)

# Write one fasta per genotype
for (i in 1:length(split)){
  # Store the scaffold names and corresponding genotype
  tmp <- split[[i]] %>% select(qseqid, genotype) %>% 
    # If the same scaffold has multiple blast hits to the same reference
    # it will cover multiple lines
    distinct()
  # Subset the scaffolds
  geno_fa <- scaffolds_fa[tmp$qseqid]
  # Write to file
  write.fasta(sequences = geno_fa, names = names(geno_fa), file.out = paste0(prefix, ".", tmp$genotype[1], "_scaffolds.fa"))
}

# Identify major and minor reference

# Simply take the best blast hit as determined by Blastn - the top hit - as the majority reference
major_name <- scaf %>% 
  head(n=1) %>% 
  pull(sseqid)

# write the fasta file
write.fasta(sequences = fasta[major_name], names = major_name, file.out = paste0(prefix, ".", major_name, "_major.fa"))

# For the minor, take the second best blast hit that belongs to a different genotype
major_geno <- str_sub(major_name, 1, 1)
minor_name <- scaf %>% 
  filter(str_detect(genotype, paste0("^", major_geno), negate = T)) %>% # Genotype cannot begin with the major_geno
  head(n=1) %>% 
  pull(sseqid)

# write the fasta file
write.fasta(sequences = fasta[minor_name], names = minor_name, file.out = paste0(prefix, ".", minor_name, "_minor.fa"))

# Create a csv file with info on major and minor references and blast hit length
df <- matrix(nrow = 1, ncol = 5)
colnames(df) <- c("sample", "major_ref", "major_length", "minor_ref", "minor_length")
df <- as.data.frame(df)
df$sample <- prefix
df$major_ref <- major_name
df$minor_ref <- minor_name

major_length <- scaf %>% 
  filter(sseqid == major_name) %>% 
  arrange(desc(length)) %>% 
  head(n=1) %>% 
  pull(length)
df$major_length <- major_length

minor_length <- scaf %>% 
  filter(sseqid == minor_name) %>% 
  arrange(desc(length)) %>% 
  head(n=1) %>% 
  pull(length)
df$minor_length <- minor_length

write_csv(df, file = paste0(prefix, ".blastparse.csv"))

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")