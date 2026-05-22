#!/usr/bin/env Rscript
# consensus_distance.R
# Compare a consensus FASTA (from iVar) against its mapping reference FASTA.
# Outputs a TSV with: sample, reference, similarity_pct, n_differences, alignment_length, consensus_length
#
# Usage: consensus_distance.R <consensus.fa> <reference.fa> <output.tsv>

library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: consensus_distance.R <consensus.fa> <reference.fa> <output.tsv>")
}

consensus_file <- args[1]
reference_file <- args[2]
output_file    <- args[3]

# Read sequences
consensus_seqs <- read.fasta(consensus_file, seqtype = "DNA", forceDNAtolower = TRUE)
reference_seqs <- read.fasta(reference_file, seqtype = "DNA", forceDNAtolower = TRUE)

# Take the first sequence from each file
cons_seq <- consensus_seqs[[1]]
ref_seq  <- reference_seqs[[1]]

# Get sequence names
cons_name <- names(consensus_seqs)[1]
ref_name  <- names(reference_seqs)[1]

# Convert to character vectors
cons_chars <- as.character(cons_seq)
ref_chars  <- as.character(ref_seq)

# iVar consensus sequences are the same length as the reference (position-by-position).
# Positions where coverage was too low are filled with 'n'.
# We compare only positions where BOTH sequences have a called base (not 'n' or '-').

# Create masks for callable positions
cons_callable <- !(cons_chars %in% c("n", "-"))
ref_callable  <- !(ref_chars %in% c("n", "-"))
both_callable <- cons_callable & ref_callable

# Number of comparable positions
alignment_length <- sum(both_callable)

if (alignment_length == 0) {
  # No comparable positions — cannot compute distance
  similarity_pct <- NA_real_
  n_differences  <- NA_integer_
} else {
  # Count differences at callable positions
  matches <- cons_chars[both_callable] == ref_chars[both_callable]
  n_differences  <- sum(!matches)
  similarity_pct <- round(sum(matches) / alignment_length * 100, 4)
}

# Consensus length (non-N bases)
consensus_length <- sum(cons_callable)

# Write output
result <- data.frame(
  sample           = cons_name,
  reference        = ref_name,
  similarity_pct   = similarity_pct,
  n_differences    = n_differences,
  alignment_length = alignment_length,
  consensus_length = consensus_length,
  stringsAsFactors = FALSE
)

write.table(result, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
