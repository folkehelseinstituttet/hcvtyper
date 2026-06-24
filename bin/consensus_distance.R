#!/usr/bin/env Rscript
# consensus_distance.R
# Compare a consensus FASTA (from iVar) against its mapping reference FASTA.
# Outputs a TSV with: sample, reference, similarity_pct, n_differences, alignment_length, consensus_length
#
# iVar is run with -aa (all reference positions) and -n N (mask uncovered positions),
# so the consensus is in reference coordinates: N at each uncovered position, one
# character per reference position.  The consensus may be shorter than the reference
# by a few bases when the 3'/5' tail has zero coverage and was trimmed before output.
#
# We compare only positions where BOTH sequences have a called base (not 'n' or '-').
# Lengths are capped at min(cons_len, ref_len) before masking so that R never recycles
# the shorter logical vector into a longer one — recycling silently produces out-of-bounds
# NA values that propagate through sum() and make the result NA even when there are many
# valid callable pairs.
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

# Consensus length (non-N bases) — computed over the full consensus before truncation
cons_callable_full <- !(cons_chars %in% c("n", "-"))
consensus_length   <- sum(cons_callable_full)

# Cap at the shorter length before building the callable masks.
# Without this, `cons_callable & ref_callable` recycles the shorter vector, placing
# recycled TRUE values beyond the end of the shorter vector; indexing with those
# positions returns NA, which propagates through sum() to produce a spurious NA result.
compare_len <- min(length(cons_chars), length(ref_chars))
cons_cmp    <- cons_chars[seq_len(compare_len)]
ref_cmp     <- ref_chars[seq_len(compare_len)]

cons_callable <- !(cons_cmp %in% c("n", "-"))
ref_callable  <- !(ref_cmp  %in% c("n", "-"))
both_callable <- cons_callable & ref_callable   # same length — no recycling

alignment_length <- sum(both_callable)

if (alignment_length == 0) {
  similarity_pct <- NA_real_
  n_differences  <- NA_integer_
} else {
  matches        <- cons_cmp[both_callable] == ref_cmp[both_callable]
  n_differences  <- sum(!matches)
  similarity_pct <- round(sum(matches) / alignment_length * 100, 4)
}

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
