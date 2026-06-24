#!/usr/bin/env Rscript
# consensus_distance.R
# Compare a consensus FASTA (from iVar) against its mapping reference FASTA.
# Outputs a TSV with: sample, reference, similarity_pct, n_differences, alignment_length, consensus_length
#
# iVar is run with -aa (all reference positions) and -n N (mask uncovered positions),
# so the consensus is in reference coordinates: N at each uncovered position, one
# character per reference position.  The consensus may be shorter than the reference
# by a few bases when the 3'/5' tail has zero coverage and was trimmed before output.
# Some iVar versions also output '-' at positions with a confirmed consensus deletion.
#
# A global Needleman-Wunsch alignment (Biostrings/pwalign pairwiseAlignment) is used
# so that internal insertions and deletions are properly placed rather than causing a
# frame-shift in a position-by-position comparison.  Existing '-' characters are
# stripped from the consensus before alignment (pairwiseAlignment requires ungapped
# input); the aligner re-places them at the optimal positions.  Gap columns introduced
# by the aligner count as differences.  N / zero-coverage columns are excluded from
# the denominator (alignment_length) exactly as before.
#
# Usage: consensus_distance.R <consensus.fa> <reference.fa> <output.tsv>

suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: consensus_distance.R <consensus.fa> <reference.fa> <output.tsv>")
}

consensus_file <- args[1]
reference_file <- args[2]
output_file    <- args[3]

# Read sequences
cons_raw <- readDNAStringSet(consensus_file)[[1]]
ref_raw  <- readDNAStringSet(reference_file)[[1]]

cons_name <- names(readDNAStringSet(consensus_file))[1]
ref_name  <- names(readDNAStringSet(reference_file))[1]

# Consensus length: non-N, non-gap called bases in the full consensus (before alignment)
cons_str_full    <- toupper(as.character(cons_raw))
cons_chars_full  <- strsplit(cons_str_full, "")[[1]]
consensus_length <- sum(!(cons_chars_full %in% c("N", "-")))

# Strip existing '-' (confirmed iVar deletions) before alignment — pairwiseAlignment
# requires ungapped input; the aligner will re-place gaps at optimal positions.
cons_seq <- DNAString(gsub("-", "", cons_str_full))
ref_seq  <- DNAString(toupper(as.character(ref_raw)))

# Global pairwise alignment (Needleman-Wunsch)
submat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE, type = "DNA")
aln <- pairwiseAlignment(
  cons_seq, ref_seq,
  type               = "global",
  substitutionMatrix = submat,
  gapOpening         = 10,
  gapExtension       = 0.5
)

# as.character() on AlignedXStringSet returns the aligned string with internal gap
# characters ('-') included; terminal overhangs (positions beyond the shorter sequence)
# are not represented and are therefore not counted — consistent with the prior
# min-length cap behaviour.
cons_aln <- strsplit(as.character(pattern(aln)), "")[[1]]
ref_aln  <- strsplit(as.character(subject(aln)), "")[[1]]

# Callable mask: exclude positions where either sequence has N (uncovered / zero-coverage)
both_callable    <- cons_aln != "N" & ref_aln != "N"
alignment_length <- sum(both_callable)

if (alignment_length == 0) {
  similarity_pct <- NA_real_
  n_differences  <- NA_integer_
} else {
  # Gap characters ('-') at callable positions count as differences (indels)
  n_differences  <- sum(cons_aln[both_callable] != ref_aln[both_callable])
  similarity_pct <- round((alignment_length - n_differences) / alignment_length * 100, 4)
}

# Write output — column names are identical to the previous implementation
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
