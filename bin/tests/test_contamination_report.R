#!/usr/bin/env Rscript

# test_contamination_report.R ---------------------------------------------
# Quick task 261016-m5n. bin/contamination_report.R compares SPAdes contig k-mer
# coverage against the read depth expected from index hopping. SPAdes k-mer
# coverage relates to read depth as
#     kmer_cov = read_depth * (L - k + 1) / L
# and the script used to multiply by (L - k + 1) only, with L hard-coded to 150 —
# a 150-fold overstatement that flagged nearly every pair as index hopping.
#
# The script is a CLI, so it is run end to end in a temp dir on a two-sample
# fixture (all-vs-all BLAST hits + fastp JSONs) and its pairs TSV is checked:
#   CR-1  read-depth threshold uses the measured read length
#   CR-2  flag is TRUE below the corrected k-mer threshold, FALSE above it
#         (the pre-fix script flags both)
#   CR-3  read length is the reads-weighted mean across samples
#   CR-4  L < k  =>  threshold cannot be converted, flag is NA (not TRUE/FALSE)
#   CR-5  no fastp JSONs  =>  flag and threshold stay NA
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "contamination_report.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Two cross-sample pairs. S1 is the high-coverage source in both; the recipient
# contig in S2 sits at 2.0x (pair A) and 3.0x (pair B) k-mer coverage.
blast_lines <- c(
  "S1|NODE_1_length_5000_cov_500.0\tS2|NODE_2_length_5000_cov_2.0\t99.9\t5000\t5\t0\t1\t5000\t1\t5000\t0\t9000",
  "S1|NODE_3_length_4000_cov_400.0\tS2|NODE_4_length_4000_cov_3.0\t99.9\t4000\t4\t0\t1\t4000\t1\t4000\t0\t7000"
)

fastp_json <- function(total_reads, read1_mean_length) {
  af <- list(total_reads = total_reads)
  if (!is.na(read1_mean_length)) af$read1_mean_length <- read1_mean_length
  jsonlite::toJSON(list(summary = list(after_filtering = af)), auto_unbox = TRUE)
}

# Run the script in a fresh temp dir; fastp = named list sample -> c(reads, len), or NULL.
run_report <- function(fastp = NULL, kmer_size = 127) {
  wd <- tempfile("contam_test_")
  dir.create(wd)
  writeLines(blast_lines, file.path(wd, "hits.txt"))
  if (!is.null(fastp)) {
    dir.create(file.path(wd, "fastp"))
    for (s in names(fastp)) {
      writeLines(fastp_json(fastp[[s]][1], fastp[[s]][2]),
                 file.path(wd, "fastp", paste0(s, ".fastp.json")))
    }
  }
  # hop_rate 0.001, min_dir_ratio 10, genome_size 9500, kmer_size, min_aln_cov 0.9
  cmd_args <- c(script, "hits.txt", "contamination", "0.001", "10", "9500", kmer_size, "0.9")
  old <- setwd(wd)
  on.exit(setwd(old))
  out <- system2("Rscript", cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(out, sep = "\n")
    fail(sprintf("contamination_report.R exited %s", status))
  }
  read_tsv(file.path(wd, "contamination.contamination_pairs.tsv"),
           show_col_types = FALSE, col_types = cols(.default = col_character()))
}

pair <- function(tsv, recipient_contig) {
  r <- tsv %>% filter(contig1 == recipient_contig | contig2 == recipient_contig)
  if (nrow(r) != 1) fail(sprintf("expected exactly one pair for %s, got %d", recipient_contig, nrow(r)))
  r
}
near <- function(a, b, tol = 1e-3) abs(as.numeric(a) - b) < tol

# ── CR-1 / CR-2: 150 bp reads, k = 127 ────────────────────────────────────────
# run_reads = 2e6, n = 2 -> depth threshold = 2e6 * 0.001 * 150 / (2 * 9500) = 15.7895
# k-mer threshold = 15.7895 * (150 - 127 + 1) / 150 = 2.5263
tsv <- run_report(list(S1 = c(1e6, 150), S2 = c(1e6, 150)))
a <- pair(tsv, "NODE_2_length_5000_cov_2.0")
b <- pair(tsv, "NODE_4_length_4000_cov_3.0")
if (!near(a$index_hop_threshold, 15.7895)) fail(sprintf("CR-1 depth threshold %s != 15.7895", a$index_hop_threshold))
ok("CR-1 read-depth threshold = 15.7895x for 150 bp reads")
if (a$flag_index_hopping != "TRUE") fail("CR-2 recipient at 2.0x (below 2.53x k-mer threshold) must be flagged")
if (b$flag_index_hopping != "FALSE") fail("CR-2 recipient at 3.0x (above 2.53x k-mer threshold) must NOT be flagged (pre-fix threshold was 379x)")
ok("CR-2 flag TRUE at 2.0x, FALSE at 3.0x against the corrected 2.53x k-mer threshold")

# ── CR-3: reads-weighted mean read length ─────────────────────────────────────
# L = (3e6 * 300 + 1e6 * 100) / 4e6 = 250 -> depth = 4e6 * 0.001 * 250 / (2 * 9500) = 52.6316
# k-mer threshold = 52.6316 * 124 / 250 = 26.1053 -> both recipients flagged
tsv <- run_report(list(S1 = c(3e6, 300), S2 = c(1e6, 100)))
a <- pair(tsv, "NODE_2_length_5000_cov_2.0")
if (!near(a$index_hop_threshold, 52.6316)) fail(sprintf("CR-3 depth threshold %s != 52.6316 (weighted L = 250)", a$index_hop_threshold))
ok("CR-3 read length is the reads-weighted mean (250 bp)")

# ── CR-4: reads shorter than k ────────────────────────────────────────────────
tsv <- run_report(list(S1 = c(1e6, 100), S2 = c(1e6, 100)), kmer_size = 127)
if (!all(is.na(tsv$flag_index_hopping))) fail("CR-4 L < k must leave flag_index_hopping NA")
ok("CR-4 L (100) < k (127) -> flag NA, no spurious TRUE/FALSE")

# ── CR-5: no fastp JSONs ──────────────────────────────────────────────────────
tsv <- run_report(NULL)
if (!all(is.na(tsv$flag_index_hopping)) || !all(is.na(tsv$index_hop_threshold)))
  fail("CR-5 without fastp JSONs the hop columns must stay NA")
ok("CR-5 no fastp JSONs -> hop columns NA (unchanged behaviour)")

cat("\nAll contamination_report tests passed.\n")
