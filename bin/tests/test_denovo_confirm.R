#!/usr/bin/env Rscript

# test_denovo_confirm.R ---------------------------------------------------
# Self-contained unit tests for the pure classify_minor_denovo() helper in
# bin/denovo_confirm.R. Exercises every decision branch (confirm / refute /
# unconfirmed / NA) plus the empty-input guard, sc_length (not alignment
# length) substantiality, and genotype-level (not subtype) matching.
#
# Run from the repo root:
#   Rscript bin/tests/test_denovo_confirm.R
# Exits 0 and prints "ALL PASS" when all assertions hold; stop()s otherwise.
# Phase 4's TEST-01 will fold these fixtures into the formal test harness.
# -------------------------------------------------------------------------

suppressMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
bin_dir <- if (length(this_file) == 1) dirname(dirname(normalizePath(this_file))) else "bin"

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "denovo_confirm.R"))

# Test 1 (CONF-01 confirm): substantial minor-genotype contig.
df_minor <- tibble(qseqid = "N1", subtype = "2b", pident = 91.32, length = 300,
                   evalue = 0, bitscore = 100, sc_length = 6426, kmer_cov = 4.79)
stopifnot(classify_minor_denovo(df_minor, "1", "2") == "confirmed_by_denovo")

# Test 2 (CONF-02 refute): substantial major-genotype-1 contig, no minor.
df_major <- tibble(qseqid = "N1", subtype = "1a", pident = 99, length = 9000,
                   evalue = 0, bitscore = 999, sc_length = 9189, kmer_cov = 15056)
stopifnot(classify_minor_denovo(df_major, "1", "4") == "refuted")

# Test 3 (CONF-03): major contig fails length floor BUT substantial minor present -> confirm.
df_shortmaj_minor <- tibble(
  qseqid   = c("N_maj", "N_min"),
  subtype  = c("3a",    "1a"),
  pident   = c(99,      95),
  length   = c(800,     8000),
  evalue   = c(0,       0),
  bitscore = c(900,     900),
  sc_length = c(972,    8500),
  kmer_cov  = c(15056,  10)
)
stopifnot(classify_minor_denovo(df_shortmaj_minor, "3", "1") == "confirmed_by_denovo")

# Test 3b: no substantial major AND no substantial minor -> unconfirmed (never refute).
df_allshort <- tibble(qseqid = "N1", subtype = "1a", pident = 99, length = 300,
                      evalue = 0, bitscore = 100, sc_length = 300, kmer_cov = 1.0)
stopifnot(classify_minor_denovo(df_allshort, "3", "1") == "unconfirmed")

# Test 4 (empty/NULL DoS guard): NULL and zero-row both -> unconfirmed, never error.
stopifnot(classify_minor_denovo(NULL, "1", "2") == "unconfirmed")
df_empty <- df_minor[0, ]
stopifnot(classify_minor_denovo(df_empty, "1", "2") == "unconfirmed")

# Test 5 (no minor candidate): minor_geno = NA -> NA_character_.
stopifnot(is.na(classify_minor_denovo(df_minor, "1", NA_character_)))

# Test 6 (substantial uses sc_length, not alignment length):
# sc_length=6426 with tiny alignment length=300 still substantial.
df_bigsc_smallaln <- tibble(qseqid = "N1", subtype = "2b", pident = 91, length = 300,
                            evalue = 0, bitscore = 100, sc_length = 6426, kmer_cov = 4.79)
stopifnot(classify_minor_denovo(df_bigsc_smallaln, "1", "2") == "confirmed_by_denovo")
# sc_length=300 with huge alignment length fails regardless.
df_smallsc_bigaln <- tibble(qseqid = "N1", subtype = "2b", pident = 91, length = 9000,
                            evalue = 0, bitscore = 999, sc_length = 300, kmer_cov = 4.79)
stopifnot(classify_minor_denovo(df_smallsc_bigaln, "1", "2") == "unconfirmed")

# Test 7 (genotype-level match, CONF-05): subtype "3a" matches minor_geno "3".
df_3a <- tibble(qseqid = "N1", subtype = "3a", pident = 95, length = 8000,
                evalue = 0, bitscore = 900, sc_length = 8000, kmer_cov = 10)
stopifnot(classify_minor_denovo(df_3a, "1", "3") == "confirmed_by_denovo")
# And does NOT confirm at genotype level when genotypes differ ("3a" vs minor "2").
stopifnot(classify_minor_denovo(df_3a, "3", "2") == "refuted")

cat("ALL PASS\n")
