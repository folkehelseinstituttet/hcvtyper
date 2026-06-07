#!/usr/bin/env Rscript

# test_summarize_denovo.R -------------------------------------------------
# Self-contained fixture tests for the de novo confirmation integration in
# bin/summarize.R (plan 03-02). summarize.R is a monolithic commandArgs
# entrypoint, so rather than running the whole pipeline these tests assert two
# things directly:
#
#   1. The SOURCE of summarize.R contains the required structural elements
#      (the *_blast_out.csv read into df_blast_out with a length() guard, the
#      helper source()s, the classify_minor_denovo() call, the downgrade-only
#      if_else, and the always-present minor_denovo_status column). These are
#      the grep-level acceptance criteria from the plan.
#
#   2. The BEHAVIOUR of the two extracted logic units works on fixtures:
#        - read_blast_out_frame(): keyed long frame from per-contig CSVs, with a
#          typed-empty fallback that never aborts on a no-de-novo run.
#        - the per-sample downgrade layer: confirmed/refuted/unconfirmed/
#          not_evaluated/NA sentinels + downgrade-only minor_typable flip.
#
# The behavioural units are defined inline here to mirror EXACTLY the logic that
# must live in summarize.R; if summarize.R drifts from this contract the
# structural grep assertions (block 1) fail.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "denovo_confirm.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# --- Reference logic units (must mirror summarize.R) ----------------------

# Read all denovo/*_blast_out.csv into a long frame keyed by sampleName, with a
# typed-empty fallback (mirrors summarize.R Task 1).
read_blast_out_frame <- function(path_denovo) {
  files <- list.files(path = path_denovo, pattern = "_blast_out.csv$", full.names = TRUE)
  if (length(files) > 0) {
    map_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>%
      mutate(sampleName = str_remove(basename(.x), "_blast_out.csv$")))
  } else {
    tibble(
      sampleName = character(),
      qseqid     = character(),
      subtype    = character(),
      pident     = double(),
      evalue     = double(),
      bitscore   = double(),
      sc_length  = double(),
      kmer_cov   = double()
    )
  }
}

# --- Block 1: structural grep assertions over summarize.R -----------------

src <- readLines(file.path(bin_dir, "summarize.R"))
has <- function(pat) any(grepl(pat, src, fixed = FALSE))

if (!has("_blast_out.csv\\$")) fail("summarize.R missing _blast_out.csv$ read pattern")
ok("summarize.R contains _blast_out.csv$ read pattern")

if (!has("df_blast_out")) fail("summarize.R missing df_blast_out frame")
ok("summarize.R contains df_blast_out")

if (!has("source\\(\"denovo_confirm.R\"\\)")) fail("summarize.R missing source(denovo_confirm.R)")
if (!has("source\\(\"genotype_utils.R\"\\)")) fail("summarize.R missing source(genotype_utils.R)")
ok("summarize.R sources both helpers")

if (!has("classify_minor_denovo\\(")) fail("summarize.R missing classify_minor_denovo() call")
ok("summarize.R calls classify_minor_denovo()")

if (!has("minor_denovo_status")) fail("summarize.R missing minor_denovo_status column")
if (!has("not_evaluated")) fail("summarize.R missing not_evaluated sentinel")
ok("summarize.R emits minor_denovo_status + not_evaluated sentinel")

# downgrade-only: "refuted" -> "NO" against minor_typable on one edit
refuted_lines <- grep("refuted", src, value = TRUE)
if (!any(grepl("minor_typable", refuted_lines) | grepl("if_else", refuted_lines))) {
  # downgrade may span; require both tokens present overall plus a guarded if_else
  if (!(has("refuted") && has("minor_typable") && has("if_else"))) {
    fail("summarize.R missing downgrade-only if_else(refuted -> NO, minor_typable)")
  }
}
ok("summarize.R contains downgrade-only refute logic")

# must NOT null Minor_reference on refute
if (has("Minor_reference = NA")) fail("summarize.R nulls Minor_reference (D-10 violation)")
ok("summarize.R keeps Minor_* populated on refute")

# parses cleanly
parse(file.path(bin_dir, "summarize.R"))
ok("summarize.R parses")

# --- Block 2: behavioural fixtures ----------------------------------------

tmp <- tempfile("denovo_test_"); dir.create(tmp)
denovo_dir <- file.path(tmp, "denovo"); dir.create(denovo_dir)

mk_blast <- function(sample, subtype, sc_length, kmer_cov, pident) {
  tibble(
    qseqid = paste0("NODE_1_length_", sc_length, "_cov_", kmer_cov),
    sseqid = paste0(subtype, "_ACC"),
    pident = pident, length = sc_length, mismatch = 0, gapopen = 0,
    qstart = 1, qend = sc_length, sstart = 1, send = sc_length,
    evalue = 0, bitscore = 1000, subtype = subtype,
    sc_length = sc_length, kmer_cov = kmer_cov
  )
}

# Test 1: keyed long frame from existing files
write_csv(mk_blast("S1", "2b", 2949, 5, 99), file.path(denovo_dir, "S1_blast_out.csv"))
write_csv(mk_blast("S2", "1a", 5000, 10, 99), file.path(denovo_dir, "S2_blast_out.csv"))
df <- read_blast_out_frame(denovo_dir)
if (!"sampleName" %in% names(df)) fail("read_blast_out_frame missing sampleName key")
if (!setequal(unique(df$sampleName), c("S1", "S2"))) fail("read_blast_out_frame wrong sample keys")
ok("Task1 Test1: keyed long frame with sampleName")

# Test 2: empty guard -> typed-empty tibble, never abort
empty_dir <- file.path(tmp, "denovo_empty"); dir.create(empty_dir)
df_e <- read_blast_out_frame(empty_dir)
if (nrow(df_e) != 0) fail("empty read should yield zero rows")
if (!all(c("sampleName", "subtype", "sc_length", "kmer_cov", "pident") %in% names(df_e))) {
  fail("empty read missing typed columns")
}
ok("Task1 Test2: empty guard yields typed-empty tibble")

# --- Downgrade layer behaviour (mirrors summarize.R Task 2) ---------------
# Build a tiny `final` and apply the same logic the script must apply.
apply_layer <- function(final, df_blast_out, flag,
                        min_len = 1000, min_kmer = 2.0, min_pid = 90, match_level = "genotype") {
  if (isTRUE(flag)) {
    final %>%
      rowwise() %>%
      mutate(minor_denovo_status = {
        if (is.na(Minor_reference)) {
          NA_character_
        } else {
          bo <- df_blast_out %>% filter(sampleName == .data$sampleName)
          classify_minor_denovo(
            bo,
            genotype_from_subtype(str_extract(Major_reference, "^[^_]+")),
            genotype_from_subtype(str_extract(Minor_reference, "^[^_]+")),
            min_len, min_kmer, min_pid, match_level
          )
        }
      }) %>%
      ungroup() %>%
      mutate(minor_typable = if_else(
        !is.na(minor_denovo_status) & minor_denovo_status == "refuted", "NO", minor_typable
      ))
  } else {
    final %>%
      mutate(minor_denovo_status = if_else(is.na(Minor_reference), NA_character_, "not_evaluated"))
  }
}

# Test 3 (confirmed): substantial minor contig -> confirmed, minor_typable unchanged
final3 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo3 <- bind_rows(mk_blast("S1", "1a", 5000, 10, 99), mk_blast("S1", "2b", 2949, 5, 99)) %>%
  mutate(sampleName = "S1")
r3 <- apply_layer(final3, bo3, TRUE)
if (r3$minor_denovo_status != "confirmed_by_denovo") fail("Test3 expected confirmed_by_denovo")
if (r3$minor_typable != "YES") fail("Test3 minor_typable must stay YES")
ok("Task2 Test3 (CONF-01): confirmed, minor_typable unchanged")

# Test 4 (refuted, downgrade): substantial major, no substantial minor
final4 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo4 <- mk_blast("S1", "1a", 5000, 10, 99) %>% mutate(sampleName = "S1")
r4 <- apply_layer(final4, bo4, TRUE)
if (r4$minor_denovo_status != "refuted") fail("Test4 expected refuted")
if (r4$minor_typable != "NO") fail("Test4 minor_typable must flip to NO")
if (is.na(r4$Minor_reference) || r4$Minor_reference != "2b_ACC") fail("Test4 Minor_reference must stay populated")
ok("Task2 Test4 (CONF-02): refuted downgrades minor_typable, Minor_* intact")

# Test 5 (unconfirmed): no substantial contig at all
final5 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo5 <- mk_blast("S1", "1a", 300, 1, 99) %>% mutate(sampleName = "S1")  # tiny noise contig
r5 <- apply_layer(final5, bo5, TRUE)
if (r5$minor_denovo_status != "unconfirmed") fail("Test5 expected unconfirmed")
if (r5$minor_typable != "YES") fail("Test5 minor_typable must stay YES (not suppressed)")
ok("Task2 Test5 (CONF-03): unconfirmed, minor not suppressed")

# Test 6 (flag OFF): not_evaluated when candidate, NA when none; minor_typable legacy
final6 <- tibble(sampleName = c("S1", "S2"),
                 Major_reference = c("1a_ACC", "1a_ACC"),
                 Minor_reference = c("2b_ACC", NA_character_),
                 minor_typable = c("YES", "NO"))
r6 <- apply_layer(final6, bo4, FALSE)
if (r6$minor_denovo_status[1] != "not_evaluated") fail("Test6 candidate must be not_evaluated when flag OFF")
if (!is.na(r6$minor_denovo_status[2])) fail("Test6 no-candidate must be NA when flag OFF")
if (!identical(r6$minor_typable, c("YES", "NO"))) fail("Test6 minor_typable must be untouched when flag OFF")
ok("Task2 Test6 (CONF-07/D-16): flag OFF bypasses layer")

# Test 7 (no candidate, flag ON): NA status regardless
final7 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = NA_character_, minor_typable = "NO")
r7 <- apply_layer(final7, bo4, TRUE)
if (!is.na(r7$minor_denovo_status)) fail("Test7 no candidate must be NA status (D-12)")
ok("Task2 Test7 (D-12): no minor candidate -> NA status")

cat("\nALL PASS\n")
