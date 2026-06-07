#!/usr/bin/env Rscript

# test_summarize_denovo.R -------------------------------------------------
# Fixture tests for the de novo confirmation downgrade layer used by
# bin/summarize.R. As of plan 04-01 the layer lives in bin/denovo_layer.R as a
# sourceable apply_denovo_layer(), so these tests exercise the REAL function
# rather than an inline re-implementation. This kills the test-drift weakness:
# if summarize.R's downgrade behaviour changes, it changes here (same code path).
#
# The test asserts three things:
#
#   1. summarize.R still parses (a single smoke check — the brittle source-text
#      greps were removed in 04-01 now that the real function is under test).
#
#   2. read_blast_out_frame() behaviour: keyed long frame from per-contig CSVs,
#      with a typed-empty fallback that never aborts on a no-de-novo run.
#
#   3. apply_denovo_layer() behaviour on fixtures: confirmed/refuted/unconfirmed/
#      not_evaluated/NA sentinels + downgrade-only minor_typable flip, PLUS a
#      flag-OFF legacy-reproduction differential against a committed golden CSV
#      (D-05/D-06) — comparing the pre-existing column subset, EXCLUDING the
#      additive minor_denovo_status column.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "denovo_confirm.R"))
source(file.path(bin_dir, "denovo_layer.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

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

# --- Block 1: summarize.R smoke check -------------------------------------
# The downgrade logic is now exercised through the REAL apply_denovo_layer()
# below, so the brittle source-text greps were dropped (04-01). We keep only a
# parse() smoke check to catch a syntactically broken script.
parse(file.path(bin_dir, "summarize.R"))
ok("summarize.R parses")

# Sanity: summarize.R sources and calls the extracted layer (key_links contract).
src <- readLines(file.path(bin_dir, "summarize.R"))
if (!any(grepl("source\\(\"denovo_layer.R\"\\)", src))) fail("summarize.R missing source(denovo_layer.R)")
if (!any(grepl("apply_denovo_layer\\(", src))) fail("summarize.R missing apply_denovo_layer() call")
ok("summarize.R sources + calls apply_denovo_layer()")

# --- Block 2: read_blast_out_frame() fixtures -----------------------------

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
ok("Test1: keyed long frame with sampleName")

# Test 2: empty guard -> typed-empty tibble, never abort
empty_dir <- file.path(tmp, "denovo_empty"); dir.create(empty_dir)
df_e <- read_blast_out_frame(empty_dir)
if (nrow(df_e) != 0) fail("empty read should yield zero rows")
if (!all(c("sampleName", "subtype", "sc_length", "kmer_cov", "pident") %in% names(df_e))) {
  fail("empty read missing typed columns")
}
ok("Test2: empty guard yields typed-empty tibble")

# --- Block 3: apply_denovo_layer() behaviour (the REAL function) ----------

# Test 3 (confirmed): substantial minor contig -> confirmed, minor_typable unchanged
final3 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo3 <- bind_rows(mk_blast("S1", "1a", 5000, 10, 99), mk_blast("S1", "2b", 2949, 5, 99)) %>%
  mutate(sampleName = "S1")
r3 <- apply_denovo_layer(final3, bo3, TRUE)
if (r3$minor_denovo_status != "confirmed_by_denovo") fail("Test3 expected confirmed_by_denovo")
if (r3$minor_typable != "YES") fail("Test3 minor_typable must stay YES")
ok("Test3 (CONF-01): confirmed, minor_typable unchanged")

# Test 4 (refuted, downgrade): substantial major, no substantial minor
final4 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo4 <- mk_blast("S1", "1a", 5000, 10, 99) %>% mutate(sampleName = "S1")
r4 <- apply_denovo_layer(final4, bo4, TRUE)
if (r4$minor_denovo_status != "refuted") fail("Test4 expected refuted")
if (r4$minor_typable != "NO") fail("Test4 minor_typable must flip to NO")
if (is.na(r4$Minor_reference) || r4$Minor_reference != "2b_ACC") fail("Test4 Minor_reference must stay populated")
ok("Test4 (CONF-02): refuted downgrades minor_typable, Minor_* intact")

# Test 5 (unconfirmed): no substantial contig at all
final5 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = "2b_ACC", minor_typable = "YES")
bo5 <- mk_blast("S1", "1a", 300, 1, 99) %>% mutate(sampleName = "S1")  # tiny noise contig
r5 <- apply_denovo_layer(final5, bo5, TRUE)
if (r5$minor_denovo_status != "unconfirmed") fail("Test5 expected unconfirmed")
if (r5$minor_typable != "YES") fail("Test5 minor_typable must stay YES (not suppressed)")
ok("Test5 (CONF-03): unconfirmed, minor not suppressed")

# Test 6 (flag OFF): not_evaluated when candidate, NA when none; minor_typable legacy
final6 <- tibble(sampleName = c("S1", "S2"),
                 Major_reference = c("1a_ACC", "1a_ACC"),
                 Minor_reference = c("2b_ACC", NA_character_),
                 minor_typable = c("YES", "NO"))
r6 <- apply_denovo_layer(final6, bo4, FALSE)
if (r6$minor_denovo_status[1] != "not_evaluated") fail("Test6 candidate must be not_evaluated when flag OFF")
if (!is.na(r6$minor_denovo_status[2])) fail("Test6 no-candidate must be NA when flag OFF")
if (!identical(r6$minor_typable, c("YES", "NO"))) fail("Test6 minor_typable must be untouched when flag OFF")
ok("Test6 (CONF-07/D-16): flag OFF bypasses layer")

# Test 7 (no candidate, flag ON): NA status regardless
final7 <- tibble(sampleName = "S1", Major_reference = "1a_ACC",
                 Minor_reference = NA_character_, minor_typable = "NO")
r7 <- apply_denovo_layer(final7, bo4, TRUE)
if (!is.na(r7$minor_denovo_status)) fail("Test7 no candidate must be NA status (D-12)")
ok("Test7 (D-12): no minor candidate -> NA status")

# --- Block 4: flag-OFF legacy-reproduction differential (D-05/D-06) -------
# Run the layer with flag OFF on a small fixture (one candidate-minor row, one
# no-candidate row), select the PRE-EXISTING (legacy) columns ONLY — i.e. every
# column EXCEPT the additive minor_denovo_status — and assert byte-equality to a
# committed golden baseline. This proves OFF reproduces legacy output on the
# legacy column subset; the additive status column never breaks the differential.
golden_path <- file.path(bin_dir, "tests", "fixtures", "flagoff_golden.csv")

final_off <- tibble(
  sampleName      = c("S1", "S2"),
  Major_reference = c("1a_ACC", "3a_ACC"),
  Minor_reference = c("2b_ACC", NA_character_),
  major_typable   = c("YES", "YES"),
  minor_typable   = c("YES", "NO")
)
r_off <- apply_denovo_layer(final_off, bo4, FALSE)

# Legacy column subset = all columns except the additive minor_denovo_status.
legacy_cols <- setdiff(names(r_off), "minor_denovo_status")
r_off_legacy <- r_off %>% select(all_of(legacy_cols))

# minor_typable must retain its pure legacy value under flag-OFF.
if (!identical(r_off_legacy$minor_typable, c("YES", "NO"))) {
  fail("flag-OFF minor_typable must equal pure legacy value")
}

if (!file.exists(golden_path)) {
  # First-run bootstrap: write the golden baseline once for the reviewer to
  # inspect and commit (D-06). Subsequent runs assert against it.
  dir.create(dirname(golden_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(r_off_legacy, golden_path)
  fail(paste0("golden baseline written to ", golden_path,
              " — inspect and commit, then re-run (D-06 bootstrap)"))
}

golden <- read_csv(golden_path, show_col_types = FALSE)
# Compare the legacy subset on both sides, EXCLUDING minor_denovo_status by
# construction (it is not in legacy_cols and not written to the golden file).
golden_legacy <- golden %>% select(any_of(legacy_cols))
if (!isTRUE(all.equal(as.data.frame(r_off_legacy), as.data.frame(golden_legacy)))) {
  fail("flag-OFF legacy-column output diverged from committed golden baseline")
}
if ("minor_denovo_status" %in% names(golden)) {
  fail("golden baseline must NOT contain minor_denovo_status (additive column excluded)")
}
ok("Block4 (D-05/D-06): flag-OFF reproduces committed legacy golden baseline")

cat("\nALL PASS\n")
