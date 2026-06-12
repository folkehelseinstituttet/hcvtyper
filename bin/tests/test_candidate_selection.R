#!/usr/bin/env Rscript

# test_candidate_selection.R ---------------------------------------------
# Subprocess-contract test (Phase 6, REFSEL-01) for the REAL neutral-selection
# entrypoint bin/summarize_mapping_to_all_references.R. Modelled on
# test_major_gate.R: we invoke the script via system2("Rscript", ...) on small
# synthetic idxstats / depth / refs.fa fixtures built inline in a tempdir, and
# assert on what it writes:
#   * the NEW long-format candidate table (<sample>.candidates.csv): one row per
#     selected candidate, columns
#       sample, candidate_rank, candidate_ref, candidate_subtype,
#       candidate_genotype, candidate_reads, candidate_cov, confirmation_status
#   * the reconstructed legacy 10-column wide CSV (<sample>.parsefirstmapping.csv)
#     with rank1 -> major_*, rank2 -> minor_*, all 10 columns always present.
#
# Behaviours asserted (06-CONTEXT D-01/D-02/D-03/D-05):
#   (a) Neutral top-N-by-reads ranking + distinct-subtype dedup: given subtypes
#       {1a:1000+800 reads (two refs), 1b:600, 3a:50}, candidate_rank 1 is the
#       1a ref with the most reads, candidate_rank 2 is the 1b ref. The
#       same-subtype second 1a ref is NOT selected (distinct-subtype dedup, D-02).
#   (b) 1a/1b preservation: a 1a-dominant + 1b-second sample yields TWO candidates
#       (both genotype 1, different subtypes) — dedup does not collapse them.
#   (c) No validity filtering (D-05): a 1a-dominant + 2k1b-second sample still
#       selects a 2nd candidate (the old is_valid_minor() 2k1b suppression is gone).
#   (d) Long-format contract: the candidate table has exactly the column set above,
#       with n_candidates rows (or fewer if fewer distinct subtypes exist).
#   (e) Shim reconstruction: the legacy 10-column wide CSV has rank1 in major_*,
#       rank2 in minor_*; all 10 columns present even on a single-candidate sample
#       (minor_* NA-filled), gate_flag defaults to "no_mapping" on empty input.
#
# The script sources genotype_utils.R cwd-relative and writes its outputs into the
# cwd, so each case runs in its own tempdir with genotype_utils.R copied in and
# absolute fixture paths passed as args. genotype_utils.R is located relative to
# this test file (the staged-path contract) — never via an absolute hardcode.
#
# RED until Task 3 lands the script rewrite (long-format emit + neutral ranking).
#
# Run from any cwd:
#   Rscript bin/tests/test_candidate_selection.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "summarize_mapping_to_all_references.R")
geno_src <- file.path(bin_dir, "genotype_utils.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Build idxstats / depth / refs.fa fixtures inline and invoke the REAL script in
# an isolated per-case tempdir with genotype_utils.R copied in. `refs` is a named
# list ref-name -> integer reads; `covs` is ref-name -> per-position coverage
# value (constant across `genome_len` positions). Returns
# list(exit, candidates, legacy) read back from the two output CSVs (or NULL).
run_select <- function(case, sampleName, refs, covs, n_candidates,
                       minRead = 100, minCov = 10, genome_len = 100) {
  wd <- tempfile(paste0("cand_", case, "_")); dir.create(wd)
  file.copy(geno_src, file.path(wd, "genotype_utils.R"), overwrite = TRUE)

  # idxstats: <name> <length> <mapped_reads> <unmapped>; '*' unmapped sentinel row.
  idx_rows <- map_chr(names(refs), function(r) {
    paste(r, genome_len, refs[[r]], 0, sep = "\t")
  })
  idx_rows <- c(idx_rows, paste("*", 0, 0, 1000, sep = "\t"))
  idx_path <- file.path(wd, paste0(case, "_idxstats.tsv"))
  writeLines(idx_rows, idx_path)

  # depth: <name> <pos> <coverage>; constant coverage per ref over genome_len.
  depth_lines <- character(0)
  for (r in names(covs)) {
    depth_lines <- c(depth_lines,
      paste(r, seq_len(genome_len), covs[[r]], sep = "\t"))
  }
  depth_path <- file.path(wd, paste0(case, "_depth.tsv"))
  if (length(depth_lines) == 0) depth_lines <- character(0)
  writeLines(depth_lines, depth_path)

  # refs.fa: a trivial sequence per reference name.
  fa_lines <- unlist(lapply(names(refs), function(r) {
    c(paste0(">", r), strrep("A", 60))
  }))
  refs_path <- file.path(wd, paste0(case, "_refs.fa"))
  writeLines(fa_lines, refs_path)

  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(idx_path), shQuote(depth_path),
      shQuote(sampleName), shQuote(refs_path), minRead, minCov, n_candidates),
    stdout = FALSE, stderr = FALSE
  )

  cand_path   <- file.path(wd, paste0(sampleName, ".candidates.csv"))
  legacy_path <- file.path(wd, paste0(sampleName, ".parsefirstmapping.csv"))
  candidates <- if (file.exists(cand_path))   read_csv(cand_path,   show_col_types = FALSE) else NULL
  legacy     <- if (file.exists(legacy_path)) read_csv(legacy_path, show_col_types = FALSE) else NULL
  list(exit = exit, candidates = candidates, legacy = legacy)
}

LONG_COLS <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
               "candidate_genotype", "candidate_reads", "candidate_cov",
               "confirmation_status")
LEGACY_COLS <- c("sample", "total_mapped_reads", "major_ref", "major_reads",
                 "major_cov", "minor_ref", "minor_reads", "minor_cov",
                 "minor_call", "gate_flag")

# --- Case (a) neutral ranking + distinct-subtype dedup (D-01/D-02/D-03) -----
# Subtypes: 1a has two refs (1000 + 800), 1b one ref (600), 3a one ref (50).
# Expected (n=2): rank1 = 1a_TOP (most reads in the top subtype), rank2 = 1b_ONE
# (top per distinct subtype, top-2 subtypes by total reads). The second 1a ref is
# NOT selected.
r_a <- run_select(
  "rank",
  "RANK",
  refs = list("1a_TOP" = 1000, "1a_SEC" = 800, "1b_ONE" = 600, "3a_LOW" = 50),
  covs = list("1a_TOP" = 50, "1a_SEC" = 50, "1b_ONE" = 50, "3a_LOW" = 50),
  n_candidates = 2
)
if (is.null(r_a$candidates)) fail("rank: no long-format candidates CSV written")
if (!all(LONG_COLS %in% colnames(r_a$candidates)))
  fail(paste("rank: long-format columns missing; got:",
             paste(colnames(r_a$candidates), collapse = ",")))
if (nrow(r_a$candidates) != 2)
  fail(paste("rank: expected 2 candidate rows, got", nrow(r_a$candidates)))
cand1 <- r_a$candidates %>% filter(candidate_rank == 1)
cand2 <- r_a$candidates %>% filter(candidate_rank == 2)
if (nrow(cand1) != 1 || cand1$candidate_ref[1] != "1a_TOP")
  fail(paste("rank: candidate_rank 1 should be 1a_TOP, got",
             paste(cand1$candidate_ref, collapse = ",")))
if (nrow(cand2) != 1 || cand2$candidate_ref[1] != "1b_ONE")
  fail(paste("rank: candidate_rank 2 should be 1b_ONE (distinct-subtype dedup), got",
             paste(cand2$candidate_ref, collapse = ",")))
if ("1a_SEC" %in% r_a$candidates$candidate_ref)
  fail("rank: same-subtype second 1a ref must NOT be selected (D-02 dedup)")
ok("rank -> rank1=1a_TOP, rank2=1b_ONE, distinct-subtype dedup drops 1a_SEC")

# --- Case (b) 1a/1b preservation (both genotype 1, different subtypes) -------
r_b <- run_select(
  "oneab",
  "ONEAB",
  refs = list("1a_DOM" = 900, "1b_SUB" = 400),
  covs = list("1a_DOM" = 50, "1b_SUB" = 50),
  n_candidates = 2
)
if (is.null(r_b$candidates)) fail("oneab: no candidates CSV written")
if (nrow(r_b$candidates) != 2)
  fail(paste("oneab: 1a/1b co-infection must yield 2 candidates, got",
             nrow(r_b$candidates)))
gts <- sort(unique(as.character(r_b$candidates$candidate_genotype)))
if (!identical(gts, c("1")))
  fail(paste("oneab: both candidates should be genotype 1, got",
             paste(gts, collapse = ",")))
sts <- sort(unique(as.character(r_b$candidates$candidate_subtype)))
if (!identical(sts, c("1a", "1b")))
  fail(paste("oneab: candidates should be subtypes 1a and 1b, got",
             paste(sts, collapse = ",")))
ok("oneab -> 1a + 1b both selected (genotype-1 dedup does not collapse them)")

# --- Case (c) no validity filtering: 2k1b second still selected (D-05) -------
r_c <- run_select(
  "novalid",
  "NOVALID",
  refs = list("1a_DOM" = 900, "2k1b_SUB" = 400),
  covs = list("1a_DOM" = 50, "2k1b_SUB" = 50),
  n_candidates = 2
)
if (is.null(r_c$candidates)) fail("novalid: no candidates CSV written")
if (nrow(r_c$candidates) != 2)
  fail(paste("novalid: 2k1b second candidate must NOT be suppressed (D-05), got",
             nrow(r_c$candidates), "candidate rows"))
if (!("2k1b_SUB" %in% r_c$candidates$candidate_ref))
  fail("novalid: 2k1b_SUB must be selectable (is_valid_minor removed, D-05)")
ok("novalid -> 2k1b second candidate selected (no validity filtering)")

# --- Case (d) shim reconstruction: rank1->major_*, rank2->minor_* ------------
if (is.null(r_a$legacy)) fail("shim: no legacy parsefirstmapping CSV written")
if (!all(LEGACY_COLS %in% colnames(r_a$legacy)))
  fail(paste("shim: legacy 10 columns missing; got:",
             paste(colnames(r_a$legacy), collapse = ",")))
if (r_a$legacy$major_ref[1] != "1a_TOP")
  fail(paste("shim: major_ref should be rank1 1a_TOP, got", r_a$legacy$major_ref[1]))
if (r_a$legacy$minor_ref[1] != "1b_ONE")
  fail(paste("shim: minor_ref should be rank2 1b_ONE, got", r_a$legacy$minor_ref[1]))
ok("shim -> rank1->major_ref, rank2->minor_ref, all 10 legacy columns present")

# --- Case (e) single-candidate NA-fill + empty-input no_mapping default ------
r_single <- run_select(
  "single",
  "SINGLE",
  refs = list("3a_ONLY" = 700),
  covs = list("3a_ONLY" = 50),
  n_candidates = 2
)
if (is.null(r_single$candidates)) fail("single: no candidates CSV written")
if (nrow(r_single$candidates) != 1)
  fail(paste("single: one distinct subtype must yield exactly 1 candidate, got",
             nrow(r_single$candidates)))
if (is.null(r_single$legacy)) fail("single: no legacy CSV written")
if (!all(LEGACY_COLS %in% colnames(r_single$legacy)))
  fail("single: all 10 legacy columns must be present even on a single candidate")
if (!is.na(r_single$legacy$minor_ref[1]))
  fail(paste("single: minor_ref must be NA-filled on a single-candidate sample, got",
             r_single$legacy$minor_ref[1]))
ok("single -> 1 candidate, minor_* NA-filled, all 10 legacy columns present")

# Empty input (only the '*' unmapped sentinel): default gate_flag = no_mapping,
# no crash on the long-format / legacy reconstruction.
r_empty <- run_select(
  "empty",
  "EMPTY",
  refs = list(),
  covs = list(),
  n_candidates = 2
)
if (is.null(r_empty$legacy)) fail("empty: legacy CSV must still be written on no-mapping input")
if (r_empty$legacy$gate_flag[1] != "no_mapping")
  fail(paste("empty: gate_flag must default to no_mapping, got",
             r_empty$legacy$gate_flag[1]))
ok("empty -> gate_flag=no_mapping default, no crash on empty frame")

cat("\nALL PASS\n")
