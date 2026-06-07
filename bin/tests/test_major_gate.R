#!/usr/bin/env Rscript

# test_major_gate.R -------------------------------------------------------
# Subprocess contract test (D-03) for the REAL major-gate entrypoint
# bin/summarize_mapping_to_all_references.R. We do NOT modify Phase-1 code:
# we invoke the script unchanged via `system2("Rscript", ...)` on static
# committed idxstats / depth / refs.fa fixtures and assert on the columns it
# writes to <sampleName>.parsefirstmapping.csv — gate_flag + minor_call.
#
# Three gate states are exercised:
#   (a) gate_pass        -> gate_flag == "ok",                    minor_call == "yes"
#   (b) gate_majorfail   -> gate_flag == "major_below_threshold", minor_call == "no"
#       (ERR1810469-class: the read-count winner itself is sub-threshold; ALL
#        refs are below threshold so the script cannot re-select a passing major)
#   (c) gate_nomapping   -> gate_flag == "no_mapping",            minor_call == "no"
#       The script writes the correct CSV and THEN crashes non-zero on the
#       write.fasta(fasta[major_ref]) line because major_ref is unset. This is
#       a documented Phase-1 LATENT BUG that D-03 forbids fixing here, so the
#       test tolerates the non-zero exit and asserts on the CSV instead.
#
# The script sources genotype_utils.R via a cwd-relative path and writes its
# output CSV/FASTA into the cwd, so each case runs in its own tempdir with
# genotype_utils.R copied in and absolute fixture paths passed as args. This
# isolates the three subprocess outputs and keeps Phase-1 code untouched.
#
# Run from any cwd:
#   Rscript bin/tests/test_major_gate.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd
# (scaffold copied from test_summarize_denovo.R:28-41).
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
fixt_dir <- file.path(this_dir, "fixtures")
script   <- file.path(bin_dir, "summarize_mapping_to_all_references.R")
geno_src <- file.path(bin_dir, "genotype_utils.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Invoke the REAL script in an isolated per-case tempdir. The script sources
# "genotype_utils.R" cwd-relative and writes "<sample>.parsefirstmapping.csv"
# cwd-relative, so we run with the subprocess cwd set to a fresh tempdir that
# has genotype_utils.R copied in, and pass absolute fixture paths as args.
# Returns list(exit, csv) where csv is the read-back data frame (or NULL).
run_gate <- function(case, sampleName, minRead, minCov) {
  idx   <- normalizePath(file.path(fixt_dir, paste0(case, "_idxstats.tsv")))
  depth <- normalizePath(file.path(fixt_dir, paste0(case, "_depth.tsv")))
  refs  <- normalizePath(file.path(fixt_dir, paste0(case, "_refs.fa")))

  wd <- tempfile(paste0("gate_", case, "_")); dir.create(wd)
  file.copy(geno_src, file.path(wd, "genotype_utils.R"), overwrite = TRUE)

  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(idx), shQuote(depth), shQuote(sampleName),
      shQuote(refs), minRead, minCov),
    stdout = FALSE, stderr = FALSE
  )

  csv_path <- file.path(wd, paste0(sampleName, ".parsefirstmapping.csv"))
  csv <- if (file.exists(csv_path)) read_csv(csv_path, show_col_types = FALSE) else NULL
  list(exit = exit, csv = csv)
}

# --- Case (a) gate_pass: both refs above threshold, genotypes 1 vs 3 -------
r_pass <- run_gate("gate_pass", "GP", 500, 30)
if (is.null(r_pass$csv)) fail("gate_pass: no CSV written")
if (r_pass$csv$gate_flag[1] != "ok") fail(paste("gate_pass gate_flag:", r_pass$csv$gate_flag[1]))
if (r_pass$csv$minor_call[1] != "yes") fail(paste("gate_pass minor_call:", r_pass$csv$minor_call[1]))
ok("gate_pass -> gate_flag=ok, minor_call=yes")

# --- Case (b) gate_majorfail (ERR1810469-class): all refs sub-threshold ----
r_fail <- run_gate("gate_majorfail", "GF", 500, 30)
if (is.null(r_fail$csv)) fail("gate_majorfail: no CSV written")
if (r_fail$csv$gate_flag[1] != "major_below_threshold")
  fail(paste("gate_majorfail gate_flag:", r_fail$csv$gate_flag[1]))
if (r_fail$csv$minor_call[1] != "no") fail(paste("gate_majorfail minor_call:", r_fail$csv$minor_call[1]))
ok("gate_majorfail -> gate_flag=major_below_threshold, minor_call=no")

# --- Case (c) gate_nomapping: only the '*' unmapped idxstats row -----------
# KNOWN PHASE-1 LATENT BUG (D-03 forbids fixing it here): the script writes the
# correct CSV (no_mapping / no) then exits NON-ZERO on write.fasta(fasta[major_ref])
# because major_ref is never assigned on the no-mapping branch. We therefore do
# NOT require exit 0; we assert on the written CSV. stderr noise (e.g. the
# cur_data() deprecation warning) is non-fatal and ignored — we assert the CSV.
r_nomap <- run_gate("gate_nomapping", "GN", 500, 30)
if (is.null(r_nomap$csv)) fail("gate_nomapping: CSV must still be written before the latent FASTA crash")
if (r_nomap$csv$gate_flag[1] != "no_mapping")
  fail(paste("gate_nomapping gate_flag:", r_nomap$csv$gate_flag[1]))
if (r_nomap$csv$minor_call[1] != "no") fail(paste("gate_nomapping minor_call:", r_nomap$csv$minor_call[1]))
ok("gate_nomapping -> CSV gate_flag=no_mapping, minor_call=no (non-zero exit tolerated, latent Phase-1 bug)")

cat("\nALL PASS\n")
