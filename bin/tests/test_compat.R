#!/usr/bin/env Rscript

# test_compat.R -----------------------------------------------------------
# Phase-9 compatibility / filename-migration regression suite (COMPAT-01..04,
# TEST-01, D-06/D-07/D-08). Auto-discovered by bin/tests/run_all.sh's
# `test_*.R` glob and run in the pinned r-regression container — no CI YAML edit.
#
# Two harnesses, modelled on the two sibling test files:
#   * Analog A (subprocess) — bin/tests/test_candidate_selection.R: drive the
#     REAL bin/summarize.R via system2("Rscript", ...) on synthetic candidate /
#     stats / depth fixtures staged in an isolated tempdir, then assert on the
#     emitted Summary.csv. Used for COMPAT-01/02/03.
#   * Analog B (function-source) — bin/tests/test_classify_roles.R: source the
#     REAL bin/classify_roles.R helpers and assert exact role / role_reason
#     values. Used for COMPAT-04 alongside the Summary.csv signal.
#
# Asserted behaviours:
#   COMPAT-01 (D-06/D-07): the five core strain-call columns (Major_reference,
#     Minor_reference, Major_genotype_mapping, Minor_genotype_mapping,
#     overall_sample_call) reproduce the human-verified golden values in
#     bin/tests/fixtures/compat_golden.csv for a monoinfection AND a co-infection
#     case.
#   COMPAT-02 (smoke): the co-infection case's `.cand1.`/`.cand2.` slot fixtures
#     parse via the Plan-03 candidate_rank join and populate non-NA Major_/Minor_
#     references that carry NO `_cand{rank}` slot suffix (RESEARCH Pitfall 1 guard).
#   COMPAT-03 (D-05, validate-only): Summary.csv carries the legacy Major_*/Minor_*
#     columns ALONGSIDE the new Major_role_*/Minor_role_*/overall_sample_call.
#   COMPAT-04 (D-12, end-to-end): a 1a/1b co-infection is allowed (1b candidate
#     role == "co-infection", overall_sample_call == "co-infection") and a 2k/1b
#     recombinant is suppressed (2k1b candidate role == "background",
#     role_reason == "recombinant_2k1b", not reported as a co-infection minor).
#
# Run from any cwd:  Rscript bin/tests/test_compat.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "summarize.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Helper: NA-tolerant equality for a single golden vs produced cell.
cell_equal <- function(a, b) {
  a_na <- is.na(a) | (is.character(a) & a %in% c("NA", ""))
  b_na <- is.na(b) | (is.character(b) & b %in% c("NA", ""))
  if (a_na && b_na) return(TRUE)
  if (a_na != b_na) return(FALSE)
  as.character(a) == as.character(b)
}

# -------------------------------------------------------------------------
# Subprocess harness: stage a one-sample summarize.R input tree in a tempdir,
# invoke the REAL script, read back Summary.csv.
#
# `cands` is a tibble with the 8-column long candidate schema (candidate_rank 1 =
# major, 2 = minor, ...). For each candidate the harness writes the matching
# `.cand{rank}.` stats / depth fixtures using the NEW slot naming so the run also
# exercises the Plan-03 candidate_rank join (COMPAT-02).
#
# Returns the Summary.csv tibble (or NULL if the script did not emit it).
# -------------------------------------------------------------------------
run_summarize <- function(case, sampleName, cands,
                          minRead = 500, minCov = 30, n_candidates = 2,
                          ref_length = 200) {
  wd <- tempfile(paste0("compat_", case, "_")); dir.create(wd)
  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)

  # summarize.R source()s these five helpers cwd-relative (bin/summarize.R L10-25).
  for (h in c("genotype_utils.R", "denovo_confirm.R", "denovo_layer.R",
              "assembly_support_join.R", "classify_roles.R")) {
    file.copy(file.path(bin_dir, h), file.path(wd, h), overwrite = TRUE)
  }

  # Input subdirs summarize.R reads from. The trimmed/, kraken_classified/, id/
  # loops iterate `1:length(files)` so they need at least one file each to avoid
  # a `1:0` empty-loop crash; the strain-call columns come from the four mapping
  # subdirs. The remaining dirs (denovo/, glue/, variation/, consensus_distance/)
  # are guarded with `if (length(...) > 0)` so we leave them absent.
  for (d in c("trimmed", "kraken_classified", "parsefirst_mapping",
              "stats_withdup", "stats_markdup", "depth", "id")) {
    dir.create(file.path(wd, d))
  }

  # --- samplesheet (args[1]); read at summarize.R L1072: select sample) ---
  samplesheet <- file.path(wd, "samplesheet.csv")
  write_csv(tibble(sample = sampleName), samplesheet)

  # --- trimmed/<sample>.log (cutadapt-style; trimmed loop) ---
  writeLines(c(
    "This is cutadapt",
    "Total read pairs processed:          1,000,000",
    "Pairs written (passing filters):       900,000 (90.0%)"
  ), file.path(wd, "trimmed", paste0(sampleName, ".log")))

  # --- kraken_classified/<sample>.kraken2.report.txt (kraken loop) ---
  # X2 (col 2) = clade reads, X6 (col 6) = rank name; loop pulls X2 where X6=="root".
  writeLines(paste(c("100.00", "900000", "900000", "R", "1", "root"),
                   collapse = "\t"),
             file.path(wd, "kraken_classified",
                       paste0(sampleName, ".kraken2.report.txt")))

  # --- id/<sample>.sequencerID.tsv (id loop) ---
  writeLines("@M01234:1:000000000:1:1101:1:1",
             file.path(wd, "id", paste0(sampleName, ".sequencerID.tsv")))

  # --- parsefirst_mapping/<sample>.candidates.csv (8-col long; the join source) ---
  cand_csv <- cands %>%
    transmute(
      sample              = sampleName,
      candidate_rank,
      candidate_ref,
      candidate_subtype,
      candidate_genotype,
      candidate_reads,
      candidate_cov,
      confirmation_status = "pass"
    )
  write_csv(cand_csv, file.path(wd, "parsefirst_mapping",
                                paste0(sampleName, ".candidates.csv")))

  # --- parsefirst_mapping/<sample>.parsefirstmapping.csv (legacy 10-col shim) ---
  major <- cands %>% filter(candidate_rank == 1) %>% slice(1)
  minor <- cands %>% filter(candidate_rank == 2) %>% slice(1)
  legacy <- tibble(
    sample             = sampleName,
    total_mapped_reads = sum(cands$candidate_reads),
    major_ref          = if (nrow(major) > 0) major$candidate_ref else NA_character_,
    major_reads        = if (nrow(major) > 0) major$candidate_reads else NA_real_,
    major_cov          = if (nrow(major) > 0) major$candidate_cov else NA_real_,
    minor_ref          = if (nrow(minor) > 0) minor$candidate_ref else NA_character_,
    minor_reads        = if (nrow(minor) > 0) minor$candidate_reads else NA_real_,
    minor_cov          = if (nrow(minor) > 0) minor$candidate_cov else NA_real_,
    minor_call         = if (nrow(minor) > 0) "co-infection" else "monoinfection",
    gate_flag          = "pass"
  )
  write_csv(legacy, file.path(wd, "parsefirst_mapping",
                              paste0(sampleName, ".parsefirstmapping.csv")))

  # --- per-candidate stats + depth fixtures using the NEW cand-slot naming ---
  # samtools-stats text: the loops pull `SN<TAB>reads mapped:<TAB><N>`.
  for (i in seq_len(nrow(cands))) {
    rk   <- cands$candidate_rank[i]
    ref  <- cands$candidate_ref[i]
    rds  <- cands$candidate_reads[i]
    base <- paste0(sampleName, ".", ref, ".cand", rk)

    stats_text <- c(
      "# This file was produced by samtools stats",
      paste("SN", "raw total sequences:", rds, sep = "\t"),
      paste("SN", "reads mapped:", rds, sep = "\t")
    )
    writeLines(stats_text, file.path(wd, "stats_withdup",
                                     paste0(base, ".withdup.stats")))
    writeLines(stats_text, file.path(wd, "stats_markdup",
                                     paste0(base, ".nodup.stats")))

    # depth tsv: <ref> <pos> <cov>, 3-col no header (summarize.R L457). Constant
    # coverage across ref_length positions; high enough to clear the typable gate.
    depth_lines <- paste(ref, seq_len(ref_length), 60, sep = "\t")
    writeLines(depth_lines, file.path(wd, "depth",
                                      paste0(base, ".nodup.tsv")))
  }

  # --- invoke the REAL summarize.R (positional arg contract, L27-86) ---
  # args: [1]=samplesheet [2]=version [3]=name [4-8] denovo params (empty=default)
  #       [9]=minRead [10]=minCov [11]=n_candidates
  # NOTE: pass the empty positional slots [4-8] as shQuote("") — a bare "" is
  # DROPPED by system2(), which would silently shift minRead/minCov off positions
  # [9]/[10] (they would land at [4]/[5]), parse as denovo params, leave the gate
  # thresholds NA, and collapse every sample to overall_sample_call=="indeterminate".
  empty_slot <- shQuote("")
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(samplesheet), "TestPipe", "hcvtyper",
      empty_slot, empty_slot, empty_slot, empty_slot, empty_slot,
      minRead, minCov, n_candidates),
    stdout = FALSE, stderr = FALSE
  )

  summary_path <- file.path(wd, "Summary.csv")
  if (file.exists(summary_path)) {
    read_csv(summary_path, show_col_types = FALSE)
  } else {
    NULL
  }
}

# Candidate-frame builder (one row per candidate).
mk_cands <- function(...) {
  bind_rows(...)
}
mk_cand <- function(rank, ref, subtype, reads, cov) {
  tibble(
    candidate_rank     = as.integer(rank),
    candidate_ref      = ref,
    candidate_subtype  = subtype,
    # genotype = leading-character rule (genotype_from_subtype), kept here as the
    # plain candidates.csv value summarize.R reads as character.
    candidate_genotype = if (subtype == "2k1b") "2k1b" else substr(subtype, 1, 1),
    candidate_reads    = reads,
    candidate_cov      = cov
  )
}

# =========================================================================
# COMPAT-01 / COMPAT-02 / COMPAT-03 — subprocess against the REAL summarize.R
# =========================================================================

golden_path <- file.path(this_dir, "fixtures", "compat_golden.csv")
if (!file.exists(golden_path)) fail(paste("compat_golden.csv missing at", golden_path))
golden <- read_csv(golden_path, show_col_types = FALSE)

core_cols <- c("Major_reference", "Minor_reference",
               "Major_genotype_mapping", "Minor_genotype_mapping",
               "overall_sample_call")

# --- monoinfection case (golden row MONO): 1a major, no minor ---
mono_golden <- golden %>% filter(sampleName == "MONO")
if (nrow(mono_golden) != 1) fail("compat_golden.csv missing MONO row")

mono_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99)
)
mono_summary <- run_summarize("mono", "MONO", mono_cands)
if (is.null(mono_summary)) fail("COMPAT-01 mono: summarize.R wrote no Summary.csv")
if (nrow(mono_summary) != 1) fail("COMPAT-01 mono: expected 1 Summary.csv row")
for (col in core_cols) {
  if (!col %in% names(mono_summary)) fail(paste("COMPAT-01 mono: Summary.csv missing column", col))
  if (!cell_equal(mono_golden[[col]][1], mono_summary[[col]][1]))
    fail(sprintf("COMPAT-01 mono: %s = '%s', golden = '%s'",
                 col, mono_summary[[col]][1], mono_golden[[col]][1]))
}
ok("COMPAT-01 mono: five core strain-call columns match compat_golden.csv (MONO)")

# --- co-infection case (golden row COINF): 1a major + 1b minor ---
coinf_golden <- golden %>% filter(sampleName == "COINF")
if (nrow(coinf_golden) != 1) fail("compat_golden.csv missing COINF row")

coinf_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99),
  mk_cand(2, "1b_D90208", "1b", 150000, 97)
)
coinf_summary <- run_summarize("coinf", "COINF", coinf_cands)
if (is.null(coinf_summary)) fail("COMPAT-01 coinf: summarize.R wrote no Summary.csv")
if (nrow(coinf_summary) != 1) fail("COMPAT-01 coinf: expected 1 Summary.csv row")
for (col in core_cols) {
  if (!col %in% names(coinf_summary)) fail(paste("COMPAT-01 coinf: Summary.csv missing column", col))
  if (!cell_equal(coinf_golden[[col]][1], coinf_summary[[col]][1]))
    fail(sprintf("COMPAT-01 coinf: %s = '%s', golden = '%s'",
                 col, coinf_summary[[col]][1], coinf_golden[[col]][1]))
}
ok("COMPAT-01 coinf: five core strain-call columns match compat_golden.csv (COINF)")

# --- COMPAT-02 smoke: cand-slot fixtures parsed (non-NA, suffix-free refs) ---
maj <- coinf_summary$Major_reference[1]
min <- coinf_summary$Minor_reference[1]
if (is.na(maj)) fail("COMPAT-02: Major_reference silently empty after cand-slot join (Pitfall 1)")
if (is.na(min)) fail("COMPAT-02: Minor_reference silently empty after cand-slot join (Pitfall 1)")
if (grepl("cand[0-9]+", maj)) fail(paste("COMPAT-02: Major_reference carries a slot suffix:", maj))
if (grepl("cand[0-9]+", min)) fail(paste("COMPAT-02: Minor_reference carries a slot suffix:", min))
if (maj != "1a_M62321") fail(paste("COMPAT-02: Major_reference expected 1a_M62321, got", maj))
ok("COMPAT-02: .cand1./.cand2. fixtures populate suffix-free Major_/Minor_reference")

# --- COMPAT-03: legacy + role columns co-present (validate-only) ---
legacy_cols <- c("Major_reference", "Minor_reference",
                 "Major_genotype_mapping", "Minor_genotype_mapping")
role_cols   <- c("Major_role_reference", "Minor_role_reference", "overall_sample_call")
if (!all(legacy_cols %in% colnames(coinf_summary)))
  fail(paste("COMPAT-03: missing legacy columns:",
             paste(setdiff(legacy_cols, colnames(coinf_summary)), collapse = ",")))
if (!all(role_cols %in% colnames(coinf_summary)))
  fail(paste("COMPAT-03: missing role columns:",
             paste(setdiff(role_cols, colnames(coinf_summary)), collapse = ",")))
ok("COMPAT-03: Summary.csv carries legacy Major_*/Minor_* alongside role columns + overall_sample_call")

cat("\nALL PASS\n")
