#!/usr/bin/env Rscript

# test_rescue_evaluation.R -----------------------------------------------
# Subprocess-contract test (Phase 10, denovo-subtype-rescue) for the REAL
# de-novo subtype-rescue evaluator in bin/rescue_evaluation.R. Modelled on
# bin/tests/test_assembly_support.R: we invoke the script via
# system2("Rscript", ...) on small synthetic candidates / blastparse /
# assembly_support / references fixtures built inline in a tempdir, and assert
# on what it writes to <prefix>.candidates.csv (+ the rescued *_cand{rank}.fa).
#
# This test is the executable specification (TDD RED) for the rescue contract.
# It encodes the six rescue behaviours locked in 10-CONTEXT / 10-PATTERNS
# (D-02..D-10):
#
#   (1) rescue-fires : candidate_subtype != denovo top-hit subtype AND all four
#       floors pass -> rescued_from == original candidate ref, rescue_trigger
#       non-NA, confirmation_status forced to "pass", and a FASTA whose basename
#       contains "_cand{rank}." is written for the rescued reference.
#   (2) floor-fail   : any single floor fails -> rescued_from is NA and
#       confirmation_status is unchanged (no rescue).
#   (3) 2k1b-rule    : a 2k1b candidate (or 2k1b denovo hit) plus a genotype-2
#       contig meeting the floors triggers rescue to the genotype-2 reference
#       (D-03 special recombinant rule).
#   (4) 1a1b-floor   : the 1a/1b boundary uses the higher rescue_1a1b_length
#       floor — a 1a-vs-1b mismatch with a 3500bp contig does NOT rescue, a
#       5200bp one does (D-04).
#   (5) skip-assembly: empty blastparse/support tables (zero data-row input
#       CSVs) -> candidates pass through unchanged with rescued_from /
#       rescue_trigger NA-filled and the columns still present (D-10 DoS guard).
#   (6) forced-pass  : the rescued FASTA basename matches "_cand{rank}." and
#       confirmation_status is forced to "pass" even when the original candidate
#       row was "below_threshold" (D-05/D-06).
#
# The thresholds are passed positionally (PARSEFIRSTMAPPING convention), matching
# the D-02 defaults: length 3000, pident 85, aln_length 3000, kmer_cov 2,
# 1a1b_length 5000.
#
# RED until Plan 02 lands bin/rescue_evaluation.R (the script does not yet exist,
# so every subtest's system2() call returns non-zero / writes nothing). Do NOT
# create bin/rescue_evaluation.R in this plan — the failure is the intended state.
#
# Run from any cwd:
#   Rscript bin/tests/test_rescue_evaluation.R
# Exits 0 and prints "ALL PASS" when all assertions hold (i.e. once Plan 02 lands
# the script); exits non-zero (RED) until then.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "rescue_evaluation.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Guard: distinguish "script missing" (intended RED) from a syntax error in the
# test itself. When the script is absent every subtest fails on the contract
# assertions below; we surface a clear RED banner so the failure is unambiguous.
if (!file.exists(script)) {
  cat("RED: bin/rescue_evaluation.R does not exist yet — this is the intended",
      "TDD RED state for Phase 10 Plan 01.\n")
  cat("FAIL: rescue_evaluation.R is not yet implemented (expected in Plan 02)\n")
  quit(status = 1)
}

# The 8-column long candidates schema (modules/local/parsefirstmapping/main.nf:77)
# that the rescue script reads, PLUS the two audit columns (D-08) it must add.
CAND_IN_COLS <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
                  "candidate_genotype", "candidate_reads", "candidate_cov",
                  "confirmation_status")
RESCUE_COLS  <- c("rescued_from", "rescue_trigger")

# assembly_support.csv schema (bin/blast_parse.R §4b) — the four-floor source.
SUPPORT_COLS <- c("sample", "subtype", "best_contig_length", "best_contig_pident",
                  "best_contig_aln_length", "best_contig_kmer_cov")

# D-02 default thresholds, passed positionally to the script.
TH_LENGTH    <- 3000
TH_PIDENT    <- 85
TH_ALN       <- 3000
TH_KMER      <- 2
TH_1A1B_LEN  <- 5000

# -------------------------------------------------------------------------
# Subprocess harness: stage candidates / blastparse / support / references
# fixtures in a tempdir, invoke the REAL rescue_evaluation.R, read back the
# rewritten <prefix>.candidates.csv and list any written FASTA files.
#
#   cands   : tibble with the 8-column long candidate schema.
#   support : tibble with the 6-column assembly_support schema (one row/subtype).
#   bparse  : tibble with the de-novo blastparse top-hit schema
#             (sample, major_ref, major_contig_length, minor_ref,
#              minor_contig_length).
#   refs    : named character vector subtype-tagged accession -> dummy sequence;
#             FASTA headers use the {subtype}_{accession} convention.
#
# Returns list(exit, cands_out, fastas, path).
# -------------------------------------------------------------------------
run_rescue <- function(case, prefix, cands, support, bparse, refs) {
  wd <- tempfile(paste0("rescue_", case, "_")); dir.create(wd)

  cand_path <- file.path(wd, paste0(prefix, ".candidates.in.csv"))
  write_csv(cands, cand_path)

  support_path <- file.path(wd, paste0(prefix, ".assembly_support.csv"))
  write_csv(support, support_path)

  bparse_path <- file.path(wd, paste0(prefix, ".blastparse.csv"))
  write_csv(bparse, bparse_path)

  # references FASTA: {subtype}_{accession} headers, dummy sequence per record.
  ref_lines <- unlist(lapply(names(refs), function(h) c(paste0(">", h), refs[[h]])))
  refs_path <- file.path(wd, paste0(prefix, "_references.fa"))
  writeLines(ref_lines, refs_path)

  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(prefix),
      shQuote(cand_path), shQuote(bparse_path), shQuote(support_path),
      shQuote(refs_path),
      TH_LENGTH, TH_PIDENT, TH_ALN, TH_KMER, TH_1A1B_LEN),
    stdout = FALSE, stderr = FALSE
  )

  out_path <- file.path(wd, paste0(prefix, ".candidates.csv"))
  cands_out <- if (file.exists(out_path)) {
    read_csv(out_path, show_col_types = FALSE)
  } else {
    NULL
  }
  fastas <- list.files(wd, pattern = "\\.fa$", full.names = FALSE)
  list(exit = exit, cands_out = cands_out, fastas = fastas, path = out_path)
}

# Candidate-row builder (8-column long schema).
mk_cand <- function(rank, ref, subtype, genotype, reads, cov, status) {
  tibble(
    candidate_rank      = as.integer(rank),
    candidate_ref       = ref,
    candidate_subtype   = subtype,
    candidate_genotype  = genotype,
    candidate_reads     = reads,
    candidate_cov       = cov,
    confirmation_status = status
  )
}
mk_cands <- function(sample, ...) {
  bind_rows(...) %>% mutate(sample = sample) %>% select(all_of(CAND_IN_COLS))
}

mk_support <- function(sample, ...) {
  bind_rows(...) %>% mutate(sample = sample) %>% select(all_of(SUPPORT_COLS))
}
mk_support_row <- function(subtype, length, pident, aln, kmer) {
  tibble(subtype = subtype, best_contig_length = length,
         best_contig_pident = pident, best_contig_aln_length = aln,
         best_contig_kmer_cov = kmer)
}

# A tiny multi-record reference set covering the subtypes the subtests reach for.
REFS <- list(
  "2b_AY232748" = strrep("A", 60),
  "2k1b_AF177036" = strrep("C", 60),
  "1a_M62321"   = strrep("G", 60),
  "1b_D90208"   = strrep("T", 60),
  "2a_AB047639" = strrep("A", 60)
)

# Empty/typed blastparse + support frames for the skip-assembly case.
empty_support <- tibble(
  sample = character(0), subtype = character(0),
  best_contig_length = double(0), best_contig_pident = double(0),
  best_contig_aln_length = double(0), best_contig_kmer_cov = double(0)
)
empty_bparse <- tibble(
  sample = character(0), major_ref = character(0),
  major_contig_length = double(0), minor_ref = character(0),
  minor_contig_length = double(0)
)

assert_schema <- function(case, cands_out) {
  if (is.null(cands_out))
    fail(paste(case, ": no <prefix>.candidates.csv written"))
  missing <- setdiff(c(CAND_IN_COLS, RESCUE_COLS), colnames(cands_out))
  if (length(missing) > 0)
    fail(paste(case, ": output candidates.csv missing columns:",
               paste(missing, collapse = ",")))
}

# =========================================================================
# Subtest 1: rescue-fires (mismatch + all four floors pass)
# =========================================================================
# Candidate rank-1 is mapped to a 2k1b/1b-flavoured ref but the de-novo top hit
# is subtype 2b with a contig clearing all four floors -> rescue to the 2b ref.
r1 <- run_rescue(
  "fires", "FIRES",
  cands = mk_cands("FIRES",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("FIRES",
    mk_support_row("2b", 4120, 96.2, 3840, 4.1)),
  bparse = tibble(sample = "FIRES", major_ref = "2b_AY232748",
                  major_contig_length = 4120, minor_ref = NA_character_,
                  minor_contig_length = NA_real_),
  refs = REFS
)
assert_schema("fires", r1$cands_out)
row1 <- r1$cands_out %>% filter(candidate_rank == 1)
if (is.na(row1$rescued_from[1]))
  fail("fires: rescued_from must record the original candidate ref when rescue fires")
if (row1$rescued_from[1] != "1a_M62321")
  fail(paste("fires: rescued_from must equal the original candidate ref 1a_M62321, got",
             row1$rescued_from[1]))
if (is.na(row1$rescue_trigger[1]))
  fail("fires: rescue_trigger must be a non-NA description when rescue fires")
if (row1$confirmation_status[1] != "pass")
  fail(paste("fires: confirmation_status must be forced to 'pass' on rescue, got",
             row1$confirmation_status[1]))
if (!any(grepl("_cand1\\.", r1$fastas)))
  fail(paste("fires: a rescued FASTA with a '_cand1.' basename must be written; got:",
             paste(r1$fastas, collapse = ",")))
ok("rescue-fires -> rescued_from set, rescue_trigger non-NA, status forced pass, _cand1. FASTA written")

# =========================================================================
# Subtest 2: floor-fail (one floor fails -> no rescue)
# =========================================================================
# Same mismatch as subtest 1 but the de-novo contig's pident (80) is below the
# 85 floor -> no rescue: rescued_from NA, confirmation_status unchanged.
r2 <- run_rescue(
  "floorfail", "FLOORFAIL",
  cands = mk_cands("FLOORFAIL",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("FLOORFAIL",
    mk_support_row("2b", 4120, 80.0, 3840, 4.1)),   # pident 80 < 85 floor
  bparse = tibble(sample = "FLOORFAIL", major_ref = "2b_AY232748",
                  major_contig_length = 4120, minor_ref = NA_character_,
                  minor_contig_length = NA_real_),
  refs = REFS
)
assert_schema("floorfail", r2$cands_out)
f2 <- r2$cands_out %>% filter(candidate_rank == 1)
if (!is.na(f2$rescued_from[1]))
  fail(paste("floorfail: rescued_from must be NA when a floor fails, got",
             f2$rescued_from[1]))
if (f2$confirmation_status[1] != "pass")
  fail("floorfail: confirmation_status must be unchanged ('pass') when no rescue fires")
ok("floor-fail -> rescued_from NA, confirmation_status unchanged (no rescue on sub-floor pident)")

# =========================================================================
# Subtest 3: 2k1b special rule (D-03)
# =========================================================================
# A 2k1b candidate plus a genotype-2 contig meeting the floors triggers rescue to
# the genotype-2 reference (the de-novo 2a top hit).
r3 <- run_rescue(
  "k2b", "K2B",
  cands = mk_cands("K2B",
    mk_cand(1, "2k1b_AF177036", "2k1b", "2k1b", 60000, 92, "pass")),
  support = mk_support("K2B",
    mk_support_row("2a", 4500, 97.0, 4200, 6.0)),
  bparse = tibble(sample = "K2B", major_ref = "2a_AB047639",
                  major_contig_length = 4500, minor_ref = NA_character_,
                  minor_contig_length = NA_real_),
  refs = REFS
)
assert_schema("k2b", r3$cands_out)
k3 <- r3$cands_out %>% filter(candidate_rank == 1)
if (is.na(k3$rescued_from[1]))
  fail("2k1b-rule: a 2k1b candidate over a genotype-2 contig meeting floors must rescue")
if (k3$rescued_from[1] != "2k1b_AF177036")
  fail(paste("2k1b-rule: rescued_from must record the original 2k1b ref, got",
             k3$rescued_from[1]))
if (is.na(k3$rescue_trigger[1]))
  fail("2k1b-rule: rescue_trigger must be non-NA for the 2k1b special-rule rescue")
ok("2k1b-rule -> genotype-2 contig over a 2k1b candidate triggers rescue (D-03)")

# =========================================================================
# Subtest 4: 1a/1b boundary uses the higher length floor (D-04)
# =========================================================================
# A 1a-vs-1b mismatch with a 3500bp contig does NOT rescue (below the 5000bp
# 1a1b floor); a 5200bp one DOES.
r4a <- run_rescue(
  "ab_short", "ABSHORT",
  cands = mk_cands("ABSHORT",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("ABSHORT",
    mk_support_row("1b", 3500, 99.0, 3400, 8.0)),    # 3500 < 5000 (1a1b floor)
  bparse = tibble(sample = "ABSHORT", major_ref = "1b_D90208",
                  major_contig_length = 3500, minor_ref = NA_character_,
                  minor_contig_length = NA_real_),
  refs = REFS
)
assert_schema("ab_short", r4a$cands_out)
a4 <- r4a$cands_out %>% filter(candidate_rank == 1)
if (!is.na(a4$rescued_from[1]))
  fail(paste("1a1b-floor: a 3500bp 1a/1b contig must NOT rescue (below the 5000bp floor), got rescued_from",
             a4$rescued_from[1]))

r4b <- run_rescue(
  "ab_long", "ABLONG",
  cands = mk_cands("ABLONG",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("ABLONG",
    mk_support_row("1b", 5200, 99.0, 5100, 8.0)),    # 5200 >= 5000 floor
  bparse = tibble(sample = "ABLONG", major_ref = "1b_D90208",
                  major_contig_length = 5200, minor_ref = NA_character_,
                  minor_contig_length = NA_real_),
  refs = REFS
)
assert_schema("ab_long", r4b$cands_out)
b4 <- r4b$cands_out %>% filter(candidate_rank == 1)
if (is.na(b4$rescued_from[1]))
  fail("1a1b-floor: a 5200bp 1a/1b contig (>= the 5000bp floor) must rescue")
ok("1a1b-floor -> 3500bp no rescue, 5200bp rescues (higher rescue_1a1b_length floor, D-04)")

# =========================================================================
# Subtest 5: skip-assembly empty tables -> pass-through with NA columns (D-10)
# =========================================================================
# Zero-row blastparse + support inputs: candidates pass through unchanged, the two
# rescue columns are present and NA-filled, the script does NOT abort.
r5 <- run_rescue(
  "skip", "SKIP",
  cands = mk_cands("SKIP",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass"),
    mk_cand(2, "1b_D90208", "1b", "1", 40000, 70, "below_threshold")),
  support = empty_support,
  bparse = empty_bparse,
  refs = REFS
)
if (r5$exit != 0)
  fail(paste("skip-assembly: script must exit 0 on empty blastparse/support input, got exit",
             r5$exit))
assert_schema("skip", r5$cands_out)
if (nrow(r5$cands_out) != 2)
  fail(paste("skip-assembly: both candidates must pass through, got",
             nrow(r5$cands_out), "rows"))
if (any(!is.na(r5$cands_out$rescued_from)))
  fail("skip-assembly: rescued_from must be NA for every row when there is no assembly evidence")
if (any(!is.na(r5$cands_out$rescue_trigger)))
  fail("skip-assembly: rescue_trigger must be NA for every row when there is no assembly evidence")
ok("skip-assembly -> candidates pass through unchanged, rescue columns present and NA (D-10 DoS guard)")

# =========================================================================
# Subtest 6: rescued FASTA naming + forced pass on a below_threshold row (D-05/D-06)
# =========================================================================
# A rank-2 candidate that was "below_threshold" gets rescued: the written FASTA
# basename matches "_cand2." and confirmation_status is forced to "pass".
r6 <- run_rescue(
  "forced", "FORCED",
  cands = mk_cands("FORCED",
    mk_cand(1, "2a_AB047639", "2a", "2", 90000, 96, "pass"),
    mk_cand(2, "1a_M62321",   "1a", "1", 20000, 40, "below_threshold")),
  support = mk_support("FORCED",
    mk_support_row("2b", 4120, 96.2, 3840, 4.1)),
  bparse = tibble(sample = "FORCED", major_ref = "2a_AB047639",
                  major_contig_length = 4500, minor_ref = "2b_AY232748",
                  minor_contig_length = 4120),
  refs = REFS
)
assert_schema("forced", r6$cands_out)
f6 <- r6$cands_out %>% filter(candidate_rank == 2)
if (is.na(f6$rescued_from[1]))
  fail("forced-pass: the below_threshold rank-2 candidate must be rescued (rescued_from set)")
if (f6$confirmation_status[1] != "pass")
  fail(paste("forced-pass: confirmation_status must be forced to 'pass' on a rescued below_threshold row, got",
             f6$confirmation_status[1]))
if (!any(grepl("_cand2\\.", r6$fastas)))
  fail(paste("forced-pass: the rescued FASTA basename must match '_cand2.'; got:",
             paste(r6$fastas, collapse = ",")))
ok("forced-pass -> below_threshold row rescued, status forced pass, _cand2. FASTA written (D-05/D-06)")

cat("\nALL PASS\n")
