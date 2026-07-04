#!/usr/bin/env Rscript

# test_evidence_report.R --------------------------------------------------
# Phase 13 Plan 01 Task 2 (EVID-05). Unit tests for the pure
# build_evidence_summary() helper in bin/classify_roles.R: every reachable
# (role, role_reason, evidence_state) tuple from classify_one_sample() must map
# to ONE non-empty, contig-language sentence that names the candidate
# ("candidate <rank> (<subtype>_<ref>)"), states the evidence_state, and states
# contig corroboration/conflict. The weak-but-matching case must state the
# actual measured contig length/identity numbers (D-02 pattern 3 / D-07), never
# editorial words. "own assembly"/"own support" phrasing is forbidden (D-01).
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Call build_evidence_summary() for one synthetic candidate. Contig metrics
# default to NA (no contig); pass them for the present-contig cases.
summ <- function(role, role_reason, evidence_state, rank = 2L,
                 ref = "2c_JX227949", subtype = "2c", assembly_exists = TRUE,
                 len = NA_real_, pid = NA_real_, kmer = NA_real_,
                 concordance_status = NA_character_) {
  build_evidence_summary(
    role = role, role_reason = role_reason, evidence_state = evidence_state,
    concordance_status = concordance_status,
    candidate_rank = rank, candidate_ref = ref, candidate_subtype = subtype,
    assembly_exists = assembly_exists,
    best_contig_length = len, best_contig_pident = pid, best_contig_kmer_cov = kmer
  )
}

# The full reachable role_reason vocabulary from classify_one_sample() L606-734,
# each paired with a plausible (role, evidence_state) and contig-metric context.
tuples <- list(
  list(role = "dominant",     rr = "dominant",                         st = "confirmed", len = 9189, pid = 100,    kmer = 40),
  list(role = "co-infection", rr = "corroborated",                     st = "confirmed", len = 9479, pid = 89.0,   kmer = 514.2),
  list(role = "co-infection", rr = "corroborated",                     st = "probable",  len = 3000, pid = 80.0,   kmer = 5),
  list(role = "background",   rr = "refuted_denovo",                   st = "refuted",   len = 600,  pid = 85.0,   kmer = 3),
  list(role = "background",   rr = "discordant_identity",              st = "refuted",   len = 700,  pid = 84.0,   kmer = 4),
  list(role = "background",   rr = "no_own_assembly",                  st = "weak",      len = NA,   pid = NA,     kmer = NA, assembly_exists = FALSE),
  list(role = "background",   rr = "weak_own_assembly_below_floor",    st = "weak",      len = 200,  pid = 50.0,   kmer = 0.5),
  list(role = "background",   rr = "same_genotype_as_dominant",        st = "confirmed", len = 9000, pid = 94.0,   kmer = 30),
  list(role = "background",   rr = "recombinant_2k1b",                 st = "confirmed", len = 9000, pid = 95.0,   kmer = 30),
  list(role = "indeterminate", rr = "indeterminate_dominance_conflict", st = "confirmed", len = 8000, pid = 92.0,  kmer = 20)
)

for (t in tuples) {
  ae <- if (!is.null(t$assembly_exists)) t$assembly_exists else TRUE
  s <- summ(t$role, t$rr, t$st, len = t$len, pid = t$pid, kmer = t$kmer, assembly_exists = ae)
  if (!is.character(s) || length(s) != 1 || is.na(s) || !nzchar(s))
    fail(sprintf("build_evidence_summary(%s) must return one non-empty string", t$rr))
  # Names the candidate: rank token + ref/subtype token (D-05).
  if (!grepl("candidate 2", s, fixed = TRUE))
    fail(sprintf("[%s] summary must name the candidate rank, got: %s", t$rr, s))
  if (!grepl("2c_JX227949", s, fixed = TRUE))
    fail(sprintf("[%s] summary must carry the ref/subtype token, got: %s", t$rr, s))
  # States the evidence_state so a reader can reconstruct it (success criterion 2).
  if (!grepl(t$st, s, fixed = TRUE))
    fail(sprintf("[%s] summary must state the evidence_state '%s', got: %s", t$rr, t$st, s))
  # D-01: contig language only — never "own assembly"/"own support".
  if (grepl("own assembly", s, fixed = TRUE) || grepl("own support", s, fixed = TRUE))
    fail(sprintf("[%s] summary must NOT use 'own assembly'/'own support' (D-01), got: %s", t$rr, s))
}
ok("EvidenceReport-1 (EVID-05/D-01): every reachable role_reason yields a candidate-named, contig-language, evidence_state-bearing sentence")

# D-02 pattern 1: no_own_assembly states the no-matching-contig phrasing plainly,
# with no fabricated numbers and no "own assembly".
s_no <- summ("background", "no_own_assembly", "weak", assembly_exists = FALSE)
if (!grepl("no contig", s_no, fixed = TRUE))
  fail(sprintf("no_own_assembly must use the no-matching-contig phrasing, got: %s", s_no))
if (grepl("NA", s_no, fixed = TRUE))
  fail(sprintf("no_own_assembly must not render literal 'NA', got: %s", s_no))
ok("EvidenceReport-2 (D-02 pattern 1): no_own_assembly -> 'no contig with the same genotype/subtype', no NA leak")

# D-02 pattern 3 / D-07: weak-but-matching states the ACTUAL measured numbers
# (length 200 bp, identity 50%), not editorial words like "essentially no support".
s_weak <- summ("background", "weak_own_assembly_below_floor", "weak", len = 200, pid = 50, kmer = 0.5)
if (!grepl("200", s_weak, fixed = TRUE))
  fail(sprintf("weak_own_assembly must state the measured contig length 200, got: %s", s_weak))
if (!grepl("50.0", s_weak, fixed = TRUE))
  fail(sprintf("weak_own_assembly must state the measured identity 50.0%%, got: %s", s_weak))
if (grepl("essentially no support", s_weak, fixed = TRUE))
  fail("weak_own_assembly must NOT editorialize ('essentially no support') (D-07)")
ok("EvidenceReport-3 (D-02 pattern 3 / D-07): weak-but-matching contig states measured length/identity, no editorializing")

# T-13-01/T-13-02: NA-heavy input tolerated — a defined sentence, no stop(), no 'NA' token.
s_na <- build_evidence_summary(
  role = NA_character_, role_reason = NA_character_, evidence_state = NA_character_,
  concordance_status = NA_character_, candidate_rank = 1L,
  candidate_ref = "1a_M62321", candidate_subtype = "1a", assembly_exists = NA)
if (!is.character(s_na) || length(s_na) != 1 || !nzchar(s_na))
  fail("build_evidence_summary must return a defined sentence on all-NA input, never stop()")
if (grepl("NA", s_na, fixed = TRUE))
  fail(sprintf("all-NA input must not render a literal 'NA' token, got: %s", s_na))
ok("EvidenceReport-4 (T-13-01/02): NA-heavy input yields a defined sentence, no stop(), no literal NA")

# Vectorizable via pmap_chr over a candidate frame (the production call shape).
frame <- tibble(
  role = c("dominant", "co-infection", "background"),
  role_reason = c("dominant", "corroborated", "no_own_assembly"),
  evidence_state = c("confirmed", "probable", "weak"),
  concordance_status = NA_character_,
  candidate_rank = c(1L, 2L, 3L),
  candidate_ref = c("1a_M62321", "2c_JX227949", "4d_DQ418786"),
  candidate_subtype = c("1a", "2c", "4d"),
  assembly_exists = c(TRUE, TRUE, FALSE),
  best_contig_length = c(9189, 3000, NA),
  best_contig_pident = c(100, 80, NA),
  best_contig_kmer_cov = c(40, 5, NA)
)
vec <- pmap_chr(frame, build_evidence_summary)
if (length(vec) != 3 || any(is.na(vec)) || any(!nzchar(vec)))
  fail("build_evidence_summary must vectorize cleanly via pmap_chr over a candidate frame")
if (!grepl("candidate 3", vec[3], fixed = TRUE) || !grepl("4d_DQ418786", vec[3], fixed = TRUE))
  fail("vectorized summary must name each candidate independently by rank + ref")
ok("EvidenceReport-5: build_evidence_summary vectorizes via pmap_chr over a candidate frame (production call shape)")

cat("ALL PASS\n")
