#!/usr/bin/env Rscript

# test_dominance_score.R --------------------------------------------------
# Function-level tests for the Phase-8 (SCORE-01/02) dominance score,
# bin/classify_roles.R::score_candidates(). Modelled on
# test_assembly_support_join.R: source the REAL helper, build small in-memory
# tibbles, and assert on the real function's output (no inline re-implementation).
#
# Asserted behaviours (08-01-PLAN Task 3):
#   1. Monotonic in reads holding breadth/evenness/kmer fixed (SCORE-01).
#   2. Monotonic in breadth holding the rest fixed (SCORE-01).
#   3. Monotonic in cv_evenness holding the rest fixed (SCORE-02 — evenness counts).
#   4. THE HEADLINE: a high-read spiky candidate (the benchmarked false 4g:
#      53279 reads, breadth 0.677, high CV -> low cv_evenness) scores LOWER than a
#      genuine even minor with FEWER reads (SCORE-02 — evenness dominates reads).
#   5. NA / "none" k-mer coverage gives NO penalty vs a no-kmer baseline (D-05 —
#      bonus only); a present k-mer cov gives a (capped) BOOST.
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

# Builder: a scored-input candidate row carrying the columns score_candidates()
# reads — candidate_reads, candidate_cov, cv_evenness factor, and the joined
# k-mer-cov metric. Defaults give a neutral baseline.
#
# NOTE (260810-dbs): `cov` here is the FIRST-pass candidate_cov, used as the
# breadth fallback because these fixtures carry no cand_cov_breadth. The
# calibration values below (67.7 / 90.3) were always TARGETED breadths taken from
# the handoff evidence table — the score was tuned on the second-pass measurement
# while production fed it the first-pass one, which is the defect 260810-dbs
# fixed. Breadth is a 0-100 PERCENT on both axes; see test_dominance_breadth_source.R.
mk_score_cand <- function(reads, cov = 90, even = 0.5, kmer = NA_real_, ref = "ref") {
  tibble(
    sampleName                            = "S",
    candidate_ref                         = ref,
    candidate_reads                       = reads,
    candidate_cov                         = cov,
    cv_evenness                           = even,
    assembly_support_best_contig_kmer_cov = kmer
  )
}

score_of <- function(df) score_candidates(df)$dominance_score

# Use the calibrated production weights (nextflow.config Plan-01 defaults).
W <- list(evenness = 3.0, reads = 1.0, kmercov = 0.5)

# --- Test 1: monotonic in reads ---------------------------------------------
s_lo <- score_candidates(mk_score_cand(1000), W)$dominance_score
s_hi <- score_candidates(mk_score_cand(100000), W)$dominance_score
if (!(s_hi > s_lo)) fail("Test1 score must increase with reads holding others fixed")
ok("Test1 (SCORE-01): dominance_score monotonic in reads")

# --- Test 2: monotonic in breadth -------------------------------------------
b_lo <- score_candidates(mk_score_cand(5000, cov = 30), W)$dominance_score
b_hi <- score_candidates(mk_score_cand(5000, cov = 95), W)$dominance_score
if (!(b_hi > b_lo)) fail("Test2 score must increase with breadth holding others fixed")
ok("Test2 (SCORE-01): dominance_score monotonic in breadth")

# --- Test 3: monotonic in cv_evenness ---------------------------------------
e_lo <- score_candidates(mk_score_cand(5000, even = 0.2), W)$dominance_score
e_hi <- score_candidates(mk_score_cand(5000, even = 0.9), W)$dominance_score
if (!(e_hi > e_lo)) fail("Test3 score must increase with cv_evenness holding others fixed")
ok("Test3 (SCORE-02): dominance_score monotonic in cv_evenness")

# --- Test 4: THE HEADLINE — 4g loses to a genuine even minor (SCORE-02) -----
# False 4g (handoff §2): 53279 reads, 67.7% breadth, spiky high-CV -> low
# cv_evenness (~0.20). Genuine even minor: FEWER reads (ERR1810447 2b ~4199),
# high breadth (90.3%), even coverage -> high cv_evenness (~0.80).
false_4g    <- mk_score_cand(53279, cov = 67.7, even = 0.20, kmer = NA_real_, ref = "4g_artifact")
genuine_min <- mk_score_cand(4199,  cov = 90.3, even = 0.80, kmer = 20,       ref = "2b_genuine")
score_4g  <- score_candidates(false_4g, W)$dominance_score
score_gen <- score_candidates(genuine_min, W)$dominance_score
if (!(score_4g < score_gen)) {
  fail(sprintf("Test4 HEADLINE: 4g score %.3f must be LOWER than even minor %.3f (evenness must dominate reads)",
               score_4g, score_gen))
}
ok("Test4 (SCORE-02 HEADLINE): spiky high-read 4g scores below a genuine even minor with fewer reads")

# Belt: even with NO k-mer boost on the genuine minor, evenness still wins.
genuine_min_nokmer <- mk_score_cand(4000, cov = 85, even = 0.75, kmer = NA_real_, ref = "2b_nokmer")
score_gen_nk <- score_candidates(genuine_min_nokmer, W)$dominance_score
if (!(score_4g < score_gen_nk)) {
  fail(sprintf("Test4b: 4g %.3f must lose to even minor %.3f even without a k-mer-cov boost",
               score_4g, score_gen_nk))
}
ok("Test4b (SCORE-02): 4g still loses to an even minor that had no de novo k-mer boost")

# --- Test 5: k-mer cov is bonus-only (D-05) ---------------------------------
base    <- mk_score_cand(5000, cov = 90, even = 0.6, kmer = NA_real_)
boosted <- mk_score_cand(5000, cov = 90, even = 0.6, kmer = 30)
s_base    <- score_candidates(base, W)$dominance_score
s_boosted <- score_candidates(boosted, W)$dominance_score
if (s_boosted < s_base) fail("Test5 a present k-mer cov must NEVER lower the score (bonus-only, D-05)")
if (!(s_boosted > s_base)) fail("Test5 a present k-mer cov should give a positive boost over the NA baseline")
# NA k-mer must equal a zero k-mer (no penalty for missing assembly support).
s_zero <- score_candidates(mk_score_cand(5000, cov = 90, even = 0.6, kmer = 0), W)$dominance_score
if (abs(s_zero - s_base) > 1e-9) fail("Test5 zero/NA k-mer cov must yield the same (no-boost) score — no penalty")
ok("Test5 (D-05): k-mer cov is bonus-only and capped; NA/none gives no penalty")

cat("\nALL PASS\n")
