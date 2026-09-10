#!/usr/bin/env Rscript

# test_dominance_breadth_source.R -----------------------------------------
# Quick task 260810-dbs. Function-level tests for the two-axis coverage
# resolution in bin/classify_roles.R:
#   .prefer_targeted()   -> which measurement the SCORE and the FLOOR read
#   .coverage_evidence() -> whether a candidate has ANY coverage (eligibility)
#
# Background. Every candidate carries up to two breadth@>=5x measurements:
#   candidate_cov     first pass, all-reference, with duplicates. NA for a
#                     rescued/nominated candidate (the number described the
#                     DISPLACED reference).
#   cand_cov_breadth  second pass, joint-mapping, deduplicated. NA for a
#                     candidate that was never targeted-mapped.
# D-08 (08-CONTEXT.md:33) makes the second-pass value authoritative for the
# score and the floor. Before 260810-dbs nothing emitted cand_cov_breadth, so
# score_candidates() always fell through to the first-pass column, and the
# rescued/nominated candidates (NA there) scored as if NOTHING of their
# reference were covered.
#
# Asserted behaviours:
#   DBS-1  score reads cand_cov_breadth when present (ERR1810469 3a anchor)
#   DBS-2  TRAP: preference is per ROW — a never-targeted candidate falls back
#          to candidate_cov, it does NOT drop to 0
#   DBS-3  TRAP: breadth is a PERCENT unconditionally — 1% must not read as 100%
#   DBS-4  a rescued candidate (candidate_cov NA) earns real breadth points
#   DBS-5  a rescued candidate is ELIGIBLE and can be selected dominant
#   DBS-6  a sole rescued candidate is not stranded at `untypable`
#   DBS-7  eligibility only ever ADDS: targeted 0 + first-pass >0 stays eligible
#   DBS-8  genuine no-coverage (both axes 0/NA) still yields `untypable`
#   DBS-9  below_floor reads the targeted axes, per its min_targeted_* params
#
# Verified against the pre-fix code (HEAD 1a6519d): DBS-2/3/5/6/9 all FAIL, and
# DBS-6 reproduces the shipped contradiction verbatim — role "co-infection"
# together with overall_sample_call "untypable".
#
# DBS-1, DBS-4, DBS-7 and DBS-8 pass both before and after BY DESIGN. They are
# contract and non-change pins, not regression tests: score_candidates() always
# honoured cand_cov_breadth when the column was present in the frame, and the
# actual defect was that no module ever PUT it there. That wiring cannot be
# reached from a function-level test — the guard for it is the
# `summary/candidates.csv` header assertion in tests/default.nf.test. If that
# assertion is ever removed, this whole file goes green against a pipeline that
# has silently reverted to first-pass breadth.
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

near <- function(a, b, tol = 1e-6) isTRUE(abs(a - b) < tol)

# Production weights (nextflow.config Plan-01 defaults).
W <- list(evenness = 3.0, reads = 1.0, kmercov = 0.5)

# Builder. `cov` is the FIRST-pass candidate_cov; `tgt` the SECOND-pass
# cand_cov_breadth. Pass NA to either to model the real absence cases. The
# cand_cov_breadth column is only attached when `tgt` is supplied, so the
# "column absent entirely" path (pre-260810-dbs frames, and every caller that
# never joins the cov loop) is exercised too.
mk <- function(cov, reads, even = 0.5, kmer = NA_real_, tgt = NULL,
               ref = "ref", subtype = "1a", sample = "S") {
  d <- tibble(
    sampleName                            = sample,
    candidate_ref                         = ref,
    candidate_subtype                     = subtype,
    candidate_genotype                    = genotype_from_subtype(subtype),
    candidate_reads                       = reads,
    targeted_reads_nodup                  = reads,
    candidate_cov                         = cov,
    cv_evenness                           = even,
    assembly_support_best_contig_kmer_cov = kmer
  )
  if (!is.null(tgt)) d$cand_cov_breadth <- tgt
  d
}

score_of <- function(df) score_candidates(df, W)$dominance_score

# --- DBS-1: the score reads the TARGETED breadth when it is present ---------
# ERR1810469's 3a candidate (Thomson 2016), the pipeline's own worked dominance
# example: 46% first-pass breadth against 18.17% targeted, 190 targeted nodup
# reads, cv_evenness 0.4785, best contig k-mer cov 51.65 (capped at 50).
err469_3a_firstpass <- score_of(mk(46, 190, 0.4785, 51.65))
err469_3a_targeted  <- score_of(mk(46, 190, 0.4785, 51.65, tgt = 18.17))

# The shipped (pre-fix) value, reproduced from candidates.csv: 5.947897.
if (!near(err469_3a_firstpass, 5.947897, tol = 5e-4))
  fail(sprintf("DBS-1 first-pass-only frame must reproduce the shipped 5.947897, got %.6f",
               err469_3a_firstpass))
# With the targeted breadth the 3a candidate loses exactly we*(0.46 - 0.1817).
expected_drop <- 3.0 * (0.46 - 0.1817)
if (!near(err469_3a_firstpass - err469_3a_targeted, expected_drop, tol = 1e-6))
  fail(sprintf("DBS-1 targeted breadth must move the score by %.6f, moved %.6f",
               expected_drop, err469_3a_firstpass - err469_3a_targeted))
ok("DBS-1 (D-08): score reads cand_cov_breadth over candidate_cov when present")

# --- DBS-2: THE TRAP — preference is per ROW, not per column ----------------
# Two candidates in one frame: rank 1 was targeted-mapped, rank 2 was
# below-threshold so it has no depth file and NA targeted breadth. A per-COLUMN
# preference (the shape of the pre-fix code) sees the column present and sends
# rank 2's breadth to 0, wiping a real 46% first-pass measurement.
mixed <- bind_rows(
  mk(69, 1017, 0.6922, 5.40, tgt = 96.02,      ref = "1a_HQ850279", subtype = "1a"),
  mk(46,  190, 0.4785, 51.65, tgt = NA_real_,  ref = "3a_D17763",   subtype = "3a")
)
mixed_scores <- score_candidates(mixed, W)$dominance_score
# Rank 2 must score exactly as a first-pass-only frame does.
if (!near(mixed_scores[2], err469_3a_firstpass, tol = 1e-9))
  fail(sprintf("DBS-2 a never-targeted candidate must fall back to candidate_cov (%.6f), got %.6f",
               err469_3a_firstpass, mixed_scores[2]))
# Belt: it must NOT have collapsed to the zero-breadth score.
zero_breadth <- score_of(mk(NA_real_, 190, 0.4785, 51.65, tgt = NA_real_))
if (near(mixed_scores[2], zero_breadth, tol = 1e-9))
  fail("DBS-2 per-COLUMN preference regression: never-targeted candidate dropped to breadth 0")
ok("DBS-2 (per-row trap): a never-targeted candidate falls back to candidate_cov, not to 0")

# --- DBS-3: THE UNITS TRAP — breadth is a percent, unconditionally ----------
# The old `ifelse(breadth_src > 1, /100, as-is)` heuristic read any value in
# [0,1] as an already-fractional breadth. percent_gt_4_int is round()ed, so
# 0.5-1.5% first-pass breadth landed on exactly 1 and collected the FULL 3.0
# breadth points: 1% scored 6.599 against 2%'s 3.659, an inversion of 2.94.
# cov_breadth_min_5 carries two decimals, which makes the (0,1) band routine.
s_1pct   <- score_of(mk(1, 500, 0.3))
s_2pct   <- score_of(mk(2, 500, 0.3))
s_frac   <- score_of(mk(NA_real_, 500, 0.3, tgt = 0.42))   # 0.42% targeted breadth
s_0pct   <- score_of(mk(0, 500, 0.3))
if (!(s_1pct < s_2pct))
  fail(sprintf("DBS-3 1%% breadth (%.4f) must score BELOW 2%% breadth (%.4f) — units inversion",
               s_1pct, s_2pct))
if (!near(s_1pct - s_0pct, 3.0 * 0.01, tol = 1e-9))
  fail(sprintf("DBS-3 1%% breadth must contribute we*0.01 = 0.03 points, contributed %.6f",
               s_1pct - s_0pct))
if (!near(s_frac - s_0pct, 3.0 * 0.0042, tol = 1e-9))
  fail(sprintf("DBS-3 a 0.42%% targeted breadth must contribute we*0.0042, contributed %.6f",
               s_frac - s_0pct))
ok("DBS-3 (units trap): breadth is coerced as a 0-100 percent unconditionally, no [0,1] heuristic")

# --- DBS-4: a rescued candidate earns real breadth points -------------------
# rescue_evaluation.R blanks candidate_cov after a replacement (:549-553) and
# creates nominations with NA from the start (:482). ERR1810447's 2b: 95.9%
# targeted breadth, which used to contribute nothing at all.
resc_blind  <- score_of(mk(NA_real_, 4199, 0.80, 20))
resc_fixed  <- score_of(mk(NA_real_, 4199, 0.80, 20, tgt = 95.9))
if (!near(resc_fixed - resc_blind, 3.0 * 0.959, tol = 1e-9))
  fail(sprintf("DBS-4 a rescued candidate must earn we*0.959 = 2.877 breadth points, earned %.6f",
               resc_fixed - resc_blind))
ok("DBS-4: a rescued/nominated candidate (candidate_cov NA) is scored on its targeted breadth")

# --- Role-level fixtures ----------------------------------------------------
classify <- function(df, minRead = 500, minCov = 30) {
  classify_roles(score_candidates(df, W), minRead = minRead, minCov = minCov,
                 denovo_min_contig_length = 500, denovo_min_kmer_cov = 2.0,
                 denovo_min_blast_identity = 90, match_level = "genotype")
}
with_support <- function(d, len, pid, kmer) {
  d %>% mutate(
    assembly_support_best_contig_length   = len,
    assembly_support_best_contig_pident   = pid,
    assembly_support_best_contig_kmer_cov = kmer
  )
}

# --- DBS-5: a rescued candidate is eligible and can WIN dominance -----------
# Pre-fix, `eligible` tested `!is.na(candidate_cov) & candidate_cov > 0`, so a
# rescued/nominated candidate was structurally barred from the dominant pool no
# matter how strong its targeted evidence.
duel <- bind_rows(
  # Rescued 2b: no first-pass numbers at all, but 95.9% targeted breadth.
  with_support(mk(NA_real_, 40000, 0.90, 30, tgt = 95.9,
                  ref = "2b_D10988", subtype = "2b", sample = "DUEL"), 9207, 96, 30),
  # First-pass 1b that never assembled anything.
  with_support(mk(60, 900, 0.35, NA_real_, tgt = 22.0,
                  ref = "1b_ref", subtype = "1b", sample = "DUEL"), NA_real_, NA_real_, NA_real_)
)
r_duel <- classify(duel)
dom_duel <- r_duel %>% filter(role == "dominant") %>% pull(candidate_ref)
if (length(dom_duel) != 1 || dom_duel != "2b_D10988")
  fail(sprintf("DBS-5 the rescued 2b must be selectable as dominant, dominant was: %s",
               paste(dom_duel, collapse = ", ")))
ok("DBS-5: a rescued candidate (candidate_cov NA) is eligible and can be selected dominant")

# --- DBS-6: a SOLE rescued candidate is not stranded at `untypable` ---------
# A rank-1 replacement is reachable: the dominant_protect_cov guard
# (rescue_evaluation.R:239) only holds when first-pass cov >= 90 AND the own
# contig clears the quality floors, and the 2k1b rule bypasses the own-subtype
# check entirely. Pre-fix this sample reported role "co-infection" together
# with overall_sample_call "untypable" — two mutually contradictory statements,
# and `untypable` is documented as "no usable coverage on any candidate".
sole <- with_support(
  mk(NA_real_, 4199, 0.80, 20, tgt = 95.9,
     ref = "2b_D10988", subtype = "2b", sample = "SOLE"), 9207, 96, 30)
r_sole <- classify(sole)
if (r_sole$overall_sample_call[1] == "untypable")
  fail("DBS-6 a sole rescued candidate with 95.9% targeted breadth must not be `untypable`")
if (r_sole$role[1] != "dominant")
  fail(sprintf("DBS-6 a sole rescued candidate must be the dominant, got role %s", r_sole$role[1]))
if (r_sole$overall_sample_call[1] != "monoinfection")
  fail(sprintf("DBS-6 a sole rescued dominant must call monoinfection, got %s",
               r_sole$overall_sample_call[1]))
ok("DBS-6: a sole rescued candidate is typable — no `co-infection` + `untypable` contradiction")

# --- DBS-7: eligibility only ever ADDS ------------------------------------
# .coverage_evidence() takes the MAX over the axes, deliberately NOT the
# targeted-first preference the score uses. A candidate that lost every read to
# its joint-mapping competitor (targeted 0) but has a real first-pass
# measurement stays eligible, exactly as before 260810-dbs. This pins the
# deliberate NON-change: switching eligibility to targeted-only can newly
# strand a sample at `untypable` and needs a cohort re-run to justify.
lost_all <- with_support(
  mk(46, 5009, 0.15, NA_real_, tgt = 0,
     ref = "3a_D17763", subtype = "3a", sample = "LOST"), NA_real_, NA_real_, NA_real_)
r_lost <- classify(lost_all)
if (r_lost$overall_sample_call[1] == "untypable")
  fail("DBS-7 targeted breadth 0 with a real first-pass 46% must stay eligible (max, not prefer)")
if (r_lost$role[1] != "dominant")
  fail(sprintf("DBS-7 the only eligible candidate must be dominant, got %s", r_lost$role[1]))
ok("DBS-7 (non-change pin): eligibility takes the MAX over axes, so it never demotes")

# --- DBS-8: genuine no-coverage still yields `untypable` --------------------
# The repair must not make `untypable` unreachable. Both axes zero, and both
# axes NA, must still produce it.
r_zero <- classify(bind_rows(
  mk(0, 248,  0, NA_real_, tgt = 0, ref = "3a_ref", subtype = "3a", sample = "ZERO"),
  mk(0, 1030, 0, NA_real_, tgt = 0, ref = "1a_ref", subtype = "1a", sample = "ZERO")
))
if (unique(r_zero$overall_sample_call) != "untypable")
  fail("DBS-8 an all-zero-coverage sample must still be untypable")
r_nacov <- classify(bind_rows(
  mk(NA_real_, 248,  0, NA_real_, tgt = NA_real_, ref = "3a_ref", subtype = "3a", sample = "NAC"),
  mk(NA_real_, 1030, 0, NA_real_, tgt = NA_real_, ref = "1a_ref", subtype = "1a", sample = "NAC")
))
if (unique(r_nacov$overall_sample_call) != "untypable")
  fail("DBS-8 a sample with NO coverage measurement on either axis must still be untypable")
ok("DBS-8: `untypable` stays reachable — both-zero and both-NA coverage still produce it")

# --- DBS-9: below_floor reads the TARGETED axes ----------------------------
# The thresholds are named min_targeted_read / min_targeted_cov at the
# summarize.R call site and D-08 says the floor and the score share one
# coverage source. Pre-fix this annotation read candidate_reads (first-pass,
# WITH duplicates) and candidate_cov, so every rescued candidate reported FALSE
# against a floor it clears by a wide margin.
r_floor <- classify(sole, minRead = 500, minCov = 30)
if (!isTRUE(r_floor$below_floor[1]))
  fail("DBS-9 a rescued candidate at 4199 targeted reads / 95.9% targeted breadth must clear the floor")
# And the targeted value governs when the two axes disagree across the floor:
# first-pass 95% would pass, targeted 12% must not.
r_floor_disagree <- classify(
  with_support(mk(95, 4199, 0.80, 20, tgt = 12.0,
                  ref = "2b_D10988", subtype = "2b", sample = "DIS"), 9207, 96, 30),
  minRead = 500, minCov = 30)
if (isTRUE(r_floor_disagree$below_floor[1]))
  fail("DBS-9 targeted breadth 12% must fail a 30% floor even when first-pass reads 95%")
ok("DBS-9 (D-08): below_floor reads targeted_reads_nodup + cand_cov_breadth, per its parameter names")

cat("\nALL PASS\n")
