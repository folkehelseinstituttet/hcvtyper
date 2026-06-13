#!/usr/bin/env Rscript

# test_classify_roles.R ---------------------------------------------------
# Function-level tests for the Phase-8 (CLASS-01..04, COMPAT-04) strain-role
# classifier, bin/classify_roles.R::classify_roles() (run after
# score_candidates()). Modelled on test_assembly_support_join.R: source the REAL
# helpers, build small in-memory tibbles anchored to the handoff §2 evidence
# table, and assert on the real functions (no inline re-implementation).
#
# Asserted behaviours (08-01-PLAN Task 3):
#   - false 4g (no 4g contig, dominant 1a HAS a substantial contig)
#       -> role=="background" && role_reason=="refuted_denovo"          (D-11)
#   - ERR1810447 full 9207bp 2b + ERR1810453 2949bp partial 2b (kmer ~5)
#       -> role=="co-infection" && role_reason=="corroborated"          (CLASS-02)
#   - sim1 1a:1b + sim2 2a:3a true co-infections preserved as co-infection
#   - IVT extreme-ratio genuine minor, de novo failed for BOTH
#       -> role=="co-infection" && role_reason=="uncorroborated_kept"   (D-11)
#   - same-genotype non-1a/1b candidate -> background/same_genotype_as_dominant (D-12)
#   - 2k1b pair -> background/recombinant_2k1b                          (D-12)
#   - no candidate passes the gate -> overall_sample_call=="indeterminate" (D-01/D-14)
#   - every candidate gets exactly one role                            (CLASS-01)
#   - zero-row input -> typed frame, no abort                          (CLASS-03)
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

# Builder: one candidate row carrying identity, abundance, the cv_evenness factor,
# and the joined assembly-support metrics (NA => no support). subtype drives the
# precomputed genotype via the real genotype_from_subtype().
mk_cand <- function(sample, ref, subtype, reads, cov, even,
                    sup_len = NA_real_, sup_kmer = NA_real_, sup_pid = NA_real_) {
  tibble(
    sampleName                              = sample,
    candidate_ref                           = ref,
    candidate_subtype                       = subtype,
    candidate_genotype                      = genotype_from_subtype(subtype),
    candidate_reads                         = reads,
    candidate_cov                           = cov,
    cv_evenness                             = even,
    assembly_support_best_contig_length     = sup_len,
    assembly_support_best_contig_kmer_cov   = sup_kmer,
    assembly_support_best_contig_pident     = sup_pid
  )
}

# Run the full pipeline (score then classify) at the validated 1000/2.0/90 floors
# and the default minRead/minCov (500/30) from the handoff.
classify <- function(df, minRead = 500, minCov = 30) {
  classify_roles(score_candidates(df), minRead = minRead, minCov = minCov,
                 denovo_min_contig_length = 1000, denovo_min_kmer_cov = 2.0,
                 denovo_min_blast_identity = 90, match_level = "genotype")
}

role_of   <- function(r, the_ref) r %>% filter(candidate_ref == the_ref) %>% pull(role)
reason_of <- function(r, the_ref) r %>% filter(candidate_ref == the_ref) %>% pull(role_reason)

# --- Test 1: false 4g refuted (D-11 asymmetric) -----------------------------
# Dominant 1a (genuine, full contig). 4g minor clears the floor on reads/cov but
# has NO 4g contig, while de novo demonstrably WORKED for the dominant 1a.
false_4g <- bind_rows(
  mk_cand("sim1_1a", "1a_M62321", "1a", 200000, 99,   0.90, sup_len = 9189, sup_kmer = 40, sup_pid = 100),
  mk_cand("sim1_1a", "4g_artif",  "4g", 53279,  67.7, 0.20)  # no support
)
r1 <- classify(false_4g)
if (role_of(r1, "4g_artif") != "background") fail("Test1 false 4g must be background")
if (reason_of(r1, "4g_artif") != "refuted_denovo") fail("Test1 false 4g reason must be refuted_denovo")
if (role_of(r1, "1a_M62321") != "dominant") fail("Test1 genuine 1a must be dominant")
if ((r1 %>% pull(overall_sample_call) %>% unique()) != "monoinfection")
  fail("Test1 a refuted-only sample must be monoinfection")
ok("Test1 (D-11): false 4g -> background/refuted_denovo; sample monoinfection")

# --- Test 2: genuine 2b co-infections corroborated (CLASS-02) ---------------
# ERR1810447: full 9207bp 2b. ERR1810453: partial 2949bp 2b at k-mer cov ~5
# (passes 1000/2.0/90 — the reconciled floor; would FAIL the old 10.0 floor).
err447 <- bind_rows(
  mk_cand("ERR1810447", "1b_ref", "1b", 300000, 98, 0.92, sup_len = 8900, sup_kmer = 45, sup_pid = 99),
  mk_cand("ERR1810447", "2b_ref", "2b", 4199,   90.3, 0.78, sup_len = 9207, sup_kmer = 30, sup_pid = 91)
)
r2 <- classify(err447)
if (role_of(r2, "2b_ref") != "co-infection") fail("Test2 ERR1810447 full 2b must be co-infection")
if (reason_of(r2, "2b_ref") != "corroborated") fail("Test2 ERR1810447 2b reason must be corroborated")

err453 <- bind_rows(
  mk_cand("ERR1810453", "1a_ref", "1a", 250000, 97, 0.90, sup_len = 9100, sup_kmer = 50, sup_pid = 99),
  mk_cand("ERR1810453", "2b_ref", "2b", 3263,   85,  0.74, sup_len = 2949, sup_kmer = 5,  sup_pid = 93)
)
r3 <- classify(err453)
if (role_of(r3, "2b_ref") != "co-infection")
  fail("Test2 ERR1810453 partial 2b (kmer ~5) must be co-infection at the reconciled 2.0 floor")
if (reason_of(r3, "2b_ref") != "corroborated")
  fail("Test2 ERR1810453 partial 2b reason must be corroborated")
if ((r3 %>% pull(overall_sample_call) %>% unique()) != "co-infection")
  fail("Test2 a corroborated minor sample must be co-infection")
ok("Test2 (CLASS-02): full + partial genuine 2b both -> co-infection/corroborated")

# --- Test 2b: the OLD 10.0 k-mer floor WOULD have refuted ERR1810453 --------
# Locks the Pitfall-1 reconciliation: at the stricter 10.0 floor the kmer-5
# partial 2b is no longer substantial and (dominant assembled) gets refuted.
r3_strict <- classify_roles(score_candidates(err453), minRead = 500, minCov = 30,
                            denovo_min_contig_length = 1000, denovo_min_kmer_cov = 10.0,
                            denovo_min_blast_identity = 90, match_level = "genotype")
if (role_of(r3_strict, "2b_ref") != "background")
  fail("Test2b precondition: at the old 10.0 floor the kmer-5 2b should be refuted (motivates the 2.0 reconciliation)")
ok("Test2b: the validated 2.0 floor (not 10.0) is what preserves the genuine partial 2b")

# --- Test 3: true co-infections preserved (sim1 1a:1b, sim2 2a:3a) ----------
sim1 <- bind_rows(
  mk_cand("sim1", "1a_ref", "1a", 150000, 99, 0.90, sup_len = 9076, sup_kmer = 40, sup_pid = 99),
  mk_cand("sim1", "1b_ref", "1b", 120000, 97, 0.88, sup_len = 9339, sup_kmer = 38, sup_pid = 99)
)
r4 <- classify(sim1)
# 1a/1b is the explicit D-12 exception: cross-subtype within gt1 is ALLOWED.
nondom1 <- r4 %>% filter(role != "dominant")
if (nrow(nondom1) != 1 || nondom1$role != "co-infection")
  fail("Test3 sim1 1a:1b must be preserved as a co-infection (D-12 1a/1b allowance)")
if ((r4 %>% pull(overall_sample_call) %>% unique()) != "co-infection")
  fail("Test3 sim1 must be a co-infection sample")

sim2 <- bind_rows(
  mk_cand("sim2", "2a_ref", "2a", 140000, 98, 0.89, sup_len = 9500, sup_kmer = 42, sup_pid = 99),
  mk_cand("sim2", "3a_ref", "3a", 110000, 96, 0.87, sup_len = 9450, sup_kmer = 36, sup_pid = 99)
)
r5 <- classify(sim2)
nondom2 <- r5 %>% filter(role != "dominant")
if (nrow(nondom2) != 1 || nondom2$role != "co-infection")
  fail("Test3 sim2 2a:3a must be preserved as a co-infection")
ok("Test3: sim1 1a:1b and sim2 2a:3a true co-infections preserved")

# --- Test 4: IVT extreme-ratio genuine minor, de novo failed for BOTH -------
# Genuine low-yield minor clears the floor but assembled nothing; the DOMINANT
# also assembled nothing substantial -> de novo inconclusive -> keep, do not refute.
ivt <- bind_rows(
  mk_cand("IVT", "1a_ref", "1a", 500000, 99, 0.95),   # dominant, NO support
  mk_cand("IVT", "3a_ref", "3a", 800,    35, 0.40)    # minor clears floor, NO support
)
r6 <- classify(ivt)
if (role_of(r6, "3a_ref") != "co-infection")
  fail("Test4 IVT minor with both-de-novo-failed must be kept as co-infection")
if (reason_of(r6, "3a_ref") != "uncorroborated_kept")
  fail("Test4 IVT minor reason must be uncorroborated_kept (D-11 never suppress on inconclusive de novo)")
ok("Test4 (D-11): de-novo-failed-for-both genuine minor -> co-infection/uncorroborated_kept")

# --- Test 5: D-12 same-genotype (non-1a/1b) demotion ------------------------
# A 3b candidate against a dominant 3a: same genotype 3, NOT a 1a/1b pair ->
# is_valid_minor() returns FALSE -> demote to background. Even WITH corroborating
# 3-genotype support, the exception may only demote (never promote).
same_gt <- bind_rows(
  mk_cand("S_3", "3a_ref", "3a", 200000, 98, 0.90, sup_len = 9000, sup_kmer = 40, sup_pid = 99),
  mk_cand("S_3", "3b_ref", "3b", 5000,   80, 0.70, sup_len = 8500, sup_kmer = 30, sup_pid = 95)
)
r7 <- classify(same_gt)
if (role_of(r7, "3b_ref") != "background")
  fail("Test5 same-genotype 3b vs 3a must be demoted to background")
if (reason_of(r7, "3b_ref") != "same_genotype_as_dominant")
  fail("Test5 same-genotype demotion reason must be same_genotype_as_dominant")
ok("Test5 (D-12): same-genotype non-1a/1b candidate -> background/same_genotype_as_dominant")

# --- Test 6: D-12 2k1b pair demotion ----------------------------------------
# 2k1b recombinant paired with a genotype-2 dominant -> blocked by is_valid_minor.
pair_2k1b <- bind_rows(
  mk_cand("S_2k", "2a_ref",  "2a",   200000, 98, 0.90, sup_len = 9000, sup_kmer = 40, sup_pid = 99),
  mk_cand("S_2k", "2k1b_ref","2k1b", 6000,   82, 0.72, sup_len = 8800, sup_kmer = 33, sup_pid = 96)
)
r8 <- classify(pair_2k1b)
if (role_of(r8, "2k1b_ref") != "background")
  fail("Test6 2k1b vs genotype-2 dominant must be demoted to background")
if (reason_of(r8, "2k1b_ref") != "recombinant_2k1b")
  fail("Test6 2k1b demotion reason must be recombinant_2k1b")
ok("Test6 (D-12): 2k1b pair -> background/recombinant_2k1b")

# --- Test 7: no candidate passes the gate -> indeterminate (D-01/D-14) -------
# ERR1810469-style: the 3a major failed coverage; nothing clears minRead/minCov.
no_gate <- bind_rows(
  mk_cand("ERR1810469", "3a_ref", "3a", 248,  8.8, 0.15),
  mk_cand("ERR1810469", "1a_ref", "1a", 1030, 28,  0.40)  # cov below 30 -> fails gate
)
r9 <- classify(no_gate)
if ((r9 %>% pull(overall_sample_call) %>% unique()) != "indeterminate")
  fail("Test7 no candidate passing the gate must yield overall_sample_call=indeterminate")
if (any(r9$role == "dominant")) fail("Test7 there must be NO dominant when none passes the gate")
ok("Test7 (D-01/D-14): no-gate-pass -> no dominant, overall_sample_call=indeterminate")

# --- Test 8: every candidate gets exactly one role (CLASS-01) ---------------
all_rows <- bind_rows(r1, r2, r3, r4, r5, r6, r7, r8, r9)
if (any(is.na(all_rows$role))) fail("Test8 every candidate must carry a non-NA role (CLASS-01)")
if (!all(all_rows$role %in% c("dominant", "co-infection", "background")))
  fail("Test8 role must be exactly one of dominant/co-infection/background")
if (!all(all_rows$overall_sample_call %in% c("monoinfection", "co-infection", "indeterminate")))
  fail("Test8 overall_sample_call must be one of monoinfection/co-infection/indeterminate")
ok("Test8 (CLASS-01/CLASS-04): every candidate has exactly one valid role + 3-value sample call")

# --- Test 9: zero-row input -> typed frame, no abort (CLASS-03) -------------
empty_in <- mk_cand("S", "x", "1a", 1, 1, 0)[0, ]
r10 <- classify(empty_in)
if (nrow(r10) != 0) fail("Test9 zero-row input must yield zero rows")
for (col in c("role", "dominance_score", "role_reason", "overall_sample_call")) {
  if (!col %in% names(r10)) fail(sprintf("Test9 zero-row frame must still carry the %s column", col))
}
# NULL input must also not abort.
r11 <- classify_roles(NULL, minRead = 500, minCov = 30)
if (nrow(r11) != 0) fail("Test9 NULL input must yield a zero-row frame, not an abort")
ok("Test9 (CLASS-03/T-08-01): zero-row and NULL input return typed frames, never stop()")

cat("\nALL PASS\n")
