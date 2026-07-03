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
                 denovo_min_contig_length = 500, denovo_min_kmer_cov = 2.0,
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

# --- Test 3b: dominance reads term uses targeted_reads_nodup, not candidate_reads ---
# Regression for the sim2 dominance-score inversion. In a real 2a:3a 70:30
# co-infection the neutral all-reference first-mapping (candidate_reads, WITH
# duplicates) MIS-RECRUITS reads and INVERTS the true abundance: it ranks the 30%
# minor (3a, 624520) above the 70% major (2a, 273886). The deduplicated targeted
# count (targeted_reads_nodup) is the truth: 2a 507708 > 3a 231132. The dominant
# must be selected on the targeted count, so 2a (the true major) wins.
# Values are the verbatim sim2 candidates.csv numbers.
sim2_inverted <- bind_rows(
  mk_cand("sim2inv", "3a_D17763", "3a", 624520, 100, 0.8114415381340713,
          sup_len = 9446, sup_kmer = 4184.394999, sup_pid = 100) %>%
    mutate(targeted_reads_nodup = 231132),
  mk_cand("sim2inv", "2a_D00944", "2a", 273886, 100, 0.8315985489939417,
          sup_len = 9695, sup_kmer = 9918.197952, sup_pid = 96.041) %>%
    mutate(targeted_reads_nodup = 507708)
)
r5b <- classify(sim2_inverted)
if (role_of(r5b, "2a_D00944") != "dominant")
  fail("Test3b: the true major 2a (higher targeted_reads_nodup) must be dominant, NOT the mis-recruited 3a candidate_reads leader")
if (role_of(r5b, "3a_D17763") != "co-infection")
  fail("Test3b: the true minor 3a must be co-infection, not dominant")
# Direct score check: scoring on targeted_reads_nodup must rank 2a above 3a.
sc5b <- score_candidates(sim2_inverted)
sc_2a <- sc5b %>% filter(candidate_ref == "2a_D00944") %>% pull(dominance_score)
sc_3a <- sc5b %>% filter(candidate_ref == "3a_D17763") %>% pull(dominance_score)
if (!(sc_2a > sc_3a))
  fail(sprintf("Test3b: 2a score (%.4f) must exceed 3a score (%.4f) when scoring on targeted_reads_nodup", sc_2a, sc_3a))
ok("Test3b: dominance reads term uses targeted_reads_nodup -> true major 2a wins despite inverted candidate_reads")

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

# --- Test 7 (D1/D3): concordance gate and fallback calls --------------------

# Scenario A: discordant candidate demoted regardless of abundance (D1).
# No concordance_status on 1a -> treated as eligible -> becomes dominant.
disc_a <- bind_rows(
  mk_cand("SA", "1a_ref", "1a", 200000, 90, 0.85) %>%
    mutate(concordance_status = "unconfirmed"),
  mk_cand("SA", "4g_ref", "4g",  80000, 70, 0.60) %>%
    mutate(concordance_status = "discordant")
)
ra <- classify(disc_a)
if (role_of(ra, "4g_ref") != "background")
  fail("Test7A discordant candidate must be background regardless of reads/cov")
if (reason_of(ra, "4g_ref") != "discordant_identity")
  fail("Test7A discordant reason must be discordant_identity")
if (role_of(ra, "1a_ref") != "dominant")
  fail("Test7A unconfirmed eligible 1a must become dominant (no floor gate)")

# Scenario B: all discordant -> overall_sample_call == "indeterminate" (D-14 new fallback).
disc_b <- bind_rows(
  mk_cand("SB", "3a_ref", "3a", 248,  8.8, 0.15) %>% mutate(concordance_status = "discordant"),
  mk_cand("SB", "1a_ref", "1a", 1030, 28,  0.40) %>% mutate(concordance_status = "discordant")
)
rb <- classify(disc_b)
if ((rb %>% pull(overall_sample_call) %>% unique()) != "indeterminate")
  fail("Test7B all-discordant sample must yield indeterminate")

# Scenario C: all candidates have cov == 0 -> "untypable" (D-14 new fallback).
no_cov <- bind_rows(
  mk_cand("SC", "3a_ref", "3a", 248,  0, 0),
  mk_cand("SC", "1a_ref", "1a", 1030, 0, 0)
)
rc <- classify(no_cov)
if ((rc %>% pull(overall_sample_call) %>% unique()) != "untypable")
  fail("Test7C zero-cov sample must yield untypable")

ok("Test7 (D1/D3): discordant->background/discordant_identity; all-discordant->indeterminate; no-cov->untypable")

# --- Test 12: D2 indeterminate dominance trigger (§7.6 Q2) ------------------

# Scenario A — trigger fires (ERR1810469-class).
# targeted reads favour 1a (1030 > 248); k-mer favours 3a (60x > 5.4x).
err469 <- bind_rows(
  mk_cand("E469", "1a_ref", "1a", 1069, 69, 0.78, sup_len = 4503, sup_kmer =  5.4, sup_pid = 93) %>%
    mutate(concordance_status = "confirmed", targeted_reads_nodup = 1030L),
  mk_cand("E469", "3a_ref", "3a", 5009, 46, 0.15, sup_len =  972, sup_kmer = 60.0, sup_pid = 95) %>%
    mutate(concordance_status = "confirmed", targeted_reads_nodup =  248L)
)
r12a <- classify(err469)
if (!all(r12a$role == "indeterminate"))
  fail("Test12A both candidates must be indeterminate when reads and k-mer rankings disagree")
if (!all(r12a$role_reason == "indeterminate_dominance_conflict"))
  fail("Test12A role_reason must be indeterminate_dominance_conflict for both")
if ((r12a %>% pull(overall_sample_call) %>% unique()) != "co-infection (indeterminate dominance)")
  fail("Test12A overall_sample_call must be 'co-infection (indeterminate dominance)'")

# Scenario B — trigger does NOT fire (ERR1810447-class: reads and k-mer agree).
# 1b leads on both reads (300000 >> 4199) and k-mer (45x >> 30x).
err447_12 <- bind_rows(
  mk_cand("E447", "1b_ref", "1b", 300000, 98, 0.92, sup_len = 8900, sup_kmer = 45, sup_pid = 99) %>%
    mutate(concordance_status = "confirmed", targeted_reads_nodup = 200000L),
  mk_cand("E447", "2b_ref", "2b",   4199, 90, 0.78, sup_len = 9207, sup_kmer = 30, sup_pid = 91) %>%
    mutate(concordance_status = "confirmed", targeted_reads_nodup =   3000L)
)
r12b <- classify(err447_12)
if (role_of(r12b, "1b_ref") != "dominant")
  fail("Test12B 1b must be dominant when reads and k-mer agree")
if (role_of(r12b, "2b_ref") != "co-infection")
  fail("Test12B 2b must be co-infection; trigger must NOT fire when rankings agree")
if ((r12b %>% pull(overall_sample_call) %>% unique()) != "co-infection")
  fail("Test12B overall_sample_call must be plain 'co-infection' when trigger is silent")

ok("Test12 (D2/§7.6 Q2): reads-vs-kmer disagreement -> indeterminate; agreement -> dominant/co-infection")

# --- Test 8: every candidate gets exactly one role (CLASS-01) ---------------
all_rows <- bind_rows(r1, r2, r3, r4, r5, r5b, r6, r7, r8, ra, rb, rc, r12a, r12b)
if (any(is.na(all_rows$role))) fail("Test8 every candidate must carry a non-NA role (CLASS-01)")
if (!all(all_rows$role %in% c("dominant", "co-infection", "background", "indeterminate")))
  fail("Test8 role must be one of dominant/co-infection/background/indeterminate")
if (!all(all_rows$overall_sample_call %in%
         c("monoinfection", "co-infection", "indeterminate", "untypable",
           "co-infection (indeterminate dominance)")))
  fail("Test8 overall_sample_call must be one of the five valid values")
ok("Test8 (CLASS-01/CLASS-04): every candidate has exactly one valid role + 5-value sample call")

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

# --- Test 10: 500 bp floor keeps ERR1810507-class 829bp minor contig (D5) -------
# Dominant 1a has a full contig. Minor 3a has a genuine but short contig: 829 bp /
# 30x / 92%. At 500: own_substantial=TRUE -> co-infection/corroborated.
# At 1000: own_substantial=FALSE (829 < 1000) AND dominant assembled -> refuted_denovo.
err507 <- bind_rows(
  mk_cand("ERR1810507", "1a_ref", "1a", 300000, 98, 0.92, sup_len = 9000, sup_kmer = 40, sup_pid = 99),
  mk_cand("ERR1810507", "3a_ref", "3a", 5000,   80, 0.70, sup_len = 829,  sup_kmer = 30, sup_pid = 92)
)
r_err507 <- classify(err507)   # uses the 500-floor classify() helper
if (role_of(r_err507, "3a_ref") != "co-infection")
  fail("Test10: ERR1810507 829bp 3a contig must be co-infection at the 500bp floor")
if (reason_of(r_err507, "3a_ref") != "corroborated")
  fail("Test10: ERR1810507 829bp 3a contig reason must be corroborated")
# Precondition: at the old 1000bp floor the same contig would be refuted.
r_err507_strict <- classify_roles(score_candidates(err507), minRead = 500, minCov = 30,
                                  denovo_min_contig_length = 1000, denovo_min_kmer_cov = 2.0,
                                  denovo_min_blast_identity = 90, match_level = "genotype")
if (role_of(r_err507_strict, "3a_ref") != "background")
  fail("Test10 precondition: at the old 1000bp floor the 829bp contig should be refuted")
ok("Test10 (D5): 829bp minor contig -> co-infection/corroborated at 500bp; refuted_denovo at old 1000bp floor")

# --- Test 11: apply_concordance() three-way check (D8) -----------------------
suppressPackageStartupMessages(library(tidyverse))
source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))

mk_conc <- function(sample, subtype, glue_gt, supp, supp_subtype) {
  tibble(
    sampleName              = sample,
    candidate_ref           = paste0(subtype, "_ref"),
    candidate_subtype       = subtype,
    candidate_genotype      = genotype_from_subtype(subtype),
    candidate_glue_genotype = glue_gt,
    assembly_support        = supp,
    assembly_support_subtype = supp_subtype
  )
}

conc_df <- bind_rows(
  # confirmed: mapping 1a + GLUE gt1 + de novo 1a -> confirmed
  mk_conc("S1", "1a", "1", "supported", "1a"),
  # discordant (mapping vs GLUE): mapping 4g + GLUE gt1 -> discordant
  mk_conc("S2", "4g", "1", "none",       NA),
  # discordant (mapping vs de novo): mapping 1a + de novo 3a -> discordant
  mk_conc("S3", "1a", NA,  "supported", "3a"),
  # unconfirmed: mapping 1a + no GLUE + no de novo
  mk_conc("S4", "1a", NA,  "none",       NA),
  # unconfirmed: mapping 1a + GLUE agrees + no de novo
  mk_conc("S5", "1a", "1", "none",       NA)
)
conc_out <- apply_concordance(conc_df)

status_of <- function(r, s) r %>% filter(sampleName == s) %>% pull(concordance_status)
reason_of_c <- function(r, s) r %>% filter(sampleName == s) %>% pull(concordance_reason)

if (status_of(conc_out, "S1") != "confirmed")
  fail(sprintf("Test11: S1 (all legs agree) must be confirmed, got '%s'", status_of(conc_out, "S1")))
if (status_of(conc_out, "S2") != "discordant")
  fail(sprintf("Test11: S2 (4g vs GLUE gt1) must be discordant, got '%s'", status_of(conc_out, "S2")))
if (reason_of_c(conc_out, "S2") != "discordant_mapping_vs_glue")
  fail(sprintf("Test11: S2 reason must be discordant_mapping_vs_glue, got '%s'", reason_of_c(conc_out, "S2")))
if (status_of(conc_out, "S3") != "discordant")
  fail(sprintf("Test11: S3 (mapping 1a vs de novo 3a) must be discordant, got '%s'", status_of(conc_out, "S3")))
if (reason_of_c(conc_out, "S3") != "discordant_mapping_vs_denovo")
  fail(sprintf("Test11: S3 reason must be discordant_mapping_vs_denovo, got '%s'", reason_of_c(conc_out, "S3")))
if (status_of(conc_out, "S4") != "unconfirmed")
  fail(sprintf("Test11: S4 (no corroborating legs) must be unconfirmed, got '%s'", status_of(conc_out, "S4")))
if (reason_of_c(conc_out, "S4") != "no_corroborating_legs")
  fail(sprintf("Test11: S4 reason must be no_corroborating_legs, got '%s'", reason_of_c(conc_out, "S4")))
if (status_of(conc_out, "S5") != "unconfirmed")
  fail(sprintf("Test11: S5 (GLUE only, no de novo) must be unconfirmed, got '%s'", status_of(conc_out, "S5")))
if (reason_of_c(conc_out, "S5") != "two_legs_glue_only")
  fail(sprintf("Test11: S5 reason must be two_legs_glue_only, got '%s'", reason_of_c(conc_out, "S5")))

# NULL input must not abort.
conc_null <- apply_concordance(NULL)
if (nrow(conc_null) != 0) fail("Test11: NULL input must yield zero-row frame")
ok("Test11 (D8): apply_concordance() correctly classifies confirmed/unconfirmed/discordant cases")

# --- Test 13: apply_concordance() 2k1b GLUE exemption (COMPAT-04) ------------
# HCV-GLUE has no CRF_02k/1b clade, so a genuine 2k/1b recombinant is always
# reported by GLUE as genotype 1 or 2, never 2k1b. apply_concordance() must NOT
# flag that as discordant. Reuses mk_conc()/status_of()/reason_of_c() from Test 11.
# NOTE: labelled Test 13 (not 12) — the file already reuses the "Test 12" label
# for the D2 indeterminate-dominance test above.
conc_2k1b <- bind_rows(
  # Exact shape of the three failing TEST-run samples: mapping 2k1b + GLUE gt1 +
  # de novo 2k1b -> confirmed / confirmed_2k1b_recombinant.
  mk_conc("S6", "2k1b", "1", "supported", "2k1b"),
  # GLUE gt2 still exempt; no de novo leg -> unconfirmed / two_legs_2k1b_recombinant_glue_only.
  mk_conc("S7", "2k1b", "2", "none",       NA),
  # GLUE gt3 is NOT in the exemption list -> still discordant (proves narrow scope).
  mk_conc("S8", "2k1b", "3", "supported", "2k1b")
)
conc_2k1b_out <- apply_concordance(conc_2k1b)

if (status_of(conc_2k1b_out, "S6") != "confirmed")
  fail(sprintf("Test13: S6 (2k1b + GLUE gt1 + de novo 2k1b) must be confirmed, got '%s'", status_of(conc_2k1b_out, "S6")))
if (reason_of_c(conc_2k1b_out, "S6") != "confirmed_2k1b_recombinant")
  fail(sprintf("Test13: S6 reason must be confirmed_2k1b_recombinant, got '%s'", reason_of_c(conc_2k1b_out, "S6")))
if (status_of(conc_2k1b_out, "S7") != "unconfirmed")
  fail(sprintf("Test13: S7 (2k1b + GLUE gt2, no de novo) must be unconfirmed, got '%s'", status_of(conc_2k1b_out, "S7")))
if (reason_of_c(conc_2k1b_out, "S7") != "two_legs_2k1b_recombinant_glue_only")
  fail(sprintf("Test13: S7 reason must be two_legs_2k1b_recombinant_glue_only, got '%s'", reason_of_c(conc_2k1b_out, "S7")))
if (status_of(conc_2k1b_out, "S8") != "discordant")
  fail(sprintf("Test13: S8 (2k1b vs GLUE gt3) must STILL be discordant, got '%s'", status_of(conc_2k1b_out, "S8")))
if (reason_of_c(conc_2k1b_out, "S8") != "discordant_mapping_vs_glue")
  fail(sprintf("Test13: S8 reason must be discordant_mapping_vs_glue, got '%s'", reason_of_c(conc_2k1b_out, "S8")))
ok("Test13 (COMPAT-04): apply_concordance() 2k1b GLUE exemption — S6 confirmed, S7 unconfirmed-glue-only, S8 gt3 still discordant")

# --- Test 14: end-to-end 2k1b -> dominant/monoinfection ----------------------
# The regression that shipped once: a lone 2k1b candidate (mapping 2k1b + GLUE
# gt1 + de novo 2k1b) was flagged discordant, zeroing its eligibility and
# yielding overall_sample_call = indeterminate / role = background. Prove the
# full apply_concordance() -> classify() chain now admits it as dominant.
# mk_cand() does not emit the concordance legs, so add them via mutate(); and
# classify_roles() does not call apply_concordance() internally (Test 7 relies on
# concordance_status being pre-set), so run apply_concordance() FIRST.
cand_2k1b <- mk_cand("S14", "2k1b_ref", "2k1b", 200000, 98, 0.90,
                     sup_len = 9000, sup_kmer = 40, sup_pid = 99) %>%
  mutate(candidate_glue_genotype  = "1",
         assembly_support         = "supported",
         assembly_support_subtype = "2k1b")
r14 <- classify(apply_concordance(cand_2k1b))
if (role_of(r14, "2k1b_ref") != "dominant")
  fail(sprintf("Test14: lone 2k1b candidate must become dominant, got '%s'", role_of(r14, "2k1b_ref")))
if ((r14 %>% pull(overall_sample_call) %>% unique()) != "monoinfection")
  fail(sprintf("Test14: lone 2k1b sample must be monoinfection, got '%s'",
               paste(r14 %>% pull(overall_sample_call) %>% unique(), collapse = ",")))
ok("Test14 (COMPAT-04 e2e): lone 2k1b candidate reaches dominant / monoinfection through apply_concordance() -> classify()")

cat("\nALL PASS\n")
