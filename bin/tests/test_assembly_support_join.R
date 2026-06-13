#!/usr/bin/env Rscript

# test_assembly_support_join.R --------------------------------------------
# Function-level tests for the Phase-7 (ASUP-02) genotype-level assembly-support
# join, bin/assembly_support_join.R::join_assembly_support(). Modelled on
# test_summarize_denovo.R: source the real helper, build small in-memory
# tibbles, and assert on the REAL function's output (no inline re-implementation).
#
# Asserted behaviours (07-02-PLAN criteria #2/#3/#4):
#   1. Genotype match (criterion #2): a 3b candidate is corroborated by 3a
#      support at match_level="genotype" (both genotype "3"), inheriting the 3a
#      contig's metrics (length 2949).
#   2. Subtype non-match: the same inputs at match_level="subtype" leave the 3b
#      candidate at assembly_support="none" (no 3b support row).
#   3. No-row-loss (criterion #3): a candidate whose match key has no support row
#      keeps its row with assembly_support="none" + NA metrics; output row count
#      == input candidate row count.
#   4. Criterion #4 reproduction: on a fixture mirroring a regression sample at
#      default flags + N=2, the cand_2 (minor-slot) assembly_support_best_contig_length
#      equals the best DIFFERENT-genotype contig length that today's blast_parse.R
#      §7 minor logic would put in denovo_minor_contig_length.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "assembly_support_join.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Builder for a Phase-6 long-format candidate row (8 columns, precomputed keys).
mk_cand <- function(sample, rank, ref, subtype, reads, cov, status = "pass") {
  tibble(
    sample              = sample,
    sampleName          = sample,
    candidate_rank      = rank,
    candidate_ref       = ref,
    candidate_subtype   = subtype,
    candidate_genotype  = genotype_from_subtype(subtype),
    candidate_reads     = reads,
    candidate_cov       = cov,
    confirmation_status = status
  )
}

# Builder for a per-subtype assembly-support row (the six Plan-01 columns +
# the sampleName key that summarize.R attaches on read).
mk_support <- function(sample, subtype, length, pident, aln_length, kmer_cov) {
  tibble(
    sampleName             = sample,
    subtype                = subtype,
    best_contig_length     = length,
    best_contig_pident     = pident,
    best_contig_aln_length = aln_length,
    best_contig_kmer_cov   = kmer_cov
  )
}

# --- Test 1: genotype match (criterion #2) — 3b candidate <- 3a support -----
cand1 <- bind_rows(
  mk_cand("S1", 1, "1a_ACC", "1a", 8000, 95),
  mk_cand("S1", 2, "3b_ACC", "3b", 200, 8)
)
support1 <- bind_rows(
  mk_support("S1", "1a", 5000, 99, 4900, 30),
  mk_support("S1", "3a", 2949, 97, 2800, 6)
)
r1 <- join_assembly_support(cand1, support1, match_level = "genotype")
if (nrow(r1) != 2) fail("Test1 row count must equal 2 candidates")
cand2_row <- r1 %>% filter(candidate_rank == 2)
if (cand2_row$assembly_support != "supported") fail("Test1 cand_2 (3b) must be supported by 3a at genotype level")
if (cand2_row$assembly_support_best_contig_length != 2949) fail("Test1 cand_2 must inherit 3a contig length 2949")
if (cand2_row$assembly_support_subtype != "3a") fail("Test1 cand_2 support subtype must be 3a")
cand1_row <- r1 %>% filter(candidate_rank == 1)
if (cand1_row$assembly_support_best_contig_length != 5000) fail("Test1 cand_1 (1a) must inherit 1a contig length 5000")
ok("Test1 (criterion #2): 3b candidate corroborated by 3a support at genotype level")

# --- Test 2: subtype non-match -> cand_2 stays "none" -----------------------
r2 <- join_assembly_support(cand1, support1, match_level = "subtype")
cand2_row2 <- r2 %>% filter(candidate_rank == 2)
if (cand2_row2$assembly_support != "none") fail("Test2 cand_2 (3b) must be none at subtype level (no 3b support)")
if (!is.na(cand2_row2$assembly_support_best_contig_length)) fail("Test2 cand_2 metrics must be NA at subtype level")
# cand_1 (1a) still matches a 1a support row at subtype level
cand1_row2 <- r2 %>% filter(candidate_rank == 1)
if (cand1_row2$assembly_support != "supported") fail("Test2 cand_1 (1a) must still match 1a support at subtype level")
ok("Test2: subtype-level non-match leaves 3b candidate as none")

# --- Test 3: no-row-loss + NA-fill (criterion #3) ---------------------------
cand3 <- bind_rows(
  mk_cand("S2", 1, "2a_ACC", "2a", 6000, 90),
  mk_cand("S2", 2, "5a_ACC", "5a", 100, 3)
)
# Support only for genotype 2; genotype 5 has NO support row.
support3 <- mk_support("S2", "2a", 4000, 98, 3900, 20)
r3 <- join_assembly_support(cand3, support3, match_level = "genotype")
if (nrow(r3) != nrow(cand3)) fail("Test3 output row count must equal input candidate row count")
unsupported <- r3 %>% filter(candidate_rank == 2)
if (unsupported$assembly_support != "none") fail("Test3 unmatched candidate must be assembly_support='none'")
if (!is.na(unsupported$assembly_support_best_contig_length)) fail("Test3 unmatched candidate metrics must be NA")
if (!is.na(unsupported$assembly_support_subtype)) fail("Test3 unmatched candidate subtype must be NA")
ok("Test3 (criterion #3): no row loss, explicit none + NA metrics for unsupported candidate")

# --- Test 3b: empty support_df (skip-assembly) -> all none, no abort --------
r3b <- join_assembly_support(cand3, support3[0, ], match_level = "genotype")
if (nrow(r3b) != nrow(cand3)) fail("Test3b empty-support row count must equal candidates")
if (!all(r3b$assembly_support == "none")) fail("Test3b empty support must yield all none")
ok("Test3b: empty support_df NA-fills every candidate, no abort")

# --- Test 3c: empty candidates_df -> typed zero-row frame, no abort ---------
r3c <- join_assembly_support(cand3[0, ], support3, match_level = "genotype")
if (nrow(r3c) != 0) fail("Test3c empty candidates must yield zero rows")
if (!"assembly_support" %in% names(r3c)) fail("Test3c must still carry the assembly_support column")
ok("Test3c: empty candidates_df returns typed zero-row frame")

# --- Test 4: criterion #4 reproduction (cand_2 == denovo_minor) -------------
# Mirror a regression sample: cand_1 maps to genotype 1 (major), cand_2 to a
# different genotype (minor slot). Today's blast_parse.R §7 selects the longest
# contig of a DIFFERENT genotype than the major as the minor evidence. Build a
# support_df whose best different-genotype contig is the 2b 2949 contig (the
# ERR1810453 anchor) plus shorter noise contigs; the cand_2 (genotype 2)
# assembly_support_best_contig_length must equal that 2949 length.
cand4 <- bind_rows(
  mk_cand("ERR", 1, "1a_M62321", "1a", 50000, 99),
  mk_cand("ERR", 2, "2b_ACC",    "2b", 300, 4)
)
# Best 2b/genotype-2 contig is 2949 (the genuine partial), with a 300bp noise 2b
# hit that must NOT win the collapse (single-best-by-length, D-03).
support4 <- bind_rows(
  mk_support("ERR", "1a", 9000, 99, 8900, 50),  # major-genotype contig
  mk_support("ERR", "2b", 2949, 97, 2800, 5)    # the minor-genotype best contig (denovo_minor)
)
# denovo_minor_contig_length today = longest different-genotype-than-major contig.
denovo_minor_contig_length <- support4 %>%
  mutate(geno = genotype_from_subtype(subtype)) %>%
  filter(geno != "1") %>%
  slice_max(best_contig_length, n = 1, with_ties = FALSE) %>%
  pull(best_contig_length)
r4 <- join_assembly_support(cand4, support4, match_level = "genotype")
cand2_len <- r4 %>% filter(candidate_rank == 2) %>% pull(assembly_support_best_contig_length)
if (length(cand2_len) != 1) fail("Test4 expected exactly one cand_2 row")
if (cand2_len != denovo_minor_contig_length) {
  fail(sprintf("Test4 cand_2 length %s must equal denovo_minor_contig_length %s",
               cand2_len, denovo_minor_contig_length))
}
if (cand2_len != 2949) fail("Test4 cand_2 length must be the 2949 best 2b contig, not the 300bp noise")
ok("Test4 (criterion #4): cand_2 assembly support equals today's denovo_minor_contig_length")

cat("\nALL PASS\n")
