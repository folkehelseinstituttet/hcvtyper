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

# --- Test 5: production read_csv/map_dfr typing path (locks CR-01 + CR-02) ----
# The in-memory tests above never exercise the readr type-inference that breaks
# the real summarize.R path (IN-01). Here we WRITE candidates + support to temp
# CSVs and READ them back exactly as summarize.R does — including a header-only
# (0-row) file mixed with a populated one in the map_dfr combine — then run the
# join. This reproduces CR-01 (numeric candidate_genotype) and CR-02 (header-only
# + populated bind_rows type clash) and asserts the join still succeeds.
tmp <- tempfile("asup_join_csv_"); dir.create(tmp)

# Candidate CSV written WITHOUT pinned types: genotype values are purely-digit
# ("3", "1"), so a naive read_csv infers <double> for candidate_genotype — the
# CR-01 trigger. The helper must coerce to character and still join at genotype.
cand_csv <- bind_rows(
  mk_cand("S5", 1, "1a_ACC", "1a", 8000, 95),
  mk_cand("S5", 2, "3b_ACC", "3b", 200, 8)
)
cand_path <- file.path(tmp, "S5.candidates.csv")
# The real *.candidates.csv carries `sample` (not `sampleName`); summarize.R adds
# sampleName via rename on read. mk_cand() carries both for the in-memory tests,
# so drop sampleName before writing to mirror the production CSV schema.
write_csv(cand_csv %>% select(-sampleName), cand_path)
# Round-trip exactly like summarize.R (col_types pinned there for CR-01/CR-02).
candidates_long5 <- map_dfr(cand_path, ~ read_csv(.x, col_types = cols(
  sample              = col_character(),
  candidate_rank      = col_integer(),
  candidate_ref       = col_character(),
  candidate_subtype   = col_character(),
  candidate_genotype  = col_character(),
  candidate_reads     = col_double(),
  candidate_cov       = col_double(),
  confirmation_status = col_character()
))) %>% rename(sampleName = sample)

# Two support CSVs: one populated, one HEADER-ONLY (the skip-assembly sample).
# Reading them with map_dfr is the CR-02 trigger when types are not pinned. The
# real *.assembly_support.csv carries `sample` (blast_parse.R writes it); mirror
# that schema so the rename(sampleName = sample) round-trips like production.
support_pop <- mk_support("S5", "3a", 2949, 97, 2800, 6) %>%
  rename(sample = sampleName)
support_pop_path <- file.path(tmp, "S5.assembly_support.csv")
write_csv(support_pop, support_pop_path)

support_empty_path <- file.path(tmp, "S6.assembly_support.csv")
# A genuinely header-only CSV: write a zero-row tibble with the six columns.
write_csv(support_pop[0, ], support_empty_path)

support_df5 <- map_dfr(c(support_pop_path, support_empty_path), ~ read_csv(.x, col_types = cols(
  sample                 = col_character(),
  subtype                = col_character(),
  best_contig_length     = col_double(),
  best_contig_pident     = col_double(),
  best_contig_aln_length = col_double(),
  best_contig_kmer_cov   = col_double()
))) %>% rename(sampleName = sample)

# The join must succeed at the DEFAULT genotype level despite the numeric-looking
# genotype column and the mixed header-only/populated support combine.
r5 <- join_assembly_support(candidates_long5, support_df5, match_level = "genotype")
if (nrow(r5) != 2) fail("Test5 row count must equal 2 candidates (no row loss through CSV path)")
cand2_row5 <- r5 %>% filter(candidate_rank == 2)
if (cand2_row5$assembly_support != "supported")
  fail("Test5 cand_2 (3b) must be supported by 3a at genotype level through the read_csv path")
if (cand2_row5$assembly_support_best_contig_length != 2949)
  fail("Test5 cand_2 must inherit the 3a 2949 contig length through the CSV round-trip")
ok("Test5 (CR-01/CR-02): join survives numeric genotype + header-only/populated map_dfr combine")

# --- Test 5b: helper coerces a <double> candidate_genotype itself (CR-01 belt) -
# Independently of summarize.R's pinned col_types, the helper must be robust to a
# caller that handed it a numeric candidate_genotype (the readr default for a
# purely-digit column). Read the SAME candidates CSV WITHOUT col_types so readr
# infers <double> for candidate_genotype, then join at the default genotype
# level. An un-coerced left_join would abort on <double> vs <character> keys.
candidates_double <- read_csv(cand_path, show_col_types = FALSE) %>%
  rename(sampleName = sample)
if (!is.numeric(candidates_double$candidate_genotype))
  fail("Test5b precondition: readr must infer numeric candidate_genotype here")
r5b <- join_assembly_support(candidates_double, support_df5, match_level = "genotype")
if (nrow(r5b) != 2) fail("Test5b row count must equal 2 candidates")
if ((r5b %>% filter(candidate_rank == 2) %>% pull(assembly_support)) != "supported")
  fail("Test5b helper must coerce numeric candidate_genotype and still match 3a at genotype level")
ok("Test5b (CR-01): helper coerces a numeric candidate_genotype key and joins cleanly")

# --- Test 6: match_level validation (WR-04) ---------------------------------
bad <- tryCatch({
  join_assembly_support(cand1, support1, match_level = "geneotype")  # typo
  FALSE
}, error = function(e) TRUE)
if (!bad) fail("Test6 an invalid match_level must error, not silently mean genotype")
ok("Test6 (WR-04): invalid match_level is rejected")

cat("\nALL PASS\n")
