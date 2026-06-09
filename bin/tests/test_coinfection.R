#!/usr/bin/env Rscript

# test_coinfection.R ------------------------------------------------------
# Scenario-specific co-infection NON-SUPPRESSION test (D-07 / D-08 / D-09).
#
# Every assertion calls the REAL pure classify_minor_denovo() (bin/denovo_confirm.R)
# — NEVER an inline re-implementation — with raw subtypes routed through the
# canonical genotype_from_subtype() helper for the major/minor genotype args.
# Each named sample is built as a synthetic per-contig BLAST frame using the
# 8-column schema classify_minor_denovo reads (qseqid, subtype, pident, length,
# evalue, bitscore, sc_length, kmer_cov) and asserted against its EXACT verified
# outcome from 04-RESEARCH §"Co-Infection Fixture Outcomes".
#
# D-09 INTEGRITY RULE (the single most important constraint in this phase):
# the worst-case extreme-ratio IVT mixture GENUINELY returns "refuted" under the
# shipped asymmetric refute rule. This test asserts that REAL behaviour and its
# existence DOCUMENTS the known limitation. It does NOT lower thresholds, does NOT
# fabricate a convenient substantial-minor contig, and contains NO blanket
# "never refuted" assertion. Papering over that refute would be a correctness lie.
#
# Run from any cwd:
#   Rscript bin/tests/test_coinfection.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd
# (verbatim self-location pattern from test_denovo_confirm.R:17-20).
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
bin_dir <- if (length(this_file) == 1) dirname(dirname(normalizePath(this_file))) else "bin"

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "denovo_confirm.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Build a per-contig frame (8-col schema). subtype/sc_length/kmer_cov/pident are
# the discriminating columns; length/evalue/bitscore are filler the top-hit
# selection (arrange(evalue, desc(bitscore))) uses but the substantiality floor
# ignores (substantiality is on sc_length, NOT alignment length).
mk <- function(subtype, sc_length, kmer_cov, pident, qseqid = NULL) {
  if (is.null(qseqid)) qseqid <- paste0("NODE_", subtype, "_", sc_length)
  tibble(
    qseqid    = qseqid,
    subtype   = subtype,
    pident    = pident,
    length    = sc_length,
    evalue    = 0,
    bitscore  = 1000,
    sc_length = sc_length,
    kmer_cov  = kmer_cov
  )
}

# Assert one named scenario. major_sub/minor_sub are RAW subtypes routed through
# genotype_from_subtype() — we never hardcode "3a" -> "3".
assert_scenario <- function(name, contigs, major_sub, minor_sub, expect) {
  major_geno <- genotype_from_subtype(major_sub)
  minor_geno <- if (is.na(minor_sub)) NA_character_ else genotype_from_subtype(minor_sub)
  got <- classify_minor_denovo(contigs, major_geno, minor_geno,
                               min_len = 1000, min_kmer = 2.0, min_pid = 90,
                               match_level = "genotype")
  if (!identical(got, expect)) {
    fail(sprintf("%s: expected %s, got %s", name, expect, got))
  }
  ok(sprintf("%s -> %s", name, expect))
}

# --- Genuine co-infections that MUST be confirmed -------------------------

# sim1 (1a:1b). NOTE: 1a and 1b are the SAME genotype 1, so minor_geno == "1".
# has_minor_sub is satisfied at GENOTYPE level by the 1a MAJOR contig (both
# contigs key to genotype "1"). This is correct for the genuine 1a/1b
# co-infection — it is a genotype-level match, NOT a separate minor-contig match
# (Pitfall 4). The shipped genotype-level rule deliberately treats this as
# confirmed.
assert_scenario("sim1 (1a:1b)",
  bind_rows(mk("1b", 9339, 100, 95), mk("1a", 9076, 100, 95)),
  "1a", "1b", "confirmed_by_denovo")

# sim2 (2a:3a -> geno 2:3): both substantial, distinct genotypes.
assert_scenario("sim2 (2a:3a)",
  bind_rows(mk("2a", 9500, 100, 95), mk("3a", 9500, 100, 95)),
  "2a", "3a", "confirmed_by_denovo")

# ERR1810505 (1a:3a -> 1:3): both substantial.
assert_scenario("ERR1810505 (1a:3a)",
  bind_rows(mk("1a", 9100, 100, 95), mk("3a", 9100, 100, 95)),
  "1a", "3a", "confirmed_by_denovo")

# ERR1810511 (1a:2a -> 1:2): both substantial (kmer_cov 50).
assert_scenario("ERR1810511 (1a:2a)",
  bind_rows(mk("1a", 8177, 50, 95), mk("2a", 7517, 50, 95)),
  "1a", "2a", "confirmed_by_denovo")

# ERR1810453 (1a:2b -> 1:2): the PARTIAL-CONTIG calibration anchor. The 2b minor
# is only 2,949 bp / 5x but MUST confirm (>=1000 bp, >=2.0x, >=90% pid all pass).
assert_scenario("ERR1810453 (1a:2b partial-contig anchor)",
  bind_rows(mk("1a", 9000, 100, 95), mk("2b", 2949, 5, 95)),
  "1a", "2b", "confirmed_by_denovo")

# Qiu SRR1762352-55 (1b:3a -> 1:3): both substantial.
assert_scenario("Qiu SRR1762352-55 (1b:3a)",
  bind_rows(mk("1b", 9300, 120, 95), mk("3a", 9400, 90, 95)),
  "1b", "3a", "confirmed_by_denovo")

# ivt_minor_assembles (1a:3a -> 1:3): minor barely assembles (1200 bp / 3.0x) but
# clears all floors -> protected, confirmed.
assert_scenario("ivt_minor_assembles (1a:3a, minor barely assembles)",
  bind_rows(mk("1a", 9000, 5000, 95), mk("3a", 1200, 3.0, 95)),
  "1a", "3a", "confirmed_by_denovo")

# --- Fall-back / refute branches ------------------------------------------

# ivt_denovo_failed_overall (1a:3a -> 1:3): nothing substantial (both contigs
# tiny / low cov) -> de novo failed overall -> unconfirmed (NEVER refute), the
# fall-back protects the minor.
assert_scenario("ivt_denovo_failed_overall (1a:3a, de novo failed)",
  bind_rows(mk("1a", 400, 1.5, 95), mk("3a", 300, 1.0, 95)),
  "1a", "3a", "unconfirmed")

# sim1_1a_single (1a:4g -> 1:4): the HEADLINE refute artefact. A huge substantial
# 1a major contig and NO 4g contig at all -> refuted (the cross-mapping artefact
# the whole feature exists to suppress).
assert_scenario("sim1_1a_single (1a:4g, headline refute artefact)",
  mk("1a", 9189, 15056, 99),
  "1a", "4g", "refuted")

# ivt_extreme_ratio_minor_unassembled (1a:3a -> 1:3): KNOWN LIMITATION (D-09).
# The genuine minor exists biologically but at an extreme ratio it assembles only
# as a 300 bp / 1.0x noise contig, below every floor, while the major is huge and
# substantial. The shipped asymmetric refute rule therefore returns "refuted",
# SUPPRESSING a genuine extreme-ratio minor. This test asserts the REAL refuted
# outcome on purpose: its existence DOCUMENTS the gap. Do NOT lower thresholds,
# do NOT fabricate a substantial 3a contig, do NOT weaken this to a blanket
# "never refuted" — that would paper over the genuine limitation.
assert_scenario("ivt_extreme_ratio_minor_unassembled (1a:3a) [KNOWN LIMITATION, D-09]",
  bind_rows(mk("1a", 9000, 5000, 95), mk("3a", 300, 1.0, 95)),
  "1a", "3a", "refuted")

cat("\nALL PASS\n")
