#!/usr/bin/env Rscript

# test_assembly_support.R ------------------------------------------------
# Subprocess-contract test (Phase 7, ASUP-01) for the REAL neutral per-subtype
# assembly-support roll-up in bin/blast_parse.R (§4b). Modelled on
# test_candidate_selection.R: we invoke the script via system2("Rscript", ...)
# on small synthetic blast_out / contigs / references fixtures built inline in a
# tempdir, and assert on what it writes to <prefix>.assembly_support.csv.
#
# Behaviours asserted (07-CONTEXT D-01/D-02/D-03, 07-PATTERNS):
#   Case A (single-best-by-sc_length, D-03): blast_out with 3a contigs of length
#     2949 and 300 plus a 2b contig of length 1500 -> exactly 2 rows; the 3a row's
#     best_contig_length == 2949 (NOT 300), and its pident / aln_length / kmer_cov
#     come from the 2949 contig (one coherent contig backs all four metrics).
#   Case B (multi-hit de-dup): a single contig with two BLAST hits to the same
#     subtype reference -> one row for that subtype, length not duplicated.
#   Case C (empty / no hits): empty blast_out -> a header-only assembly_support.csv
#     is written, the script exits 0, and the header equals the seven-column contract
#     (T-07-01 DoS guard: never abort on skip-assembly).
#
# blast_parse.R sources nothing test-local; it reads the references FASTA (must be
# non-empty or it aborts) and the contigs FASTA. Each case runs in its own tempdir
# with absolute fixture paths passed as args. The script is located relative to
# this test file via --file= (the staged-path contract) — never an absolute hardcode.
#
# RED until Task 1 lands the §4b roll-up (the script currently writes no
# *.assembly_support.csv).
#
# Run from any cwd:
#   Rscript bin/tests/test_assembly_support.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "blast_parse.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Seven-column contract: `best_ref` (the winning contig's BLAST subject / reference)
# was added after `subtype` mid-project so the rescue step (rescue_evaluation.R) can
# rank alternate subtypes by their best_ref; blast_parse.R §4b emits it in this order.
SUPPORT_COLS <- c("sample", "subtype", "best_ref", "best_contig_length",
                  "best_contig_pident", "best_contig_aln_length", "best_contig_kmer_cov")

# The 12 standard outfmt6 columns. `hits` is a list of named lists each carrying
# the fields we vary (qseqid, sseqid, pident, length); the rest are filled with
# plausible constants. blast_out is tab-separated, no header (blast_parse.R reads
# col_names = FALSE). Returns the assembly-support tibble read back (or NULL) plus
# the subprocess exit status.
run_blast_parse <- function(case, prefix, hits) {
  wd <- tempfile(paste0("asup_", case, "_")); dir.create(wd)

  # blast_out outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend
  #                    sstart send evalue bitscore
  blast_rows <- map_chr(hits, function(h) {
    paste(h$qseqid, h$sseqid, h$pident, h$length,
          0, 0, 1, h$length, 1, h$length, "1e-100", h$length,
          sep = "\t")
  })
  blast_path <- file.path(wd, paste0(case, "_blast_out.tsv"))
  writeLines(blast_rows, blast_path)

  # contigs FASTA: one record per distinct qseqid (names must match qseqids).
  qseqids <- unique(map_chr(hits, "qseqid"))
  contig_lines <- unlist(lapply(qseqids, function(q) c(paste0(">", q), strrep("A", 60))))
  if (is.null(contig_lines)) contig_lines <- c(">dummy_NODE", strrep("A", 60))
  contigs_path <- file.path(wd, paste0(case, "_contigs.fa"))
  writeLines(contig_lines, contigs_path)

  # references FASTA: one record per distinct sseqid (must be non-empty or
  # blast_parse.R aborts on the reference read). Keep at least one record.
  sseqids <- unique(map_chr(hits, "sseqid"))
  if (length(sseqids) == 0) sseqids <- "3a_DUMMY"
  ref_lines <- unlist(lapply(sseqids, function(s) c(paste0(">", s), strrep("C", 60))))
  refs_path <- file.path(wd, paste0(case, "_refs.fa"))
  writeLines(ref_lines, refs_path)

  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(prefix), shQuote(blast_path),
      shQuote(contigs_path), shQuote(refs_path), "HCV"),
    stdout = FALSE, stderr = FALSE
  )

  support_path <- file.path(wd, paste0(prefix, ".assembly_support.csv"))
  support <- if (file.exists(support_path)) {
    read_csv(support_path, show_col_types = FALSE)
  } else {
    NULL
  }
  # 260803-ogc: blastparse.csv is read back too, so the minor-slot coherence
  # contract (ref / contig / length from ONE row) can be asserted directly.
  bp_path <- file.path(wd, paste0(prefix, ".blastparse.csv"))
  blastparse <- if (file.exists(bp_path)) read_csv(bp_path, show_col_types = FALSE) else NULL
  # Also read the raw header line so the empty case can assert the contract even
  # when read_csv yields a zero-row frame.
  header <- if (file.exists(support_path)) readLines(support_path, n = 1) else NA_character_
  list(exit = exit, support = support, header = header, path = support_path,
       blastparse = blastparse)
}

# Helper to build an outfmt6 hit. qseqid carries NODE_<n>_length_<len>_cov_<cov>
# so blast_parse.R extracts sc_length and kmer_cov from the header.
node <- function(n, len, cov, sseqid, pident, aln_length) {
  list(
    qseqid = paste0("NODE_", n, "_length_", len, "_cov_", cov),
    sseqid = sseqid, pident = pident, length = aln_length
  )
}

# --- Case A: single-best-by-sc_length across multi-subtype (D-03) ------------
# 3a has two contigs (2949 and 300); 2b has one (1500). Expect 2 rows; the 3a
# row carries the 2949 contig's metrics, NOT the 300 one.
r_a <- run_blast_parse(
  "best",
  "BEST",
  hits = list(
    node(1, 2949, 12.0, "3a_REF1", 99.0, 2900),
    node(2,  300,  4.0, "3a_REF1", 95.0,  290),
    node(3, 1500,  8.0, "2b_REF1", 97.0, 1480)
  )
)
if (is.null(r_a$support)) fail("best: no assembly_support.csv written")
if (!identical(colnames(r_a$support), SUPPORT_COLS))
  fail(paste("best: header must be the seven-column contract; got:",
             paste(colnames(r_a$support), collapse = ",")))
if (nrow(r_a$support) != 2)
  fail(paste("best: expected exactly 2 rows (one per subtype), got", nrow(r_a$support)))
row3a <- r_a$support %>% filter(subtype == "3a")
row2b <- r_a$support %>% filter(subtype == "2b")
if (nrow(row3a) != 1) fail("best: expected exactly one 3a row")
if (row3a$best_contig_length[1] != 2949)
  fail(paste("best: 3a best_contig_length must be 2949 (single-best-by-sc_length), got",
             row3a$best_contig_length[1]))
# All four metrics must come from the SAME (2949) contig, not the 300 one.
if (row3a$best_contig_pident[1] != 99.0)
  fail(paste("best: 3a pident must come from the 2949 contig (99.0), got",
             row3a$best_contig_pident[1]))
if (row3a$best_contig_aln_length[1] != 2900)
  fail(paste("best: 3a aln_length must come from the 2949 contig (2900), got",
             row3a$best_contig_aln_length[1]))
if (row3a$best_contig_kmer_cov[1] != 12.0)
  fail(paste("best: 3a kmer_cov must come from the 2949 contig (12.0), got",
             row3a$best_contig_kmer_cov[1]))
if (nrow(row2b) != 1 || row2b$best_contig_length[1] != 1500)
  fail("best: 2b row must carry its 1500-length contig")
if (any(r_a$support$sample != "BEST"))
  fail("best: sample column must equal the prefix 'BEST'")
ok("best -> 2 rows; 3a carries the 2949 contig's four metrics (D-03 single-best-by-length)")

# --- Case B: multi-hit de-dup (one contig, two hits to same subtype ref) -----
# A single 2b contig (length 1500) has TWO BLAST hits to the same reference.
# Expect exactly one row for 2b; the length is NOT duplicated.
r_b <- run_blast_parse(
  "dedup",
  "DEDUP",
  hits = list(
    node(1, 1500, 8.0, "2b_REF1", 97.0, 1480),
    node(1, 1500, 8.0, "2b_REF1", 96.0,  700)
  )
)
if (is.null(r_b$support)) fail("dedup: no assembly_support.csv written")
b2b <- r_b$support %>% filter(subtype == "2b")
if (nrow(b2b) != 1)
  fail(paste("dedup: a single contig with two hits must yield ONE 2b row, got",
             nrow(b2b)))
if (b2b$best_contig_length[1] != 1500)
  fail(paste("dedup: 2b best_contig_length must be 1500 (not duplicated), got",
             b2b$best_contig_length[1]))
ok("dedup -> one 2b row, length not duplicated across the two same-ref hits")

# --- Case C: empty / no hits -> header-only CSV, exit 0, seven-column contract --
r_c <- run_blast_parse("empty", "EMPTY", hits = list())
if (r_c$exit != 0)
  fail(paste("empty: script must exit 0 on no-hit input, got exit", r_c$exit))
if (is.na(r_c$header))
  fail("empty: a header-only assembly_support.csv must still be written")
got_cols <- strsplit(r_c$header, ",", fixed = TRUE)[[1]]
if (!identical(got_cols, SUPPORT_COLS))
  fail(paste("empty: header-only CSV must carry the seven-column contract; got:",
             r_c$header))
if (!is.null(r_c$support) && nrow(r_c$support) != 0)
  fail(paste("empty: header-only CSV must have zero data rows, got",
             nrow(r_c$support)))
ok("empty -> header-only CSV, exit 0, seven-column contract (T-07-01 DoS guard)")

# --- Case D: minor_ref / minor_contig / minor_contig_length name ONE contig ----
# Regression for the 2633901 defect (260803-ogc). Shape reproduced exactly:
#   NODE_1  5159 bp, top hit 1a, aln 5160  -> the major (best overall hit)
#   NODE_2  3232 bp, top hit 1a, aln 3088  -> a 1a contig that ALSO hits the 6i
#                                             reference at aln 879 (5'UTR/core)
#   NODE_3  1620 bp, only hit  6i, aln 69  -> the actual off-genotype contig
# bitscore == aln length in this harness, so NODE_2's 879 bp hit to 6i outscores
# NODE_3's 69 bp hit to the same reference by an order of magnitude. That is what
# used to make summarize.R name NODE_2 while the reference and the length described
# NODE_3.
r_d <- run_blast_parse("minorcoherence", "MINCOH", hits = list(
  node(1, 5159, 9828.93, "1a_HQ850279", 93.9, 5160),
  node(2, 3232, 6285.47, "1a_HQ850279", 93.0, 3088),
  node(2, 3232, 6285.47, "6i_DQ835770", 88.7,  879),
  node(3, 1620,    1.02, "6i_DQ835770", 91.3,   69)
))
if (r_d$exit != 0) fail(paste("minor-coherence: exit", r_d$exit))
bp <- r_d$blastparse
if (is.null(bp)) fail("minor-coherence: blastparse.csv must be written")
if (!"minor_contig" %in% names(bp))
  fail(paste("minor-coherence: blastparse.csv must carry minor_contig; got:",
             paste(names(bp), collapse = ",")))
if (!identical(bp$minor_ref[1], "6i_DQ835770"))
  fail(paste("minor-coherence: minor_ref must be the 6i reference, got", bp$minor_ref[1]))
# The whole point: the NAMED contig must be NODE_3, never NODE_2.
if (!grepl("^NODE_3_", bp$minor_contig[1]))
  fail(paste("minor-coherence: minor_contig must be NODE_3 (the off-genotype contig),",
             "not the 1a contig that outscores it on the same reference; got",
             bp$minor_contig[1]))
if (bp$minor_contig_length[1] != 1620)
  fail(paste("minor-coherence: minor_contig_length must be NODE_3's 1620, got",
             bp$minor_contig_length[1]))
# And the three fields must be mutually consistent: the length embedded in the
# contig name must equal the reported length.
name_len <- as.numeric(sub(".*_length_([0-9]+)_.*", "\\1", bp$minor_contig[1]))
if (!identical(name_len, as.numeric(bp$minor_contig_length[1])))
  fail(sprintf("minor-coherence: contig name says %s bp but minor_contig_length says %s",
               name_len, bp$minor_contig_length[1]))
ok("minor-coherence: minor_ref / minor_contig / minor_contig_length all describe ONE contig (2633901 regression)")

# --- Case E: several contigs share the minor reference as their top hit --------
# The changed length semantics, pinned. Two genotype-2 contigs both top-hit the same
# 2b reference: NODE_5 is LONGER (3000 bp) but aligns over only 100 bp; NODE_6 is
# shorter (1000 bp) but aligns over 900 and therefore wins the selection. The
# reported contig and length must be the SELECTED one, not the longest of the group —
# otherwise the pair goes back out of sync.
r_e <- run_blast_parse("minorgroup", "MINGRP", hits = list(
  node(4, 9000, 500.0, "1a_HQ850279", 99.0, 8900),
  node(5, 3000,   2.0, "2b_ACC",      90.0,  100),
  node(6, 1000,  30.0, "2b_ACC",      96.0,  900)
))
if (r_e$exit != 0) fail(paste("minor-group: exit", r_e$exit))
bpe <- r_e$blastparse
if (!grepl("^NODE_6_", bpe$minor_contig[1]))
  fail(paste("minor-group: minor_contig must be the selected NODE_6, got", bpe$minor_contig[1]))
if (bpe$minor_contig_length[1] != 1000)
  fail(paste("minor-group: minor_contig_length must be NODE_6's 1000 (the selected",
             "contig), not NODE_5's 3000 (the longest); got", bpe$minor_contig_length[1]))
ok("minor-group: with several contigs on the minor reference, the SELECTED contig is reported, not the longest")

cat("\nALL PASS\n")
