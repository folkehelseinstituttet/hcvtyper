#!/usr/bin/env Rscript

# test_compat.R -----------------------------------------------------------
# Phase-9 compatibility / filename-migration regression suite (COMPAT-01..04,
# TEST-01, D-06/D-07/D-08). Auto-discovered by bin/tests/run_all.sh's
# `test_*.R` glob and run in the pinned r-regression container — no CI YAML edit.
#
# Two harnesses, modelled on the two sibling test files:
#   * Analog A (subprocess) — bin/tests/test_candidate_selection.R: drive the
#     REAL bin/summarize.R via system2("Rscript", ...) on synthetic candidate /
#     stats / depth fixtures staged in an isolated tempdir, then assert on the
#     emitted Summary.csv. Used for COMPAT-01/02/03.
#   * Analog B (function-source) — bin/tests/test_classify_roles.R: source the
#     REAL bin/classify_roles.R helpers and assert exact role / role_reason
#     values. Used for COMPAT-04 alongside the Summary.csv signal.
#
# Asserted behaviours:
#   COMPAT-01 (D-06/D-07): the five core strain-call columns (Major_reference,
#     Minor_reference, Major_genotype_mapping, Minor_genotype_mapping,
#     overall_sample_call) reproduce the human-verified golden values in
#     bin/tests/fixtures/compat_golden.csv for a monoinfection AND a co-infection
#     case.
#   COMPAT-02 (smoke): the co-infection case's `.cand1.`/`.cand2.` slot fixtures
#     parse via the Plan-03 candidate_rank join and populate non-NA Major_/Minor_
#     references that carry NO `_cand{rank}` slot suffix (RESEARCH Pitfall 1 guard).
#   COMPAT-03 (D-05, validate-only): Summary.csv carries the legacy Major_*/Minor_*
#     columns ALONGSIDE the new Major_role_*/Minor_role_*/overall_sample_call.
#   COMPAT-04 (D-12, end-to-end): a 1a/1b co-infection is allowed (1b candidate
#     role == "co-infection", overall_sample_call == "co-infection") and a 2k/1b
#     recombinant is suppressed (2k1b candidate role == "background",
#     role_reason == "recombinant_2k1b", not reported as a co-infection minor).
#
# Run from any cwd:  Rscript bin/tests/test_compat.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Resolve bin/ relative to this test file so it runs from any cwd.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "summarize.R")

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

# Helper: NA-tolerant equality for a single golden vs produced cell.
cell_equal <- function(a, b) {
  a_na <- is.na(a) | (is.character(a) & a %in% c("NA", ""))
  b_na <- is.na(b) | (is.character(b) & b %in% c("NA", ""))
  if (a_na && b_na) return(TRUE)
  if (a_na != b_na) return(FALSE)
  as.character(a) == as.character(b)
}

# -------------------------------------------------------------------------
# Subprocess harness: stage a one-sample summarize.R input tree in a tempdir,
# invoke the REAL script, read back Summary.csv.
#
# `cands` is a tibble with the 8-column long candidate schema (candidate_rank 1 =
# major, 2 = minor, ...). For each candidate the harness writes the matching
# `.cand{rank}.` stats / depth fixtures using the NEW slot naming so the run also
# exercises the Plan-03 candidate_rank join (COMPAT-02).
#
# Returns the Summary.csv tibble (or NULL if the script did not emit it).
# -------------------------------------------------------------------------
run_summarize <- function(case, sampleName, cands,
                          minRead = 500, minCov = 30, n_candidates = 2,
                          ref_length = 200, denovo = NULL, return_mqc = FALSE,
                          assembly_support = NULL) {
  wd <- tempfile(paste0("compat_", case, "_")); dir.create(wd)
  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)

  # summarize.R source()s these five helpers cwd-relative (bin/summarize.R L10-25).
  for (h in c("genotype_utils.R", "denovo_confirm.R", "denovo_layer.R",
              "assembly_support_join.R", "classify_roles.R")) {
    file.copy(file.path(bin_dir, h), file.path(wd, h), overwrite = TRUE)
  }

  # Input subdirs summarize.R reads from. The trimmed/, kraken_classified/, id/
  # loops iterate `1:length(files)` so they need at least one file each to avoid
  # a `1:0` empty-loop crash; the strain-call columns come from the four mapping
  # subdirs. The remaining dirs (denovo/, glue/, variation/, consensus_distance/)
  # are guarded with `if (length(...) > 0)` so we leave them absent.
  for (d in c("trimmed", "kraken_classified", "parsefirst_mapping",
              "stats_withdup", "stats_markdup", "depth", "id")) {
    dir.create(file.path(wd, d))
  }

  # --- samplesheet (args[1]); read at summarize.R L1072: select sample) ---
  samplesheet <- file.path(wd, "samplesheet.csv")
  write_csv(tibble(sample = sampleName), samplesheet)

  # --- trimmed/<sample>.log (cutadapt-style; trimmed loop) ---
  writeLines(c(
    "This is cutadapt",
    "Total read pairs processed:          1,000,000",
    "Pairs written (passing filters):       900,000 (90.0%)"
  ), file.path(wd, "trimmed", paste0(sampleName, ".log")))

  # --- kraken_classified/<sample>.kraken2.report.txt (kraken loop) ---
  # X2 (col 2) = clade reads, X6 (col 6) = rank name; loop pulls X2 where X6=="root".
  writeLines(paste(c("100.00", "900000", "900000", "R", "1", "root"),
                   collapse = "\t"),
             file.path(wd, "kraken_classified",
                       paste0(sampleName, ".kraken2.report.txt")))

  # --- id/<sample>.sequencerID.tsv (id loop) ---
  writeLines("@M01234:1:000000000:1:1101:1:1",
             file.path(wd, "id", paste0(sampleName, ".sequencerID.tsv")))

  # --- parsefirst_mapping/<sample>.candidates.csv (8-col long; the join source) ---
  cand_csv <- cands %>%
    transmute(
      sample              = sampleName,
      candidate_rank,
      candidate_ref,
      candidate_subtype,
      candidate_genotype,
      candidate_reads,
      candidate_cov,
      confirmation_status = "pass"
    )
  write_csv(cand_csv, file.path(wd, "parsefirst_mapping",
                                paste0(sampleName, ".candidates.csv")))

  # --- parsefirst_mapping/<sample>.parsefirstmapping.csv (legacy 10-col shim) ---
  major <- cands %>% filter(candidate_rank == 1) %>% slice(1)
  minor <- cands %>% filter(candidate_rank == 2) %>% slice(1)
  legacy <- tibble(
    sample             = sampleName,
    total_mapped_reads = sum(cands$candidate_reads),
    major_ref          = if (nrow(major) > 0) major$candidate_ref else NA_character_,
    major_reads        = if (nrow(major) > 0) major$candidate_reads else NA_real_,
    major_cov          = if (nrow(major) > 0) major$candidate_cov else NA_real_,
    minor_ref          = if (nrow(minor) > 0) minor$candidate_ref else NA_character_,
    minor_reads        = if (nrow(minor) > 0) minor$candidate_reads else NA_real_,
    minor_cov          = if (nrow(minor) > 0) minor$candidate_cov else NA_real_,
    minor_call         = if (nrow(minor) > 0) "yes" else "no",
    gate_flag          = "ok"
  )
  write_csv(legacy, file.path(wd, "parsefirst_mapping",
                              paste0(sampleName, ".parsefirstmapping.csv")))

  # --- read-count + depth fixtures (Plan 03 / JMAP-03: idxstats format) -------
  # Read counts now come from SAMTOOLS_IDXSTATS on the COMBINED BAM (pre/post
  # dedup), staged as ONE file PER SAMPLE: <sample>.withdup.idxstats /
  # <sample>.nodup.idxstats. Each is a 4-column no-header TSV
  # (refname, seqlen, mapped, unmapped) with ONE data row per candidate ref plus
  # a trailing `*` unmapped row. The idxstats refnames are the BARE candidate
  # reference (no `.cand{rank}.` filename slot — idxstats is per-sample, so rank
  # is recovered via the candidate_rank_lookup join in summarize.R, NOT a
  # filename parse). `withdup_reads`/`nodup_reads` optionally override the read
  # count written per candidate (defaults to candidate_reads).
  idx_seqlen <- 9456L  # representative HCV genome length; value is unused downstream

  withdup_rows <- character(0)
  nodup_rows   <- character(0)
  for (i in seq_len(nrow(cands))) {
    ref      <- cands$candidate_ref[i]
    wd_reads <- if (!is.null(cands$withdup_reads)) cands$withdup_reads[i] else cands$candidate_reads[i]
    nd_reads <- if (!is.null(cands$nodup_reads))   cands$nodup_reads[i]   else cands$candidate_reads[i]
    withdup_rows <- c(withdup_rows, paste(ref, idx_seqlen, wd_reads, 0, sep = "\t"))
    nodup_rows   <- c(nodup_rows,   paste(ref, idx_seqlen, nd_reads, 0, sep = "\t"))
  }
  # Trailing `*` unmapped row (must be dropped by the loop's candidate_ref != "*").
  withdup_rows <- c(withdup_rows, paste("*", 0, 0, 190, sep = "\t"))
  nodup_rows   <- c(nodup_rows,   paste("*", 0, 0, 190, sep = "\t"))
  writeLines(withdup_rows, file.path(wd, "stats_withdup",
                                     paste0(sampleName, ".withdup.idxstats")))
  writeLines(nodup_rows, file.path(wd, "stats_markdup",
                                   paste0(sampleName, ".nodup.idxstats")))

  # depth tsv: <ref> <pos> <cov>, 3-col no header (summarize.R coverage loop).
  # Depth is still per-candidate (post-split BAM) so it keeps the `.cand{rank}.`
  # slot the coverage loop recovers rank from via the lookup join.
  for (i in seq_len(nrow(cands))) {
    rk   <- cands$candidate_rank[i]
    ref  <- cands$candidate_ref[i]
    base <- paste0(sampleName, ".", ref, ".cand", rk)
    # Constant coverage across ref_length positions; high enough to clear gate.
    depth_lines <- paste(ref, seq_len(ref_length), 60, sep = "\t")
    writeLines(depth_lines, file.path(wd, "depth",
                                      paste0(base, ".nodup.tsv")))
  }

  # --- OPTIONAL de novo fixtures (RPT-CONTIG): stage denovo/ so summarize.R
  # resolves denovo_*_ref (from *.blastparse.csv) AND denovo_*_contig (top-bitscore
  # qseqid whose sseqid == that ref, from *_blast_out.csv). Absent by default, so
  # every existing caller stays a no-denovo run (four denovo_* columns NA-fill).
  # `denovo` is a named list: major_ref/minor_ref (blastparse) + major_contig/
  # minor_contig (the intended top-bitscore contig for each ref). Decoy lower-bitscore
  # hits on each ref prove slice_max picks the named contig, not just the only row.
  if (!is.null(denovo)) {
    dir.create(file.path(wd, "denovo"))
    # *.blastparse.csv: exact col-types summarize.R pins at L630-636.
    blastparse <- tibble(
      sample              = sampleName,
      major_ref           = denovo$major_ref,
      major_contig_length = 3000L,
      minor_ref           = denovo$minor_ref,
      # 260803-ogc: blast_parse.R now emits the minor contig NAME from the same
      # scaf_top row as minor_ref / minor_contig_length, and summarize.R takes it
      # from here rather than re-resolving it against the full hit table.
      minor_contig        = denovo$minor_contig,
      minor_contig_length = 2500L
    )
    write_csv(blastparse, file.path(wd, "denovo",
                                    paste0(sampleName, ".blastparse.csv")))

    # <sample>_blast_out.csv: outfmt6+ schema summarize.R df_blast_out reads
    # (L676-692). One row per (contig, ref) hit. Bitscore chosen so the top hit per
    # (sampleName, sseqid) carries the intended contig; other columns constant/dummy.
    blast_out <- tribble(
      ~qseqid,               ~sseqid,            ~bitscore,
      denovo$major_contig,   denovo$major_ref,   5000,   # top hit on major_ref
      "NODE_99_len_400",     denovo$major_ref,   900,    # decoy on major_ref (lower)
      denovo$minor_contig,   denovo$minor_ref,   4200,   # the SELECTED minor contig
      # 260803-ogc: this decoy deliberately OUTSCORES the selected minor contig on
      # the minor reference, reproducing sample 2633901 — where NODE_2 (an on-genotype
      # contig whose conserved 5'UTR/core region hit the off-genotype reference at
      # bitscore 1074) beat NODE_3 (the actual off-genotype contig, bitscore 97) in
      # summarize.R's top-bitscore-per-reference lookup. The old re-resolution would
      # name THIS contig; the fixed code must report denovo$minor_contig, taken from
      # blastparse.csv. The major decoy stays lower because the major slot legitimately
      # keeps the re-resolution (its ref is the globally best hit, so the grains agree).
      "NODE_98_len_350",     denovo$minor_ref,   9000    # decoy that must NOT win
    ) %>%
      mutate(
        subtype   = NA_character_, pident = 99, length = 3000, mismatch = 5,
        gapopen   = 1, qstart = 1, qend = 3000, sstart = 1, send = 3000,
        evalue    = 0, sc_length = 3000, kmer_cov = 40
      ) %>%
      select(qseqid, sseqid, subtype, pident, length, mismatch, gapopen,
             qstart, qend, sstart, send, evalue, bitscore, sc_length, kmer_cov)
    write_csv(blast_out, file.path(wd, "denovo",
                                   paste0(sampleName, "_blast_out.csv")))
  }

  # --- OPTIONAL Phase-7 assembly-support fixture (CR-01/CR-02/WR-05, 12-REVIEW):
  # stage <sample>.assembly_support.csv so join_assembly_support() (called from
  # summarize.R L799) populates real assembly_support_* metrics on candidate_support,
  # which score_assembly_support() -> classify_roles() then turns into a per-candidate
  # evidence_state (confirmed/probable/weak/refuted). Absent by default (assembly_exists
  # stays FALSE / evidence_state stays "weak" for every candidate, matching every
  # existing caller's fixtures, which never set this up — this is a distinct fixture
  # from the RPT-CONTIG `denovo` block above, which only staged the rescue-path
  # blastparse/blast_out files, never the Phase-7 per-subtype support table).
  # `assembly_support` is a tibble: subtype, best_contig_length, best_contig_pident,
  # best_contig_kmer_cov (best_contig_aln_length NA-fills — unused by scoring).
  # Exact column types pinned to match summarize.R's col_types (L773-780).
  if (!is.null(assembly_support)) {
    if (!dir.exists(file.path(wd, "denovo"))) dir.create(file.path(wd, "denovo"))
    support_csv <- assembly_support %>%
      transmute(
        sample                  = sampleName,
        subtype                 = subtype,
        best_contig_length      = as.double(best_contig_length),
        best_contig_pident      = as.double(best_contig_pident),
        best_contig_aln_length  = as.double(if ("best_contig_aln_length" %in% names(assembly_support)) best_contig_aln_length else NA_real_),
        best_contig_kmer_cov    = as.double(if ("best_contig_kmer_cov" %in% names(assembly_support)) best_contig_kmer_cov else NA_real_)
      )
    write_csv(support_csv, file.path(wd, "denovo",
                                     paste0(sampleName, ".assembly_support.csv")))
  }

  # --- invoke the REAL summarize.R (positional arg contract, L27-86) ---
  # args: [1]=samplesheet [2]=version [3]=name [4-8] denovo params (empty=default)
  #       [9]=minRead [10]=minCov [11]=n_candidates
  # NOTE: pass the empty positional slots [4-8] as shQuote("") — a bare "" is
  # DROPPED by system2(), which would silently shift minRead/minCov off positions
  # [9]/[10] (they would land at [4]/[5]), parse as denovo params, leave the gate
  # thresholds NA, and collapse every sample to overall_sample_call=="indeterminate".
  empty_slot <- shQuote("")
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(samplesheet), "TestPipe", "hcvtyper",
      empty_slot, empty_slot, empty_slot, empty_slot, empty_slot,
      minRead, minCov, n_candidates),
    stdout = FALSE, stderr = FALSE
  )

  summary_path <- file.path(wd, "Summary.csv")
  summary_tbl <- if (file.exists(summary_path)) {
    read_csv(summary_path, show_col_types = FALSE)
  } else {
    NULL
  }

  if (return_mqc) {
    # summary_mqc.tsv = manual colname header line + write_tsv data rows (no second
    # header). read_tsv treats the first line as the header (bin/summarize.R L1993-1995).
    mqc_path <- file.path(wd, "summary_mqc.tsv")
    mqc_tbl <- if (file.exists(mqc_path)) {
      read_tsv(mqc_path, show_col_types = FALSE)
    } else {
      NULL
    }
    return(list(summary = summary_tbl, mqc = mqc_tbl))
  }

  summary_tbl
}

# Candidate-frame builder (one row per candidate).
mk_cands <- function(...) {
  bind_rows(...)
}
mk_cand <- function(rank, ref, subtype, reads, cov,
                    withdup_reads = NA_real_, nodup_reads = NA_real_) {
  tibble(
    candidate_rank     = as.integer(rank),
    candidate_ref      = ref,
    candidate_subtype  = subtype,
    # genotype = leading-character rule (genotype_from_subtype), kept here as the
    # plain candidates.csv value summarize.R reads as character.
    candidate_genotype = if (subtype == "2k1b") "2k1b" else substr(subtype, 1, 1),
    candidate_reads    = reads,
    candidate_cov      = cov,
    # JMAP-03: optional per-rank read-count overrides written into the idxstats
    # fixtures (default to candidate_reads). Lets a test assert that rank 1 -> major
    # and rank 2 -> minor pick up the right idxstats column-3 value independently.
    withdup_reads      = ifelse(is.na(withdup_reads), reads, withdup_reads),
    nodup_reads        = ifelse(is.na(nodup_reads),   reads, nodup_reads)
  )
}

# =========================================================================
# COMPAT-01 / COMPAT-02 / COMPAT-03 — subprocess against the REAL summarize.R
# =========================================================================

golden_path <- file.path(this_dir, "fixtures", "compat_golden.csv")
if (!file.exists(golden_path)) fail(paste("compat_golden.csv missing at", golden_path))
golden <- read_csv(golden_path, show_col_types = FALSE)

core_cols <- c("Major_reference", "Minor_reference",
               "Major_genotype_mapping", "Minor_genotype_mapping",
               "overall_sample_call")

# --- monoinfection case (golden row MONO): 1a major, no minor ---
mono_golden <- golden %>% filter(sampleName == "MONO")
if (nrow(mono_golden) != 1) fail("compat_golden.csv missing MONO row")

mono_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99)
)
mono_summary <- run_summarize("mono", "MONO", mono_cands)
if (is.null(mono_summary)) fail("COMPAT-01 mono: summarize.R wrote no Summary.csv")
if (nrow(mono_summary) != 1) fail("COMPAT-01 mono: expected 1 Summary.csv row")
for (col in core_cols) {
  if (!col %in% names(mono_summary)) fail(paste("COMPAT-01 mono: Summary.csv missing column", col))
  if (!cell_equal(mono_golden[[col]][1], mono_summary[[col]][1]))
    fail(sprintf("COMPAT-01 mono: %s = '%s', golden = '%s'",
                 col, mono_summary[[col]][1], mono_golden[[col]][1]))
}
ok("COMPAT-01 mono: five core strain-call columns match compat_golden.csv (MONO)")

# --- co-infection case (golden row COINF): 1a major + 1b minor ---
coinf_golden <- golden %>% filter(sampleName == "COINF")
if (nrow(coinf_golden) != 1) fail("compat_golden.csv missing COINF row")

coinf_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99),
  mk_cand(2, "1b_D90208", "1b", 150000, 97)
)
# Plan 12-05: skip_assembly was removed, so every real run now always produces
# assembly_support.csv. Under the continuous-evidence model a candidate with NO
# own de novo support scores "weak" and can never be a co-infection minor (D-06),
# so an unstaged coinf run would collapse the golden co-infection to a false
# monoinfection. Stage realistic per-subtype assembly support so BOTH genotype-1
# candidates (1a major, 1b minor) land at evidence_state "confirmed" and the
# golden co-infection outcome (compat_golden.csv COINF) is reproduced. The join
# is genotype-level (assembly_support_join.R default), so both genotype-1
# candidates share the single best genotype-1 contig — this is concordant (the
# own-denovo conflict check is genotype-level: both map and assemble to gt1).
coinf_support <- tibble(
  subtype              = c("1a", "1b"),
  best_contig_length   = c(9076, 9339),
  best_contig_pident   = c(99, 99),
  best_contig_kmer_cov = c(40, 38)
)
coinf_summary <- run_summarize("coinf", "COINF", coinf_cands,
                               assembly_support = coinf_support)
if (is.null(coinf_summary)) fail("COMPAT-01 coinf: summarize.R wrote no Summary.csv")
if (nrow(coinf_summary) != 1) fail("COMPAT-01 coinf: expected 1 Summary.csv row")
for (col in core_cols) {
  if (!col %in% names(coinf_summary)) fail(paste("COMPAT-01 coinf: Summary.csv missing column", col))
  if (!cell_equal(coinf_golden[[col]][1], coinf_summary[[col]][1]))
    fail(sprintf("COMPAT-01 coinf: %s = '%s', golden = '%s'",
                 col, coinf_summary[[col]][1], coinf_golden[[col]][1]))
}
ok("COMPAT-01 coinf: five core strain-call columns match compat_golden.csv (COINF)")

# --- COMPAT-02 smoke: cand-slot fixtures parsed (non-NA, suffix-free refs) ---
maj <- coinf_summary$Major_reference[1]
min <- coinf_summary$Minor_reference[1]
if (is.na(maj)) fail("COMPAT-02: Major_reference silently empty after cand-slot join (Pitfall 1)")
if (is.na(min)) fail("COMPAT-02: Minor_reference silently empty after cand-slot join (Pitfall 1)")
if (grepl("cand[0-9]+", maj)) fail(paste("COMPAT-02: Major_reference carries a slot suffix:", maj))
if (grepl("cand[0-9]+", min)) fail(paste("COMPAT-02: Minor_reference carries a slot suffix:", min))
if (maj != "1a_M62321") fail(paste("COMPAT-02: Major_reference expected 1a_M62321, got", maj))
ok("COMPAT-02: .cand1./.cand2. fixtures populate suffix-free Major_/Minor_reference")

# --- COMPAT-03: legacy + role columns co-present (validate-only) ---
legacy_cols <- c("Major_reference", "Minor_reference",
                 "Major_genotype_mapping", "Minor_genotype_mapping")
role_cols   <- c("Major_role_reference", "Minor_role_reference", "overall_sample_call")
if (!all(legacy_cols %in% colnames(coinf_summary)))
  fail(paste("COMPAT-03: missing legacy columns:",
             paste(setdiff(legacy_cols, colnames(coinf_summary)), collapse = ",")))
if (!all(role_cols %in% colnames(coinf_summary)))
  fail(paste("COMPAT-03: missing role columns:",
             paste(setdiff(role_cols, colnames(coinf_summary)), collapse = ",")))
ok("COMPAT-03: Summary.csv carries legacy Major_*/Minor_* alongside role columns + overall_sample_call")

# =========================================================================
# COMPAT-04 — 1a/1b allowed, 2k/1b suppressed (D-12), exercised end-to-end
# through the FULL summarize.R path AND locked at the helper level.
# =========================================================================

# --- Part A: observable Summary.csv signal via the run_summarize() harness ---
# A 1a-dominant + 1b-corroborated sample must surface as a co-infection with the
# 1b ref in the Minor_role slot. We reuse the COMPAT-01 co-infection run (1a/1b),
# which already exercises the production score_candidates -> classify_roles call
# site inside summarize.R, and assert the role-level outcome here.
if (is.na(coinf_summary$overall_sample_call[1]) ||
    coinf_summary$overall_sample_call[1] != "co-infection")
  fail(sprintf("COMPAT-04 (1a/1b): overall_sample_call must be co-infection, got '%s'",
               coinf_summary$overall_sample_call[1]))
if (is.na(coinf_summary$Minor_role_reference[1]) ||
    coinf_summary$Minor_role_reference[1] != "1b_D90208")
  fail(sprintf("COMPAT-04 (1a/1b): Minor_role_reference must be the 1b ref, got '%s'",
               coinf_summary$Minor_role_reference[1]))
ok("COMPAT-04 (1a/1b): end-to-end summarize.R yields co-infection with 1b in the minor role slot")

# A genotype-2-dominant + 2k1b-second sample must NOT be reported as a co-infection
# minor — the 2k1b recombinant is suppressed to background, so its ref never reaches
# the Minor_role slot and the sample is not a co-infection.
suppress_cands <- mk_cands(
  mk_cand(1, "2a_ref",   "2a",   200000, 98),
  mk_cand(2, "2k1b_ref", "2k1b", 150000, 96)
)
suppress_summary <- run_summarize("suppress", "SUPPRESS", suppress_cands)
if (is.null(suppress_summary)) fail("COMPAT-04 (2k1b): summarize.R wrote no Summary.csv")
if (!is.na(suppress_summary$overall_sample_call[1]) &&
    suppress_summary$overall_sample_call[1] == "co-infection")
  fail("COMPAT-04 (2k1b): a suppressed 2k1b recombinant must NOT yield a co-infection sample")
if (!is.na(suppress_summary$Minor_role_reference[1]) &&
    suppress_summary$Minor_role_reference[1] == "2k1b_ref")
  fail("COMPAT-04 (2k1b): the suppressed 2k1b ref must NOT appear in Minor_role_reference")
ok("COMPAT-04 (2k1b): end-to-end summarize.R suppresses the 2k1b recombinant (not a co-infection minor)")

# --- Part B: lock the exact role / role_reason vocabulary at the helper level ---
# The wide Summary.csv only carries the dominant + corroborated-minor role slots, so
# it cannot surface the per-candidate role_reason. Source the REAL helpers and mirror
# the test_classify_roles.R sim1 (1a/1b) + 2k1b cases to assert the precise role and
# role_reason values the D-12 exception preserves. Scoped to ONLY these two cases
# (COMPAT-04) — the broader classifier behaviours live in test_classify_roles.R.
source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))

mk_helper_cand <- function(sample, ref, subtype, reads, cov, even,
                           sup_len = NA_real_, sup_kmer = NA_real_, sup_pid = NA_real_) {
  tibble(
    sampleName                            = sample,
    candidate_ref                         = ref,
    candidate_subtype                     = subtype,
    candidate_genotype                    = genotype_from_subtype(subtype),
    candidate_reads                       = reads,
    candidate_cov                         = cov,
    cv_evenness                           = even,
    assembly_support_best_contig_length   = sup_len,
    assembly_support_best_contig_kmer_cov = sup_kmer,
    assembly_support_best_contig_pident   = sup_pid
  )
}
classify_helper <- function(df) {
  classify_roles(score_candidates(df), minRead = 0, minCov = 0,
                 denovo_min_contig_length = 500, denovo_min_kmer_cov = 2.0,
                 denovo_min_blast_identity = 90, match_level = "genotype")
}
role_of   <- function(r, ref) r %>% filter(candidate_ref == ref) %>% pull(role)
reason_of <- function(r, ref) r %>% filter(candidate_ref == ref) %>% pull(role_reason)

# 1a/1b allowance (D-12): the 1b candidate is preserved as co-infection.
ab <- bind_rows(
  mk_helper_cand("ab", "1a_ref", "1a", 150000, 99, 0.90, sup_len = 9076, sup_kmer = 40, sup_pid = 99),
  mk_helper_cand("ab", "1b_ref", "1b", 120000, 97, 0.88, sup_len = 9339, sup_kmer = 38, sup_pid = 99)
)
r_ab <- classify_helper(ab)
if (!identical(role_of(r_ab, "1b_ref"), "co-infection"))
  fail("COMPAT-04 helper: 1a/1b cross-subtype-within-gt1 must keep the 1b candidate as co-infection")
if ((r_ab %>% pull(overall_sample_call) %>% unique()) != "co-infection")
  fail("COMPAT-04 helper: the 1a/1b sample must be a co-infection")
ok("COMPAT-04 helper (D-12): 1a/1b -> 1b candidate role == co-infection")

# 2k1b suppression (D-12): the 2k1b candidate is demoted to background/recombinant_2k1b.
k2 <- bind_rows(
  mk_helper_cand("k2", "2a_ref",   "2a",   200000, 98, 0.90, sup_len = 9000, sup_kmer = 40, sup_pid = 99),
  mk_helper_cand("k2", "2k1b_ref", "2k1b", 6000,   82, 0.72, sup_len = 8800, sup_kmer = 33, sup_pid = 96)
)
r_k2 <- classify_helper(k2)
if (!identical(role_of(r_k2, "2k1b_ref"), "background"))
  fail("COMPAT-04 helper: a 2k1b recombinant paired with a genotype-2 dominant must be background")
if (!identical(reason_of(r_k2, "2k1b_ref"), "recombinant_2k1b"))
  fail("COMPAT-04 helper: the 2k1b suppression reason must be recombinant_2k1b")
ok("COMPAT-04 helper (D-12): 2k1b -> background/recombinant_2k1b")

# --- GATE03-REMOVE: removing GATE-03 must NOT suppress minor when major has low nodup reads ---
# ERR1810469-class: rank-1 candidate has only 248 reads — below minRead=500.
# Old GATE-03 would flip minor_typable from YES to NO because Reads_nodup_mapped_major
# (248) <= minRead (500). After D4 removal it must stay YES (minor has good own coverage).
gate03_cands <- mk_cands(
  mk_cand(1, "3a_D17763", "3a", 248,  99),
  mk_cand(2, "1a_M62321", "1a", 1030, 99)
)
# Plan 12-05: skip_assembly removed => assembly always runs, so this fixture must
# stage assembly_support for BOTH candidates (as a real run would) or the new
# evidence model scores the 1a minor "weak" and flips minor_typable to NO for a
# reason unrelated to GATE-03. Major 3a (genotype 3) and minor 1a (genotype 1) are
# distinct genotypes, so each gets its own confirmed contig (no genotype-collapse
# cross-contamination). With the minor confirmed, minor_typable staying YES
# genuinely isolates the GATE-03-removal behaviour this test guards.
gate03_support <- tibble(
  subtype              = c("3a", "1a"),
  best_contig_length   = c(9000, 9076),
  best_contig_pident   = c(99, 99),
  best_contig_kmer_cov = c(40, 40)
)
gate03_summary <- run_summarize("gate03remove", "G03", gate03_cands,
                                minRead = 500, minCov = 30,
                                assembly_support = gate03_support)
if (is.null(gate03_summary)) fail("GATE03-REMOVE: summarize.R wrote no Summary.csv")
if (!identical(gate03_summary$minor_typable[1], "YES"))
  fail(sprintf("GATE03-REMOVE: minor_typable must be YES after GATE-03 removal, got '%s'",
               gate03_summary$minor_typable[1]))
ok("GATE03-REMOVE: minor_typable stays YES despite low-nodup-read major (GATE-03 removed per D4)")

# =========================================================================
# RESCUE-AUDIT (Phase 10, D-08): Summary.csv carries the rescue audit columns,
# NA-filled on a no-rescue input. Drives the REAL summarize.R via the same
# run_summarize() harness. The candidates fixture written by run_summarize()
# has NO rescued_from/rescue_trigger columns (a no-rescue / pre-rescue input),
# so summarize.R must (a) surface the wide audit columns regardless and
# (b) leave them NA and rescue_flag FALSE when nothing was rescued.
#
# Column spelling uses the in-file cand_{rank}_{value} pivot convention
# (bin/summarize.R candidate_support_wide, names_glue "cand_{candidate_rank}_{.value}"),
# locked in Plan 03 — NOT the cand2_ literal from CONTEXT D-12.
#
# RED until Plan 03 adds rescued_from/rescue_trigger to the candidate_support
# pivot and the rescue_flag rollup in summarize.R.
# -------------------------------------------------------------------------
rescue_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99),
  mk_cand(2, "1b_D90208", "1b", 150000, 97)
)
rescue_summary <- run_summarize("rescueaudit", "RESCUEAUDIT", rescue_cands)
if (is.null(rescue_summary)) fail("RESCUE-AUDIT: summarize.R wrote no Summary.csv")

rescue_audit_cols <- c("rescue_flag",
                       "cand_1_rescued_from", "cand_2_rescued_from",
                       "cand_1_rescue_trigger", "cand_2_rescue_trigger")
missing_rescue <- setdiff(rescue_audit_cols, colnames(rescue_summary))
if (length(missing_rescue) > 0)
  fail(paste("RESCUE-AUDIT: Summary.csv missing rescue audit columns:",
             paste(missing_rescue, collapse = ",")))

# On a no-rescue input, rescue_flag must be FALSE and every cand_*_rescued_from NA.
if (!cell_equal(rescue_summary$rescue_flag[1], FALSE) &&
    !identical(as.logical(rescue_summary$rescue_flag[1]), FALSE))
  fail(sprintf("RESCUE-AUDIT: rescue_flag must be FALSE on a no-rescue input, got '%s'",
               rescue_summary$rescue_flag[1]))
for (col in c("cand_1_rescued_from", "cand_2_rescued_from")) {
  v <- rescue_summary[[col]][1]
  if (!(is.na(v) || (is.character(v) && v %in% c("NA", ""))))
    fail(sprintf("RESCUE-AUDIT: %s must be NA on a no-rescue input, got '%s'", col, v))
}
ok("RESCUE-AUDIT (D-08): Summary.csv carries rescue_flag + cand_{rank}_rescued_from/rescue_trigger, NA on no-rescue")

# =========================================================================
# JMAP-03 (Phase 11, D-15): read-count columns are sourced from SAMTOOLS_IDXSTATS
# (combined BAM, one file per sample) instead of the removed SAMTOOLS_STATS.
# The run_summarize() harness now stages <sample>.withdup.idxstats /
# <sample>.nodup.idxstats 4-column TSV fixtures (one row per candidate + a `*`
# trailer). Assert the four read-count columns populate from idxstats column 3,
# rank 1 -> *_major, rank 2 -> *_minor, with the withdup/nodup values distinct.
#
# RED until Plan 03 migrates summarize.R's path_4/path_5 loops to parse idxstats.
# -------------------------------------------------------------------------
jmap_cands <- mk_cands(
  mk_cand(1, "3a_D17763",   "3a", 200000, 99, withdup_reads = 17257, nodup_reads = 1007),
  mk_cand(2, "3i_JX227955", "3i", 150000, 97, withdup_reads = 227,   nodup_reads = 47)
)
jmap_summary <- run_summarize("jmapidx", "JMAPIDX", jmap_cands)
if (is.null(jmap_summary)) fail("JMAP-03: summarize.R wrote no Summary.csv")
if (nrow(jmap_summary) != 1) fail("JMAP-03: expected exactly 1 Summary.csv row")

read_cols <- c("Reads_withdup_mapped_major", "Reads_withdup_mapped_minor",
               "Reads_nodup_mapped_major",   "Reads_nodup_mapped_minor")
for (col in read_cols) {
  if (!col %in% names(jmap_summary))
    fail(paste("JMAP-03: Summary.csv missing read-count column", col))
  v <- jmap_summary[[col]][1]
  if (is.na(v))
    fail(paste("JMAP-03:", col, "is NA — idxstats read counts did not populate"))
}
# Exact idxstats column-3 values, rank 1 -> major, rank 2 -> minor.
if (as.numeric(jmap_summary$Reads_withdup_mapped_major[1]) != 17257)
  fail(sprintf("JMAP-03: Reads_withdup_mapped_major = %s, expected 17257 (idxstats col 3, rank 1)",
               jmap_summary$Reads_withdup_mapped_major[1]))
if (as.numeric(jmap_summary$Reads_withdup_mapped_minor[1]) != 227)
  fail(sprintf("JMAP-03: Reads_withdup_mapped_minor = %s, expected 227 (idxstats col 3, rank 2)",
               jmap_summary$Reads_withdup_mapped_minor[1]))
if (as.numeric(jmap_summary$Reads_nodup_mapped_major[1]) != 1007)
  fail(sprintf("JMAP-03: Reads_nodup_mapped_major = %s, expected 1007 (idxstats col 3, rank 1)",
               jmap_summary$Reads_nodup_mapped_major[1]))
if (as.numeric(jmap_summary$Reads_nodup_mapped_minor[1]) != 47)
  fail(sprintf("JMAP-03: Reads_nodup_mapped_minor = %s, expected 47 (idxstats col 3, rank 2)",
               jmap_summary$Reads_nodup_mapped_minor[1]))
ok("JMAP-03: Reads_withdup/nodup_mapped_major/minor sourced from idxstats column 3 (rank 1->major, rank 2->minor)")

# =========================================================================
# JMAP-03 shared-reference guard (Pitfall 4 / commit 1d1a051): when two
# candidates share ONE reference name, the combined-BAM idxstats has a single
# row for that ref. The candidate_rank_lookup join must still disambiguate the
# two ranks WITHOUT exploding into a many-to-many join — exactly ONE Summary.csv
# row per sample (no duplicated sample rows in df_mapped_reads).
#
# RED until Plan 03 builds the lookup so both ranks of a shared ref survive and
# the df_with_dups/df_nodups full_join collapses to one row per sample.
# -------------------------------------------------------------------------
shared_cands <- mk_cands(
  mk_cand(1, "3a_D17763", "3a", 200000, 99, withdup_reads = 5000, nodup_reads = 400),
  mk_cand(2, "3a_D17763", "3a", 150000, 97, withdup_reads = 5000, nodup_reads = 400)
)
shared_summary <- run_summarize("jmapshared", "JMAPSHARED", shared_cands)
if (is.null(shared_summary)) fail("JMAP-03 shared-ref: summarize.R wrote no Summary.csv")
if (nrow(shared_summary) != 1)
  fail(sprintf("JMAP-03 shared-ref: expected exactly 1 Summary.csv row, got %d (many-to-many join, 1d1a051 regression)",
               nrow(shared_summary)))
ok("JMAP-03 shared-ref: two candidates sharing one reference yield exactly one Summary.csv row (many-to-many guard)")

# =========================================================================
# RPT-CONTIG (260703-dpj): de novo contig NAME + its BLAST ref surfaced per
# strain. Summary.csv gains denovo_major_contig / denovo_minor_contig (the
# top-bitscore contig qseqid backing each strain's ref) alongside the already-
# present denovo_major_ref / denovo_minor_ref. summary_mqc.tsv carries the generic
# denovo_best_contig / denovo_best_contig_ref, with the [minor] row promoted from
# the denovo_minor_* counterparts (col_major/col_minor pairing) — so a [minor] row
# shows the MINOR strain's contig+ref, not the major's (T-dpj-03).
#
# Drives the REAL summarize.R via run_summarize(denovo=, return_mqc=TRUE), which
# stages denovo/<sample>.blastparse.csv + <sample>_blast_out.csv fixtures.
# -------------------------------------------------------------------------
contig_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99),
  mk_cand(2, "1b_D90208", "1b", 150000, 97)
)
contig_denovo <- list(
  major_ref    = "1a_M62321", major_contig = "NODE_1_len_3000",
  minor_ref    = "1b_D90208", minor_contig = "NODE_7_len_2500"
)
# Plan 12-05: the [minor] summary_mqc.tsv row only exists for a co-infection, so
# (like COMPAT-01 coinf) this run must stage assembly_support to confirm the 1b
# minor now that skip_assembly is gone and assembly is load-bearing. The `denovo`
# fixture above feeds the orthogonal denovo_*_ref/contig reporting columns, NOT
# the assembly_support_* evidence columns that drive evidence_state.
contig_support <- tibble(
  subtype              = c("1a", "1b"),
  best_contig_length   = c(9076, 9339),
  best_contig_pident   = c(99, 99),
  best_contig_kmer_cov = c(40, 38)
)
contig_out <- run_summarize("rptcontig", "RPTCONTIG", contig_cands,
                            denovo = contig_denovo,
                            assembly_support = contig_support,
                            return_mqc = TRUE)
cs <- contig_out$summary
cm <- contig_out$mqc
if (is.null(cs)) fail("RPT-CONTIG: summarize.R wrote no Summary.csv")
if (is.null(cm)) fail("RPT-CONTIG: summarize.R wrote no summary_mqc.tsv")

# --- Summary.csv: all four denovo_*_ref / denovo_*_contig present AND populated ---
contig_cols <- c("denovo_major_ref", "denovo_major_contig",
                 "denovo_minor_ref", "denovo_minor_contig")
missing_contig <- setdiff(contig_cols, colnames(cs))
if (length(missing_contig) > 0)
  fail(paste("RPT-CONTIG: Summary.csv missing columns:",
             paste(missing_contig, collapse = ",")))
if (!identical(as.character(cs$denovo_major_ref[1]), "1a_M62321"))
  fail(sprintf("RPT-CONTIG: denovo_major_ref = '%s', expected 1a_M62321", cs$denovo_major_ref[1]))
if (!identical(as.character(cs$denovo_major_contig[1]), "NODE_1_len_3000"))
  fail(sprintf("RPT-CONTIG: denovo_major_contig = '%s', expected NODE_1_len_3000 (top-bitscore contig on 1a_M62321)", cs$denovo_major_contig[1]))
if (!identical(as.character(cs$denovo_minor_ref[1]), "1b_D90208"))
  fail(sprintf("RPT-CONTIG: denovo_minor_ref = '%s', expected 1b_D90208", cs$denovo_minor_ref[1]))
if (!identical(as.character(cs$denovo_minor_contig[1]), "NODE_7_len_2500"))
  fail(sprintf(paste("RPT-CONTIG: denovo_minor_contig = '%s', expected NODE_7_len_2500 —",
                     "the contig blast_parse.R SELECTED, not NODE_98_len_350 which",
                     "outscores it on 1b_D90208 (the 2633901 defect: a re-resolution",
                     "against the full hit table names a contig the minor selection",
                     "had excluded)"), cs$denovo_minor_contig[1]))
ok("RPT-CONTIG: Summary.csv denovo_major/minor_contig + _ref present and populated per strain")

# --- summary_mqc.tsv: generic columns exist; [major] row = major strain's
#     contig+ref, [minor] row = minor strain's (per-strain correctness, T-dpj-03) ---
mqc_contig_cols <- c("denovo_best_contig", "denovo_best_contig_ref")
missing_mqc <- setdiff(mqc_contig_cols, colnames(cm))
if (length(missing_mqc) > 0)
  fail(paste("RPT-CONTIG: summary_mqc.tsv missing columns:",
             paste(missing_mqc, collapse = ",")))
major_row <- cm %>% filter(sampleName == "RPTCONTIG [major]")
minor_row <- cm %>% filter(sampleName == "RPTCONTIG [minor]")
if (nrow(major_row) != 1) fail("RPT-CONTIG: summary_mqc.tsv missing the 'RPTCONTIG [major]' row")
if (nrow(minor_row) != 1) fail("RPT-CONTIG: summary_mqc.tsv missing the 'RPTCONTIG [minor]' row")
if (!identical(as.character(major_row$denovo_best_contig[1]), "NODE_1_len_3000"))
  fail(sprintf("RPT-CONTIG: [major] denovo_best_contig = '%s', expected NODE_1_len_3000", major_row$denovo_best_contig[1]))
if (!identical(as.character(major_row$denovo_best_contig_ref[1]), "1a_M62321"))
  fail(sprintf("RPT-CONTIG: [major] denovo_best_contig_ref = '%s', expected 1a_M62321", major_row$denovo_best_contig_ref[1]))
if (!identical(as.character(minor_row$denovo_best_contig[1]), "NODE_7_len_2500"))
  fail(sprintf("RPT-CONTIG: [minor] denovo_best_contig = '%s', expected NODE_7_len_2500 (minor strain's contig, not the major's)", minor_row$denovo_best_contig[1]))
if (!identical(as.character(minor_row$denovo_best_contig_ref[1]), "1b_D90208"))
  fail(sprintf("RPT-CONTIG: [minor] denovo_best_contig_ref = '%s', expected 1b_D90208 (minor strain's ref, not the major's)", minor_row$denovo_best_contig_ref[1]))
ok("RPT-CONTIG: summary_mqc.tsv [major]/[minor] rows carry the correct per-strain contig + ref via col_major/col_minor pairing")

# =========================================================================
# EVID-02 (Phase 12, folded todo 2026-06-22): a co-infection whose MINOR
# candidate's combined-BAM idxstats row reports mapped=0 (competitive joint
# mapping assigned all reads to the dominant, Phase 11) must NOT be silently
# dropped. Before the fix, summarize.R's `filter(mapped > 0)` in both idxstats
# parse loops discarded the minor row, Minor_reference collapsed to NA, and the
# sample was reported as a false monoinfection. After relaxing the guard to
# `filter(candidate_ref != "*")`, the zero-read minor survives to Summary.csv so
# the downstream evidence engine — not a silent upstream filter — decides its role.
#
# The minor is a genuine trace co-infection: selected during first-pass with a
# low-but-above-floor candidate_reads count (> minRead), yet its withdup/nodup
# idxstats mapped count is 0 because competitive joint mapping assigned all its
# reads to the dominant. Its first-pass count stays trace-level so the dominant
# (candidate 1) remains the role-dominant and the zero-read minor genuinely lands
# in the Minor slot. Assert Minor_reference is retained (populated, not NA).
# -------------------------------------------------------------------------
zeroread_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99, withdup_reads = 200000, nodup_reads = 15000),
  mk_cand(2, "1b_D90208", "1b", 1000,   97, withdup_reads = 0,      nodup_reads = 0)
)
zeroread_summary <- run_summarize("zeroread", "ZEROREAD", zeroread_cands)
if (is.null(zeroread_summary)) fail("EVID-02 zero-read: summarize.R wrote no Summary.csv")
if (nrow(zeroread_summary) != 1)
  fail(sprintf("EVID-02 zero-read: expected exactly 1 Summary.csv row, got %d", nrow(zeroread_summary)))
# The minor reference must be retained (populated, not silently dropped to NA).
zr_minor <- zeroread_summary$Minor_reference[1]
if (is.na(zr_minor) || (is.character(zr_minor) && zr_minor %in% c("NA", "")))
  fail("EVID-02 zero-read: Minor_reference was silently dropped to NA for a 0-mapped-read minor (mapped>0 filter regression)")
if (zr_minor != "1b_D90208")
  fail(sprintf("EVID-02 zero-read: Minor_reference = '%s', expected 1b_D90208 (the retained zero-read minor)", zr_minor))
# The retained minor carries a numeric 0 read count (flows through, not NA-dropped).
if (!cell_equal(zeroread_summary$Reads_nodup_mapped_minor[1], 0))
  fail(sprintf("EVID-02 zero-read: Reads_nodup_mapped_minor = '%s', expected 0 (the zero-read minor's idxstats col 3)",
               zeroread_summary$Reads_nodup_mapped_minor[1]))
ok("EVID-02 zero-read: a 0-mapped-read co-infection minor is retained in Summary.csv (not silently dropped)")

# =========================================================================
# CR-01 / WR-05 (12-REVIEW): any_probable_only review_flag + call_confidence
# trigger, exercised end-to-end through the REAL summarize.R (not just at the
# classify_roles() helper level, per the Fix note: "add a test_compat.R case
# that exercises the review_flag/call_confidence text end-to-end so a future
# role_reason rename cannot silently break it again"). A 1a-dominant + 2c-minor
# sample where the 2c candidate's own de novo contig scores in the marginal
# "probable" band (assembly_support_score ~0.577, in [0.50, 0.72)) must surface
# the co-infection minor AND fire the any_probable_only review sentence /
# demote call_confidence to "provisional" — the exact signal the retired
# any_uncorroborated trigger (built on the no-longer-reachable role_reason ==
# "uncorroborated_kept") could never produce.
# -------------------------------------------------------------------------
probable_cands <- mk_cands(
  mk_cand(1, "1a_M62321", "1a", 200000, 99),
  mk_cand(2, "2c_JX227949", "2c", 5000, 50)
)
probable_support <- tibble(
  subtype                = "2c",
  best_contig_length     = 3000,
  best_contig_pident     = 87,
  best_contig_kmer_cov   = NA_real_
)
probable_summary <- run_summarize("probable", "PROBABLE", probable_cands,
                                  assembly_support = probable_support)
if (is.null(probable_summary)) fail("CR-01 probable-band: summarize.R wrote no Summary.csv")
if (is.na(probable_summary$overall_sample_call[1]) ||
    probable_summary$overall_sample_call[1] != "co-infection")
  fail(sprintf("CR-01 probable-band: overall_sample_call must be co-infection, got '%s'",
               probable_summary$overall_sample_call[1]))
if (is.na(probable_summary$Minor_role_reference[1]) ||
    probable_summary$Minor_role_reference[1] != "2c_JX227949")
  fail(sprintf("CR-01 probable-band: Minor_role_reference must be the 2c ref, got '%s'",
               probable_summary$Minor_role_reference[1]))
pb_flag <- probable_summary$review_flag[1]
if (is.na(pb_flag) || !grepl("marginally corroborated", pb_flag))
  fail(sprintf("CR-01 probable-band: review_flag must carry the any_probable_only sentence, got '%s'",
               pb_flag))
if (is.na(probable_summary$call_confidence[1]) ||
    probable_summary$call_confidence[1] != "provisional")
  fail(sprintf("CR-01 probable-band: call_confidence must be 'provisional', got '%s'",
               probable_summary$call_confidence[1]))
# WR-05: the confirmed/probable/weak/refuted band must also be surfaced as its
# own wide Summary.csv column, not just folded into the review_flag prose.
if (is.na(probable_summary$Minor_evidence_state[1]) ||
    probable_summary$Minor_evidence_state[1] != "probable")
  fail(sprintf("WR-05: Minor_evidence_state must be 'probable', got '%s'",
               probable_summary$Minor_evidence_state[1]))
if (is.na(probable_summary$Major_evidence_state[1]) ||
    probable_summary$Major_evidence_state[1] != "weak")
  fail(sprintf("WR-05: Major_evidence_state must be 'weak' (no assembly_support fixture staged for the major), got '%s'",
               probable_summary$Major_evidence_state[1]))
ok("CR-01/WR-05: a probable-band (evidence_state=probable) co-infection minor fires the any_probable_only review_flag sentence, demotes call_confidence to provisional, and surfaces Major_evidence_state/Minor_evidence_state in Summary.csv, end-to-end via the real summarize.R")

cat("\nALL PASS\n")
