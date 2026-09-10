#!/usr/bin/env Rscript

# assembly_support_join.R -------------------------------------------------
# Genotype-level join of the per-subtype de novo assembly-support table
# (Plan 07-01's *.assembly_support.csv) onto the Phase-6 long-format candidate
# set (*.candidates.csv). This is the ASUP-02 join side (Phase 7).
#
# This is a SOURCED helper, not an entrypoint: it defines a single function
# (join_assembly_support) and has no commandArgs parser and no top-level file
# I/O. Consumers source() it from their task workdir (the file is staged there
# as a declared `path` process input by the calling Nextflow module), AFTER
# first sourcing genotype_utils.R so that genotype_from_subtype() is already in
# scope. This helper does NOT re-source genotype_utils.R (mirrors the
# denovo_layer.R / denovo_confirm.R sourced-helper convention).
#
# Design (07-CONTEXT D-02/D-03, criterion #3; 07-PATTERNS lines 76-131):
#   - Candidates are the LEFT side of the join => NO row loss (criterion #3).
#   - Match level is parameterized (default "genotype"): a 3b candidate is
#     corroborated by 3a support at genotype level (both genotype "3"), but only
#     by a 3b support row at subtype level (criterion #2 / the subtype non-match).
#   - The candidate-side match key reuses the precomputed Phase-6 columns
#     (candidate_genotype / candidate_subtype) rather than recomputing.
#   - The support-side match key is derived with the exact denovo_confirm.R line-64
#     branch: `if (match_level == "subtype") subtype else genotype_from_subtype(subtype)`.
#   - When several subtypes roll up to the same genotype, the genotype winner is
#     the single best contig by best_contig_length (D-03 carried through the
#     collapse): slice_max(best_contig_length, n = 1, with_ties = FALSE).
#   - A candidate with no support at its match level resolves to an explicit
#     assembly_support = "none" with NA metric columns (left_join NA-fill).
#   - Empty support_df (skip-assembly) => every candidate gets "none" + NA, no
#     abort. Empty candidates_df => a typed zero-row frame is returned, no abort.
#
# The legacy minor-coupled path (apply_denovo_layer / denovo_minor_*) is NOT
# touched here; it runs in parallel this phase (D-04).
# -------------------------------------------------------------------------

# Defensive: consumers already load tidyverse, which provides the pipe,
# group_by()/slice_max()/mutate()/left_join()/if_else(). This guard only fires
# if sourced into a session that has not.
if (!exists("group_by")) {
  library(tidyverse)
}

join_assembly_support <- function(candidates_df, support_df, match_level = "genotype") {
  # WR-04: reject any match_level other than the two supported values up front.
  # The downstream `if (match_level == "subtype") ... else ...` branch otherwise
  # treats every non-"subtype" value (a typo, NA, or empty string from a
  # mis-parsed arg) as the genotype path silently — a correctness risk for a
  # clinical genotyping tool.
  stopifnot(match_level %in% c("genotype", "subtype"))

  # The per-candidate support columns this function attaches. Declared once so
  # the typed zero-row path and the NA-fill path stay in lockstep.
  supported_cols <- c(
    "assembly_support",
    "assembly_support_subtype",
    "assembly_support_best_contig_length",
    "assembly_support_best_contig_pident",
    "assembly_support_best_contig_aln_length",
    "assembly_support_best_contig_kmer_cov"
  )

  # Zero-row candidate set: return a typed frame carrying the candidate columns
  # plus the new support columns, never abort (criterion #3 generalized).
  if (is.null(candidates_df) || nrow(candidates_df) == 0) {
    out <- candidates_df
    if (is.null(out)) {
      out <- tibble()
    }
    out <- out %>%
      mutate(
        assembly_support                        = character(),
        assembly_support_subtype                = character(),
        assembly_support_best_contig_length     = double(),
        assembly_support_best_contig_pident     = double(),
        assembly_support_best_contig_aln_length = double(),
        assembly_support_best_contig_kmer_cov   = double()
      )
    return(out)
  }

  # Candidate-side match key: reuse the precomputed Phase-6 columns
  # (07-PATTERNS line 118 — prefer reusing over recomputing). CR-01: coerce to
  # character so the key type is deterministic regardless of how the candidates
  # CSV was typed by readr. HCV genotypes 1–7 are purely-digit, so readr infers
  # <double> for candidate_genotype on a real run; the support side is always
  # character (genotype_from_subtype), and an un-coerced left_join would abort on
  # incompatible key types at the default match_level="genotype".
  cand <- candidates_df %>%
    mutate(.match_key = as.character(
      if (match_level == "subtype") candidate_subtype else candidate_genotype
    ))

  # Support-side match key: derive with the exact denovo_confirm.R line-64 branch
  # so both sides compute the key identically. Collapse to one row per
  # (sampleName, match_key) by the single best contig (D-03).
  if (is.null(support_df) || nrow(support_df) == 0) {
    # Empty support: typed zero-row collapsed frame so every candidate NA-fills.
    support_collapsed <- tibble(
      sampleName                              = character(),
      .match_key                              = character(),
      assembly_support_subtype                = character(),
      assembly_support_best_contig_length     = double(),
      assembly_support_best_contig_pident     = double(),
      assembly_support_best_contig_aln_length = double(),
      assembly_support_best_contig_kmer_cov   = double()
    )
  } else {
    support_collapsed <- support_df %>%
      mutate(.match_key = as.character(
        if (match_level == "subtype") subtype else genotype_from_subtype(subtype)
      )) %>%
      group_by(sampleName, .match_key) %>%
      slice_max(best_contig_length, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      # Carry the winning contig's subtype + four metrics under the assembly_support_* names.
      transmute(
        sampleName,
        .match_key,
        assembly_support_subtype                = subtype,
        assembly_support_best_contig_length     = best_contig_length,
        assembly_support_best_contig_pident     = best_contig_pident,
        assembly_support_best_contig_aln_length = best_contig_aln_length,
        assembly_support_best_contig_kmer_cov   = best_contig_kmer_cov
      )
  }

  # left_join candidates (LEFT) onto collapsed support by (sampleName, match_key)
  # so candidates anchor the left side — no candidate row is ever dropped.
  out <- cand %>%
    left_join(support_collapsed, by = c("sampleName", ".match_key")) %>%
    mutate(
      assembly_support = if_else(is.na(assembly_support_best_contig_length), "none", "supported")
    ) %>%
    select(-.match_key)

  out
}
