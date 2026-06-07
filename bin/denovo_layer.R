#!/usr/bin/env Rscript

# denovo_layer.R ----------------------------------------------------------
# De novo confirmation downgrade layer for Change 2 (D-02, D-13: DOWNGRADE-ONLY).
#
# This is a SOURCED helper, not an entrypoint: it defines a single function
# (apply_denovo_layer) and has no commandArgs parser and no top-level file I/O.
# Consumers source() it from their task workdir (the file is staged there as a
# declared `path` process input by the calling Nextflow module), AFTER first
# sourcing genotype_utils.R and denovo_confirm.R so that genotype_from_subtype()
# and classify_minor_denovo() are already in scope. This helper does NOT
# re-source those files (mirrors the denovo_confirm.R sourced-helper convention).
#
# apply_denovo_layer() is the behaviour-preserving extraction of the inline
# downgrade block formerly at summarize.R:778-806. It runs IMMEDIATELY after the
# minor_typable case_when so it can only ever flip a YES->NO (refute), never
# resurrect a suppressed minor (CONF-06, by construction). The 1a/1b allowance
# and the upstream 2k1b suppression are untouched.
#
# Genotypes are derived UNCONDITIONALLY from the mapping reference names
# (`<subtype>_<acc>` -> leading subtype token -> genotype_from_subtype()), so the
# layer is robust to GLUE absence (does NOT depend on Major_subtype/Minor_subtype).
# minor_denovo_status is ALWAYS present afterwards (D-12), on both branches.
#
# Defaults 1000 / 2.0 / 90 / "genotype" are calibration-VALIDATED (03-RESEARCH
# "Calibration Evidence") and mirror classify_minor_denovo() — do not change.
# -------------------------------------------------------------------------

# Defensive: consumers already load tidyverse, which provides the pipe,
# rowwise()/mutate()/ungroup()/if_else()/filter(). This guard only fires if
# sourced into a session that has not.
if (!exists("group_by")) {
  library(tidyverse)
}

apply_denovo_layer <- function(final, df_blast_out, denovo_confirm_minor,
                               min_len = 1000, min_kmer = 2.0, min_pid = 90,
                               match_level = "genotype") {
  if (isTRUE(denovo_confirm_minor)) {
    final %>%
      rowwise() %>%
      mutate(minor_denovo_status = {
        if (is.na(Minor_reference)) {
          NA_character_                                   # no minor candidate (D-12)
        } else {
          bo <- df_blast_out %>% filter(sampleName == .data$sampleName)
          classify_minor_denovo(
            bo,
            genotype_from_subtype(str_extract(Major_reference, "^[^_]+")),
            genotype_from_subtype(str_extract(Minor_reference, "^[^_]+")),
            min_len, min_kmer,
            min_pid, match_level
          )
        }
      }) %>%
      ungroup() %>%
      # Downgrade-only (D-13): refute flips minor_typable YES->NO; Minor_* columns
      # stay populated (D-10 — never null a Minor_* field on refute).
      mutate(minor_typable = if_else(
        !is.na(minor_denovo_status) & minor_denovo_status == "refuted", "NO", minor_typable
      ))
  } else {
    # Flag OFF (CONF-07 / D-16): bypass the layer entirely. minor_typable keeps its
    # pure legacy value; status = not_evaluated when a candidate exists, NA otherwise.
    final %>%
      mutate(minor_denovo_status = if_else(is.na(Minor_reference), NA_character_, "not_evaluated"))
  }
}
