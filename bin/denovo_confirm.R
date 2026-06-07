#!/usr/bin/env Rscript

# denovo_confirm.R --------------------------------------------------------
# Pure de novo confirmation core for Change 2 (CONF-01..05).
#
# This is a SOURCED helper, not an entrypoint: it defines a single function
# and has no commandArgs parser. Consumers source() it from their task workdir
# (the file is staged there as a declared `path` process input by the calling
# Nextflow module), after first sourcing genotype_utils.R so that
# genotype_from_subtype() is already in scope.
#
# classify_minor_denovo() takes a per-contig BLAST data frame (the in-memory
# *_blast_out.csv schema produced by blast_parse.R) plus the mapping major/minor
# genotypes, and returns exactly one of:
#   "confirmed_by_denovo" — a substantial minor-genotype contig exists (CONF-01)
#   "refuted"             — a substantial major contig but no substantial minor (CONF-02)
#   "unconfirmed"         — de novo failed overall; never refute (CONF-03)
#   NA_character_         — no minor candidate to evaluate
#
# Purity is a firm requirement (D-08): no file I/O, no commandArgs, no global
# mutation. This lets Phase 4's TEST-01 exercise every branch on synthetic input
# without a pipeline run.
#
# Defaults 1000 / 2.0 / 90 / "genotype" are calibration-VALIDATED (03-RESEARCH
# "Calibration Evidence") — do not change. Substantiality uses the full contig
# length (sc_length), NOT the BLAST alignment length column.
# -------------------------------------------------------------------------

# Defensive: consumers already load tidyverse, which provides the pipe,
# arrange()/group_by()/slice()/mutate() and if_else(). This guard only fires if
# sourced into a session that has not.
if (!exists("group_by")) {
  library(tidyverse)
}

classify_minor_denovo <- function(blast_out_df, major_geno, minor_geno,
                                  min_len = 1000, min_kmer = 2.0, min_pid = 90,
                                  match_level = "genotype") {
  # 1. No minor candidate -> NA (caller also guards; defensive here).
  if (is.na(minor_geno)) {
    return(NA_character_)
  }

  # 2. No de novo evidence at all -> de novo failed -> unconfirmed (CONF-03,
  #    NEVER refute). This is the empty/NULL DoS guard (T-03-01): return a value,
  #    never stop().
  if (is.null(blast_out_df) || nrow(blast_out_df) == 0) {
    return("unconfirmed")
  }

  # 3. Best hit per contig — mirror blast_parse.R scaf_top (lines 148-152).
  top <- blast_out_df %>%
    arrange(evalue, desc(bitscore)) %>%
    group_by(qseqid) %>%
    slice(1) %>%
    ungroup()

  # 4. Match key: raw subtype if match_level == "subtype", else genotype level
  #    via genotype_from_subtype() (default; CONF-05 — genotype not subtype).
  #    5. substantial = ANDed floors on sc_length (full contig length), kmer_cov
  #       and pident — NOT the BLAST alignment length column.
  top <- top %>%
    mutate(
      match_key = if (match_level == "subtype") subtype else genotype_from_subtype(subtype),
      substantial = sc_length >= min_len & kmer_cov >= min_kmer & pident >= min_pid
    )

  # 6. Presence of a substantial contig at each genotype.
  has_minor_sub <- any(top$substantial & top$match_key == minor_geno, na.rm = TRUE)
  has_major_sub <- any(top$substantial & top$match_key == major_geno, na.rm = TRUE)

  # 7. Asymmetric return order (D-05): minor present wins, then major-present
  #    refutes, else de novo failed overall -> unconfirmed.
  if (has_minor_sub) {
    return("confirmed_by_denovo")  # CONF-01 (D-04)
  }
  if (has_major_sub) {
    return("refuted")              # CONF-02 (D-05): major present, minor absent
  }
  "unconfirmed"                    # CONF-03 (D-06): de novo failed overall
}
