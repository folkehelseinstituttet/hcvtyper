#!/usr/bin/env Rscript

library(tidyverse)
library(seqinr)

# Source the canonical 2k1b-aware genotype-from-subtype helper.
# Relative path: the file is staged into the task workdir as a declared
# `path(genotype_utils)` process input (PARSEFIRSTMAPPING). Do NOT use an
# absolute or projectDir path — that would break container portability.
source("genotype_utils.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Usage: summarize_mapping_to_all_references.R <idxstats file> <depth file> <sample name> <references> <minRead> <minCov> [n_candidates]", call.=FALSE)
}

idxstats   <- args[1]
depth      <- args[2]
sampleName <- args[3]
references <- args[4]
# V5 input validation: as.numeric coerces; a non-numeric arg yields NA, which
# fails the `>` gate comparisons safely (no crash, no minor_call='yes') rather
# than producing a NumberFormatException. Values originate from tracked config
# (conf/modules_hcv.config minRead/minCov), not external input.
minRead    <- as.numeric(args[5])
minCov     <- as.numeric(args[6])
# Phase 6 (REFSEL-02): number of neutrally-ranked candidates to select. New
# positional arg after minCov; as.integer mirrors the minRead/minCov coercion.
# Default 2 when absent (reproduces the legacy two-slot topology under the shim).
# A non-integer arg coerces to NA -> fall back to the safe default 2 (T-06-02).
n_candidates <- if (length(args) >= 7) as.integer(args[7]) else 2L
if (is.na(n_candidates) || n_candidates < 1L) n_candidates <- 2L

# First calculate coverage for all references
# Read the depth file from the first mapping.
# The file can be empty (no reads mapped) -> a 0-column tibble with no X1, which
# crashes group_by(X1). Guard it (T-06-01): on an empty/X1-less depth frame emit
# a well-typed 0-row cov frame so the no-mapping branch reaches its gate_flag=
# "no_mapping" default instead of erroring before any output is written.
depth_raw <- read_tsv(depth, col_names = FALSE)
if (nrow(depth_raw) == 0 || !("X1" %in% colnames(depth_raw))) {
  cov <- tibble(
    X1               = character(0),
    total_rows       = integer(0),
    count_gt_4       = integer(0),
    percent_gt_4     = numeric(0),
    percent_gt_4_int = numeric(0)
  )
} else {
  cov <- depth_raw %>%
    group_by(X1) %>% # Group by name of the reference
    summarise(
      total_rows = n(), # Get the total number of positions for the reference (genome length)
      count_gt_4 = sum(X3 > 4), # Get the number of positions with coverage >= 5
      percent_gt_4 = (count_gt_4 / total_rows) * 100,
      percent_gt_4_int = round(percent_gt_4, digits = 0) # Round to nearest integer. Need integer for groovy/nextflow filtering later
    )
}

# Then read mapped reads from the first mapping
# Read the idsxtats output of the first mapping
df <- read_table(idxstats, col_names = FALSE) %>%
  # Separate subtype and reference from the sequence names
  separate(X1, into = c("Subtype", "Reference"), sep = "_", remove = FALSE) %>%
  # Discard the unmapped reads marked by an * (more precisely these are unmapped reads without coordinates)
  filter(X1 != "*") %>%
  # Separate the genotype from the subtype.
  # For 2k1b we use the whole name for genotype also (see bin/genotype_utils.R).
  mutate(Genotype = genotype_from_subtype(Subtype))

# Join the percent coverage to the mapping statistics
df <- left_join(df, cov, by = c("X1" = "X1"))

# ---------------------------------------------------------------------------
# Phase 6 long-format candidate table (D-01). One row per neutrally-ranked
# candidate. This is the new downstream contract (consumed by Plans 02/03).
# Built first as an empty 0-row frame so the no-mapping branch still writes a
# valid (header-only) candidates CSV without crashing.
# ---------------------------------------------------------------------------
candidates_long <- tibble(
  sample              = character(0),
  candidate_rank      = integer(0),
  candidate_ref       = character(0),
  candidate_subtype   = character(0),
  candidate_genotype  = character(0),
  candidate_reads     = numeric(0),
  candidate_cov       = numeric(0),
  confirmation_status = character(0)
)

# Create empty final (legacy, shim) dataframe to populate (D-06). All 10 columns
# are ALWAYS present so bin/summarize.R can pull them unconditionally (NA tolerated).
df_final <- as.data.frame(matrix(nrow = 1, ncol = 10))
colnames(df_final) <- c("sample", "total_mapped_reads", "major_ref", "major_reads", "major_cov", "minor_ref", "minor_reads", "minor_cov", "minor_call", "gate_flag")

# Add sample name
df_final$sample[1] <- sampleName

# Gate-decision defaults. These always carry a value so a row is never emitted
# with an unexplained empty gate state (D-07). The empty-df / no-mapping branch
# leaves these defaults in place; the populated branch overwrites them below.
df_final$minor_call[1] <- "no"
df_final$gate_flag[1]  <- "no_mapping"

# Track the selected candidate references (rank-keyed) for the FASTA write below.
# Empty by default so the no-mapping branch writes no FASTA.
selected_refs <- character(0)

# Sometimes the mappings stats are completely empty
if (nrow(df) > 0) {

  # First get the total number of mapped reads to all references
  df_final$total_mapped_reads[1] <- sum(df$X3, na.rm = TRUE)

  # ---- Neutral candidate ranking (REFSEL-01; D-01/D-02/D-03) --------------
  # Total reads per subtype (for ordering subtypes), genotype retained for the
  # long-format table and the legacy shim columns.
  subtype_reads <- df %>%
    group_by(Subtype, Genotype) %>%
    summarise(reads = sum(X3), .groups = "drop") %>%
    arrange(desc(reads))

  # Top reference per distinct subtype (most reads within the subtype). This is
  # the distinct-subtype de-duplication (D-02): one ref per subtype, so a
  # same-subtype second reference is never selected. The old validity-filtering
  # function is DELETED here: the validity rules (different-genotype / 1a-1b
  # allow / 2k1b block) moved OUT of selection and INTO Phase 8 classification
  # (D-05). Selection is now mechanical (read-recruitment ranking only).
  top_ref_per_subtype <- df %>%
    group_by(Subtype) %>%
    arrange(desc(X3), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(Subtype, candidate_ref = X1, percent_gt_4_int)

  # Order the per-subtype top refs by the subtype's total reads, then take the
  # top-N subtypes (head(n = n_candidates)). Fewer rows if fewer distinct
  # subtypes exist (single-candidate / sparse samples).
  ranked <- subtype_reads %>%
    left_join(top_ref_per_subtype, by = "Subtype") %>%
    arrange(desc(reads)) %>%
    head(n = n_candidates) %>%
    mutate(candidate_rank = row_number())

  # Per-candidate confirmation_status: the existing reads>minRead && cov>minCov
  # threshold comparison generalized per candidate (Open Q2 resolution). Keeps
  # the legacy major-pass behaviour observable; Phase 8 redefines the vocabulary.
  candidates_long <- ranked %>%
    transmute(
      sample              = sampleName,
      candidate_rank      = as.integer(candidate_rank),
      candidate_ref       = candidate_ref,
      candidate_subtype   = Subtype,
      candidate_genotype  = Genotype,
      candidate_reads     = reads,
      candidate_cov       = percent_gt_4_int,
      confirmation_status = if_else(
        reads > minRead & percent_gt_4_int > minCov, "pass", "below_threshold"
      )
    )

  selected_refs <- candidates_long$candidate_ref

  # ---- Legacy 10-column shim reconstruction (D-06) ------------------------
  # rank 1 -> major_*, rank 2 -> minor_*. All columns always present; minor_*
  # stays NA when there is no 2nd candidate (single-subtype sample).
  rank1 <- candidates_long %>% filter(candidate_rank == 1)
  rank2 <- candidates_long %>% filter(candidate_rank == 2)

  if (nrow(rank1) == 1) {
    df_final$major_ref[1]   <- rank1$candidate_ref[1]
    df_final$major_reads[1] <- rank1$candidate_reads[1]
    df_final$major_cov[1]   <- rank1$candidate_cov[1]
  }

  if (nrow(rank2) == 1) {
    df_final$minor_ref[1]   <- rank2$candidate_ref[1]
    df_final$minor_reads[1] <- rank2$candidate_reads[1]
    df_final$minor_cov[1]   <- rank2$candidate_cov[1]
  }

  # ---- Gate decision (legacy minor_call / gate_flag for the shim) ---------
  # Reconstructed from the same per-candidate threshold comparison so non-gated
  # samples keep their existing minor_call/gate_flag values. The major must pass
  # BOTH thresholds before any minor can be reported.
  major_pass <- nrow(rank1) == 1 &&
                !is.na(df_final$major_reads[1]) &&
                (df_final$major_reads[1] > minRead) &&
                (df_final$major_cov[1] > minCov)
  minor_pass <- nrow(rank2) == 1 &&
                !is.na(df_final$minor_reads[1]) &&
                (df_final$minor_reads[1] > minRead) &&
                (df_final$minor_cov[1] > minCov)

  df_final$minor_call[1] <- if (isTRUE(major_pass) && isTRUE(minor_pass)) "yes" else "no"
  df_final$gate_flag[1]  <- if (!isTRUE(major_pass)) "major_below_threshold" else "ok"
}

# Write the new long-format candidate table (distinct glob from the legacy wide
# CSV so they never collide). Always written, even header-only on no-mapping.
write_csv(candidates_long, file = paste0(sampleName, ".candidates.csv"))

# Write the reconstructed legacy wide CSV (shim).
write_csv(df_final, file = paste0(sampleName, ".parsefirstmapping.csv"))

# ---- FASTA write — uniform per-rank cand{k} slot, guarded (D-06; Pitfall 3) --
# Writes one <sample>.<ref>_cand{k}.fa per selected candidate (k = 1..N), so the
# downstream meta.reference enrichment (fasta basename split) yields <ref>_cand{k}
# consistently with the new config slot. Only write when candidates exist, so the
# no-mapping branch (selected_refs empty) writes nothing and never crashes.
if (length(selected_refs) > 0) {
  fasta <- read.fasta(file = references)
  for (k in seq_along(selected_refs)) {
    ref <- selected_refs[k]
    write.fasta(sequences = fasta[ref], names = ref,
                file.out = paste0(sampleName, ".", ref, "_cand", k, ".fa"))
  }
}
