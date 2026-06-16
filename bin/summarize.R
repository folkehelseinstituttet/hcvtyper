#!/usr/bin/env Rscript

library(tidyverse)

# De novo confirmation helpers (Phase 3). Sourced from the task workdir, where
# they are staged as declared `path` process inputs (plan 03-03), mirroring how
# summarize_mapping_to_all_references.R sources genotype_utils.R relatively.
# genotype_utils.R MUST be sourced first so genotype_from_subtype() is in scope
# before denovo_confirm.R / the confirmation layer use it.
source("genotype_utils.R")
source("denovo_confirm.R")
source("denovo_layer.R")
# Phase 7 (ASUP-02): genotype-level assembly-support join helper. Sourced AFTER
# genotype_utils.R so genotype_from_subtype() is already in scope (the helper
# does not re-source it). Pure sourced helper; defines join_assembly_support().
source("assembly_support_join.R")
# Phase 8 (SCORE-01/02, CLASS-01..04): dominance scoring + strain-role classifier.
# Sourced AFTER genotype_utils.R (relies on genotype_from_subtype() being in scope;
# does not re-source it) and AFTER assembly_support_join.R. Pure sourced helper;
# defines score_candidates() + classify_roles() + is_valid_minor(). This REPLACES
# the consumption of denovo_layer.R / denovo_confirm.R below (D-15 — one
# confirmation system, not two). The source() lines for the legacy helpers are
# retained so the staged files still load cleanly, but apply_denovo_layer() is no
# longer called.
source("classify_roles.R")

args = commandArgs(trailingOnly=TRUE)

# Define variables --------------------------------------------------------

samplesheet      <- args[1]
pipeline_version <- args[2]
pipeline_name    <- args[3]

# De novo confirmation params (D-05): parsed-but-unused in Phase 2. Read defensively
# (cf. contamination_report.R optional-arg pattern). These are NEVER branched on this
# phase; Phase 3 consumes them. Defaults mirror plan 01's nextflow.config defaults.
denovo_min_contig_length  <- if (length(args) >= 4 && nchar(args[4]) > 0) as.numeric(args[4]) else 1000
denovo_min_kmer_cov       <- if (length(args) >= 5 && nchar(args[5]) > 0) as.numeric(args[5]) else 2.0
denovo_min_blast_identity <- if (length(args) >= 6 && nchar(args[6]) > 0) as.numeric(args[6]) else 90
denovo_match_level        <- if (length(args) >= 7 && nchar(args[7]) > 0) args[7] else "genotype"
# WR-04: reject a malformed match level at the source rather than silently
# misinterpreting any non-"subtype" value as "genotype" downstream.
stopifnot(denovo_match_level %in% c("genotype", "subtype"))
denovo_confirm_minor      <- if (length(args) >= 8 && nchar(args[8]) > 0) as.logical(args[8]) else TRUE
# First-mapping candidate prune thresholds; reused as the classify_roles()
# dominance-gate floor (D-07/D-09). GATE-03 (the v1.2.0 dedup-targeted minor
# suppressor) removed per D4 — these are no longer applied to minor_typable here.
min_targeted_read <- if (length(args) >= 9  && nchar(args[9])  > 0) as.numeric(args[9])  else NA_real_
min_targeted_cov  <- if (length(args) >= 10 && nchar(args[10]) > 0) as.numeric(args[10]) else NA_real_
# Number of neutrally-ranked candidate slots (cand_1..cand_n). Used to give the
# Phase-7 assembly-support wide block a FIXED column set (WR-01/WR-02) so the
# Summary.csv schema is stable regardless of which ranks happen to populate in a
# given batch. Default 2 mirrors nextflow.config params.n_candidates.
n_candidates <- if (length(args) >= 11 && nchar(args[11]) > 0) as.integer(args[11]) else 2L

# Phase 8 dominance-score weights (D-04). APPENDED at args[12]+ — strictly AFTER
# n_candidates (args[11]) — because the SUMMARIZE ext.args string is positional
# (conf/modules_hcv.config) and the R side reads args by index; inserting these
# mid-string would silently re-map every later position (Pitfall 2 / T-08-04).
# Read with the same defensive index-guarded form; defaults mirror the
# nextflow.config values declared in Plan 01 (evenness 3.0 > reads 1.0 so
# breadth-evenness dominates raw read count, SCORE-02).
score_weight_evenness <- if (length(args) >= 12 && nchar(args[12]) > 0) as.numeric(args[12]) else 3.0
score_weight_reads    <- if (length(args) >= 13 && nchar(args[13]) > 0) as.numeric(args[13]) else 1.0
score_weight_kmercov  <- if (length(args) >= 14 && nchar(args[14]) > 0) as.numeric(args[14]) else 0.5
# Evenness-transform constant (factor = 1/(1 + k*CV)); plumbed end-to-end so a
# caller supplying a raw CV instead of a precomputed factor stays configurable.
score_evenness_k      <- if (length(args) >= 15 && nchar(args[15]) > 0) as.numeric(args[15]) else 1.0

script_name_version <- if (!is.na(pipeline_version) && nzchar(trimws(pipeline_version))) {
  paste(pipeline_name, pipeline_version)
} else {
  paste(pipeline_name, "(version unknown)")
}

path_1 <- "trimmed/"
path_2 <- "kraken_classified/"
path_3 <- "parsefirst_mapping/"
path_4 <- "stats_withdup/"
path_5 <- "stats_markdup/"
path_6 <- "depth/"
path_denovo <- "denovo/"
path_8 <- "glue/"
path_9 <- "id/"
path_10 <- "variation/"
path_11 <- "consensus_distance/"



# Trimmed ----------------------------------------------------------------

# List files
trimmed_files <- list.files(path = path_1, pattern = ".log$", full.names = TRUE)

# Empty df
trimmed_df <- tibble(
  sampleName = rep(NA_character_, length(trimmed_files)),
  total_raw_reads = rep(NA_real_, length(trimmed_files)),
  total_trimmed_reads = rep(NA_real_, length(trimmed_files))
)

for (i in seq_along(trimmed_files)) {
  f <- trimmed_files[i]
  # Get the sampleName
  trimmed_df$sampleName[i] <- str_split(basename(f), "\\.")[[1]][1]

  if (!file.exists(f) || file.size(f) == 0) {
    warning(glue::glue("Log file missing or empty: {f}"))
    next
  }

  try(rm(trimmed_stats), silent = TRUE)

  # Read the log file from either cutadapt or fastp
  trimmed_stats <- read_lines(f) %>% tibble(value = .)

  # Check if the string fastp is found in the log file. Then process accordingly
  if (any(str_detect(trimmed_stats$value, "fastp v"))) {
    ## ---- fastp logs ----

    # Raw reads = Read1 before + Read2 before
    raw1_idx <- which(str_detect(trimmed_stats$value, "Read1 before filtering:")) # Finds the line number
    raw2_idx <- which(str_detect(trimmed_stats$value, "Read2 before filtering:"))

    raw1_val <- if (length(raw1_idx) > 0) {
      as.numeric(str_extract(trimmed_stats$value[raw1_idx + 1], "\\d+")) # Pull out the number on the line after
    } else NA
    raw2_val <- if (length(raw2_idx) > 0) {
      as.numeric(str_extract(trimmed_stats$value[raw2_idx + 1], "\\d+"))
    } else NA

    if (!is.na(raw1_val) && !is.na(raw2_val)) {
      trimmed_df$total_raw_reads[i] <- raw1_val + raw2_val
    }

    # Trimmed reads = "reads passed filter:"
    passed <- trimmed_stats %>%
      filter(str_detect(value, "^reads passed filter:"))
    if (nrow(passed) > 0) {
      tmp <- as.numeric(str_extract(passed$value, "\\d+"))
      trimmed_df$total_trimmed_reads[i] <- tmp
    }

  } else {
    ## ---- Cutadapt logs ----
    raw <- filter(trimmed_stats, str_detect(value, "Total read pairs processed:"))
    if (nrow(raw) > 0) {
      tmp <- str_split(raw$value, "\\s+")[[1]] # Split on white space and get the list content
      tmp <- as.numeric(str_remove_all(tmp[length(tmp)], ",")) # extract the last element which contains the pair number, remove commas and create a numeric
      tmp <- tmp * 2 # Double to get reads
      trimmed_df$total_raw_reads[i] <- tmp
    }

    trimmed <- filter(trimmed_stats, str_detect(value, "Pairs written \\(passing filters\\)"))
    if (nrow(trimmed) > 0) {
      tmp <- str_split(trimmed$value, "\\s+")[[1]]
      tmp <- as.numeric(str_remove_all(tmp[length(tmp)-1], ","))
      tmp <- tmp * 2
      trimmed_df$total_trimmed_reads[i] <- tmp
    }
  }
}

trimmed_df <- as_tibble(trimmed_df)

# Kraken ------------------------------------------------------------------
# List files
kraken_files <- list.files(path = path_2, pattern = "kraken2.report.txt$", full.names = TRUE)

# Empty df
kraken_df <- as.data.frame(matrix(nrow = length(kraken_files), ncol = 2))
colnames(kraken_df) <- c("sampleName", "total_classified_reads")

for (i in 1:length(kraken_files)) {
  try(rm(kraken_stats))
  # Get sample name
  kraken_df$sampleName[i] <- str_split(basename(kraken_files[i]), "\\.")[[1]][1]

  # Get total of trimmed sequences put in to the mapping. Adding try() if no root sequences
  try(kraken_stats <- read_tsv(kraken_files[i], col_names = FALSE) %>% filter(X6 == "root") %>% pull(X2))
  if (exists("kraken_stats") & length(kraken_stats) > 0) {
    kraken_df$total_classified_reads[i] <- kraken_stats*2 # Kraken reports read pairs (fragments)
  }
}
kraken_df <- as_tibble(kraken_df)

# Total mapped reads to all references, with duplicates -------------------
first_mapping_files <- list.files(path = path_3, pattern = "parsefirstmapping.csv$", full.names = TRUE)

# Empty df
parsefirstmapping_df <- tibble(
  sampleName = rep(NA_character_, length(first_mapping_files)),
  total_mapped_reads = rep(NA_real_, length(first_mapping_files)),
  major_mapped_reads = rep(NA_real_, length(first_mapping_files)),
  minor_mapped_reads = rep(NA_real_, length(first_mapping_files)),
  major_cov_firstmapping = rep(NA_real_, length(first_mapping_files)),
  major_ref_firstmapping = rep(NA_character_, length(first_mapping_files)),
  gate_flag = rep(NA_character_, length(first_mapping_files))
)

# If the length of parsefirstmapping_files is non-zero
if (length(first_mapping_files) > 0) {
  for (i in 1:length(first_mapping_files)) {
    try(rm(sample_parsefirstmapping))
    # Get sample name
    parsefirstmapping_df$sampleName[i] <- str_split(basename(first_mapping_files[i]), "\\.")[[1]][1]

    # Read the mapping stats
    sample_parsefirstmapping <- read_csv(first_mapping_files[i])

    # Get number of mapped reads before duplicate removal
    parsefirstmapping_df$total_mapped_reads[i] <- sample_parsefirstmapping %>% pull(total_mapped_reads)

    # Get the number of mapped reads against all major references belonging to the major subtype
    parsefirstmapping_df$major_mapped_reads[i] <- sample_parsefirstmapping %>% pull(major_reads)

    # Get the number of mapped reads against all minor references belonging to the minor subtype
    parsefirstmapping_df$minor_mapped_reads[i] <- sample_parsefirstmapping %>% pull(minor_reads)

    # If present, capture coverage, reference and gate_flag from the first-mapping report
    if ("major_cov" %in% colnames(sample_parsefirstmapping)) {
      parsefirstmapping_df$major_cov_firstmapping[i] <- sample_parsefirstmapping %>% pull(major_cov)
    }
    if ("major_ref" %in% colnames(sample_parsefirstmapping)) {
      parsefirstmapping_df$major_ref_firstmapping[i] <- sample_parsefirstmapping %>% pull(major_ref)
    }
    if ("gate_flag" %in% colnames(sample_parsefirstmapping)) {
      parsefirstmapping_df$gate_flag[i] <- sample_parsefirstmapping %>% pull(gate_flag)
    }
  }
}

parsefirstmapping_df <- as_tibble(parsefirstmapping_df) %>%
  # Calculate the median of total mapped reads accross all samples, and then for each sample the fraction of mapped reads compared to the median
  # In addition calculate the fraction of mapped reads against major and minor references of the total mapped reads
  mutate(
    median_mapped = median(total_mapped_reads),
    fraction_mapped_reads_vs_median = total_mapped_reads / median_mapped,
    percent_mapped_reads_major_firstmapping = round(major_mapped_reads / total_mapped_reads * 100, digits = 2),
    percent_mapped_reads_minor_firstmapping = round(minor_mapped_reads / total_mapped_reads * 100, digits = 2)
  ) %>%
  select(sampleName, total_mapped_reads, fraction_mapped_reads_vs_median, percent_mapped_reads_major_firstmapping, percent_mapped_reads_minor_firstmapping, major_cov_firstmapping, major_ref_firstmapping, gate_flag)

# Phase-9 (COMPAT-02 / D-02): load the long-format candidate table BEFORE the three
# stats loops so each loop can recover its candidate rank by joining its per-file
# reference token to candidate_rank, instead of parsing a major/minor slot from
# filename position 3. The load depends only on path_3 (set near the top), so it is
# safe to hoist here ahead of the stats loops (originally loaded near the Phase-7
# assembly-support join). The Phase-7 assembly-support join below reuses this frame.
candidates_files <- list.files(path = path_3, pattern = "\\.candidates.csv$", full.names = TRUE)

if (length(candidates_files) > 0) {
  # CR-01/CR-02: pin column types so (a) a header-only candidates.csv (no-mapping
  # branch) and a populated one combine cleanly under map_dfr/bind_rows, and (b)
  # candidate_genotype is read as character — HCV genotypes 1–7 are purely-digit,
  # which readr would otherwise infer as <double>, breaking the genotype-level
  # join. col_types declares only the columns we depend on; the rest are inferred.
  candidates_long <- map_dfr(candidates_files, ~ read_csv(.x, col_types = cols(
    sample              = col_character(),
    candidate_rank      = col_integer(),
    candidate_ref       = col_character(),
    candidate_subtype   = col_character(),
    candidate_genotype  = col_character(),
    candidate_reads     = col_double(),
    candidate_cov       = col_double(),
    confirmation_status = col_character()
  ))) %>%
    rename(sampleName = sample)
} else {
  # Declare all eight Phase-6 candidate columns with their types so a no-candidate
  # run yields a typed zero-row frame (the join then returns a typed zero-row frame).
  candidates_long <- tibble(
    sampleName          = character(),
    candidate_rank      = integer(),
    candidate_ref       = character(),
    candidate_subtype   = character(),
    candidate_genotype  = character(),
    candidate_reads     = double(),
    candidate_cov       = double(),
    confirmation_status = character()
  )
}

# Phase-9 (COMPAT-02 / D-02): per-sample (candidate_ref -> candidate_rank) lookup
# used by all three stats loops to recover the slot by join. Built once from the
# hoisted candidates_long frame.
candidate_rank_lookup <- candidates_long %>%
  select(sampleName, candidate_ref, candidate_rank) %>%
  distinct(sampleName, candidate_ref, .keep_all = TRUE)

# Second mapping, reads mapped with duplicates ----------------------------
# List files
stats_files <- list.files(path = path_4, pattern = "\\withdup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 3))
colnames(tmp_df) <- c("sampleName", "reference", "trimmed_reads_withdups_mapped")

for (i in 1:length(stats_files)) {
  try(rm(map_stats))
  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][2]

  # Phase-9 (COMPAT-02 / D-02): the candidate rank is no longer parsed from
  # filename position 3 (.major./.minor./.cand{rank}.). It is recovered by joining
  # the cleaned reference token to candidate_rank below.

  # Read the mapping stats
  map_stats <- read_tsv(stats_files[i], col_names = FALSE, comment = "#")

  # Get number of mapped reads before duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)

  mapped_reads <- as.numeric(mapped_reads)
  tmp_df$trimmed_reads_withdups_mapped[i] <- mapped_reads
}
tmp_df <- as_tibble(tmp_df)

# Phase-9 (COMPAT-02 / D-02): recover candidate_rank by join. Strip the new
# `_cand{rank}` slot suffix off the reference token to get the cleaned candidate_ref
# (e.g. `3a_D17763`), then join to the per-sample candidate_rank_lookup.
tmp_df <- tmp_df %>%
  mutate(candidate_ref = str_remove(reference, "_cand[0-9]+$")) %>%
  left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref"))

# Add number of raw and trimmed reads - needed for calculation of percentages
tmp_df <- left_join(tmp_df, trimmed_df, by = "sampleName")

# Add number of classified reads from Kraken2 -  - needed for calculation of percentages
tmp_df <- left_join(tmp_df, kraken_df, by = "sampleName")

df_with_dups <- tmp_df %>%
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_withdup_mapped_major = case_when(candidate_rank == 1 ~ trimmed_reads_withdups_mapped)) %>%
  mutate(Reads_withdup_mapped_minor = case_when(candidate_rank == 2 ~ trimmed_reads_withdups_mapped)) %>%
  # Don't include number of reads mapped in the first mapping. Info must be taken from another process if we should include
  #mutate(Reads_withdup_mapped_first_mapping = case_when(reference == "first_mapping" ~ trimmed_reads_withdups_mapped)) %>%
  select(-trimmed_reads_withdups_mapped) %>%
  # Create columns for the major and minor references (cleaned candidate_ref — no
  # slot suffix — so it stays byte-identical to the legacy join key, COMPAT-02 D-02).
  mutate(Major_reference = case_when(candidate_rank == 1 ~ candidate_ref)) %>%
  mutate(Minor_reference = case_when(candidate_rank == 2 ~ candidate_ref)) %>%
  select(-reference, -candidate_ref) %>%
  # Calculate percent of the trimmed reads mapped
  #mutate(total_trimmed_reads_with_dups = as.integer(total_trimmed_reads_with_dups),
  #       Reads_withdup_mapped_major = as.integer(Reads_withdup_mapped_major),
  #       Reads_withdup_mapped_minor = as.integer(Reads_withdup_mapped_minor)) %>%
         #Reads_withdup_mapped_first_mapping = as.integer(Reads_withdup_mapped_first_mapping)) %>%
  mutate(Percent_reads_mapped_of_trimmed_with_dups_major = Reads_withdup_mapped_major / total_trimmed_reads * 100,
         Percent_reads_mapped_of_trimmed_with_dups_minor = Reads_withdup_mapped_minor / total_trimmed_reads * 100) %>%
         #Percent_reads_mapped_with_dups_first_mapping = Reads_withdup_mapped_first_mapping / total_trimmed_reads_with_dups * 100) %>%
  # Create one row per sample
  select(-candidate_rank) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Reads mapped no duplicates ----------------------------------------------
# List files
stats_files <- list.files(path = path_5, pattern = "nodup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 3))
colnames(tmp_df) <- c("sampleName", "reference", "trimmed_reads_nodups_mapped")

for (i in 1:length(stats_files)) {
  try(rm(mapped_reads))

  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][2]

  # Phase-9 (COMPAT-02 / D-02): candidate rank recovered by join below, not from
  # filename position 3.

  # Read the mapping stats
  map_stats <- read_tsv(stats_files[i], col_names = FALSE, comment = "#")

  # Get number of mapped reads after duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)

  mapped_reads <- as.numeric(mapped_reads)
  tmp_df$trimmed_reads_nodups_mapped[i] <- mapped_reads

}
tmp_df <- as_tibble(tmp_df)

# Phase-9 (COMPAT-02 / D-02): recover candidate_rank by join. Strip the `_cand{rank}`
# slot suffix to get candidate_ref; the `first_mapping` reference has no suffix and
# no candidate_rank, so it survives the join with candidate_rank == NA and is matched
# below by `reference == "first_mapping"`.
tmp_df <- tmp_df %>%
  mutate(candidate_ref = str_remove(reference, "_cand[0-9]+$")) %>%
  left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref"))

df_nodups <- tmp_df %>%
  # Create columns for major and minor
  separate(candidate_ref, into = c("genotype", NA), sep = "_", remove = F) %>%
  mutate(Major_genotype_mapping = case_when(candidate_rank == 1 ~ genotype)) %>%
  mutate(Minor_genotype_mapping = case_when(candidate_rank == 2 ~ genotype)) %>%
  select(-genotype) %>%
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_nodup_mapped_major = case_when(candidate_rank == 1 ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_minor = case_when(candidate_rank == 2 ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_first_mapping = case_when(reference == "first_mapping" ~ trimmed_reads_nodups_mapped)) %>%
  select(-trimmed_reads_nodups_mapped) %>%
  #mutate(Percent_mapped_major = case_when(candidate_rank == 1 ~ Percent_trimmed_reads_mapped)) %>%
  #mutate(Percent_mapped_minor = case_when(candidate_rank == 2 ~ Percent_trimmed_reads_mapped)) %>%
  # Create columns for the major and minor references (cleaned candidate_ref — no
  # slot suffix — byte-identical to the legacy join key, COMPAT-02 D-02).
  mutate(Major_reference = case_when(candidate_rank == 1 ~ candidate_ref)) %>%
  mutate(Minor_reference = case_when(candidate_rank == 2 ~ candidate_ref)) %>%
  # Create one row per sample
  select(-reference, -candidate_ref, -candidate_rank) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Combine mapped reads data
df_mapped_reads <- full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Remove columns for total_raw_reads, total_trimmed_reads and total_classified_reads.
  # These will be added later to ensure info is kept for samples that were filtered out before the second mapping
  select(-total_raw_reads, -total_trimmed_reads, -total_classified_reads)

# Coverage ----------------------------------------------------------------

# Add both breadth (in percent) and depth (average depth)
# All this is without duplicates

# List files
cov_files <- list.files(path = path_6, pattern = "tsv$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(cov_files), ncol = 7))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_1", "cov_breadth_min_5", "cov_breadth_min_10", "avg_depth", "cv_evenness")

for (i in 1:length(cov_files)) {
  try(rm(cov))

  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][2]

  # Phase-9 (COMPAT-02 / D-02): candidate rank recovered by join below, not from
  # filename position 3.

  # Read the depth per position
  cov <- read_tsv(cov_files[i], col_names = FALSE)

  # Reference length
  ref_length <- nrow(cov)

  # Average depth
  tmp_df$avg_depth[i] <- mean(cov$X3)

  # CV-evenness factor (Phase 8, D-03). SAMTOOLS_DEPTH runs with `-aa`
  # (conf/modules_hcv.config:259) so cov$X3 spans EVERY reference position incl.
  # zeros — the coefficient of variation (sd/mean) over this full vector captures
  # how spiky/uneven the per-position depth is. Map it to a 0–1 evenness factor
  # 1/(1+CV): a perfectly flat pileup -> 1, a spiky cross-mapping artefact -> ~0.
  # The factor is computed HERE because cov$X3 exists only inside this loop and is
  # discarded once the loop iterates. Guard the zero-mean / zero-length edge
  # (no reads mapped, or empty depth file) -> cv_evenness = 0, never NaN/Inf
  # (Pitfall 3 / T-08-06).
  cv_raw <- if (ref_length > 0 && mean(cov$X3) > 0) sd(cov$X3) / mean(cov$X3) else NA_real_
  tmp_df$cv_evenness[i] <- if (!is.na(cv_raw)) 1 / (1 + cv_raw) else 0

  # Nr. of positions with coverage >=1, >= 5 and > 9
  # If ref_length is zero it means that no reads were mapped. Set coverage to zero.
  # Coverage may also be zero if there are reads mapped, but never more than 5 per position
  if (ref_length > 0) {
    pos_1 <- nrow(
      cov %>%
        filter(X3 >= 1)
    )

    pos_5 <- nrow(
      cov %>%
        filter(X3 >= 5)
    )

    pos_10 <- nrow(
      cov %>%
        filter(X3 > 9)
    )
    # Coverage breadth
    breadth_1 <- round(pos_1 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_1[i] <- breadth_1

    breadth_5 <- round(pos_5 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_5[i] <- breadth_5

    breadth_10 <- round(pos_10 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_10[i] <- breadth_10
  } else if (ref_length == 0) {
    tmp_df$cov_breadth_min_1[i] <- 0
    tmp_df$cov_breadth_min_5[i] <- 0
    tmp_df$cov_breadth_min_10[i] <- 0
  }
}

# Create column for subtype and Sample_ref
tmp_df <- as_tibble(tmp_df)

# Phase-9 (COMPAT-02 / D-02): recover candidate_rank by join. Strip the `_cand{rank}`
# slot suffix to get candidate_ref; `first_mapping` has no suffix / no candidate_rank
# and is filtered out below before any rank-based logic runs.
tmp_df <- tmp_df %>%
  mutate(candidate_ref = str_remove(reference, "_cand[0-9]+$")) %>%
  left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref"))

# Phase 8 (D-03): per-reference cv_evenness lookup for the role classifier.
# df_coverage (below) collapses to one row/sample with Major_/Minor_ slots, which
# loses the per-candidate granularity the classifier needs. Build a long lookup
# keyed by sampleName + the cleaned mapping-reference name (the cov-loop `reference`
# carries a `_major`/`_minor` targeted-mapping suffix; strip it so it matches the
# Phase-6 `candidate_ref` token like `3a_D17763`). This is left_joined into the
# long candidate_support frame so score_candidates() can read cv_evenness per
# candidate. A candidate whose reference was never targeted-mapped (no depth file)
# NA-fills here and score_candidates() treats the missing factor as neutral 0.
cv_by_ref <- tmp_df %>%
  filter(reference != "first_mapping") %>%
  # Phase-9 (COMPAT-02 / D-02): the cov-loop `reference` carries the new `_cand{rank}`
  # targeted-mapping suffix; strip it (identical regex to the stats-loop joins above)
  # so candidate_ref matches the Phase-6 token like `3a_D17763`. tmp_df already holds
  # this column from the coverage candidate_rank join, but re-derive it here so the
  # strip is explicit and cannot drift from the other sites.
  mutate(candidate_ref = str_remove(reference, "_cand[0-9]+$")) %>%
  select(sampleName, candidate_ref, cv_evenness) %>%
  filter(!is.na(candidate_ref)) %>%
  distinct(sampleName, candidate_ref, .keep_all = TRUE)

df_coverage <- tmp_df %>%
  # Don't need first mapping data
  filter(reference != "first_mapping") %>%
  # Create columns for major and minor coverage
  mutate(Major_cov_breadth_min_1 = case_when(candidate_rank == 1 ~ cov_breadth_min_1)) %>%
  mutate(Minor_cov_breadth_min_1 = case_when(candidate_rank == 2 ~ cov_breadth_min_1)) %>%
  mutate(Major_cov_breadth_min_5 = case_when(candidate_rank == 1 ~ cov_breadth_min_5)) %>%
  mutate(Minor_cov_breadth_min_5 = case_when(candidate_rank == 2 ~ cov_breadth_min_5)) %>%
  mutate(Major_cov_breadth_min_10 = case_when(candidate_rank == 1 ~ cov_breadth_min_10)) %>%
  mutate(Minor_cov_breadth_min_10 = case_when(candidate_rank == 2 ~ cov_breadth_min_10)) %>%
  # Create columns for major and minor average depth
  mutate(Major_avg_depth = case_when(candidate_rank == 1 ~ avg_depth)) %>%
  mutate(Minor_avg_depth = case_when(candidate_rank == 2 ~ avg_depth)) %>%
  # Create columns for the major and minor references (cleaned candidate_ref — no
  # slot suffix — byte-identical to the legacy join key, COMPAT-02 D-02).
  mutate(Major_reference = case_when(candidate_rank == 1 ~ candidate_ref)) %>%
  mutate(Minor_reference = case_when(candidate_rank == 2 ~ candidate_ref)) %>%
  # Create one row per sample
  select(-reference, -candidate_ref, -candidate_rank, -cov_breadth_min_1, -cov_breadth_min_5, -cov_breadth_min_10, -avg_depth) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)


# De novo evidence --------------------------------------------------------

# Read the parsed de novo / BLAST evidence emitted by BLASTPARSE (*.blastparse.csv).
# Columns (cf. blast_parse.R summary CSV): sample, major_ref, major_contig_length,
# minor_ref, minor_contig_length. The `sample` column == prefix == meta.id; trust it
# as the join key (renamed to sampleName). Guarded on length(...) > 0 so a
# skip-assembly / no-de-novo run yields a typed empty tibble -> NA fields, never abort
# (T-02-02). The *_blast_out.csv per-contig table is plumbed THROUGH the denovo/
# staging dir for Phase 3 only (Open Q2 / D-05) and is NOT read or aggregated here.
blastparse_files <- list.files(path = path_denovo, pattern = "blastparse.csv$", full.names = TRUE)

if (length(blastparse_files) > 0) {
  df_denovo <- map_dfr(blastparse_files, read_csv) %>%
    rename(
      sampleName                 = sample,
      denovo_major_ref           = major_ref,
      denovo_major_contig_length = major_contig_length,
      denovo_minor_ref           = minor_ref,
      denovo_minor_contig_length = minor_contig_length
    )
} else {
  # PLUMB-02: a skip-assembly / no-de-novo run must still yield the four denovo_*
  # columns NA-filled (never drop them). An empty tibble carrying ONLY sampleName
  # contributes no columns to the left_join, so the denovo_* columns would vanish
  # from Summary.csv on skip-assembly. Declare all four columns with their
  # blast_parse.R types (ref = character, contig_length = integer) so the join
  # always emits them; with zero rows here every existing sample row NA-fills.
  df_denovo <- tibble(
    sampleName                 = character(),
    denovo_major_ref           = character(),
    denovo_major_contig_length = integer(),
    denovo_minor_ref           = character(),
    denovo_minor_contig_length = integer()
  )
}

# Per-contig de novo BLAST table (Phase 3 confirmation, CONF-01/02/03). Read ALL
# denovo/*_blast_out.csv ONCE into a long frame keyed by sampleName (basename
# prefix before "_blast_out.csv"); the per-sample confirmation scan filters this
# by sampleName at the chokepoint. Schema mirrors blast_parse.R (qseqid, subtype,
# pident, evalue, bitscore, sc_length, kmer_cov). Guarded with length(...) > 0 so
# a skip-assembly / no-de-novo run yields a typed-empty tibble and NEVER aborts
# (T-03-01 DoS guard); a sample absent from this frame classifies as "unconfirmed",
# never "refuted". This is read-only staging — NOT joined into `final` here.
blast_out_files <- list.files(path = path_denovo, pattern = "_blast_out.csv$", full.names = TRUE)

if (length(blast_out_files) > 0) {
  df_blast_out <- map_dfr(blast_out_files, ~ read_csv(.x, show_col_types = FALSE) %>%
    mutate(sampleName = str_remove(basename(.x), "_blast_out.csv$")))
} else {
  df_blast_out <- tibble(
    sampleName = character(),
    qseqid     = character(),
    subtype    = character(),
    pident     = double(),
    evalue     = double(),
    bitscore   = double(),
    sc_length  = double(),
    kmer_cov   = double()
  )
}

# Phase-7 per-subtype assembly support (ASUP-02). The long-format candidate table
# (candidates_long) is loaded earlier — before the three stats loops — so those
# loops can recover their candidate rank by join (Phase-9 / COMPAT-02 / D-02). The
# typed-empty-tibble guard there mirrors the PLUMB-02 pattern so a skip-assembly /
# no-candidate run NA-fills the new support columns and never aborts (T-07-03 DoS
# guard). candidates_long is the LEFT side of the genotype-level join (criterion #3,
# no row loss); support_df is the per-subtype RIGHT side.
support_files <- list.files(path = path_denovo, pattern = "\\.assembly_support.csv$", full.names = TRUE)

if (length(support_files) > 0) {
  # CR-02: pin column types so a header-only assembly_support.csv (skip-assembly /
  # no-contig sample) and a populated one from another sample combine cleanly
  # under map_dfr/bind_rows. readr types every column of a header-only CSV as
  # <character>, while a populated CSV types the metric columns as <double>;
  # without col_types the bind_rows across them aborts SUMMARIZE.
  support_df <- map_dfr(support_files, ~ read_csv(.x, col_types = cols(
    sample                 = col_character(),
    subtype                = col_character(),
    best_contig_length     = col_double(),
    best_contig_pident     = col_double(),
    best_contig_aln_length = col_double(),
    best_contig_kmer_cov   = col_double()
  ))) %>%
    rename(sampleName = sample)
} else {
  # Declare the six Plan-01 assembly-support columns (+ sampleName) with their
  # types so a skip-assembly run yields a zero-row frame -> every candidate NA-fills.
  support_df <- tibble(
    sampleName             = character(),
    subtype                = character(),
    best_contig_length     = double(),
    best_contig_pident     = double(),
    best_contig_aln_length = double(),
    best_contig_kmer_cov   = double()
  )
}

# Join assembly support to candidates at the parameterized denovo_match_level
# (default genotype). Candidates anchor the LEFT side (no row loss). The per-
# candidate support columns are pivoted to wide cand_<rank> slots below so they
# can left_join onto `final` (one row/sample) without exploding rows.
candidate_support <- join_assembly_support(candidates_long, support_df, denovo_match_level)

# Phase 8 dominance scoring + role classification (SCORE-01/02, CLASS-01..04, D-15).
# Attach the per-candidate cv_evenness factor computed in the cov loop (joined on
# sampleName + candidate_ref; a candidate whose reference was never targeted-mapped
# NA-fills and score_candidates() treats it as neutral 0), then run the pure
# classifier: score_candidates() emits dominance_score (breadth-evenness dominating
# raw reads), classify_roles() emits role / role_reason / overall_sample_call per
# candidate. This is the N-candidate role model that REPLACES the legacy
# apply_denovo_layer() / minor_denovo_status / coinfection_flag path retired below.
# left_join (not inner) so a no-cov / no-candidate batch keeps every candidate row
# (the classifier's typed zero-row guard handles the empty frame, T-08-01/CLASS-03).
candidate_support <- candidate_support %>%
  left_join(cv_by_ref, by = c("sampleName", "candidate_ref"))

candidate_support <- score_candidates(
  candidate_support,
  score_weights  = list(
    evenness = score_weight_evenness,
    reads    = score_weight_reads,
    kmercov  = score_weight_kmercov
  ),
  evenness_const = score_evenness_k
)

candidate_support <- classify_roles(
  candidate_support,
  # D-07/D-09: the dominant major-gate and the co-infection floor share ONE
  # threshold set. The SUMMARIZE ext.args passes ${params.minRead}/${params.minCov}
  # at positions 6/7 -> args[9]/args[10], parsed above as min_targeted_read /
  # min_targeted_cov. classify_roles() applies `reads > minRead & cov > minCov`,
  # so the gate floor is exactly the configured minRead/minCov (NA -> treated as a
  # non-passing gate via the NA guard inside classify_roles()).
  minRead = min_targeted_read,
  minCov  = min_targeted_cov,
  denovo_min_contig_length  = denovo_min_contig_length,
  denovo_min_kmer_cov       = denovo_min_kmer_cov,
  denovo_min_blast_identity = denovo_min_blast_identity,
  match_level               = denovo_match_level
)

# Enriched long *.candidates.csv (D-16, CLASS-03). Write EVERY candidate — incl.
# background / refuted — carrying the new role / dominance_score / role_reason +
# overall_sample_call alongside the original Phase-6 candidate columns and the
# Phase-7 assembly-support join. Backgrounds are surfaced here, never dropped; the
# wide Summary.csv below only carries the dominant + corroborated co-infection
# slots. Empty-batch safe: classify_roles() returns a typed zero-row frame, so a
# no-candidate run writes a header-only candidates.csv (never aborts, T-08-01).
write_csv(candidate_support, file = "candidates.csv")

# Map the per-candidate roles into the wide one-row-per-sample layout (D-16). The
# `dominant` candidate fills the role_* major slot and the first corroborated
# `co-infection` candidate fills the role_* minor slot; `overall_sample_call` is a
# new wide column. These role_* columns sit ALONGSIDE the legacy Major_*/Minor_*
# mapping/coverage slots (the legacy Major_*/Minor_* column ALIASING onto the role
# model is Phase 9 / COMPAT-03 — Phase 8 only retires the legacy confirmation
# LOGIC, D-15). overall_sample_call is taken per-sample (constant within a sample).
if (nrow(candidate_support) > 0) {
  role_dominant <- candidate_support %>%
    filter(role == "dominant") %>%
    group_by(sampleName) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      sampleName,
      Major_role_reference     = candidate_ref,
      Major_role_subtype       = candidate_subtype,
      Major_dominance_score    = dominance_score,
      Major_role_reason        = role_reason
    )

  role_minor <- candidate_support %>%
    filter(role == "co-infection") %>%
    group_by(sampleName) %>%
    # Deterministic: the highest-scoring corroborated co-infection fills the minor
    # slot (ties already broken inside classify_roles()'s dominant selection; here
    # we just take the top remaining co-infection by dominance_score then ref name).
    arrange(desc(dominance_score), candidate_ref, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      sampleName,
      Minor_role_reference     = candidate_ref,
      Minor_role_subtype       = candidate_subtype,
      Minor_dominance_score    = dominance_score,
      Minor_role_reason        = role_reason
    )

  overall_call <- candidate_support %>%
    distinct(sampleName, overall_sample_call)

  candidate_roles_wide <- overall_call %>%
    left_join(role_dominant, by = "sampleName") %>%
    left_join(role_minor, by = "sampleName")
} else {
  # No-candidate batch: typed zero-row wide frame so the left_join onto `final`
  # below still emits every role_* column (NA-filled per sample), never drops them.
  candidate_roles_wide <- tibble(
    sampleName               = character(),
    overall_sample_call      = character(),
    Major_role_reference     = character(),
    Major_role_subtype       = character(),
    Major_dominance_score    = double(),
    Major_role_reason        = character(),
    Minor_role_reference     = character(),
    Minor_role_subtype       = character(),
    Minor_dominance_score    = double(),
    Minor_role_reason        = character()
  )
}

# WR-01/WR-02: the wide assembly-support block must carry a FIXED column set —
# cand_1..cand_{n_candidates} × the six support values — regardless of which
# ranks actually populate in a given batch (otherwise Summary.csv's schema drifts
# run to run, and a no-candidate batch silently drops the whole block). Build the
# full schema explicitly, then fill it from whatever candidate_support holds.
support_value_cols <- c(
  "assembly_support",
  "assembly_support_subtype",
  "assembly_support_best_contig_length",
  "assembly_support_best_contig_pident",
  "assembly_support_best_contig_aln_length",
  "assembly_support_best_contig_kmer_cov"
)
# Column types per support value, in support_value_cols order: the two status/
# subtype columns are character, the four metrics are double.
support_value_is_character <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)

# The complete, deterministic set of wide column names (cand_<rank>_<value>) for
# ranks 1..n_candidates. names_glue below emits "cand_{rank}_{value}".
expected_wide_cols <- as.vector(t(outer(
  seq_len(n_candidates),
  support_value_cols,
  function(rk, val) paste0("cand_", rk, "_", val)
)))
expected_wide_is_character <- as.vector(t(outer(
  seq_len(n_candidates),
  support_value_is_character,
  function(rk, is_chr) is_chr
)))

# A zero-row tibble carrying sampleName + every expected wide column at its type.
# Used both as the empty-branch frame and to back-fill any rank columns missing
# from the pivot (e.g. a batch where no sample reached rank 2).
empty_wide <- tibble(sampleName = character())
for (i in seq_along(expected_wide_cols)) {
  empty_wide[[expected_wide_cols[i]]] <-
    if (expected_wide_is_character[i]) character() else double()
}

if (nrow(candidate_support) > 0) {
  candidate_support_wide <- candidate_support %>%
    select(
      sampleName,
      candidate_rank,
      all_of(support_value_cols)
    ) %>%
    pivot_wider(
      id_cols     = sampleName,
      names_from  = candidate_rank,
      names_glue  = "cand_{candidate_rank}_{.value}",
      values_from = all_of(support_value_cols)
    )
  # Add any expected cand_<rank>_* columns that this batch's ranks did not produce
  # (WR-02), typed to match empty_wide so the schema is identical run to run.
  missing_wide <- setdiff(expected_wide_cols, names(candidate_support_wide))
  for (col in missing_wide) {
    candidate_support_wide[[col]] <- empty_wide[[col]][NA_integer_][seq_len(nrow(candidate_support_wide))]
  }
  # Reorder to the canonical sampleName + expected-column order.
  candidate_support_wide <- candidate_support_wide %>%
    select(sampleName, all_of(expected_wide_cols))
} else {
  # WR-01: a no-candidate batch must still contribute the full (zero-row) column
  # block so the left_join below never drops the assembly-support columns.
  candidate_support_wide <- empty_wide
}

# GLUE --------------------------------------------------------------------

glue_file <- list.files(path = path_8, pattern = "GLUE_collected_report_major.tsv$", full.names = TRUE)
# Guard the read so an empty glue/ dir (ch_glue -> [], D-07 caveat) does not abort.
# nrow(tibble()) == 0 reproduces the GLUE-absent branch exactly, leaving the existing
# `if (nrow(glue_report) > 0)` guards inert when GLUE is present (PLUMB-04, T-02-03).
glue_report <- if (length(glue_file) > 0) read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character())) else tibble()

# Collect also the minor GLUE report
glue_file_minor <- list.files(path = path_8, pattern = "GLUE_collected_report_minor.tsv$", full.names = TRUE)
glue_report_minor <- if (length(glue_file_minor) > 0) read_tsv(glue_file_minor, col_types = cols(GLUE_subtype = col_character())) else tibble()

# Extract the GLUE genotypes and subtypes for major and minor and compare them

if (nrow(glue_report) > 0) {
  major_gt <- glue_report %>%
    select(Sample, GLUE_genotype, GLUE_subtype) %>%
    rename(Major_genotype = GLUE_genotype,
           Major_subtype = GLUE_subtype)
}

if (nrow(glue_report_minor) > 0) {
  minor_gt <- glue_report_minor %>%
    select(Sample, GLUE_genotype, GLUE_subtype) %>%
    rename(Minor_genotype = GLUE_genotype,
           Minor_subtype = GLUE_subtype)
}

# Check if the same genotype has been called for major and minor. If Yes, then minor is not typable
if (exists("major_gt") & exists("minor_gt")) {
  gt_check <- major_gt %>%
    left_join(minor_gt, by = "Sample") %>%
    mutate(
      identical_geno = case_when(
        Major_genotype == Minor_genotype ~ "YES",
        is.na(Minor_genotype) ~ NA,
        .default = "NO"
      ),
      identical_subgeno = case_when(
        Major_subtype == Minor_subtype ~ "YES",
        is.na(Minor_subtype) ~ NA,
        .default = "NO"
      )
    )
}

# Sequencer ID ------------------------------------------------------------
id_files <- list.files(path = path_9, pattern = "sequencerID.tsv$", full.names = TRUE)

# Initialize id_df unconditionally so the tibble conversion at line 931 never
# references an undefined variable when id_files is empty.
# (WR-02: "1:0" anti-pattern yielded c(1L,0L) and crashed with "object id_df not found")
id_df <- as.data.frame(matrix(nrow = length(id_files), ncol = 2))
colnames(id_df) <- c("sampleName", "sequencer_id")

for (i in seq_along(id_files)) {
  try(rm(id))

  # Get sample name
  id_df$sampleName[i] <- str_split(basename(id_files[i]), "\\.")[[1]][1]

  # Get the sequencer id
  id <- read_tsv(id_files[i], col_names = FALSE)

  # Check if file is empty (i.e., there were no FASTQ reads for this sample)
  if (nrow(id) == 0 || ncol(id) == 0) {
    # Set sequencer ID to NA
    id_df$sequencer_id[i] <- NA
  } else if (str_detect(id$X1, "^@SRR")) {
    # Extract the first field if header start with '@SRR'
    id_df$sequencer_id[i] <- id %>% pull(X1)
  } else {
    id_df$sequencer_id[i] <- id %>%
      pull(X1) %>%
      # Extract string up to the first ":".
      # The "?" means a "lazy", or non-greedy, match to get the shortest string that satisfies the criteria.
      # This is useful because there are several ":"
      str_extract("^.*?:") %>%
      # Remove the leading "@" and the trailing ":"
      str_remove_all("^@|:$")
  }
}

id_df <- as_tibble(id_df)

# Variation ---------------------------------------------------------------
# Read both the major and minor variation plots
variation_plot_files <- list.files(path = path_10, pattern = ".*variation_plot.*\\.png$", full.names = TRUE)

# Create R-code that will gather all the variation plots, then plot them as a grid with four columns and as many rows as needed.
# Make separate grids for files containing the string "major" and "minor"
if (length(variation_plot_files) > 0) {
  # Create a grid of plots for major (cand1) and minor (cand2).
  # Phase-9 naming: filenames contain "_cand1." / "_cand2." (not "major"/"minor").
  major_plots <- variation_plot_files[grepl("_cand1\\.", variation_plot_files)]
  minor_plots <- variation_plot_files[grepl("_cand2\\.", variation_plot_files)]

  # Create a grid of plots for major
  if (length(major_plots) > 0) {
    # Read images and convert to grobs
    image_list <- lapply(major_plots, function(file) {
      img <- png::readPNG(file)
      grid::rasterGrob(img, interpolate = TRUE)
    })

    # Create a list of ggplot objects containing the images
    plot_list <- lapply(image_list, function(g) {
      ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void()
    })

    # Arrange the plots in a grid with 4 columns
    grid_arrange_major <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

    # Save the arranged grid to a PNG file with a white background
    ggsave(filename = "Variation_plot_major.png", plot = grid_arrange_major, width = 12, height = 8, dpi = 300, bg = "white")
  }

  # Create a grid of plots for minor
  if (length(minor_plots) > 0) {    # Read images and convert to grobs
    image_list <- lapply(minor_plots, function(file) {
      img <- png::readPNG(file)
      grid::rasterGrob(img, interpolate = TRUE)
    })

    # Create a list of ggplot objects containing the images
    plot_list <- lapply(image_list, function(g) {
      ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void()
    })

    # Arrange the plots in a grid with 4 columns
    grid_arrange_minor <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

    # Save the arranged grid to a PNG file with a white background
    ggsave(filename = "Variation_plot_minor.png", plot = grid_arrange_minor, width = 12, height = 8, dpi = 300, bg = "white")
  }
}

# Consensus distance to reference -----------------------------------------
# Read the TSV files produced by CONSENSUS_DISTANCE for each ranked candidate.
# Each file has columns: sample, reference, similarity_pct, n_differences, alignment_length, consensus_length

distance_files <- list.files(path = path_11, pattern = "consensus_distance\\.tsv$", full.names = TRUE)

df_consensus_distance <- tibble(
  sampleName                    = character(),
  candidate_rank                = integer(),
  consensus_similarity_pct      = numeric(),
  consensus_n_differences       = integer(),
  consensus_alignment_length    = integer(),
  consensus_length              = integer()
)

if (length(distance_files) > 0) {
  for (df_file in distance_files) {
    # Phase-9 (COMPAT-02 / D-02): parse the sample name and candidate slot from the
    # filename. The CONSENSUS_DISTANCE prefix is `${meta.id}.cand{rank}` (no reference
    # field), so the slot lives at filename position 2 and carries the integer rank.
    # Expected pattern: <sampleName>.cand{rank}.consensus_distance.tsv
    fname_parts <- str_split(basename(df_file), "\\.")[[1]]
    sample_name <- fname_parts[1]
    slot        <- fname_parts[2]
    cand_rank   <- as.integer(str_remove(slot, "^cand"))

    dat <- tryCatch(
      read_tsv(df_file, col_types = cols(
        sample = col_character(),
        reference = col_character(),
        similarity_pct = col_double(),
        n_differences = col_integer(),
        alignment_length = col_integer(),
        consensus_length = col_integer()
      )),
      error = function(e) NULL
    )

    if (!is.null(dat) && nrow(dat) > 0) {
      sample_col <- dat$sample[1]
      # Extract sampleName and slot from sample_col for validation (the iVar consensus
      # header carries the `${meta.id}.cand{rank}` prefix written by IVAR_CONSENSUS).
      file_sample_name <- str_extract(sample_col, "(?<=Consensus_)[^\\.]+")
      file_slot        <- str_extract(sample_col, "cand[0-9]+")
      if (!identical(sample_name, file_sample_name) || !identical(slot, file_slot)) {
        warning(glue::glue(
          "Consensus distance file {basename(df_file)}: filename sample/slot ({sample_name}, {slot}) does not match file content ({file_sample_name}, {file_slot})"
        ))
      }
      df_consensus_distance <- bind_rows(
        df_consensus_distance,
        tibble(
          sampleName                 = sample_name,
          candidate_rank             = cand_rank,
          consensus_similarity_pct   = dat$similarity_pct[1],
          consensus_n_differences    = dat$n_differences[1],
          consensus_alignment_length = dat$alignment_length[1],
          consensus_length           = dat$consensus_length[1]
        )
      )
    }
  }
}

# Pivot to wide format: separate columns for major (rank 1) and minor (rank 2)
df_distance_wide <- df_consensus_distance %>%
  mutate(
    Major_consensus_similarity_pct   = case_when(candidate_rank == 1 ~ consensus_similarity_pct),
    Major_consensus_n_differences    = case_when(candidate_rank == 1 ~ consensus_n_differences),
    Minor_consensus_similarity_pct   = case_when(candidate_rank == 2 ~ consensus_similarity_pct),
    Minor_consensus_n_differences    = case_when(candidate_rank == 2 ~ consensus_n_differences)
  ) %>%
  select(sampleName,
         Major_consensus_similarity_pct, Major_consensus_n_differences,
         Minor_consensus_similarity_pct, Minor_consensus_n_differences) %>%
  group_by(sampleName) %>%
  fill(everything(), .direction = "downup") %>%
  slice(1) %>%
  ungroup()


# Join dataframes ---------------------------------------------------------

# Start with the original input samplesheeet and extract only the sample names. This is to ensure that all samples are included in the final summary, even if they have no data.
input_samplesheet <- read_csv(samplesheet) %>%
  select("sampleName" = sample)

final <- input_samplesheet %>%
  # Add sequencer id
  left_join(id_df, join_by(sampleName)) %>%
  # Add number of raw and trimmed reads
  left_join(trimmed_df, by = "sampleName") %>%
  # Add number of classified reads from Kraken2
  left_join(kraken_df, by = "sampleName") %>%
  # Add total mapped reads from first mapping
  left_join(parsefirstmapping_df, join_by(sampleName)) %>%
  # Add mapped reads stats
  left_join(df_mapped_reads, join_by(sampleName)) %>%
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Add consensus distance to reference
  left_join(df_distance_wide, join_by(sampleName)) %>%
  # Add de novo / BLAST evidence (PLUMB-01). Samplesheet anchors the left side so a
  # sample with no de novo output keeps its row with NA de novo fields (PLUMB-02).
  # The downstream select(..., everything()) carries the four denovo_ columns through.
  left_join(df_denovo, join_by(sampleName)) %>%
  # Phase 7 (ASUP-02): per-candidate assembly support, pivoted to wide cand_<rank>
  # slots. Samplesheet still anchors the left side so a sample with no candidate /
  # no de novo support keeps its row with NA support columns (criterion #3, no row
  # loss). everything() in the select reorder below carries these columns through.
  left_join(candidate_support_wide, join_by(sampleName)) %>%
  # Phase 8 (D-16): per-sample role slots + overall_sample_call. Dominant ->
  # Major_role_*, corroborated co-infection -> Minor_role_*, plus the new wide
  # overall_sample_call. Samplesheet anchors the left side so a no-candidate
  # sample keeps its row with NA role columns. everything() in the select reorder
  # below carries the role columns through.
  left_join(candidate_roles_wide, join_by(sampleName))

if (nrow(glue_report) > 0) {
  final <- final %>%
    # Add glue result. Only Majority currently
    left_join(glue_report, by = c("sampleName" = "Sample"))
}

# Add script name and version
final <- final %>%
  add_column("pipeline_version" = script_name_version)


# Decide if a sample is "typable" or not
final <- final %>%
  mutate(major_typable = case_when(
    Major_cov_breadth_min_1 >= 10 & Major_avg_depth >= 2 ~ "YES",
    .default = "NO"
  )) %>%
  mutate(minor_typable = case_when(
    Minor_cov_breadth_min_1 >= 10 & Minor_avg_depth >= 2 ~ "YES",
    .default = "NO"
  ))

# If minor genotype is the same as major, then not typable. But only possible if there are minor glue reports available
# But allow for the co-infection of 1a and 1b even though these belong to the same genotype
if (nrow(glue_report) > 0 & exists("gt_check")) {
  final <- final %>%
    left_join(gt_check, by = c("sampleName" = "Sample")) %>%
    mutate(minor_typable = case_when(
      identical_geno == "NO" ~ "YES",                              # Different genotypes, so minor is typable
      identical_geno == "YES" & identical_subgeno == "NO" &
        ((Major_subtype == "1a" & Minor_subtype == "1b") |
         (Major_subtype == "1b" & Minor_subtype == "1a")) ~ "YES", # If the genotype is the same and subtypes are different, but must be 1a and 1b combination. Then allow typable Minor
      identical_geno == "YES" & identical_subgeno == "NO" ~ "NO",  # If the genotype is the same and subtypes are different, but not 1a and 1b combination. Then not typable Minor
      identical_geno == "YES" & identical_subgeno == "YES" ~ "NO", # Same genotype & same subtype → not typable
      is.na(identical_geno) ~ "UNKNOWN"
    ))
}

# De novo confirmation of the reported minor — RETIRED (D-15). The legacy
# apply_denovo_layer() / minor_denovo_status / coinfection_flag chokepoint has been
# REPLACED by the Phase-8 N-candidate role classifier (score_candidates() +
# classify_roles() above, run over the long candidate_support frame). There is now
# ONE confirmation system, not two: the corroboration verdict + asymmetric refute
# (D-10/D-11) + the HCV exceptions (is_valid_minor, D-12) live in classify_roles(),
# and review_flag (below) is rewired onto role / role_reason / overall_sample_call.
# bin/denovo_layer.R / bin/denovo_confirm.R are still SOURCED above so the staged
# files load cleanly and the unit suite can exercise them directly, but their
# consumption here is removed. Phase 9 handles the Major_*/Minor_* column aliasing
# (COMPAT-03); Phase 8 only retires the legacy LOGIC.

# De novo subtype comparison columns (ODH-01). Extract the leading subtype token
# from both the mapping reference names (Major_reference / Minor_reference) and the
# de novo BLAST top-hit reference names (denovo_major_ref / denovo_minor_ref), then
# cross-compare them. All four columns are additive; existing columns are unchanged.
# str_extract returns NA for NA/NULL inputs (safe; see T-odh-01 in threat model).
final <- final %>%
  mutate(
    denovo_major_subtype = str_extract(denovo_major_ref, "^[^_]+"),
    denovo_minor_subtype = str_extract(denovo_minor_ref, "^[^_]+"),
    denovo_major_subtype_match = case_when(
      is.na(denovo_major_subtype) | is.na(Major_reference) ~ NA_character_,
      str_extract(Major_reference, "^[^_]+") == denovo_major_subtype ~ "YES",
      .default = "NO"
    ),
    denovo_minor_subtype_match = case_when(
      is.na(denovo_minor_subtype) | is.na(Minor_reference) ~ NA_character_,
      str_extract(Minor_reference, "^[^_]+") == denovo_minor_subtype ~ "YES",
      .default = "NO"
    )
  )

# If the GLUE report is missing, and GLUE columns with NAs
if (!"GLUE_genotype" %in% colnames(final)) {
  final <- final %>%
    add_column("Reference" = NA_character_,
               "GLUE_genotype" = NA_character_,
               "GLUE_subtype" = NA_character_,
               "glecaprevir" = NA_character_,
               "glecaprevir_mut" = NA_character_,
               "glecaprevir_mut_short" = NA_character_,
               "grazoprevir" = NA_character_,
               "grazoprevir_mut" = NA_character_,
               "grazoprevir_mut_short" = NA_character_,
               "paritaprevir" = NA_character_,
               "paritaprevir_mut" = NA_character_,
               "paritaprevir_mut_short" = NA_character_,
               "voxilaprevir" = NA_character_,
               "voxilaprevir_mut" = NA_character_,
               "voxilaprevir_mut_short" = NA_character_,
               "NS34A" = NA_character_,
               "NS34A_short" = NA_character_,
               "daclatasvir" = NA_character_,
               "daclatasvir_mut" = NA_character_,
               "daclatasvir_mut_short" = NA_character_,
               "elbasvir" = NA_character_,
               "elbasvir_mut" = NA_character_,
               "elbasvir_mut_short" = NA_character_,
               "ledipasvir" = NA_character_,
               "ledipasvir_mut" = NA_character_,
               "ledipasvir_mut_short" = NA_character_,
               "ombitasvir" = NA_character_,
               "ombitasvir_mut" = NA_character_,
               "ombitasvir_mut_short" = NA_character_,
               "pibrentasvir" = NA_character_,
               "pibrentasvir_mut" = NA_character_,
               "pibrentasvir_mut_short" = NA_character_,
               "velpatasvir" = NA_character_,
               "velpatasvir_mut" = NA_character_,
               "velpatasvir_mut_short" = NA_character_,
               "NS5A" = NA_character_,
               "NS5A_short" = NA_character_,
               "dasabuvir" = NA_character_,
               "dasabuvir_mut" = NA_character_,
               "dasabuvir_mut_short" = NA_character_,
               "sofosbuvir" = NA_character_,
               "sofosbuvir_mut" = NA_character_,
               "sofosbuvir_mut_short" = NA_character_,
               "NS5B" = NA_character_,
               "NS5B_short" = NA_character_,
               "HCV project version" = NA_character_,
               "GLUE engine version" = NA_character_,
               "PHE drug resistance extension version" = NA_character_,
               )
}

# Ensure de novo subtype comparison columns are always present in the schema,
# even when GLUE is absent and the above add_column() block runs but does not
# include them. Since the columns are derived unconditionally above, they already
# exist at this point; this guard is a no-op in normal execution and exists only
# as a safety net for any future refactor that moves the derivation block.
if (!"denovo_major_subtype" %in% colnames(final)) {
  final <- final %>%
    add_column(
      "denovo_major_subtype"       = NA_character_,
      "denovo_minor_subtype"       = NA_character_,
      "denovo_major_subtype_match" = NA_character_,
      "denovo_minor_subtype_match" = NA_character_
    )
}

# Per-sample role-reason summary for the rewired review_flag (D-13/D-15). The
# review sentences are now derived from the N-candidate role model, not the retired
# minor_denovo_status / coinfection_flag. Roll the long classified candidate_support
# frame up to one row per sample, capturing whether ANY candidate was refuted by de
# novo (background/refuted_denovo) or kept as an uncorroborated co-infection
# (co-infection/uncorroborated_kept). These per-sample booleans feed the pmap_chr
# below alongside overall_sample_call + gate_flag + the subtype-match columns.
if (nrow(candidate_support) > 0) {
  role_review <- candidate_support %>%
    group_by(sampleName) %>%
    summarise(
      any_refuted_denovo     = any(role_reason == "refuted_denovo", na.rm = TRUE),
      any_uncorroborated     = any(role_reason == "uncorroborated_kept", na.rm = TRUE),
      .groups = "drop"
    )
} else {
  role_review <- tibble(
    sampleName         = character(),
    any_refuted_denovo = logical(),
    any_uncorroborated = logical()
  )
}

final <- final %>%
  left_join(role_review, join_by(sampleName))

# Review flag (REVIEW-01), rewired onto the Phase-8 roles (D-13/D-15). Human-readable
# inspection prompts for samples that warrant manual review, joined with " | ". NA
# when no reasons fire. The verbatim message text lives in the pmap_chr() below; the
# triggers, in order, are now driven by role / role_reason / overall_sample_call
# instead of the retired minor_denovo_status / coinfection_flag:
#   1. overall_sample_call == "co-infection" AND a major or minor subtype mismatch
#      (denovo_*_subtype_match == "NO") — co-infection called but major/minor
#      assignment uncertain (de novo and mapping disagree on which strain is dominant)
#   2. monoinfection AND denovo_major_subtype_match == "NO" — major subtype conflict
#      between de novo assembly and mapping
#   3. any_refuted_denovo — a minor candidate refuted by de novo; likely single infection
#   4. any_uncorroborated — a co-infection kept without de novo corroboration (de novo
#      inconclusive for both strains); warrants analyst review of QC plots / contigs
#   5. overall_sample_call == "indeterminate" — no candidate passed the major-gate
#   6. gate_flag != "ok" — major failed the first-mapping quality thresholds
# (Earlier versions emitted semicolon-separated reason codes; rewritten to full
# sentences in commit ff12009; rewired onto roles in Phase 8 / D-15.)
#
# MultiQC orange-highlight note: in assets/multiqc_config.yml the results_summary
# custom_data block includes a cond_formatting_rules entry for this column that
# colours any non-NA value orange (warn class). See the pconfig.cond_formatting_rules
# key added there. If that config is absent (older deployments), MultiQC falls back
# to plain text — the column is still useful as a text summary.
final <- final %>%
  mutate(review_flag = {
    pmap_chr(
      list(
        denovo_major_subtype_match,
        denovo_minor_subtype_match,
        overall_sample_call,
        any_refuted_denovo,
        any_uncorroborated,
        gate_flag
      ),
      function(maj_match, min_match, sample_call, refuted, uncorr, gflag) {
        msgs        <- character(0)
        is_coinf    <- !is.na(sample_call) && sample_call == "co-infection"
        is_mono     <- !is.na(sample_call) && sample_call == "monoinfection"
        is_indet    <- !is.na(sample_call) && sample_call == "indeterminate"
        subtype_dis <- (!is.na(maj_match) && maj_match == "NO") || (!is.na(min_match) && min_match == "NO")
        if (is_coinf && subtype_dis)
          msgs <- c(msgs, "Co-infection confirmed, but major/minor assignment uncertain — de novo and mapping disagree on which strain is dominant. Please review.")
        if (is_mono && !is.na(maj_match) && maj_match == "NO")
          msgs <- c(msgs, "Major subtype conflict between de novo assembly and mapping — possible reference mismatch or highly divergent strain. Please review.")
        if (isTRUE(refuted))
          msgs <- c(msgs, "Minor strain candidate refuted by de novo assembly — likely single infection.")
        if (isTRUE(uncorr))
          msgs <- c(msgs, "Co-infection kept without de novo corroboration — de novo inconclusive for both strains; minor strain may be a genuine low-yield co-infection. Please review.")
        if (is_indet)
          msgs <- c(msgs, "No candidate passed the major-gate — overall sample call indeterminate.")
        if (!is.na(gflag) && gflag != "ok")
          msgs <- c(msgs, "Major strain failed mapping quality thresholds — genotype call uncertain.")
        if (length(msgs) == 0) NA_character_ else paste(msgs, collapse = " | ")
      }
    )
  }) %>%
  # Drop the per-sample role-review helper booleans now they have been consumed.
  select(-any_refuted_denovo, -any_uncorroborated)

# Reorder columns
final <- final %>%
  select(sampleName,
         total_raw_reads,
         total_trimmed_reads,
         total_classified_reads,
         total_mapped_reads,
         fraction_mapped_reads_vs_median,
         Major_genotype_mapping,
         Major_reference,
         Minor_genotype_mapping,
         Minor_reference,
         major_typable,
         minor_typable,
         # Phase 8 (D-16): overall sample call from the N-candidate role model,
         # placed where the retired minor_denovo_status / coinfection_flag sat.
         overall_sample_call,
         Major_role_reference,
         Major_role_subtype,
         Major_dominance_score,
         Major_role_reason,
         Minor_role_reference,
         Minor_role_subtype,
         Minor_dominance_score,
         Minor_role_reason,
         denovo_major_subtype,
         denovo_minor_subtype,
         denovo_major_subtype_match,
         denovo_minor_subtype_match,
         review_flag,
         Reads_withdup_mapped_major,
         Reads_nodup_mapped_major,
         Percent_reads_mapped_of_trimmed_with_dups_major,
         Major_cov_breadth_min_5,
         Major_cov_breadth_min_10,
         percent_mapped_reads_major_firstmapping,
         any_of(c("Major_consensus_similarity_pct", "Major_consensus_n_differences")),
         Reads_withdup_mapped_minor,
         Reads_nodup_mapped_minor,
         Percent_reads_mapped_of_trimmed_with_dups_minor,
         Minor_cov_breadth_min_5,
         Minor_cov_breadth_min_10,
         percent_mapped_reads_minor_firstmapping,
         any_of(c("Minor_consensus_similarity_pct", "Minor_consensus_n_differences")),
         everything()) %>%
  distinct() %>% # Remove any duplicated rows from the different joins
  # If there are no minor genotype reports, then the identical_geno and identical_subgeno columns will not exist. Therefore use any_of in case they are not there
  select(-any_of(c("Major_minor", "identical_geno", "identical_subgeno")))

# Write file
write_csv(final, file = "Summary.csv")

# Write file for MultiQC
# Add MultiQC info lines
header <- c("# id: 'summary'",
            "# section_name: 'Summary'",
            "# description: 'These statistics are generated from the process SUMMARIZE and the R script summarize.R",
            "# format: 'csv'")

# Convert final data to data frame
tt <- as.data.frame(final)

# TSV avoids quoting issues when field values contain commas (e.g. review_flag sentences).
# MultiQC config must match: file_format: tsv, fn: "*/summary_mqc.tsv"
file <- "summary_mqc.tsv"

# Add the column names to file
tt %>% colnames() %>% paste0(collapse = "\t") %>% write_lines(file, append = TRUE)

# Write the data to file
write_tsv(tt, file, append = TRUE) # colnames will not be included

