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
denovo_min_contig_length  <- if (length(args) >= 4 && nchar(args[4]) > 0) as.numeric(args[4]) else 500
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
# Contig-length floor for the monoinfection different-genotype-contig REVIEW sentence
# (260803-ogc). APPENDED at args[16] — strictly after score_evenness_k — for the same
# positional reason as the block above (T-08-04): inserting mid-string silently
# re-maps every later index. Default mirrors nextflow.config.
review_min_offgenotype_contig_length <-
  if (length(args) >= 16 && nchar(args[16]) > 0) as.numeric(args[16]) else 1000

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
    confirmation_status = col_character(),
    # Phase-10 rescue audit columns (RESCUE_EVALUATION emits them on every candidates
    # CSV; declare as character so an all-NA column is not inferred to <logical>).
    rescued_from        = col_character(),
    rescue_trigger      = col_character()
  ))) %>%
    rename(sampleName = sample)
  # A pre-Phase-10 / skip-assembly candidates CSV may lack the rescue columns entirely;
  # add them NA-filled so the wide pivot + rescue_flag rollup below always find them.
  if (!"rescued_from" %in% names(candidates_long)) {
    candidates_long <- candidates_long %>% mutate(rescued_from = NA_character_)
  }
  if (!"rescue_trigger" %in% names(candidates_long)) {
    candidates_long <- candidates_long %>% mutate(rescue_trigger = NA_character_)
  }
} else {
  # Declare all eight Phase-6 candidate columns + the two Phase-10 rescue audit
  # columns with their types so a no-candidate run yields a typed zero-row frame
  # (the join then returns a typed zero-row frame; rescue columns never vanish).
  candidates_long <- tibble(
    sampleName          = character(),
    candidate_rank      = integer(),
    candidate_ref       = character(),
    candidate_subtype   = character(),
    candidate_genotype  = character(),
    candidate_reads     = double(),
    candidate_cov       = double(),
    confirmation_status = character(),
    rescued_from        = character(),
    rescue_trigger      = character()
  )
}

# Phase-9 (COMPAT-02 / D-02): per-sample (candidate_ref -> candidate_rank) lookup
# used by the idxstats read-count loops + the coverage loop to recover the slot
# by join. Built once from the hoisted candidates_long frame.
#
# Phase-11 (JMAP-03 / Pitfall 4): the read-count loops now parse SAMTOOLS_IDXSTATS
# on the COMBINED BAM — ONE file per sample, keyed only on the bare reference name
# (no `.cand{rank}.` filename slot to carry the rank). Rank is therefore recovered
# ENTIRELY by joining the idxstats refname to this lookup. The lookup must keep
# EVERY (sampleName, candidate_rank) row — NOT collapse duplicate refs — so that
# when two candidates share one reference (the commit-1d1a051 scenario) BOTH ranks
# survive the join. A previous `distinct(sampleName, candidate_ref, .keep_all)`
# collapsed a shared ref to a single rank, silently dropping the minor candidate.
# distinct() over ALL three key columns only de-duplicates true duplicate rows.
candidate_rank_lookup <- candidates_long %>%
  select(sampleName, candidate_ref, candidate_rank) %>%
  filter(!is.na(candidate_rank)) %>%
  distinct(sampleName, candidate_ref, candidate_rank)

# Second mapping, reads mapped with duplicates ----------------------------
# Phase-11 (JMAP-03 / D-15): read counts now come from SAMTOOLS_IDXSTATS on the
# COMBINED withdup BAM — ONE file per sample (<sample>.withdup.idxstats), a
# 4-column no-header TSV (refname, seqlen, mapped, unmapped) with one DATA row per
# candidate reference plus a trailing `*` unmapped row. This replaces the former
# per-candidate SAMTOOLS_STATS text files parsed via `X2 == "reads mapped:"`.
# The withdup loop therefore changes from one-file-per-candidate to
# one-file-per-sample, emitting one row per candidate.
#
# The glob pins `\.withdup\.idxstats$` so the PARSEFIRSTMAPPING
# `*.firstmapping.withdup.idxstats` files (if ever co-staged) do NOT match — only
# JOINT_MAPPING's per-sample `<sample>.withdup.idxstats` is consumed here.
# List files
stats_files <- list.files(path = path_4, pattern = "\\.withdup\\.idxstats$", full.names = TRUE)
# Exclude any PARSEFIRSTMAPPING `*.firstmapping.withdup.idxstats` that may share the
# dir — only JOINT_MAPPING's per-sample `<sample>.withdup.idxstats` is consumed here.
stats_files <- stats_files[!str_detect(basename(stats_files), "\\.firstmapping\\.")]

if (length(stats_files) > 0) {
  tmp_df <- purrr::map_dfr(stats_files, function(f) {
    # Sample name = first dot-delimited basename token.
    sampleName <- str_split(basename(f), "\\.")[[1]][1]
    # 4-col no-header idxstats TSV: refname, seqlen, mapped, unmapped.
    read_tsv(f, col_names = c("candidate_ref", "seqlen", "mapped", "unmapped"),
             comment = "", show_col_types = FALSE) %>%
      # Drop only the `*` unmapped trailer row. A real candidate reference whose
      # mapped count is 0 under competitive joint mapping (Phase 11 assigned all
      # reads to the dominant) is KEPT — the Plan 01/02 evidence engine scores it
      # on its own assembly evidence rather than relying on a silent upstream drop
      # (EVID-02 / folded todo 2026-06-22). The 0 read count flows through as a
      # numeric 0 (not NA), so the downstream case_when/percentage/fill logic
      # stays NA-tolerant and the candidate_rank_lookup join stays one-row-per-rank.
      filter(candidate_ref != "*") %>%
      # idxstats refnames are already the bare candidate reference; strip a
      # `_cand{rank}` slot defensively (no-op for combined-BAM idxstats).
      mutate(candidate_ref = str_remove(candidate_ref, "_cand[0-9]+$")) %>%
      transmute(sampleName, candidate_ref,
                trimmed_reads_withdups_mapped = as.numeric(mapped))
  })
} else {
  tmp_df <- tibble(sampleName = character(), candidate_ref = character(),
                   trimmed_reads_withdups_mapped = numeric())
}
# Phase-11 (JMAP-03 / Pitfall 4): recover candidate_rank by joining the idxstats
# refname to candidate_rank_lookup. idxstats is per-sample so the filename carries
# NO per-candidate slot — the lookup is the ONLY rank source. When two candidates
# share a reference, the lookup holds both ranks for that ref, so the join yields
# one row per rank (each with the same per-reference mapped count); the per-sample
# group/fill/slice below collapses them to a single Summary.csv row without a
# many-to-many explosion (the full_join keys on Major_/Minor_reference, identical
# across df_with_dups/df_nodups).
tmp_df <- tmp_df %>%
  left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref")) %>%
  mutate(candidate_rank = suppressWarnings(as.integer(candidate_rank)))

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
  select(-candidate_ref) %>%
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
# Phase-11 (JMAP-03 / D-15): symmetric idxstats migration of the nodup loop.
# Read counts come from SAMTOOLS_IDXSTATS on the COMBINED dedup BAM — ONE file per
# sample (<sample>.nodup.idxstats), same 4-column no-header TSV (refname, seqlen,
# mapped, unmapped) with one row per candidate + a `*` trailer. Replaces the
# per-candidate SAMTOOLS_STATS `X2 == "reads mapped:"` parse. The glob pins
# `\.nodup\.idxstats$` and excludes any `*.firstmapping.nodup.idxstats` so only
# JOINT_MAPPING's per-sample nodup idxstats is consumed.
# List files
stats_files <- list.files(path = path_5, pattern = "\\.nodup\\.idxstats$", full.names = TRUE)
# Exclude any PARSEFIRSTMAPPING `*.firstmapping.nodup.idxstats` that may share the dir.
stats_files <- stats_files[!str_detect(basename(stats_files), "\\.firstmapping\\.")]

if (length(stats_files) > 0) {
  tmp_df <- purrr::map_dfr(stats_files, function(f) {
    sampleName <- str_split(basename(f), "\\.")[[1]][1]
    read_tsv(f, col_names = c("candidate_ref", "seqlen", "mapped", "unmapped"),
             comment = "", show_col_types = FALSE) %>%
      # Symmetric with the withdup loop: drop only the `*` trailer, keep a real
      # candidate whose mapped count is 0 so a zero-read co-infection minor is not
      # silently dropped (EVID-02 / folded todo 2026-06-22).
      filter(candidate_ref != "*") %>%
      mutate(candidate_ref = str_remove(candidate_ref, "_cand[0-9]+$")) %>%
      transmute(sampleName, candidate_ref,
                trimmed_reads_nodups_mapped = as.numeric(mapped))
  })
} else {
  tmp_df <- tibble(sampleName = character(), candidate_ref = character(),
                   trimmed_reads_nodups_mapped = numeric())
}
# Phase-11 (JMAP-03 / Pitfall 4): rank recovered via candidate_rank_lookup (the
# only rank source — idxstats is per-sample, no filename slot). Both ranks of a
# shared reference survive (the lookup keeps every (ref, rank) pair), so
# targeted_nodup_per_cand below stays keyed on candidate_rank with no duplicate
# (sampleName, candidate_rank) rows and no many-to-many join downstream.
tmp_df <- tmp_df %>%
  left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref")) %>%
  mutate(candidate_rank = suppressWarnings(as.integer(candidate_rank)))

# Phase-10/11: targeted_nodup_per_cand joins by candidate_rank (not candidate_ref)
# so rescued samples where two slots share a reference do not create a many-to-many
# join. With the idxstats migration the rank comes from candidate_rank_lookup.
targeted_nodup_per_cand <- tmp_df %>%
  filter(!is.na(candidate_rank)) %>%
  distinct(sampleName, candidate_rank, .keep_all = TRUE) %>%
  select(sampleName, candidate_rank, targeted_reads_nodup = trimmed_reads_nodups_mapped)

df_nodups <- tmp_df %>%
  # Create columns for major and minor
  separate(candidate_ref, into = c("genotype", NA), sep = "_", remove = F) %>%
  mutate(Major_genotype_mapping = case_when(candidate_rank == 1 ~ genotype)) %>%
  mutate(Minor_genotype_mapping = case_when(candidate_rank == 2 ~ genotype)) %>%
  select(-genotype) %>%
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_nodup_mapped_major = case_when(candidate_rank == 1 ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_minor = case_when(candidate_rank == 2 ~ trimmed_reads_nodups_mapped)) %>%
  # Phase-11 (JMAP-03): the combined-BAM nodup idxstats carries ONLY candidate
  # references (no `first_mapping` row), so this vestigial column is always NA —
  # kept NA-filled to preserve the Summary.csv header contract. The legacy
  # `reference == "first_mapping"` source no longer exists post-migration.
  mutate(Reads_nodup_mapped_first_mapping = NA_real_) %>%
  select(-trimmed_reads_nodups_mapped) %>%
  #mutate(Percent_mapped_major = case_when(candidate_rank == 1 ~ Percent_trimmed_reads_mapped)) %>%
  #mutate(Percent_mapped_minor = case_when(candidate_rank == 2 ~ Percent_trimmed_reads_mapped)) %>%
  # Create columns for the major and minor references (cleaned candidate_ref — no
  # slot suffix — byte-identical to the legacy join key, COMPAT-02 D-02).
  mutate(Major_reference = case_when(candidate_rank == 1 ~ candidate_ref)) %>%
  mutate(Minor_reference = case_when(candidate_rank == 2 ~ candidate_ref)) %>%
  # Create one row per sample
  select(-candidate_ref, -candidate_rank) %>%
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
  # Pin col_types so per-file readr inference can't disagree across samples: a
  # sample whose major/minor contig length is empty (NA) infers logical/character,
  # while a populated one infers double -> map_dfr/bind_rows aborts ("Can't combine
  # <double> and <character>"). Same header-only/empty-cell trap guarded below for
  # candidates/assembly_support. Schema mirrors blast_parse.R summary CSV.
  df_denovo <- map_dfr(blastparse_files, ~ read_csv(.x, col_types = cols(
    sample              = col_character(),
    major_ref           = col_character(),
    major_contig_length = col_double(),
    minor_ref           = col_character(),
    minor_contig_length = col_double()
  ))) %>%
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
  # Pin col_types to the full blast_parse.R _blast_out.csv schema: a sample whose
  # contigs were all dropped yields a header-only file, which readr would type as
  # all-logical and abort map_dfr/bind_rows against the populated files (same trap
  # that crashed the blastparse read above). Coercing here keeps a header-only file
  # contributing zero rows (NA-fill, "unconfirmed") instead of failing the run.
  df_blast_out <- map_dfr(blast_out_files, ~ read_csv(.x, col_types = cols(
    qseqid    = col_character(),
    sseqid    = col_character(),
    subtype   = col_character(),
    pident    = col_double(),
    length    = col_double(),
    mismatch  = col_double(),
    gapopen   = col_double(),
    qstart    = col_double(),
    qend      = col_double(),
    sstart    = col_double(),
    send      = col_double(),
    evalue    = col_double(),
    bitscore  = col_double(),
    sc_length = col_double(),
    kmer_cov  = col_double()
  )) %>%
    mutate(sampleName = str_remove(basename(.x), "_blast_out.csv$")))
} else {
  df_blast_out <- tibble(
    sampleName = character(),
    qseqid     = character(),
    sseqid     = character(),
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
# GLUE per-candidate (D8): read the major/minor collected reports early so the
# per-candidate GLUE genotype is available for apply_concordance() after the
# assembly-support join. These objects are also reused by the same-genotype check
# and the final Summary.csv join further below (the blocks there remain in place;
# they read glue_report / glue_report_minor which are now defined here).
glue_file <- list.files(path = path_8, pattern = "GLUE_collected_report_major.tsv$", full.names = TRUE)
# Read every column as character (.default): the GLUE report schema is entirely
# textual (genotype/subtype, drug-resistance status, mutation strings, *_short
# codes, versions). A resistance column that is all-NA in one report infers as
# logical/double while the other report has mutation strings -> bind_rows of the
# major/minor reports (glue_by_rank, glue_per_cand) aborts ("Can't combine
# <character> and <double>"). Forcing character pins the whole schema in one shot.
glue_report <- if (length(glue_file) > 0) read_tsv(glue_file, col_types = cols(.default = col_character())) else tibble()
glue_file_minor <- list.files(path = path_8, pattern = "GLUE_collected_report_minor.tsv$", full.names = TRUE)
glue_report_minor <- if (length(glue_file_minor) > 0) read_tsv(glue_file_minor, col_types = cols(.default = col_character())) else tibble()

# Build the per-candidate GLUE frame (rank 1 = major report, rank 2 = minor report)
# and join candidate_glue_genotype + candidate_glue_subtype onto candidates_long so
# assembly_support_join() passes them through into candidate_support.
glue_cand1 <- if (nrow(glue_report) > 0) {
  glue_report %>% transmute(sampleName = Sample, candidate_rank = 1L,
                             candidate_glue_genotype = GLUE_genotype,
                             candidate_glue_subtype  = GLUE_subtype)
} else {
  tibble(sampleName = character(), candidate_rank = integer(),
         candidate_glue_genotype = character(), candidate_glue_subtype = character())
}

glue_cand2 <- if (nrow(glue_report_minor) > 0) {
  glue_report_minor %>% transmute(sampleName = Sample, candidate_rank = 2L,
                                   candidate_glue_genotype = GLUE_genotype,
                                   candidate_glue_subtype  = GLUE_subtype)
} else {
  tibble(sampleName = character(), candidate_rank = integer(),
         candidate_glue_genotype = character(), candidate_glue_subtype = character())
}

glue_per_cand <- bind_rows(glue_cand1, glue_cand2)
candidates_long <- candidates_long %>%
  left_join(glue_per_cand, by = c("sampleName", "candidate_rank"))

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
  left_join(cv_by_ref,               by = c("sampleName", "candidate_ref")) %>%
  left_join(targeted_nodup_per_cand, by = c("sampleName", "candidate_rank"))

# D8 concordance pre-annotation: annotate each candidate with concordance_status
# (confirmed/unconfirmed/discordant) + concordance_reason BEFORE scoring and role
# classification. classify_roles() will consume concordance_status in Task 4.
candidate_support <- apply_concordance(candidate_support)

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

# Per-sample candidate_rank of the role-dominant. Used below to swap Major_*/Minor_*
# stat/coverage/consensus/GLUE columns when the role classifier's dominant is at
# candidate_rank 2 (i.e. the mapping-minor carried the true dominant strain, as in
# a co-infection where first-mapping read mis-recruitment inverted the abundance order).
role_dominant_rank <- if (nrow(candidate_support) > 0) {
  candidate_support %>%
    filter(role == "dominant") %>%
    group_by(sampleName) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(sampleName, dominant_cand_rank = as.integer(candidate_rank))
} else {
  tibble(sampleName = character(), dominant_cand_rank = integer())
}

# EVID-06/D-11: per-sample REFERENCE name of the role-dominant candidate. role_dominant_rank
# (above) carries only the rank; sample_review_message() also needs the ref to NAME the
# dominant candidate in its D-11 triggers (dominant_unconfirmed / major_ref_changed) and the
# enriched monoinfection subtype-conflict messages (RESEARCH Pattern 3 / Pitfall 3: enrich the
# sample-level trigger with the candidate rank+ref from the L850 lookup). Transient — joined
# into final just before the review_flag mutate and dropped in the same select() afterwards.
role_dominant_ref <- if (nrow(candidate_support) > 0) {
  candidate_support %>%
    filter(role == "dominant") %>%
    group_by(sampleName) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(sampleName, dominant_cand_ref = as.character(candidate_ref))
} else {
  tibble(sampleName = character(), dominant_cand_ref = character())
}

# EVID-05: one contig-language evidence_summary sentence per candidate, built by the
# pure Plan-01 helper build_evidence_summary() (defined in classify_roles.R). Threads
# the candidate's role/state + its OWN best-contig metrics so a reader can reconstruct
# the candidate's evidence_state from candidates.csv alone (success criteria 1+2). The
# three contig_*_contribution columns already ride along from score_assembly_support()
# (Plan 01), so the un-select()ed raw write_csv below surfaces all four automatically.
# pmap_chr over a zero-row frame returns character(0), so this is empty-batch safe.
candidate_support <- candidate_support %>%
  mutate(evidence_summary = pmap_chr(
    list(role, role_reason, evidence_state, concordance_status,
         candidate_rank, candidate_ref, candidate_subtype, assembly_exists,
         assembly_support_best_contig_length,
         assembly_support_best_contig_pident,
         assembly_support_best_contig_kmer_cov),
    build_evidence_summary
  ))

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
      Major_role_reason        = role_reason,
      # WR-05 (12-REVIEW): surface the calibrated confirmed-vs-probable evidence
      # band on the wide, human-facing Summary.csv — before this, evidence_state
      # survived only in the long candidates.csv, so a marginal 0.50-score
      # co-infection call was indistinguishable from a strong 0.95 one at a glance.
      Major_evidence_state     = evidence_state,
      # EVID-05/D-04: carry this slot's OWN best-contig metrics into the wide frame
      # so build_evidence() can append a contig-corroboration token to Major_evidence.
      # Transient — consumed by build_evidence() then dropped by the reorder select.
      Major_best_contig_length   = assembly_support_best_contig_length,
      Major_best_contig_pident   = assembly_support_best_contig_pident,
      Major_best_contig_kmer_cov = assembly_support_best_contig_kmer_cov
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
      Minor_role_reason        = role_reason,
      # WR-05: see the Major_evidence_state comment above.
      Minor_evidence_state     = evidence_state,
      # EVID-05/D-04: this slot's OWN best-contig metrics for the Minor_evidence token.
      Minor_best_contig_length   = assembly_support_best_contig_length,
      Minor_best_contig_pident   = assembly_support_best_contig_pident,
      Minor_best_contig_kmer_cov = assembly_support_best_contig_kmer_cov
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
    Major_evidence_state     = character(),
    Major_best_contig_length   = double(),
    Major_best_contig_pident   = double(),
    Major_best_contig_kmer_cov = double(),
    Minor_role_reference     = character(),
    Minor_role_subtype       = character(),
    Minor_dominance_score    = double(),
    Minor_role_reason        = character(),
    Minor_evidence_state     = character(),
    Minor_best_contig_length   = double(),
    Minor_best_contig_pident   = double(),
    Minor_best_contig_kmer_cov = double()
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
  "assembly_support_best_contig_kmer_cov",
  # Phase-10 rescue audit columns: pivot to cand_{rank}_rescued_from /
  # cand_{rank}_rescue_trigger (D-12 intent, in-file cand_{rank}_ underscore spelling).
  "rescued_from",
  "rescue_trigger"
)
# Column types per support value, in support_value_cols order: the two status/
# subtype columns are character, the four metrics are double, and the two rescue
# audit columns are character.
support_value_is_character <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)

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
# glue_report / glue_report_minor are read earlier (before the assembly-support
# section) so apply_concordance() has the GLUE leg available. The same-genotype
# check and the final Summary.csv join below consume these already-defined objects.

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

# Attach the per-sample dominant_cand_rank so the swap and GLUE join below can
# use it. Samples with no classified candidate (no-mapping) get NA → no swap.
final <- final %>%
  left_join(role_dominant_rank, by = "sampleName")

# Re-key Major_*/Minor_* stat/coverage/consensus columns to the role-dominant strain.
# All stats loops above key by candidate_rank (1=mapping-major, 2=mapping-minor).
# When classify_roles() assigns dominance to candidate_rank 2 (e.g. sim2 70:30 2a:3a
# where first-mapping mis-recruits reads to 3a), every "Major_*" column currently
# describes the wrong strain. Swap the pairs for those samples so the Summary.csv
# row is internally consistent: Major_* ↔ the role-dominant, Minor_* ↔ the co-infection.
swap_col_pairs <- list(
  c("Major_reference",                                 "Minor_reference"),
  c("Major_genotype_mapping",                          "Minor_genotype_mapping"),
  c("Reads_withdup_mapped_major",                      "Reads_withdup_mapped_minor"),
  c("Reads_nodup_mapped_major",                        "Reads_nodup_mapped_minor"),
  c("Percent_reads_mapped_of_trimmed_with_dups_major", "Percent_reads_mapped_of_trimmed_with_dups_minor"),
  c("percent_mapped_reads_major_firstmapping",         "percent_mapped_reads_minor_firstmapping"),
  c("Major_cov_breadth_min_1",                         "Minor_cov_breadth_min_1"),
  c("Major_cov_breadth_min_5",                         "Minor_cov_breadth_min_5"),
  c("Major_cov_breadth_min_10",                        "Minor_cov_breadth_min_10"),
  c("Major_avg_depth",                                 "Minor_avg_depth"),
  c("Major_consensus_similarity_pct",                  "Minor_consensus_similarity_pct"),
  c("Major_consensus_n_differences",                   "Minor_consensus_n_differences")
)

needs_swap <- !is.na(final$dominant_cand_rank) & final$dominant_cand_rank == 2L

if (any(needs_swap, na.rm = TRUE)) {
  for (pair in swap_col_pairs) {
    maj_col <- pair[1]; min_col <- pair[2]
    if (!maj_col %in% names(final) || !min_col %in% names(final)) next
    tmp <- final[[maj_col]]
    final[[maj_col]][needs_swap] <- final[[min_col]][needs_swap]
    final[[min_col]][needs_swap] <- tmp[needs_swap]
  }
}

# Join GLUE results keyed to the ROLE-DOMINANT candidate, not always candidate_rank 1
# (the mapping-major). When classify_roles() assigns dominance to candidate_rank 2,
# glue_report_minor describes the dominant strain's resistance profile.
# Build a rank-keyed combined GLUE frame; join via dominant_cand_rank so
# GLUE_subtype and drug-resistance columns follow the role assignment.
# Samples with no role result fall back to rank-1 GLUE (dominant_cand_rank NA → 1L).
#
# Design decision — resistance reporting for co-infections (D-GLUE-COINFECTION):
# HCVGLUE is run on ALL candidate BAMs (cand1 + cand2), so a resistance profile
# exists for both strains when a co-infection is confirmed. However, the primary
# output columns (GLUE_genotype, GLUE_subtype, drug-resistance columns) carry
# only the role-DOMINANT strain's profile. The role-MINOR strain's resistance
# data lives in the opposing GLUE report file but is not surfaced as separate
# columns in Summary.csv.  Rationale: clinical guidance is anchored to the
# dominant strain; adding a second resistance column set would double the column
# count and complicate downstream parsing for the common (monoinfection) case.
# If per-strain resistance for co-infections is needed in future, expose the
# non-dominant GLUE report as supplementary output rather than widening Summary.csv.
if (nrow(glue_report) > 0 || nrow(glue_report_minor) > 0) {
  glue_by_rank <- bind_rows(
    if (nrow(glue_report) > 0)
      glue_report %>% rename(sampleName = Sample) %>% mutate(.dom_rank = 1L)
    else
      tibble(sampleName = character(), .dom_rank = integer()),
    if (nrow(glue_report_minor) > 0)
      glue_report_minor %>% rename(sampleName = Sample) %>% mutate(.dom_rank = 2L)
    else
      tibble(sampleName = character(), .dom_rank = integer())
  )
  final <- final %>%
    mutate(.dom_rank = coalesce(dominant_cand_rank, 1L)) %>%
    left_join(glue_by_rank, by = c("sampleName", ".dom_rank")) %>%
    select(-.dom_rank)
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

# Consistency override: monoinfection samples always get minor_typable = NO.
# The gt_check block above produces UNKNOWN when identical_geno is NA (the minor
# GLUE subtype is absent for that sample). Whether gt_check ran or not depends on
# whether glue_report_minor was present — so refuted spurious minors ended up as
# UNKNOWN or NO depending on unrelated pipeline state. Both cases represent the
# same biology: no real minor strain. "UNKNOWN" implies ambiguity that doesn't
# exist once the role classifier has called monoinfection.
final <- final %>%
  mutate(minor_typable = if_else(
    !is.na(overall_sample_call) & overall_sample_call == "monoinfection",
    "NO",
    minor_typable
  ))

# gt_check derives Major_subtype / Minor_subtype from glue_report / glue_report_minor,
# which map to cand1 / cand2 by file-system slot name. When the dominant candidate is
# not in cand1 (rescue or neutral-ranking flip), Major_subtype carries the wrong strain.
# Override here with the role-based subtypes, which are correctly assigned by the
# Phase-8 role classifier regardless of slot order.
# Guard: when GLUE reports are absent (no-GLUE batches) gt_check never runs and
# Major_subtype / Minor_subtype are never added to final. NA-fill them here so
# the coalesce below never errors (BM2-01; guard moved from post-review_flag).
if (!"Major_subtype" %in% colnames(final))
  final <- final %>% add_column("Major_subtype" = NA_character_)
if (!"Minor_subtype" %in% colnames(final))
  final <- final %>% add_column("Minor_subtype" = NA_character_)
final <- final %>%
  mutate(
    Major_subtype = coalesce(Major_role_subtype, Major_subtype),
    Minor_subtype = coalesce(Minor_role_subtype, Minor_subtype)
  )

# Major_genotype / Minor_genotype — re-derived from the role-corrected subtypes
# immediately above (260803-ogc; comparison report §6). Two defects are fixed here:
#
#   (1) MISSING COLUMNS. They were only ever created by the gt_check block, which is
#       itself guarded on `nrow(glue_report) > 0`. A batch where HCV-GLUE could not be
#       run never got the columns at all — observed as 134 vs 136 columns across five
#       runs of the SAME pipeline version, so a downstream consumer reading
#       Major_genotype breaks on some runs and silently reports nothing on others.
#       mutate() creates them unconditionally, so the emitted schema is now stable.
#
#   (2) STALE ORDERING. They came from glue_report / glue_report_minor, which map to
#       cand1 / cand2 by file-system SLOT. Major_subtype / Minor_subtype get corrected
#       to the role-based assignment just above; the genotype columns never did. On the
#       one cohort sample whose roles are reversed relative to first-mapping (Sample51K)
#       they disagreed with the subtype columns: Major_genotype=2 alongside
#       Major_subtype=3a. Deriving them from the corrected subtypes keeps the two in
#       lockstep by construction.
#
# genotype_from_subtype() (not substr) so the 2k1b recombinant keeps its full name as
# its genotype, matching every other genotype derivation in the pipeline.
#
# ORDERING IS LOAD-BEARING: this must stay AFTER the coalesce above and AFTER the
# gt_check join that feeds minor_typable. gt_check's identical_geno / identical_subgeno
# are computed from the RAW GLUE slot values inside gt_check itself and are consumed by
# the minor_typable case_when earlier — re-deriving there would change typability,
# which IS a call. Here it changes only two reported columns.
final <- final %>%
  mutate(
    Major_genotype = genotype_from_subtype(Major_subtype),
    Minor_genotype = genotype_from_subtype(Minor_subtype)
  )

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

# De novo contig NAME per strain (RPT-CONTIG-01/02). We already carry the best
# BLAST-hit reference for each strain (denovo_major_ref / denovo_minor_ref, ingested
# from *.blastparse.csv above). Here we surface the contig (qseqid) that produced
# that hit, by resolving the per-contig long BLAST table (df_blast_out) to the
# major/minor grain — a pure summarize.R join, no change to blast_parse.R.
#
# Pairing rationale: for a given (sample, reference), the contig backing that call is
# the top-bitscore hit whose sseqid == that reference. So column (1) denovo_*_contig
# and column (2) denovo_*_ref come from ONE BLAST-hit row — a coherent same-row pair
# (the named contig's best hit landed on that reference).
#
# slice_max(..., with_ties = FALSE) collapses to one contig per (sampleName, sseqid)
# BEFORE the left_join, so a reference shared by several contigs cannot explode rows
# (T-dpj-02 many-to-many guard; mirrors the JMAP-03 1d1a051 fix). An empty
# df_blast_out (skip-assembly / no-denovo batch — the sseqid empty-tibble fix above
# is what lets this build) yields a zero-row contig_by_ref, so both left_joins still
# ADD the denovo_*_contig columns NA-filled (always present downstream).
contig_by_ref <- df_blast_out %>%
  filter(!is.na(sseqid)) %>%
  group_by(sampleName, sseqid) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sampleName, .ref = sseqid, .contig = qseqid)

final <- final %>%
  left_join(
    contig_by_ref %>% rename(denovo_major_ref = .ref, denovo_major_contig = .contig),
    by = c("sampleName", "denovo_major_ref")
  ) %>%
  left_join(
    contig_by_ref %>% rename(denovo_minor_ref = .ref, denovo_minor_contig = .contig),
    by = c("sampleName", "denovo_minor_ref")
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
# novo (background/refuted_denovo). These per-sample booleans feed the pmap_chr
# below alongside overall_sample_call + gate_flag + the subtype-match columns.
#
# CR-01/WR-05 (12-REVIEW): the old any_uncorroborated trigger checked
# role_reason == "uncorroborated_kept", a value the Plan-03 evidence_state redesign
# retired entirely (no code path can produce it any more — see the state->role map
# in classify_roles.R::classify_one_sample()), so the old trigger was permanently
# FALSE and its review_flag sentence / call_confidence "provisional" contribution
# could never fire. any_probable_only is the closest new-model equivalent of "kept
# without strong corroboration": a co-infection candidate whose OWN evidence only
# cleared the marginal "probable" band (assembly_support_score in [0.50, 0.72)),
# not the stronger "confirmed" band (>=0.72) — surfacing exactly the confidence
# gradient the continuous-score redesign was meant to preserve.
#
# CR-02 (12-REVIEW): any_refuted_denovo checks role_reason == "refuted_denovo",
# whose ONLY production trigger (denovo_contradicts & quality_fails_state,
# classify_roles.R L482-507) is preempted in every real run by the
# discordant_identity hard gate in classify_one_sample() (L565-570), which runs
# first and reads apply_concordance()'s identical denovo_conflict predicate
# (L106) — apply_concordance() ALWAYS runs before classify_roles() in this file
# (L818/L830). any_refuted_denovo is therefore effectively dead on real data;
# any_discordant_identity is the trigger that ACTUALLY fires for a genuine own-
# assembly-vs-mapping (or GLUE) identity conflict, and is added here as an
# ADDITIONAL review signal alongside (not replacing) the pre-existing
# any_refuted_denovo — see test_classify_roles.R Test20 for the integrated
# (apply_concordance -> score_candidates -> classify_roles) call-order proof.
if (nrow(candidate_support) > 0) {
  role_review <- candidate_support %>%
    group_by(sampleName) %>%
    summarise(
      any_refuted_denovo      = any(role_reason == "refuted_denovo",       na.rm = TRUE),
      any_discordant_identity = any(role_reason == "discordant_identity",  na.rm = TRUE),
      any_probable_only       = any(role == "co-infection" &
                                    !is.na(evidence_state) &
                                    evidence_state == "probable",           na.rm = TRUE),
      dominant_unconfirmed    = any(role == "dominant" &
                                    !is.na(concordance_status) &
                                    concordance_status == "unconfirmed",    na.rm = TRUE),
      .groups = "drop"
    )
} else {
  role_review <- tibble(
    sampleName               = character(),
    any_refuted_denovo       = logical(),
    any_discordant_identity  = logical(),
    any_probable_only        = logical(),
    dominant_unconfirmed     = logical()
  )
}

final <- final %>%
  left_join(role_review, join_by(sampleName))

# EVID-06 hybrid split (RESEARCH Pattern 3): the PER-CANDIDATE review reasons
# (weak/probable/refuted/discordant evidence states + the D-08 demotions) are built
# row-wise on the long candidate_support frame via the pure Plan-01 helper
# candidate_review_fragment(), then collapsed to ONE fragment string per sample —
# arranged by candidate_rank first for a DETERMINISTIC order (RESEARCH anti-pattern:
# non-deterministic collapse), NA-filtered before paste (Pitfall 2: never render a
# literal "NA"), joined with " | ". A clean sample (every fragment NA) yields NA. This
# collapsed column is transient: sample_review_message() merges it into review_flag
# below, then it is dropped in the same select() as the role_review helper booleans.
if (nrow(candidate_support) > 0) {
  candidate_flag_review <- candidate_support %>%
    mutate(.frag = pmap_chr(
      list(role, role_reason, evidence_state, concordance_status,
           candidate_rank, candidate_ref, candidate_subtype,
           assembly_support_best_contig_length,
           assembly_support_best_contig_pident,
           assembly_support_best_contig_kmer_cov),
      candidate_review_fragment
    )) %>%
    arrange(sampleName, candidate_rank) %>%
    group_by(sampleName) %>%
    summarise(
      candidate_flag_fragment = {
        f <- .frag[!is.na(.frag)]
        if (length(f) == 0) NA_character_ else paste(f, collapse = " | ")
      },
      .groups = "drop"
    )
} else {
  candidate_flag_review <- tibble(
    sampleName              = character(),
    candidate_flag_fragment = character()
  )
}

final <- final %>%
  left_join(candidate_flag_review, join_by(sampleName)) %>%
  # EVID-06/D-11: name of the dominant candidate's reference, consumed by
  # sample_review_message() below and dropped in the same select() as the helper booleans.
  left_join(role_dominant_ref, join_by(sampleName))

# Phase-10 rescue_flag (D-08): a per-sample boolean, NEW and INDEPENDENT of review_flag.
# TRUE iff ANY candidate of the sample had its mapped reference REPLACED by a de-novo
# rescue (rescued_from non-NA). Copies the role_review group_by/summarise/left_join
# structure; placed ALONGSIDE review_flag (NOT coupled into its pmap_chr logic, and
# classify_roles.R / apply_concordance() are untouched, per D-08). NA-safe: a no-rescue
# / skip-assembly sample yields FALSE; a no-candidate batch yields a typed zero-row frame.
if (nrow(candidate_support) > 0) {
  rescue_review <- candidate_support %>%
    group_by(sampleName) %>%
    summarise(
      rescue_flag = any(!is.na(rescued_from), na.rm = TRUE),
      .groups = "drop"
    )
} else {
  rescue_review <- tibble(
    sampleName  = character(),
    rescue_flag = logical()
  )
}

final <- final %>%
  left_join(rescue_review, join_by(sampleName)) %>%
  # A sample with no candidate_support row (left_join NA) is, by definition, not rescued.
  mutate(rescue_flag = if_else(is.na(rescue_flag), FALSE, rescue_flag))

# Phase-10 rescue_effect (handoff §4b): a per-sample categorical that reports WHERE a
# surviving de-novo rescue landed — on the dominant/Major slot vs a minor slot vs
# nowhere. Derived ONLY from data summarize.R already holds (candidate_support role +
# rescued_from), so no new channel input is added (avoids SUMMARIZE channel-arity
# fragility). NOTE: rescue_effect reflects ONLY rescues that SURVIVED into the final
# call; the standalone {prefix}.rescue_audit.csv (rescue_evaluation.R) is authoritative
# and additionally captures dropped_collapse / dropped_cap events. Mirrors the
# rescue_review group_by/summarise/left_join structure.
if (nrow(candidate_support) > 0) {
  rescue_effect_review <- candidate_support %>%
    group_by(sampleName) %>%
    summarise(
      major_changed = any(role == "dominant" & !is.na(rescued_from), na.rm = TRUE),
      minor_changed = any(!is.na(rescued_from), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(rescue_effect = case_when(
      major_changed ~ "major_ref_changed",
      minor_changed ~ "minor_ref_changed",
      TRUE          ~ "none"
    )) %>%
    select(sampleName, rescue_effect)
} else {
  rescue_effect_review <- tibble(
    sampleName    = character(),
    rescue_effect = character()
  )
}

final <- final %>%
  left_join(rescue_effect_review, join_by(sampleName)) %>%
  # A sample with no candidate_support row (left_join NA) had no surviving rescue.
  mutate(rescue_effect = if_else(is.na(rescue_effect), "none", rescue_effect))

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
#   3. monoinfection AND denovo_minor_subtype is non-NA AND its genotype differs from
#      the major genotype — de novo found evidence of a second, different-genotype strain
#      that the role classifier demoted to background (possible missed co-infection)
#   4. any_refuted_denovo — a minor candidate refuted by de novo; likely single infection
#      (CR-02: on real data this specific trigger is effectively dead — see #4a)
#   4a. any_discordant_identity — a candidate's mapping identity genuinely conflicts
#      with its own GLUE and/or de novo assembly identity (role_reason ==
#      "discordant_identity"). This is the trigger that ACTUALLY fires in production
#      for a genuine own-assembly-vs-mapping contradiction: apply_concordance()
#      always runs before classify_roles() in this file, so its discordant_identity
#      hard gate pre-empts the refuted_denovo path above before it can ever be
#      reached (CR-02, 12-REVIEW; see test_classify_roles.R Test20 for the proof).
#      Added as an ADDITIONAL trigger alongside any_refuted_denovo, not a replacement.
#   5. any_probable_only — a co-infection candidate's own evidence only cleared the
#      marginal "probable" band (assembly_support_score 0.50-0.72), not "confirmed"
#      (>=0.72); warrants analyst review of QC plots / contigs (CR-01/WR-05, replaces
#      the pre-Phase-12 any_uncorroborated trigger, which checked for the retired
#      role_reason == "uncorroborated_kept" and could never fire)
#   6. overall_sample_call == "indeterminate" — no candidate passed the major-gate
#   7. gate_flag != "ok" — major failed the first-mapping quality thresholds
#   8. rescue_effect == "major_ref_changed" — the de-novo rescue OVERRODE the
#      dominant/Major reference chosen by first-mapping. This is the single
#      highest-stakes automated decision in the pipeline (it silently replaces the
#      primary call for the dominant strain), so it ALWAYS warrants human
#      adjudication — see rescue_audit.csv for the from/to and the trigger. Note
#      rescue_flag (D-08) stays independent; only a rescue that landed on the MAJOR
#      slot is coupled into review_flag here (a minor-slot rescue is surfaced via
#      rescue_effect but does not on its own force review).
# (Earlier versions emitted semicolon-separated reason codes; rewritten to full
# sentences in commit ff12009; rewired onto roles in Phase 8 / D-15.)
#
# MultiQC orange-highlight note: in assets/multiqc_config.yml the results_summary
# custom_data block includes a cond_formatting_rules entry for this column that
# colours any non-NA value orange (warn class). See the pconfig.cond_formatting_rules
# key added there. If that config is absent (older deployments), MultiQC falls back
# to plain text — the column is still useful as a text summary.
# EVID-06 (D-05/D-09/D-10/D-11): the sample-level review_flag is now built by the pure,
# unit-tested Plan-01 helper sample_review_message() (classify_roles.R) instead of an inline
# anonymous closure (RESEARCH Pitfall 4: do not widen an 11-arg positional closure). It keeps
# the D-10 triggers generic (is_indet / gate_flag / is_indet_dom), enriches the monoinfection
# subtype-conflict + different-genotype-contig triggers with the NAMED dominant candidate
# (rank+ref) and the actual conflicting subtype values, resolves D-11 by naming the dominant
# candidate for dominant_unconfirmed and rescue_effect == "major_ref_changed", and merges the
# pre-collapsed per-candidate fragment (candidate_flag_fragment, built above via
# candidate_review_fragment) with " | ". It returns NA_character_ when nothing fires so a clean
# sample stays NA (the MultiQC cond_formatting_rules colours any non-NA review_flag orange).
# Args are threaded POSITIONALLY to match sample_review_message()'s signature.

# 260803-ogc: the CONTIG-LENGTH leg of the different-genotype-contig review trigger.
# Masked into a TRANSIENT column rather than applied to denovo_minor_subtype itself:
# that column is emitted in Summary.csv, feeds Minor_evidence via build_evidence()
# below, and already produced denovo_minor_subtype_match upstream — overwriting it
# would corrupt three reported values to change one sentence. The genotype-difference
# and 2k1b legs live in offgenotype_contig_reviewable() (classify_roles.R) and run
# inside sample_review_message(); only the length leg needs the param and the length
# column, both of which are in scope here. denovo_minor_contig_length is joined into
# `final` at the df_denovo left_join above, so no new read or join is introduced.
final <- final %>%
  mutate(
    denovo_minor_subtype_reviewable = if_else(
      coalesce(denovo_minor_contig_length, Inf) >= review_min_offgenotype_contig_length,
      denovo_minor_subtype,
      NA_character_
    )
  )

# 260803-ogc option C: the MEASURED EVIDENCE behind that sentence. Until now it named
# a subtype and nothing else, then asked a human to review it — while the aligned
# length, identity and k-mer coverage of the very contig it was talking about existed
# ONLY in blastparse/<sample>.assembly_support.csv. They never reached Summary.csv,
# because join_assembly_support() attaches support at CANDIDATE grain and the
# off-genotype contig is by definition not a candidate. So the flag that asks for a
# judgement withheld every number needed to make it.
#
# The fix is one join of the SAME support frame at the off-genotype subtype's grain.
# support_df is already read and typed above (line ~772) and is untouched since, so
# no new file read is introduced.
#
# na_matches = "never": denovo_minor_subtype is NA on most samples, and dplyr's
# default would happily match those against any NA subtype on the support side.
# slice_max(with_ties = FALSE) is a defensive many-to-many guard in the JMAP-03 /
# T-dpj-02 idiom — blast_parse.R emits one row per (sample, subtype) today, and this
# keeps a future duplicate from multiplying rows in `final`.
#
# NOTE the deliberate asymmetry: the length GATE above uses denovo_minor_contig_length
# (from blastparse.csv, per-reference, the validated sweep column), while the DISPLAYED
# metrics come from this per-subtype assembly_support row so that all four numbers are
# a coherent same-row pair. The two length definitions were identical on 68 of 72
# flagged samples (max difference 176 bp), which is immaterial in a human-readable
# sentence but would matter if the gate were switched to this column.
offgeno_support <- support_df %>%
  group_by(sampleName, subtype) %>%
  slice_max(best_contig_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    sampleName,
    denovo_minor_subtype  = subtype,
    offgeno_contig_length = best_contig_length,
    offgeno_contig_aln    = best_contig_aln_length,
    offgeno_contig_pident = best_contig_pident,
    offgeno_contig_kmer   = best_contig_kmer_cov
  )

final <- final %>%
  left_join(offgeno_support, by = c("sampleName", "denovo_minor_subtype"),
            na_matches = "never") %>%
  mutate(offgeno_note = pmap_chr(
    list(offgeno_contig_length, offgeno_contig_aln, offgeno_contig_pident, offgeno_contig_kmer),
    offgenotype_contig_note
  ))

final <- final %>%
  mutate(review_flag = pmap_chr(
    list(
      overall_sample_call,
      denovo_major_subtype_match,
      denovo_minor_subtype_match,
      gate_flag,
      denovo_minor_subtype_reviewable,
      denovo_major_subtype,
      Major_subtype,
      rescue_effect,
      dominant_unconfirmed,
      dominant_cand_rank,
      dominant_cand_ref,
      candidate_flag_fragment,
      offgeno_note
    ),
    sample_review_message
  )) %>%
  # call_confidence (evidence-tier axis) — an EXPLICIT rollup of the confidence
  # signals that are otherwise scattered across review_flag / gate_flag /
  # rescue_effect / subtype-match / role columns. It does NOT introduce new science:
  # every input below is already computed upstream. The point is to separate "what
  # the data indicate" (overall_sample_call — kept concrete) from "how much to trust
  # it" (this column), so a clean call reads differently from a marginal one at a
  # glance without a reader reverse-engineering the flag columns. Deliberately a
  # small controlled vocabulary; tiers are ordered most-severe-first (first match
  # wins), so a single hard signal caps the tier regardless of softer ones.
  #
  #   indeterminate — no actionable call (no candidate passed the gate / untypable).
  #   review        — a HARD conflict that should block auto-reporting: major failed
  #                   mapping QC, major subtype conflict (de novo vs mapping), the
  #                   rescue overrode the Major reference, or dominance is ambiguous.
  #   provisional   — a SOFT caveat worth noting but not blocking: dominant identity
  #                   uncorroborated, co-infection kept without corroboration, a minor
  #                   refuted by de novo, a minor-subtype conflict, or a minor-slot
  #                   rescue.
  #   high          — none of the above fired AND no review note; the call stands
  #                   on clean evidence.
  #
  # Invariant (one-directional, enforced by construction): a sample with a non-empty
  # review_flag is NEVER `high` — the final `!is.na(review_flag)` clause below
  # demotes any remaining review-noted sample to at least `provisional`. This catches
  # review_flag triggers that have no dedicated column signal here (e.g. the
  # "monoinfection but de novo found a different-genotype contig" trigger, where the
  # Minor slot — and thus denovo_minor_subtype_match — is NA). So `high` always has
  # an empty review_flag. The converse does NOT hold: some `provisional` samples
  # carry no prose note (a minor-slot rescue, an uncorroborated co-infection
  # dominant), which is intended — `provisional` is a softer bucket than a sentence.
  mutate(
    call_confidence = case_when(
      is.na(overall_sample_call) |
        overall_sample_call %in% c("indeterminate", "untypable")        ~ "indeterminate",
      (!is.na(gate_flag) & gate_flag != "ok") |
        (!is.na(denovo_major_subtype_match) & denovo_major_subtype_match == "NO") |
        (!is.na(rescue_effect) & rescue_effect == "major_ref_changed") |
        overall_sample_call == "co-infection (indeterminate dominance)"  ~ "review",
      coalesce(dominant_unconfirmed, FALSE) |
        coalesce(any_probable_only, FALSE) |
        coalesce(any_refuted_denovo, FALSE) |
        coalesce(any_discordant_identity, FALSE) |
        (!is.na(denovo_minor_subtype_match) & denovo_minor_subtype_match == "NO") |
        (!is.na(rescue_effect) & rescue_effect == "minor_ref_changed") |
        !is.na(review_flag)                                              ~ "provisional",
      TRUE                                                               ~ "high"
    )
  ) %>%
  # Drop the per-sample role-review helper booleans now they have been consumed.
  # Also drop the two EVID-06 transients: candidate_flag_fragment (per-sample collapsed
  # per-candidate fragment) and dominant_cand_ref (dominant reference name) — both merged
  # into review_flag by sample_review_message() above and not part of the emitted schema.
  select(-any_refuted_denovo, -any_discordant_identity, -any_probable_only, -dominant_unconfirmed,
         -candidate_flag_fragment, -dominant_cand_ref,
         # 260803-ogc transients: the length-masked copy of denovo_minor_subtype and the
         # off-genotype contig metrics + their rendered clause, all consumed by
         # sample_review_message() above. Dropped here for the same CR-02 reason as the
         # two transients beside them — not part of the emitted schema. The metrics
         # themselves remain available per sample in blastparse/*.assembly_support.csv,
         # and are now surfaced in prose inside review_flag.
         -denovo_minor_subtype_reviewable, -offgeno_note,
         -offgeno_contig_length, -offgeno_contig_aln,
         -offgeno_contig_pident, -offgeno_contig_kmer)

# Shorthand aliases surfaced near the front of Summary.csv for at-a-glance reading.
# Pure verbatim copies of the existing, buried Major_subtype / Minor_subtype — no
# subtype is recomputed or re-derived here (BM2-01).
final <- final %>%
  mutate(
    Major = Major_subtype,
    Minor = Minor_subtype
  )

# Per-strain evidence basis (Major_evidence / Minor_evidence). A single compact
# string per strain answering the reviewer's real question — "on what evidence is
# THIS strain being called?" — so the basis for a call can be read in one cell
# instead of cross-referencing the ~dozen scattered Major_*/Minor_* metric columns.
# Pure presentation: every token is a value already computed upstream; nothing is
# re-derived. Optional columns (consensus identity, present only when the distance
# leg ran) are NA-filled first so the row-wise builder never errors, and each token
# is emitted only when its source value is non-NA. A strain with no reference (e.g.
# the Minor slot of a monoinfection) yields NA — no empty scaffold string.
evidence_cols <- c(
  "Major_reference", "Reads_nodup_mapped_major", "Major_cov_breadth_min_10",
  "denovo_major_subtype", "Major_consensus_similarity_pct",
  "Minor_reference", "Reads_nodup_mapped_minor", "Minor_cov_breadth_min_10",
  "denovo_minor_subtype", "Minor_consensus_similarity_pct",
  # EVID-05/D-04: contig-corroboration source columns for the extended build_evidence().
  # NA-filled here so a batch missing any of them NA-fills instead of erroring row-wise.
  "Major_evidence_state", "Major_best_contig_length", "Major_best_contig_pident",
  "Major_best_contig_kmer_cov",
  "Minor_evidence_state", "Minor_best_contig_length", "Minor_best_contig_pident",
  "Minor_best_contig_kmer_cov"
)
for (col in evidence_cols) {
  if (!col %in% colnames(final)) final <- final %>% add_column(!!col := NA)
}

# Build one strain's evidence string from its (already-computed) parts. `rescued`
# is TRUE when the de-novo rescue reassigned THIS slot's reference (from rescue_effect).
build_evidence <- function(ref, reads, breadth10, dn_sub, cons_pct, rescued,
                           ev_state = NA_character_,
                           contig_len = NA_real_, contig_pid = NA_real_,
                           contig_kmer = NA_real_) {
  if (is.na(ref)) return(NA_character_)
  toks <- ref
  if (!is.na(reads))     toks <- c(toks, paste0(format(round(reads), big.mark = ",", trim = TRUE, scientific = FALSE), " reads (nodup)"))
  if (!is.na(breadth10)) toks <- c(toks, paste0(round(breadth10), "% breadth@10x"))  # cov_breadth_min_10 is already 0-100
  if (!is.na(dn_sub))    toks <- c(toks, paste0("de novo ", dn_sub))
  if (!is.na(cons_pct))  toks <- c(toks, paste0(round(cons_pct, 1), "% consensus id"))
  toks <- c(toks, if (isTRUE(rescued)) "ref REASSIGNED by de-novo rescue" else "ref from first-mapping")
  # EVID-05/D-04: append the contig-corroboration basis (the own best-contig
  # identity/length/k-mer, in contig language) plus the calibrated evidence_state, so
  # Major_evidence / Minor_evidence answer "on what CONTIG evidence?" in the same cell.
  # Purely additive — no existing token above is changed. Each metric emits only when
  # its source value is non-NA (never renders a literal "NA").
  contig_toks <- character(0)
  if (length(contig_len)  == 1 && !is.na(contig_len))  contig_toks <- c(contig_toks, paste0(round(contig_len), " bp contig"))
  if (length(contig_pid)  == 1 && !is.na(contig_pid))  contig_toks <- c(contig_toks, paste0(round(contig_pid, 1), "% contig identity"))
  if (length(contig_kmer) == 1 && !is.na(contig_kmer)) contig_toks <- c(contig_toks, paste0(round(contig_kmer, 2), " k-mer cov"))
  if (length(contig_toks) > 0) toks <- c(toks, paste0("contig support: ", paste(contig_toks, collapse = ", ")))
  if (length(ev_state) == 1 && !is.na(ev_state)) toks <- c(toks, paste0("evidence_state=", ev_state))
  paste(toks, collapse = " | ")
}

final <- final %>%
  mutate(
    Major_evidence = pmap_chr(
      list(Major_reference, Reads_nodup_mapped_major, Major_cov_breadth_min_10,
           denovo_major_subtype, Major_consensus_similarity_pct,
           !is.na(rescue_effect) & rescue_effect == "major_ref_changed",
           Major_evidence_state,
           Major_best_contig_length, Major_best_contig_pident, Major_best_contig_kmer_cov),
      build_evidence
    ),
    Minor_evidence = pmap_chr(
      list(Minor_reference, Reads_nodup_mapped_minor, Minor_cov_breadth_min_10,
           denovo_minor_subtype, Minor_consensus_similarity_pct,
           !is.na(rescue_effect) & rescue_effect == "minor_ref_changed",
           Minor_evidence_state,
           Minor_best_contig_length, Minor_best_contig_pident, Minor_best_contig_kmer_cov),
      build_evidence
    )
  )

# Save dominant-rank lookup before the reorder select drops it (used for
# resistance MQC minor rows — the non-dominant GLUE report for co-infections).
dom_rank_lookup <- if ("dominant_cand_rank" %in% colnames(final)) {
  select(final, sampleName, dominant_cand_rank)
} else {
  tibble(sampleName = character(), dominant_cand_rank = integer())
}

# Reorder columns
final <- final %>%
  select(sampleName,
         Major,
         Minor,
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
         # Explicit evidence-tier axis paired with the call (see mutate above).
         call_confidence,
         Major_role_reference,
         Major_role_subtype,
         Major_dominance_score,
         Major_role_reason,
         # WR-05 (12-REVIEW): the confirmed/probable/weak/refuted evidence band
         # behind Major/Minor_role_reason's "corroborated" bucket.
         Major_evidence_state,
         Minor_role_reference,
         Minor_role_subtype,
         Minor_dominance_score,
         Minor_role_reason,
         Minor_evidence_state,
         denovo_major_subtype,
         denovo_minor_subtype,
         denovo_major_subtype_match,
         denovo_minor_subtype_match,
         review_flag,
         # Per-strain evidence basis strings (built above) — the "why this call".
         Major_evidence,
         Minor_evidence,
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
  # CR-02 (13-REVIEW): Major_/Minor_best_contig_length/pident/kmer_cov are transient
  # scratch inputs consumed by build_evidence() to build the Major_evidence/Minor_evidence
  # prose token (EVID-05/D-04) — they were never meant to be their own visible columns
  # and must be dropped here alongside the other transients, or they leak into Summary.csv.
  select(-any_of(c("Major_minor", "identical_geno", "identical_subgeno", "dominant_cand_rank",
                    "Major_best_contig_length", "Major_best_contig_pident", "Major_best_contig_kmer_cov",
                    "Minor_best_contig_length", "Minor_best_contig_pident", "Minor_best_contig_kmer_cov")))

# Write file
write_csv(final, file = "Summary.csv")

# Write triage table for MultiQC (summary_mqc.tsv)
# Contains only the 6 problem-signal columns plus key context. All columns in the
# resistance section are split out to a separate file (glue_resistance_mqc.tsv) so
# MultiQC can render them as a distinct section with appropriate colour config.
# TSV format avoids quoting issues with commas in review_flag sentences.
# MultiQC config must match: file_format: tsv, fn: "*/summary_mqc.tsv"

triage <- final %>%
  mutate(
    # Signal 1: subtype conflict (mapping vs. de novo)
    subtype_conflict = case_when(
      denovo_major_subtype_match == "NO"  ~ "CONFLICT",
      denovo_major_subtype_match == "YES" ~ "OK",
      TRUE                               ~ NA_character_
    ),
    # Per-row genotype: defaults to the dominant strain; overridden to Minor
    # for the [minor] co-infection rows created below.
    Genotype = Major
  ) %>%
  select(
    sampleName,
    overall_sample_call,                             # placement 2
    call_confidence,                                 # placement 2b
    review_flag,                                     # placement 3
    Genotype,                                        # placement 4
    Major_avg_depth,                                 # placement 6
    subtype_conflict,                                # placement 10
    any_of(c("denovo_major_contig",                  # placement 11 — de novo contig NAME (major strain)
             "denovo_major_ref")),                   # placement 12 — that contig's best BLAST-hit ref
    rescue_flag,                                     # placement 20
    rescue_effect,                                   # placement 20b — where a surviving rescue landed
    total_trimmed_reads,                             # placement 25
    Major_cov_breadth_min_10,                        # placement 30
    Major_cov_breadth_min_5,                         # placement 35
    Percent_reads_mapped_of_trimmed_with_dups_major, # placement 38
    Reads_nodup_mapped_major,                        # placement 40
    percent_mapped_reads_major_firstmapping,         # placement 50
    any_of(c("Major_consensus_similarity_pct", "Major_consensus_n_differences")), # placement 150, 160
    # Minor columns kept temporarily to set Genotype on [minor] rows
    Minor,
    any_of(c("Minor_avg_depth",
             "Reads_nodup_mapped_minor",
             "percent_mapped_reads_minor_firstmapping",
             "Minor_cov_breadth_min_10",
             "Minor_cov_breadth_min_5",
             "Percent_reads_mapped_of_trimmed_with_dups_minor",
             "Minor_consensus_similarity_pct",
             "Minor_consensus_n_differences",
             "denovo_minor_contig",                   # promoted onto the [minor] row's denovo_best_contig
             "denovo_minor_ref"))                     # promoted onto the [minor] row's denovo_best_contig_ref
  )

# Paired metric columns.  major → the column present in triage;
# minor → the counterpart column to promote for the [minor] row;
# generic → the final output column name used in both rows.
col_major   <- c("Major_avg_depth",
                 "Reads_nodup_mapped_major",
                 "percent_mapped_reads_major_firstmapping",
                 "Major_cov_breadth_min_10",
                 "Major_cov_breadth_min_5",
                 "Percent_reads_mapped_of_trimmed_with_dups_major",
                 "Major_consensus_similarity_pct",
                 "Major_consensus_n_differences",
                 "denovo_major_contig",
                 "denovo_major_ref")
col_minor   <- c("Minor_avg_depth",
                 "Reads_nodup_mapped_minor",
                 "percent_mapped_reads_minor_firstmapping",
                 "Minor_cov_breadth_min_10",
                 "Minor_cov_breadth_min_5",
                 "Percent_reads_mapped_of_trimmed_with_dups_minor",
                 "Minor_consensus_similarity_pct",
                 "Minor_consensus_n_differences",
                 "denovo_minor_contig",
                 "denovo_minor_ref")
col_generic <- c("average_depth_0",
                 "Reads_nodup_mapped",
                 "percent_mapped_reads_firstmapping",
                 "cov_breadth_min_10",
                 "cov_breadth_min_5",
                 "Percent_reads_mapped_of_trimmed_with_dups",
                 "pct_similarity_to_nearest_reference",
                 "n_differences_to_nearest_reference",
                 "denovo_best_contig",
                 "denovo_best_contig_ref")

# Restrict to pairs where the major column actually exists in triage
present     <- col_major %in% colnames(triage)
col_major   <- col_major[present]
col_minor   <- col_minor[present]
col_generic <- col_generic[present]

# Expand co-infection samples: add a [major] row and a [minor] row.
# Monoinfection / indeterminate samples: one row, no brackets.
coinf_rows <- triage$overall_sample_call %in% c("co-infection", "co-infection (indeterminate dominance)")

if (any(coinf_rows)) {
  triage_mono <- triage[!coinf_rows, , drop = FALSE]

  # [major] row: identical to the original row but with the bracket label
  triage_coinf_major <- triage[coinf_rows, , drop = FALSE] %>%
    mutate(sampleName = paste0(sampleName, " [major]"))

  # [minor] row: Genotype set to the minor strain's subtype; metric columns
  # populated from the Minor_ counterparts.
  triage_coinf_minor <- triage[coinf_rows, , drop = FALSE] %>%
    mutate(sampleName = paste0(sampleName, " [minor]"),
           Genotype = Minor)
  for (i in seq_along(col_major)) {
    mc <- col_major[i]; nc <- col_minor[i]
    triage_coinf_minor[[mc]] <- if (nc %in% colnames(triage_coinf_minor)) triage_coinf_minor[[nc]] else NA
  }

  # Interleave: sort by base sample name, then [major] before [minor]
  triage <- bind_rows(triage_mono, triage_coinf_major, triage_coinf_minor) %>%
    arrange(str_remove(sampleName, " \\[(major|minor)\\]$"), sampleName)
}

# Drop minor-specific columns and the temporary Minor helper (now consumed by Genotype)
triage <- triage %>% select(-any_of(c(col_minor, "Minor")))

# Rename Major_* columns to final generic names
rename_map <- setNames(col_major, col_generic)
rename_map <- rename_map[rename_map %in% colnames(triage)]
triage <- triage %>% rename(all_of(rename_map)) %>% as.data.frame()

triage_file <- "summary_mqc.tsv"
triage %>% colnames() %>% paste0(collapse = "\t") %>% write_lines(triage_file, append = TRUE)
write_tsv(triage, triage_file, append = TRUE)

# Write resistance table for MultiQC (glue_resistance_mqc.tsv)
# Co-infection samples: two rows — sampleName (major) and "sampleName [minor]".
# Monoinfection samples: one row — sampleName.
# Columns follow the order in the design spec (NS3/4A → NS5A → NS5B drug groups).
# MultiQC config must match: file_format: tsv, fn: "*/glue_resistance_mqc.tsv"

# NS class-summary columns lead each drug group; _short columns are omitted.
resistance_col_order <- c(
  "NS34A",
  "glecaprevir", "glecaprevir_mut",
  "grazoprevir", "grazoprevir_mut",
  "paritaprevir", "paritaprevir_mut",
  "voxilaprevir", "voxilaprevir_mut",
  "NS5A",
  "daclatasvir", "daclatasvir_mut",
  "elbasvir",    "elbasvir_mut",
  "ledipasvir",  "ledipasvir_mut",
  "ombitasvir",  "ombitasvir_mut",
  "pibrentasvir","pibrentasvir_mut",
  "velpatasvir", "velpatasvir_mut",
  "NS5B",
  "dasabuvir",   "dasabuvir_mut",
  "sofosbuvir",  "sofosbuvir_mut"
)
available_res_cols <- resistance_col_order[resistance_col_order %in% colnames(final)]

if (length(available_res_cols) > 0) {
  # Major rows: dominant-strain resistance is already in final after the GLUE join
  major_res <- final %>%
    select(sampleName, all_of(available_res_cols)) %>%
    rename(Sample = sampleName)

  # Minor rows for co-infection samples: use the non-dominant GLUE report.
  # dom_rank_lookup gives the dominant rank; minor rank = 3 - dominant_rank.
  # glue_by_rank (built during the GLUE join section) holds both ranks with .dom_rank.
  minor_res <- tibble()
  if (exists("glue_by_rank") && nrow(glue_by_rank) > 0 && nrow(dom_rank_lookup) > 0) {
    coinf_names <- final %>%
      filter(overall_sample_call == "co-infection") %>%
      pull(sampleName)

    if (length(coinf_names) > 0) {
      minor_rank_df <- dom_rank_lookup %>%
        filter(sampleName %in% coinf_names) %>%
        mutate(minor_rank = 3L - coalesce(dominant_cand_rank, 1L))

      minor_res <- glue_by_rank %>%
        inner_join(minor_rank_df, by = c("sampleName", ".dom_rank" = "minor_rank")) %>%
        select(sampleName, any_of(available_res_cols)) %>%
        mutate(Sample = paste0(sampleName, " [minor]")) %>%
        select(Sample, everything(), -sampleName)
    }
  }

  # Combine; sort so "[minor]" rows appear immediately after their major row
  if (nrow(minor_res) > 0) {
    res_data <- bind_rows(major_res, minor_res) %>%
      arrange(Sample) %>%
      as.data.frame()
  } else {
    res_data <- major_res %>% as.data.frame()
  }

  # Overall Resistance summary column: strongest signal across all status columns.
  # Priority: "Resistance" > "Probable/Possible resistance" > "No resistance" > NA.
  resistance_status_cols <- c(
    "NS34A", "glecaprevir", "grazoprevir", "paritaprevir", "voxilaprevir",
    "NS5A",  "daclatasvir", "elbasvir",    "ledipasvir",   "ombitasvir",
             "pibrentasvir","velpatasvir",
    "NS5B",  "dasabuvir",   "sofosbuvir"
  )
  status_cols_present <- intersect(resistance_status_cols, colnames(res_data))
  if (length(status_cols_present) > 0) {
    status_matrix <- res_data[, status_cols_present, drop = FALSE]
    overall_resistance <- apply(status_matrix, 1, function(vals) {
      vals <- as.character(vals)
      vals <- vals[!is.na(vals) & vals != "NA"]
      if (length(vals) == 0) return(NA_character_)
      if (any(vals == "Resistance")) return("Resistance")
      if (any(grepl("Probable resistance|Possible resistance", vals, ignore.case = TRUE))) return("Probable/Possible resistance")
      if (any(vals == "No resistance")) return("No resistance")
      return(NA_character_)
    })
    res_data <- res_data %>%
      mutate(Resistance = overall_resistance) %>%
      select(Sample, Resistance, everything())
  }

  res_file <- "glue_resistance_mqc.tsv"
  res_data %>% colnames() %>% paste0(collapse = "\t") %>% write_lines(res_file, append = TRUE)
  write_tsv(res_data, res_file, append = TRUE)
}

