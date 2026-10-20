# background_signals.R ----------------------------------------------------
# Two sample-level signals that a co-infection minor may be background reads
# rather than a second strain. Neither needs any other sample of the run, so
# replicate samples, mixed runs and sources outside the run do not matter.
#
#   1. Read fraction. Background bleed (index hopping, demultiplexing bleed,
#      carry-over) is roughly a constant FRACTION of reads, while probe-capture
#      depth makes absolute read floors easy to pass. A minor's share of the
#      sample's candidate reads is therefore the depth-invariant quantity.
#
#   2. Duplicate discordance. Within one library every molecule is amplified and
#      sequenced alike, so a genuine minor shares the library's PCR-duplicate rate.
#      Reads that arrived from ANOTHER library after its PCR are a sparse draw from
#      a large pool and are almost all unique. A minor whose contigs are far less
#      duplicated than the dominant genotype's is therefore foreign reads. The test
#      is one-directional: contamination that entered BEFORE library PCR is
#      amplified like the sample's own and looks concordant.
#
# Measured on the Thomson 2016 MiSeq run: background minors at 0.17-1.7 % of
# reads with contigs 0-15 % duplicated vs 76-96 % for the sample's own genome
# (hcvtyper_handoff_dedup_before_assembly_2026-10-16.md).
#
# Neither signal changes a call. summarize.R adds the numbers as columns and a
# review sentence. Pure functions; typed-empty in, typed-empty out.
# genotype_from_subtype() comes from genotype_utils.R, sourced before this file.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Share of the sample's candidate reads that each candidate holds, on the
# competitive joint-mapping BAM with duplicates included (raw abundance, which is
# what a bleed rate describes). `reads` has sampleName, candidate_rank,
# reads_withdup. Fractions sum to 1 within a sample.
candidate_read_fractions <- function(reads) {
  if (nrow(reads) == 0) {
    return(tibble(sampleName = character(), candidate_rank = integer(),
                  candidate_read_fraction = double()))
  }
  reads %>%
    filter(!is.na(candidate_rank)) %>%
    group_by(sampleName) %>%
    mutate(candidate_read_fraction = if (sum(reads_withdup, na.rm = TRUE) > 0)
                                       reads_withdup / sum(reads_withdup, na.rm = TRUE)
                                     else NA_real_) %>%
    ungroup() %>%
    select(sampleName, candidate_rank, candidate_read_fraction)
}

# Read every <sample>.contig_dups.tsv in `dir` (CONTIG_DUPLICATES output).
read_contig_dups <- function(dir) {
  empty <- tibble(sampleName = character(), qseqid = character(),
                  fragments = double(), duplicates = double())
  f <- list.files(dir, pattern = "\\.contig_dups\\.tsv$", full.names = TRUE)
  if (length(f) == 0) return(empty)
  map_dfr(f, function(x) {
    d <- read_tsv(x, col_types = cols(contig = col_character(),
                                      fragments = col_double(),
                                      duplicates = col_double()),
                  show_col_types = FALSE)
    if (nrow(d) == 0) return(empty)
    d %>% transmute(sampleName = str_remove(basename(x), "\\.contig_dups\\.tsv$"),
                    qseqid = contig, fragments, duplicates)
  })
}

# Pool contig duplicate counts by the genotype of each contig's top BLAST hit
# (lowest e-value, then highest bitscore — the blast_parse.R rule). Contigs with
# no HCV hit are left out. Returns one row per (sampleName, genotype).
genotype_dup_rates <- function(contig_dups, blast_out) {
  empty <- tibble(sampleName = character(), genotype = character(),
                  dup_fragments = double(), dup_duplicates = double(),
                  dup_rate = double())
  if (nrow(contig_dups) == 0 || nrow(blast_out) == 0) return(empty)
  top <- blast_out %>%
    filter(!is.na(subtype)) %>%
    arrange(evalue, desc(bitscore)) %>%
    group_by(sampleName, qseqid) %>% slice(1) %>% ungroup() %>%
    transmute(sampleName, qseqid, genotype = genotype_from_subtype(subtype))
  contig_dups %>%
    inner_join(top, by = c("sampleName", "qseqid")) %>%
    group_by(sampleName, genotype) %>%
    summarise(dup_fragments  = sum(fragments),
              dup_duplicates = sum(duplicates),
              .groups = "drop") %>%
    mutate(dup_rate = if_else(dup_fragments > 0, dup_duplicates / dup_fragments, NA_real_))
}

# Per candidate: its genotype's duplicate rate, the library rate (the dominant
# candidate's genotype), their ratio, and whether the candidate is discordant.
#
# `candidates` has sampleName, candidate_rank, candidate_ref and role; `dominant`
# has sampleName and dominant_ref (the Major-slot candidate). dup_discordant is
#   TRUE  when the library is saturated enough to show duplication
#         (library rate >= min_library_dup), the candidate has >= min_fragments,
#         and its rate is <= max_ratio x the library rate;
#   FALSE when those first two hold and the rate is above that;
#   NA    when the test has no power (unsaturated library, too few fragments,
#         no contigs) or the candidate IS the dominant genotype.
candidate_dup_signal <- function(candidates, dominant, geno_dups,
                                 min_library_dup = 0.5, min_fragments = 20,
                                 max_ratio = 0.5) {
  if (nrow(candidates) == 0) {
    return(tibble(sampleName = character(), candidate_rank = integer(),
                  candidate_dup_rate = double(), candidate_dup_fragments = double(),
                  library_dup_rate = double(), dup_ratio = double(),
                  dup_discordant = logical()))
  }
  geno_of <- function(ref) genotype_from_subtype(str_extract(ref, "^[^_]+"))
  lib <- dominant %>%
    transmute(sampleName, .lib_geno = geno_of(dominant_ref)) %>%
    left_join(geno_dups %>% select(sampleName, .lib_geno = genotype,
                                   library_dup_rate = dup_rate),
              by = c("sampleName", ".lib_geno"))
  candidates %>%
    transmute(sampleName, candidate_rank, .geno = geno_of(candidate_ref)) %>%
    left_join(geno_dups %>% select(sampleName, .geno = genotype,
                                   candidate_dup_fragments = dup_fragments,
                                   candidate_dup_rate = dup_rate),
              by = c("sampleName", ".geno")) %>%
    left_join(lib, by = "sampleName") %>%
    mutate(
      dup_ratio = if_else(!is.na(library_dup_rate) & library_dup_rate > 0,
                          candidate_dup_rate / library_dup_rate, NA_real_),
      .testable = !is.na(.lib_geno) & .geno != .lib_geno &
                  !is.na(library_dup_rate) & library_dup_rate >= min_library_dup &
                  !is.na(candidate_dup_fragments) & candidate_dup_fragments >= min_fragments,
      dup_discordant = if_else(.testable, dup_ratio <= max_ratio, NA)
    ) %>%
    select(sampleName, candidate_rank, candidate_dup_rate, candidate_dup_fragments,
           library_dup_rate, dup_ratio, dup_discordant)
}

# The review sentence for one co-infection candidate, or NA when neither signal
# fires. Carries the measured numbers so it can be judged without a re-run.
background_note <- function(candidate_rank, candidate_ref, read_fraction,
                            dup_rate, library_dup_rate, dup_discordant,
                            min_read_fraction) {
  low_frac <- !is.na(read_fraction) && read_fraction < min_read_fraction
  discord  <- isTRUE(dup_discordant)
  if (!low_frac && !discord) return(NA_character_)
  pct <- function(x) paste0(formatC(100 * x, format = "f", digits = if (x < 0.1) 2 else 1), "%")
  parts <- character(0)
  if (low_frac) {
    parts <- c(parts, sprintf("holds %s of the sample's HCV reads (below %s)",
                              pct(read_fraction), pct(min_read_fraction)))
  }
  if (discord) {
    parts <- c(parts, sprintf(
      "its contigs are %s PCR duplicates against %s for the dominant genotype, so its reads did not go through this library's PCR",
      pct(dup_rate), pct(library_dup_rate)))
  }
  sprintf(
    "Co-infection candidate %s (%s) %s — consistent with background reads from another sample (index hopping or carry-over); check contamination_check/ before reporting.",
    candidate_rank, candidate_ref, paste(parts, collapse = ", and "))
}
