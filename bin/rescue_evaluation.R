#!/usr/bin/env Rscript

# rescue_evaluation.R ----------------------------------------------------------
# De-novo subtype-rescue evaluator (Phase 10, denovo-subtype-rescue).
#
# Reads the long-format candidates table (one row per ranked candidate) together
# with the de-novo BLAST top-hit summary (blastparse.csv) and the four-floor
# assembly_support.csv, and decides — per candidate slot — whether the reference
# the reads were mapped to should be REPLACED by the de-novo-derived reference
# (a "rescue"). A rescue fires when the candidate subtype disagrees with the
# de-novo top-hit subtype AND the de-novo contig clears all four quality floors,
# with two special cases:
#
#   * 2k1b special rule (D-03): a 2k1b candidate (or 2k1b de-novo hit) over a
#     genotype-2 contig meeting the floors rescues to that genotype-2 reference,
#     regardless of the normal mismatch check.
#   * 1a/1b boundary (D-04): when both the candidate and de-novo subtypes are in
#     {1a, 1b}, the contig-length floor is the stricter rescue_1a1b_length.
#
# On rescue: rescued_from records the ORIGINAL candidate ref, candidate_ref is
# overwritten with the rescue reference, rescue_trigger gets a human-readable
# evidence string, confirmation_status is forced to "pass" (D-05; the downstream
# hcvtyper.nf filter drops non-pass rows), and the rescue reference is re-written
# as {prefix}.{ref}_cand{rank}.fa via the verbatim write_ref_fasta() membership
# guard copied from blast_parse.R (V5 / T-10-02: a ref absent from
# params.references writes NO FASTA → candidate guarded out, never a wrong map).
#
# Empty / skip-assembly inputs (zero-row blastparse/support) are handled by typed
# -empty tibbles so the joins yield zero matches → candidates pass through with
# rescued_from / rescue_trigger NA-filled, and the script ALWAYS exits 0
# (T-10-03 DoS guard). The output preserves the 8 original candidate columns and
# adds exactly rescued_from + rescue_trigger.
#
# Positional args (PARSEFIRSTMAPPING convention):
#   prefix candidates_csv blastparse_csv support_csv references \
#   rescue_min_length rescue_min_pident rescue_min_aln_length \
#   rescue_min_kmer_cov rescue_1a1b_length
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(seqinr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: rescue_evaluation.R <prefix> <candidates.csv> <blastparse.csv> ",
       "<assembly_support.csv> <references.fa> <rescue_min_length> ",
       "<rescue_min_pident> <rescue_min_aln_length> <rescue_min_kmer_cov> ",
       "<rescue_1a1b_length>")
}

prefix          <- args[1]
candidates_csv  <- args[2]
blastparse_csv  <- args[3]
support_csv     <- args[4]
references      <- args[5]
rescue_min_length     <- as.numeric(args[6])
rescue_min_pident     <- as.numeric(args[7])
rescue_min_aln_length <- as.numeric(args[8])
rescue_min_kmer_cov   <- as.numeric(args[9])
rescue_1a1b_length    <- as.numeric(args[10])

# --- Subtype-token helper -----------------------------------------------------
# Extract the subtype token from a {subtype}_{accession} reference name
# (bin/summarize.R:1203 convention). NA-safe.
subtype_of <- function(ref) {
  ifelse(is.na(ref), NA_character_, str_extract(ref, "^[^_]+"))
}

# --- Typed-empty guards (T-10-03 DoS guard) -----------------------------------
# Build typed-empty frames mirroring blast_parse.R:184-197 so downstream joins
# never abort on zero-row / unreadable input — candidates simply pass through.
empty_candidates <- tibble(
  sample              = character(0),
  candidate_rank      = integer(0),
  candidate_ref       = character(0),
  candidate_subtype   = character(0),
  candidate_genotype  = character(0),
  candidate_reads     = double(0),
  candidate_cov       = double(0),
  confirmation_status = character(0)
)
empty_blastparse <- tibble(
  sample              = character(0),
  major_ref           = character(0),
  major_contig_length = double(0),
  minor_ref           = character(0),
  minor_contig_length = double(0)
)
empty_support <- tibble(
  sample                 = character(0),
  subtype                = character(0),
  best_contig_length     = double(0),
  best_contig_pident     = double(0),
  best_contig_aln_length = double(0),
  best_contig_kmer_cov   = double(0)
)

read_csv_guarded <- function(path, empty_tbl, col_types) {
  if (is.na(path) || !file.exists(path)) return(empty_tbl)
  tbl <- tryCatch(
    suppressWarnings(read_csv(path, col_types = col_types, progress = FALSE)),
    error = function(e) NULL
  )
  if (is.null(tbl) || nrow(tbl) == 0) return(empty_tbl)
  # Ensure all expected columns are present; missing ones become typed NA.
  missing <- setdiff(colnames(empty_tbl), colnames(tbl))
  for (m in missing) tbl[[m]] <- empty_tbl[[m]][NA_integer_][seq_len(nrow(tbl))]
  tbl
}

candidates <- read_csv_guarded(
  candidates_csv, empty_candidates,
  cols(
    sample              = col_character(),
    candidate_rank      = col_integer(),
    candidate_ref       = col_character(),
    candidate_subtype   = col_character(),
    candidate_genotype  = col_character(),
    candidate_reads     = col_double(),
    candidate_cov       = col_double(),
    confirmation_status = col_character(),
    .default            = col_character()
  )
)
blastparse <- read_csv_guarded(
  blastparse_csv, empty_blastparse,
  cols(
    sample              = col_character(),
    major_ref           = col_character(),
    major_contig_length = col_double(),
    minor_ref           = col_character(),
    minor_contig_length = col_double(),
    .default            = col_character()
  )
)
support <- read_csv_guarded(
  support_csv, empty_support,
  cols(
    sample                 = col_character(),
    subtype                = col_character(),
    best_contig_length     = col_double(),
    best_contig_pident     = col_double(),
    best_contig_aln_length = col_double(),
    best_contig_kmer_cov   = col_double(),
    .default               = col_character()
  )
)

# --- Reference FASTA re-extraction (COPY VERBATIM from blast_parse.R:319-327) --
# The `ref %in% names(ref_fa)` membership guard is the V5 / T-10-02 integrity
# control: a reference absent from params.references writes NO FASTA.
ref_fa <- if (!is.na(references) && file.exists(references)) {
  tryCatch(read.fasta(file = references), error = function(e) list())
} else {
  list()
}
write_ref_fasta <- function(ref_name, tag) {
  if (!is.na(ref_name) && ref_name %in% names(ref_fa)) {
    write.fasta(
      sequences = ref_fa[ref_name],
      names     = ref_name,
      file.out  = paste0(prefix, ".", ref_name, "_", tag, ".fa")
    )
  }
}

# --- Floor evaluation ---------------------------------------------------------
# >  for length / pident, >= for aln_length / kmer_cov, exactly per D-02.
floors_ok <- function(srow, length_floor) {
  if (nrow(srow) == 0) return(FALSE)
  isTRUE(
    srow$best_contig_length[1]     >  length_floor          &
    srow$best_contig_pident[1]     >  rescue_min_pident     &
    srow$best_contig_aln_length[1] >= rescue_min_aln_length &
    srow$best_contig_kmer_cov[1]   >= rescue_min_kmer_cov
  )
}

# Find the genotype of a subtype token (leading character, e.g. "2a" -> "2").
genotype_of <- function(subtype) {
  ifelse(is.na(subtype), NA_character_, str_sub(subtype, 1, 1))
}

# Per-sample set of subtypes that look 2k1b (candidate OR de-novo hit) for D-03.
sample_has_2k1b <- function(sample_id) {
  cand_subs   <- candidates %>% filter(sample == sample_id) %>% pull(candidate_subtype)
  bp          <- blastparse %>% filter(sample == sample_id)
  denovo_subs <- c(subtype_of(bp$major_ref), subtype_of(bp$minor_ref))
  any(c(cand_subs, denovo_subs) == "2k1b", na.rm = TRUE)
}

# --- Per-candidate rescue evaluation -----------------------------------------
evaluate_row <- function(row) {
  sample_id <- row$sample
  rank      <- row$candidate_rank
  orig_ref  <- row$candidate_ref
  cand_sub  <- subtype_of(orig_ref)

  bp <- blastparse %>% filter(sample == sample_id)

  # De-novo top-hit ref for this slot: major for rank 1, minor for rank 2.
  denovo_ref <- NA_character_
  if (nrow(bp) > 0) {
    denovo_ref <- if (rank == 1L) bp$major_ref[1]
                  else if (rank == 2L) bp$minor_ref[1]
                  else NA_character_
  }
  denovo_sub <- subtype_of(denovo_ref)

  no_rescue <- list(rescued_from = NA_character_, rescue_ref = NA_character_,
                    rescue_trigger = NA_character_)

  # ---- 2k1b special rule (D-03) ---------------------------------------------
  # If 2k1b appears anywhere for the sample and a genotype-2 contig meets the
  # floors, rescue this slot to the corresponding genotype-2 reference.
  if (sample_has_2k1b(sample_id)) {
    g2_support <- support %>%
      filter(sample == sample_id, genotype_of(subtype) == "2") %>%
      arrange(desc(best_contig_length))
    if (nrow(g2_support) > 0) {
      g2_row <- g2_support[1, ]
      if (floors_ok(g2_row, rescue_min_length)) {
        # Pick the de-novo ref whose subtype matches this genotype-2 contig.
        g2_sub <- g2_row$subtype[1]
        rescue_ref <- if (!is.na(denovo_sub) && denovo_sub == g2_sub) denovo_ref
                      else {
                        cand_match <- c(bp$major_ref, bp$minor_ref)
                        cand_match <- cand_match[subtype_of(cand_match) == g2_sub]
                        cand_match <- cand_match[!is.na(cand_match)]
                        if (length(cand_match) > 0) cand_match[1] else denovo_ref
                      }
        if (!is.na(rescue_ref) && rescue_ref != orig_ref) {
          trig <- sprintf(
            "2k1b-rule denovo %s contig %gbp pident=%g aln=%gbp kmer_cov=%g (replaced %s)",
            rescue_ref, g2_row$best_contig_length[1], g2_row$best_contig_pident[1],
            g2_row$best_contig_aln_length[1], g2_row$best_contig_kmer_cov[1], orig_ref)
          return(list(rescued_from = orig_ref, rescue_ref = rescue_ref,
                      rescue_trigger = trig))
        }
      }
    }
  }

  # ---- Normal mismatch + four-floor rescue ----------------------------------
  if (is.na(denovo_ref) || is.na(denovo_sub) || is.na(cand_sub)) return(no_rescue)
  if (denovo_sub == cand_sub) return(no_rescue)  # no mismatch → no rescue

  # 1a/1b boundary: stricter length floor when both subtypes are 1a/1b (D-04).
  length_floor <- rescue_min_length
  if (cand_sub %in% c("1a", "1b") && denovo_sub %in% c("1a", "1b")) {
    length_floor <- rescue_1a1b_length
  }

  srow <- support %>% filter(sample == sample_id, subtype == denovo_sub)
  if (nrow(srow) == 0) return(no_rescue)
  srow <- srow %>% arrange(desc(best_contig_length)) %>% slice(1)

  if (!floors_ok(srow, length_floor)) return(no_rescue)

  trig <- sprintf(
    "denovo %s contig %gbp pident=%g aln=%gbp kmer_cov=%g (replaced %s)",
    denovo_ref, srow$best_contig_length[1], srow$best_contig_pident[1],
    srow$best_contig_aln_length[1], srow$best_contig_kmer_cov[1], orig_ref)
  list(rescued_from = orig_ref, rescue_ref = denovo_ref, rescue_trigger = trig)
}

# --- Drive evaluation over every candidate row --------------------------------
out <- candidates %>%
  mutate(rescued_from = NA_character_, rescue_trigger = NA_character_)

if (nrow(out) > 0) {
  for (i in seq_len(nrow(out))) {
    res <- evaluate_row(out[i, ])
    if (!is.na(res$rescued_from) && !is.na(res$rescue_ref) &&
        res$rescue_ref %in% names(ref_fa)) {
      # Rescue fires AND the rescue ref exists in the references FASTA (guard).
      out$rescued_from[i]        <- res$rescued_from
      out$candidate_ref[i]       <- res$rescue_ref
      out$candidate_subtype[i]   <- subtype_of(res$rescue_ref)
      out$candidate_genotype[i]  <- genotype_of(subtype_of(res$rescue_ref))
      out$rescue_trigger[i]      <- res$rescue_trigger
      out$confirmation_status[i] <- "pass"  # D-05 force pass
      write_ref_fasta(res$rescue_ref, paste0("cand", out$candidate_rank[i]))
    }
    # else: membership guard failed or no rescue → leave NA, status unchanged.
  }
}

# --- Write the corrected candidates CSV --------------------------------------
# Preserve the 8 original columns + the 2 rescue audit columns.
out_cols <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
              "candidate_genotype", "candidate_reads", "candidate_cov",
              "confirmation_status", "rescued_from", "rescue_trigger")
out <- out %>% select(all_of(out_cols))
# Output name MUST differ from the input candidates CSV ({prefix}.candidates.csv,
# the PARSEFIRSTMAPPING emit staged as our input): Nextflow excludes input-named
# files from output matching, so an identically-named output is reported MISSING.
# Use a distinct {prefix}.rescued.candidates.csv that still matches the downstream
# `\\.candidates.csv$` glob in summarize.R and the module's *.candidates.csv emit.
write_csv(out, paste0(prefix, ".rescued.candidates.csv"))
