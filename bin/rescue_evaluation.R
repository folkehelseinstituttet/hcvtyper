#!/usr/bin/env Rscript

# rescue_evaluation.R ----------------------------------------------------------
# De-novo subtype-rescue evaluator (Phase 10, denovo-subtype-rescue).
#
# Reads the long-format candidates table (one row per ranked candidate) together
# with the per-subtype assembly_support.csv (blastparse.R §4b), and decides —
# per candidate slot — whether the reference the reads were mapped to should be
# REPLACED by the de-novo-derived reference (a "rescue").
#
# Rescue logic uses subtype-level matching from assembly_support.csv rather than
# the legacy major/minor rank-slot assignment from blastparse.csv. A rescue fires
# when:
#   (a) the candidate's own subtype has no or weak assembly support (fails the
#       standard four floors), AND
#   (b) a different subtype has strong support (passes the floors).
# The rescue reference is taken from the best_ref column of assembly_support.csv
# (the closest database reference for that contig), keeping the decision fully
# within the assembly evidence without rank-to-slot assumptions.
#
# Two special cases apply on top of this:
#   * 2k1b special rule (D-03): a 2k1b candidate (or a sample where the best
#     alternative contig is 2k1b) over a genotype-2 contig meeting the floors
#     rescues to that genotype-2 reference regardless of the own-support guard.
#   * 1a/1b boundary (D-04): when both the candidate and alternative subtypes are
#     in {1a, 1b}, the contig-length floor is the stricter rescue_1a1b_length.
#
# Collapse guard: a rescue is blocked if its target reference is already the
# original reference of another candidate slot in the same sample, OR has already
# been committed as a rescue target for a previous slot. This prevents two
# distinct biological slots collapsing into the same reference and silently
# discarding a genuine strain signal.
#
# On rescue: rescued_from records the ORIGINAL candidate ref, candidate_ref is
# overwritten with the rescue reference, rescue_trigger gets a human-readable
# evidence string, confirmation_status is forced to "pass" (D-05; the downstream
# hcvtyper.nf filter drops non-pass rows), and the rescue reference is re-written
# as {prefix}.{ref}_cand{rank}.fa via the verbatim write_ref_fasta() membership
# guard (V5 / T-10-02: a ref absent from params.references writes NO FASTA).
#
# Empty / skip-assembly inputs (zero-row support) are handled by typed-empty
# tibbles so joins yield zero matches → candidates pass through with rescued_from
# / rescue_trigger NA-filled, and the script ALWAYS exits 0 (T-10-03 DoS guard).
#
# Positional args:
#   prefix candidates_csv support_csv references \
#   rescue_min_length rescue_min_pident rescue_min_aln_length \
#   rescue_min_kmer_cov rescue_1a1b_length
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(seqinr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: rescue_evaluation.R <prefix> <candidates.csv> ",
       "<assembly_support.csv> <references.fa> <rescue_min_length> ",
       "<rescue_min_pident> <rescue_min_aln_length> <rescue_min_kmer_cov> ",
       "<rescue_1a1b_length>")
}

prefix          <- args[1]
candidates_csv  <- args[2]
support_csv     <- args[3]
references      <- args[4]
rescue_min_length     <- as.numeric(args[5])
rescue_min_pident     <- as.numeric(args[6])
rescue_min_aln_length <- as.numeric(args[7])
rescue_min_kmer_cov   <- as.numeric(args[8])
rescue_1a1b_length    <- as.numeric(args[9])

# --- Subtype-token helper -----------------------------------------------------
subtype_of <- function(ref) {
  ifelse(is.na(ref), NA_character_, str_extract(ref, "^[^_]+"))
}

# --- Typed-empty guards (T-10-03 DoS guard) -----------------------------------
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
empty_support <- tibble(
  sample                 = character(0),
  subtype                = character(0),
  best_ref               = character(0),
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
support <- read_csv_guarded(
  support_csv, empty_support,
  cols(
    sample                 = col_character(),
    subtype                = col_character(),
    best_ref               = col_character(),
    best_contig_length     = col_double(),
    best_contig_pident     = col_double(),
    best_contig_aln_length = col_double(),
    best_contig_kmer_cov   = col_double(),
    .default               = col_character()
  )
)

# --- Reference FASTA re-extraction -------------------------------------------
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
floors_ok <- function(srow, length_floor) {
  if (nrow(srow) == 0) return(FALSE)
  isTRUE(
    srow$best_contig_length[1]     >  length_floor          &
    srow$best_contig_pident[1]     >  rescue_min_pident     &
    srow$best_contig_aln_length[1] >= rescue_min_aln_length &
    srow$best_contig_kmer_cov[1]   >= rescue_min_kmer_cov
  )
}

genotype_of <- function(subtype) {
  ifelse(is.na(subtype), NA_character_, str_sub(subtype, 1, 1))
}

# --- Per-candidate rescue evaluation -----------------------------------------
# Uses assembly_support.csv with subtype-level matching. For each candidate:
#   1. If the candidate's own subtype has strong assembly support (passes
#      floors), the de novo confirms the reference — no rescue (unless the
#      candidate is 2k1b, where the special rule takes priority).
#   2. Find the best-supported alternative subtype (longest contig, different
#      subtype). Apply 2k1b special rule if applicable, otherwise attempt the
#      normal four-floor rescue.
evaluate_row <- function(row) {
  sample_id <- row$sample
  orig_ref  <- row$candidate_ref
  cand_sub  <- subtype_of(orig_ref)

  no_rescue <- list(rescued_from = NA_character_, rescue_ref = NA_character_,
                    rescue_trigger = NA_character_)

  if (is.na(cand_sub)) return(no_rescue)

  sample_support <- support %>% filter(sample == sample_id)
  if (nrow(sample_support) == 0) return(no_rescue)

  # Own-subtype floor check: if de novo confirms the candidate's own subtype,
  # no rescue needed. Skipped for 2k1b candidates — the 2k1b special rule
  # (D-03) fires regardless of own-support quality.
  if (!isTRUE(cand_sub == "2k1b")) {
    own_row <- sample_support %>% filter(subtype == cand_sub) %>%
      arrange(desc(best_contig_length)) %>% slice(1)
    if (nrow(own_row) > 0 && floors_ok(own_row, rescue_min_length)) return(no_rescue)
  }

  # Find the best-supported alternative subtype (longest contig, different
  # from the candidate's own subtype). This is the assembly's primary
  # disagreement signal.
  alt_support <- sample_support %>%
    filter(subtype != cand_sub) %>%
    arrange(desc(best_contig_length)) %>%
    slice(1)

  if (nrow(alt_support) == 0) return(no_rescue)

  denovo_sub <- alt_support$subtype[1]
  denovo_ref <- alt_support$best_ref[1]

  # ---- 2k1b special rule (D-03) ---------------------------------------------
  slot_is_2k1b <- isTRUE(cand_sub == "2k1b") || isTRUE(denovo_sub == "2k1b")
  if (slot_is_2k1b) {
    g2_support <- sample_support %>%
      filter(genotype_of(subtype) == "2") %>%
      arrange(desc(best_contig_length))
    if (nrow(g2_support) > 0) {
      g2_row <- g2_support[1, ]
      if (floors_ok(g2_row, rescue_min_length)) {
        rescue_ref <- g2_row$best_ref[1]
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

  # ---- Normal four-floor rescue ---------------------------------------------
  if (is.na(denovo_ref) || is.na(denovo_sub)) return(no_rescue)

  # 1a/1b boundary: stricter length floor when both subtypes are 1a/1b (D-04).
  length_floor <- rescue_min_length
  if (cand_sub %in% c("1a", "1b") && denovo_sub %in% c("1a", "1b")) {
    length_floor <- rescue_1a1b_length
  }

  if (!floors_ok(alt_support, length_floor)) return(no_rescue)

  trig <- sprintf(
    "denovo %s contig %gbp pident=%g aln=%gbp kmer_cov=%g (replaced %s)",
    denovo_ref, alt_support$best_contig_length[1], alt_support$best_contig_pident[1],
    alt_support$best_contig_aln_length[1], alt_support$best_contig_kmer_cov[1], orig_ref)
  list(rescued_from = orig_ref, rescue_ref = denovo_ref, rescue_trigger = trig)
}

# --- Drive evaluation over every candidate row --------------------------------
out <- candidates %>%
  mutate(rescued_from = NA_character_, rescue_trigger = NA_character_)

if (nrow(out) > 0) {
  # Snapshot original candidate refs before any rescue mutations.
  # The collapse guard checks against both the original refs (ordering-
  # independent detection of original slot conflicts) and the current
  # post-rescue state (preventing two candidates being rescued to the same
  # new reference within the same sample).
  orig_refs_snap <- candidates %>% select(sample, candidate_rank, candidate_ref)

  for (i in seq_len(nrow(out))) {
    res <- evaluate_row(out[i, ])

    # Collapse guard: block rescue if rescue_ref is already held by another
    # candidate slot — either as its original ref or as an already-committed
    # rescue target in this run.
    rescue_would_collapse <- if (!is.na(res$rescue_ref)) {
      other_orig <- orig_refs_snap %>%
        filter(sample == out$sample[i], candidate_rank != out$candidate_rank[i]) %>%
        pull(candidate_ref)
      current_other <- out %>%
        filter(sample == out$sample[i], candidate_rank != out$candidate_rank[i]) %>%
        pull(candidate_ref)
      res$rescue_ref %in% union(other_orig, current_other)
    } else {
      FALSE
    }

    if (!is.na(res$rescued_from) && !is.na(res$rescue_ref) &&
        res$rescue_ref %in% names(ref_fa) &&
        !rescue_would_collapse) {
      out$rescued_from[i]        <- res$rescued_from
      out$candidate_ref[i]       <- res$rescue_ref
      out$candidate_subtype[i]   <- subtype_of(res$rescue_ref)
      out$candidate_genotype[i]  <- genotype_of(subtype_of(res$rescue_ref))
      out$rescue_trigger[i]      <- res$rescue_trigger
      out$confirmation_status[i] <- "pass"  # D-05 force pass
      write_ref_fasta(res$rescue_ref, paste0("cand", out$candidate_rank[i]))
      # Remove the now-stale pass-through FASTA for the REPLACED ref at this rank.
      stale_fa <- paste0(prefix, ".", res$rescued_from, "_cand",
                         out$candidate_rank[i], ".fa")
      if (!identical(res$rescued_from, res$rescue_ref) && file.exists(stale_fa)) {
        file.remove(stale_fa)
      }
    }
    # else: membership guard failed, collapse guard blocked, or no rescue → leave NA, status unchanged.
  }
}

# --- Write the corrected candidates CSV --------------------------------------
out_cols <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
              "candidate_genotype", "candidate_reads", "candidate_cov",
              "confirmation_status", "rescued_from", "rescue_trigger")
out <- out %>% select(all_of(out_cols))
write_csv(out, paste0(prefix, ".rescued.candidates.csv"))
