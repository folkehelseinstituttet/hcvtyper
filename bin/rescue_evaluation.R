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
# Same-genotype guard: a rescue is blocked if the rescue target's genotype
# already matches another candidate slot's genotype in the same sample.
# Exception: 1a and 1b are treated as a permitted cross-subtype co-infection pair
# (matching classify_roles.R is_valid_minor() logic), so a 1a+1b pair is allowed.
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
# Optional 10th positional arg: the genotype-diverse candidate cap (params.n_candidates).
# Bounds the final candidate set to at most this many DISTINCT-genotype references
# after de-novo nomination + genotype collapse. Absent / non-integer => no cap (Inf),
# which preserves the legacy behaviour for callers (and unit tests) that omit it.
n_candidates_cap <- if (length(args) >= 10) suppressWarnings(as.integer(args[10])) else NA_integer_
if (is.na(n_candidates_cap) || n_candidates_cap < 1L) n_candidates_cap <- Inf

# Optional 11th/12th positional args (Fix #2/#3). Absent => guard disabled, so
# the 9-arg subprocess test and any legacy caller preserve prior behaviour.
rescue_kmer_cov_ratio <- if (length(args) >= 11) suppressWarnings(as.numeric(args[11])) else NA_real_
if (is.na(rescue_kmer_cov_ratio) || rescue_kmer_cov_ratio <= 0) rescue_kmer_cov_ratio <- Inf
dominant_protect_cov  <- if (length(args) >= 12) suppressWarnings(as.numeric(args[12])) else NA_real_
if (is.na(dominant_protect_cov)) dominant_protect_cov <- Inf

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

  # Best own-subtype contig depth — reference point for Fix #2 / #3.
  own_rows <- sample_support %>% filter(subtype == cand_sub)
  own_kmer <- if (nrow(own_rows) > 0) max(own_rows$best_contig_kmer_cov, na.rm = TRUE) else NA_real_

  # Fix #3 — relative k-mer-coverage guard (used at both rescue-return points).
  # A REPLACE target whose k-mer depth is dwarfed by the candidate's own assembled
  # subtype is cross-mapping noise sitting on a real strain; refuse. Inert when the
  # own subtype has no contig (own_kmer NA) — preserves subtest 1.
  kmer_ratio_blocks <- function(target_kmer) {
    !is.na(own_kmer) && is.finite(own_kmer) && !is.na(target_kmer) &&
      target_kmer > 0 && (own_kmer / target_kmer) >= rescue_kmer_cov_ratio
  }

  # Own-subtype confirmation. Skipped for 2k1b (D-03 fires regardless).
  if (!isTRUE(cand_sub == "2k1b") && nrow(own_rows) > 0) {
    # Standard: longest own-subtype contig passes all four floors.
    own_longest <- own_rows %>% arrange(desc(best_contig_length)) %>% slice(1)
    if (floors_ok(own_longest, rescue_min_length)) return(no_rescue)
    # Fix #2 — mapping-aware confirmation. A candidate already well covered by
    # first-mapping reads is a real dominant strain; its de-novo contig merely
    # assembling SHORT must not make it eligible for replacement. Confirm if the
    # own subtype clears the QUALITY floors (identity + k-mer depth) even when it
    # fails the LENGTH / ALN floors. (Subtest 1 stays a replace: there the own
    # subtype has NO contig, so this cannot fire.)
    own_best_q <- own_rows %>% arrange(desc(best_contig_kmer_cov)) %>% slice(1)
    own_quality_ok <- isTRUE(
      own_best_q$best_contig_pident[1]   >  rescue_min_pident &
      own_best_q$best_contig_kmer_cov[1] >= rescue_min_kmer_cov
    )
    if (isTRUE(row$candidate_cov >= dominant_protect_cov) && own_quality_ok)
      return(no_rescue)
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
        if (!is.na(rescue_ref) && rescue_ref != orig_ref &&
            !kmer_ratio_blocks(g2_row$best_contig_kmer_cov[1])) {
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
  if (kmer_ratio_blocks(alt_support$best_contig_kmer_cov[1])) return(no_rescue)  # Fix #3

  trig <- sprintf(
    "denovo %s contig %gbp pident=%g aln=%gbp kmer_cov=%g (replaced %s)",
    denovo_ref, alt_support$best_contig_length[1], alt_support$best_contig_pident[1],
    alt_support$best_contig_aln_length[1], alt_support$best_contig_kmer_cov[1], orig_ref)
  list(rescued_from = orig_ref, rescue_ref = denovo_ref, rescue_trigger = trig)
}

# --- Drive evaluation over every candidate row --------------------------------
out <- candidates %>%
  mutate(rescued_from = NA_character_, rescue_trigger = NA_character_)

# --- Rescue audit ledger (Fix #4) --------------------------------------------
# Records EVERY rescue/nomination/block decision so a fired-then-dropped rescue
# (collapsed or capped) remains traceable in the published output.
audit_rows <- list()
add_audit <- function(sample, event, orig_ref, target_ref, evidence) {
  audit_rows[[length(audit_rows) + 1]] <<- tibble(
    sample = sample, event = event,
    original_ref = orig_ref, original_subtype = subtype_of(orig_ref),
    target_ref = target_ref, target_subtype = subtype_of(target_ref),
    target_genotype = genotype_of(subtype_of(target_ref)),
    evidence = evidence)
}

# Pre-cap survivor set (Fix #4). Declared in outer scope so it exists even when
# there are no candidates (empty / skip-assembly inputs).
kept_refs_precap <- character(0)

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

    rescue_would_dup_genotype <- if (!is.na(res$rescue_ref)) {
      rescue_sub  <- subtype_of(res$rescue_ref)
      rescue_geno <- genotype_of(rescue_sub)
      other_refs <- out %>%
        filter(sample == out$sample[i], candidate_rank != out$candidate_rank[i]) %>%
        pull(candidate_ref)
      any(vapply(other_refs, function(r) {
        other_sub  <- subtype_of(r)
        other_geno <- genotype_of(other_sub)
        if (rescue_geno != other_geno) return(FALSE)
        !( rescue_sub %in% c("1a","1b") &&
           other_sub  %in% c("1a","1b") &&
           rescue_sub != other_sub )
      }, logical(1)))
    } else {
      FALSE
    }

    if (!is.na(res$rescued_from) && !is.na(res$rescue_ref) &&
        res$rescue_ref %in% names(ref_fa) &&
        !rescue_would_collapse &&
        !rescue_would_dup_genotype) {
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
      add_audit(out$sample[i], "replace", res$rescued_from, res$rescue_ref, res$rescue_trigger)
    } else if (!is.na(res$rescue_ref) && rescue_would_collapse) {
      add_audit(out$sample[i], "blocked_collapse", out$candidate_ref[i], res$rescue_ref, res$rescue_trigger)
    } else if (!is.na(res$rescue_ref) && rescue_would_dup_genotype) {
      add_audit(out$sample[i], "blocked_dup_genotype", out$candidate_ref[i], res$rescue_ref, res$rescue_trigger)
    }
    # else: membership guard failed or no rescue → leave NA, status unchanged.
  }
}

# =============================================================================
# Genotype-diverse candidate finalization -------------------------------------
# The pipeline CANNOT resolve within-genotype (same-genotype) co-infections, so
# at most ONE reference per genotype may enter JOINT_MAPPING — the sole exception
# being the 1a/1b cross-subtype pair (mirrors classify_roles::is_valid_minor()).
# Two steps run AFTER the per-slot rescue above:
#   (1) De-novo nomination (additive): a genuine DIFFERENT-genotype strain that
#       first-mapping under-ranked but de novo assembled strongly (e.g. the
#       ERR1810475 1a) is APPENDED as a new candidate from assembly_support,
#       gated on the same four rescue floors. First-mapping ranks only the top-N
#       subtypes by read recruitment, so a low-recruitment second genotype that
#       only survives as a de-novo contig is otherwise unmappable; this is the
#       only path that surfaces it in the mapping minor call.
#   (2) Genotype collapse: the candidate set is reduced to one reference per
#       genotype (best first-mapping read recruiter kept), 1a/1b kept as a pair,
#       the 2k1b recombinant blocked vs genotype 1/2. This removes redundant
#       same-genotype slots (the cosmetic same-genotype "Minor", and the divergent
#       genotype-4 second reference) at SELECTION rather than only labelling them
#       downstream — a sub-optimal-but-correct-genotype reference still yields the
#       right consensus under Bowtie2 mismatch tolerance.
# Candidates are then re-ranked 1..N by read recruitment (dominant first; nominated
# NA-read candidates last), capped at n_candidates_cap, and the per-rank candidate
# FASTAs are reconciled to the final set (stale ones removed, survivors written) so
# exactly one *_cand{rank}.fa exists per surviving candidate (dup-@SQ guard).
# -----------------------------------------------------------------------------

# Local genotype + validity helpers. genotype_utils.R is NOT staged into this
# module, so replicate the minimal 2k1b-aware logic of genotype_from_subtype()
# and the verbatim three-rule is_valid_minor() (classify_roles.R:137).
geno_of_subtype <- function(subtype) {
  ifelse(is.na(subtype), NA_character_,
         ifelse(str_starts(subtype, "2k1b"), "2k1b", str_sub(subtype, 1, 1)))
}
valid_minor_pair <- function(minor_sub, minor_geno, major_sub, major_geno) {
  if (is.na(minor_sub) || is.na(major_sub)) return(FALSE)
  # Rule 1: allow 1a/1b cross-subtype co-infection.
  if (minor_sub %in% c("1a", "1b") && major_sub %in% c("1a", "1b") &&
      minor_sub != major_sub) return(TRUE)
  # Rule 2: block 2k1b paired with genotype 1 / 2 / 2k1b.
  if ((major_geno == "2k1b" && minor_geno %in% c("1", "2", "2k1b")) ||
      (minor_geno == "2k1b" && major_geno %in% c("1", "2", "2k1b"))) return(FALSE)
  # Rule 3: otherwise require a different genotype.
  major_geno != minor_geno
}

if (nrow(out) > 0) {
  sample_id <- out$sample[1]

  # ---- (1) De-novo additive nomination --------------------------------------
  present_genos <- unique(out$candidate_genotype)
  # Dominant (for the validity check) = highest first-mapping read recruiter.
  dom_row  <- out %>% arrange(desc(candidate_reads)) %>% slice(1)
  dom_sub  <- dom_row$candidate_subtype[1]
  dom_geno <- dom_row$candidate_genotype[1]

  nominable <- support %>%
    filter(!is.na(best_ref), best_ref %in% names(ref_fa)) %>%
    mutate(.geno = geno_of_subtype(subtype)) %>%
    filter(
      !is.na(.geno), !(.geno %in% present_genos),
      best_contig_length     >  rescue_min_length,
      best_contig_pident     >  rescue_min_pident,
      best_contig_aln_length >= rescue_min_aln_length,
      best_contig_kmer_cov   >= rescue_min_kmer_cov
    ) %>%
    group_by(.geno) %>%
    arrange(desc(best_contig_length), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

  if (nrow(nominable) > 0) {
    next_rank <- max(out$candidate_rank, na.rm = TRUE)
    for (j in seq_len(nrow(nominable))) {
      nsub  <- nominable$subtype[j]
      ngeno <- nominable$.geno[j]
      nref  <- nominable$best_ref[j]
      # Validity vs the dominant (e.g. block a 2k1b nomination over a genotype-2
      # dominant). The collapse below would drop an invalid pair anyway; this just
      # avoids creating then deleting its FASTA.
      if (!valid_minor_pair(nsub, ngeno, dom_sub, dom_geno)) next
      next_rank <- next_rank + 1L
      trig <- sprintf(
        "denovo-nomination %s contig %gbp pident=%g aln=%gbp kmer_cov=%g (different-genotype strain absent from first-mapping candidates)",
        nref, nominable$best_contig_length[j], nominable$best_contig_pident[j],
        nominable$best_contig_aln_length[j], nominable$best_contig_kmer_cov[j])
      out <- bind_rows(out, tibble(
        sample              = sample_id,
        candidate_rank      = as.integer(next_rank),
        candidate_ref       = nref,
        candidate_subtype   = nsub,
        candidate_genotype  = ngeno,
        candidate_reads     = NA_real_,
        candidate_cov       = NA_real_,
        confirmation_status = "pass",
        rescued_from        = NA_character_,
        rescue_trigger      = trig
      ))
      add_audit(sample_id, "nominate", NA_character_, nref, trig)
    }
  }

  # ---- (2) Genotype collapse to distinct genotypes (1a/1b pair exempt) -------
  ord <- out %>% arrange(desc(candidate_reads))
  dom_sub2  <- ord$candidate_subtype[1]
  dom_geno2 <- ord$candidate_genotype[1]
  keep_idx   <- integer(0)
  kept_subs  <- character(0)
  kept_genos <- character(0)
  for (i in seq_len(nrow(ord))) {
    csub  <- ord$candidate_subtype[i]
    cgeno <- ord$candidate_genotype[i]
    if (length(keep_idx) == 0) {            # the dominant always survives
      keep_idx <- i; kept_subs <- csub; kept_genos <- cgeno; next
    }
    # Must be a valid minor w.r.t. the dominant ...
    if (!valid_minor_pair(csub, cgeno, dom_sub2, dom_geno2)) next
    # ... and a distinct genotype from EVERY already-kept candidate (1a/1b pair ok).
    dup <- any(vapply(seq_along(kept_subs), function(k) {
      if (cgeno != kept_genos[k]) return(FALSE)
      !(csub %in% c("1a", "1b") && kept_subs[k] %in% c("1a", "1b") &&
        csub != kept_subs[k])
    }, logical(1)))
    if (dup) next
    keep_idx   <- c(keep_idx, i)
    kept_subs  <- c(kept_subs, csub)
    kept_genos <- c(kept_genos, cgeno)
  }

  # ---- Re-rank 1..N (dominant first) and cap at n_candidates_cap -------------
  # Capture the pre-cap survivor set so the audit can distinguish a rescue that
  # was dropped by genotype-collapse from one dropped by the n_candidates cap.
  kept_refs_precap <- ord[keep_idx, , drop = FALSE]$candidate_ref
  out <- ord[keep_idx, , drop = FALSE] %>%
    arrange(desc(candidate_reads)) %>%
    slice(seq_len(min(n(), n_candidates_cap))) %>%
    mutate(candidate_rank = row_number())

  # ---- Reconcile per-rank candidate FASTAs to the final set -----------------
  # Exactly one {prefix}.{ref}_cand{rank}.fa must survive per final candidate:
  # remove any stale / collapsed-away / wrong-rank per-rank FASTA, then (re)write
  # the survivors. A duplicated or stale per-rank FASTA would duplicate an @SQ line
  # in the combined per-sample reference and crash BOWTIE2_BUILD.
  final_fastas <- paste0(prefix, ".", out$candidate_ref, "_cand", out$candidate_rank, ".fa")
  existing_fastas <- list.files(".", pattern = "_cand[0-9]+\\.fa$")
  existing_fastas <- existing_fastas[startsWith(existing_fastas, paste0(prefix, "."))]
  for (f in setdiff(existing_fastas, final_fastas)) {
    if (file.exists(f)) file.remove(f)
  }
  for (i in seq_len(nrow(out))) {
    write_ref_fasta(out$candidate_ref[i], paste0("cand", out$candidate_rank[i]))
  }
}

# --- Stale first-mapping-stat guard ------------------------------------------
# A REPLACED candidate's first-mapping read count / coverage described the
# DISPLACED reference, not the reference now shown in candidate_ref. Blank them
# (as the additive-nomination path already does) so no downstream table reports
# the old ref's reads / first-mapping % against the new ref. Targeted mapping
# recomputes the real per-candidate numbers.
out <- out %>%
  mutate(
    candidate_reads = if_else(!is.na(rescued_from), NA_real_, candidate_reads),
    candidate_cov   = if_else(!is.na(rescued_from), NA_real_, candidate_cov)
  )

# --- Write the corrected candidates CSV --------------------------------------
out_cols <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
              "candidate_genotype", "candidate_reads", "candidate_cov",
              "confirmation_status", "rescued_from", "rescue_trigger")
out <- out %>% select(all_of(out_cols))

# Resolve each audit row's fate against the FINAL candidate set.
audit <- if (length(audit_rows) > 0) bind_rows(audit_rows) else tibble(
  sample = character(), event = character(),
  original_ref = character(), original_subtype = character(),
  target_ref = character(), target_subtype = character(),
  target_genotype = character(), evidence = character())
final_rank <- out %>% select(target_ref = candidate_ref, .rank = candidate_rank)
audit <- audit %>%
  left_join(final_rank, by = "target_ref") %>%
  mutate(disposition = case_when(
    event %in% c("blocked_collapse", "blocked_dup_genotype") ~ "blocked",
    !is.na(.rank)                                            ~ paste0("retained_rank_", .rank),
    target_ref %in% kept_refs_precap                         ~ "dropped_cap",
    TRUE                                                     ~ "dropped_collapse"
  )) %>%
  select(-.rank)
write_csv(audit, paste0(prefix, ".rescue_audit.csv"))

write_csv(out, paste0(prefix, ".rescued.candidates.csv"))
