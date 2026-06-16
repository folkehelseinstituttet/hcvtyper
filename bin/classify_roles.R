#!/usr/bin/env Rscript

# classify_roles.R --------------------------------------------------------
# Phase-8 (SCORE-01/02, CLASS-01..04, COMPAT-04) dominance scoring + strain-role
# classification. This hosts the pure combination logic that turns the Phase-6
# neutral candidate set + Phase-7 genotype-level assembly-support join + the
# per-candidate targeted-mapping coverage (cv_evenness, breadth) into:
#   - a per-candidate numeric `dominance_score`   (score_candidates(); SCORE-01/02)
#   - a per-candidate `role` + `role_reason`       (classify_roles(); CLASS-01..03)
#   - one `overall_sample_call` per sample         (classify_roles(); CLASS-04)
# plus the verbatim-recovered HCV exception predicate is_valid_minor() (D-12).
#
# This is a SOURCED helper, not an entrypoint: it defines functions ONLY and has
# no commandArgs parser and no top-level file I/O and no global mutation (the
# bin/denovo_confirm.R / bin/assembly_support_join.R house pattern). Consumers
# source() it from their task workdir (the file is staged there as a declared
# `path` process input by the calling Nextflow module), AFTER first sourcing
# genotype_utils.R so that genotype_from_subtype() is already in scope. This
# helper does NOT re-source genotype_utils.R. The summarize.R wiring (arg parse,
# cov-loop cv_evenness, review_flag rewire, candidate CSV emit) is Plan 02 — this
# file is the standalone, unit-testable decision core with NO pipeline run.
#
# Decision spec (08-CONTEXT D-01..D-16; calibration anchor: handoff §2 evidence
# table). Key behaviours:
#   D-01  dominant = highest dominance_score among candidates that ALSO pass the
#         minRead/minCov major-gate; none pass => no dominant, call="indeterminate".
#   D-02  score = weighted sum of log10(reads), breadth fraction (0-1), the
#         CV-evenness factor (0-1) and a bonus-only log10(1+kmer_cov) term, with
#         breadth-evenness dominating raw read count (SCORE-02). Breadth and the
#         CV-evenness factor share the breadth-evenness weight (score_weight_evenness);
#         reads carry score_weight_reads; the k-mer-cov bonus carries score_weight_kmercov.
#   D-05  k-mer-cov term is bonus-ONLY and capped: NA / "none" support => no boost,
#         never a penalty.
#   D-06  deterministic tie-break: score, then reads, then ref name.
#   D-10  corroboration verdict = ANDed denovo floors (length/kmer_cov/pident) on
#         the candidate's joined assembly_support_* metrics at match_level.
#   D-11  asymmetric refute: cleared-floor candidate with no own support is
#         `background`/refuted_denovo ONLY if the DOMINANT assembled substantially;
#         if the dominant also failed de novo, keep it co-infection/uncorroborated_kept.
#   D-12  HCV exceptions (verbatim is_valid_minor): same-genotype (non-1a/1b) or a
#         2k1b pair => DEMOTE to background (exceptions never promote).
#   D-14  overall call: >=1 co-infection role => "co-infection"; 1 dominant + only
#         background => "monoinfection"; no candidate passes the gate => "indeterminate".
#   T-08-01 / CLASS-03  zero-row / NULL input => typed zero-row frame, never stop().
#
# The defaults 500 / 2.0 / 90 / "genotype" are calibration-VALIDATED (03-RESEARCH);
# the runtime ext.args (reconciled to 500/2.0/90 in nextflow.config, Plan 01 Task 1)
# overrides them.
# -------------------------------------------------------------------------

# Defensive: consumers already load tidyverse, which provides the pipe,
# group_by()/slice_max()/mutate()/arrange()/if_else(). This guard only fires if
# sourced into a session that has not.
if (!exists("group_by")) {
  library(tidyverse)
}

# is_valid_minor() — D-12, recovered VERBATIM from git 43904de~1:
# bin/summarize_mapping_to_all_references.R (removed in 43904de "feat(06-01)").
# Reconstructed as a PURE function taking explicit candidate/dominant subtype +
# genotype args instead of the old major_* closure vars. The three rules are
# preserved EXACTLY (COMPAT-04):
#   (1) allow 1a/1b cross-subtype as co-infection;
#   (2) block 2k1b paired with genotype {1, 2, 2k1b};
#   (3) otherwise require a different genotype.
# Returns TRUE when the candidate is a VALID minor (i.e. allowed as co-infection),
# FALSE when it must be demoted to background. Genotype comparisons use the
# already-in-scope genotype_from_subtype() — no hand-rolled substr.
is_valid_minor <- function(cand_subtype, cand_genotype, dom_subtype, dom_genotype) {
  major_subtype  <- dom_subtype
  major_genotype <- dom_genotype
  minor_subtype  <- cand_subtype
  minor_genotype <- cand_genotype

  # Rule: allow 1a and 1b co-infection
  if ((major_subtype %in% c("1a", "1b")) & (minor_subtype %in% c("1a", "1b")) & (major_subtype != minor_subtype)) {
    return(TRUE)
  }

  # Rule: block 2k1b co-infections with any genotype 1 or 2 (and itself)
  if ((major_genotype == "2k1b" & minor_genotype %in% c("1", "2", "2k1b")) |
      (minor_genotype == "2k1b" & major_genotype %in% c("1", "2", "2k1b"))) {
    return(FALSE)
  }

  # Rule: allow only different genotypes
  return(major_genotype != minor_genotype)
}

# Default dominance-score weights. Mirror the nextflow.config defaults
# (Plan 01 Task 1): breadth-evenness DOMINATES raw read count (SCORE-02). Consumers
# pass the runtime-configured weights; these defaults keep the helper self-contained
# for unit tests.
.default_score_weights <- function() {
  list(evenness = 3.0, reads = 1.0, kmercov = 0.5)
}

# score_candidates(df, score_weights, evenness_const, kmercov_cap)
#   df             : candidate frame. Expected columns (NA-tolerant):
#                    candidate_reads (numeric), the per-candidate breadth fraction
#                    (cand_cov_breadth as a 0-100 percent OR candidate_cov; coerced
#                    to a 0-1 fraction), cv_evenness (0-1 factor, supplied by the
#                    Plan-02 cov loop), and assembly_support_best_contig_kmer_cov
#                    (numeric; NA/"none" => no boost).
#   score_weights  : list(evenness=, reads=, kmercov=). evenness weights BOTH the
#                    breadth fraction and the cv_evenness factor (the breadth-evenness
#                    headline). Defaults to .default_score_weights().
#   evenness_const : transform constant for any caller that supplies a raw CV instead
#                    of a precomputed factor (factor = 1/(1 + k*CV)); unused when
#                    cv_evenness is already a 0-1 factor (the production path). Kept
#                    in the signature so the param is plumbed end-to-end (D-04).
#   kmercov_cap    : cap on the k-mer-cov bonus (D-05).
# Returns df with an added numeric `dominance_score`. Pure; no file I/O.
score_candidates <- function(df, score_weights = .default_score_weights(),
                             evenness_const = 1.0, kmercov_cap = 50) {
  if (is.null(df) || nrow(df) == 0) {
    out <- if (is.null(df)) tibble() else df
    return(out %>% mutate(dominance_score = double()))
  }

  we <- score_weights$evenness %||% 3.0
  wr <- score_weights$reads    %||% 1.0
  wk <- score_weights$kmercov  %||% 0.5

  # Breadth fraction (0-1). Prefer an explicit breadth column; fall back to
  # candidate_cov (the targeted-mapping coverage percent). Coerce a 0-100 percent
  # to a 0-1 fraction; an already-fractional value (<=1) is left as-is.
  breadth_src <- if ("cand_cov_breadth" %in% names(df)) {
    df$cand_cov_breadth
  } else if ("candidate_cov" %in% names(df)) {
    df$candidate_cov
  } else {
    rep(NA_real_, nrow(df))
  }
  breadth_frac <- ifelse(is.na(breadth_src), 0,
                         ifelse(breadth_src > 1, breadth_src / 100, breadth_src))
  breadth_frac <- pmax(0, pmin(1, breadth_frac))

  # CV-evenness factor (0-1). If a precomputed cv_evenness column is present, use
  # it (production path). Otherwise, if a raw cv column is present, transform it
  # via 1/(1 + k*CV). Otherwise neutral 0.
  if ("cv_evenness" %in% names(df)) {
    even_fac <- ifelse(is.na(df$cv_evenness), 0, pmax(0, pmin(1, df$cv_evenness)))
  } else if ("cv" %in% names(df)) {
    even_fac <- ifelse(is.na(df$cv), 0, 1 / (1 + evenness_const * df$cv))
  } else {
    even_fac <- rep(0, nrow(df))
  }

  reads <- df$candidate_reads
  reads_term <- ifelse(is.na(reads) | reads <= 0, 0, log10(reads))

  # k-mer-cov: bonus-ONLY, capped (D-05). NA / none => 0 boost, never a penalty.
  kmer <- if ("assembly_support_best_contig_kmer_cov" %in% names(df)) {
    df$assembly_support_best_contig_kmer_cov
  } else {
    rep(NA_real_, nrow(df))
  }
  kmer_capped <- ifelse(is.na(kmer) | kmer <= 0, 0, pmin(kmer, kmercov_cap))
  kmer_term   <- ifelse(kmer_capped > 0, log10(1 + kmer_capped), 0)

  df %>%
    mutate(
      dominance_score = wr * reads_term +
        we * breadth_frac +
        we * even_fac +
        wk * kmer_term
    )
}

# classify_roles(scored_df, minRead, minCov, denovo_*, match_level)
#   scored_df : candidate frame ALREADY carrying dominance_score (run
#               score_candidates() first) plus per-candidate identity
#               (candidate_subtype, candidate_genotype, candidate_ref), the
#               abundance columns the floor reads (candidate_reads, candidate_cov),
#               and the joined assembly_support_* metrics. One row per candidate;
#               may span multiple samples (grouped by sampleName).
#   minRead / minCov : the shared major-gate / co-infection floor (D-07/D-09).
#   denovo_min_contig_length / denovo_min_kmer_cov / denovo_min_blast_identity :
#               the ANDed corroboration floors (D-10), defaults validated 1000/2.0/90.
#   match_level : "genotype" (default) or "subtype" — for the D-12 exceptions; the
#               assembly-support join already collapsed support to this level.
# Returns scored_df + `role` (dominant/co-infection/background), `role_reason`
# (coded vocabulary), and `overall_sample_call` (monoinfection/co-infection/
# indeterminate). Pure; never stop() on empty input (T-08-01 / CLASS-03).
classify_roles <- function(scored_df, minRead, minCov,
                           denovo_min_contig_length = 500,
                           denovo_min_kmer_cov = 2.0,
                           denovo_min_blast_identity = 90,
                           match_level = "genotype") {
  # WR-04 (ported from assembly_support_join.R:50): reject any unsupported
  # match_level up front rather than silently treating a typo as genotype.
  stopifnot(match_level %in% c("genotype", "subtype"))

  # Typed zero-row / NULL guard (generalized from assembly_support_join.R:63-80):
  # return a typed frame carrying the new columns, never abort.
  if (is.null(scored_df) || nrow(scored_df) == 0) {
    out <- if (is.null(scored_df)) tibble() else scored_df
    out <- out %>%
      mutate(
        dominance_score     = if ("dominance_score" %in% names(.)) dominance_score else double(),
        role                = character(),
        role_reason         = character(),
        overall_sample_call = character()
      )
    return(out)
  }

  # Ensure a dominance_score column exists (defensive — caller normally scores first).
  if (!"dominance_score" %in% names(scored_df)) {
    scored_df <- score_candidates(scored_df)
  }

  # Per-candidate substantiality of OWN assembly support (D-10 ANDed floors on the
  # joined metrics). Missing metrics (assembly_support="none" / NA) => FALSE, never NA.
  has_len  <- if ("assembly_support_best_contig_length"  %in% names(scored_df)) scored_df$assembly_support_best_contig_length  else rep(NA_real_, nrow(scored_df))
  has_kmer <- if ("assembly_support_best_contig_kmer_cov" %in% names(scored_df)) scored_df$assembly_support_best_contig_kmer_cov else rep(NA_real_, nrow(scored_df))
  has_pid  <- if ("assembly_support_best_contig_pident"  %in% names(scored_df)) scored_df$assembly_support_best_contig_pident  else rep(NA_real_, nrow(scored_df))
  own_substantial <- !is.na(has_len) & !is.na(has_kmer) & !is.na(has_pid) &
    has_len  >= denovo_min_contig_length &
    has_kmer >= denovo_min_kmer_cov &
    has_pid  >= denovo_min_blast_identity

  # Per-candidate floor pass (D-07/D-09: same minRead/minCov as the dominant gate).
  reads <- scored_df$candidate_reads
  cov   <- if ("candidate_cov" %in% names(scored_df)) scored_df$candidate_cov else rep(NA_real_, nrow(scored_df))
  clears_floor <- !is.na(reads) & !is.na(cov) & reads > minRead & cov > minCov

  scored_df <- scored_df %>%
    mutate(
      .own_substantial = own_substantial,
      .clears_floor    = clears_floor,
      .row_order       = row_number()
    )

  # Group by sample so dominant determination + the asymmetric refute are per-sample.
  sample_key <- if ("sampleName" %in% names(scored_df)) "sampleName" else NULL

  classify_one_sample <- function(g) {
    n <- nrow(g)
    g$role <- NA_character_
    g$role_reason <- NA_character_

    # D-01: dominant = highest dominance_score among floor-passing candidates,
    # deterministic tie-break score -> reads -> ref name (D-06).
    gated <- which(g$.clears_floor)
    dom_idx <- NA_integer_
    if (length(gated) > 0) {
      ord <- order(
        -g$dominance_score[gated],
        -g$candidate_reads[gated],
        as.character(g$candidate_ref[gated])
      )
      dom_idx <- gated[ord[1]]
    }

    # "Did de novo work for the dominant?" (D-11) — the dominant's OWN substantiality.
    dom_substantial <- if (!is.na(dom_idx)) isTRUE(g$.own_substantial[dom_idx]) else FALSE
    dom_subtype  <- if (!is.na(dom_idx)) as.character(g$candidate_subtype[dom_idx])  else NA_character_
    dom_genotype <- if (!is.na(dom_idx)) as.character(g$candidate_genotype[dom_idx]) else NA_character_

    for (i in seq_len(n)) {
      if (!is.na(dom_idx) && i == dom_idx) {
        g$role[i] <- "dominant"
        g$role_reason[i] <- "dominant"
        next
      }

      # Below floor => background/below_floor (never promote).
      if (!isTRUE(g$.clears_floor[i])) {
        g$role[i] <- "background"
        g$role_reason[i] <- "below_floor"
        next
      }

      # Cleared the floor. Provisional corroboration verdict (D-10/D-11):
      if (isTRUE(g$.own_substantial[i])) {
        role <- "co-infection"
        reason <- "corroborated"
      } else if (dom_substantial) {
        # Candidate has no own support BUT de novo demonstrably worked for the
        # dominant => absence is real evidence (the 4g case). Asymmetric refute.
        role <- "background"
        reason <- "refuted_denovo"
      } else {
        # Neither the candidate nor the dominant assembled substantially => de novo
        # inconclusive; never suppress a genuine low-yield minor (D-11).
        role <- "co-infection"
        reason <- "uncorroborated_kept"
      }

      # D-12 HCV exceptions: applied to candidates that would be co-infection. They
      # may only DEMOTE to background, never promote (RESEARCH anti-pattern). When
      # there is no dominant, there is no genotype to compare against — skip.
      if (role == "co-infection" && !is.na(dom_idx)) {
        cand_subtype  <- as.character(g$candidate_subtype[i])
        cand_genotype <- as.character(g$candidate_genotype[i])
        valid <- is_valid_minor(cand_subtype, cand_genotype, dom_subtype, dom_genotype)
        if (!isTRUE(valid)) {
          role <- "background"
          # Distinguish the 2k1b pair from the generic same-genotype demotion.
          if (cand_genotype == "2k1b" || dom_genotype == "2k1b") {
            reason <- "recombinant_2k1b"
          } else {
            reason <- "same_genotype_as_dominant"
          }
        }
      }

      g$role[i] <- role
      g$role_reason[i] <- reason
    }

    # D-14 overall sample call.
    if (is.na(dom_idx)) {
      call <- "indeterminate"
    } else if (any(g$role == "co-infection")) {
      call <- "co-infection"
    } else {
      call <- "monoinfection"
    }
    g$overall_sample_call <- call
    g
  }

  if (!is.null(sample_key)) {
    out <- scored_df %>%
      group_split(.data[[sample_key]]) %>%
      map_dfr(classify_one_sample)
  } else {
    out <- classify_one_sample(scored_df)
  }

  out %>%
    arrange(.row_order) %>%
    select(-.own_substantial, -.clears_floor, -.row_order)
}
