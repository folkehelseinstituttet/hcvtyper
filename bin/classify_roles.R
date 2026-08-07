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
#   D-10/D-11 (WR-01, 12-REVIEW: superseded by Phase-12 Plan-03 — see EVID-01..04
#         below; kept here for the ID numbering, not for the OLD behaviour they used
#         to describe). Role / role_reason are now DERIVED from a per-candidate,
#         continuous `evidence_state` (confirmed/probable/weak; `refuted` removed
#         260805 as structurally unreachable) computed by
#         score_assembly_support() + classify_roles()'s evidence-band derivation
#         (see the D-01/D-09/D-10/D-05 comment above evidence_hi_cut/evidence_lo_cut,
#         and the EVID-02/EVID-04 comments in classify_one_sample()). This REPLACED
#         the old binary ANDed-floor `own_substantial` corroboration check and its
#         dominance-DEPENDENT asymmetric refute (a no-own-support candidate used to
#         be `background`/refuted_denovo only if the dominant ALSO assembled
#         substantially, else kept as `co-infection`/uncorroborated_kept). Under the
#         new model a candidate's evidence_state — and hence its role — is identical
#         regardless of which candidate wins dominance (EVID-02); `uncorroborated_kept`
#         no longer exists in the role_reason vocabulary at all.
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

# has_own_denovo(assembly_support, assembly_support_subtype)
#   Does this candidate have a de novo contig attributed to it? Presence only.
#
#   260805 (refuted-unreachable): this replaces own_denovo_conflict(), which also
#   returned a `conflict` flag — "the candidate's OWN de novo genotype differs from
#   its OWN mapping genotype". That predicate was STRUCTURALLY UNSATISFIABLE and has
#   been removed along with its two consumers.
#
#   Why it could never fire: assembly_support_join() attributes a contig to a
#   candidate only when the contig's genotype EQUALS the candidate's — the join key
#   is `candidate_genotype` on one side and `genotype_from_subtype(subtype)` on the
#   other (assembly_support_join.R:90-131), and assembly_support_subtype is then the
#   winning contig's subtype from within that key group. So
#   genotype_from_subtype(assembly_support_subtype) == candidate_genotype by
#   construction, and a predicate testing those two for INEQUALITY is always FALSE.
#   Confirmed on 358 candidate rows across 10 result directories: zero rows where the
#   attributed contig's genotype differs from the candidate's.
#
#   WR-04 had unified the two copies of that predicate so they could not drift apart.
#   They could not — they were identically dead.
#
#   Genuine mapping-vs-assembly contradiction is owned by RESCUE_EVALUATION
#   (rescue_evaluation.R), which filters the support table by SAMPLE only and then
#   explicitly seeks the best contig of a DIFFERENT subtype (L207, L246-249). Because
#   it is not keyed on genotype equality it can see what this function never could.
#
#   Pure; vectorised. NA/blank subtype counts as absent.
has_own_denovo <- function(assembly_support, assembly_support_subtype, n = NULL) {
  if (is.null(n)) n <- max(length(assembly_support), length(assembly_support_subtype))
  asup_subtype <- if (is.null(assembly_support_subtype)) rep(NA_character_, n) else as.character(assembly_support_subtype)
  supported <- if (is.null(assembly_support)) rep(FALSE, n) else (!is.na(assembly_support) & assembly_support == "supported")
  supported & !is.na(asup_subtype) & nzchar(asup_subtype)
}

# apply_concordance(df) — D8 pre-annotation helper.
# Compares three identity legs at genotype level and annotates each candidate with:
#   concordance_status : "confirmed" (all legs present and agree),
#                        "unconfirmed" (no conflict, but a leg is absent),
#                        "discordant" (mapping vs GLUE or mapping vs de novo conflict)
#   concordance_reason : short coded reason string
# Inputs required in df:
#   candidate_genotype          : mapping-derived genotype (character, always present)
#   candidate_glue_genotype     : GLUE-derived genotype (character, NA when GLUE absent)
#   assembly_support            : "supported" or "none" (from assembly_support_join)
#   assembly_support_subtype    : subtype of best de novo contig (character, NA when none)
# Pure — no file I/O, no side effects. Safe on NULL/zero-row input.
apply_concordance <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    out <- if (is.null(df)) tibble() else df
    return(out %>% mutate(concordance_status = character(), concordance_reason = character()))
  }

  # Ensure expected columns exist with NA defaults when absent (robustness for unit tests
  # that may not supply all three legs).
  if (!"candidate_glue_genotype"  %in% names(df)) df$candidate_glue_genotype  <- NA_character_
  if (!"assembly_support"         %in% names(df)) df$assembly_support         <- "none"
  if (!"assembly_support_subtype" %in% names(df)) df$assembly_support_subtype <- NA_character_

  has_glue <- !is.na(df$candidate_glue_genotype) & nzchar(as.character(df$candidate_glue_genotype))

  # 260805: the de novo leg is PRESENCE-ONLY. It contributes to confirmed vs
  # unconfirmed (a de novo contig corroborates the mapping identity) but can no
  # longer contribute to `discordant`, because an attributed contig always shares the
  # candidate's genotype — see has_own_denovo()'s docstring. The GLUE leg below is
  # unaffected and remains the live source of `discordant`.
  has_denovo <- has_own_denovo(df$assembly_support, df$assembly_support_subtype, nrow(df))

  map_gt    <- as.character(df$candidate_genotype)
  glue_gt   <- as.character(df$candidate_glue_genotype)

  # 2k1b structural exception (CLAUDE.md Constraints): HCV-GLUE's clade-placement
  # tree has no CRF_02k/1b category, so a genuine 2k/1b recombinant is ALWAYS
  # reported by GLUE as genotype 1 or 2 (whichever region/majority-length portion
  # dominates the consensus), never "2k1b". This is expected GLUE behaviour, not
  # evidence of a wrong mapping/de novo call. Mirrors the 2k1b-aware exception
  # already applied to co-infection pairing in is_valid_minor() (below) and to
  # the D-03 rescue rule in rescue_evaluation.R. GLUE leg only — the de novo leg
  # already agrees natively via genotype_from_subtype()'s 2k1b-aware rule.
  glue_2k1b_exempt <- map_gt == "2k1b" & glue_gt %in% c("1", "2")

  glue_conflict   <- has_glue & !is.na(glue_gt) & map_gt != glue_gt & !glue_2k1b_exempt

  status <- character(nrow(df))
  reason <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    # 260805: `discordant` now has exactly one source, the GLUE leg. The former
    # de novo leg was structurally unsatisfiable, so the two reason strings that
    # depended on it ("discordant_all_legs", "discordant_mapping_vs_denovo") were
    # unreachable and have been removed with it. `discordant` itself and the
    # downstream discordant_identity role gate are LIVE and untouched.
    if (glue_conflict[i]) {
      status[i] <- "discordant"
      reason[i] <- "discordant_mapping_vs_glue"
    } else if (has_glue[i] && has_denovo[i]) {
      status[i] <- "confirmed"
      reason[i] <- if (glue_2k1b_exempt[i]) "confirmed_2k1b_recombinant" else "all_legs_concordant"
    } else if (has_glue[i] || has_denovo[i]) {
      status[i] <- "unconfirmed"
      reason[i] <- if (has_glue[i]) {
        if (glue_2k1b_exempt[i]) "two_legs_2k1b_recombinant_glue_only" else "two_legs_glue_only"
      } else {
        "two_legs_denovo_only"
      }
    } else {
      status[i] <- "unconfirmed"
      reason[i] <- "no_corroborating_legs"
    }
  }

  df %>% mutate(concordance_status = status, concordance_reason = reason)
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

# Default assembly-support-score weights (EVID-01 / D-11). Parallels
# .default_score_weights() in house style, but carries ONLY the identity/length/
# kmer_bonus keys — deliberately NO reads/breadth/evenness keys, so the score can
# never read candidate_reads/candidate_cov/cv_evenness (D-11: assembly_support_score
# stays purely about the candidate's OWN de novo/assembly evidence, orthogonal to
# dominance_score). Constants are the Phase-12 calibration-VALIDATED defaults
# (12-RESEARCH §Real Data Calibration; see score_assembly_support() header).
.default_assembly_weights <- function() {
  list(identity = 0.65, length = 0.35, kmer_bonus = 0.10)
}

# score_assembly_support(df, w, id_center, id_slope, len_ref, kmer_cap)
#   Adds a bounded continuous `assembly_support_score` (0-1) and an `assembly_exists`
#   boolean to df, replacing the meaning of the binary ANDed `own_substantial`
#   collapse with a smooth, calibrated evidence signal (EVID-01). The score reads
#   ONLY the candidate's own assembly metrics (D-11):
#     assembly_support_best_contig_pident   (identity %),
#     assembly_support_best_contig_length   (bp),
#     assembly_support_best_contig_kmer_cov (k-mer coverage; bonus-only).
#   It NEVER reads candidate_reads / candidate_cov / cv_evenness / targeted_reads_nodup.
#
#   Formula (bounded to [0,1], D-08):
#     id_term  = logistic((pident - id_center) / id_slope)         -- smooth, NO 90% cliff
#     len_term = min(length / len_ref, 1)                          -- saturating partial credit
#     base     = w$identity * id_term + w$length * len_term
#     kbonus   = w$kmer_bonus * log10(1 + min(kmer, kmer_cap)) / log10(1 + kmer_cap)  -- bonus-only (D-10)
#     score    = min(base + kbonus, 1)                             -- floored to 0 when no assembly (D-09)
#
#   Calibration (12-RESEARCH §Real Data Calibration, D-12/D-13): the identity
#   logistic is centered at 88 — BELOW 90 — because real dominant assemblies reach
#   81.34% identity and the genuine corroborated co-infection band clusters at
#   90.4-94.1%; centering a steep curve at 90 would re-create the very cliff EVID-01
#   removes. Named-anchor scores under these constants: Sample51K-2c (89.009%, 9479bp,
#   k-mer 514) -> 0.874; Sample61K-2c (88.987%, 9477bp, k-mer 103) -> 0.872;
#   2714372 1a (90.996%, 6811bp, k-mer 1.93 — k-mer no longer penalizes) -> 0.941;
#   2768856 4d (no assembly) -> 0.000. The 15 genuine corroborated minors all land
#   >= 0.87, cleanly separable from the score-0 no-assembly floor.
#
#   D-09: a candidate with NO own assembly (assembly_exists FALSE, or NA identity/
#   length) returns score EXACTLY 0, not NA. `assembly_exists` is carried as a
#   SEPARATE boolean (derived from assembly_support == "supported" OR non-NA
#   identity+length) so no_own_assembly stays distinguishable from a weak-but-
#   present assembly (D-06).
#
#   Pure; no file I/O. NA-tolerant; never stop() on empty/NULL input (CLASS-03/
#   T-08-01) — the typed zero-row/NULL guard carries the new columns.
#
#   WR-02 (12-REVIEW) — CALIBRATION-LOCKED BY DESIGN, not runtime-configurable:
#   unlike score_candidates()'s score_weight_* args (threaded end-to-end from
#   nextflow.config -> conf/modules_hcv.config ext.args -> summarize.R argv ->
#   score_candidates(score_weights=...)), score_assembly_support()'s w/id_center/
#   id_slope/len_ref/kmer_cap args and classify_roles()'s local evidence_hi_cut/
#   evidence_lo_cut band cutpoints (defined just above the evidence_state_col
#   ifelse() cascade) are ALWAYS called at their compiled-in defaults —
#   classify_roles() invokes score_assembly_support(scored_df) with zero extra
#   args. This is DELIBERATE, not an oversight: these six-plus numbers are
#   calibration anchors validated against the real 203-candidate dataset
#   (12-RESEARCH §Real Data Calibration) and are meant to move only via a new
#   calibration pass with fresh named-anchor evidence, not via a per-run
#   ext.args knob an operator could silently mistune. If a future need for
#   runtime tuning arises, thread them through classify_roles()'s signature and
#   summarize.R's arg list mirroring the score_weight_* pattern — but that is an
#   intentional escalation, not the current design.
score_assembly_support <- function(df, w = .default_assembly_weights(),
                                   id_center = 88, id_slope = 1.6,
                                   len_ref = 3000, kmer_cap = 50) {
  if (is.null(df) || nrow(df) == 0) {
    out <- if (is.null(df)) tibble() else df
    # EVID-05/D-03: the zero-row / NULL path declares the SAME schema as the
    # populated path so a header-only candidates.csv stays byte-stable (T-08-01).
    # Mirrors the assembly_support_join.R:52-79 declare-once lockstep — the three
    # contribution columns are carried here as length-0 doubles in step with the
    # populated mutate below.
    return(out %>% mutate(
      assembly_support_score       = double(),
      assembly_exists              = logical(),
      contig_identity_contribution = double(),
      contig_length_contribution   = double(),
      contig_kmer_contribution     = double()
    ))
  }

  wi <- w$identity   %||% 0.65
  wl <- w$length     %||% 0.35
  wk <- w$kmer_bonus %||% 0.10

  # NA-tolerant, presence-checked extraction (score_candidates() idiom). Reads ONLY
  # the three OWN-assembly metric columns — never reads/breadth/evenness (D-11).
  pid <- if ("assembly_support_best_contig_pident"   %in% names(df)) df$assembly_support_best_contig_pident   else rep(NA_real_, nrow(df))
  len <- if ("assembly_support_best_contig_length"   %in% names(df)) df$assembly_support_best_contig_length   else rep(NA_real_, nrow(df))
  kmer <- if ("assembly_support_best_contig_kmer_cov" %in% names(df)) df$assembly_support_best_contig_kmer_cov else rep(NA_real_, nrow(df))

  # assembly_exists (D-06/D-09): a contig is present when the join flags it
  # "supported" OR when both identity and length metrics are non-NA. Kept SEPARATE
  # from the score so no_own_assembly (exists FALSE) vs weak_own_assembly (exists
  # TRUE, low score) stays distinguishable downstream.
  metric_present <- !is.na(pid) & !is.na(len)
  if ("assembly_support" %in% names(df)) {
    supported <- !is.na(df$assembly_support) & df$assembly_support == "supported"
    assembly_exists <- supported | metric_present
  } else {
    assembly_exists <- metric_present
  }

  # Identity: logistic centered BELOW 90 (D-12) — smooth, no cliff (EVID-01).
  id_term  <- 1 / (1 + exp(-(pid - id_center) / id_slope))
  # Length: saturating partial credit, clamped to [0,1] (D-08).
  len_term <- pmax(0, pmin(len / len_ref, 1))
  base     <- wi * id_term + wl * len_term

  # k-mer: BONUS-ONLY, capped + log-compressed (D-10) — NA / <=0 => 0 boost, never a
  # penalty. Reuses the kmercov_cap precedent from score_candidates().
  kmer_capped <- ifelse(is.na(kmer) | kmer <= 0, 0, pmin(kmer, kmer_cap))
  kbonus      <- ifelse(kmer_capped > 0, wk * (log10(1 + kmer_capped) / log10(1 + kmer_cap)), 0)

  raw   <- base + kbonus
  # D-09: floor to 0 for no assembly / NA identity or length; D-08: bound to [0,1].
  floored <- !assembly_exists | is.na(pid) | is.na(len)
  score <- ifelse(floored, 0, pmin(raw, 1))
  score <- pmax(0, pmin(1, score))

  # EVID-05/D-03: surface the already-computed weighted per-metric contribution
  # terms (wi*id_term, wl*len_term, kbonus) so a reader can see each metric's push
  # into the score — ATTACHED, never recomputed. On a floored row (no contig / NA
  # identity or length) all three report 0 in lockstep with the floored 0 score, so
  # the reported contributions never leak NA on a scored row and never sum above a
  # zeroed score. On a scored row the three sum to the PRE-CLAMP score (base + kbonus
  # = raw), which is what "contribution to the score" means.
  contig_identity_contribution <- ifelse(floored, 0, wi * id_term)
  contig_length_contribution   <- ifelse(floored, 0, wl * len_term)
  contig_kmer_contribution     <- ifelse(floored, 0, kbonus)

  df %>%
    mutate(
      assembly_support_score       = score,
      assembly_exists              = assembly_exists,
      contig_identity_contribution = contig_identity_contribution,
      contig_length_contribution   = contig_length_contribution,
      contig_kmer_contribution     = contig_kmer_contribution
    )
}

# score_candidates(df, score_weights, evenness_const, kmercov_cap)
#   df             : candidate frame. Expected columns (NA-tolerant):
#                    candidate_reads (numeric) AND/OR targeted_reads_nodup (numeric,
#                    preferred for the reads term — the deduplicated targeted count;
#                    candidate_reads is the neutral first-mapping count which can be
#                    inverted by co-infection read mis-recruitment), the per-candidate breadth fraction
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

  # Reads term: prefer the deduplicated TARGETED read count over candidate_reads.
  # candidate_reads is the neutral all-reference first-mapping count (with
  # duplicates). In a co-infection that competitive mapping mis-recruits reads
  # between similar references, inverting the true abundance (sim2 70:30 case:
  # candidate_reads ranks the 30% minor above the 70% major). targeted_reads_nodup
  # is the deduplicated targeted re-mapping count and reflects true abundance, so
  # use it whenever it is present and valid. This mirrors the D2 tie-break in
  # classify_roles(), which already prefers targeted_reads_nodup over candidate_reads.
  # Fall back to candidate_reads per-row only where the targeted count is absent.
  cand_reads <- if ("candidate_reads" %in% names(df)) df$candidate_reads else rep(NA_real_, nrow(df))
  if ("targeted_reads_nodup" %in% names(df)) {
    tnodup <- df$targeted_reads_nodup
    reads <- ifelse(!is.na(tnodup) & tnodup > 0, tnodup, cand_reads)
  } else {
    reads <- cand_reads
  }
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
        dominance_score        = if ("dominance_score" %in% names(.)) dominance_score else double(),
        assembly_support_score = if ("assembly_support_score" %in% names(.)) assembly_support_score else double(),
        assembly_exists        = if ("assembly_exists" %in% names(.)) assembly_exists else logical(),
        # EVID-05/D-03: carry the three contribution columns on the zero-row path so a
        # header-only candidates.csv keeps the same schema as a populated run.
        contig_identity_contribution = if ("contig_identity_contribution" %in% names(.)) contig_identity_contribution else double(),
        contig_length_contribution   = if ("contig_length_contribution"   %in% names(.)) contig_length_contribution   else double(),
        contig_kmer_contribution     = if ("contig_kmer_contribution"     %in% names(.)) contig_kmer_contribution     else double(),
        evidence_state         = character(),
        role                   = character(),
        role_reason            = character(),
        overall_sample_call    = character()
      )
    return(out)
  }

  # Ensure a dominance_score column exists (defensive — caller normally scores first).
  if (!"dominance_score" %in% names(scored_df)) {
    scored_df <- score_candidates(scored_df)
  }

  # D-14: emit the continuous assembly_support_score + assembly_exists columns on
  # every output row (EVID-01). As of Plan 03, role / role_reason / overall_sample_call
  # are DERIVED from the per-candidate evidence_state (below). The binary ANDed
  # own_substantial floor is gone entirely as of 260805 — its last consumer was the
  # refuted band's quality re-check, and that band was unreachable.
  # score_assembly_support() reads only the candidate's OWN assembly metrics
  # (identity/length/k-mer, D-11).
  # WR-03 (12-REVIEW): score_assembly_support() always sets assembly_support_score
  # AND assembly_exists TOGETHER, so the recompute guard must check BOTH columns,
  # not just one — a caller-supplied frame carrying assembly_support_score without
  # assembly_exists (or vice versa) would otherwise skip recomputation and leave
  # assembly_exists NULL, silently zero-length-vectoring the ifelse() below into a
  # length-mismatch error instead of a clear diagnostic.
  if (!all(c("assembly_support_score", "assembly_exists") %in% names(scored_df))) {
    scored_df <- score_assembly_support(scored_df)
  }

  # 260805: the per-candidate `own_substantial` ANDed-floor computation was removed
  # here. Its last remaining consumer was the refuted band's quality re-check
  # (quality_fails_state), and that band is gone. The three denovo_min_* parameters
  # are RETAINED in the signature — callers pass them positionally, and
  # score_assembly_support() still uses the same metrics to build
  # assembly_support_score — but classify_roles() itself no longer applies them as a
  # floor. Do not mistake their presence in the signature for an active gate here.

  # --- EVID-02 / EVID-03 / EVID-04: per-candidate evidence_state ---------------
  # Computed for EACH candidate from its OWN assembly_support_score + assembly_exists
  # (Plan 01) and its OWN de novo contradiction — NEVER from another candidate's
  # dominance (EVID-02: the state is identical whether or not a stronger candidate
  # shares the sample; this vectorised outer computation cannot see dom_idx, which is
  # only determined per-sample inside classify_one_sample()). As of Plan 03, this
  # evidence_state is the SOLE driver of a non-dominant candidate's role /
  # role_reason (see the state -> role map in classify_one_sample()), which is what
  # decouples the sample-level co-infection call from dominance ordering (EVID-04).
  #
  # Band cutpoints are calibration-VALIDATED against the real 203-candidate dataset
  # (12-RESEARCH §4, D-13): the 15 genuine corroborated minors + the Sample51K/61K 2c
  # (~0.87) + 2714372 1a (~0.94) anchors all score >= hi_cut (confirmed); the
  # no-assembly floor scores 0 (weak). GLUE agreement is NOT read here — it can only
  # help reach confirmed, never gate it (D-05).
  evidence_hi_cut <- 0.72   # confirmed when assembly_support_score >= hi_cut (D-13)
  evidence_lo_cut <- 0.50   # probable in [lo_cut, hi_cut); weak below lo_cut (D-13)

  asup_score  <- scored_df$assembly_support_score
  asup_exists <- scored_df$assembly_exists

  # Bands (D-09/D-10/D-05): assembly_exists==FALSE forces weak regardless of score
  # (D-09); otherwise the assembly_support_score decides.
  #
  # 260805 (refuted-unreachable): the fourth band, `refuted`, has been REMOVED. It
  # required `denovo_contradicts & quality_fails_state`, where denovo_contradicts came
  # from the structurally unsatisfiable own_denovo_conflict() predicate — see
  # has_own_denovo()'s docstring. The band was dead code and evidence_state is now
  # three-valued: confirmed / probable / weak.
  #
  # It is not being repaired, and repairing it would be harmful. Note what the band
  # actually required: a contradicting contig that ALSO FAILS the substantiality
  # floors (quality_fails_state = !own_substantial). That is the profile of assembly
  # noise, not of a real second strain. The IVT dilution series demonstrates this
  # directly — nine controlled 1a/2a mixtures carry six off-genotype contigs at
  # 226-320 bp / 0.79-1.22x k-mer, every one failing own_substantial. Had attribution
  # been widened so these could be seen, each would have REFUTED a correct dominant
  # candidate on the strength of noise.
  #
  # The case genuinely worth catching — a SUBSTANTIAL contradicting contig — is owned
  # by RESCUE_EVALUATION (floors 3000 bp / 85% / 3000 bp aln / 2.0x k-mer), which is
  # unconditional and sample-scoped. A candidate whose only contig is off-genotype
  # still gets assembly_support = "none" here and lands in `weak`, so it is demoted
  # either way; only the label differs.
  evidence_state_col <- ifelse(
    !asup_exists,
    "weak",
    ifelse(asup_score >= evidence_hi_cut, "confirmed",
           ifelse(asup_score >= evidence_lo_cut, "probable", "weak"))
  )

  # Per-candidate floor pass (D-07/D-09: now informational annotation only, not a hard gate).
  reads <- scored_df$candidate_reads
  cov   <- if ("candidate_cov" %in% names(scored_df)) scored_df$candidate_cov else rep(NA_real_, nrow(scored_df))
  clears_floor <- !is.na(reads) & !is.na(cov) & reads > minRead & cov > minCov

  # New eligible pool: not discordant + has any coverage (breadth@>=1x sanity).
  has_concordance_outer <- "concordance_status" %in% names(scored_df)
  concordance_ok_outer <- if (!has_concordance_outer) rep(TRUE, nrow(scored_df)) else
    (is.na(scored_df$concordance_status) | scored_df$concordance_status != "discordant")
  eligible <- concordance_ok_outer & !is.na(cov) & cov > 0

  scored_df <- scored_df %>%
    mutate(
      below_floor      = clears_floor,
      .eligible        = eligible,
      .row_order       = row_number(),
      evidence_state   = evidence_state_col
    )

  # Group by sample so dominant determination + the asymmetric refute are per-sample.
  sample_key <- if ("sampleName" %in% names(scored_df)) "sampleName" else NULL

  classify_one_sample <- function(g) {
    n <- nrow(g)
    g$role <- NA_character_
    g$role_reason <- NA_character_

    has_concordance <- "concordance_status" %in% names(g)
    cov_vec <- if ("candidate_cov" %in% names(g)) g$candidate_cov else rep(NA_real_, nrow(g))
    concordance_ok <- if (!has_concordance) rep(TRUE, nrow(g)) else
      (is.na(g$concordance_status) | g$concordance_status != "discordant")
    eligible <- concordance_ok & !is.na(cov_vec) & cov_vec > 0

    # D-01: dominant = highest dominance_score among eligible candidates,
    # deterministic tie-break score -> reads -> ref name (D-06).
    gated <- which(eligible)
    dom_idx <- NA_integer_
    if (length(gated) > 0) {
      ord <- order(
        -g$dominance_score[gated],
        -g$candidate_reads[gated],
        as.character(g$candidate_ref[gated])
      )
      dom_idx <- gated[ord[1]]
    }

    dom_subtype  <- if (!is.na(dom_idx)) as.character(g$candidate_subtype[dom_idx])  else NA_character_
    dom_genotype <- if (!is.na(dom_idx)) as.character(g$candidate_genotype[dom_idx]) else NA_character_

    for (i in seq_len(n)) {
      if (!is.na(dom_idx) && i == dom_idx) {
        g$role[i] <- "dominant"
        g$role_reason[i] <- "dominant"
        next
      }

      if (has_concordance && !is.na(g$concordance_status[i]) &&
          g$concordance_status[i] == "discordant") {
        g$role[i]        <- "background"
        g$role_reason[i] <- "discordant_identity"
        next
      }

      # state -> role (EVID-04 / D-16): a non-dominant candidate's role is DERIVED
      # from its OWN per-candidate evidence_state (computed at outer scope, never
      # from dom_substantial), replacing the old own_substantial + dominance-
      # dependent asymmetric-refute branch. This is what lets a strong non-dominant
      # candidate (Sample51K-2c) surface as a co-infection MEMBER independent of
      # which candidate wins dominance, while a no-assembly candidate nets out
      # background regardless of whether the dominant assembled.
      #   confirmed / probable -> co-infection (subject to the UNCHANGED
      #                           is_valid_minor() demotion below; D-17)
      #   weak                 -> background, with the D-06/D-18 reason split:
      #                           no_own_assembly (assembly_exists FALSE) vs
      #                           weak_own_assembly_below_floor (present, below cut)
      #
      # 260805: the `refuted -> background / refuted_denovo` arm has been removed
      # with the unreachable evidence_state band that fed it. A candidate whose only
      # contig is off-genotype gets assembly_support = "none" -> weak -> background /
      # no_own_assembly, which is the same role by a differently-worded reason.
      st <- g$evidence_state[i]
      if (st %in% c("confirmed", "probable")) {
        role <- "co-infection"
        reason <- "corroborated"
      } else {
        role <- "background"
        reason <- if (isTRUE(g$assembly_exists[i])) "weak_own_assembly_below_floor" else "no_own_assembly"
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

    # D-14 updated (D3 fallback ladder).
    if (is.na(dom_idx)) {
      has_any_cov <- any(!is.na(cov_vec) & cov_vec > 0)
      call <- if (!has_any_cov) "untypable" else "indeterminate"
    } else if (any(g$role == "co-infection")) {
      call <- "co-infection"
    } else {
      call <- "monoinfection"
    }

    # D2: per-candidate indeterminate role trigger (§7.6 Q2).
    # Among non-background candidates, rank by (a) targeted nodup reads and
    # (b) best-contig k-mer coverage. If the #1 candidate differs -> both become
    # "indeterminate"; overall_sample_call -> "co-infection (indeterminate dominance)".
    non_bg <- which(!is.na(g$role) & g$role != "background")
    if (length(non_bg) >= 2 && !is.na(dom_idx)) {
      reads_col  <- if ("targeted_reads_nodup" %in% names(g)) "targeted_reads_nodup" else "candidate_reads"
      reads_vals <- g[[reads_col]][non_bg]
      kmer_vals  <- if ("assembly_support_best_contig_kmer_cov" %in% names(g))
                      g$assembly_support_best_contig_kmer_cov[non_bg]
                    else rep(NA_real_, length(non_bg))
      any_kmer   <- any(!is.na(kmer_vals) & kmer_vals > 0)
      if (any_kmer && !all(is.na(reads_vals))) {
        top_reads_local <- which.max(ifelse(is.na(reads_vals), -Inf, reads_vals))
        top_kmer_local  <- which.max(ifelse(!is.na(kmer_vals) & kmer_vals > 0, kmer_vals, -Inf))
        if (top_reads_local != top_kmer_local) {
          top_reads_idx <- non_bg[top_reads_local]
          top_kmer_idx  <- non_bg[top_kmer_local]
          g$role[top_reads_idx]        <- "indeterminate"
          g$role_reason[top_reads_idx] <- "indeterminate_dominance_conflict"
          g$role[top_kmer_idx]         <- "indeterminate"
          g$role_reason[top_kmer_idx]  <- "indeterminate_dominance_conflict"
          call <- "co-infection (indeterminate dominance)"
        }
      }
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
    select(-.eligible, -.row_order)
}

# =============================================================================
# Phase 13 interpretability helpers (EVID-05 / EVID-06)
#
# Pure, scalar-per-candidate functions that turn the already-computed
# role / role_reason / evidence_state / concordance_status + measured contig
# metrics into human-readable, contig-language text. They RECOMPUTE nothing
# (D-03) and read only fields already present on candidate_support after
# classify_roles(). Language rule D-01: always "contig", never "own
# assembly"/"own support". Severity is conveyed by concrete measured values,
# never editorial words (D-07).
# =============================================================================

# Candidate identity token per D-05: "candidate <rank> (<subtype>_<ref>)".
# Production reference names already carry the subtype prefix (e.g. "2c_JX227949"),
# so we avoid doubling it when ref already starts with "<subtype>_"; if ref is the
# bare accession we prepend the subtype. NA-tolerant (T-13-01): never stop().
.candidate_label <- function(rank, ref, subtype) {
  ref_s <- if (length(ref) == 0 || is.na(ref)) NA_character_ else as.character(ref)
  sub_s <- if (length(subtype) == 0 || is.na(subtype)) NA_character_ else as.character(subtype)
  ident <- if (!is.na(ref_s) && !is.na(sub_s)) {
    if (startsWith(ref_s, paste0(sub_s, "_")) || identical(ref_s, sub_s)) ref_s
    else paste0(sub_s, "_", ref_s)
  } else if (!is.na(ref_s)) {
    ref_s
  } else if (!is.na(sub_s)) {
    sub_s
  } else {
    "unknown"
  }
  rank_s <- if (length(rank) == 0 || is.na(rank)) "?" else as.character(rank)
  paste0("candidate ", rank_s, " (", ident, ")")
}

# Measured contig metrics as a plain, factual phrase (D-02 pattern 3 / D-07): emits
# ONLY the tokens whose source value is non-NA — never renders a literal "NA"
# (T-13-02). Returns "" when no metric is present so callers can test length.
.contig_metrics_phrase <- function(len, pid, kmer) {
  toks <- character(0)
  if (length(len)  == 1 && !is.na(len))  toks <- c(toks, paste0("contig length ", round(len), " bp"))
  if (length(pid)  == 1 && !is.na(pid))  toks <- c(toks, paste0("identity ", sprintf("%.1f", pid), "%"))
  if (length(kmer) == 1 && !is.na(kmer)) toks <- c(toks, paste0("k-mer coverage ", sprintf("%.2f", kmer)))
  paste(toks, collapse = ", ")
}

# build_evidence_summary(): one contig-language sentence per candidate stating its
# role + evidence_state + contig corroboration/conflict (EVID-05/D-01/D-02). Scalar;
# apply row-wise via purrr::pmap_chr(). Reads ONLY contig support (never
# reads/cov/evenness — D-11 spirit). Every reachable role_reason from
# classify_one_sample() maps to a non-empty, candidate-named sentence; an unmapped
# reason falls through to a defined generic sentence (never empty, never stop()).
build_evidence_summary <- function(role, role_reason, evidence_state,
                                   concordance_status,
                                   candidate_rank, candidate_ref, candidate_subtype,
                                   assembly_exists,
                                   best_contig_length = NA_real_,
                                   best_contig_pident = NA_real_,
                                   best_contig_kmer_cov = NA_real_) {
  label <- .candidate_label(candidate_rank, candidate_ref, candidate_subtype)
  metrics <- .contig_metrics_phrase(best_contig_length, best_contig_pident, best_contig_kmer_cov)
  st  <- if (length(evidence_state) == 1 && !is.na(evidence_state)) as.character(evidence_state) else NA_character_
  rr  <- if (length(role_reason)    == 1 && !is.na(role_reason))    as.character(role_reason)    else NA_character_
  state_tok <- if (!is.na(st)) paste0(" (evidence_state=", st, ")") else ""
  # A trailing "; contig <metrics>" clause when a matching contig is present.
  contig_clause <- if (nzchar(metrics)) paste0("; contig ", metrics) else ""

  body <- if (identical(rr, "dominant")) {
    paste0("dominant candidate", if (nzchar(metrics)) paste0(", contig ", metrics) else "")
  } else if (identical(rr, "corroborated")) {
    paste0("co-infection member corroborated by a matching-genotype contig",
           if (nzchar(metrics)) paste0(" (", metrics, ")") else "")
  } else if (identical(rr, "discordant_identity")) {
    paste0("background — mapping identity conflicts with the contig/GLUE identity", contig_clause)
  } else if (identical(rr, "no_own_assembly")) {
    "background — no contig with the same genotype/subtype as the candidate"
  } else if (identical(rr, "weak_own_assembly_below_floor")) {
    paste0("background — a matching contig was found but is weak: ",
           if (nzchar(metrics)) metrics else "no measurable contig support")
  } else if (identical(rr, "same_genotype_as_dominant")) {
    paste0("demoted to background — same genotype as the dominant candidate",
           if (nzchar(metrics)) paste0(", though its own contig evidence was ", metrics) else "")
  } else if (identical(rr, "recombinant_2k1b")) {
    "demoted to background — 2k/1b recombinant pair (reported as genotype 1 or 2, no 2k1b GLUE clade)"
  } else if (identical(rr, "indeterminate_dominance_conflict")) {
    "indeterminate — read-count and k-mer-coverage dominance rankings disagree"
  } else {
    role_s <- if (length(role) == 1 && !is.na(role)) as.character(role) else "unclassified"
    paste0(role_s, if (nzchar(metrics)) paste0("; contig ", metrics) else "")
  }

  paste0(label, ": ", body, state_tok, ".")
}

# candidate_review_fragment(): per-candidate review_flag fragment (EVID-06). Returns
# NA_character_ for a clean dominant / clean confirmed / clean background candidate,
# else "candidate <rank> (<subtype>_<ref>): <reason with the concrete driving value>".
# Fires for EVERY non-dominant, non-clean candidate — weak/probable evidence
# states and the D-08 demotions (same_genotype_as_dominant / recombinant_2k1b) that
# fire even though the candidate's own contig evidence was good. NO numeric flag-gating
# threshold is introduced (D-06); severity is conveyed by the concrete measured value in
# the wording (D-07). Scalar; apply row-wise via pmap_chr, then collapse per sample.
candidate_review_fragment <- function(role, role_reason, evidence_state, concordance_status,
                                      candidate_rank, candidate_ref, candidate_subtype,
                                      best_contig_length = NA_real_,
                                      best_contig_pident = NA_real_,
                                      best_contig_kmer_cov = NA_real_) {
  rr     <- if (length(role_reason)    == 1 && !is.na(role_reason))    as.character(role_reason)    else NA_character_
  st     <- if (length(evidence_state) == 1 && !is.na(evidence_state)) as.character(evidence_state) else NA_character_
  role_s <- if (length(role)           == 1 && !is.na(role))           as.character(role)           else NA_character_

  # The dominant candidate's review is surfaced at the SAMPLE level
  # (dominant_unconfirmed), and the dominance-ordering conflict is a generic sample
  # message (D-10) — neither yields a per-candidate fragment here.
  if (identical(role_s, "dominant") || identical(rr, "dominant")) return(NA_character_)
  if (identical(rr, "indeterminate_dominance_conflict")) return(NA_character_)

  is_demotion  <- !is.na(rr) && rr %in% c("same_genotype_as_dominant", "recombinant_2k1b")
  # 260805: "refuted_denovo" removed from is_conflict and "refuted" from is_weakstate —
  # both were unreachable evidence_state/role_reason values. discordant_identity stays:
  # it is LIVE via apply_concordance()'s GLUE leg.
  is_conflict  <- !is.na(rr) && rr %in% c("discordant_identity")
  is_weakstate <- !is.na(st) && st %in% c("weak", "probable")

  # D-06: flag every non-dominant, non-clean candidate. A clean confirmed co-infection
  # / clean background (no flaggable state, not a demotion or conflict) => NA.
  if (!(is_demotion || is_conflict || is_weakstate)) return(NA_character_)

  label   <- .candidate_label(candidate_rank, candidate_ref, candidate_subtype)
  metrics <- .contig_metrics_phrase(best_contig_length, best_contig_pident, best_contig_kmer_cov)
  m_paren <- if (nzchar(metrics)) paste0(" (", metrics, ")") else ""

  reason <- if (identical(rr, "weak_own_assembly_below_floor")) {
    paste0("matching contig but weak — ", if (nzchar(metrics)) metrics else "no measurable contig support")
  } else if (identical(rr, "no_own_assembly")) {
    "no contig with the same genotype/subtype as the candidate"
  } else if (identical(rr, "corroborated") && identical(st, "probable")) {
    paste0("co-infection only marginally corroborated by contig", m_paren)
  } else if (identical(rr, "discordant_identity")) {
    paste0("mapping identity conflicts with the contig/GLUE identity", m_paren)
  } else if (identical(rr, "same_genotype_as_dominant")) {
    paste0("demoted — same genotype as the dominant candidate",
           if (nzchar(metrics)) paste0(", though its own contig evidence was ", metrics) else "")
  } else if (identical(rr, "recombinant_2k1b")) {
    paste0("demoted — 2k/1b recombinant pair",
           if (nzchar(metrics)) paste0(", though its own contig evidence was ", metrics) else "")
  } else {
    rr_tok <- if (!is.na(rr)) rr else "flagged"
    paste0(rr_tok, if (nzchar(metrics)) paste0(" — ", metrics) else "")
  }

  # Every fragment carries the evidence_state token so severity is reconstructable.
  suffix <- if (!is.na(st) && !grepl("evidence_state", reason, fixed = TRUE)) paste0(" [evidence_state=", st, "]") else ""
  paste0(label, ": ", reason, suffix)
}

# offgenotype_contig_reviewable(): gate for the monoinfection "de novo assembly
# found a different-genotype contig" review trigger (quick task 260803-ogc).
#
# The trigger used to fire on ANY off-genotype contig with no substantiality floor
# whatsoever: 72 of 140 samples (51%) across the five v1.3.0-g28a568d routine runs,
# median triggering contig 606 bp, shortest 142 bp — below even
# denovo_min_contig_length (500). It was the ONLY sentence on all 72, so it alone
# accounted for 92% of the cohort's `provisional` calls: the one signal meant to mark
# a possibly missed co-infection was firing mostly on assembly noise.
#
# Two independent gates, both settled empirically by the 2026-08-03 threshold
# sweep over that cohort (bin/tests/offgeno_flag_sweep.R re-derives it; the sweep
# report itself is not in the repo, as it quotes diagnostic sample identifiers):
#
#   (1) CONTIG LENGTH >= min_length. Length is the ONLY usable leg. The three
#       legacy-typable minors this build demotes to monoinfection — the cases that
#       most deserve the flag — measure 4467/2787/2706 bp but sit at k-mer coverage
#       1.42/1.97/1.62, i.e. BELOW denovo_min_kmer_cov (2.0). A k-mer leg would
#       suppress exactly the samples worth reviewing, and the artefacts run the other
#       way (2k1b fragments co-assembled off an abundant major reach k-mer 2797-3257).
#       Applied by the CALLER (summarize.R), which owns the param and the length
#       column; passing min_length here keeps the helper self-contained for tests.
#
#   (2) NOT a 2k/1b recombinant pair. 2k1b is a 2k/1b recombinant reference: its 1b
#       portion is not a different genotype from a 1b major, so a 2k1b contig against
#       a genotype-1 or -2 major is a taxonomy artefact rather than a co-infection.
#       This is the SAME policy is_valid_minor() rule 2 (above) already applies to
#       candidate promotion — the review trigger simply never consulted it. 14 of the
#       72 flags are this class.
#
# WHY NOT genotype_from_subtype() ALONE: it maps "2k1b" -> "2k1b" (the 2k1b-aware
# rule), so "2k1b" != "1" and a 2k1b-vs-1b pair STILL differs. Swapping the old
# substr(x, 1, 1) comparison for genotype_from_subtype() removes zero flags, and for
# a 2k1b contig against a 2c major it ADDS one that substr never fired ("2" == "2").
# The explicit pair test below is what does the work. Using genotype_from_subtype()
# for the difference test is still correct — combined with the pair test it is
# behaviourally identical to substr for every 2k1b case — and it retires the
# hand-rolled idiom the codebase replaced everywhere else.
#
# WHY NOT is_valid_minor() WHOLESALE: rule 1 returns TRUE for 1a/1b pairs, so
# delegating the whole test to it would newly fire this trigger on the
# within-genotype-1 cross-mapping artefacts (29 samples) that the role model already
# discards correctly.
#
# NA contig_length fails OPEN (keeps the flag): an unmeasurable contig is not
# evidence of insubstantiality, and silently dropping a possible co-infection is the
# worse error. A batch with no de novo leg never reaches here — contig_subtype is NA
# and the difference test below is FALSE.
#
# Vectorised (called from a mutate over the whole sample frame). Returns a bare
# logical with no NAs.
offgenotype_contig_reviewable <- function(contig_subtype, major_subtype,
                                          contig_length = NA_real_,
                                          min_length = 0) {
  cg <- genotype_from_subtype(contig_subtype)
  mg <- genotype_from_subtype(major_subtype)

  different_genotype <- !is.na(cg) & !is.na(mg) & cg != mg

  # is_valid_minor() rule 2, restated over genotypes only.
  pair_2k1b <- (cg == "2k1b" & mg %in% c("1", "2", "2k1b")) |
               (mg == "2k1b" & cg %in% c("1", "2", "2k1b"))
  pair_2k1b <- !is.na(pair_2k1b) & pair_2k1b

  long_enough <- is.na(contig_length) | contig_length >= min_length

  different_genotype & !pair_2k1b & long_enough
}

# contig_evidence_note(): the measured-evidence clause for a review sentence that
# makes a claim about a de novo contig (quick task 260803-ogc, option C).
#
# Used by BOTH monoinfection contig triggers, via `context`:
#   "offgenotype"    — "de novo found a different-genotype contig (X)"
#   "major_conflict" — "Major subtype conflict — mapping (X) vs contig (Y)"
# The metrics render identically; only the interpretation clause differs, because a
# poorly-aligned contig means something different in each case (a phantom second
# strain vs a phantom disagreement).
#
# The sentence used to name a subtype and nothing else — "de novo assembly found a
# different-genotype contig (6i)" — and then ask a human to review it. None of the
# four numbers needed to judge that claim were anywhere the analyst would look: the
# aligned length lives ONLY in blastparse/<sample>.assembly_support.csv, because the
# candidate-grain join never reaches a subtype that is not itself a candidate.
#
# Worked example (2633901, run 20251212-01, major 1a). A 1,620 bp contig was reported
# as a genotype-6i off-genotype contig. Its BLAST alignment is 69 bp — 4% of the
# contig. Confirmed independently against nt: the contig's only HCV-like region is a
# ~212 bp tail (13% of it), closest to 1a — the SAME genotype as the major. So the
# contig is ~87% non-HCV, the "6i" label is an artefact of a short anchor, and there
# is no second strain. A reader given "6i" alone cannot possibly reach that
# conclusion; a reader given "69 bp aligned (4%)" reaches it immediately.
#
# Why the aligned FRACTION is the discriminator, not the contig length: this contig
# is 1,620 bp, comfortably over any length floor. Length says the contig is real;
# the aligned fraction says how much of it is actually HCV. The three genuine minors
# this build demotes align over 99-100% of their contigs (2705/2706, 2762/2787,
# 4448/4467) against 4% here — a ~25x separation, so the min_aln_frac cut point is
# not sensitive.
#
# min_aln_frac only changes the WORDING, never whether the sentence fires. Kept as a
# plain default rather than a plumbed param: it is presentational, and the sweep
# output (aln_frac_pct) should inform any tuning before it earns an arg position.
#
# Returns NA_character_ when no metric is known, so the caller falls back to the
# original wording. Scalar (called via pmap_chr), mirroring build_evidence() and
# candidate_review_fragment().
contig_evidence_note <- function(contig_length, aln_length, pident, kmer_cov,
                                 context = "offgenotype", min_aln_frac = 0.5) {
  one <- function(x) if (length(x) == 0) NA else x[[1]]
  len <- suppressWarnings(as.numeric(one(contig_length)))
  aln <- suppressWarnings(as.numeric(one(aln_length)))
  pid <- suppressWarnings(as.numeric(one(pident)))
  km  <- suppressWarnings(as.numeric(one(kmer_cov)))

  frac <- if (!is.na(len) && !is.na(aln) && len > 0) aln / len else NA_real_

  toks <- character(0)
  if (!is.na(len)) toks <- c(toks, paste0(round(len), " bp contig"))
  if (!is.na(aln)) toks <- c(toks, paste0(round(aln), " bp aligned",
                                          if (!is.na(frac)) paste0(" (", round(100 * frac), "%)") else ""))
  if (!is.na(pid)) toks <- c(toks, paste0(formatC(pid, format = "f", digits = 1), "% identity"))
  if (!is.na(km))  toks <- c(toks, paste0("k-mer cov ", formatC(km, format = "f", digits = 1)))
  if (length(toks) == 0) return(NA_character_)

  weak    <- !is.na(frac) && frac < min_aln_frac
  is_conf <- identical(context, "major_conflict")
  interp <- if (weak && is_conf) {
    # A short anchor cannot support a subtype call, so it cannot support a
    # DISAGREEMENT with one either. Say that, rather than implying a real conflict.
    paste0("only ", round(100 * frac), "% of the contig aligns to any reference in the panel, so the ",
           "contig subtype is weakly supported and the conflict may be an artefact of a short anchor ",
           "rather than a real discrepancy")
  } else if (weak) {
    paste0("only ", round(100 * frac), "% of the contig aligns to any reference in the panel, so the ",
           "subtype assignment is weakly supported — the contig may be largely non-HCV, chimeric, ",
           "or too divergent to type")
  } else if (is_conf) {
    "possible reference mismatch or highly divergent strain"
  } else {
    "possible missed co-infection or contamination"
  }
  paste0("— ", paste(toks, collapse = ", "), "; ", interp, ".")
}

# sample_review_message(): the rewritten sample-level review_flag builder (EVID-06).
# Keeps the D-10 triggers (is_indet, gate_flag != "ok", is_indet_dom) GENERIC with no
# candidate name; enriches the subtype-conflict / different-genotype-contig triggers
# with the named candidate + the actual conflicting subtype values (D-02 pattern 2,
# sample-level per RESEARCH Pitfall 3); resolves D-11 by NAMING the candidate for
# dominant_unconfirmed and rescue_effect == "major_ref_changed"; accepts the
# pre-collapsed per-candidate fragment string (from candidate_review_fragment) as one
# argument and merges everything with " | "; returns NA_character_ when nothing fires
# (the MultiQC NA sentinel — Pitfall 2: filter !is.na before paste, never render "NA").
sample_review_message <- function(overall_sample_call,
                                  denovo_major_subtype_match,
                                  denovo_minor_subtype_match,
                                  gate_flag,
                                  denovo_minor_subtype,
                                  denovo_major_subtype,
                                  major_subtype,
                                  rescue_effect,
                                  dominant_unconfirmed,
                                  dominant_rank = NA_integer_,
                                  dominant_ref = NA_character_,
                                  candidate_fragment = NA_character_,
                                  # 260803-ogc option C: pre-built measured-evidence clause for the
                                  # different-genotype-contig sentence (offgenotype_contig_note()).
                                  # Trailing arg WITH a default so every existing positional and
                                  # named caller keeps working unchanged; NA falls back to the
                                  # original wording.
                                  offgeno_note = NA_character_,
                                  # 260803-ogc follow-up 1: the same clause for the
                                  # MAJOR subtype-conflict sentence, describing the
                                  # contig behind denovo_major_subtype. Trailing arg
                                  # with a default for the same back-compatibility
                                  # reason as offgeno_note.
                                  majconf_note = NA_character_) {
  one <- function(x) if (length(x) == 0) NA else x[[1]]
  sc        <- one(overall_sample_call)
  maj_match <- one(denovo_major_subtype_match)
  # 260804-dnrank: denovo_minor_subtype_match no longer has a reader here — its only
  # consumer was the removed co-infection dominance message. The PARAMETER is kept
  # (it is the 3rd positional arg; dropping it would silently shift every positional
  # caller) but nothing binds it locally any more.
  gflag     <- one(gate_flag)
  resc      <- one(rescue_effect)
  dv_minor  <- one(denovo_minor_subtype)
  dv_major  <- one(denovo_major_subtype)
  maj_sub   <- one(major_subtype)
  dom_unconf <- one(dominant_unconfirmed)
  cf        <- one(candidate_fragment)
  ogn       <- one(offgeno_note)
  mcn       <- one(majconf_note)
  tok <- function(x) if (length(x) == 0 || is.na(x)) "unknown" else as.character(x)

  msgs <- character(0)
  is_mono      <- !is.na(sc) && sc == "monoinfection"
  is_indet     <- !is.na(sc) && sc == "indeterminate"
  is_indet_dom <- !is.na(sc) && sc == "co-infection (indeterminate dominance)"
  dom_label    <- .candidate_label(dominant_rank, dominant_ref, maj_sub)

  # D-10 generic (no candidate name).
  if (is_indet_dom)
    msgs <- c(msgs, "Dominance ordering uncertain — read-count and k-mer-coverage rankings disagree. Both genotypes reported as present; co-infection vs contamination agnostic. Please review.")
  # Co-infection subtype conflict (sample-level de novo vs mapping).
  #
  # 260804-dnrank: REMOVED. This asserted a disagreement about DOMINANCE, but
  # denovo_major_ref / denovo_minor_ref are ordered by BLAST bitscore
  # (blast_parse.R:301, frame sorted at :150-153), which ranks strains by proximity
  # to the reference PANEL, not by abundance. On the designed 7:3 sim2 mixture it
  # handed "major" to 3a on a 100% panel match (bitscore 17,444) despite 2.4x less
  # k-mer coverage and 2.2x fewer mapped reads than 2a, and raised
  # call_confidence = review on a clean call.
  #
  # Dominance is already owned by the D2 trigger above, which compares
  # top-by-targeted-reads against top-by-best-contig-k-mer-coverage -- both
  # abundance measures -- and emits "Dominance ordering uncertain". D2 correctly
  # stays silent on sim2. Keeping two dominance checks that measure different
  # things guaranteed they would disagree.
  #
  # denovo_*_subtype_match is retained for what it CAN support: the monoinfection
  # major-subtype cross-check immediately below, where a single strain means no
  # dominance question arises.
  # Monoinfection major subtype conflict: name the candidate + the actual subtype values (D-02 pattern 2).
  # 260803-ogc follow-up 1: carry the conflicting contig's measured evidence. This
  # trigger sets call_confidence = "review" on its own, so it is the HARDEST signal
  # the sample-level builder emits — and until now it named two subtypes and no
  # numbers. The metrics it needs are not in Major_best_contig_*, which describes the
  # CANDIDATE's genotype group; the contig that caused the disagreement belongs to
  # denovo_major_subtype, a different group entirely.
  if (is_mono && !is.na(maj_match) && maj_match == "NO")
    msgs <- c(msgs, if (!is.na(mcn) && nzchar(mcn)) paste0(
      "Major subtype conflict for ", dom_label, " — mapping (", tok(maj_sub),
      ") vs contig (", tok(dv_major), ") ", mcn, " Please review."
    ) else paste0(
      "Major subtype conflict for ", dom_label, " — mapping (", tok(maj_sub),
      ") vs contig (", tok(dv_major), "). Possible reference mismatch or highly divergent strain. Please review."))
  # Monoinfection different-genotype contig: names the contig subtype value.
  # 260803-ogc: the genotype-difference + 2k1b policy now lives in the pure
  # offgenotype_contig_reviewable() helper above (replacing a hand-rolled
  # substr(x, 1, 1)). The CONTIG-LENGTH leg is applied by the caller, which masks
  # dv_minor to NA when the contig is below review_min_offgenotype_contig_length —
  # so a sub-floor contig short-circuits on the is.na() guard here.
  if (is_mono && !is.na(dv_minor) && !is.na(maj_sub) &&
      offgenotype_contig_reviewable(dv_minor, maj_sub))
    msgs <- c(msgs, paste0(
      "Monoinfection called for ", dom_label, ", but de novo assembly found a different-genotype contig (",
      dv_minor, ") ",
      # Measured evidence when the metrics resolved; the original bare wording otherwise.
      if (!is.na(ogn) && nzchar(ogn)) ogn else "— possible missed co-infection or contamination.",
      " Please review."))
  # D-11: name the dominant candidate for a provisional (uncorroborated) call.
  if (isTRUE(dom_unconf) && is_mono)
    msgs <- c(msgs, paste0("Genotype call for ", dom_label,
      " is provisional — identity not corroborated (mapping evidence only; no GLUE or de novo confirmation). Please review."))
  # D-10 generic (no candidate name).
  if (is_indet)
    msgs <- c(msgs, "No candidate passed the major-gate — overall sample call indeterminate.")
  if (!is.na(gflag) && gflag != "ok")
    msgs <- c(msgs, "Major strain failed mapping quality thresholds — genotype call uncertain.")
  # D-11: name the reassigned dominant candidate for a Major-slot rescue.
  if (!is.na(resc) && resc == "major_ref_changed")
    msgs <- c(msgs, paste0("De novo rescue overrode the dominant/Major reference for ", dom_label,
      " — the primary call was reassigned automatically. Confirm against rescue_audit.csv (from/to + trigger) before reporting. Please review."))
  # Merge the pre-collapsed per-candidate fragment(s), NA-filtered (Pitfall 2).
  if (!is.na(cf) && nzchar(cf)) msgs <- c(msgs, cf)

  if (length(msgs) == 0) NA_character_ else paste(msgs, collapse = " | ")
}
