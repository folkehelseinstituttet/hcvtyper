#!/usr/bin/env Rscript

# offgeno_flag_sweep.R ----------------------------------------------------
# Sweep the contig-length floor for the monoinfection "de novo assembly found a
# different-genotype contig" review trigger (quick task 260803-ogc), on the column
# the gate actually reads (denovo_minor_contig_length), with the comparison
# report's assembly_support best_contig_length arm reproduced alongside.
#
# This is the tool that chose review_min_offgenotype_contig_length = 1000. Kept in
# the repo so the threshold can be re-derived on a new cohort rather than argued
# from memory. READ-ONLY: touches nothing but its own output directory.
#
# Deliberately not named test_*.R — run_all.sh globs that pattern and runs each
# file with no arguments; this one needs a data root.
#
# Usage:
#   Rscript bin/tests/offgeno_flag_sweep.R <root> [outdir] [runs] [expected_n]
#
#   root        directory containing the run folders (each with summary/Summary.csv)
#   outdir      default ./offgeno_sweep_out
#   runs        comma-separated run directory names. Default: the five
#               v1.3.0-g28a568d routine runs. Pass "" to accept every run under
#               <root> (and then set expected_n to match, or 0 to disable the check).
#   expected_n  expected sample count; the script ABORTS on a mismatch. Default 140,
#               0 disables.
#
# Two bugs in the first version of this script, both fixed here, both of which
# produced plausible output and exit 0:
#
#   BUG 1 — the recursive Summary.csv glob picked up every result directory under
#   the root, not just the cohort: 10 runs / 437 samples instead of 5 / 140,
#   including benchmark dirs built from different pipeline versions. Fixed by the
#   explicit `runs` allowlist AND a hard abort on the cohort-size check, which
#   previously only printed the expected value as a comment.
#
#   BUG 2 — sampleName in Summary.csv carries the "-HCV" suffix
#   (summarize.R derives it as the first dot-separated token of the trimmed-read
#   basename), while the assembly_support prefix does too but was stripped on only
#   ONE side of the join. The comparison arm silently resolved to zero rows for the
#   whole cohort, and MUST_KEEP (bare ids) never matched, so must_keep_kept read
#   0/3 in every sweep row — the disqualifying criterion the exercise hangs on never
#   evaluated. Fixed by deriving a bare `sample_id` once and keying everything on it,
#   plus a hard abort when the must-keep samples are not all found.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args       <- commandArgs(trailingOnly = TRUE)
root       <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
outdir     <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "offgeno_sweep_out"
runs_arg   <- if (length(args) >= 3) args[3] else NA_character_
expected_n <- if (length(args) >= 4 && nzchar(args[4])) as.integer(args[4]) else 140L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

DEFAULT_RUNS <- c("NGS_SEQ-20251113-02", "NGS_SEQ-20251212-01", "NGS_SEQ-20260326-01",
                  "NGS_SEQ-20260521-01", "NGS_SEQ-20260625-02")
runs_wanted <- if (is.na(runs_arg)) DEFAULT_RUNS else
  if (nzchar(runs_arg)) trimws(str_split(runs_arg, ",")[[1]]) else character(0)

say  <- function(...) cat(..., "\n", sep = "")
rule <- function(t) say("\n", strrep("=", 72), "\n", t, "\n", strrep("=", 72))
die  <- function(...) { say("\n*** ABORT: ", ..., " ***"); quit(status = 1) }

MUST_KEEP <- c("2743986", "2726018", "2714375")
FLOORS    <- c(0, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000)
TRIG_TEXT <- "different-genotype contig"

## -- 1. Load the cohort ---------------------------------------------------
# All-character read: Major_genotype/Minor_genotype were absent from 2 of the 5
# runs (fixed separately in this quick task) and per-file type inference disagrees
# across runs — either would abort a typed bind_rows. Coerce numerics after.
sum_files <- list.files(root, pattern = "^Summary\\.csv$", recursive = TRUE, full.names = TRUE)
sum_files <- sum_files[basename(dirname(sum_files)) == "summary"]
if (length(sum_files) == 0) die("no summary/Summary.csv found under: ", root)

run_of <- function(f) basename(dirname(dirname(f)))
if (length(runs_wanted) > 0) {
  keep <- run_of(sum_files) %in% runs_wanted     # BUG 1 fix
  missing_runs <- setdiff(runs_wanted, run_of(sum_files))
  if (length(missing_runs)) die("requested run(s) not found under ", root, ": ",
                                paste(missing_runs, collapse = ", "))
  say("Excluded ", sum(!keep), " run dir(s) not in the cohort allowlist: ",
      paste(sort(unique(run_of(sum_files)[!keep])), collapse = ", "))
  sum_files <- sum_files[keep]
}

sm <- map_dfr(sum_files, function(f) {
  read_csv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    mutate(run = run_of(f), .before = 1)
})

rule("1. Cohort")
say("Runs  : ", paste(sort(unique(sm$run)), collapse = ", "))
say("Samples: ", nrow(sm))
if (expected_n > 0 && nrow(sm) != expected_n)
  die("loaded ", nrow(sm), " samples but expected ", expected_n,
      ". Fix the run allowlist (arg 3) or the expected count (arg 4) before trusting any number below.")

if (!"Major_subtype" %in% names(sm) && "Major" %in% names(sm)) {
  sm <- mutate(sm, Major_subtype = Major)
  say("NOTE: Major_subtype absent; using Major (identical by construction).")
}
need <- c("sampleName", "overall_sample_call", "call_confidence", "review_flag",
          "Major_subtype", "denovo_minor_subtype", "denovo_minor_ref",
          "denovo_minor_contig_length")
missing <- setdiff(need, names(sm))
if (length(missing)) die("Summary.csv is missing required columns: ",
                         paste(missing, collapse = ", "))

# BUG 2 fix: one bare id, derived once, used for every join and every lookup.
sm <- sm %>% mutate(
  L_summary = suppressWarnings(as.numeric(denovo_minor_contig_length)),
  sample_id = str_remove(sampleName, "-HCV.*$")
)

## -- 2. Reconstruct the trigger, and check it against what actually fired --
geno1 <- function(x) substr(x, 1, 1)
sm <- sm %>% mutate(
  pred_fires = !is.na(overall_sample_call) & overall_sample_call == "monoinfection" &
               !is.na(denovo_minor_subtype) & !is.na(Major_subtype) &
               geno1(denovo_minor_subtype) != geno1(Major_subtype),
  text_fires = !is.na(review_flag) & str_detect(review_flag, fixed(TRIG_TEXT)),
  n_sentences = if_else(is.na(review_flag), 0L, str_count(review_flag, fixed(" | ")) + 1L)
)

rule("2. Trigger reconstruction")
say("Predicate (classify_roles.R, pre-260803-ogc) fires on : ", sum(sm$pred_fires))
say("review_flag text actually contains the sentence       : ", sum(sm$text_fires))
disagree <- sm %>% filter(pred_fires != text_fires)
if (nrow(disagree) > 0) {
  say("\n*** WARNING: predicate and text disagree on ", nrow(disagree), " sample(s) —",
      " the reconstruction is not exact. ***")
  print(disagree %>% select(run, sampleName, overall_sample_call, Major_subtype,
                            denovo_minor_subtype, pred_fires, text_fires))
  write_csv(disagree, file.path(outdir, "predicate_text_disagreement.csv"))
} else {
  say("Predicate and text agree on all ", nrow(sm), " samples. Reconstruction is exact.")
}

flagged <- sm %>% filter(text_fires)
say("\nSole sentence on every flagged sample? ",
    if (all(flagged$n_sentences == 1)) "YES" else
      paste0("NO — ", sum(flagged$n_sentences > 1), " carry additional sentences"))

## -- 3. Comparison arm: assembly_support best_contig_* --------------------
asup_files <- list.files(root, pattern = "\\.assembly_support\\.csv$",
                         recursive = TRUE, full.names = TRUE)
if (length(runs_wanted) > 0) asup_files <- asup_files[run_of(asup_files) %in% runs_wanted]
have_asup <- length(asup_files) > 0

if (have_asup) {
  asup <- map_dfr(asup_files, function(f) {
    read_csv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
      mutate(run = run_of(f), .before = 1)
  }) %>%
    mutate(sample_id = str_remove(sample, "-HCV.*$"),          # BUG 2 fix
           L_asup    = suppressWarnings(as.numeric(best_contig_length)),
           pident    = suppressWarnings(as.numeric(best_contig_pident)),
           kmer_cov  = suppressWarnings(as.numeric(best_contig_kmer_cov)),
           # Carried through per the results write-up: 2633901 has a 1620 bp contig
           # with a 69 bp alignment (4.3%), which a length-only floor lets through.
           aln_len   = suppressWarnings(as.numeric(best_contig_aln_length))) %>%
    group_by(run, sample_id, subtype) %>%
    slice_max(L_asup, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(run, sample_id, subtype, L_asup, pident, kmer_cov, aln_len)

  flagged <- flagged %>%
    left_join(asup, by = c("run", "sample_id", "denovo_minor_subtype" = "subtype")) %>%
    mutate(aln_frac_pct = round(100 * aln_len / L_asup, 1))
} else {
  flagged <- flagged %>% mutate(L_asup = NA_real_, pident = NA_real_,
                                kmer_cov = NA_real_, aln_len = NA_real_,
                                aln_frac_pct = NA_real_)
}

rule("3. The two length definitions")
say("assembly_support.csv files: ", length(asup_files),
    if (!have_asup) "  -> comparison arm UNAVAILABLE" else "")
say("Flagged with L_summary : ", sum(!is.na(flagged$L_summary)), " / ", nrow(flagged))
say("Flagged with L_asup    : ", sum(!is.na(flagged$L_asup)),    " / ", nrow(flagged))

both <- flagged %>% filter(!is.na(L_summary), !is.na(L_asup))
if (nrow(both) > 0) {
  say("\nWhere both resolve (", nrow(both), " samples):")
  say("  identical                : ", sum(both$L_summary == both$L_asup))
  say("  differ                   : ", sum(both$L_summary != both$L_asup))
  say("  median |difference| (bp) : ", median(abs(both$L_summary - both$L_asup)))
  say("  max    |difference| (bp) : ", max(abs(both$L_summary - both$L_asup)))
}
say("\nL_summary distribution across the flagged set:"); print(summary(flagged$L_summary))
say("\nTriggering subtype frequency:")
print(flagged %>% count(denovo_minor_subtype, sort = TRUE) %>% as.data.frame())

## -- 4. 2k1b pairing rule (is_valid_minor rule 2) -------------------------
geno_of <- function(s) if_else(s == "2k1b", s, substr(s, 1, 1))
flagged <- flagged %>%
  mutate(ga = geno_of(denovo_minor_subtype), gb = geno_of(Major_subtype),
         pair_2k1b = coalesce((ga == "2k1b" & gb %in% c("1", "2", "2k1b")) |
                              (gb == "2k1b" & ga %in% c("1", "2", "2k1b")), FALSE)) %>%
  select(-ga, -gb)
say("\n2k1b/genotype-{1,2} pairs among the flagged set: ", sum(flagged$pair_2k1b))

## -- 5. Must-keep samples -------------------------------------------------
rule("4. Must-keep samples (the demoted legacy typable=YES minors)")
mk <- flagged %>% filter(sample_id %in% MUST_KEEP) %>%
  select(run, sample_id, Major_subtype, denovo_minor_subtype,
         L_summary, L_asup, pident, kmer_cov, aln_len, aln_frac_pct, pair_2k1b)
if (nrow(mk) < length(MUST_KEEP))
  die("only ", nrow(mk), " of ", length(MUST_KEEP), " must-keep samples are in the ",
      "flagged set (missing: ", paste(setdiff(MUST_KEEP, mk$sample_id), collapse = ", "),
      "). must_keep_kept would be meaningless — fix the cohort or the id keying first.")
print(as.data.frame(mk), row.names = FALSE)

## -- 6. The sweep ---------------------------------------------------------
n_cohort <- nrow(sm)
sweep_one <- function(floor_bp, len_col, excl_2k1b) {
  L <- flagged[[len_col]]
  keep <- coalesce(L, Inf) >= floor_bp     # NA fails OPEN, matching the pipeline gate
  if (excl_2k1b) keep <- keep & !flagged$pair_2k1b
  tibble(length_source = len_col, exclude_2k1b = excl_2k1b, floor_bp = floor_bp,
         flags_kept = sum(keep), suppressed = sum(!keep),
         pct_of_cohort = round(100 * sum(keep) / n_cohort, 1),
         must_keep_kept = paste0(sum(keep & flagged$sample_id %in% MUST_KEEP),
                                 "/", length(MUST_KEEP)),
         to_high = sum(!keep & flagged$n_sentences == 1))
}
sources <- c("L_summary", if (have_asup && any(!is.na(flagged$L_asup))) "L_asup")
sweep <- expand_grid(floor_bp = FLOORS, len_col = sources, excl = c(FALSE, TRUE)) %>%
  pmap_dfr(function(floor_bp, len_col, excl) sweep_one(floor_bp, len_col, excl)) %>%
  arrange(length_source, exclude_2k1b, floor_bp)

rule("5. Sweep")
say("Cohort ", n_cohort, ". to_high = samples that would move provisional -> high.\n")
print(as.data.frame(sweep), row.names = FALSE)

## -- 7. Per-floor disagreement between the two definitions ----------------
cross <- NULL
if (nrow(both) > 0) {
  rule("6. Where the two definitions disagree, by floor")
  cross <- map_dfr(FLOORS, function(f) {
    a <- both$L_summary >= f; b <- both$L_asup >= f
    tibble(floor_bp = f, agree = sum(a == b), disagree = sum(a != b),
           only_summary_keeps = sum(a & !b), only_asup_keeps = sum(b & !a))
  })
  print(as.data.frame(cross), row.names = FALSE)
}

## -- 8. Outputs -----------------------------------------------------------
per_sample <- flagged %>%
  mutate(decision_1000_2k1b = if_else(coalesce(L_summary, Inf) >= 1000 & !pair_2k1b,
                                      "KEEP", "SUPPRESS"),
         is_must_keep = sample_id %in% MUST_KEEP) %>%
  select(run, sample_id, sampleName, call_confidence, overall_sample_call,
         Major_subtype, denovo_minor_subtype, denovo_minor_ref,
         L_summary, L_asup, pident, kmer_cov, aln_len, aln_frac_pct,
         pair_2k1b, n_sentences, is_must_keep, decision_1000_2k1b, review_flag) %>%
  arrange(desc(L_summary))

write_csv(per_sample, file.path(outdir, "flagged_samples.csv"))
write_csv(sweep,      file.path(outdir, "sweep.csv"))
if (!is.null(cross)) write_csv(cross, file.path(outdir, "definition_disagreement.csv"))

rule("Done")
say("Wrote ", file.path(outdir, "flagged_samples.csv"), ", sweep.csv",
    if (!is.null(cross)) ", definition_disagreement.csv" else "")
