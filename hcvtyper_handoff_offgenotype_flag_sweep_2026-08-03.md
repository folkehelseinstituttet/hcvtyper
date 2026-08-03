# Handoff: re-derive the off-genotype-contig review-flag sweep

**Date:** 2026-08-03
**Author:** analysis prepared against hcvtyper dev @ `28a568d` (the build that produced the results)
**Run on:** the machine where `HCV/2026/HCV_paper_revisjon/` and the per-run `summary/` + `blastparse/` directories are mounted.
**Access needed:** read-only. This procedure writes nothing except its own output directory.

---

## 1. Why this needs re-deriving

The comparison report `hcvtyper_v1.3.0-g28a568d_vs_v1.1.x_comparison.html` (§5.3) established that the
review sentence

> "Monoinfection called for candidate N, but de novo assembly found a different-genotype contig (X) — possible missed co-infection or contamination. Please review."

fires on **72 of 140 samples (51%)** with no contig-length floor, and proposed gating it on contig length.
It measured each triggering contig's length from `blastparse/<sample>.assembly_support.csv`, matched on
the subtype named in the flag — i.e. the column **`best_contig_length`**, which is *the longest contig
per subtype*.

**That is not the number the fix would gate on.** The value already available in `Summary.csv` is
**`denovo_minor_contig_length`**, produced at `bin/blast_parse.R:328-332` as the longest contig hitting
the chosen *reference*, excluding the major contig. The two are related but not identical:

| | `best_contig_length` (what the report measured) | `denovo_minor_contig_length` (what the gate would use) |
|---|---|---|
| Grain | per **subtype** | per **reference** (`denovo_minor_ref`) |
| Source | `blastparse/*.assembly_support.csv` | `Summary.csv`, already joined at `bin/summarize.R:1301` |
| Excludes the major contig? | no | yes (`filter(qseqid != major_contig)`) |
| Selection | `slice_max(sc_length)` within subtype | `slice_max(sc_length)` within reference |

There is a further wrinkle: `denovo_minor_ref` itself is chosen by **BLAST rank**, not by contig size
(`bin/blast_parse.R:308-312`, `slice_head(n = 1)` over the different-genotype hits). So the reference that
names the flag and the contig whose length is reported can, in principle, be different contigs. Any
threshold picked from the report's numbers may therefore land differently once applied to the real column.

**Deliverable:** the same sweep, computed on `denovo_minor_contig_length`, with the report's
`best_contig_length` arm reproduced alongside it for comparison, so a defensible default can be locked in.

---

## 2. Inputs

Five runs, 140 samples total, all under the new-build results root:

```
HCV/2026/HCV_paper_revisjon/
  20251113-02/   (n=23)
  20251212-01/   (n=22)
  20260326-01/   (n=28)
  20260521-01/   (n=36)
  20260625-02/   (n=31)
```

Per run the script needs:

- `<run>/summary/Summary.csv` — required.
- `<run>/blastparse/*.assembly_support.csv` — optional; only the comparison arm needs it. If absent, the
  script still produces the primary result and reports the comparison arm as unavailable.

Legacy v1.1.x results are **not** needed. This is entirely a within-v1.3.0 analysis.

### Columns consumed from `Summary.csv`

`sampleName`, `overall_sample_call`, `call_confidence`, `review_flag`, `Major_subtype` (or `Major` as
fallback — they are identical), `denovo_minor_subtype`, `denovo_minor_ref`, `denovo_minor_contig_length`.

> **Note.** `Major_genotype` / `Minor_genotype` are missing entirely from runs 20251113-02 and
> 20260625-02 (§6 of the report, a separate known defect). Nothing here reads them, and the script reads
> every column as character before coercing, so the ragged schema across runs cannot break the load.

### Columns consumed from `assembly_support.csv`

`sample`, `subtype`, `best_contig_length`, `best_contig_pident`, `best_contig_kmer_cov`.

---

## 3. What the script computes

**Step 1 — reconstruct the trigger and validate against reality.**
The predicate is transcribed verbatim from `bin/classify_roles.R:962-966`:

```r
overall_sample_call == "monoinfection" &
  !is.na(denovo_minor_subtype) & !is.na(Major_subtype) &
  substr(denovo_minor_subtype, 1, 1) != substr(Major_subtype, 1, 1)
```

It is cross-checked against which rows *actually* carry the sentence in `review_flag`.
**If those two sets disagree, stop and report** — every number downstream is built on the reconstruction
being exact. The text-matched set is treated as authoritative for the sweep.

Expected: **72 samples**, and on every one of them it is the only sentence in `review_flag`. Both facts
are re-verified rather than assumed.

**Step 2 — attach both length definitions** to those samples: `denovo_minor_contig_length` from
`Summary.csv`, and `best_contig_length` from `assembly_support.csv` joined on
`(run, sample, subtype = denovo_minor_subtype)`.

**Step 3 — quantify how far apart the two definitions are.** Exact agreement count, median absolute
difference, and — the number that actually matters — **how many samples land on opposite sides of each
candidate floor** depending on which definition is used.

**Step 4 — the 2k1b suppression rule**, mirroring `is_valid_minor()` rule 2 at
`bin/classify_roles.R:200-203`, using the 2k1b-aware genotype rule from `bin/genotype_utils.R:24`:

```r
geno_of <- function(s) ifelse(s == "2k1b", "2k1b", substr(s, 1, 1))
# suppress when one side is 2k1b and the other is genotype 1, 2, or 2k1b
```

Note this is deliberately *not* a wholesale swap to `is_valid_minor()`: that function returns `TRUE` for
1a/1b pairs, which would newly fire the trigger on the within-genotype-1 artefact class the build already
eliminated correctly.

**Step 5 — the sweep.** Floors `0, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000` bp × both length
definitions × {2k1b exclusion off, on}. For each cell: flags kept, flags suppressed, % of the 140-sample
cohort, and whether all three must-keep samples survive.

**Must-keep set** — the three legacy-typable=YES minors that v1.3.0 demotes to monoinfection (§4.1 of the
report). These are the cases that most deserve a flag; any floor that drops one is disqualified:

| Sample | Run | Report's `best_contig_length` | pident | k-mer cov |
|---|---|---|---|---|
| 2743986 | 20260521-01 | 4,467 bp (1b) | 91.9 | 1.42 |
| 2726018 | 20260326-01 | 2,787 bp (1a) | 93.3 | 1.97 |
| 2714375 | 20260326-01 | 2,706 bp (1b) | 92.6 | 1.62 |

The script prints their `denovo_minor_contig_length` values explicitly — **if these three shift materially
under the new definition, the recommended floor has to move with them.**

**Step 6 — confidence-tier impact.** For each configuration, how many samples would move
`provisional → high` (the trigger is their only `review_flag` sentence, so suppressing it empties the
field, and `!is.na(review_flag)` at `bin/summarize.R:1863` is the only thing currently holding them at
`provisional`).

---

## 4. The script

Requires R with `tidyverse` — any of the pipeline's own R containers will do (e.g.
`community.wave.seqera.io/library/r-seqinr_r-tidyverse`). Save as `offgeno_flag_sweep.R`, then:

```bash
Rscript offgeno_flag_sweep.R /path/to/HCV/2026/HCV_paper_revisjon ./offgeno_sweep_out
```

Both arguments are optional (default: current directory, `./offgeno_sweep_out`). The root can be any
directory containing the five run folders; the script finds `summary/Summary.csv` recursively.

```r
#!/usr/bin/env Rscript
# offgeno_flag_sweep.R -----------------------------------------------------
# Re-derive the off-genotype-contig review-flag sweep on the column the fix
# would actually gate on (denovo_minor_contig_length), with the comparison
# report's assembly_support best_contig_length arm reproduced alongside.
# READ-ONLY: touches nothing but its own output directory.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
root   <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
outdir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "offgeno_sweep_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(t) say("\n", strrep("=", 72), "\n", t, "\n", strrep("=", 72))

MUST_KEEP <- c("2743986", "2726018", "2714375")
FLOORS    <- c(0, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000)
TRIG_TEXT <- "different-genotype contig"

## -- 1. Load every Summary.csv --------------------------------------------
# All-character read: the Major_genotype/Minor_genotype columns are absent from
# 2 of the 5 runs and per-file type inference disagrees across runs, either of
# which would abort a typed bind_rows. Coerce the few numerics explicitly after.
sum_files <- list.files(root, pattern = "^Summary\\.csv$",
                        recursive = TRUE, full.names = TRUE)
sum_files <- sum_files[basename(dirname(sum_files)) == "summary"]
if (length(sum_files) == 0) stop("No summary/Summary.csv found under: ", root)

read_one <- function(f) {
  read_csv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
    mutate(run = basename(dirname(dirname(f))), .before = 1)
}
sm <- map_dfr(sum_files, read_one)

rule("1. Cohort")
say("Summary.csv files found : ", length(sum_files))
say("Runs                    : ", paste(sort(unique(sm$run)), collapse = ", "))
say("Rows (samples)          : ", nrow(sm), "   [report: 140]")

# Major_subtype and Major are verbatim copies; accept either.
if (!"Major_subtype" %in% names(sm) && "Major" %in% names(sm)) {
  sm <- mutate(sm, Major_subtype = Major)
  say("NOTE: Major_subtype absent; using Major (identical by construction).")
}
need <- c("sampleName", "overall_sample_call", "call_confidence", "review_flag",
          "Major_subtype", "denovo_minor_subtype", "denovo_minor_ref",
          "denovo_minor_contig_length")
missing <- setdiff(need, names(sm))
if (length(missing)) stop("Summary.csv is missing required columns: ",
                          paste(missing, collapse = ", "))

sm <- sm %>% mutate(L_summary = suppressWarnings(as.numeric(denovo_minor_contig_length)))

## -- 2. Reconstruct the trigger, and check it against what actually fired --
geno1 <- function(x) substr(x, 1, 1)

sm <- sm %>%
  mutate(
    pred_fires = !is.na(overall_sample_call) & overall_sample_call == "monoinfection" &
                 !is.na(denovo_minor_subtype) & !is.na(Major_subtype) &
                 geno1(denovo_minor_subtype) != geno1(Major_subtype),
    text_fires = !is.na(review_flag) & str_detect(review_flag, fixed(TRIG_TEXT)),
    n_sentences = if_else(is.na(review_flag), 0L, str_count(review_flag, fixed(" | ")) + 1L)
  )

rule("2. Trigger reconstruction")
say("Predicate (classify_roles.R:962-966) fires on : ", sum(sm$pred_fires))
say("review_flag text actually contains sentence   : ", sum(sm$text_fires), "   [report: 72]")
disagree <- sm %>% filter(pred_fires != text_fires)
if (nrow(disagree) > 0) {
  say("\n*** WARNING: predicate and text disagree on ", nrow(disagree), " sample(s). ***")
  say("*** The reconstruction is not exact — inspect before trusting the sweep. ***")
  print(disagree %>% select(run, sampleName, overall_sample_call, Major_subtype,
                            denovo_minor_subtype, pred_fires, text_fires))
  write_csv(disagree, file.path(outdir, "predicate_text_disagreement.csv"))
} else {
  say("Predicate and text agree on all ", nrow(sm), " samples. Reconstruction is exact.")
}

flagged <- sm %>% filter(text_fires)
say("\nIs it the sole sentence on every flagged sample? ",
    if (all(flagged$n_sentences == 1)) "YES (as reported)" else
      paste0("NO — ", sum(flagged$n_sentences > 1), " carry additional sentences"))

## -- 3. The comparison arm: assembly_support best_contig_length -----------
asup_files <- list.files(root, pattern = "\\.assembly_support\\.csv$",
                         recursive = TRUE, full.names = TRUE)
have_asup <- length(asup_files) > 0

if (have_asup) {
  asup <- map_dfr(asup_files, function(f) {
    read_csv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
      mutate(run = basename(dirname(dirname(f))), .before = 1)
  })
  # Summary keys on sampleName; assembly_support keys on the BLASTPARSE prefix
  # (e.g. "2714375-HCV"). Strip at "-HCV", matching how the report paired them.
  asup <- asup %>%
    mutate(sampleName = str_remove(sample, "-HCV.*$"),
           L_asup     = suppressWarnings(as.numeric(best_contig_length)),
           pident     = suppressWarnings(as.numeric(best_contig_pident)),
           kmer_cov   = suppressWarnings(as.numeric(best_contig_kmer_cov))) %>%
    group_by(run, sampleName, subtype) %>%
    slice_max(L_asup, n = 1, with_ties = FALSE) %>%   # guard against dup rows
    ungroup() %>%
    select(run, sampleName, subtype, L_asup, pident, kmer_cov)

  flagged <- flagged %>%
    left_join(asup, by = c("run", "sampleName", "denovo_minor_subtype" = "subtype"))
} else {
  flagged <- flagged %>% mutate(L_asup = NA_real_, pident = NA_real_, kmer_cov = NA_real_)
}

rule("3. The two length definitions")
say("assembly_support.csv files found: ", length(asup_files),
    if (!have_asup) "  -> comparison arm UNAVAILABLE" else "")
say("Flagged samples with L_summary  : ", sum(!is.na(flagged$L_summary)), " / ", nrow(flagged))
say("Flagged samples with L_asup     : ", sum(!is.na(flagged$L_asup)),    " / ", nrow(flagged))

both <- flagged %>% filter(!is.na(L_summary), !is.na(L_asup))
if (nrow(both) > 0) {
  say("\nOn the ", nrow(both), " samples where both are resolvable:")
  say("  identical                     : ", sum(both$L_summary == both$L_asup))
  say("  differ                        : ", sum(both$L_summary != both$L_asup))
  say("  median |difference| (bp)      : ", median(abs(both$L_summary - both$L_asup)))
  say("  max    |difference| (bp)      : ", max(abs(both$L_summary - both$L_asup)))
}

say("\nDistribution of L_summary across the flagged set:")
print(summary(flagged$L_summary))
say("\nTriggering subtype frequency:")
print(flagged %>% count(denovo_minor_subtype, sort = TRUE) %>% as.data.frame())

## -- 4. 2k1b pairing rule (is_valid_minor rule 2) -------------------------
geno_of <- function(s) if_else(s == "2k1b", s, substr(s, 1, 1))
is_2k1b_pair <- function(a, b) {
  ga <- geno_of(a); gb <- geno_of(b)
  coalesce((ga == "2k1b" & gb %in% c("1", "2", "2k1b")) |
           (gb == "2k1b" & ga %in% c("1", "2", "2k1b")), FALSE)
}
flagged <- flagged %>%
  mutate(pair_2k1b = is_2k1b_pair(denovo_minor_subtype, Major_subtype))

say("\n2k1b/genotype-{1,2} pairs among the flagged set: ", sum(flagged$pair_2k1b),
    "   [report: 13]")

## -- 5. Must-keep samples -------------------------------------------------
rule("4. Must-keep samples (the three demoted legacy typable=YES minors)")
mk <- flagged %>% filter(sampleName %in% MUST_KEEP) %>%
  select(run, sampleName, Major_subtype, denovo_minor_subtype,
         L_summary, L_asup, pident, kmer_cov, pair_2k1b)
if (nrow(mk) < length(MUST_KEEP)) {
  say("*** WARNING: only ", nrow(mk), " of ", length(MUST_KEEP),
      " must-keep samples are in the flagged set. ***")
  say("*** Missing: ", paste(setdiff(MUST_KEEP, mk$sampleName), collapse = ", "), " ***")
}
print(as.data.frame(mk))
say("\nReport's best_contig_length for these: 2743986=4467, 2726018=2787, 2714375=2706")
say("If L_summary differs materially from those, the recommended floor must move with it.")

## -- 6. The sweep ---------------------------------------------------------
n_cohort <- nrow(sm)
sweep_one <- function(floor_bp, len_col, excl_2k1b) {
  L <- flagged[[len_col]]
  keep <- coalesce(L, 0) >= floor_bp
  if (excl_2k1b) keep <- keep & !flagged$pair_2k1b
  tibble(
    length_source  = len_col,
    exclude_2k1b   = excl_2k1b,
    floor_bp       = floor_bp,
    flags_kept     = sum(keep),
    suppressed     = sum(!keep),
    pct_of_cohort  = round(100 * sum(keep) / n_cohort, 1),
    must_keep_kept = paste0(sum(keep & flagged$sampleName %in% MUST_KEEP),
                            "/", length(MUST_KEEP)),
    # sole-trigger samples that lose their only sentence -> would become `high`
    to_high        = sum(!keep & flagged$n_sentences == 1)
  )
}

sources <- c("L_summary", if (have_asup && any(!is.na(flagged$L_asup))) "L_asup")
sweep <- expand_grid(floor_bp = FLOORS, len_col = sources, excl = c(FALSE, TRUE)) %>%
  pmap_dfr(function(floor_bp, len_col, excl) sweep_one(floor_bp, len_col, excl)) %>%
  arrange(length_source, exclude_2k1b, floor_bp)

rule("5. Sweep")
say("Cohort size: ", n_cohort, ". `to_high` = samples that would move provisional -> high.\n")
print(as.data.frame(sweep))

## -- 7. Side-by-side disagreement at each floor ---------------------------
if (nrow(both) > 0) {
  rule("6. Where the two definitions disagree, by floor")
  cross <- map_dfr(FLOORS, function(f) {
    a <- both$L_summary >= f; b <- both$L_asup >= f
    tibble(floor_bp = f, agree = sum(a == b), disagree = sum(a != b),
           only_summary_keeps = sum(a & !b), only_asup_keeps = sum(b & !a))
  })
  print(as.data.frame(cross))
  say("\n'disagree' > 0 means the floor chosen from the report's numbers would")
  say("behave differently once wired to denovo_minor_contig_length.")
}

## -- 8. Outputs -----------------------------------------------------------
per_sample <- flagged %>%
  mutate(
    decision_1000_no2k1b = if_else(coalesce(L_summary, 0) >= 1000, "KEEP", "SUPPRESS"),
    decision_1000_2k1b   = if_else(coalesce(L_summary, 0) >= 1000 & !pair_2k1b,
                                   "KEEP", "SUPPRESS"),
    decision_2500_2k1b   = if_else(coalesce(L_summary, 0) >= 2500 & !pair_2k1b,
                                   "KEEP", "SUPPRESS"),
    is_must_keep = sampleName %in% MUST_KEEP
  ) %>%
  select(run, sampleName, call_confidence, overall_sample_call,
         Major_subtype, denovo_minor_subtype, denovo_minor_ref,
         L_summary, L_asup, pident, kmer_cov, pair_2k1b, n_sentences, is_must_keep,
         starts_with("decision_"), review_flag) %>%
  arrange(desc(L_summary))

write_csv(per_sample, file.path(outdir, "flagged_samples.csv"))
write_csv(sweep,      file.path(outdir, "sweep.csv"))
if (nrow(both) > 0) write_csv(cross, file.path(outdir, "definition_disagreement.csv"))

rule("Done")
say("Wrote:")
say("  ", file.path(outdir, "flagged_samples.csv"), "  (one row per flagged sample, both lengths)")
say("  ", file.path(outdir, "sweep.csv"))
if (nrow(both) > 0) say("  ", file.path(outdir, "definition_disagreement.csv"))
```

---

## 5. What to send back

The three CSVs, plus the console output (`Rscript ... 2>&1 | tee offgeno_sweep_out/console.log`).

The specific questions the output needs to answer:

1. **Does the reconstruction hold?** Predicate count == text count == 72. If not, everything else is
   suspect and the predicate needs re-deriving before any threshold is chosen.
2. **How far apart are the two length definitions?** Specifically the `definition_disagreement.csv` row
   at 1000 bp — if `disagree` is 0 or 1, the report's sweep transfers directly and the 1,000 bp
   recommendation stands as-is. If it is large, the floor has to be picked from `L_summary` alone.
3. **Where do the three must-keep samples sit under `L_summary`?** They are at 4,467 / 2,787 / 2,706
   under the report's definition. The chosen floor needs a comfortable margin below the smallest of them —
   this is the constraint that argues for 1,000 bp over the report's 2,500 bp.
4. **`flags_kept` and `to_high` at floor = 1000, `exclude_2k1b = TRUE`, `length_source = L_summary`.**
   That is the candidate default. The expectation from the report's numbers is roughly 15–21 flags kept
   (~11–15% of the cohort) and ~50–57 samples moving to `high`.

If `to_high` is uncomfortably large, that is the signal to take the two-tier option — emit a short
informational token on the sub-threshold contigs into `Major_evidence` instead of dropping them entirely,
so `call_confidence` stays clean without losing the observation.

---

## 6. Scope note

This procedure settles the **threshold only**. Two things it deliberately does not cover:

- The 2k1b exclusion is evaluated here but is not really a tuning question — the codebase already holds
  that position in `is_valid_minor()` rule 2 (`bin/classify_roles.R:200-203`) and the review trigger simply
  fails to consult it, using a hand-rolled `substr(x, 1, 1)` instead of the 2k1b-aware
  `genotype_from_subtype()`. That is a consistency fix regardless of what the sweep says.
- The `Major_genotype` / `Minor_genotype` defect (§6 of the report — columns absent from 2 of 5 runs, and
  swapped against the subtype columns on Sample51K) is unrelated and needs fixing on its own.
