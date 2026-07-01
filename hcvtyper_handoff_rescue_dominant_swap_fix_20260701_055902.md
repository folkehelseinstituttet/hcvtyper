# HCVTyper Handoff — De-novo Rescue Swaps Out the Dominant Strain + Rescue Traceability

- **Date:** 2026-07-01 05:59:02
- **Branch:** dev
- **Author:** jon.brate@fhi.no (analysis via Claude Code)
- **Status:** Bug diagnosed; fix **drafted, NOT yet applied** — implement via `/gsd-quick`.

> ## How to use this handoff (READ FIRST)
>
> **This handoff is self-contained. You do NOT need any results files, mounts, or
> `/mnt/N/...` paths to implement or test the fix.** Everything below works against
> the repository alone (dev branch). The observed run data that motivated the fix is
> embedded inline as worked examples, and all unit-test fixtures are synthetic and
> defined in-file.
>
> Files you will edit (all in the repo):
> - `bin/rescue_evaluation.R`
> - `bin/tests/test_rescue_evaluation.R`
> - `nextflow.config`
> - `conf/modules_hcv.config`
> - `modules/local/rescueevaluation/main.nf`
> - `bin/summarize.R`  *(traceability column; small, optional-to-defer)*
>
> Acceptance = the R contract test and the module nf-test pass (Section 4). Neither
> needs external data. A full end-to-end re-run on real data is *optional validation*
> for whoever has the analysis host — it is NOT required to land the fix.
>
> Appendix A (run-wide audit) is **reference/context only** and depends on results
> that live on the analysis host; skip it for implementation.

---

## 1. Problem summary (worked example, embedded — no files needed)

On a real HCV sample the de-novo "rescue" step **discarded the genuine dominant
strain** and replaced it with a thin cross-mapping artefact, and it carried the
displaced reference's read/coverage numbers onto the new label.

### 1a. The candidate table before vs after rescue (embedded)

**First mapping produced two candidates** (both genotype 1):

| rank | ref | subtype | geno | first-map reads | cov |
|---|---|---|---|---|---|
| 1 | `1a_HQ850279` | 1a | 1 | **6,534,184** | 100 |
| 2 | `1_KJ439780`  | 1  | 1 | 55,438 | 91 |

**Assembly support for that sample** (`assembly_support.csv`, embedded):

| subtype | best_ref | contig len | pident | aln | **kmer_cov** |
|---|---|---|---|---|---|
| **1a** (own) | `1a_HQ850279` | 2118 | 92.30 | 2118 | **31647×** |
| 1b | `1b_EU781827` | 1824 | 90.74 | 1804 | 1.26 |
| 3a (alt) | `3a_X76918` | 4189 | 92.36 | 4187 | **2.73×** |

**After rescue** the rank-1 slot was overwritten `1a → 3a`, keeping 1a's numbers:

| rank | ref | rescued_from | candidate_reads (STALE) | cov (STALE) | targeted_reads_nodup (real) | role |
|---|---|---|---|---|---|---|
| 1 | `3a_X76918`  | `1a_HQ850279` | 6,534,184 | 100 | **3,180** | co-infection |
| 2 | `1_KJ439780` | NA | 55,438 | 91 | **97,242** | dominant |

Net effect on the final call: the subtype-resolved `1a_HQ850279` Major was lost;
the Major degraded to the generic genotype-1 `1_KJ439780` (subtype "1"); a spurious
3a co-infection was invented from a 2.73× k-mer contig; and 3a displayed 1a's
"95.3% of mapped reads" (that figure is `6,534,184 / 6,853,211`, i.e. 1a's share).
GLUE on the surviving Major consensus still returns genotype 1 / subtype 1a, so the
clinical **genotype** was recoverable, but the reported reference/subtype and one
summary genotype cell were wrong.

---

## 2. Root cause (code-grounded — `bin/rescue_evaluation.R`)

Configured floors (`nextflow.config`): `rescue_min_length=3000`,
`rescue_min_pident=85`, `rescue_min_aln_length=3000`, `rescue_min_kmer_cov=2.0`,
`rescue_1a1b_length=5000`.

`evaluate_row()` (≈lines 190-263) decides rescues from **`assembly_support.csv`
only** — it never reads `candidate_reads`/`candidate_cov`, so mapping dominance is
invisible to the rescue. For the example above:

1. **Own-subtype confirmation FAILS on length** (≈lines 206-209): 1a's contig is
   2118 bp < 3000 → `floors_ok` FALSE → 1a treated as "unconfirmed", despite
   kmer_cov 31,647× (clearly the dominant strain).
2. **Best alternative = the long, thin 3a** (≈lines 215-218): 4189 bp.
3. **3a clears all four floors** (≈line 256): 4189>3000, 92.36>85, 4187≥3000,
   kmer_cov 2.73≥2.0. Rescue fires; rank-1 rewritten 1a→3a with 1a's reads/cov
   carried over (≈lines 317-322).

**Why 1a can't return:** additive nomination (≈lines 383-433) only surfaces
genotypes *not already present*; genotype 1 was still "present" via `1_KJ439780`, so
nothing re-nominated the proper 1a reference.

**The stale-stats bug:** on REPLACE, `candidate_ref` is overwritten but
`candidate_reads`/`candidate_cov` are not reset (≈lines 317-322), so they keep
describing the displaced reference.

**Key constraint discovered in the test suite** (`bin/tests/test_rescue_evaluation.R`):
Subtest 1 ("rescue-fires") deliberately encodes a cross-genotype REPLACE as
*correct* — a lone `1a` candidate with **no** `1a` assembly support, de novo says
`2b`, is *supposed* to be remapped `1a → 2b`. So a blanket "block all cross-genotype
replacement" would break intended behaviour. The real discriminator between that
(replace is right) and the example above (replace is wrong) is **whether the
candidate's own subtype is assembled**: subtest 1 has no `1a` contig; the example
has one (kmer 31,647×), merely fragmented below the length floor. All guards below
are scoped around that discriminator so they are **non-breaking** against the 13
existing subtests (guards are inert when the new args are absent or the own subtype
is unassembled).

---

## 3. Fixes to implement

### 3.0 Stale first-mapping stats — reset on REPLACE

`bin/rescue_evaluation.R`, insert immediately before the `out_cols <- ...` line
(≈485), i.e. after collapse/re-rank so ranking still uses real first-map reads:

```r
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
```

### 3.1 New parameters + wiring

`nextflow.config`, after `rescue_1a1b_length` (≈line 67):

```groovy
    rescue_kmer_cov_ratio       = 10        // Fix #3: block a REPLACE when the candidate's own-subtype contig k-mer coverage exceeds the rescue target's by >= this factor (cross-mapping-noise guard). 0/blank => disabled.
    rescue_dominant_protect_cov = 90        // Fix #2: a first-mapping candidate at >= this coverage whose own subtype IS assembled (identity+k-mer floors met) is self-confirming and never replaced. Blank => disabled.
```

`conf/modules_hcv.config`, RESCUE_EVALUATION `ext.args` (≈line 108) — append the
11th/12th positional args after `n_candidates`:

```groovy
        ext.args   = { "${params.rescue_min_length} ${params.rescue_min_pident} ${params.rescue_min_aln_length} ${params.rescue_min_kmer_cov} ${params.rescue_1a1b_length} ${params.n_candidates} ${params.rescue_kmer_cov_ratio} ${params.rescue_dominant_protect_cov}" }
```

`bin/rescue_evaluation.R`, after the `n_candidates_cap` block (≈line 83):

```r
# Optional 11th/12th positional args (Fix #2/#3). Absent => guard disabled, so
# the 9-arg subprocess test and any legacy caller preserve prior behaviour.
rescue_kmer_cov_ratio <- if (length(args) >= 11) suppressWarnings(as.numeric(args[11])) else NA_real_
if (is.na(rescue_kmer_cov_ratio) || rescue_kmer_cov_ratio <= 0) rescue_kmer_cov_ratio <- Inf
dominant_protect_cov  <- if (length(args) >= 12) suppressWarnings(as.numeric(args[12])) else NA_real_
if (is.na(dominant_protect_cov)) dominant_protect_cov <- Inf
```

### 3.2 Fix #2 (primary) — mapping-aware own-subtype confirmation

`bin/rescue_evaluation.R`, inside `evaluate_row()`, replace the block that currently
reads (≈lines 200-210):

```r
  sample_support <- support %>% filter(sample == sample_id)
  if (nrow(sample_support) == 0) return(no_rescue)

  # Own-subtype floor check ...
  if (!isTRUE(cand_sub == "2k1b")) {
    own_row <- sample_support %>% filter(subtype == cand_sub) %>%
      arrange(desc(best_contig_length)) %>% slice(1)
    if (nrow(own_row) > 0 && floors_ok(own_row, rescue_min_length)) return(no_rescue)
  }
```

with:

```r
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
```

Effect: the example's 1a (cov 100, own contig pident 92.3 / kmer 31,647×) is
confirmed → not replaced. Genotype 3 is then still absent, so the **existing
nomination path** re-surfaces 3a additively (NA reads), and collapse keeps
`{1a_HQ850279 (Major), 3a_X76918 (minor)}`. **1a preserved as Major.** This + the
existing nomination path also achieve the original "route different-genotype
evidence through nomination, not replace" goal without breaking subtest 1, so no
standalone cross-genotype block is added.

### 3.3 Fix #3 (backstop) — apply the k-mer guard at the rescue-return points

`kmer_ratio_blocks()` is defined in 3.2. Wire it into both rescues in
`evaluate_row()`:

2k1b path (≈line 235) — add to the condition:

```r
        if (!is.na(rescue_ref) && rescue_ref != orig_ref &&
            !kmer_ratio_blocks(g2_row$best_contig_kmer_cov[1])) {
```

Normal four-floor path (≈after line 256):

```r
  if (!floors_ok(alt_support, length_floor)) return(no_rescue)
  if (kmer_ratio_blocks(alt_support$best_contig_kmer_cov[1])) return(no_rescue)  # Fix #3
```

For the example: own 1a 31,647× vs target 3a 2.73× → ratio ≈ 11,600 ≥ 10 → blocked.
**REPLACE path only — never nomination** (a relative-kmer guard on nomination would
break subtest 12 / ERR1810475, which nominates a thin 7.1× contig under a 30,307×
dominant).

### 3.4 Fix #4 (NEW) — make the rescue traceable in the output files

**Goal:** the pipeline output must let a reader determine, per sample, whether the
de-novo rescue took effect — **including rescues that fired and were then discarded
by genotype-collapse or the `n_candidates` cap**, which currently leave no trace in
`candidates.csv` (the observability gap). Two parts:

#### 4a. Per-sample ledger: `{prefix}.rescue_audit.csv` (authoritative trace)

`bin/rescue_evaluation.R` — accumulate every rescue/nomination/block decision and
its final disposition, and ALWAYS write the ledger (header-only when nothing fired,
so downstream joins never miss a file).

Add near the top of the driver section (before the per-slot loop, ≈line 266):

```r
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
```

In the per-slot REPLACE loop (≈277-332), record each outcome. Where a rescue is
committed (inside the `if (... !rescue_would_collapse && !rescue_would_dup_genotype)`
block), add:

```r
      add_audit(out$sample[i], "replace", res$rescued_from, res$rescue_ref, res$rescue_trigger)
```

and in the `else`/guard branches, record the block (add explicit branches so both
guards are logged):

```r
    } else if (!is.na(res$rescue_ref) && rescue_would_collapse) {
      add_audit(out$sample[i], "blocked_collapse", out$candidate_ref[i], res$rescue_ref, res$rescue_trigger)
    } else if (!is.na(res$rescue_ref) && rescue_would_dup_genotype) {
      add_audit(out$sample[i], "blocked_dup_genotype", out$candidate_ref[i], res$rescue_ref, res$rescue_trigger)
    }
```

In the nomination loop (≈after the `bind_rows(out, tibble(...))` append, ≈line 431):

```r
      add_audit(sample_id, "nominate", NA_character_, nref, trig)
```

Capture the pre-cap survivor set so disposition can distinguish collapse vs cap.
Just before the cap `slice()` (≈line 465), record:

```r
  kept_refs_precap <- ord[keep_idx, , drop = FALSE]$candidate_ref
```

(Declare `kept_refs_precap <- character(0)` once in outer scope so it exists when no
candidates.) Finally, before writing `rescued.candidates.csv` (≈line 489), resolve
dispositions and write the ledger:

```r
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
```

`modules/local/rescueevaluation/main.nf` — declare the new output (script path) and
a stub:

```groovy
    // in output:
    tuple val(meta), path("*.rescue_audit.csv"), emit: rescue_audit
```

```bash
    # in stub: header-only ledger so -stub-run matches the emit glob
    printf "sample,event,original_ref,original_subtype,target_ref,target_subtype,target_genotype,evidence,disposition\n" > ${prefix}.rescue_audit.csv
```

Ensure the RESCUE_EVALUATION `publishDir` (in `conf/modules_hcv.config`) publishes
`*.rescue_audit.csv` alongside the candidates CSV (directory publish or add the glob).

#### 4b. Top-level flag: `rescue_effect` column in `Summary.csv`

`bin/summarize.R` already ingests `rescued_from`/`rescue_trigger` from the candidates
CSVs and rolls up `rescue_flag`. Add a derived per-sample `rescue_effect` from data
it already has (no new channel input — avoids the SUMMARIZE channel-arity fragility):
after roles are assigned, set

- `major_ref_changed` — the candidate holding the dominant/Major role has
  `rescued_from` non-NA;
- else `minor_ref_changed` — any surviving candidate has `rescued_from` non-NA;
- else `none`.

Document in the column/comment that `rescue_effect` reflects only rescues that
**survived** into the final call; the standalone `rescue_audit.csv` (4a) is the
authoritative record that additionally captures `dropped_collapse` / `dropped_cap`
events. (Optional, phase-2: also glob the per-sample `*.rescue_audit.csv` in
summarize to emit `n_rescue_events` / `n_rescue_dropped`; defer unless needed, per
the SUMMARIZE arity caution.)

### 3.5 New unit tests (`bin/tests/test_rescue_evaluation.R`, self-contained)

Add an `extra_args` param to `run_rescue()` so a test can enable the new guards
(guards are inert unless args 10-12 are supplied):

```r
run_rescue <- function(case, prefix, cands, support, refs, prestage = character(0), extra_args = NULL) {
  ...
  exit <- system2("Rscript",
    c(shQuote(script), shQuote(prefix), shQuote(cand_path), shQuote(support_path),
      shQuote(refs_path),
      TH_LENGTH, TH_PIDENT, TH_ALN, TH_KMER, TH_1A1B_LEN,
      if (is.null(extra_args)) character(0) else extra_args),
    stdout = FALSE, stderr = FALSE)
  ...
}
```

Then append:

```r
REFS_ERR <- c(REFS, list("1a_HQ850279" = strrep("G",60),
                         "3a_X76918"   = strrep("T",60),
                         "1_KJ439780"  = strrep("A",60)))

# Subtest 14: dominant strain protected (real-sample shape).
rA <- run_rescue("err_dominant", "ERRDOM",
  cands = mk_cands("ERRDOM",
    mk_cand(1, "1a_HQ850279", "1a", "1", 6534184, 100, "pass"),
    mk_cand(2, "1_KJ439780",  "1",  "1",   55438,  91, "pass")),
  support = mk_support("ERRDOM",
    mk_support_row("1a", "1a_HQ850279", 2118, 92.304, 2118, 31647.19),
    mk_support_row("3a", "3a_X76918",   4189, 92.357, 4187, 2.728)),
  refs = REFS_ERR, extra_args = c("2", "10", "90"))
assert_schema("err_dominant", rA$cands_out)
a1 <- rA$cands_out %>% filter(candidate_ref == "1a_HQ850279")
if (nrow(a1) != 1 || !is.na(a1$rescued_from[1]))
  fail("err-dominant: dominant 1a must survive UNREPLACED (Fix #2)")
if (!any(rA$cands_out$candidate_ref == "3a_X76918"))
  fail("err-dominant: 3a must be surfaced additively via nomination")
if (any(rA$cands_out$candidate_ref == "1_KJ439780"))
  fail("err-dominant: redundant genotype-1 slot must collapse")
if ((rA$cands_out %>% arrange(candidate_rank) %>% slice(1))$candidate_ref[1] != "1a_HQ850279")
  fail("err-dominant: Major (rank 1) must remain 1a_HQ850279")
ok("dominant-protect -> well-covered assembled 1a kept as Major; thin 3a additive")

# Subtest 15: Fix #3 in isolation (below the #2 cov floor).
rB <- run_rescue("kmer_guard", "KMERG",
  cands = mk_cands("KMERG",
    mk_cand(1, "1a_HQ850279", "1a", "1", 500000, 50, "pass")),
  support = mk_support("KMERG",
    mk_support_row("1a", "1a_HQ850279", 2118, 92.3,   2118, 31647.19),
    mk_support_row("3a", "3a_X76918",   4189, 92.357, 4187, 2.728)),
  refs = REFS_ERR, extra_args = c("2", "10", "90"))
b1 <- rB$cands_out %>% filter(candidate_ref == "1a_HQ850279")
if (nrow(b1) != 1 || !is.na(b1$rescued_from[1]))
  fail("kmer-guard: 1a REPLACE must be blocked by the relative k-mer guard (Fix #3)")
ok("kmer-guard -> own-subtype depth >> target depth blocks the REPLACE (Fix #3)")

# Subtest 16: legitimate correction still fires AND stale stats are blanked.
rC <- run_rescue("stale_reads", "STALEREADS",
  cands = mk_cands("STALEREADS",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("STALEREADS",
    mk_support_row("2b", "2b_AY232748", 4120, 96.2, 3840, 4.1)),
  refs = REFS, extra_args = c("2", "10", "90"))
c1 <- rC$cands_out %>% filter(candidate_rank == 1)
if (is.na(c1$rescued_from[1]))
  fail("stale-reads: legitimate 1a->2b correction must still fire with guards enabled")
if (!is.na(c1$candidate_reads[1]) || !is.na(c1$candidate_cov[1]))
  fail("stale-reads: a REPLACED candidate must not carry the displaced ref's reads/cov")
ok("stale-reads -> replace still fires; displaced ref's first-mapping reads/cov blanked")

# Subtest 17: rescue audit ledger is written and records the replace.
audit_path <- file.path(dirname(rC$path), "STALEREADS.rescue_audit.csv")
if (!file.exists(audit_path))
  fail("rescue-audit: {prefix}.rescue_audit.csv must be written")
aud <- read_csv(audit_path, show_col_types = FALSE)
if (!any(aud$event == "replace" & aud$target_ref == "2b_AY232748"))
  fail("rescue-audit: the 1a->2b replace must appear in the ledger")
if (!any(grepl("^retained_rank_", aud$disposition)))
  fail("rescue-audit: a surviving replace must have a retained_rank_* disposition")
ok("rescue-audit -> ledger records the replace with its final disposition (Fix #4)")
```

> Note: subtest 17 assumes `run_rescue` runs in a tempdir it can read back from; the
> harness already exposes `$path` (the rescued CSV) — derive the audit path from
> `dirname($path)` as shown, or extend `run_rescue` to also return the audit path.

---

## 4. Acceptance criteria (all repo-local, no external data)

1. `Rscript bin/tests/test_rescue_evaluation.R` prints `ALL PASS` (existing 13
   subtests unchanged + new 14-17).
2. Module nf-test passes:
   `modules/local/rescueevaluation/tests/main.nf.test`
   (uses `~/.nf-test` 0.9.3, `--profile docker`; update the snapshot for the new
   `rescue_audit` emit). **Watch host disk** — nf-test/docker can fill it (see the
   project note on near-full disk); `sudo rm -rf work/ && docker volume prune` if it
   ENOSPCs.
3. `nextflow config -profile test` resolves the two new params without error.
4. Behavioural expectation encoded by the tests:

   | | Before | After |
   |---|---|---|
   | Example candidate set | {3a_X76918, 1_KJ439780} | {**1a_HQ850279** (Major), 3a_X76918 (minor)} |
   | Major subtype | "1" (unresolved) | **1a** |
   | 3a's reads/cov | 1a's stale 6.5M / 100 | **NA** (never first-mapped) |
   | rescue traceability | surviving rescues only | full ledger incl. dropped_collapse/dropped_cap |
   | existing subtests 1-13 | pass | **pass** |

---

## 5. Notes / open questions

- `rescue_dominant_protect_cov = 90` and `rescue_kmer_cov_ratio = 10` are starting
  values; tune against the benchmark panel so the milestone's closed failure modes
  don't reopen. Both default to "disabled" when the args are absent, so nothing
  changes for callers that don't pass them.
- The spurious 3a **minor** nomination is a separate calibration question (same thin-
  contig shape subtest 12 legitimately keeps) and is what already raises the
  `review_flag`; out of scope here.
- Distinct summarize-side bug seen on the example: `Summary.csv` `Major_genotype`
  read `3` while `Major_reference`/GLUE were genotype 1 — the genotype/subtype summary
  cells followed candidate RANK (stale-reads-inflated) while the reference cells
  followed ROLE. The fixes above (1a kept as real rank-1 Major) make this cell
  compute as 1; if it persists after the fix, treat as a separate summarize ticket.

---

## Appendix A — Run-wide audit evidence (REFERENCE ONLY; needs the analysis host)

> This appendix documents *why* the fix is low-risk and *how broad* the impact was.
> It was derived from result files on the analysis host and is **NOT required to
> implement or test the fix**. Do not block on reproducing it.

Across three revision runs (68 samples, 81 candidate rows), the REPLACE rescue fired
only **twice**, and only one touched a Major:

| Sample | REPLACE hit | Original Major | Post-fix-era Major | Major affected? |
|---|---|---|---|---|
| 2672173 | Minor slot (`1c_AY651061 → 3a_D17763`) | 1a_HQ850279 (1a) | 1a_HQ850279 (1a) | No — minor was *corrected* to match the original |
| 2679090 | Dominant slot (`1a_HQ850279 → 3a_X76918`) | 1a_HQ850279 (1a) | 1_KJ439780 ("1") | Yes — subtype de-resolved (genotype 1 preserved) |

- **0** dominant rows had `NA` first-mapping reads → no *nominated* candidate ever
  became a Major (nominations only land as minors).
- Stale reads/cov existed only on those 2 rescued rows and were display-only; role
  classification uses real targeted reads, so no Major flipped from stale reads.

**"Fire-then-collapse" caveat** — `candidates.csv` shows only survivors of
REPLACE → nominate → collapse → cap, so a rescue that fired then was collapsed/capped
leaves no trace (this is exactly what Fix #4 closes). Its bounded impact:
- **Major calls: provably immune.** Major = top read-recruiter (rank 1); collapse
  only drops the lower recruiter of a same-genotype pair, the cap only drops beyond
  rank N — neither can remove rank 1.
- **Genotype calls: bounded, empirically zero here.** Collapse only drops a
  genotype already represented by a survivor; the only way a *distinct* genotype is
  lost is the `n_candidates=2` cap dropping a genuine 3rd genotype. Evidence: max **2**
  distinct floor-passing genotypes per sample (8 samples had 2, none had 3+) → the
  cap never dropped a real genotype.
- **Performance:** the real residual risk is a false-negative *minor* (same-genotype
  co-infection — a known/deliberate limitation) and, mainly, the **observability gap**
  (reported rescue count is a lower bound) — which is why Fix #4 exists.

**Hygiene finding (host-side):** stale, old-schema `sim*` outputs (no `best_ref`
column) were co-mingled with a current run's `blastparse/` dir, and the
`assembly_support.csv` schema drifted mid-project (`best_ref` added). Glob-based
aggregation can pick up stale files. Worth a cleanup + a schema-version guard, but
unrelated to the code fix here.
