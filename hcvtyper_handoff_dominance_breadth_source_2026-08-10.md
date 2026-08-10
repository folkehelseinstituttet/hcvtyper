# Handoff — the dominance score's breadth term reads first-pass coverage, not the targeted mapping

**Date:** 2026-08-10
**Branch:** `dev` @ `42b4e38`
**Purpose:** report a scoring defect found while verifying Table S_DOM of the Access Microbiology
revision against a production run, and propose a minimal fix.
**Status:** **bug confirmed, reproduced against a real run; fix proposed, not applied.**

---

## 1. Summary

`score_candidates()` computes the dominance score's breadth term *B* from **`candidate_cov`, which is
the first-pass all-reference mapping coverage**, not from the targeted second-pass mapping. Every
other term in the same score is a second-pass quantity. The result is a score that mixes two
different mappings without saying so, and in the worst case credits a candidate with coverage it
does not have once reads are competitively assigned.

Two distinct manifestations, one root cause:

1. **Wrong source.** For a normal candidate, *B* is the first-pass breadth. For `ERR1810469`'s 3a
   candidate that is **46%**, where the targeted breadth is **18.2%** — the score awards 0.83 points
   of breadth the candidate has not earned.
2. **Silently zero.** For a **rescued or nominated** candidate, `rescue_evaluation.R` deliberately
   blanks `candidate_cov` (it described the *displaced* reference), so *B* falls to `NA` → `0`. The
   candidate scores as though **none of its reference were covered at all**, losing the full 3.0-point
   breadth contribution. In the Thomson 2016 run this hits **all five** of the de-novo-surfaced
   co-infections, whose real targeted breadths are 62–99%.

The intended behaviour appears to be the targeted breadth: `score_candidates()` already prefers a
column named `cand_cov_breadth` and only falls back to `candidate_cov`, and the comment beside the
fallback describes `candidate_cov` as "the targeted-mapping coverage percent" — which it is not.
**No module anywhere in the pipeline emits `cand_cov_breadth`**, so the fallback has always fired and
the preferred path has never executed in production.

I agree with the manuscript author's position that this should be the second-pass metric: the score
exists to rank candidates *after* competitive assignment, and every other input to it is already
post-competition.

---

## 2. Root cause

**`bin/classify_roles.R:389-401`** — the breadth source, inside `score_candidates()`:

```r
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
  breadth_frac <- ifelse(is.na(breadth_src), 0, ...)
```

Three facts make this misfire:

1. `cand_cov_breadth` is **never produced**. `grep -rn "cand_cov_breadth" bin/ modules/ workflows/`
   returns only the two references inside `classify_roles.R` itself (the doc comment at :366 and the
   branch at :393). The fallback is the only live path.
2. `candidate_cov` is **first-pass**. It originates in
   `bin/summarize_mapping_to_all_references.R:156` as `candidate_cov = percent_gt_4_int`, computed
   from the all-reference first mapping — the same measure as the `--minCov` entry criterion. It is
   carried unchanged through `rescue_evaluation.R` into `summary/candidates.csv`.
3. For rescued/nominated candidates it is **deliberately NA**. `bin/rescue_evaluation.R:549-553`
   blanks `candidate_reads` and `candidate_cov` after a replacement, correctly, because those numbers
   described the reference that was displaced; nominated candidates are created with `NA` from the
   start (:480-481). The comment there says "targeted mapping recomputes the real per-candidate
   numbers" — and it does, for reads (`targeted_reads_nodup`), but nothing ever recomputes breadth
   for the score.

The asymmetry with the neighbouring terms is what makes this a bug rather than a design choice.
In the same function, the reads term explicitly prefers the targeted count over the first-pass one
(`classify_roles.R:424-431`) for exactly the reason that applies here — "candidate_reads is the
neutral all-reference first-mapping count (with duplicates) ... targeted_reads_nodup ... reflects
true abundance". `cv_evenness` is computed from the targeted deduplicated depth file
(`bin/summarize.R:505-537`). Only breadth was left on the first-pass measurement.

---

## 3. Evidence — `ERR1810469` (Thomson 2016; the dev team has this run)

Run directory:
`/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon/HCVTyper_SRA_data/`

This is the pipeline's own worked dominance example — a 3a/1a sample reported as
`co-infection (indeterminate dominance)`.

| | 3a_D17763 (rank 1) | 1a_HQ850279 (rank 2) |
|---|---|---|
| `candidate_cov` (first-pass, **used by the score**) | **46%** | **69%** |
| targeted breadth ≥5× (`Major/Minor_cov_breadth_min_5`) | **18.17%** | **96.02%** |
| targeted breadth ≥10× | 4.18% | 77.73% |
| `targeted_reads_nodup` | 190 | 1,017 |
| `cv_evenness` | 0.4785 | 0.6922 |
| best contig k-mer cov | 51.65 | 5.40 |
| **`dominance_score` as shipped** | **5.9479** | **7.5571** |

The published score reconstructs exactly from the **first-pass** breadth, confirming which column was
read:

```
3a: 1.0*log10(190) + 3*0.46 + 3*0.47845 + 0.5*log10(1+50)
  = 2.2788 + 1.3800 + 1.4354 + 0.8538 = 5.94790   (candidates.csv: 5.947896923084096)
```

Substituting the targeted breadth gives 3a **5.11** (−0.83) and 1a **8.37** (+0.81) — a 1.6-point
swing in the gap. The ordering does not change here, but the table is wrong about why.

Note the direction of the error: the 3a candidate collapses from 5,009 first-pass reads to 190 after
competitive assignment and deduplication, and its targeted breadth collapses with it (46% → 18.2%).
The first-pass number is precisely the pre-competition figure the second-pass mapping exists to
correct, and it is the one the score used.

### The rescued/nominated case is worse

All five de-novo-surfaced co-infections in this run score with **`breadth_frac = 0`**:

| sample | candidate | route | `candidate_cov` | **actual targeted breadth ≥5×** | breadth points awarded |
|---|---|---|---|---|---|
| ERR1810475 | 1a_HQ850279 | nomination | NA | 99.3% | 0.00 (should be 2.98) |
| ERR1810447 | 2b_D10988 | rescue (2k1b rule) | NA | 95.9% | 0.00 (should be 2.88) |
| ERR1810453 | 2b_D10988 | rescue (four-floor) | NA | 94.9% | 0.00 (should be 2.85) |
| ERR1810449 | 4d_DQ418786 | nomination | NA | 74.5% | 0.00 (should be 2.24) |
| ERR1810485 | 3a_D17763 | nomination | NA | 62.2% | 0.00 (should be 1.87) |

These are the candidates the de novo layer exists to surface, and the dominance score is blind to the
single strongest piece of evidence that they are real. Verified for `ERR1810447`'s 2b: the shipped
`dominance_score` of 5.231 is only reproducible with the breadth term at zero.

No dominance ordering flips in this run: of 45 multi-candidate samples, the 15 in which every
candidate has a per-candidate depth file were recomputed with targeted breadth and the top-scoring
candidate was unchanged in all 15 (the other 30 have at least one candidate that was never
targeted-mapped, so no targeted breadth exists to substitute). The defect is currently latent for
call outcomes in this dataset — but it is a 3.0-point systematic penalty applied to exactly the
candidates the pipeline is least confident about, and a nominated strain that is genuinely dominant
would be ranked below a first-pass candidate that is not.

---

## 4. Suggested fix

The targeted breadth is already computed, per candidate, in the same loop that produces
`cv_evenness`, and there is already a per-candidate lookup joined into the scoring frame. Two small
edits.

**(a) `bin/summarize.R`, `cv_by_ref` (~:592-602)** — carry breadth alongside evenness and name it
what the classifier already asks for:

```r
cv_by_ref <- tmp_df %>%
  filter(reference != "first_mapping") %>%
  mutate(candidate_ref = str_remove(reference, "_cand[0-9]+$")) %>%
  select(sampleName, candidate_ref, cv_evenness,
         cand_cov_breadth = cov_breadth_min_5) %>%   # NEW: targeted breadth @>=5x, 0-100 percent
  filter(!is.na(candidate_ref)) %>%
  distinct(sampleName, candidate_ref, .keep_all = TRUE)
```

`cov_breadth_min_5` is the ≥5× breadth from the **deduplicated** targeted depth file, matching both
the `--minCov` semantics and the manuscript's definition of *B*. It is a 0-100 percent, which
`score_candidates()` already coerces (`classify_roles.R:399-400`). The existing `left_join(cv_by_ref,
...)` at `summarize.R:822` needs no change.

**(b) `bin/classify_roles.R:392-398`** — make the preference per-row rather than per-column, so a
candidate that was never targeted-mapped still falls back instead of silently scoring 0:

```r
  breadth_src <- if ("cand_cov_breadth" %in% names(df)) {
    if ("candidate_cov" %in% names(df)) {
      dplyr::coalesce(as.numeric(df$cand_cov_breadth), as.numeric(df$candidate_cov))
    } else {
      as.numeric(df$cand_cov_breadth)
    }
  } else if ("candidate_cov" %in% names(df)) {
    df$candidate_cov
  } else {
    rep(NA_real_, nrow(df))
  }
```

Without (b), adding the column in (a) would introduce a new NA path: any candidate with no depth file
would go from "first-pass breadth" straight to 0. Note the coalesce is *not* a full rescue of the
rescued/nominated case on its own — for those rows `candidate_cov` is NA too — which is why (a) is
the substantive fix and (b) is the guard.

Also fix the stale comment at `classify_roles.R:390`, which claims `candidate_cov` is "the
targeted-mapping coverage percent". That comment is how this survived review.

### Do not change these at the same time

`candidate_cov` has three other consumers, all of which legitimately want the **first-pass** value:

- `classify_roles.R:589-591` `clears_floor` — the informational `--minRead`/`--minCov` annotation.
- `classify_roles.R:616-618` `eligible` — the "has any coverage" sanity gate on dominance.
- `rescue_evaluation.R:239` `rescue_dominant_protect_cov` — the ≥90% first-pass protection guard,
  which must stay first-pass because it runs *before* the targeted mapping exists.

Scope the change to the score's breadth term only.

### Tests

- `bin/tests/test_dominance_score.R`: add a case with `cand_cov_breadth` present and `candidate_cov`
  divergent, asserting the score uses the former; and a case with `cand_cov_breadth = NA` +
  `candidate_cov` present, asserting the coalesce fallback.
- Add a regression case with `candidate_cov = NA` and `cand_cov_breadth = 95` — i.e. the rescued
  candidate — asserting a non-zero breadth contribution.
- The nf-test snapshot will move: every rescued/nominated candidate's `dominance_score` changes, and
  most others shift slightly. Expect a re-baseline, and check no `role` /
  `overall_sample_call` changes in the process (none expected from the 15-sample check above).

---

## 5. Manuscript implication

Table S_DOM in `Access_Microbiology/Revision_1/revised_manuscript_complete.md` publishes the *B* row
as "Coverage breadth (≥5×)" with 0.46 / 0.69, while the Results text for the same sample reports
18.2% and 77.7%. Both are correct measurements of different things, but as printed they read as a
contradiction, and a referee checking the supplementary table against the Results will find it.

Whichever way this is resolved, the two must agree: either the code is fixed and Table S_DOM is
regenerated with targeted breadth (preferred, and the author's stated preference), or the manuscript
must state explicitly that *B* is the first-pass breadth while *R* and *E* are second-pass. The
manuscript deadline is 2026-08-11, so if the fix cannot land and be re-run in time, the interim
action is the explicit wording in the manuscript plus this handoff as the record.
