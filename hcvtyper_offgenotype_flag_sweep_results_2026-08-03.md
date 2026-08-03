# Results: off-genotype-contig review-flag sweep

**Responds to:** `hcvtyper_handoff_offgenotype_flag_sweep_2026-08-03.md` (commit `c577b5e`)
**Date:** 2026-08-03
**Run environment:** conda env `/home/jon.brate/.conda/R_shared` (R + tidyverse)
**Root:** `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon`
**Read-only.** Nothing was written to `/mnt/N`.

**Deliverables** in `offgeno_sweep_out/`:

| file | contents |
|---|---|
| `out_fixed_console.log` | full console output of the valid run |
| `flagged_samples.csv` | one row per flagged sample, both length definitions |
| `sweep.csv` | 36-cell sweep (9 floors × 2 length sources × 2k1b on/off) |
| `definition_disagreement.csv` | per-floor cross-tab of the two definitions |
| `offgeno_flag_sweep_fixed.R` | the script with the two fixes below applied |
| `console_UNPATCHED_invalid.log` | output of the script as shipped — kept as evidence, **do not use its numbers** |

---

## 1. Two bugs in the script as shipped — both had to be fixed before the output meant anything

The script ran to completion with exit 0 and printed a plausible-looking sweep. Its numbers were wrong.

### Bug 1 — the run glob picks up every result directory under the root, not the five

```r
sum_files <- list.files(root, pattern = "^Summary\\.csv$", recursive = TRUE, full.names = TRUE)
sum_files <- sum_files[basename(dirname(sum_files)) == "summary"]
```

`HCV_paper_revisjon/` also holds `HCVTyper_SRA_data`, `Orig_HCVTyper_SRA_data`, `Orig_HCVTyper_sim_data`,
`v1.2.0-HCVTyper_SRA_data` and `NGS_SEQ-20260210-01`. The last of these was **not re-run with `28a568d`**,
and the four benchmark dirs are different builds entirely. The script loaded **10 runs / 437 samples**
instead of 5 / 140, and reported `Rows (samples): 437 [report: 140]` without stopping.

**Fix:** restrict explicitly to the five runs named in §2 of the handoff.

```r
RUNS_WANTED <- c("NGS_SEQ-20251113-02", "NGS_SEQ-20251212-01", "NGS_SEQ-20260326-01",
                 "NGS_SEQ-20260521-01", "NGS_SEQ-20260625-02")
sum_files <- sum_files[basename(dirname(dirname(sum_files))) %in% RUNS_WANTED]
```

Suggest the script also **hard-fail** when the loaded row count differs from the expected cohort size,
rather than printing the mismatch as a comment.

### Bug 2 — `sampleName` key mismatch silently breaks both the join and the must-keep check

`Summary.csv` carries the suffix (`sampleName = "2743986-HCV"`). The script strips it when building the
`assembly_support.csv` side (`str_remove(sample, "-HCV.*$")` → `"2743986"`) but then joins on the
**unstripped** `sampleName`, and compares the **unstripped** `sampleName` against a bare `MUST_KEEP`.

Consequences in the shipped run:

- `Flagged samples with L_asup: 33 / 112` — the entire comparison arm was empty for the runs of interest.
  The 33 that did resolve were SRA samples (`ERR1810442` etc.), which carry no `-HCV` suffix and so
  happened to match.
- `*** WARNING: only 0 of 3 must-keep samples are in the flagged set ***`, and `must_keep_kept = 0/3`
  in **every row of the sweep**. The disqualifying criterion the whole exercise hangs on never evaluated.

**Fix:** derive a bare id once and key everything on it.

```r
sm <- sm %>% mutate(L_summary = suppressWarnings(as.numeric(denovo_minor_contig_length)),
                    sample_id = str_remove(sampleName, "-HCV.*$"))
# ... join: by = c("run", "sample_id" = "sampleName", "denovo_minor_subtype" = "subtype")
# ... must-keep: filter(sample_id %in% MUST_KEEP)
```

Suggest the script **hard-fail** when fewer than `length(MUST_KEEP)` must-keep samples land in the
flagged set — as written it warns and carries on producing a sweep whose key column is meaningless.

Everything below is from the patched run.

---

## 2. Answers to the four questions in §5 of the handoff

### Q1 — Does the reconstruction hold? **Yes, exactly.**

```
Predicate (classify_roles.R:962-966) fires on : 72
review_flag text actually contains sentence   : 72   [report: 72]
Predicate and text agree on all 140 samples. Reconstruction is exact.
Is it the sole sentence on every flagged sample? YES (as reported)
```

I also checked the transcription against `git show 28a568d:bin/classify_roles.R` directly — the predicate
in the handoff is verbatim, including the hand-rolled `substr(x, 1, 1)`. This was the one genuine gap in
the original report, which identified the 72 by text-matching `review_flag` alone. It closes clean.

### Q2 — How far apart are the two length definitions? **Not far enough to matter.**

| | |
|---|---|
| identical | **68 / 72** |
| differ | 4 |
| median abs. difference | 0 bp |
| max abs. difference | 176 bp |
| subtype label mismatch | 0 / 72 |

The four that differ:

| run | sample | subtype | `denovo_minor_ref` | `L_summary` | `L_asup` | Δ |
|---|---|---|---|---|---|---|
| 20251212-01 | `2634413` | 1b | `1b_EU781827` | 683 | 689 | +6 |
| 20251212-01 | `2637466` | 1b | `1b_EU781827` | 662 | 699 | +37 |
| 20260521-01 | `2746037` | 1b | `1b_D90208` | 628 | 804 | +176 |
| 20260521-01 | `Sample71K` | 3a | `3a_X76918` | 529 | 532 | +3 |

Per-floor cross-tab — **the two definitions disagree on exactly one sample, at one floor**:

| floor | agree | disagree |
|---|---|---|
| 0, 500, 1000, 1250, 1500, 2000, 2500, 3000 bp | 72 | **0** |
| 750 bp | 71 | 1 |

**At 1,000 bp the disagreement is zero.** Per the handoff's own criterion — *"if `disagree` is 0 or 1, the
report's sweep transfers directly and the 1,000 bp recommendation stands as-is"* — it stands as-is.

The BLAST-rank wrinkle flagged in §1 of the handoff (`denovo_minor_ref` chosen by BLAST rank at
`blast_parse.R:308-312`, not by contig size, so the naming reference and the measured contig could
diverge) is real in principle but does not bite on this cohort: the subtype label agrees in 72/72 cases
and the lengths agree in 68/72.

### Q3 — Where do the must-keep samples sit under `L_summary`? **Identical. No shift.**

| sample | run | major | contig | `L_summary` | `L_asup` | pident | k-mer cov |
|---|---|---|---|---|---|---|---|
| `2743986` | 20260521-01 | 3a | 1b | **4,467** | 4,467 | 91.93 | 1.42 |
| `2726018` | 20260326-01 | 3a | 1a | **2,787** | 2,787 | 93.27 | 1.97 |
| `2714375` | 20260326-01 | 3a | 1b | **2,706** | 2,706 | 92.64 | 1.62 |

Byte-identical to the report's figures under both definitions. The recommended floor does not have to move.

This also re-confirms the point that matters most for parameter choice: all three sit **below**
`denovo_min_kmer_cov = 2.0`. The k-mer leg would suppress exactly the samples that most deserve the flag,
and must not be applied to this trigger. Length is the only usable leg.

### Q4 — `flags_kept` and `to_high` at floor = 1000, `exclude_2k1b = TRUE`, `L_summary`

```
flags_kept = 17    (12.1% of the 140-sample cohort)
suppressed = 55
to_high    = 55
must_keep  = 3/3
```

Inside the handoff's predicted range (15–21 flags, ~50–57 to `high`).

**`to_high = 55` is not uncomfortable — it is the intended effect.** All 72 flagged samples are currently
`provisional` and the trigger is their only sentence, so the confidence distribution moves:

| tier | now | after |
|---|---|---|
| high | 38 | **93** |
| provisional | 78 | **23** |
| review | 14 | 14 |
| indeterminate | 10 | 10 |

The 23 remaining `provisional` are the 17 that keep the flag plus 6 that are provisional for other reasons.
That looks like a healthy distribution for a routine cohort. **I would not reach for the two-tier
`Major_evidence` token option** — it adds a second reporting channel to solve a problem the floor already
solves, and the suppressed contigs are 142–999 bp fragments whose observation carries no clinical weight.

---

## 3. One correction to the handoff, and one to my own report

### The 2k1b count is 14, not 13 — your rule is right and mine was the undercount

```
2k1b/genotype-{1,2} pairs among the flagged set: 14   [report: 13]
```

The report used a narrower hand-rolled rule (2k1b contig against a `1b` or `2k` major). `is_valid_minor()`
rule 2 (`classify_roles.R:200-203`) blocks 2k1b against **any** genotype 1 or 2, which additionally catches
`2634016` (20251212-01, 2k1b contig against a **1a** major). Breakdown of the 14: 13× against a 1b major,
1× against 1a. Use the codebase rule.

Note that `2610361` (20251113-02, 2k1b contig against a **3a** major, 1,116 bp) correctly survives the
exclusion — genotype 3 is not in `{1, 2, 2k1b}`.

### The report conflated two length columns without saying so

For the record: §5.3 of `hcvtyper_v1.3.0-g28a568d_vs_v1.1.x_comparison.html` took its headline
distribution (median 606 bp, range 142–4,467) from `denovo_minor_contig_length`, but its floor-sweep tables
and the pident/k-mer analysis from `assembly_support.csv`. The handoff's §1 describes only the second.
The two arms were not distinguished in the write-up — a fair thing to have caught, even though it turned
out to change nothing.

### Parameter defaults now confirmed against the right commit

The report caveated that defaults were read at `5ff7166`, which predates the build. `28a568d` is now in the
repo; `nextflow.config` is unchanged at that commit — `denovo_min_contig_length = 500`,
`denovo_min_kmer_cov = 2.0`, `denovo_min_blast_identity = 90`, `rescue_min_length = 3000`,
`rescue_min_pident = 85`, `rescue_min_aln_length = 3000`, `rescue_min_kmer_cov = 2.0`. Caveat resolved.

---

## 4. Recommendation

**Gate the trigger at 1,000 bp on `denovo_minor_contig_length`, and consult
`genotype_from_subtype()` / `is_valid_minor()` rule 2 for the 2k1b case.**

I withdraw the report's 2,500 bp in favour of your 1,000 bp. The margin argument is decisive: 2,500 leaves
only 206 bp below the smallest must-keep sample (`2714375`, 2,706 bp), which is fragile against any future
change in assembly behaviour. 1,000 bp leaves 1,706 bp of headroom and still cuts the flag from 51% of the
cohort to 12%.

The 2k1b fix is a consistency repair independent of the threshold, as your §6 says: the trigger uses
`substr(dv_minor, 1, 1) != substr(maj_sub, 1, 1)` where `genotype_from_subtype()`
(`genotype_utils.R:24`) already encodes the 2k1b-aware rule. Switching the trigger to it alone removes
14 of the 72 flags before any length floor applies.

### Resulting review worklist — 17 samples

| # | run | sample | major | contig | length | pident | k-mer | must-keep |
|---|---|---|---|---|---|---|---|---|
| 1 | 20260521-01 | `2743986` | 3a | 1b | 4,467 | 91.93 | 1.42 | **yes** |
| 2 | 20260326-01 | `2726020` | 1a | 3a | 3,500 | 92.54 | 1.64 | |
| 3 | 20260326-01 | `2726018` | 3a | 1a | 2,787 | 93.27 | 1.97 | **yes** |
| 4 | 20260521-01 | `2753092` | 3a | 1a | 2,778 | 94.02 | 1.03 | |
| 5 | 20260326-01 | `2726029` | 3a | 1a | 2,735 | 92.91 | 1.11 | |
| 6 | 20260326-01 | `2714375` | 3a | 1b | 2,706 | 92.64 | 1.62 | **yes** |
| 7 | 20251212-01 | `2625911` | 1b | 3a | 2,337 | 92.17 | 2.21 | |
| 8 | 20260326-01 | `2726025` | 1b | 3a | 1,992 | 90.39 | 1.51 | |
| 9 | 20260521-01 | `2745085` | 1a | 3a | 1,939 | 93.46 | 1.41 | |
| 10 | 20260326-01 | `2726024` | 1a | 3a | 1,773 | 92.78 | 1.16 | |
| 11 | 20251212-01 | `2633901` | 1a | 6i | 1,620 | 91.30 | 1.02 | |
| 12 | 20260521-01 | `2740250` | 3a | 1b | 1,445 | 93.01 | 1.57 | |
| 13 | 20251212-01 | `2620329` | 1b | 3a | 1,314 | 93.91 | 1.51 | |
| 14 | 20260625-02 | `Sample62K2` | 3a | 1b | 1,246 | 91.97 | 0.99 | |
| 15 | 20260625-02 | `Sample72K2` | 1b | 3a | 1,176 | 92.69 | 1.25 | |
| 16 | 20251113-02 | `2610361` | 3a | 2k1b | 1,116 | 92.86 | 1.03 | |
| 17 | 20251212-01 | `2627883` | 1a | 3a | 1,036 | 92.55 | 1.16 | |

One to note: `2633901` (#11) has a 1,620 bp 6i contig but a BLAST **alignment length of only 69 bp**
(`best_contig_aln_length` in `blastparse/2633901-HCV.assembly_support.csv`; the script does not carry that
column through to `flagged_samples.csv` — worth adding). A length floor on the contig alone lets it through.
If you want a second leg, `rescue_min_aln_length` — not `denovo_min_kmer_cov` — is the one that would
catch this class without touching the must-keep samples (their alignment lengths are 2,705 / 2,762 / 4,448
bp, i.e. essentially full-contig).

---

## 5. Reproducing

```bash
conda activate /home/jon.brate/.conda/R_shared
Rscript offgeno_sweep_out/offgeno_flag_sweep_fixed.R \
  /mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon \
  ./offgeno_sweep_out 2>&1 | tee offgeno_sweep_out/console.log
```

Note the `assembly_support.csv` glob is still unrestricted in the patched script — harmless, because the
join carries `run`, but it scans all 262 files across all 10 result directories. Restrict it too if the
runtime ever matters.
