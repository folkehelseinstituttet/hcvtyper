# Handoff: re-running and evaluating `summarize.R` after quick task 260803-ogc

**Date:** 2026-08-03
**Branch:** `dev`, commits `0686f6d` … `58af4df` (13 commits)
**Cohort:** the same five runs / 140 samples used for
`hcvtyper_v1.3.0-g28a568d_vs_v1.1.x_comparison.html` and
`hcvtyper_offgenotype_flag_sweep_results_2026-08-03.md`
**Purpose:** confirm empirically that the change altered **only** what it was supposed to alter.

The central claim of this task is that gating the off-genotype-contig review flag changes *reporting*
and nothing else — no subtype, no typability, no `overall_sample_call`, no resistance result. That claim
was argued from the code (`review_flag` and `call_confidence` are read by nothing but MultiQC styling and
the CSV itself). This run turns the argument into a check.

---

## 1. What changed, and what should move

| Commit | Change | Observable effect |
|---|---|---|
| `0686f6d` | Off-genotype trigger gated on a 1,000 bp contig floor + the 2k/1b pair rule | `review_flag`, `call_confidence` |
| `ab7b947` | `Major_genotype`/`Minor_genotype` derived from the role-corrected subtypes | those two columns; adds them to 2 runs that lacked them |
| `2c1e5c1` | Measured contig evidence written into the off-genotype sentence | `review_flag` text |
| `78ec956` | `denovo_minor_contig` now names the contig the selection chose | `denovo_minor_contig`, `denovo_minor_contig_length` |
| `f6ca7c9` | Measured contig evidence written into the Major-conflict sentence | `review_flag` text |

Everything else in the 13 commits is tests, tooling or documentation.

**Six columns are permitted to differ.** Any other difference is a finding, not a pass.

```
review_flag  call_confidence  Major_genotype  Minor_genotype
denovo_minor_contig  denovo_minor_contig_length
```

---

## 2. Before you start

```bash
git -C /path/to/hcvtyper log --oneline -1        # expect 58af4df or later on dev
bash bin/tests/run_all.sh                        # expect: ALL R TESTS PASSED (13 files)
```

If the unit suite is red, stop — nothing below is meaningful.

**Publish to a NEW output directory.** Do not re-use the existing `HCV_paper_revisjon/<run>/` folders. The
2026-07-03 dev outputs were overwritten in place by a previous re-run, which is why the earlier comparison
had to take its "previous snapshot" column from a published HTML document rather than from files. Keep the
old `Summary.csv` intact or there is nothing to diff against.

---

## 3. Re-running — the one thing that will bite you

**`-resume` alone will crash the run.** This is not a guess; the failure was reproduced:

```
The following named parsers don't match the column names: minor_contig
ERROR: Can't rename columns that don't exist.
✖ Column `minor_contig` doesn't exist.
```

Why: Nextflow's task hash covers the process script, `ext.args` and **declared** input files. It does not
cover scripts that are merely on `PATH` from `bin/`.

- **SUMMARIZE re-runs correctly.** `classify_roles.R` is a declared `path` input
  (`workflows/hcvtyper.nf:543`), and its `ext.args` gained
  `review_min_offgenotype_contig_length`. Either alone invalidates the cache.
- **BLASTPARSE does NOT re-run.** `blast_parse.R` is invoked from `PATH` and is not a declared input, and
  the module's `ext.args` did not change. With `-resume` it stays cached and re-publishes the **old**
  `blastparse.csv` — which has no `minor_contig` column. SUMMARIZE then aborts on the rename above.

### Recommended: targeted re-run

Write `rerun-260803-ogc.config`:

```groovy
process {
    // blast_parse.R changed but is not a declared input, so Nextflow cannot see it.
    // Force BLASTPARSE to re-execute; everything upstream stays cached.
    withName: '.*:BLASTPARSE' { cache = false }
}
```

Then, per run:

```bash
nextflow run folkehelseinstituttet/hcvtyper \
    -r dev \
    -resume \
    -c rerun-260803-ogc.config \
    --input  <the run's original samplesheet> \
    --outdir <NEW output dir for this run> \
    <the same profile / params as the original run>
```

Trimming, mapping, assembly and GLUE stay cached. BLASTPARSE is a `process_single` R step and is cheap to
redo. SUMMARIZE re-runs because its declared inputs, its `ext.args` and its BLASTPARSE inputs all changed.

**Fallback:** a clean run without `-resume` is always correct, just far more expensive (re-does SPAdes and
GLUE). Use it if `-resume` behaves unexpectedly.

### Confirm the new parameter actually arrived

```bash
grep -o 'summarize.R.*' work/*/*/.command.sh | tail -1
```

The argument list is positional; `review_min_offgenotype_contig_length` is **args[16]**, the last one. If
it is absent the R side falls back to its own default of 1000, so results are unaffected — but its absence
means your config was not picked up, which is worth knowing before you interpret anything else.

---

## 4. Evaluating

### 4.1 The automated check

Per run:

```bash
Rscript bin/tests/compare_summary_regression.R \
    OLD/<run>/summary/Summary.csv \
    NEW/<run>/summary/Summary.csv
```

Exit 0 = every difference is confined to the six allowed columns. Exit 1 = it names the offending columns
and prints up to ten example samples each.

Also run the **narrow** form, which is the stronger statement:

```bash
Rscript bin/tests/compare_summary_regression.R \
    OLD/<run>/summary/Summary.csv NEW/<run>/summary/Summary.csv \
    review_flag,call_confidence
```

This excludes the genotype-column fix and the contig-coherence fix from the allowlist, so it will fail —
deliberately — and show you exactly which samples those two changes touched. That is the intended way to
inspect them, not a problem to fix.

### 4.2 The numbers to check by hand

The differ proves *what* changed; these confirm it changed by the *expected amount*.

```r
library(tidyverse)
read_all <- function(root) {
  list.files(root, "^Summary\\.csv$", recursive = TRUE, full.names = TRUE) %>%
    keep(~ basename(dirname(.x)) == "summary") %>%
    map_dfr(read_csv, col_types = cols(.default = col_character()))
}
old <- read_all("OLD"); new <- read_all("NEW")

# 1. confidence distribution
bind_rows(old = count(old, call_confidence), new = count(new, call_confidence), .id = "v") %>%
  pivot_wider(names_from = v, values_from = n)

# 2. how many samples carry the off-genotype sentence
c(old = sum(str_detect(old$review_flag, "different-genotype contig"), na.rm = TRUE),
  new = sum(str_detect(new$review_flag, "different-genotype contig"), na.rm = TRUE))

# 3. did the gate column move? (the contig-coherence fix can shift it)
sum(old$denovo_minor_contig_length != new$denovo_minor_contig_length, na.rm = TRUE)
```

---

## 5. Expected results

| Metric | Before | After | Notes |
|---|---|---|---|
| Samples carrying the off-genotype sentence | 72 | **~17** | 55 suppressed: 51 below the 1,000 bp floor, 4 more by the 2k1b pair rule |
| `call_confidence` high | 38 | **~93** | the 55 suppressed samples had no other review sentence |
| `call_confidence` provisional | 78 | **~23** | 17 that keep the flag + 6 provisional for other reasons |
| `call_confidence` review | 14 | **14** | unchanged — the Major-conflict annotation changes wording, not tier |
| `call_confidence` indeterminate | 10 | **10** | unchanged |
| Column count, runs 20251113-02 / 20260625-02 | 134 | **136** | `Major_genotype`/`Minor_genotype` now always emitted |
| Column count, other three runs | 136 | **136** | unchanged |
| `Major`, `Minor`, `*_subtype`, `minor_typable`, `overall_sample_call`, `NS34A_short`, `NS5A_short`, `NS5B_short` | — | **identical on all 140** | the whole point |
| `Major_genotype` on Sample51K | 2 | **3** | was swapped against `Major_subtype = 3a` |
| `denovo_minor_contig_length` differing | — | **a handful at most** | see §6 |

The "~" on the first three is honest: `denovo_minor_contig_length` gates the flag and the coherence fix can
move it, so 17 is the expectation, not a guarantee. A result of 15–20 needs no explanation; 5 or 40 does.

**Spot-check two samples by eye:**

- **2633901** (run 20251212-01) — should still be flagged, now reading `1620 bp contig, 69 bp aligned (4%),
  91.3% identity, k-mer cov 1.0; only 4% of the contig aligns…`. `denovo_minor_contig` should now name
  **NODE_3**, not `NODE_2_length_3232`. Major stays 1a, NS3/4A stays 122G.
- **2724348** or **2753257** (Major subtype conflict) — the sentence should now carry the conflicting
  contig's length, aligned length and fraction. If the fraction is low, it will say the conflict may be an
  artefact of a short anchor. Both stay at `call_confidence = review`.

---

## 6. Triage — what an unexpected difference means

| Unexpected diff | Reading | Action |
|---|---|---|
| `Major`, `Minor`, `Major_subtype`, `Minor_subtype` | **Serious.** The change was supposed to be display-only | Stop and report. Do not ship |
| `minor_typable`, `overall_sample_call` | **Serious.** Same | Stop and report |
| `NS34A_short`, `NS5A_short`, `NS5B_short` | **Serious.** The GLUE join moved | Stop and report |
| `denovo_minor_ref` | Unexpected — the *selection* logic was not touched, only which contig name is reported | Report; likely a BLASTPARSE re-run difference, not this change |
| `denovo_minor_contig_length` on more than a handful | The coherence fix moved more than expected. Only samples where several contigs top-hit the same reference should move | Check whether the flagged set changed as a result; re-derive the sweep if so |
| `denovo_minor_contig` on many samples | Expected on samples where an on-genotype contig outscored the off-genotype one on the same reference (the 2633901 shape). Frequency is unknown | Note the count; not a blocker |
| Flagged count far from 17 | Either the gate column moved (above), or the parameter did not arrive | Check `.command.sh` per §3, then `denovo_minor_contig_length` |
| `call_confidence` moved on a sample that **kept** its flag | Unexpected — keeping the flag should keep the tier | Report with the sample's `review_flag` |

If the sweep needs re-deriving, `bin/tests/offgeno_flag_sweep.R` does it and now also emits `aln_len` and
`aln_frac_pct` per flagged sample. It aborts rather than printing plausible numbers if the cohort size is
wrong or the three must-keep samples are missing.

---

## 7. What to report back

1. Exit status of `compare_summary_regression.R` for each of the five runs.
2. The confidence-distribution table (§4.2 item 1), old vs new.
3. The flagged-sample count, old vs new.
4. How many samples moved on `denovo_minor_contig_length` and `denovo_minor_contig`.
5. The `review_flag` text for 2633901 and for one Major-conflict sample, verbatim.
6. Anything in the §6 "Serious" rows, immediately.

A clean result is: five exit-0 runs, ~17 flagged, high ≈ 93, and the strain/resistance columns byte-identical
across all 140 samples.

---

## 8. Still open after this run

Neither is a blocker for shipping the above.

- **`-task blastn` vs the default megablast.** Measured on simulated contigs against the 224-sequence
  panel: identical below ~10% divergence, but megablast **misses 9% of queries at 15% divergence and 71% at
  22%** — it goes silent rather than wrong. The exposure is therefore false negatives on divergent strains,
  not spurious flags. Cost is ~49× runtime, and unlike this task it **can change strain calls**
  (`assembly_support` feeds `evidence_state` and `rescue_evaluation.R`). Do the cheap read-only diagnostic
  first: the distribution of top-hit alignment fraction per contig across the published
  `blastparse/*_blast_out.csv` files.
- **Reference-panel composition.** `HCVgenosubtypes_8.5.19_clean.fa` is 224 sequences: genotype 6 = 81
  (36%), genotype 4 = 47, genotype 2 = 35, genotype 1 = 35 (only **5 × 1a**), genotype 3 = 17 (only
  **4 × 3a**). One-reference-per-named-subtype, so it is dominated by the rare, finely-split genotypes and a
  short ambiguous anchor is disproportionately likely to be labelled genotype 6 — which is exactly what
  happened to 2633901.
