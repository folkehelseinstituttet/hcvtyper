# Handoff 260804-dnrank — the de novo dominance comparison measures the wrong thing

**Date:** 2026-08-04
**Base:** `dev` @ `2185e90`. Applies unchanged to `ecda63f` — the four commits in between are
documentation and `nextflow_schema.json` only, touch nothing under `bin/`, and leave every parameter
default cited here (`denovo_min_contig_length = 500`, `denovo_min_kmer_cov = 2.0`,
`denovo_min_blast_identity = 90`, `review_min_offgenotype_contig_length = 1000`) unchanged.
**Cohort:** the 93 samples in `HCV/2026/HCV_paper_revisjon/HCVTyper_SRA_data`, of two kinds:

- **86 `ERR…` samples — public SRA run accessions from Thomson *et al.* (2016)**, our real-world
  benchmark. These are patient plasma samples (single-genotype and co-infected), plus the in vitro RNA
  transcript (IVT) mixtures of H77 (1a, AF009606) and JFH-1 (2a, AB047639) at known ratios from 5000:1 to
  1:5000. Every `ERR…` accession named in this document is public and can be re-downloaded from SRA/ENA.
- **7 simulated datasets** — `sim1`, `sim2`, `sim3`, `sim11asingle`, `sim11bsingle`, `sim22asingle`,
  `sim23asingle` — generated in-house with known genotype composition and mixing ratios. Note the
  filesystem names carry no underscores: `sim11asingle` is the `sim1` 1a single-infection control,
  `sim22asingle` the `sim2` 2a single-infection control, and so on.

Both sets are available to the development team, and every claim below is reproducible from them
(see §7).

**Analysis was read-only.** Nothing was written to `/mnt/N`.

> **Independently verified 2026-08-05** by a fresh end-to-end Nextflow run on a second machine, from raw
> FASTQ, on all 7 simulated datasets plus 3 Thomson accessions. Every claim testable there reproduced —
> most to the exact digit. §9.1's root cause was additionally pinned down and is **not** where this
> document guessed; §9.1 has been corrected accordingly. Full record in §10.
>
> **The §3 patch has since been applied, measured and shipped as `2f254a3`** — two columns move on one
> sample, all 25 call-bearing columns byte-identical, and the D2 control still fires (§10.5).
>
> **§9.1 and §9.4 are now fixed too** (`f9a5591`, `b47f627`), and the §8 IVT dilution series has been
> downloaded and run — every one of its eight predicted contig length / k-mer values reproduced to the
> digit (§10.6). That run also exposed a **sensitivity cliff at 500:1 and beyond** that this document did
> not previously record: see §10.7.

---

## 1. The defect in one sentence

`summarize.R` raises

> *"Co-infection confirmed, but major/minor assignment uncertain — de novo and mapping disagree on which
> strain is dominant. Please review."* → `call_confidence = review`

whenever `denovo_major_subtype_match == "NO"`. But `denovo_major_ref` is **not** an abundance statement,
so the two sides of that comparison do not measure the same thing, and the message asserts a disagreement
the pipeline never tested.

`blast_parse.R:301` sets `major_name <- scaf_top$sseqid[1]` — "best overall hit" — from a frame sorted by
descending **bitscore** (`blast_parse.R:152`). Bitscore is alignment length × identity, so the de novo
"major" is whichever strain sits closest to the **reference panel**, regardless of how much of it is in
the sample.

### The evidence

Top BLAST hit per contig, ordered by bitscore — row 1 becomes `denovo_major_ref`:

| sample | rank 1 | rank 2 | mapping major | agree? |
|---|---|---|---|---|
| **`sim2`** | **3a** bits 17,444, pident **100%**, k-mer 4,184 | 2a bits 15,579, pident 96.04%, k-mer **9,918** | 2a — 507,708 vs 231,132 reads | **NO → false flag** |
| `sim1` | 1a bits 16,761, pident 100%, k-mer 10,626 | 1b bits 13,452, pident 93.13%, k-mer 4,169 | 1a | yes — coincidence |
| **`ERR1810505`** | **1a** bits 14,083, pident 94.45%, k-mer 214 | 3a bits 13,568, pident 93.38%, k-mer **334** | 3a — 7,472 vs 5,645 reads | **NO → false flag** |

`sim2` is a **designed 7:3 mixture**. Every abundance signal says 2a — 2.2× the mapped reads and 2.4× the
k-mer coverage. 3a wins on bitscore only because its simulation source matches `3a_D17763` exactly (100%)
while the 2a source differs from `2a_D00944` (96.04%). `sim1` escapes purely because there the abundant
strain also happens to hold the exact-matching reference; nothing about `sim1` is more correct.

### The pipeline already answers this question correctly, elsewhere

The D2 trigger (`classify_roles.R:711-733`) compares **top-by-targeted-reads** against
**top-by-best-contig-k-mer-coverage** — both abundance measures — and emits *"Dominance ordering
uncertain"*. It correctly stays silent on `sim2` and `ERR1810505`.

So there are two dominance checks in the codebase: one abundance-based and right, one bitscore-based and
wrong. The wrong one is the one that escalates to `review`.

---

## 2. An approach that does NOT work — please don't repeat it

The obvious fix is to re-rank `blast_parse.R`'s major/minor selection by k-mer coverage instead of
bitscore. **I implemented and ran this across all 92 samples. It is worse than the bug.**

The minor slot serves two incompatible purposes:

- **dominance ordering** — wants the *most abundant* strain (k-mer)
- **the off-genotype-contig review trigger** — wants the *most substantial* evidence (length)

Ranking the minor by k-mer makes short, high-coverage fragments win, and they then fall below
`review_min_offgenotype_contig_length` (1,000 bp), silently deleting the flag:

| sample | minor contig before | minor contig after | consequence |
|---|---|---|---|
| `ERR1810503` | 1,475 bp @ 11.0× | **637 bp @ 16.0×** | flag lost — a genuine genotype-3 detection at 11× |
| `ERR1810521` | 1,395 bp @ 16.7× | **845 bp @ 18.4×** | flag lost — the 1:500 IVT true minority |
| `ERR1810491` | 2,132 bp @ 1.17× | 653 bp @ 1.34× | flag lost |
| `ERR1810467` | 1,849 bp @ 1.53× | 540 bp @ 2.56× | flag lost |
| `ERR1810487` | 1,149 bp @ 1.36× | 970 bp @ 1.71× | flag lost |

Five of the seven off-genotype flags in the cohort disappeared, including two that are genuine
minor-strain detections. **The `denovo_*_ref` fields must keep their current bitscore ordering** — they
feed the review triggers, and those triggers want substantiality, not abundance.

If you ever do want an abundance-ranked pair, **add** one (`denovo_dominant_ref` / `denovo_subdominant_ref`)
rather than repurposing the existing fields.

---

## 3. Recommended patch

Fix the *comparison*, not the *ranking*. Two hunks, both reporting-layer.

### 3.1 `bin/classify_roles.R` — remove the co-infection dominance message

```diff
@@ -1109,8 +1109,23 @@
   if (is_indet_dom)
     msgs <- c(msgs, "Dominance ordering uncertain — read-count and k-mer-coverage rankings disagree. ...")
   # Co-infection subtype conflict (sample-level de novo vs mapping).
-  if (is_coinf && subtype_dis)
-    msgs <- c(msgs, "Co-infection confirmed, but major/minor assignment uncertain — de novo and mapping disagree on which strain is dominant. Please review.")
+  #
+  # 260804-dnrank: REMOVED. This asserted a disagreement about DOMINANCE, but
+  # denovo_major_ref / denovo_minor_ref are ordered by BLAST bitscore
+  # (blast_parse.R:301, frame sorted at :152), which ranks strains by proximity to
+  # the reference PANEL, not by abundance. On the designed 7:3 sim2 mixture it
+  # handed "major" to 3a on a 100% panel match despite 2.4x less k-mer coverage and
+  # 2.2x fewer mapped reads, and raised call_confidence = review on a clean call.
+  #
+  # Dominance is already owned by the D2 trigger below, which compares
+  # top-by-targeted-reads against top-by-best-contig-k-mer-coverage -- both
+  # abundance measures -- and emits "Dominance ordering uncertain". D2 correctly
+  # stays silent on sim2 and ERR1810505. Keeping two dominance checks that measure
+  # different things guaranteed they would disagree.
+  #
+  # denovo_*_subtype_match is retained for what it CAN support: the monoinfection
+  # major-subtype cross-check immediately below, where a single strain means no
+  # dominance question arises.
```

### 3.2 `bin/summarize.R` — scope both `subtype_match` legs to non-co-infection calls

```diff
@@ -1989,14 +1989,23 @@
       (!is.na(gate_flag) & gate_flag != "ok") |
-        (!is.na(denovo_major_subtype_match) & denovo_major_subtype_match == "NO") |
+        # 260804-dnrank: scoped to non-co-infection calls. On a co-infection this
+        # was a bitscore-vs-abundance artefact, not a real conflict.
+        (!is.na(denovo_major_subtype_match) & denovo_major_subtype_match == "NO" &
+           !str_detect(coalesce(overall_sample_call, ""), "^co-infection")) |
         (!is.na(rescue_effect) & rescue_effect == "major_ref_changed") |
         overall_sample_call == "co-infection (indeterminate dominance)"  ~ "review",
       coalesce(dominant_unconfirmed, FALSE) |
         ...
-        (!is.na(denovo_minor_subtype_match) & denovo_minor_subtype_match == "NO") |
+        # 260804-dnrank: same artefact, other end of the same bitscore ordering.
+        (!is.na(denovo_minor_subtype_match) & denovo_minor_subtype_match == "NO" &
+           !str_detect(coalesce(overall_sample_call, ""), "^co-infection")) |
```

The monoinfection *"Major subtype conflict"* trigger is deliberately untouched — with one strain there is
no dominance question, so `denovo_major_subtype_match` there means what it says.

**Variant A** (major leg only) was also tested. It leaves `sim2` and `ERR1810505` at `provisional` with an
**empty `review_flag`** — a confidence caveat with nothing for the analyst to act on. Variant B above
moves both to `high`, which is consistent with the documented invariant that `high` always has an empty
flag. Recommend B.

---

## 4. Exact effect on the simulated data

| sample | mapping call | before | after |
|---|---|---|---|
| `sim1` (1a:1b 7:3) | co-infection 1a/1b | `high`, no flag | **unchanged** |
| **`sim2`** (2a:3a 7:3) | co-infection 2a/3a | **`review`** + dominance message | **`high`, no flag** |
| `sim3` (4a single) | monoinfection 4a | `high`, no flag | **unchanged** |
| `sim11asingle` | monoinfection 1a | `high`, no flag | **unchanged** |
| `sim11bsingle` | monoinfection 1b | `high`, no flag | **unchanged** |
| `sim22asingle` | monoinfection 2a | `provisional` + background-candidate note | **unchanged** |
| `sim23asingle` | monoinfection 3a | `provisional` + background-candidate note | **unchanged** |

`sim2` is the only simulated sample that moves, and it moves to the correct answer: the genotypes, read
counts (507,708:231,132), coverage (99.99:99.98) and dominance ordering were all right before and after —
only the spurious escalation is gone.

`sim22asingle` / `sim23asingle` are **deliberately unaffected**. Their `provisional` comes from a separate
issue — a spurious `background`/`weak` second candidate lowering the tier — described in §9.2.

---

## 5. Exact effect on the full 93-sample cohort

**Only two columns move at all:** `review_flag` (7 samples) and `call_confidence` (2 samples).

### Confidence distribution

| tier | before | after |
|---|---|---|
| high | 52 | **54** |
| provisional | 17 | 17 |
| review | 22 | **20** |
| indeterminate | 2 | 2 |

### The 7 samples whose `review_flag` text changes

All seven lose the same sentence and keep every other sentence they carried. `ERR…` are Thomson *et al.*
(2016) SRA run accessions; `sim2` is the simulated 2a:3a 7:3 mixture.

| sample | call | confidence | what remains |
|---|---|---|---|
| `ERR1810454` | co-infection | provisional → provisional | candidate-2 marginal-corroboration note |
| `ERR1810465` | co-infection | review → review | major failed mapping QC + candidate-2 note |
| `ERR1810496` | co-infection | review → review | major failed mapping QC + candidate-2 note |
| `ERR1810497` | co-infection | review → review | major failed mapping QC + candidate-2 note |
| `ERR1810499` | co-infection | review → review | major failed mapping QC |
| `ERR1810505` | co-infection | **review → high** | (nothing — correctly clean) |
| `sim2` | co-infection | **review → high** | (nothing — correctly clean) |

Five of the seven keep `review` on their *own* merits. Only the two false positives clear entirely.

### Nothing is gained, nothing call-bearing moves

`review_flag` count 39 → 37, **0 gained**. All 25 call-bearing columns byte-identical across all 93
samples: `Major`, `Minor`, `Major_subtype`, `Minor_subtype`, `Major_genotype_mapping`,
`Minor_genotype_mapping`, `major_typable`, `minor_typable`, `overall_sample_call`, `NS34A_short`,
`NS5A_short`, `NS5B_short`, `GLUE_genotype`, `GLUE_subtype`, `Reference`, `Major_reference`,
`Minor_reference`, `Major_evidence_state`, `Minor_evidence_state`, `Major_role_reason`,
`Minor_role_reason`, `Major_cov_breadth_min_10`, `Minor_cov_breadth_min_10`,
`Reads_nodup_mapped_major`, `Reads_nodup_mapped_minor`.

All 7 off-genotype-contig flags survive: `ERR1810451`, `ERR1810467`, `ERR1810487`, `ERR1810491`,
`ERR1810493`, `ERR1810503`, `ERR1810521`.

---

## 6. Why this cannot change a genotype call

`blastparse.csv` — the file carrying `major_ref` / `minor_ref` — is consumed by **nothing but SUMMARIZE**.
`workflows/hcvtyper.nf:515` routes `BLASTPARSE.out.csv` into `ch_denovo` (SUMMARIZE staging only), while
reference selection and rescue take `BLASTPARSE.out.support` (`assembly_support.csv`) at line 386. The
recommended patch does not touch `blast_parse.R` at all, so even that indirection is moot: both hunks sit
downstream of every call-bearing decision.

---

## 7. Reproducing this

Nextflow is not required. Every input SUMMARIZE and BLASTPARSE consume is already published, so both
scripts run directly on a results directory (~150 MB staging for these 93 samples, a few minutes).

**Establish fidelity first.** Replaying with the *unpatched* scripts reproduced the published output
byte-for-byte — 184/184 `blastparse.csv` + `assembly_support.csv` files, and the whole `Summary.csv`.
Without that check, harness bugs are indistinguishable from real regressions.

Stage **by channel, not by glob**:

| staged dir | source | trap |
|---|---|---|
| `kraken_classified/` | `*.focused.kraken2.report.txt` | `*.entireDB.kraken2.report.txt` also matches |
| `stats_markdup/` | `<sample>.nodup.idxstats` | `*.firstmapping.nodup.idxstats` also matches |
| `depth/` | `<sample>.<ref>.cand*.nodup.tsv` | `*.firstmapping.withdup.sorted.tsv` also matches |
| `parsefirst_mapping/` | `parsefirstmapping/*.parsefirstmapping.csv` **+ `blastparse/*.rescued.candidates.csv`** | the second is post-rescue; `parsefirstmapping/*.candidates.csv` is pre-rescue and wrong |
| `denovo/` | `blastparse/*.{blastparse,assembly_support,_blast_out}.csv` | — |

Invocations (reference panel is on the `test-datasets` branch at
`blast_db/HCVgenosubtypes_8.5.19_clean.fa`, 224 seqs; contigs are published gzipped):

```bash
blast_parse.R <sample> blast/<sample>.txt spades/<sample>.contigs.fa \
              HCVgenosubtypes_8.5.19_clean.fa HCV
summarize.R samplesheet.csv 1.3.0 folkehelseinstituttet/hcvtyper \
            500 2.0 90 genotype true 499 29 2 3.0 1.0 0.5 1.0 1000
```

`summarize.R` sources helpers by bare relative path, so run it from a directory holding the six R scripts
and the staged input dirs. Verify with:

```bash
Rscript bin/tests/compare_summary_regression.R OLD/Summary.csv NEW/Summary.csv review_flag,call_confidence
```

Expect **exit 0** — under this patch the narrow allowlist is sufficient, which is the strongest available
statement that nothing else moved.

---

## 8. Suggested tests

1. **`sim2` as a named fixture.** A designed 7:3 mixture where one component matches the panel exactly and
   the other does not is precisely the configuration that broke this. Assert `call_confidence = high` and
   an empty `review_flag`.
2. **A dominance-disagreement fixture that must still fire**, so the removal cannot be over-applied:
   `ERR1810469` and `ERR1810510` are both `co-infection (indeterminate dominance)` via D2 and must keep
   `review`. Both are unchanged by this patch — worth locking in.
3. **Off-genotype flag count.** Assert all 7 survive; this is what the rejected approach in §2 broke.
4. **The IVT dilution series as a sensitivity fixture.** This would have caught the §2 approach
   immediately, since `ERR1810521` is the 1:500 minority. These are Thomson *et al.* SRA accessions, so
   the fixture inputs are public and re-downloadable rather than dependent on `/mnt/N` — the nine mixtures
   are `ERR1810511`/`513`/`515`/`517`/`519`/`521`/`523`/`525`/`527`.

   The series is worth building out properly: nine 1a/2a mixtures spanning 1:1 to 1:5000 in both
   directions, where the minority strain is known by construction. De novo recovers it at **every**
   dilution — 2a at 2,578 bp/326× (5:1), 3,075 bp/232× (50:1), 349 bp/1.4× (500:1), 920 bp/1.7× (5000:1);
   1a at 4,585 bp/373× (1:5), 1,563 bp/111× (1:50), 1,395 bp/16.7× (1:500), 505 bp/0.94× (1:5000). Nothing
   in `bin/tests/` currently asserts *sensitivity* — that the pipeline still detects what it should — only
   that output has not changed. This is the one dataset that can support such a test.

   **Status 2026-08-05: downloaded and run; every number above reproduced to the digit (§10.6).** The
   accession → ratio mapping is now known and tabulated there, so the fixture can be written. Note the
   run also found a sensitivity cliff at 500:1 (§10.7) — a fixture asserting only "minority recovered"
   would pass while the Summary reports `monoinfection` at `high` with no caveat, so assert the
   *reported* outcome per ratio, not just the contig's existence.

---

## 9. Related defects found in the same replay, not fixed here

These surfaced while diagnosing the above and are recorded here so they are not re-derived from scratch.
§9.1 is the most serious — a genuine regression in shipped code.

### 9.1 REGRESSION — coverage and depth silently blank when a sample has two candidates but one depth file

> **FIXED in `f9a5591` (2026-08-05).** Root cause below; the fix drops `Minor_reference` from the
> `summarize.R:1306` join key (and drops `df_coverage`'s own copy of that column, which would otherwise
> collide into `.x`/`.y` and break the 136-column schema).
>
> **This section understated the impact.** It is written as coverage and depth columns going `NA`, but the
> coverage those columns hold is what the typability gate reads — so the bug was also forcing
> **`major_typable` to `NO`**, a call-bearing column, on the affected samples. On `sim22asingle` /
> `sim23asingle` — simulated single-genotype data with ~692k reads at 100% breadth and an exact reference
> match — the pipeline was reporting the major strain as untypable. Recovered values match the
> recomputation table below exactly (cov@10× 99.99; mean depth 20,762 and 21,045).

On **15 of the 93 samples**, `Major_cov_breadth_min_1/5/10` and `Major_avg_depth` come out `NA` even
though the `cand1` depth file exists and its reference matches `Major_reference` exactly.

The trigger condition is exact:

| candidates | depth files | outcome |
|---|---|---|
| 1 | 1 | works (`sim11asingle`, `sim3`, …) |
| 2 | 2 | works (`sim1`, `sim2`, `ERR1810447`, …) |
| 2 | 1 | **all Major coverage/depth columns NA** |

Affected: `ERR1810454`, `ERR1810473`, `ERR1810477`, `ERR1810480`, `ERR1810481`, `ERR1810500`,
`ERR1810501`, `ERR1810506`, `ERR1810512`, `ERR1810514`, `ERR1810516`, `ERR1810518`, `ERR1810524`,
`sim22asingle`, `sim23asingle`.

**The `candidate_rank` join is NOT the cause — that dead end is already ruled out.** It is the obvious
suspect and it is innocent: replaying the join (`depth/` filenames → `candidate_rank_lookup`,
`summarize.R:577-580`) recovers `candidate_rank = 1` correctly for every affected sample, with **0 NA
ranks across all 91 depth rows**.

**ROOT CAUSE (identified 2026-08-05, §10).** It is the *later* join, not `df_coverage` itself.
`df_coverage` is innocent too — its `fill(.direction = "downup") %>% slice(1)` collapse was the suspect
named in the original draft of this section, and it is wrong: on an affected sample `df_coverage`
produces a perfectly correct row carrying `Major_cov_breadth_min_10` and a `Minor_reference` of `NA`.

The fault is `summarize.R:1306`:

```r
left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference))
```

`df_coverage` derives **both** reference columns from the staged `depth/` filenames, so a sample whose
rank-2 candidate never got a depth file has `Minor_reference = NA` on that side, while the main frame
carries the real candidate name from `candidates.csv`. dplyr matches `NA` to `NA` by default
(`na_matches = "na"`), which is exactly why the other two rows of the trigger table survive:

| candidates | depth files | `df_coverage$Minor_reference` | main frame | join |
|---|---|---|---|---|
| 1 | 1 | `NA` | `NA` | matches (NA-to-NA) → **works** |
| 2 | 2 | `2b_D10988` | `2b_D10988` | matches → **works** |
| 2 | 1 | `NA` | `1m_KJ439778` | **no match → whole coverage block NA** |

One mechanism explains all three rows. The fix belongs in the join key — `df_coverage` is already one row
per sample, so `Minor_reference` earns nothing as a key and only introduces the failure. Dropping it has
**not** been tested; do that before shipping, and note that joining on `sampleName` alone would collide
the two reference columns, so `join_by(sampleName, Major_reference)` is the narrower change.

**The pileups are intact; only the summary fails to carry them through.** Recomputing breadth and mean
depth directly from `samtools/<sample>.<ref>.cand1.nodup.tsv` reproduces the expected values exactly:

| sample | expected cov@10× | recomputed | expected depth | recomputed |
|---|---|---|---|---|
| `ERR1810454` | 89.00 | 89.00 | 16 | 16.1 |
| `ERR1810473` | 0.78 | 0.78 | 3 | 2.8 |
| `ERR1810477` | 98.42 | 98.42 | 55 | 55.3 |
| `ERR1810480` | 57.32 | 57.32 | 11 | 10.6 |
| `ERR1810481` | 85.58 | 85.58 | 19 | 18.7 |
| `ERR1810501` | 98.9 | 98.91 | — | 188.2 |
| `sim22asingle` | 100.0 | 99.99 | — | 20,762 |
| `sim23asingle` | 100.0 | 99.99 | — | 21,045 |

Worth a regression test asserting that a sample with a `cand1` depth file always gets a non-NA
`Major_cov_breadth_min_10`.

### 9.2 Spurious second candidates, and a correctly-rejected candidate lowering confidence

All three simulated single-infection datasets carry a cross-mapping second candidate that first mapping
correctly gated out (`minor_call = no`):

| sample | spurious candidate | first-map reads | first-map cov | retained as rank 2? | confidence |
|---|---|---|---|---|---|
| `sim11asingle` | `1i_KJ439772` | 11,664 | **17%** | **no** | `high` |
| `sim22asingle` | `1m_KJ439778` | 8,260 | 6% | **yes** | `provisional` |
| `sim23asingle` | `1n_KJ439775` | 9,688 | 6% | **yes** | `provisional` |

The candidate with the *highest* coverage is the one dropped, and the two at 6% are kept — whatever
governs retention is not monotone in the evidence. This retention is also what creates the
two-candidates-one-depth-file condition in §9.1, so the NA coverage on those two samples is downstream of
the same inconsistency.

Once retained, the candidates are handled **correctly** — `background` / `no_own_assembly` / `weak`, no
co-infection called. The issue is that a correctly-demoted background candidate drags `call_confidence`
to `provisional`, which asserts the *call* is less certain when the pipeline has in fact just done its
job. On a simulated single-genotype dataset with 692k reads at 100% coverage and an exact reference
match, `provisional` is the wrong tier. Suggest keeping the traceability note but not letting a demoted
`weak`/`background` candidate lower the tier.

### 9.3 `ERR1810505`'s `denovo_*` columns populate at HEAD where they were `NA` before

Under `28a568d` every `denovo_*` column on `ERR1810505` was `NA`; at HEAD they populate. That is a gain —
but it is what exposed the §1 defect on that sample, and it means the previous build was silently
suppressing a review flag rather than not raising one. Worth confirming the change was deliberate and
adding a fixture.

### 9.4 `major_contig_length` can describe a different contig than `major_ref`

> **FIXED in `b47f627` (2026-08-05).** The major slot now uses the same ONE ROW discipline as the minor
> slot: `scaf_top` row 1 is the row that defines `major_ref`, so reference, contig and length are all read
> off it. `sim1` 9,339 → **9,076**, `sim2` 9,695 → **9,446** (where 9,695 was literally the *minor* 2a
> contig's length sitting in the major field). Verified OK on all 9 IVT samples too (§10.6).
>
> Two notes for the reader. First, this section's risk is lower than it looks:
> `denovo_major_contig_length` has **no logic consumers** — only `denovo_minor_contig_length` feeds a
> threshold (the off-genotype trigger at `summarize.R:1867`) — so the defect could only ever mislead a
> human reading the row, never change a call. Second, `major_contig` was derived a *third* way (longest
> contig in `scaf_top`, unfiltered) and was assigned but never read; it is now correct and still internal.

`blast_parse.R`'s `major_contig_length` is derived as
`scaf %>% filter(sseqid == major_name) %>% slice_max(sc_length)`, which searches the **full** BLAST table
rather than the top-hit-per-contig frame. On `sim1` that returns 9,339 bp — the length of the **1b**
contig, which carries a secondary 78.7%-identity hit against `1a_HQ850279` — while `denovo_major_ref` is
`1a_HQ850279`, whose own contig is 9,076 bp. Two fields, two contigs: the same grain mismatch `78ec956`
fixed for the minor slot, still live in the major slot. Not touched here because it is orthogonal to the
dominance defect, but worth its own patch.

---

## 10. Independent verification, 2026-08-05

### 10.1 What was run

A **fresh end-to-end Nextflow run from raw FASTQ** on a different machine from the one that produced §1-§9
— not a replay of published outputs. This matters: §7's reproduction recipe re-runs `summarize.R` and
`blast_parse.R` over an existing results directory, so it shares that directory's staging. A full pipeline
run shares nothing, and it reproduced the same numbers.

- **Code:** `dev` @ `d64d35a`, **plus the uncommitted working-tree changes present on that machine** —
  principally `conf/modules.config` (`BOWTIE2_BUILD` resources + `publishDir` disabled; `FASTQC_TRIM`
  prefix `.trimmed` → `_trimmed`), plus `assets/multiqc_config.yml`, `nf-test.config`, `tests/`, and CI.
  **None of them touch `bin/`**, so no call-bearing logic differs from `d64d35a` — but the run is *not*
  byte-reproducible from `origin/dev` alone until that work is committed.
- **Samples (10):** all 7 simulated datasets, plus `ERR1810469`, `ERR1810451`, `ERR1810447`.
- **Profile:** `-profile docker`, GLUE enabled, `kraken_all = false`.

Only 10 of the 86 `ERR…` accessions were available locally, which set the ceiling on what §5 and §2 could
be checked against (§10.4).

### 10.2 Confirmed — reproduced to the exact digit

**§1, the core defect.** `sim2` top hits per contig, ordered by bitscore, rebuilt from
`sim2_blast_out.csv`:

| rank | sseqid | bitscore | pident | k-mer cov |
|---|---|---|---|---|
| 1 | `3a_D17763` | 17,444 | 100 | 4,184.4 |
| 2 | `2a_D00944` | 15,579 | 96.041 | 9,918.2 |

`blastparse.csv` gives `denovo_major_ref = 3a_D17763`; mapping makes 2a the major at **507,708 vs
231,132** reads. `sim1` reproduces likewise (1a: 16,761 / 100% / 10,626 — 1b: 13,452 / 93.128% / 4,168.6)
and escapes for the stated reason, not by being more correct.

`sim2`'s `candidates.csv` shows why D2 is right to stay silent: rank 1 is 2a on **both** abundance
measures (624,520 targeted reads, k-mer 9,918) and rank 2 is 3a (273,886, k-mer 4,184). The abundance
signals agree with each other and disagree only with bitscore.

**§4, all seven rows.** Reproduced exactly, including `sim2` at `review` carrying the message verbatim,
and `sim22asingle` / `sim23asingle` at `provisional` on the background-candidate note.

**§9.1.** Trigger table reproduced with no exceptions across all 10 samples; only `sim22asingle` and
`sim23asingle` are affected. Root cause then identified — see the correction in §9.1 above.

**§9.2.** `sim11asingle` carries a single candidate (the 17% `1i_KJ439772` is dropped outright, absent
from `candidates.csv`); `sim22asingle` retains `1m_KJ439778` at 8,260 reads / 6%; `sim23asingle` retains
`1n_KJ439775` at 9,688 reads / 6%. Both retained ones are `below_threshold` / `weak` / `background`. The
non-monotonicity is real.

**§9.4.** `sim1` has `major_ref = 1a_HQ850279` with `major_contig_length = 9339` — the 1b contig's length,
not the 1a contig's 9,076 bp.

**§6, statically.** `BLASTPARSE.out.csv` has exactly one consumer, `ch_denovo` at
`workflows/hcvtyper.nf:515`; rescue and reference selection take `.support` at `:386`. Confirmed by grep
over `workflows/` and `subworkflows/` — no other reference exists.

**The D2 trigger** at `classify_roles.R:709-733` does compare `targeted_reads_nodup` against
`assembly_support_best_contig_kmer_cov`, both abundance measures, as §1 claims.

### 10.3 The two Thomson controls behave as the patch requires

| sample | result | why it matters |
|---|---|---|
| `ERR1810469` | `co-infection (indeterminate dominance)`, `review`, flagged **by D2** | §8 test 2's negative control. Confirms D2 fires independently of the message §3.1 removes, so the removal cannot be over-applied. |
| `ERR1810451` | off-genotype-contig flag present — a 1,538 bp 2b contig at k-mer 2.3 against a 4a monoinfection | §8 test 3 / the flag class the §2 approach destroyed. |

`ERR1810451` also incidentally strengthens §2: its genuine flag rests on a 1,538 bp contig at k-mer 2.3 —
exactly the long-but-low-coverage profile that k-mer re-ranking displaces in favour of a short fragment.

### 10.4 Not verified — and why

Everything below needs Thomson accessions that were not on the verification machine. **None of it was
contradicted; it simply could not be reached.**

- **§5's cohort-wide tables** — the 93-sample confidence distribution and the 7 changed `review_flag`
  texts. Only 10 `ERR…` were available.
- **§2's five-flags-lost table** (`ERR1810503`, `ERR1810521`, `ERR1810491`, `ERR1810467`, `ERR1810487`)
  and **§1's second false positive `ERR1810505`**.
- **§8 test 4, the IVT dilution series** (`ERR1810511`-`527`).
- **§9.3**, which needs a `28a568d` build of `ERR1810505`.
- **Every §9.1-affected `ERR…`** — though the mechanism now established in §9.1 predicts them all.

### 10.5 The patch's "after" state — measured, and it holds

The §3 hunks were applied and `summarize.R` re-run standalone against the run's SUMMARIZE work directory.
**Fidelity was established first**, as §7 insists: the *unpatched* replay reproduced the pipeline's
`Summary.csv` and `candidates.csv` byte-for-byte, so the delta below is the patch and not the harness.

(That check earned its keep. The first replay attempt mis-staged — `trimmed/` and `variation/` were
missing — and silently `NA`-ed four read-count columns. Without the fidelity gate that would have read as
a patch effect.)

**Exactly two columns move, on exactly one sample:**

| sample | column | before | after |
|---|---|---|---|
| `sim2` | `call_confidence` | `review` | **`high`** |
| `sim2` | `review_flag` | the dominance sentence | **empty** |

- **All 25 call-bearing columns byte-identical** across all 10 samples — 0 differences, confirming §5's
  "nothing call-bearing moves" on the sample subset available here.
- **`ERR1810469` holds `review`**, flag unchanged, still raised by D2 — §8 test 2's control, and the
  evidence that the removal is not over-applied: the abundance-based dominance check still fires once the
  bitscore-based one is gone.
- `ERR1810451` keeps its off-genotype flag; `ERR1810447` holds `provisional`.
- `sim22asingle` / `sim23asingle` stay `provisional` on the background-candidate note, per §4.

Confidence distribution `high 4 / provisional 3 / review 3` → `high 5 / provisional 3 / review 2` — the
same shape as §5's cohort prediction: `review` down, `high` up, `provisional` untouched, nothing gained.

**Variant B's rationale holds.** `sim2` lands at `high` with an empty flag, consistent with the documented
invariant. Variant A would have left it `provisional` with nothing for the analyst to act on.

Shipped as `2f254a3`, which also drops the locals the removal left dead in `classify_roles.R` (`is_coinf`,
`subtype_dis`, `min_match` — the `denovo_minor_subtype_match` *parameter* is kept, since it is the 3rd
positional arg and dropping it would silently shift every positional caller). That cleanup was re-run and
is byte-neutral.

**Still open at that point:** §5's full 93-sample tables, §2's five lost flags, `ERR1810505` and §9.3.
The IVT series has since been closed — see §10.6.

### 10.6 The IVT dilution series — downloaded, run, and every §8 number reproduced

§10.4 listed the IVT series as unreachable. It is not: all nine accessions are public, and
`nf-core/fetchngs@1.12.0` retrieves them in minutes (464 MB total). They were then run end-to-end through
the pipeline at `b47f627`, i.e. with the §3 patch and both 260805 fixes in place.

Two fetchngs gotchas, since the next person will hit them: `--input` rejects a `.txt` extension outright,
and it reads a CSV **header row as an accession** (`Mixture of ids provided via --input: id`). Use a
headerless `.csv`.

**Every predicted contig length and k-mer coverage reproduced to the digit**, which also pins down the
accession → ratio mapping that §8 did not state:

| accession | ratio (1a:2a) | minority | §8 predicted | measured | call | conf |
|---|---|---|---|---|---|---|
| `ERR1810511` | 1:1 | — | — | both ≈1,550× | co-infection | high |
| `ERR1810513` | 1:5 | 1a | 4,585 bp / 373× | **4,585 / 373.03** | co-infection | high |
| `ERR1810515` | 5:1 | 2a | 2,578 bp / 326× | **2,578 / 326.40** | co-infection | high |
| `ERR1810517` | 1:50 | 1a | 1,563 bp / 111× | **1,563 / 110.96** | co-infection | high |
| `ERR1810519` | 50:1 | 2a | 3,075 bp / 232× | **3,075 / 232.43** | co-infection | high |
| `ERR1810521` | 1:500 | 1a | 1,395 bp / 16.7× | **1,395 / 16.74** | monoinfection | provisional |
| `ERR1810523` | 500:1 | 2a | 349 bp / 1.4× | **349 / 1.38** | monoinfection | high |
| `ERR1810525` | 1:5000 | 1a | 505 bp / 0.94× | **505 / 0.94** | monoinfection | high |
| `ERR1810527` | 5000:1 | 2a | 920 bp / 1.7× | **920 / 1.71** | monoinfection | high |

§8's claim holds: **de novo recovers the minority strain at every dilution, 1:1 through 1:5000.**

`ERR1810521`'s off-genotype flag fires with exactly the predicted evidence — *"1395 bp contig, 1395 bp
aligned (100%), 98.7% identity, k-mer cov 16.7"*. This is the flag §2 says k-mer re-ranking destroys, now
observed firing on its own data.

**§9.4 verified on all nine**: `major_contig_length` equals the length of `major_ref`'s own contig in
every case, checked by rebuilding `scaf_top` from the raw BLAST output rather than trusting the field.

**§9.1 was NOT exercised here.** All nine samples have matching candidate and depth-file counts (5×2/2,
4×1/1), so none reaches the 2-candidates/1-depth-file shape. This run shows the fix does no harm on
unseen data; it is not independent positive evidence. `sim22asingle` / `sim23asingle` remain the only
demonstration, and a dedicated fixture is still worth building.

One incidental observation: at 1:5000 the 1a minority's best hit is to `1a_M67463`, a different 1a
reference than the `1a_AF009606` used at every other ratio. Genotype recovery is unaffected, but a
fixture asserting an exact reference name there would be brittle.

### 10.7 NEW FINDING — a sensitivity cliff at 500:1, where a known mixture reports as clean

Not in any previous section, and the reason the IVT series is worth keeping as a permanent fixture.

Across the series the pipeline gives **three qualitatively different answers**, and the boundary is set
entirely by `review_min_offgenotype_contig_length` (1,000 bp):

| ratios | outcome | minority visible to the analyst? |
|---|---|---|
| 1:1 – 50:1 | co-infection, `high` | **yes** — reported as the Minor strain |
| 1:500 | monoinfection, `provisional` + off-genotype flag | **yes** — as a caveat (contig 1,395 bp > 1,000) |
| 500:1, 1:5000, 5000:1 | monoinfection, **`high`, empty `review_flag`** | **no** — contig 349/505/920 bp < 1,000 |

So on a sample that is a two-strain mixture *by construction*, and whose assembly demonstrably contains
the second strain, the pipeline can emit `monoinfection` at `high` confidence with no caveat whatsoever.
`ERR1810527` (5000:1) misses the flag by 80 bp — its 920 bp contig would have been flagged at 1,000.

This is the threshold behaving as documented rather than a defect, but it is a sharper sensitivity
statement than "de novo recovers it at every dilution" suggests, and it belongs in any paper claim about
co-infection detection limits. It also shows §2's mechanism in the wild: what separates a flagged
minority from a silent one is contig length crossing 1,000 bp — precisely the quantity the rejected
k-mer re-ranking would have shortened.

Worth deciding deliberately: whether `high` is the right tier for a monoinfection call that has a
sub-threshold off-genotype contig sitting in its own assembly. A `provisional` tier, or a note without a
tier change, would at least leave a trace.
