# Handoff 260804-dnrank — the de novo dominance comparison measures the wrong thing

**Date:** 2026-08-04
**Base:** `dev` @ `2185e90`. Applies unchanged to `ecda63f` — the four commits in between are
documentation and `nextflow_schema.json` only, touch nothing under `bin/`, and leave every parameter
default cited here (`denovo_min_contig_length = 500`, `denovo_min_kmer_cov = 2.0`,
`denovo_min_blast_identity = 90`, `review_min_offgenotype_contig_length = 1000`) unchanged.
**Cohort:** the 93 samples in `HCV/2026/HCV_paper_revisjon/HCVTyper_SRA_data`, of two kinds:

- **86 `ERR…` samples — public SRA run accessions from Thomson *et al.* (2016)**, the real-world
  benchmark used in the manuscript. These are patient plasma samples (single-genotype and co-infected),
  plus the in vitro RNA transcript (IVT) mixtures of H77 (1a, AF009606) and JFH-1 (2a, AB047639) at known
  ratios from 5000:1 to 1:5000. Every `ERR…` accession named in this document is public and can be
  re-downloaded from SRA/ENA.
- **7 simulated datasets** — `sim1`, `sim2`, `sim3`, `sim11asingle`, `sim11bsingle`, `sim22asingle`,
  `sim23asingle` — generated in-house with known genotype composition and mixing ratios. Note the
  filesystem names are the manuscript names with underscores stripped: `sim11asingle` = `sim1_1a_single`,
  `sim22asingle` = `sim2_2a_single`, and so on.

**Analysis was read-only.** Nothing was written to `/mnt/N`.

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
| `ERR1810503` | 1,475 bp @ 11.0× | **637 bp @ 16.0×** | flag lost — cited in the manuscript |
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
issue (a spurious `background`/`weak` second candidate lowering the tier), written up as item 8b of
`hcvtyper_backlog_from_paper_replay_2026-08-04.md`.

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
4. **The IVT dilution series as a sensitivity fixture** (backlog item 1) would have caught the §2 approach
   immediately, since `ERR1810521` is the 1:500 minority. These are Thomson *et al.* SRA accessions, so
   the fixture inputs are public and re-downloadable rather than dependent on `/mnt/N` — the nine mixtures
   are `ERR1810511`/`513`/`515`/`517`/`519`/`521`/`523`/`525`/`527`.

---

## 9. Related, not fixed here

From `hcvtyper_backlog_from_paper_replay_2026-08-04.md`:

- **item 6** — coverage/depth silently `NA` on 15 of 93 samples whenever a sample has two candidates but
  one depth file. A genuine regression; the values are fully recoverable from the published depth files.
- **item 8b** — spurious `background`/`weak` candidates dragging `call_confidence` to `provisional` even
  though the pipeline correctly rejected them. This is what makes `sim22asingle` / `sim23asingle`
  `provisional`.
- **item 5** — `ERR1810505`'s `denovo_*` columns were `NA` under `28a568d` and populate at HEAD, which is
  what exposed this defect there in the first place.

One further observation, not raised as a defect. `blast_parse.R`'s `major_contig_length` is derived as
`scaf %>% filter(sseqid == major_name) %>% slice_max(sc_length)`, which searches the **full** BLAST table
rather than the top-hit-per-contig frame. On `sim1` that returns 9,339 bp — the length of the **1b**
contig, which carries a secondary 78.7%-identity hit against `1a_HQ850279` — while `denovo_major_ref` is
`1a_HQ850279`, whose own contig is 9,076 bp. Two fields, two contigs: the same grain mismatch `78ec956`
fixed for the minor slot, still live in the major slot. Not touched here because it is orthogonal to the
dominance defect, but worth its own patch.
