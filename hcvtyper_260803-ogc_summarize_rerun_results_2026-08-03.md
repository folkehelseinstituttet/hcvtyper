# Results: re-run and evaluation of `summarize.R` after quick task 260803-ogc

**Responds to:** `hcvtyper_handoff_260803-ogc_summarize_rerun_validation_2026-08-03.md` (commit `2b80f69`)
**Date:** 2026-08-03
**Cohort:** the five runs / 140 samples under `HCV/2026/HCV_paper_revisjon/`
**Old side:** `v1.3.0-g28a568d` — the published `Summary.csv` files, and a replay from them
**New side:** `origin/dev` @ `58af4df`
**Read-only.** Nothing was written to `/mnt/N`.

## Verdict

**Clean pass.** Five exit-0 regression runs, 17 samples still flagged, `high` = 93, and every strain and
resistance column byte-identical across all 140 samples. The central claim of the task — that the change
alters reporting and nothing else — holds empirically.

Nothing in the §6 "Serious" rows fired.

---

## 1. The re-run had to be done differently than §3 describes

**`-resume` was not available.** The machine this was run on has no hcvtyper Nextflow work directory, no
write access to `/mnt/N`, and the five runs live in separate published folders. The `cache = false`
recipe in §3 needs a cache to invalidate.

**It turned out Nextflow is not needed at all.** Every input `SUMMARIZE` and `BLASTPARSE` consume is
already published to the results directory, and `summarize.R` simply globs flat directories. So both
processes can be invoked directly on the published outputs. Total staging for all five runs: **146 MB**.

This sidesteps §3 entirely — with the scripts called directly there is no task hash, so the
"`blast_parse.R` is on `PATH` and therefore invisible to the cache" problem cannot arise.

A side benefit: replaying `blast_parse.R` from the published **raw** BLAST table means BLAST itself is
never re-executed, so the only thing that can move is `blast_parse.R`. That also makes the open
megablast-vs-blastn question in §8 orthogonal to this validation.

### The fidelity check that makes the rest trustworthy

Before testing anything, both legs of the harness were replayed with the **old** scripts (`28a568d`) and
compared against the published output:

| leg | scope | result |
|---|---|---|
| `summarize.R` | all five runs | **5/5 `Summary.csv` byte-for-byte identical** |
| `blast_parse.R` | run 20251212-01, all samples | **44/44 `blastparse.csv` + `assembly_support.csv` byte-for-byte identical** |

This also settles the parameter question. The runs were launched without explicit params, so the
`nextflow.config` defaults applied, and byte-exact reproduction confirms the reconstructed `ext.args`:

```
500 2.0 90 genotype true 499 29 2 3.0 1.0 0.5 1.0        # old, args[4..15]
500 2.0 90 genotype true 499 29 2 3.0 1.0 0.5 1.0 1000   # new, args[4..16]
```

`minRead = 499` and `minCov = 29` come from `conf/modules_hcv.config:16-17`, not `nextflow.config`.

### Two traps for anyone repeating this

Both produced a full set of *false* "Serious" regressions — `Minor_subtype` differing, Sample51K
2c→3a, `denovo_minor_ref` going `NA`. The fidelity check is what exposed them as harness bugs.

**(a) Stage by channel, not by glob.** The published directories merge outputs from several processes, so
a bare glob double-counts and `summarize.R` emits one row per file (45 rows for a 23-sample run):

| staged dir | correct source | trap |
|---|---|---|
| `kraken_classified/` | `KRAKEN2_FOCUSED.out.report` | `kraken2/` also holds `*.entireDB.kraken2.report.txt` |
| `stats_markdup/` | `JOINT_MAPPING.out.idxstats_nodup` | `samtools/` also holds `*.firstmapping.nodup.idxstats` |
| `depth/` | `JOINT_MAPPING.out.depth` | `samtools/` also holds `*.firstmapping.withdup.sorted.tsv` |
| `parsefirst_mapping/` | `PARSEFIRSTMAPPING.out.csv` **+ `RESCUE_EVALUATION.out.candidates`** | the second is `blastparse/*.rescued.candidates.csv` (post-rescue), **not** `parsefirstmapping/*.candidates.csv` (pre-rescue) |

**(b) `BLASTPARSE`'s input is `blast/<sample>.txt`** — the raw headerless BLASTN table.
`blastparse/*_blast_out.csv` is BLASTPARSE's *output*; feeding it back in yields header-only
`blastparse.csv` files and a cascade of spurious diffs.

Invocation that works:

```bash
blast_parse.R <sample> blast/<sample>.txt spades/<sample>.contigs.fa \
              HCVgenosubtypes_8.5.19_clean.fa HCV      # no ext.args on this module
summarize.R   samplesheet.csv 1.3.0 folkehelseinstituttet/hcvtyper <ext.args above>
```

Contigs are published gzipped; decompress them for `seqinr::read.fasta`. The reference panel is on the
`test-datasets` branch at `blast_db/HCVgenosubtypes_8.5.19_clean.fa` (224 sequences), so no download is
needed. `summarize.R` sources its helpers by bare relative path, so it must run with the working directory
containing the six R scripts and the staged input directories.

---

## 2. The six items §7 asks for

### 2.1 `compare_summary_regression.R` exit status

```
NGS_SEQ-20251113-02  exit=0  PASS       NGS_SEQ-20260521-01  exit=0  PASS
NGS_SEQ-20251212-01  exit=0  PASS       NGS_SEQ-20260625-02  exit=0  PASS
NGS_SEQ-20260326-01  exit=0  PASS
```

Every difference confined to the six allowed columns. The **narrow** form
(`review_flag,call_confidence`) fails on all five as designed, naming only `denovo_minor_contig`,
`denovo_minor_contig_length`, `Major_genotype`, `Minor_genotype` — no unexpected column appears.

### 2.2 Confidence distribution

| tier | before | after | §5 predicted |
|---|---|---|---|
| high | 38 | **93** | ~93 |
| provisional | 78 | **23** | ~23 |
| review | 14 | **14** | 14 |
| indeterminate | 10 | **10** | 10 |

### 2.3 Flagged-sample count

**72 → 17.** 55 suppressed, **0 gained**. The 17 are the same samples the sweep predicted, at the same
lengths, with all three must-keep samples retained:

| # | run | sample | len | | # | run | sample | len |
|---|---|---|---|---|---|---|---|---|
| 1 | 20260521-01 | `2743986` **(must-keep)** | 4467 | | 10 | 20260326-01 | `2726024` | 1773 |
| 2 | 20260326-01 | `2726020` | 3500 | | 11 | 20251212-01 | `2633901` | 1620 |
| 3 | 20260326-01 | `2726018` **(must-keep)** | 2787 | | 12 | 20260521-01 | `2740250` | 1445 |
| 4 | 20260521-01 | `2753092` | 2778 | | 13 | 20251212-01 | `2620329` | 1314 |
| 5 | 20260326-01 | `2726029` | 2735 | | 14 | 20260625-02 | `Sample62K2` | 1246 |
| 6 | 20260326-01 | `2714375` **(must-keep)** | 2706 | | 15 | 20260625-02 | `Sample72K2` | 1176 |
| 7 | 20251212-01 | `2625911` | 2337 | | 16 | 20251113-02 | `2610361` | 1116 |
| 8 | 20260326-01 | `2726025` | 1992 | | 17 | 20251212-01 | `2627883` | 1036 |
| 9 | 20260521-01 | `2745085` | 1939 | | | | | |

### 2.4 Movement on the two contig columns

`denovo_minor_contig`: **43** samples. `denovo_minor_contig_length`: **8** samples.

Per §6, 8 is more than "a handful", so the follow-up was done: **three cross the 1,000 bp floor, and none
changes the flagged set.**

| run | sample | old len | new len | crosses 1000 | why the flag is unaffected |
|---|---|---|---|---|---|
| 20251113-02 | `2617257` | 1164 | 444 | yes | 2k1b vs a 1b major — already excluded by the pair rule |
| 20260625-02 | `2761293` | 1380 | 672 | yes | 2k1b vs a 1b major — already excluded by the pair rule |
| 20260326-01 | `2726019` | 704 | 9414 | yes | co-infection, so the monoinfection trigger never applies |
| 20251212-01 | `2634392` | 937 | 807 | no | below the floor either way |
| 20260326-01 | `2726026` | 543 | 516 | no | below the floor either way |
| 20260521-01 | `Sample71K` | 529 | 453 | no | below the floor either way |
| 20260521-01 | `Sample74T` | 474 | 436 | no | below the floor either way |
| 20260625-02 | `Sample74T2` | 474 | 241 | no | below the floor either way |

**No sweep re-derivation is needed.**

### 2.5 `review_flag` text, verbatim

`2633901` (20251212-01) — `denovo_minor_contig` now names `NODE_3_length_1620_cov_1.021433`, not
`NODE_2_length_3232_cov_6285.477617`. Major stays `1a`, NS3/4A stays `122G`, confidence stays `provisional`.

> Monoinfection called for candidate 1 (1a_HQ850279), but de novo assembly found a different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0; only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly supported — the contig may be largely non-HCV, chimeric, or too divergent to type. Please review.

`2724348` (20260521-01), Major subtype conflict — stays at `review`:

> Major subtype conflict for candidate 1 (3a_D28917) — mapping (3a) vs contig (1a) — 2004 bp contig, 2004 bp aligned (100%), 93.5% identity, k-mer cov 1.2; possible reference mismatch or highly divergent strain. Please review.

The annotation discriminates correctly: 2633901's 4% alignment gets the artefact caveat, 2724348's 100%
alignment does not. `2753257` behaves the same way (3540 bp, 100% aligned, 91.4%).

### 2.6 "Serious" rows — none

Byte-identical across all 140 samples: `Major`, `Minor`, `Major_subtype`, `Minor_subtype`,
`Major_genotype_mapping`, `Minor_genotype_mapping`, `major_typable`, `minor_typable`,
`overall_sample_call`, `NS34A_short`, `NS5A_short`, `NS5B_short`, `GLUE_genotype`, `GLUE_subtype`,
`Reference`, `denovo_major_ref`, **`denovo_minor_ref`**.

`denovo_minor_ref` being untouched is worth stating explicitly: the *selection* is unchanged and only the
reported contig name moved, which is exactly the intent of `78ec956`.

---

## 3. The remaining §5 expectations

| Metric | Expected | Observed |
|---|---|---|
| Column count, 20251113-02 / 20260625-02 | 134 → 136 | **134 → 136** |
| Column count, other three runs | 136 | **136** |
| `Major_genotype` on Sample51K | 2 → 3 | **2 → 3** (and `Minor_genotype` 3 → 2 — both un-swapped) |
| Row counts | unchanged | **23 / 22 / 28 / 36 / 31**, unchanged |

Two consistency checks beyond §7, both clean:

- samples whose flag was suppressed but whose confidence did **not** move: **0 of 55**
- samples that **kept** the flag but whose confidence moved: **0 of 17**

`bash bin/tests/run_all.sh` on `origin/dev`: **ALL R TESTS PASSED**.

---

## 4. Suggestions for the tooling

1. **`compare_summary_regression.R` is doing its job** — it caught the harness bugs above as loudly as it
   would have caught a real regression. No change needed.
2. **Consider recording the resolved params** in `pipeline_info/` (e.g. the Nextflow `params.json`).
   Reconstructing `ext.args` from `conf/` + `nextflow.config` worked here only because the runs used
   defaults; a run with overrides could not be replayed this way at all.
3. **A `--only-runs` guard on `bin/tests/offgeno_flag_sweep.R`.** Its `list.files(recursive = TRUE)` picks
   up every result directory under the root, and `HCV_paper_revisjon/` also contains the SRA/sim benchmark
   dirs plus `NGS_SEQ-20260210-01`, which was not re-run with this build. The abort-on-wrong-cohort-size
   check added in this task already catches it — the explicit allowlist would just make it self-documenting.
4. **Carry `best_contig_aln_length` into `flagged_samples.csv`.** The new annotation proves its worth
   (4% vs 100% is what separates 2633901 from 2724348), but the sweep's per-sample CSV still omits it.
