# Handoff: sample 2633901 — `denovo_minor_contig` names a different contig than `denovo_minor_contig_length` measures

**Date:** 2026-08-03
**Sample:** `2633901-HCV`
**Run:** `NGS_SEQ-20251212-01`, results under
`/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon/NGS_SEQ-20251212-01/`
**Build:** `folkehelseinstituttet/hcvtyper v1.3.0-g28a568d`, Nextflow 25.04.3
**Context:** this is sample #11 of the 17-sample review worklist in
`hcvtyper_offgenotype_flag_sweep_results_2026-08-03.md` §4 — the one flagged there as having a
1,620 bp contig but only a 69 bp BLAST alignment. Drilling into it surfaced a reporting defect.
**Read-only.** Nothing was written to `/mnt/N`.

---

## 1. The defect

`Summary.csv` reports three de-novo-minor fields for this sample that do not describe the same object:

| field | value |
|---|---|
| `denovo_minor_ref` | `6i_DQ835770` |
| `denovo_minor_contig_length` | **1,620** |
| `denovo_minor_contig` | **`NODE_2_length_3232_cov_6285.477617`** |

**1,620 bp is `NODE_3`. `NODE_2` is 3,232 bp.** The length reported does not belong to the contig named.

This is the wrinkle anticipated in §1 of `hcvtyper_handoff_offgenotype_flag_sweep_2026-08-03.md` —
*"`denovo_minor_ref` is chosen by BLAST rank, not by contig size (`bin/blast_parse.R:308-312`), so the
reference that names the flag and the contig whose length is reported can, in principle, be different
contigs"* — occurring in practice. The sweep found the two length definitions agreed on 68/72 samples, so
the effect is rare, but it is not hypothetical.

The two selections that disagree:

- `denovo_minor_ref` / `denovo_minor_contig` — best-scoring different-genotype BLAST hit across all contigs.
  The winner is `NODE_2`'s hit to `6i_DQ835770` (bitscore 1,074), so `NODE_2` is named.
- `denovo_minor_contig_length` (and `blastparse/*.assembly_support.csv`) — longest contig *assigned to* the
  6i subtype group. `NODE_2`'s own top hit is 1a, so it sits in the 1a group; the longest contig whose top
  hit is 6i is `NODE_3` at 1,620 bp. That length is reported against `NODE_2`'s name.

**Suggested fix:** derive `denovo_minor_contig_length` from `denovo_minor_contig` directly, so the reported
length, k-mer coverage and identity always describe the named contig. If the per-subtype grouping is wanted
for a different purpose, report it under a distinct column name rather than reusing
`denovo_minor_contig_length`.

**Impact on threshold work:** any length floor gated on `denovo_minor_contig_length` is, for samples like
this one, filtering on a contig that is not the one the flag is about. Here it happens not to change the
outcome (both candidates are spurious — §3), but the mapping between column and object should be repaired
before the floor is wired in.

---

## 2. The call itself is sound

| field | value |
|---|---|
| `Major` / `Major_subtype` | **1a** |
| `Major_reference` | `1a_HQ850279` |
| `overall_sample_call` | monoinfection |
| `call_confidence` | **provisional** |
| `major_typable` / `minor_typable` | YES / NO |
| `Major_role_reason` | dominant |
| `Major_evidence_state` | **confirmed** |
| `Major_dominance_score` | 9.41638 |
| `gate_flag` | ok |
| `cv_evenness` | 0.71405 |
| `rescue_flag` / `rescue_effect` | FALSE / none |
| `GLUE_genotype` / `GLUE_subtype` | 1 / 1a |
| `denovo_major_subtype_match` | YES |

Read budget: 6,966,158 raw → 6,062,762 trimmed → 2,022,682 classified → 2,014,667 mapped
(`fraction_mapped_reads_vs_median` = 2.538).

Resistance: NS3/4A **S122G** — *probable resistance* to grazoprevir; glecaprevir, paritaprevir and
voxilaprevir no resistance. NS5A and NS5B: no resistance across all agents.

`review_flag` (the only sentence emitted):

> Monoinfection called for candidate 1 (1a_HQ850279), but de novo assembly found a different-genotype contig (6i) — possible missed co-infection or contamination. Please review.

`Major_evidence`:

> 1a_HQ850279 | 2,633 reads (nodup) | 100% breadth@10x | de novo 1a | 93.7% consensus id | ref from first-mapping | contig support: 5159 bp contig, 93.9% contig identity, 9828.93 k-mer cov | evidence_state=confirmed

---

## 3. The 6i flag is an artefact, twice over

SPAdes produced 107 contigs. Three carry the story.

### The genuine major contig

| | |
|---|---|
| contig | `NODE_1_length_5159_cov_9828.934817` |
| length | **5,159 bp** |
| k-mer coverage | **9,828.93** |
| best hit | `1a_HQ850279`, **93.915%** identity over **5,160 bp**, bitscore 7,788, e = 0 |
| subject span | 9,159 → 4,001 (reverse) — ~5.2 kb of the 3′ half |
| `cand_1_assembly_support` | **supported** |

Next-best hits are all 1a (92.1–92.8%), then a clean drop to 84.8% for genotype-1 unassigned references.
Unambiguous.

### Candidate A for the flag — `NODE_3`, assembly noise

| | |
|---|---|
| contig | `NODE_3_length_1620_cov_1.021433` |
| length | 1,620 bp |
| k-mer coverage | **1.02** |
| BLAST hits, total | **5** |

All five hits fall inside a ~133 bp window of the contig (q 1483–1615, **8% of its length**); the remaining
~1.5 kb matches nothing in the reference set.

| subject | subtype | pident | aln len | q span | s span | bitscore |
|---|---|---|---|---|---|---|
| `6i_DQ835770` | 6i | 91.304 | **69** | 1539–1607 | 1629–1561 | 97.1 |
| `2_JF735116` | 2 | 91.176 | 68 | 1532–1599 | 1638–1571 | 93.5 |
| `2c_JX227949` | 2c | 88.312 | 77 | 1539–1615 | 1546–1471 | 91.6 |
| `6v_EU798760` | 6v | 80.702 | 114 | 1483–1596 | 1685–1572 | 89.8 |
| `6w_EU643836` | 6w | 88.235 | 51 | 1539–1589 | 1353–1303 | 62.1 |

A 1.02× k-mer contig with a 69 bp top alignment is assembly noise. It is what
`blastparse/2633901-HCV.assembly_support.csv` records for the 6i row:

```
2633901-HCV,6i,6i_DQ835770,1620,91.304,69,1.021433
```

### Candidate B for the flag — `NODE_2`, a 1a contig hitting a conserved region

| | |
|---|---|
| contig | `NODE_2_length_3232_cov_6285.477617` |
| length | 3,232 bp |
| k-mer coverage | 6,285.48 |
| **own top hit** | `1a_HQ850279`, **93.038%** over **3,088 bp**, q 1–3088, s 3087–1, bitscore **4,510** |
| 6i hit | `6i_DQ835770`, 88.737% over 879 bp, q 2347–3224, s **890–14**, bitscore 1,074 |
| 6i hit (second HSP) | 71.503% over 1,351 bp, q 682–1995, s 2560–1242, bitscore 302 |

`NODE_2` is a 1a contig by a factor of four in bitscore. Its 879 bp hit to the genotype-6 reference lands on
subject positions **14–890** of a 9,447 bp genome — the 5′UTR and the start of core. 88.7% identity there,
against ~70% genome-wide 1a↔genotype-6 divergence, is consistent with a conserved-region alignment rather
than genotype-6 material. The 1a alignment corroborates the placement: `NODE_2` maps to positions 1–3087 of
`1a_HQ850279` in reverse, so the 6i-hitting segment (q 2347–3224) sits at the genome 5′ end on both.

### Consequence for the length-floor work

This sample is why a contig-length floor alone is not sufficient. The 1,000 bp floor recommended in
`hcvtyper_offgenotype_flag_sweep_results_2026-08-03.md` §4 admits `NODE_3` on its 1,620 bp, even though only
69 bp of it aligns. **`rescue_min_aln_length` is the second leg that would catch this class** without
touching the three must-keep samples, whose alignments are essentially full-contig (2,705 / 2,762 / 4,448 bp).
`denovo_min_kmer_cov` would also catch it (1.02 ≪ 2.0) but must not be used, for the reason established in
the sweep results: the must-keep samples all sit below 2.0 themselves.

---

## 4. First mapping

`parsefirstmapping/2633901-HCV.parsefirstmapping.csv`:

```
total_mapped_reads  2,014,667
major_ref  1a_HQ850279   major_reads  1,937,589   major_cov  100
minor_ref  1_AJ851228    minor_reads     16,170   minor_cov   35
minor_call  yes          gate_flag  ok
```

Share of first-mapping reads: major **96.17%**, minor **0.80%**
(`percent_mapped_reads_major_firstmapping` / `..._minor_firstmapping`).

Both candidates cleared the floor — `parsefirstmapping/2633901-HCV.candidates.csv`:

| rank | reference | subtype | genotype | reads | cov | status |
|---|---|---|---|---|---|---|
| 1 | `1a_HQ850279` | 1a | 1 | 1,937,589 | 100 | pass |
| 2 | `1_AJ851228` | 1 | 1 | 16,170 | 35 | pass |

Read counts above are **with duplicates**. Deduplicated, the first-mapping BAM tells a much starker story —
reads spread across **199 references** with ≥1 nodup read (12,197 nodup total), the usual capture-panel
cross-recruitment:

| reference | nodup reads |
|---|---|
| `1a_HQ850279` | 2,561 |
| `1a_AF009606` | 1,089 |
| `1a_M62321` | 940 |
| `1a_EF407457` | 753 |
| `1a_M67463` | 577 |
| … | … |
| `1_AJ851228` (candidate 2) | **163** |

### Candidate 2 was correctly dropped

`blastparse/2633901-HCV.rescue_audit.csv`:

```
event=blocked_collapse
original: 1_AJ851228 (subtype 1)  ->  target: 1a_HQ850279 (subtype 1a, genotype 1)
evidence: denovo 1a_HQ850279 contig 5159bp pident=93.915 aln=5160bp kmer_cov=9828.93 (replaced 1_AJ851228)
disposition=blocked
```

The de novo evidence wanted to replace the genotype-1 unassigned reference with `1a_HQ850279`, but that
reference was already held by candidate 1, so the same-genotype collapse guard blocked the rescue and
`rescued.candidates.csv` retains a single candidate. Correct behaviour — `1_AJ851228` was 1a
cross-recruitment, not a second strain.

---

## 5. Final mapping and consensus

Reference `1a_HQ850279`, 9,191 bp.

| metric | value |
|---|---|
| `Reads_withdup_mapped_major` | 1,984,732 |
| `Reads_nodup_mapped_major` | **2,633** |
| duplication rate | **≈99.87%** |
| `Percent_reads_mapped_of_trimmed_with_dups_major` | 32.736% |
| `Major_cov_breadth_min_1` | 100% |
| `Major_cov_breadth_min_5` | 99.99% |
| `Major_cov_breadth_min_10` | **99.97%** |
| `Major_avg_depth` | **64.94×** |
| positions ≥100× (from depth tsv) | 11.38% |
| `Major_consensus_similarity_pct` | **93.6888%** |
| `Major_consensus_n_differences` | 580 (over 9,190 bp) |

The ~99.9% duplication is expected for this capture protocol and is why the with-dup and nodup read counts
differ by three orders of magnitude; it is not specific to this sample.

---

## 6. Summary

- **1a monoinfection, well supported**: 100% first-mapping coverage, 99.97% breadth at 10×, 64.9× mean depth,
  a 5.16 kb 1a contig at 93.9% identity with k-mer coverage ~9,829, and GLUE independently calling 1a.
- **The only thing holding it at `provisional` is the 6i review flag**, which is spurious under either
  interpretation of which contig triggered it.
- **One reporting defect to fix**: `denovo_minor_contig_length` does not describe `denovo_minor_contig`.
- **One reinforcement for the threshold work**: a contig-length floor alone admits this sample;
  `rescue_min_aln_length` is the second leg that would exclude it safely.

### Files used

```
NGS_SEQ-20251212-01/summary/Summary.csv
NGS_SEQ-20251212-01/parsefirstmapping/2633901-HCV.parsefirstmapping.csv
NGS_SEQ-20251212-01/parsefirstmapping/2633901-HCV.candidates.csv
NGS_SEQ-20251212-01/blastparse/2633901-HCV_blast_out.csv          (1,137 hits)
NGS_SEQ-20251212-01/blastparse/2633901-HCV.blastparse.csv
NGS_SEQ-20251212-01/blastparse/2633901-HCV.assembly_support.csv
NGS_SEQ-20251212-01/blastparse/2633901-HCV.rescue_audit.csv
NGS_SEQ-20251212-01/blastparse/2633901-HCV.rescued.candidates.csv
NGS_SEQ-20251212-01/samtools/2633901-HCV.firstmapping.nodup.idxstats
NGS_SEQ-20251212-01/samtools/2633901-HCV.firstmapping.withdup.stats
NGS_SEQ-20251212-01/samtools/2633901-HCV.1a_HQ850279.cand1.nodup.tsv
NGS_SEQ-20251212-01/consensus/2633901-HCV.cand1.consensus_distance.tsv
NGS_SEQ-20251212-01/spades/2633901-HCV.contigs.fa.gz               (107 contigs)
```
