# Handoff: missing avg-depth / coverage-breadth values in Summary.csv (2026-06-27 run)

**Date:** 2026-06-27
**Author:** Jon Bråte (FHI)
**Audience:** HCVTyper development team
**Scope:** Documents a data-completeness gap in `summary/Summary.csv` of the genotype-diverse run
(`folkehelseinstituttet/hcvtyper 1.3.0-gfe2db6e`, SRA dir `…/2026/HCVTyper_SRA_data/`, Summary.csv
2026-06-27 06:43): for 15 samples the targeted-mapping **average depth and coverage-breadth columns are
empty**, even though the samples were genotyped, a consensus was built, and the per-position depth files
exist. The question this handoff answers: **are these real missing values?** Short answer: **the values
are missing from Summary.csv, but the underlying data is not lost** — it is a summarisation/join gap, not a
coverage, biology, or data-loss problem, and every value is recoverable from files already on disk.

---

## 1. Symptom

For 15 samples, all of the following columns are blank (`NA`/empty) **as a block**:

- `Major_avg_depth`
- `Major_cov_breadth_min_1`
- `Major_cov_breadth_min_5`
- `Major_cov_breadth_min_10`

It is **always the whole block together** — there is no sample with `Major_avg_depth` present but breadth
blank, or vice versa. All other columns for these samples are populated normally (genotype call,
`Reads_nodup_mapped_major`, `Major_consensus_similarity_pct`, de novo fields, etc.).

### Affected samples (n = 15)

| Sample | Major | Major_reference | reads (nodup) | call |
|---|---|---|---|---|
| ERR1810454 | 4d | 4d_EU392172 | 1,530 | monoinfection |
| ERR1810473 | 3a | 3a_D17763 | 198 | monoinfection |
| ERR1810477 | 3a | 3a_D17763 | 3,666 | monoinfection |
| ERR1810480 | 3a | 3a_D17763 | 1,032 | monoinfection |
| ERR1810481 | 3a | 3a_D17763 | 1,293 | monoinfection |
| ERR1810500 | 3a | 3a_D17763 | 813 | monoinfection |
| ERR1810501 | 3a | 3a_D17763 | 12,660 | co-infection |
| ERR1810506 | 1a | 1a_EF407457 | 7,855 | monoinfection |
| ERR1810512 | 2a | 2a_AB047639 | 1,047 | co-infection |
| ERR1810514 | 1a | 1a_AF009606 | 652 | co-infection |
| ERR1810516 | 2a | 2a_AB047639 | 1,372 | monoinfection |
| ERR1810518 | 1a | 1a_AF009606 | 3,382 | co-infection |
| ERR1810524 | 2a | 2a_AB047639 | 956 | monoinfection |
| sim22asingle | 2a | 2a_D00944 | 692,053 | monoinfection |
| sim23asingle | 3a | 3a_D17763 | 695,614 | monoinfection |

---

## 2. Are these real missing values? — diagnosis

**The data is NOT lost — it is a Summary-population gap.** Four independent observations:

1. **The per-position depth files exist for all 15.** Each affected sample has its
   `samtools/<sample>.<ref>.cand1.nodup.tsv` (`samtools depth -a` output, one row per reference position).
   The breadth and mean depth recompute trivially from these files. For the five single-infection rows the
   recomputed values **exactly match** both the values the Summary reports for unaffected samples and the
   previous run's published values (e.g. ERR1810473 = 0.78 % @≥10×, 2.8× depth; ERR1810454 = 89.00 %,
   16.1×; ERR1810480 = 57.32 %, 10.6×). So nothing is wrong with the depth computation itself.

2. **It is not a coverage / biology effect.** The affected set includes **sim22asingle (692,053 reads)** and
   **sim23asingle (695,614 reads)** — two of the most deeply and completely covered samples in the whole run
   (~100 % breadth). A low-coverage explanation is therefore excluded.

3. **It is not reference-specific.** Every reference in the affected set is populated normally for *other*
   samples that used the same reference: `3a_D17763` is blank for ERR1810473/477/480/481/500/501/sim23a but
   populated for ERR1810469 (4.18 %) and ERR1810505 (98.91 %); `2a_AB047639`, `1a_AF009606`, `1a_EF407457`,
   `4d_EU392172`, `2a_D00944` all likewise appear in both groups. So it is not a per-reference lookup
   failure (e.g. a missing reference-length entry).

4. **It is not call-type-specific.** Affected samples span monoinfection and co-infection, real and
   simulated, and 4 genotypes (1a, 2a, 3a, 4d).

**Where the break is.** The execution trace shows `SAMTOOLS_DEPTH` ran for the affected samples (e.g.
`JOINT_MAPPING:SAMTOOLS_DEPTH (ERR1810480 / ERR1810501 / sim22asingle)`), and the `.cand1.nodup.tsv`
outputs are present. The depth is therefore computed but **not summarised/joined into Summary.csv** for these
15 rows. The failure is downstream of `SAMTOOLS_DEPTH`, in the depth-aggregation → Summary join (the step
that turns the per-position `.nodup.tsv` into the `*_avg_depth` / `*_cov_breadth_min_*` columns).

---

## 3. Likely cause (hypotheses — not yet proven)

Because the gap is per-sample, all-or-nothing on the depth block, and independent of reference / coverage /
call type, the most probable explanations are:

- **An output-completeness / copy race.** The 2026-06-25 handoff already noted that "the result files were
  still being copied to the output dir … per-position depth files mostly present." If the depth-summary
  consumes the `.nodup.tsv` while results are still being staged/copied (or Summary.csv is assembled before
  every per-sample depth-summary has landed), a subset of rows would be written with the block empty while
  the raw `.nodup.tsv` still arrives. This fits the random-looking per-sample pattern.
- **A silent per-sample failure in the depth-parse/aggregation step** (empty/short read of the `.nodup.tsv`,
  a transient I/O error, or a parse that yields zero rows and is written as `NA` rather than retried),
  which would also produce an all-or-nothing empty block without affecting any other column.

What it is **not**: a logic error keyed on genotype/reference/coverage (ruled out in §2), and not data loss.

### Suggested checks for the dev team

1. Pull the **2026-06-27 run's** nextflow log / `execution_trace` (the traces currently in `pipeline_info/`
   are dated 2026-06-23/25 and may predate this Summary) and inspect the **depth-summarisation / coverage-parse
   process** (the one consuming `*.cand1.nodup.tsv`) for the 15 samples: did it run, exit 0, and emit a
   non-empty per-sample depth-summary?
2. Confirm whether Summary.csv assembly can begin **before** all per-sample depth-summaries are present
   (ordering / `collect()` completeness) — i.e. test the copy-race hypothesis.
3. Make the join **fail loudly** rather than silently writing `NA`: if a per-sample depth-summary is empty
   or missing when a consensus exists and reads > 0, raise an error (or at minimum a `review_flag` /
   log warning) so an incomplete Summary is visible at runtime instead of discovered downstream.
4. Check the **minor** depth/breadth block (`Minor_avg_depth`, `Minor_cov_breadth_min_*`) for the same gap
   on co-infection samples; this handoff characterised the major block.

---

## 4. Recovery / workaround (already applied for the manuscript)

The values are fully recoverable. Breadth at depth *d* = (positions with depth ≥ *d*) / (reference length),
mean depth = mean of column 3, computed directly from `samtools/<sample>.<ref>.cand1.nodup.tsv`:

```
awk -F'\t' '{n++; s+=$3; if($3>=10)b10++} END{printf "br10=%.2f avgdepth=%.1f\n",100*b10/n,s/n}' \
  samtools/ERR1810480.3a_D17763.cand1.nodup.tsv
```

For the manuscript, the five affected **Table 2** single-infection rows (454, 473, 477, 480, 481) were
recomputed this way and inserted; their recomputed values match the prior run exactly. The remaining
affected samples are in Table 5 (co-infection/IVT), where the affected value was either not displayed
(minor reported as "—") or independently re-verified (e.g. ERR1810501 major 3a breadth recomputed to
98.9 %, matching the published Table 5 value), and the two simulated singles (sim22a/sim23a) report
MAFFT-based accuracy in Table 1, which does not use these columns. **Manuscript impact is therefore already
mitigated**; the fix is needed for the pipeline output itself, so that Summary.csv is complete without a
manual recompute step.
