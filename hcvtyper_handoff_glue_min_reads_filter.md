# Handoff: `glue_min_reads` — pre-GLUE candidate BAM filter

**Branch:** `dev`
**Type:** Efficiency optimization (no correctness change)
**Author of spec:** (analysis session, 2026-06-25)

---

## 1. Problem / rationale

JOINT_MAPPING is now **competitive**: all candidate references for a sample are
concatenated into one combined reference, one Bowtie2 index is built, and each read
is aligned **once** and assigned to its single best location
(`subworkflows/local/joint_mapping/main.nf:71-91`). Reads are therefore *partitioned*
across candidates instead of double-counted as they were under the old independent
per-candidate TARGETED_MAPPING.

A direct consequence: a candidate that looked supported on first-mapping evidence can
be out-competed for nearly all its reads and end up with very few post-dedup reads.

Today **every** candidate that clears the upstream `confirmation_status == 'pass'`
gate (set in PARSEFIRSTMAPPING/RESCUE_EVALUATION, i.e. *before* competitive mapping)
emits a per-candidate BAM, and **all** of those BAMs are sent to HCV-GLUE with no
intervening read-count / coverage filter:

```groovy
// workflows/hcvtyper.nf:472-475
HCVGLUE (
    JOINT_MAPPING.out.aligned.collect({it[1]}).collect(), // every passing candidate BAM
    params.hcvglue_threshold
)
```

`run_hcvglue.sh` then loops over every staged `*.bam` and runs **two** serial docker
invocations each — `reportBam` (JSON) + `reportBamAsHtml` — inside a single
gluetools-mysql container (`bin/run_hcvglue.sh:135-170`). So each near-empty BAM
costs two full container runs of wall-clock time and produces nothing useful (GLUE's
internal `hcvglue_threshold=15` coverage cutoff means a starved BAM reports nothing,
but the BAM is still loaded and processed — the threshold does not save the compute).

**Goal:** drop candidate BAMs that obviously cannot yield a resistance call before
they reach GLUE.

---

## 2. Governing constraint (READ FIRST)

This must remain a **pure efficiency optimization**. The role classifier in
`bin/summarize.R` / `bin/classify_roles.R` runs **after** GLUE and can legitimately
keep a low-yield minor as a reported `co-infection` (the `uncorroborated_kept` branch,
`classify_roles.R:384-388`). If the pre-GLUE filter prunes such a BAM, that genuine
co-infection silently loses its resistance call — a correctness regression.

Therefore:

- The filter threshold **must stay well below** the role-classification floor
  (production `minRead = 499`, `minCov = 29` in `conf/modules_hcv.config:16-17`).
  `glue_min_reads = 100` is ~5× below `minRead`, so no reportable strain is at risk.
- The filter applies **only to the GLUE leg**. The channels feeding SUMMARIZE
  (`idxstats_withdup`, `idxstats_nodup`, `depth`, `consensus`, `consensus_distance`,
  `variation`) must continue to see **all** candidates — do not filter those.

Rationale for value 100: HCV genome ≈ 9.6 kb; ~150 bp paired reads ≈ 300 bp/pair, so
~64 reads is mean 1× *with perfect spread*. Below ~100 reads a candidate cannot
sustain ≥15% depth across the NS3/NS5A/NS5B drug-target genes, so GLUE returns empty.

---

## 3. Changes

### 3.1 `nextflow.config` — new param

Add to the `params { ... }` block (near `hcvglue_threshold`, currently line 21):

```groovy
glue_min_reads = 100 // Minimum post-dedup mapped reads for a candidate BAM to be submitted to HCV-GLUE. Candidates starved by competitive mapping below this are skipped (resistance analysis only — does NOT affect SUMMARIZE / strain classification).
```

If the pipeline has an nf-schema (`nextflow_schema.json`), add a matching entry so
`nf-core lint` / schema validation passes (integer, default 100, minimum 0). See the
existing `hcvglue_threshold` schema entry for the pattern.

### 3.2 `subworkflows/local/joint_mapping/main.nf` — attach per-candidate nodup read count to each BAM's meta

**Requirement:** every emitted `aligned` BAM (= `ch_percand`, one per **candidate**,
multiple per sample) must carry its OWN post-dedup mapped read count in the meta
(propose `meta.candidate_nodup_reads`). Per-sample is NOT sufficient.

The counts already exist: `IDXSTATS_NODUP` runs on the **combined** dedup BAM and
emits one idxstats file per sample with **one row per candidate reference**
(`joint_mapping/main.nf:134-137`). idxstats TSV columns: **col 1 = reference name,
col 2 = ref length, col 3 = mapped reads, col 4 = unmapped reads** (`*` row = unmapped
bin; ignore it).

Each per-candidate element in `ch_split` / `ch_percand` already carries
`meta.candidate_ref` (the reference name for that candidate). Join each candidate to
its own count by matching `candidate_ref` against idxstats col 1.

Suggested implementation (place after `ch_percand` is defined, ~line 189, before the
`emit:` block; adapt names to local style):

```groovy
// Build a per-sample lookup of nodup mapped reads keyed by reference name, then
// attach each candidate's OWN count to its meta (multiple candidates per sample).
// IDXSTATS_NODUP emits one file per sample (meta keyed by sampleName/id), one row per
// candidate reference: col1 = ref name, col3 = mapped reads. Drop the '*' unmapped bin.
ch_percand_counted = ch_percand
    .map { meta, bam -> tuple(meta.subMap('id'), meta, bam) }            // key by sample id
    .combine(
        IDXSTATS_NODUP.out.idxstats.map { meta, idx -> tuple(meta.subMap('id'), idx) },
        by: 0
    )
    .map { _key, meta, bam, idx ->
        def counts = [:]
        idx.readLines().each { line ->
            def f = line.split('\t')
            if (f.size() >= 3 && f[0] != '*') counts[f[0]] = (f[2] as long)
        }
        def n = counts.getOrDefault(meta.candidate_ref as String, 0L)
        tuple(meta + [candidate_nodup_reads: n], bam)
    }
```

> Implementation notes:
> - Join by **sample key only** (`meta.subMap('id')`), because `IDXSTATS_NODUP`'s meta
>   is the sample-level meta whereas `ch_percand`'s meta is the per-candidate meta —
>   they will not be equal, so a plain `.join()` on full meta would drop everything.
>   `.combine(..., by: 0)` on the sample-id sub-map fans the single per-sample idxstats
>   file out to each of that sample's candidate BAMs.
> - `candidate_ref` must match idxstats col 1 **exactly**. Confirm whether the combined
>   FASTA headers (and therefore the BAM `@SQ` / idxstats names) are the raw reference
>   names or were altered by `REHEADER_CANDIDATE` / `CAT_CANDIDATES`. If they differ,
>   normalize before the lookup. **Verify with a real run's `.nodup.idxstats`.**
> - Read the idxstats from the **combined** BAM (post-dedup), which is what
>   `IDXSTATS_NODUP` already produces — do not re-derive per-candidate.

Then change the emit to use the counted channel:

```groovy
emit:
aligned            = ch_percand_counted          // per-candidate dedup BAM, meta carries candidate_nodup_reads
...
```

Leave all SUMMARIZE-bound emits (`idxstats_withdup`, `idxstats_nodup`, `depth`,
`consensus`, `consensus_distance`, `variation`) **unchanged**.

### 3.3 `workflows/hcvtyper.nf` — filter the GLUE leg only

Replace the HCVGLUE input (currently lines 472-475) so the collected BAM list is
filtered by the new meta field:

```groovy
if (!params.skip_hcvglue) {
    ch_glue_bams = JOINT_MAPPING.out.aligned
        .filter { meta, _bam -> (meta.candidate_nodup_reads ?: 0) >= params.glue_min_reads }
        .collect({ it[1] })
        .collect()

    HCVGLUE (
        ch_glue_bams,
        params.hcvglue_threshold
    )
    ...
}
```

Do **not** touch the SUMMARIZE staging block (`hcvtyper.nf:485+`) — it must keep
seeing every candidate via the idxstats/depth/consensus channels.

> Edge case: if a sample's candidates are ALL below `glue_min_reads`, the collected
> list can be empty. `HCVGLUE` already tolerates a no-BAM input (`run_hcvglue.sh`
> handles the empty-glob case and the outputs are `optional: true`), and GLUE_PARSER
> reads collected JSONs, so an empty set is safe — but confirm the
> `HCVGLUE.out.GLUE_json.collect()` → GLUE_PARSER path does not stall when zero JSONs
> are produced (should be fine since outputs are optional, but worth a sanity check).

---

## 4. Verification

1. **Unit/channel sanity:** run the `test` profile (`minRead=100, minCov=5`) and a
   `test_full` run; confirm pipeline completes.
2. **Meta correctness:** dump `meta.candidate_nodup_reads` for a multi-candidate
   sample and cross-check each value against the corresponding row in that sample's
   `*.nodup.idxstats` (col 3). They must match per candidate, not per sample.
3. **Filter behaviour:** pick a sample with a known starved candidate (few competitive
   reads); confirm its BAM is absent from the GLUE work dir / produces no GLUE JSON,
   while candidates ≥100 reads still do.
4. **No SUMMARIZE regression:** confirm the final summary TSV still has rows /
   coverage / role classification for the filtered-out candidate (it must still appear
   there via the idxstats/depth channels) — only its GLUE resistance columns should be
   empty/NA.
5. **nf-test:** if there is a JOINT_MAPPING or workflow-level nf-test snapshot,
   regenerate it; the new meta key and the filtered GLUE input will change snapshots.
   (Note from project memory: nf-test fills the host disk — clean `work/` + prune
   docker volumes if you hit ENOSPC.)

---

## 5. Files touched

- `nextflow.config` (+ `nextflow_schema.json` if present)
- `subworkflows/local/joint_mapping/main.nf`
- `workflows/hcvtyper.nf`
- nf-test snapshots (regenerate if affected)

No changes to `bin/summarize.R`, `bin/classify_roles.R`, `run_hcvglue.sh`, or the
`confirmation_status` gate.
