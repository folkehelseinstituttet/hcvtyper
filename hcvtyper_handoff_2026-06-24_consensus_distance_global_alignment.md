# Handoff: Global alignment for `consensus_distance.R` (`similarity_pct`)

**Date:** 2026-06-24
**Author:** jon.brate (via Claude Code investigation)
**Intended entry point:** `/gsd-quick` (run on the server/other machine — local host disk is near-full, container rebuild + nf-test fill it)
**Status:** investigated + decisions locked; NOT yet implemented

---

## TL;DR task for `/gsd-quick`

Replace the ungapped position-by-position comparison in `bin/consensus_distance.R` with a **Biostrings global pairwise alignment** (`pairwiseAlignment(type="global")`, Needleman-Wunsch). Count **indels (gap columns) as differences**; keep excluding N / zero-coverage columns from the denominator. Keep the output TSV columns **identical** so `summarize.R` and the summary schema are unaffected. Add `bioconductor-biostrings` to the `CONSENSUS_DISTANCE` container and regenerate the `consensus_distance` nf-test snapshot.

---

## Why (the problem)

`similarity_pct` is produced by the `CONSENSUS_DISTANCE` process and reported (per-candidate, Major=cand1 / Minor=cand2) in the final summary. It is currently **report-only** — it feeds no gating/confirmation/rescue logic.

The current script (`bin/consensus_distance.R`) does an **ungapped, position-by-position** identity comparison. Its correctness depends on the claim that the iVar consensus is *"in reference coordinates: one character per reference position."*

**That assumption is unsafe.** iVar is invoked as (`conf/modules_hcv.config:313-314`, JOINT_MAPPING:IVAR_CONSENSUS):

```
samtools mpileup --count-orphans --no-BAQ --max-depth 0 --min-BQ 0 -aa \
  | ivar consensus -t 0 -q 1 -m 10 -n N
```

- `-aa` only guarantees every reference position appears in the **pileup** — it does not stop iVar from emitting indels.
- `-t 0` (frequency threshold 0) means iVar will call insertions supported by ~any reads and represent deletions (`*`). iVar 1.4.x incorporates both into the consensus.

So the consensus can be **longer** (insertion) or internally **shorter** (deletion) than the reference. After the first internal indel, the ungapped loop is **frame-shifted for every downstream position** → ~75% of post-shift bases read as mismatches. A single true 3-bp indel can turn a real ~98% identity into a reported ~85%, with inflated `n_differences`. The existing `min()`-truncation only fixes R vector-recycling; it does nothing for internal indels. The code comment only anticipates *tail* trimming, not internal indels.

Impact today is "misleading the analyst," not mis-called strains (report-only). But it must be correct before `similarity_pct` could ever be used as confirmation evidence.

---

## Decisions (locked)

1. **Aligner:** Biostrings `pairwiseAlignment(type = "global")` (exact Needleman-Wunsch, stays in the existing R script, one container dependency add, no new process). Chosen over MAFFT (`--globalpair`) and EMBOSS `needle` — both add a tool/process/container/parsing for no benefit in the 2-sequence case.
2. **Indel semantics:** Gap columns count **as differences** (standard for genotyping distance). `n_differences` = substitutions + indel gap-columns.
3. **N / zero-coverage handling:** Still **exclude** columns where the consensus base is `N`/uncovered (and where the reference base is missing) from the denominator, so low-coverage samples are not penalized as divergent.
4. **Output contract unchanged:** Same TSV columns, so no downstream edits.

---

## Implementation guidance

### 1. `bin/consensus_distance.R`

- Swap `seqinr` for `Biostrings` (or keep `seqinr` only for FASTA reading if convenient, but `Biostrings::readDNAStringSet` is cleaner end-to-end).
- Read consensus + reference (first record of each). Preserve the consensus header name (`cons_name`) and reference name (`ref_name`) for the `sample` / `reference` output columns — `summarize.R` validates the `Consensus_<sample>...cand<rank>` pattern from the `sample` column, so keep that header intact (it comes from the IVAR consensus FASTA header; just pass it through).
- Compute `consensus_length` = number of non-`N`, non-gap called bases in the **full** consensus (before alignment), matching current semantics.
- Run a global NW alignment:
  ```r
  library(Biostrings)
  aln <- pairwiseAlignment(pattern = consensus, subject = reference, type = "global")
  cons_aln <- as.character(alignedPattern(aln))   # consensus row, with '-' gaps
  ref_aln  <- as.character(alignedSubject(aln))   # reference row, with '-' gaps
  ```
  (split to character vectors per column).
- Define per-column callability and differences:
  - A column is **excluded** (not callable) if the consensus char is `N`/`n` (uncovered) OR the reference char is `N`/`n`/missing. (Gaps are NOT excluded — they count, per decision 2.)
  - `alignment_length` = number of callable columns (i.e. total aligned columns minus excluded N/uncovered columns).
  - `n_differences` = callable columns where the two chars differ, **including** gap-vs-base columns (indels).
  - `similarity_pct = round((alignment_length - n_differences) / alignment_length * 100, 4)`.
  - Edge case: `alignment_length == 0` → `similarity_pct = NA_real_`, `n_differences = NA_integer_` (preserve current behavior).
- Decide scoring params explicitly and document them in a comment (e.g. `nucleotideSubstitutionMatrix(match=1, mismatch=-1)` or a simple match/mismatch, `gapOpening`/`gapExtension`). Pick values that don't over-fragment indels for ~10 kb HCV genomes; document the choice. Keep it deterministic.
- Keep the usage string and the 3-arg CLI (`<consensus.fa> <reference.fa> <output.tsv>`) and the same `write.table` output (tab-sep, no quotes, no row names) with columns in the same order:
  `sample, reference, similarity_pct, n_differences, alignment_length, consensus_length`.

### 2. Container — `modules/local/consensus_distance/environment.yml`

- Add `bioconductor-biostrings` (bioconda channel). Keep existing R deps as needed.
- Rebuild the Wave container and update the `container` line hash in `modules/local/consensus_distance/main.nf` (the digest URL + `community.wave.seqera.io/...` tag). Use the NEXTFLOW conda env's tooling / Wave as per project convention.
- The current container line:
  ```
  'community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab'
  ```
  becomes a new build that also includes `bioconductor-biostrings`.

### 3. nf-test snapshot

- Update/regenerate the `consensus_distance` module nf-test snapshot. The stub (`main.nf` `stub:` block) hardcodes a sample row — keep stub as-is unless the snapshot needs it.
- Use the locally pinned nf-test (`~/.nf-test`, 0.9.3) via PATH prefix + `--profile docker`.
- **Disk caution:** host disk is near-full; nf-test/docker can fill it → ENOSPC. Recovery: `sudo rm -rf work/`, `docker volume prune`. Consider running on the server.

---

## Files in scope

| File | Change |
|------|--------|
| `bin/consensus_distance.R` | Rewrite comparison → Biostrings global NW; indels = differences |
| `modules/local/consensus_distance/environment.yml` | Add `bioconductor-biostrings` |
| `modules/local/consensus_distance/main.nf` | Update Wave container hash; (stub unchanged) |
| `modules/local/consensus_distance/tests/*.snap` | Regenerate snapshot |

**Out of scope (no change needed):** `bin/summarize.R` (consumes the same columns), summary JSON/TSV schema, `subworkflows/local/joint_mapping/main.nf`, `subworkflows/local/targeted_mapping/main.nf`.

---

## Verification

1. `nf-test` for the `consensus_distance` module passes (`--profile docker`).
2. Hand-built sanity case: a consensus with a known 3-bp insertion vs reference should now report high `similarity_pct` (≈ matching identity) instead of the frame-shifted ~85%. Confirm `n_differences` ≈ indel length + true substitutions.
3. Edge cases preserved: all-N consensus → `NA` similarity; tail-trimmed consensus (shorter, no internal indel) → identity unchanged vs old behavior within rounding.
4. Output TSV column names/order unchanged; a full pipeline summary still populates `Major_consensus_similarity_pct` / `Minor_consensus_similarity_pct`.

---

## Context references

- Process: `modules/local/consensus_distance/main.nf`
- Script: `bin/consensus_distance.R` (ungapped loop at lines ~55-71)
- iVar invocation: `conf/modules_hcv.config:311-315`
- Active wiring: `subworkflows/local/joint_mapping/main.nf:239-254` → `workflows/hcvtyper.nf:537`
- Consumer: `bin/summarize.R:1089-1167` (parses `<sample>.cand{rank}.consensus_distance.tsv`)
