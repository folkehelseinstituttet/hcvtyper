# Phase 6: Neutral Candidate Selection - Pattern Map

**Mapped:** 2026-06-12
**Files analyzed:** 9 (8 modified + 1 new test)
**Analogs found:** 9 / 9 (all analogs are in-file or sibling — this is a refactor, not greenfield)

> Note: This phase is a re-wiring of existing code, not new construction. For most files the "closest analog" is the file's own current implementation (the pattern to *generalize*) or a sibling declaration in the same file. Excerpts below are the exact pattern to copy/extend, with line references verified against the repo on 2026-06-12.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `bin/summarize_mapping_to_all_references.R` | utility (R selection script) | transform (idxstats/depth → ranked CSV + FASTA) | self (current ranking/guards lines 76-189) | exact (self-generalize) |
| `bin/genotype_utils.R` | utility (sourced R helper) | transform (subtype → genotype) | reuse unchanged | reuse-only |
| `modules/local/parsefirstmapping/main.nf` | module (Nextflow process) | transform / emit | self (current emit contract lines 19-23) | exact (self-extend) |
| `workflows/hcvtyper.nf` (routing ~371-490) | route (channel topology) | event-driven fan-out | self (current major route lines 383-403) | exact (self-generalize) |
| `subworkflows/local/targeted_mapping/` | subworkflow (uniform mapping) | request-response (per-ref map+stats) | self — **no rewrite**, honor join discipline (lines 29-33, 49-56) | reuse-as-is |
| `conf/modules_hcv.config` (ext.prefix + param) | config | n/a | sibling `MAJOR/MINOR_MAPPING:*` closures (lines 125-310) + `denovo_*` param block (lines 13-21) | exact (sibling) |
| `nextflow_schema.json` (new param) | config (typed param) | n/a | `denovo_*` typed params (lines 191-216) | exact (sibling) |
| `nextflow.config` (new param) | config (param default) | n/a | `denovo_*` params block (lines 45-50) | exact (sibling) |
| `bin/tests/test_candidate_selection.R` | test (R subprocess) | n/a | `bin/tests/test_major_gate.R` (glob-picked by `run_all.sh`) | role-match |
| `modules/local/parsefirstmapping/tests/main.nf.test` (update) | test (nf-test) | n/a | self (GATE-02/03/04 cases lines 76-188) | exact (self-update) |

---

## Pattern Assignments

### `bin/summarize_mapping_to_all_references.R` (utility, transform) — PRIMARY REWRITE

**Analog:** self — generalize current logic; preserve guards and shim contract (D-01..D-03, D-05, D-06).

**Source helper pattern (reuse unchanged, lines 6-10):**
```r
# Relative path: staged into workdir as a `path(genotype_utils)` process input.
# Do NOT use absolute/projectDir — breaks container portability.
source("genotype_utils.R")
```

**Arg-parsing + NA-safe coercion pattern to extend (lines 12-26):** add a new positional `n_candidates` arg after `minCov`, coerced via `as.numeric`/`as.integer` exactly like `minRead`/`minCov`:
```r
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) { stop("Usage: ... <minRead> <minCov>", call.=FALSE) }
minRead    <- as.numeric(args[5])   # as.numeric => non-numeric yields NA, fails gate safely
minCov     <- as.numeric(args[6])
# ADD: n_candidates <- as.integer(args[7]) with default 2 if length(args) < 7
```

**Current ranking — the asymmetry being REPLACED (lines 76-104, 136-162):**
```r
summary <- df %>% group_by(Subtype, Genotype) %>%
  summarise(reads = sum(X3)) %>% arrange(desc(reads))
major_subtype <- summary$Subtype[1]                       # MAJOR = most reads
major_ref <- df %>% filter(Subtype == major_subtype) %>%
  arrange(desc(X3)) %>% head(n = 1) %>% pull(X1)           # top ref in that subtype
# MINOR by COVERAGE BREADTH among VALID different-genotype refs (D-03 removes this asymmetry):
tmp <- candidates %>% rowwise() %>%
  filter(is_valid_minor(cur_data())) %>% ungroup() %>%
  arrange(desc(percent_gt_4)) %>% slice(1)
```
**Target (neutral, D-01/D-02/D-03/D-05):** group by Subtype → top ref per subtype by reads → arrange subtypes by total reads → `head(n = n_candidates)`. **Delete `is_valid_minor()` entirely** (lines 112-129) per D-05. Emit one long-format row per selected candidate.

**Validity rule block to DELETE (lines 112-129)** — moves to Phase 8:
```r
is_valid_minor <- function(minor_row) { ... 1a/1b allow ... 2k1b block ... different-genotype ... }
```

**Empty-frame / no-mapping guards to GENERALIZE, not reinvent (lines 63-70, 136-148):**
```r
df_final$minor_call[1] <- "no"
df_final$gate_flag[1]  <- "no_mapping"            # default carried always
if (nrow(df) > 0) { ... }                          # whole populated branch guarded
if (nrow(candidates) > 0) { ... } else { tmp <- candidates %>% slice(0) }  # zero-row rowwise guard
```

**Legacy 10-column wide CSV — must be RECONSTRUCTED under shim (lines 57-58), D-06:**
```r
colnames(df_final) <- c("sample","total_mapped_reads","major_ref","major_reads",
  "major_cov","minor_ref","minor_reads","minor_cov","minor_call","gate_flag")
```
**Shim rule:** rank 1 → `major_*`; rank 2 → `minor_*`; emit the NEW long-format table alongside (suggested cols: `sample,candidate_rank,candidate_ref,candidate_subtype,candidate_genotype,candidate_reads,candidate_cov,confirmation_status`). Always emit all 10 legacy columns NA-filled when no 2nd candidate (summarize.R pulls them unconditionally — see Open Q3 in RESEARCH).

**Gate-decision threshold logic to GENERALIZE per-candidate (lines 170-177):**
```r
major_pass <- (df_final$major_reads[1] > minRead) && (df_final$major_cov[1] > minCov)
minor_pass <- length(minor_ref) > 0 && !is.na(df_final$minor_reads[1]) &&
              (df_final$minor_reads[1] > minRead) && (df_final$minor_cov[1] > minCov)
df_final$minor_call[1] <- if (isTRUE(major_pass) && isTRUE(minor_pass)) "yes" else "no"
df_final$gate_flag[1]  <- if (!isTRUE(major_pass)) "major_below_threshold" else "ok"
```
**Target:** compute a per-candidate `confirmation_status` from the same `reads>minRead && cov>minCov` comparison; STILL reconstruct `minor_call`/`gate_flag` for the shim.

**FASTA write — rank-keyed suffix, keep guard (lines 184-189):**
```r
fasta <- read.fasta(file = references)
write.fasta(sequences = fasta[major_ref], names = major_ref, file.out = paste0(sampleName, ".", major_ref, "_major.fa"))
if (length(minor_ref) > 0) {
  write.fasta(sequences = fasta[minor_ref], names = minor_ref, file.out = paste0(sampleName, ".", minor_ref, "_minor.fa"))
}
```
**Target:** cand_1 → `_major.fa`, cand_2 → `_minor.fa`; keep `length(...) > 0` guard (Pitfall 3 — never write `_minor.fa` on empty ref).

---

### `bin/genotype_utils.R` (utility, sourced helper) — REUSE UNCHANGED

**Analog:** self — read-only reuse. The canonical 2k1b-aware rule (lines 24-26):
```r
genotype_from_subtype <- function(subtype) {
  if_else(subtype == "2k1b", subtype, substr(subtype, 1, 1))
}
```
**Apply to:** subtype→genotype derivation in the new ranking/dedup. Do NOT hand-roll inline `substr`/`if` (Don't Hand-Roll). Already staged as `path(genotype_utils)` in PARSEFIRSTMAPPING.

---

### `modules/local/parsefirstmapping/main.nf` (module, transform/emit) — ADD CANDIDATES EMIT

**Analog:** self — extend the emit contract; keep all existing emits (D-06).

**Current emit contract to EXTEND (lines 19-23):**
```groovy
output:
tuple val(meta), path("*.csv"), path("*major.fa"), emit: major_mapping, optional: true
tuple val(meta), path("*.csv"), path("*minor.fa"), emit: minor_mapping, optional: true
tuple val(meta), path("*.csv"),                    emit: csv,           optional: true
path "versions.yml",                               emit: versions
```
**Target:** ADD a `candidates` emit for the long-format CSV (e.g. `path("*.candidates.csv")`), keep `major_mapping`/`minor_mapping`/`csv` for the shim. Differentiate filename globs so `*.csv` (legacy wide) and `*.candidates.csv` (long) don't collide.

**Script invocation to extend (lines 32-40):** pass the new param positionally after `${params.minCov}`:
```groovy
summarize_mapping_to_all_references.R \\
    ${idxstats} ${depth} ${prefix} ${references} \\
    ${params.minRead} ${params.minCov} \\
    $args
# ADD: ${params.n_candidates} as a new positional arg
```

**Stub contract to update (lines 57-67):** the stub must emit the new `*.candidates.csv` and full long-format header so `-stub-run` of the workflow fans out correctly. Mirror the existing 10-col legacy stub header.

---

### `workflows/hcvtyper.nf` (route, event-driven fan-out) — COLLAPSE TO ONE FAN-OUT (D-04)

**Analog:** self — generalize the current major route into a per-candidate `splitCsv` fan-out.

**Current major route — the splitCsv-into-meta pattern to GENERALIZE (lines 383-403):**
```groovy
ch_major_mapping = PARSEFIRSTMAPPING.out.major_mapping.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    .map { meta, _csv, major_fasta, _reads ->
        def elements = _csv.splitCsv( header: true, sep:',')
        def new_meta = meta + elements[0]             // lifts row[0] into meta (single-row today)
        assert new_meta.id == new_meta.sample : "Metadata mismatch: id=${new_meta.id}, sample=${new_meta.sample}"
        tuple(new_meta, major_fasta, _reads)
    }
    .filter { entry ->
        def mappedReads = entry[0]['major_reads'].toInteger()
        def majorCov = entry[0]['major_cov'].toInteger()
        mappedReads > params.minRead && majorCov > params.minCov
    }
MAJOR_MAPPING( ch_major_mapping )
```
**Current minor route — routes on R-emitted STRING, not numeric coercion (lines 419-436):**
```groovy
.map { meta, _csv, minor_fasta, _reads ->
    def elements = _csv.splitCsv( header: true, sep:',')
    def new_meta = meta + elements[0]
    assert new_meta.id == new_meta.sample : "..."
    tuple(new_meta, minor_fasta, _reads)
}
.filter { entry -> entry[0]['minor_call'] == 'yes' }   // string route avoids NA.toInteger() crash
MINOR_MAPPING ( ch_minor_mapping )
```
**Target (Pattern 1 + 2):** `splitCsv(header:true)` the **long-format** candidates CSV → **iterate ALL rows** (NOT `elements[0]`; current code reads `[0]` only because the file is single-row today — Assumption A2) → one channel element per candidate → `meta += row` carrying `candidate_rank`, `candidate_ref`, `confirmation_status` → `.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)` by `meta.id` → ONE `TARGETED_MAPPING` call. Keep the `id == sample` assert. Route/carry decisions as R-emitted STRINGS (never `NA.toInteger()`). Make per-candidate meta unique (`reference` + `candidate_rank`) before entering the subworkflow.

**Downstream feed channels to UPDATE IN LOCKSTEP (T-2, lines 443, 468-490):**
```groovy
HCVGLUE ( MAJOR_MAPPING.out.aligned.collect({it[1]}).mix(MINOR_MAPPING.out.aligned.collect({it[1]})).collect(), ... )
ch_stats_withdup = MAJOR_MAPPING.out.stats_withdup.collect({it[1]}).mix(MINOR_MAPPING.out.stats_withdup.collect({it[1]}))
ch_depth         = MAJOR_MAPPING.out.depth.collect({it[1]}).mix(MINOR_MAPPING.out.depth.collect({it[1]}))
ch_variation     = MAJOR_MAPPING.out.variation.collect().mix(MINOR_MAPPING.out.variation.collect())
ch_consensus_distance = MAJOR_MAPPING.out.consensus_distance.collect({it[1]}).mix(MINOR_MAPPING.out.consensus_distance.collect({it[1]}))
```
**Target:** replace every `.mix(MAJOR..., MINOR...)` pair with the single `TARGETED_MAPPING.out.<x>` (the fan-out already contains all candidate BAMs/stats). **Missing any one silently drops half the stats** (T-2 / Assumption A3).

---

### `subworkflows/local/targeted_mapping/` (subworkflow, request-response) — NO REWRITE, HONOR JOIN DISCIPLINE

**Analog:** self — reuse as-is; the uniform path candidates fan into.

**Meta-enrichment + metadata-key join discipline to HONOR (lines 29-33, 49-56), Pattern 2:**
```groovy
ch_input = ch_major_mapping.map { meta, fasta, reads ->
    def new_meta = meta + [ reference: fasta.getBaseName().toString().split('\\.').last() ]
    tuple(new_meta, fasta, reads)
}
...
ch_aligned_input = ch_input.reads
    .join( BOWTIE2_BUILD.out.index )   // by meta key — BUILD emits in completion order, NOT submission order
    .join( ch_input.build )            // positional join would mis-pair under parallelism
```
**Action:** ensure per-candidate meta is unique BEFORE entering this subworkflow (add `candidate_rank`). With N candidates per sample sharing `meta.id`, an id-only key cross-pairs (Pitfall 2). The `reference` enrichment here is in scope for `ext.prefix` (Assumption A1) — this is why `candidate_rank` can drive the slot literal.

---

### `conf/modules_hcv.config` (config) — SLOT-FROM-RANK ext.prefix + new param (T-1, HIGH)

**Analog:** sibling `MAJOR/MINOR_MAPPING:*` closures — collapse the two alias-keyed selectors into one `candidate_rank`-keyed closure.

**Current alias-keyed slot closures to REPLACE (T-1 root cause, lines 125-133, 199-205, 242-310):**
```groovy
withName: '...:MAJOR_MAPPING:SAMTOOLS_SORMADUP' { ext.prefix = { "${meta.id}.${meta.reference}.major.nodup" } }
withName: '...:MINOR_MAPPING:SAMTOOLS_SORMADUP' { ext.prefix = { "${meta.id}.${meta.reference}.minor.nodup" } }
withName: '...:MAJOR_MAPPING:STATS_WITHDUP'      { ext.prefix = { "${meta.id}.${meta.reference}.major.withdup" } }
withName: '...:MINOR_MAPPING:STATS_WITHDUP'      { ext.prefix = { "${meta.id}.${meta.reference}.minor.withdup" } }
withName: '...:MAJOR_MAPPING:SAMTOOLS_DEPTH'     { ext.prefix = { "${meta1.id}.${meta1.reference}.major.nodup" } }
withName: '...:MINOR_MAPPING:SAMTOOLS_DEPTH'     { ext.prefix = { "${meta1.id}.${meta1.reference}.minor.nodup" } }
withName: '...:MAJOR_MAPPING:IVAR_CONSENSUS'     { ext.prefix = { "${meta.id}.major" }; ... }
withName: '...:MINOR_MAPPING:IVAR_CONSENSUS'     { ext.prefix = { "${meta.id}.minor" }; ... }
withName: '...:MAJOR_MAPPING:CONSENSUS_DISTANCE' { ext.prefix = { "${meta.id}.major" } }
withName: '...:MINOR_MAPPING:CONSENSUS_DISTANCE' { ext.prefix = { "${meta.id}.minor" } }
```
**Target (Pattern 3) — single closure per process, slot derived from `meta.candidate_rank`:**
```groovy
withName: '...:TARGETED_MAPPING:STATS_WITHDUP' {
    ext.prefix = { "${meta.id}.${meta.reference}.${meta.candidate_rank == 1 ? 'major' : 'minor'}.withdup" }
}
```
**CRITICAL (Pitfall 1):** This is a FIRST-CLASS TASK, easy to miss because it lives in config. If omitted, mapping outputs lose the `.major.`/`.minor.` slot, `summarize.R`'s `str_split(...)[[1]][3]` (lines 226, 291) returns `withdup`/`nodup`, and `Major_reference`/`Minor_reference` go blank/swapped. Note `SAMTOOLS_DEPTH` uses `meta1.` (joined-input meta) not `meta.` — preserve that namespace when rewriting (lines 244, 249). **For N>2, cand_k≥3 has no legacy slot** — derive `.cand{k}.` and document N>2 as Phase-9-dependent (T-1 residual / Open Q1).

**New param to ADD to the `params{}` block (sibling pattern, lines 13-21):**
```groovy
params {
    minRead = 499
    minCov = 29
    denovo_match_level = 'genotype'
    denovo_confirm_minor = true
    // ADD: n_candidates = 2
}
```

---

### `nextflow_schema.json` (config, typed param) — DECLARE n_candidates

**Analog:** `denovo_*` typed params (lines 191-216).

**Pattern to mirror (lines 191-211):**
```json
"denovo_min_contig_length": { "type": "integer", "default": 500, "description": "..." },
"denovo_match_level": { "type": "string", "default": "genotype", "enum": ["genotype","subtype"], "description": "..." }
```
**Target declaration:**
```json
"n_candidates": {
  "type": "integer",
  "default": 2,
  "minimum": 1,
  "description": "Number of neutrally-ranked candidate references to select and map (cand_1..cand_n). Default 2 reproduces the legacy two-slot topology under the compatibility shim."
}
```

---

### `nextflow.config` (config, param default) — DECLARE n_candidates DEFAULT

**Analog:** `denovo_*` params block (lines 45-50).

**Pattern to mirror (lines 45-50):**
```groovy
// De novo confirmation options
denovo_min_contig_length   = 500
denovo_match_level         = 'genotype'
denovo_confirm_minor       = true
```
**Target:** add `n_candidates = 2` with an inline comment, alongside `minRead`/`minCov` or in a new "Candidate selection options" group. Declare in all three sites (schema + nextflow.config + modules_hcv.config) per the `denovo_*` triple precedent.

---

### `bin/tests/test_candidate_selection.R` (test, NEW) — AUTO-GLOBBED

**Analog:** `bin/tests/test_major_gate.R` (subprocess-contract test, picked up by `bin/tests/run_all.sh` glob `test_*.R`).

**Coverage:** REFSEL-01 ranking (top-N distinct subtypes by reads), distinct-subtype dedup (preserves 1a/1b), shim column reconstruction (rank1→major_*, rank2→minor_*), NA-fill on single-candidate. Run via `bash bin/tests/run_all.sh` (seconds, no docker) — preferred during iteration given disk-near-full.

---

### `modules/local/parsefirstmapping/tests/main.nf.test` (test, UPDATE) — T-3

**Analog:** self — existing GATE-02/03/04 cases (lines 76-188).

**Cases that CHANGE under D-05 (no validity filtering in selection):**
- GATE-03 (lines 76-106): 2k1b single-major-no-minor — a 2k1b sample may now select a 2nd candidate; assertion `minor_mapping.size() == 0` (line 97) will change.
- GATE-04 (lines 110-138): no-minor-FASTA — may now write a 2nd FASTA; update accordingly.

**Cases that LARGELY HOLD under the shim:**
- GATE-02 (lines 147-188): failing-major reports stats, no crash — the direct CSV-cell assertions (`gate_flag`, `minor_call`, `major_ref`...) still hold via the reconstructed shim columns.

**Add:** assertions for the new `candidates` emit + long-format CSV header; regenerate `.snap`. Snapshot churn is EXPECTED, not a regression (Assumption A4).

---

## Shared Patterns

### R-emits-decision, Nextflow-routes
**Source:** `workflows/hcvtyper.nf:383-403` (splitCsv→meta), `bin/summarize_mapping_to_all_references.R:170-177` (R emits gate strings)
**Apply to:** All candidate routing. R ranks & writes routing columns; Nextflow `splitCsv(header:true)` fans out per row. Route on R-emitted STRINGS, never Groovy `.toInteger()` on possibly-NA fields.

### Metadata-key joins, never positional
**Source:** `subworkflows/local/targeted_mapping/main.nf:29-33, 49-56`
**Apply to:** Every join carrying per-candidate data through TARGETED_MAPPING. Make meta unique per candidate (`reference` + `candidate_rank`) or N-per-sample candidates cross-pair (Pitfall 2).

### Filename slot derived from candidate_rank (the shim hinge)
**Source:** target Pattern 3 (replaces alias-keyed `conf/modules_hcv.config:125-310`); consumed by `bin/summarize.R:226,291`
**Apply to:** Every `TARGETED_MAPPING:*` `ext.prefix` closure that today distinguished major/minor by alias name. `candidate_rank == 1 ? 'major' : 'minor'`. Preserves `summarize.R` unchanged (D-06).

### Typed-param triple
**Source:** `nextflow_schema.json:191-211` + `nextflow.config:45-50` + `conf/modules_hcv.config:13-21` (`denovo_*` precedent)
**Apply to:** `n_candidates` — declare in all three; pass positionally to the R script (like `minRead`/`minCov`).

### NA-safe / empty-frame guards (reuse, don't reinvent)
**Source:** `bin/summarize_mapping_to_all_references.R:21-26` (as.numeric coercion), `63-70` (default gate_flag + nrow guard), `136-148` (zero-row rowwise guard)
**Apply to:** generalizing selection to a variable-length candidate set; preserve default `gate_flag="no_mapping"` and the `length(...) > 0` FASTA-write guard (Pitfall 3).

---

## No Analog Found

None. Every file in this phase modifies or sits beside existing code; all analogs are in-file (self-generalize) or sibling declarations. This is a re-wiring phase — the algorithms and idioms already exist (see RESEARCH "Don't Hand-Roll": every "new" capability is a reusable asset).

---

## Decision Tensions Carried Forward (for planner)

- **T-1 (HIGH):** D-04 (collapse aliases) destroys the alias-name signal that produces `.major.`/`.minor.`; D-06 requires those exact strings. Resolution: derive slot from `meta.candidate_rank` in single `TARGETED_MAPPING:*` `ext.prefix` closures. **Planner must make the `modules_hcv.config` ext.prefix rewrite a first-class, explicitly-listed task** (config, not code — easy to miss; silent Summary.csv corruption if omitted).
- **T-2 (MEDIUM):** HCVGLUE + SUMMARIZE feed channels `.mix(MAJOR..., MINOR...)` (hcvtyper.nf:443, 468-490) must be replaced with the single `TARGETED_MAPPING.out.*` in lockstep or stats/depth/variation/consensus drop to half.
- **T-3 (LOW):** parsefirstmapping nf-test GATE-03/04 assertions change under no-validity-filtering; snapshots regenerated (expected, not a regression).
- **D-05 CONSEQUENCE:** Do NOT treat ROADMAP success-criterion #5 (golden reproduction) as a Phase-6 gate. Phase 6 verifies topology + no-crash + shim integrity; golden reproduction is gated at Phase 9.

---

## Metadata

**Analog search scope:** `bin/`, `bin/tests/`, `modules/local/parsefirstmapping/`, `subworkflows/local/targeted_mapping/`, `workflows/hcvtyper.nf`, `conf/modules_hcv.config`, `nextflow_schema.json`, `nextflow.config`
**Files scanned (read in full or targeted range):** 9
**Pattern extraction date:** 2026-06-12
