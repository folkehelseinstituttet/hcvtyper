# Phase 6: Neutral Candidate Selection - Research

**Researched:** 2026-06-12
**Domain:** Nextflow DSL2 channel-topology refactor + R selection-script rewrite (HCV strain candidate selection)
**Confidence:** HIGH (all findings verified against the actual repo files this phase rewrites)

## Summary

Phase 6 turns reference selection from a hardwired two-slot **major/minor** asymmetry into a **neutral, parameterized top-N candidate set**, mapped through one uniform path, while keeping the v1.0 on-disk contract intact behind a shim (D-06). The rewrite touches exactly six known sites: the R selection script (`bin/summarize_mapping_to_all_references.R`), the genotype helper (`bin/genotype_utils.R`, read-only reuse), the `PARSEFIRSTMAPPING` module, the `hcvtyper.nf` routing block (~lines 371–490), the `TARGETED_MAPPING` subworkflow (mapping discipline to honor, no rewrite needed), `bin/summarize.R` (consumer, must keep working unchanged), plus the typed-param plumbing (`nextflow_schema.json` + `conf/modules_hcv.config` + `nextflow.config`).

The single most important — and under-appreciated — finding is that the **`.major.`/`.minor.` filename slot is NOT carried in candidate data**. It is produced entirely by `ext.prefix` closures in `conf/modules_hcv.config` that are keyed on the **subworkflow alias name** in the `withName:` selector (`...:MAJOR_MAPPING:STATS_WITHDUP` → `.major.`, `...:MINOR_MAPPING:...` → `.minor.`). `bin/summarize.R` then parses the 3rd dot-field of those filenames (`str_split(...)[[1]][3]`) to decide major vs minor. D-04 collapses `MAJOR_MAPPING`/`MINOR_MAPPING` into one aliased subworkflow — which **destroys the alias-name signal those `withName` selectors depend on**. Preserving the legacy `.major.`/`.minor.` filenames (D-06) therefore requires deriving the slot from candidate metadata (`candidate_rank` → slot string) inside the `ext.prefix` closures, NOT from two distinct alias names. This is a genuine Decision Tension that the planner must resolve explicitly (see Decision Tensions §).

The second key finding: the R script writes `_major.fa`/`_minor.fa` and the 10-column `major_*/minor_*/minor_call/gate_flag` CSV. Under the shim, both must survive at default N=2 (cand_1→major slot, cand_2→minor slot), with the new long-format table + per-candidate `confirmation_status` emitted **alongside** them. The validity rules (`is_valid_minor()`) are deleted from selection per D-05 — meaning at N=2 the second candidate can differ from today's minor, which is exactly why ROADMAP success-criterion #5 cannot hold standalone and must be re-scoped (D-05 CONSEQUENCE).

**Primary recommendation:** Rewrite the R script to emit BOTH (a) a new long-format candidate table (one row per candidate) and (b) the reconstructed legacy 10-column wide CSV + `_major.fa`/`_minor.fa` files, deriving the legacy slots from `candidate_rank` (1→major, 2→minor). In `PARSEFIRSTMAPPING`, add a `candidates` (long-format) emit while keeping `major_mapping`/`minor_mapping`/`csv`. In `hcvtyper.nf`, replace the two filter-and-route blocks with a single `splitCsv`-driven per-candidate fan-out into ONE `TARGETED_MAPPING` call, carrying `candidate_rank` + `confirmation_status` in the meta map. Make the filename slot a function of `candidate_rank` in `modules_hcv.config` so `summarize.R` keeps parsing `.major.`/`.minor.` unchanged. Add `params.n_candidates` (integer ≥1, default 2) modelled on the `denovo_*` precedent. **Do not delete the legacy columns/filenames in Phase 6** — that is Phase 9.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Neutral candidate ranking (top-N distinct subtypes by reads) | R selection script (`summarize_mapping_to_all_references.R`) | — | Established "R-emits-decision, Nextflow-routes" pattern; ranking is data logic |
| Subtype→genotype derivation for dedup | R helper (`genotype_utils.R`, sourced) | — | Canonical 2k1b-aware rule already single-sourced and staged as `path()` input |
| Emitting candidate set as channel-routable rows | `PARSEFIRSTMAPPING` module (CSV emit) | Nextflow `splitCsv` in workflow | R writes rows; Nextflow fans them out |
| Per-candidate fan-out to mapping | `hcvtyper.nf` workflow routing | `TARGETED_MAPPING` subworkflow | Workflow owns channel topology; subworkflow owns the uniform per-ref path |
| Per-reference mapping + stats + consensus | `TARGETED_MAPPING` subworkflow | nf-core modules (bowtie2/samtools/ivar) | Already uniform; candidates reuse it as-is |
| Output filename slot (`.major.`/`.minor.`) | `conf/modules_hcv.config` `ext.prefix` closures | `meta.reference` / `meta.candidate_rank` | Slot is config-driven, NOT data-driven today (alias-name keyed) — see Decision Tensions |
| Selection-CSV + filename consumption | `bin/summarize.R` | — | Downstream consumer; MUST stay unchanged in Phase 6 (D-06) |
| Candidate-count parameter | nf-schema + `conf/modules_hcv.config` + `nextflow.config` | — | Typed-param triple, `denovo_*` precedent |

## Standard Stack

No new external packages. This phase is pure refactor of existing project code using the already-installed stack.

### Core (already present — no install)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Nextflow DSL2 | ≥23.04.0 [CITED: nextflow.config] | Channel topology / fan-out | Project workflow engine |
| nf-schema plugin | v2.1.0 [CITED: CLAUDE.md] | Typed param validation | Already validates `denovo_*` params |
| R tidyverse | (seqera container `r-seqinr_r-tidyverse:5358395134867368`) [CITED: parsefirstmapping/main.nf] | Ranking / dedup / CSV | Already the PARSEFIRSTMAPPING container |
| R seqinr | (same container) | FASTA read/write of selected refs | Already used by selection script |

### Supporting (test infra — already present)
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| nf-test | 0.9.3 (pinned `~/.nf-test`) [VERIFIED: MEMORY.md nf-test-local-install] | Module/workflow snapshot tests | PARSEFIRSTMAPPING `.nf.test`; invoke via PATH prefix + `--profile docker` |
| Rscript test harness | `bin/tests/run_all.sh` (glob `test_*.R`) [VERIFIED: bin/tests/run_all.sh] | R subprocess-contract tests | Add a `test_candidate_selection.R`; picked up automatically by the glob |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Long-format candidate CSV (D-01, locked) | Numbered wide slots `cand1_ref,cand2_ref…` | Rejected in D-01: ragged columns, doesn't scale to N, harder fan-out |
| Slot-from-`candidate_rank` in `ext.prefix` | Two retained alias subworkflows | Keeping two aliases contradicts D-04 (uniform path); see Decision Tensions |

**Installation:** None. (No `## Package Legitimacy Audit` required — this phase installs no external packages.)

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| REFSEL-01 | Selection ranks candidates neutrally (`cand_1..cand_n`) by read recruitment, no major/minor semantics | Rewrite `summarize_mapping_to_all_references.R` ranking: group by Subtype, sum reads, top ref per subtype, top-N subtypes by reads (replaces the asymmetric major-by-reads / minor-by-coverage-breadth logic at lines 76–162). Helper `genotype_from_subtype()` reused for dedup. |
| REFSEL-02 | Candidate count is a parameter, default 2 | Add `params.n_candidates` (integer ≥1, default 2) across schema + `modules_hcv.config` + `nextflow.config`, passed as a new positional arg to the R script (precedent: `denovo_*` params, `minRead`/`minCov` already passed positionally). |
| REFSEL-03 | Each candidate independently targeted-mapped; `gate_flag`/`minor_call` → per-candidate `confirmation_status` | Collapse `MAJOR_MAPPING`/`MINOR_MAPPING` into one `TARGETED_MAPPING` fan-out over `splitCsv` of the long-format table; carry `confirmation_status` per row through the meta map with explicit metadata-key joins. |

## Architecture Patterns

### System Architecture Diagram (Phase-6 target topology)

```
                 GET_MAPPING_STATS_WITHDUP
                  (idxstats + depth, dups in)
                           |
                           v
        ch_parsefirstmapping = idxstats.join(depth).filter{tsv non-empty}
                           |
                           v
   ┌──────────────────────────────────────────────────────────┐
   │ PARSEFIRSTMAPPING                                          │
   │  summarize_mapping_to_all_references.R                     │
   │   • rank: per-subtype top ref, top-N subtypes by reads     │
   │   • NO validity filtering (D-05)                           │
   │   • emit long-format candidate table (1 row / candidate)   │
   │   • ALSO reconstruct legacy 10-col wide CSV (shim, D-06)   │
   │   • write cand_k FASTA → _major.fa (k=1) / _minor.fa (k=2) │
   └──────────────────────────────────────────────────────────┘
        | candidates (long CSV)        | csv (legacy wide, for summarize.R)
        v                              |
   splitCsv(header,row-per-candidate)  |
        |  + KRAKEN2_FOCUSED reads (join by meta.id)            |
        |  meta += {candidate_rank, candidate_ref,             |
        |           candidate_reads, candidate_cov,            |
        |           confirmation_status}                       |
        v                                                      |
   ┌──────────────────────────────────────────────┐           |
   │ TARGETED_MAPPING  (ONE call, per-candidate)   │           |
   │  bowtie2 build/align → sormadup → stats →      │           |
   │  depth → ivar consensus → distance             │           |
   │  meta.reference + meta.candidate_rank drive     │           |
   │  ext.prefix slot (.major./.minor. via rank)     │           |
   └──────────────────────────────────────────────┘           |
        | stats/depth/consensus (filenames carry .major./.minor.)
        v                                                      v
                     SUMMARIZE  (bin/summarize.R)  <-----------+
        parses filename 3rd dot-field (major/minor) UNCHANGED
        reads legacy major_reads/minor_reads/gate_flag UNCHANGED
```

### Recommended Project Structure
No new files required structurally; changes land in existing files. Optional new R test:
```
bin/
├── summarize_mapping_to_all_references.R   # PRIMARY rewrite (ranking + long-format + shim)
├── genotype_utils.R                        # reuse unchanged (sourced)
└── tests/
    └── test_candidate_selection.R          # NEW: ranking + dedup + shim contract (auto-globbed)
modules/local/parsefirstmapping/main.nf     # add `candidates` emit; keep major/minor/csv
workflows/hcvtyper.nf                        # collapse routing to one fan-out (~371–490)
conf/modules_hcv.config                      # slot-from-rank ext.prefix; n_candidates param
nextflow_schema.json + nextflow.config       # declare params.n_candidates
```

### Pattern 1: R-emits-decision, Nextflow-routes (extend, don't replace)
**What:** The selection script writes routing columns; Nextflow filters/fans out on them via `splitCsv`.
**When to use:** All candidate routing in this phase.
**Example (current, the pattern to generalize):**
```groovy
// Source: workflows/hcvtyper.nf:383-403 (current major routing)
ch_major_mapping = PARSEFIRSTMAPPING.out.major_mapping.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    .map { meta, _csv, major_fasta, _reads ->
        def elements = _csv.splitCsv( header: true, sep:',')
        def new_meta = meta + elements[0]            // lifts ALL csv cols into meta
        assert new_meta.id == new_meta.sample : "Metadata mismatch: id=${new_meta.id}, sample=${new_meta.sample}"
        tuple(new_meta, major_fasta, _reads)
    }
    .filter { entry -> entry[0]['major_reads'].toInteger() > params.minRead && entry[0]['major_cov'].toInteger() > params.minCov }
```
**Generalized target (long-format fan-out):** `splitCsv` the candidate table → one channel element per candidate row → `meta += row` (carrying `candidate_rank`, `candidate_ref`, `confirmation_status`) → join classified reads by `meta.id` → ONE `TARGETED_MAPPING`. Because `splitCsv` over a multi-row file already yields one element per row, the fan-out is natural (D-01 rationale confirmed by code).

### Pattern 2: Metadata-key joins, never positional (mandatory)
**What:** `TARGETED_MAPPING` joins reads/index/fasta by `meta` key, not by emission order.
**Why it matters here:** With N candidates per sample mapped through one subworkflow, candidates for the SAME sample share `meta.id`. The meta key must be unique per candidate or joins will cross-pair. Add `candidate_rank` (and/or `reference`) to the meta map so each candidate is a distinct join key.
```groovy
// Source: subworkflows/local/targeted_mapping/main.nf:29-33,49-56
def new_meta = meta + [ reference: fasta.getBaseName().toString().split('\\.').last() ]
...
ch_aligned_input = ch_input.reads.join( BOWTIE2_BUILD.out.index ).join( ch_input.build )  // all by meta key
```
**Action:** Ensure the per-candidate meta is unique. `meta + [reference: ...]` already gives uniqueness when candidates differ in reference; add `candidate_rank` defensively so two same-reference candidates (shouldn't happen post-dedup, but safe) can't collide.

### Pattern 3: Slot string derived from candidate rank (the shim hinge)
**What:** `ext.prefix` must emit `.major.` for rank 1 and `.minor.` for rank 2 so `summarize.R` parses unchanged.
**Example (target):**
```groovy
// conf/modules_hcv.config — one withName for the single TARGETED_MAPPING:STATS_WITHDUP
withName: 'FOLKEHELSEINSTITUTTET_HCVTYPER:HCVTYPER:TARGETED_MAPPING:STATS_WITHDUP' {
    ext.prefix = { "${meta.id}.${meta.reference}.${meta.candidate_rank == 1 ? 'major' : 'minor'}.withdup" }
}
```
**Why:** Today the slot comes from two alias names (`MAJOR_MAPPING:` vs `MINOR_MAPPING:`). With one aliased subworkflow that distinction is gone (Decision Tension D-04↔D-06). Deriving the literal slot from `meta.candidate_rank` keeps the exact `.major.`/`.minor.` strings that `summarize.R:226,291` parse — preserving the Phase-9-deferred rename boundary.

### Anti-Patterns to Avoid
- **Routing on positional channel order across N candidates** — guaranteed mis-pairing under parallelism; always join by meta key (project's documented TARGETED_MAPPING fix).
- **Deleting legacy `major_*`/`minor_*`/`minor_call`/`gate_flag` columns or `_major.fa`/`_minor.fa` files in Phase 6** — breaks `summarize.R` immediately; these are Phase 9 (COMPAT-02/03).
- **Renaming the filename slot to `.cand1.`/`.cand2.` now** — that is the lockstepped Phase-9 migration; doing it here silently drops summary rows (STATE.md blocker).
- **Hardcoding an absolute path to `genotype_utils.R`** — must stay `source("genotype_utils.R")` with the `path(genotype_utils)` staging (container portability).
- **Re-deriving `.toInteger()` on possibly-NA candidate fields in Groovy** — current minor route deliberately routes on the R-emitted string `minor_call=='yes'` to avoid `NA.toInteger()` crashes (hcvtyper.nf:431). Carry the gate/confirmation decision as a string from R.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Subtype→genotype with 2k1b special case | New inline `substr`/`if` logic | `genotype_from_subtype()` (`bin/genotype_utils.R`) | Canonical single-sourced 2k1b-aware rule; already staged + tested |
| Per-row channel fan-out | Manual index loops / flatMap juggling | Nextflow `splitCsv(header:true)` | Already the project idiom (hcvtyper.nf:388) — one element per row free |
| Empty/NA-safe selection guards | New guard logic | Generalize the existing length-zero rowwise guard + `as.numeric` coercion + default `gate_flag="no_mapping"` | Already battle-tested (summarize_mapping_to_all_references.R:21-26,63-68,136-148) |
| FASTA write of selected refs | New writer | `seqinr::write.fasta` (already used) | In-container, already imported |
| Uniform mapping/stats/consensus | New mapping subworkflow | `TARGETED_MAPPING` as-is | Already the single uniform path candidates fan into |

**Key insight:** Almost every "new" capability this phase appears to need already exists as a reusable asset. The work is *re-wiring and generalizing*, not building. The risk is in the wiring boundaries (filename slot, meta-key uniqueness, shim column reconstruction), not in any algorithm.

## Runtime State Inventory

This is a code/config refactor, not a rename of stored state. Checked all five categories:

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | **None** — no datastore keys/collections embed "major"/"minor". The `.major.`/`.minor.` strings live only in run-output filenames + transient CSV columns, regenerated each run. | None |
| Live service config | **None** — no external UI/DB-resident config references these names. (HCVGLUE feed consumes mapping BAMs by collect(), not by slot name — see Decision Tensions for the Phase-5-paused coupling.) | Verify HCVGLUE wiring still receives all candidate BAMs after fan-out (STATE.md blocker) |
| OS-registered state | **None** — no scheduler/launchd/systemd entries. | None |
| Secrets/env vars | **None** — no secret/env name references major/minor or candidate count. | None |
| Build artifacts | **None new** — R scripts are interpreted; no compiled artifacts. nf-test snapshots WILL need regeneration (see Validation Architecture). | Regenerate `parsefirstmapping/tests/main.nf.test.snap` + `tests/default.nf.test.snap` after behavior change |

**The canonical question — after every file is updated, what runtime systems still cache the old string?** Answer: only nf-test snapshot files (`.snap`) hold the old expected output; they are regenerated, not migrated. No live/stored runtime state carries major/minor.

## Common Pitfalls

### Pitfall 1: Filename slot disappears when aliases collapse
**What goes wrong:** Removing `MAJOR_MAPPING`/`MINOR_MAPPING` aliases removes the `withName:` selectors that inject `.major.`/`.minor.` into `ext.prefix`. Mapping outputs lose their slot field; `summarize.R`'s `str_split(...)[[1]][3]` returns the wrong token (e.g. `nodup`/`withdup`), silently corrupting Major/Minor columns.
**Why it happens:** The slot is config/alias-name-driven, not data-driven (the non-obvious core finding).
**How to avoid:** Derive slot from `meta.candidate_rank` in the single `TARGETED_MAPPING:*` `ext.prefix` closures (Pattern 3). Verify a real `--profile docker` run produces `<id>.<ref>.major.withdup.stats` and `.minor.` for N=2.
**Warning signs:** Summary.csv `Major_reference`/`Minor_reference` blank or swapped; `first_major_minor` column containing `withdup`/`nodup`.

### Pitfall 2: Meta-key collision across same-sample candidates
**What goes wrong:** Two candidates of one sample share `meta.id`; positional or id-only joins in `TARGETED_MAPPING` cross-pair index/fasta/reads.
**Why it happens:** N>1 candidates per sample now flow through one subworkflow simultaneously.
**How to avoid:** Make per-candidate meta unique (`reference` + `candidate_rank`). The subworkflow already joins by full meta map; ensure the map differs per candidate before entering it.
**Warning signs:** Consensus/stats files for one candidate built against another candidate's reference; nondeterministic snapshot diffs.

### Pitfall 3: NA propagation in single-candidate / no_mapping cases
**What goes wrong:** A sample with one candidate (or zero) yields NA `candidate_reads`/`cov`; Groovy `.toInteger()` on NA crashes; the legacy `_minor.fa` write fires on an empty `minor_ref`.
**Why it happens:** The script generalizes from a guaranteed-1-major/optional-1-minor shape to a variable-length set.
**How to avoid:** Reuse the existing guards: only write `_minor.fa` when `length(cand_2_ref) > 0`; route on R-emitted strings not Groovy numeric coercion; keep default `gate_flag="no_mapping"` and a per-candidate default `confirmation_status`. Preserve the documented latent no_mapping FASTA-write crash behavior (out of scope to fix — REQUIREMENTS.md Out of Scope / STATE.md T-04-03).
**Warning signs:** `NumberFormatException` on NA; FASTA-write error when no second candidate exists.

### Pitfall 4: Treating ROADMAP success-criterion #5 as a Phase-6 gate
**What goes wrong:** Verifier blocks Phase 6 because the N=2 candidate set doesn't byte-reproduce today's two-slot selection on fixtures.
**Why it happens:** D-05 removes validity filtering from selection, so the 2nd candidate legitimately differs (e.g. a same-genotype second subtype that `is_valid_minor()` excludes today).
**How to avoid:** Per D-05 CONSEQUENCE (locked), criterion #5 is re-scoped to "reproduces across Phase 6+8, verified at Phase 9". Phase 6 verification asserts topology + no-crash, not golden reproduction (see Validation Architecture / Open Questions).
**Warning signs:** Plan tasks that diff Phase-6 output against the v1.0 golden baseline as a blocking acceptance test.

## Code Examples

### Current selection ranking (the asymmetry being replaced)
```r
# Source: bin/summarize_mapping_to_all_references.R:76-104,136-162
summary <- df %>% group_by(Subtype, Genotype) %>%
  summarise(reads = sum(X3)) %>% arrange(desc(reads))
major_subtype <- summary$Subtype[1]                      # MAJOR = most reads
major_ref <- df %>% filter(Subtype == major_subtype) %>%
  arrange(desc(X3)) %>% head(n = 1) %>% pull(X1)          # top ref in that subtype
# MINOR chosen by COVERAGE BREADTH among VALID different-genotype refs:
tmp <- candidates %>% rowwise() %>%
  filter(is_valid_minor(cur_data())) %>% ungroup() %>%
  arrange(desc(percent_gt_4)) %>% slice(1)                # <-- asymmetric metric (D-03 removes)
```
**Target (neutral, by reads, distinct-subtype dedup, top-N):** group by Subtype → top ref per subtype by reads → arrange subtypes by total reads → take `head(n = n_candidates)`. NO `is_valid_minor()`. Emit one row per selected candidate.

### Current legacy CSV schema (must be reconstructable under shim)
```r
# Source: bin/summarize_mapping_to_all_references.R:57-58
colnames(df_final) <- c("sample","total_mapped_reads","major_ref","major_reads",
  "major_cov","minor_ref","minor_reads","minor_cov","minor_call","gate_flag")
```
**Shim rule:** rank 1 → `major_*`; rank 2 → `minor_*`; `minor_call`/`gate_flag` reconstructed from the existing major-pass/minor-pass threshold logic (lines 170-177) so `summarize.R:177-193` still finds `major_reads`,`minor_reads`,`gate_flag`. Add long-format columns alongside (Claude's discretion on exact names; suggested: `sample,candidate_rank,candidate_ref,candidate_subtype,candidate_genotype,candidate_reads,candidate_cov,confirmation_status`).

### Current FASTA write (slot tag origin in R)
```r
# Source: bin/summarize_mapping_to_all_references.R:186-189
write.fasta(sequences=fasta[major_ref], names=major_ref, file.out=paste0(sampleName,".",major_ref,"_major.fa"))
if (length(minor_ref) > 0)
  write.fasta(sequences=fasta[minor_ref], names=minor_ref, file.out=paste0(sampleName,".",minor_ref,"_minor.fa"))
```
**Note:** The `_major.fa`/`_minor.fa` SUFFIX here (R-side) and the `.major.`/`.minor.` SLOT (config-side `ext.prefix`) are TWO separate mechanisms. Both must be preserved. Keep this R write keyed on rank (1→`_major.fa`, 2→`_minor.fa`).

### Typed-param precedent
```json
// Source: nextflow_schema.json:191-216 (denovo_* pattern to mirror)
"denovo_min_contig_length": { "type": "integer", "default": 500, "description": "..." }
```
```groovy
// Source: conf/modules_hcv.config:16-21 + nextflow.config:46-50
minRead = 499; denovo_match_level = 'genotype'  // declared in params{} block
```
**Recommended new param:** `n_candidates` — schema: `{"type":"integer","default":2,"minimum":1,"description":"Number of neutrally-ranked candidate references to select and map (cand_1..cand_n). Default 2 reproduces the legacy two-slot topology under the compatibility shim."}`. Declare in `nextflow.config` params{} and `conf/modules_hcv.config` params{}; pass as a positional arg to the R script and reference in `PARSEFIRSTMAPPING` script block as `${params.n_candidates}`.

## State of the Art

| Old Approach | Current (Phase-6 target) | When Changed | Impact |
|--------------|--------------------------|--------------|--------|
| Major-by-reads + minor-by-coverage-breadth (asymmetric) | Uniform read-recruitment ranking, top-N distinct subtypes | Phase 6 (D-03) | 2nd candidate may differ from today's minor |
| Validity filter (`is_valid_minor`) inside selection | No filtering in selection; validity → Phase 8 | Phase 6 (D-05) | Selection is mechanical; can select same-genotype/2k1b pairs |
| Two alias subworkflows `MAJOR_MAPPING`/`MINOR_MAPPING` | One `TARGETED_MAPPING` fan-out | Phase 6 (D-04) | Slot must move from alias-name to `candidate_rank` |
| `gate_flag`/`minor_call` 2-slot plumbing | Per-candidate `confirmation_status` (+ legacy retained) | Phase 6 (D-06 shim) | New field additive; legacy stays until Phase 9 |

**Deprecated/outdated (but RETAINED until Phase 9):** `.major.`/`.minor.` filename slot, `major_*`/`minor_*`/`minor_call`/`gate_flag` columns, `_major.fa`/`_minor.fa` files. Do NOT remove in Phase 6.

## Decision Tensions

These are points where the locked decisions (D-01..D-06) collide with the code as it actually exists. Flagging per instructions rather than silently working around them.

### T-1 (HIGH severity): D-04 (collapse aliases) vs D-06 (keep `.major.`/`.minor.` filenames)
**The conflict:** The `.major.`/`.minor.` filename slot is produced ONLY by `withName: '...:MAJOR_MAPPING:...'` / `'...:MINOR_MAPPING:...'` selectors in `conf/modules_hcv.config` (lines 125-133, 199-205, 242-310). These selectors exist *because* there are two distinct alias subworkflows. D-04 deletes those aliases. Once collapsed, there is no alias-name to key the slot on — yet D-06 requires the exact `.major.`/`.minor.` strings to survive (because `summarize.R:226,291` parse them).
**Resolution (recommended, no decision reversal needed):** Derive the slot literal from `meta.candidate_rank` inside a single set of `TARGETED_MAPPING:*` `ext.prefix` closures (Pattern 3). At N=2: rank 1 → `major`, rank 2 → `minor`. This honors both D-04 (one subworkflow) and D-06 (same filenames). **The planner must explicitly include the `ext.prefix` rewrite in `modules_hcv.config` as a first-class task** — it is easy to miss because it lives in config, not code, and its omission silently corrupts Summary.csv (Pitfall 1).
**Residual risk:** For N>2 there is no legacy slot name for cand_3+. Phase 6 default is N=2 so this is latent; recommend cand_k≥3 get a non-legacy slot (e.g. `.cand3.`) OR Phase 6 documents N>2 as Phase-9-dependent. Flag for planner.

### T-2 (MEDIUM severity): D-04 fan-out vs paused-Phase-5 HCVGLUE feed
**The conflict:** STATE.md blocker notes `MAJOR_MAPPING`/`MINOR_MAPPING` aliases are consumed downstream by HCVGLUE feed wiring (hcvtyper.nf:443 `MAJOR_MAPPING.out.aligned...mix(MINOR_MAPPING.out.aligned...)`). Collapsing to one subworkflow changes that `.out.aligned` source from two channels to one.
**Resolution:** After collapse, the single `TARGETED_MAPPING.out.aligned` already contains all candidate BAMs (the fan-out maps every candidate). Replace the `.mix()` of two alias outputs with the single subworkflow's `.out.aligned.collect()`. Verify HCVGLUE still receives one BAM per candidate. Phase 5 is paused on master-fork; coordinate on resume (STATE.md). Low immediate risk because Phase 5 is not on the v3.0 branch, but the SUMMARIZE input channels (lines 468-490) `.mix(MAJOR..., MINOR...)` MUST be updated in lockstep or stats/depth/variation/consensus_distance silently drop to half.

### T-3 (LOW severity): D-05 removes filtering vs current nf-test assertions
**The conflict:** `parsefirstmapping/tests/main.nf.test` asserts validity-rule outcomes (GATE-03 2k1b single-major no-minor; GATE-04 no-minor-FASTA; GATE-02 failing-major). With validity rules removed from selection, a 2k1b sample that today yields "no valid minor" may now select a second candidate.
**Resolution:** These nf-tests must be updated to assert the NEW neutral behavior (snapshots regenerated). The GATE-02 failing-major assertions (major stats reported, no crash) largely still hold under the shim. The 2k1b "no minor" assertion (GATE-03/04) will change — expected per D-05. This is test maintenance, not a blocker. See Validation Architecture.

## Validation Architecture

> nyquist_validation not explicitly false in config — section included.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | nf-test 0.9.3 (module/workflow snapshots) + Rscript subprocess-contract tests |
| Config file | `nf-test.config` (root); R harness `bin/tests/run_all.sh` |
| Quick run command | `PATH=~/.nf-test:$PATH nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile docker` |
| R quick run | `bash bin/tests/run_all.sh` |
| Full suite command | `PATH=~/.nf-test:$PATH nf-test test --profile docker` (note: host disk near-full — see MEMORY.md; clean `work/` + `docker volume prune` if ENOSPC) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| REFSEL-01 | Neutral top-N-by-reads ranking, distinct-subtype dedup | unit (R subprocess) | `Rscript bin/tests/test_candidate_selection.R` | ❌ Wave 0 |
| REFSEL-01 | Long-format candidate CSV emitted by module | module (nf-test) | `nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile docker` | ⚠️ exists, needs new assertions + snapshot regen |
| REFSEL-02 | `params.n_candidates` validated, default 2 | schema/manual | `nextflow run . -profile test,docker --validate_params` (smoke) | ⚠️ param not yet declared |
| REFSEL-03 | Single uniform fan-out; `confirmation_status` carried; no index misalignment | workflow (nf-test) | `nf-test test tests/default.nf.test --profile docker` | ⚠️ exists, snapshot regen |
| REFSEL-03 / shim | `.major.`/`.minor.` filenames + legacy columns survive at N=2 | workflow assert | grep mapping output filenames for `.major.`/`.minor.`; assert Summary.csv has `Major_reference`/`Minor_reference` | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `bash bin/tests/run_all.sh` (R logic, seconds) + the single `parsefirstmapping` nf-test.
- **Per wave merge:** full module + `tests/default.nf.test` workflow snapshot under `--profile docker`.
- **Phase gate:** full suite green; **explicitly NOT** golden-baseline reproduction (D-05 — that gate is Phase 9).

### What Phase 6 verification SHOULD assert (per D-05, instead of criterion #5)
1. At default `n_candidates=2`, the run completes without crash on the regression fixtures (topology intact).
2. Mapping outputs still carry `.major.`/`.minor.` filename slots and `summarize.R` still produces populated `Major_reference`/`Minor_reference` columns (shim integrity).
3. The legacy `major_*`/`minor_*`/`minor_call`/`gate_flag` columns are present in the parsefirstmapping CSV.
4. A new per-candidate `confirmation_status` field exists in the long-format table.
5. For an unambiguous single-subtype-dominant fixture (where validity rules never fired anyway), the selected cand_1 == today's major. (Cases where `is_valid_minor` previously suppressed a minor are EXPECTED to differ — do not assert reproduction.)

### Wave 0 Gaps
- [ ] `bin/tests/test_candidate_selection.R` — covers REFSEL-01 (ranking + distinct-subtype dedup + shim column reconstruction). Auto-globbed by `run_all.sh`.
- [ ] New nf-test assertions in `parsefirstmapping/tests/main.nf.test` for the `candidates` emit + long-format CSV; regenerate `.snap`.
- [ ] Regenerate `tests/default.nf.test.snap` after the topology change (expected snapshot churn).
- [ ] Optional fixture: a 2-distinct-subtype sample where the 2nd candidate differs from today's minor, to lock the NEW D-03/D-05 behavior.
- [ ] Update `modules/local/parsefirstmapping/tests/main.nf.test` GATE-03/04 (2k1b/no-minor) assertions to reflect no-validity-filtering selection.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| nf-test | Module/workflow snapshot verification | ✓ | 0.9.3 (`~/.nf-test`, PATH prefix) | Run R subprocess tests + manual run |
| Rscript + tidyverse/seqinr | R logic + R tests | ✓ (in seqera container; host Rscript for `run_all.sh`) | container `r-seqinr_r-tidyverse:5358395134867368` | — |
| Docker | `--profile docker` nf-test | ✓ (project default engine) | — | singularity profile |
| Disk headroom | nf-test/docker work dirs | ⚠ near-full (~93%+, MEMORY.md) | — | `sudo rm -rf work/` + `docker volume prune` before full runs |

**Missing dependencies with no fallback:** None.
**Constraint to flag for planner:** Host disk is near-full; nf-test runs that fill `work/` can break Bash with ENOSPC (MEMORY.md). Prefer the fast R subprocess tests during iteration; gate full nf-test runs and clean between them.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Slot literal can be derived from `meta.candidate_rank` inside `ext.prefix` (meta is in scope in the closure) | Pattern 3 / T-1 | If `candidate_rank` isn't reliably in meta at the STATS process, slot reverts wrong → Summary.csv corrupted. Mitigation: it's set in the same meta-enrichment `.map` that sets `reference` (TARGETED_MAPPING:29-33), which IS in scope for `ext.prefix` per existing `meta.reference` usage (lines 126,261) — high confidence but verify on a real run. |
| A2 | `splitCsv` over a multi-row candidate CSV yields one channel element per row, enabling natural fan-out | Pattern 1 / D-01 | If consumed as `splitCsv()[0]` (first row only, as current code does), fan-out won't materialize. Current code reads `elements[0]` because the file is single-row today; the new code must iterate all rows. Behavioral, not a tool limitation. |
| A3 | HCVGLUE/SUMMARIZE `.mix(MAJOR..., MINOR...)` channels can be replaced by the single subworkflow's collected output without losing files | T-2 | If any per-alias `ext.config` difference (e.g. different publishDir) was load-bearing, collapsing could change outputs. Verified the two alias configs are symmetric except the slot literal. |
| A4 | nf-test snapshot churn is acceptable/expected (not a regression) for the behavior change | Validation / T-3 | If CI treats any snapshot change as failure without regeneration, gate blocks. Mitigation: regenerate snapshots as a planned task. |

## Open Questions (RESOLVED)

1. **Slot naming for N>2 candidates under the shim**
   - What we know: Legacy slots are only `.major.`/`.minor.` (2 names). Default N=2 maps cleanly.
   - What's unclear: What filename slot cand_3+ should get when a user sets `n_candidates>2` in Phase 6.
   - **RESOLVED:** For Phase 6, derive cand_k≥3 as `.cand{k}.` (forward-compatible with Phase 9) and document that `summarize.R` only consumes slots 1-2 until Phase 9. Default-2 (the only tested topology this phase) is unaffected either way.

2. **Exact `confirmation_status` vocabulary in Phase 6**
   - What we know: D-06 says the real gating/role logic is Phase 8; Claude's discretion (CONTEXT) on the Phase-6 vocabulary, but the current major-pass threshold behavior must stay observable.
   - What's unclear: Whether to emit a single neutral value (`selected`) or a provisional threshold outcome (`pass`/`below_threshold`) per candidate.
   - **RESOLVED:** Emit a per-candidate threshold outcome (reusing the existing `major_pass`/`minor_pass` comparison generalized to each candidate's reads/cov vs minRead/minCov) so nothing observable regresses, while keeping the legacy `minor_call`/`gate_flag` columns reconstructed for the shim. Keeps Phase 8 free to redefine the vocabulary.

3. **Does `summarize.R` require the `minor_*` columns to be non-NA, or just present?**
   - What we know: `summarize.R:183` does `pull(minor_reads)` unconditionally (will error if column ABSENT, tolerates NA value).
   - What's unclear: Nothing critical — column must be PRESENT; NA is fine (current no-minor path already leaves it NA).
   - **RESOLVED:** Always emit all 10 legacy columns (NA-filled when no 2nd candidate), exactly as today. Confirmed safe by current behavior.

## Sources

### Primary (HIGH confidence — direct repo reads)
- `bin/summarize_mapping_to_all_references.R` (full) — selection logic, schema, guards, FASTA write
- `bin/genotype_utils.R` (full) — `genotype_from_subtype()` 2k1b rule
- `modules/local/parsefirstmapping/main.nf` (full) — emit contract, staging, stub
- `workflows/hcvtyper.nf` lines 364-510, includes 25-68 — routing, aliases, SUMMARIZE/HCVGLUE feed
- `subworkflows/local/targeted_mapping/main.nf` (full) — uniform path, meta-key join discipline
- `bin/summarize.R` lines 150-334 — selection-CSV consumption + filename 3rd-dot-field parsing
- `conf/modules_hcv.config` (full) — alias-keyed `.major.`/`.minor.` `ext.prefix` (the T-1 finding)
- `nextflow_schema.json` lines 160-216, `nextflow.config` lines 46-50 — typed-param precedent
- `modules/local/parsefirstmapping/tests/main.nf.test`, `bin/tests/run_all.sh`, `bin/tests/test_major_gate.R` — test infra
- `.planning/phases/06-neutral-candidate-selection/06-CONTEXT.md`, `REQUIREMENTS.md`, `STATE.md`, `ROADMAP.md` — decisions/scope

### Secondary (MEDIUM confidence)
- `.claude/projects/.../memory/MEMORY.md` — nf-test local install, disk-near-full, commit_docs disabled, worktrees disabled

### Tertiary (LOW confidence)
- None — no WebSearch needed; this is an internal-code refactor verified entirely against the repo.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new packages; all tools verified present in repo/container.
- Architecture / routing: HIGH — exact current wiring read line-by-line; target topology follows existing idioms.
- Filename-slot tension (T-1): HIGH that it exists; MEDIUM that `meta.candidate_rank`-in-`ext.prefix` is the cleanest fix (A1 — verify on a real docker run).
- Pitfalls: HIGH — derived from documented project history (TARGETED_MAPPING join fix, NA.toInteger crash, latent no_mapping FASTA crash).
- Test impact: HIGH — fixtures and harness located and read.

**Research date:** 2026-06-12
**Valid until:** 2026-07-12 (stable internal codebase; re-verify if hcvtyper.nf routing or modules_hcv.config changes before planning)
