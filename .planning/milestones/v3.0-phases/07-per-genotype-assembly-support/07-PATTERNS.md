# Phase 7: Per-Genotype Assembly Support - Pattern Map

**Mapped:** 2026-06-13
**Files analyzed:** 4 modified (3 R scripts + 1 Nextflow module); 1 read-only helper reused on both join sides
**Analogs found:** 4 / 4 (all analogs are in-file or sibling — this is a tightly self-referential additive change)

This phase has **no new file type to introduce** — every pattern it needs already exists in the four touchpoint files. The "analog" for each new construct is an existing block in the same script (or its sibling). The planner/executor should pattern-match against those exact blocks rather than invent shape. Two roles only:

- **R selection/evidence logic** (`bin/blast_parse.R`, `bin/summarize.R`, `bin/genotype_utils.R`, `bin/denovo_confirm.R`) — emits CSV columns / computes the neutral support roll-up and the candidate join.
- **Nextflow routing** (`modules/local/blastparse/main.nf`) — declares any new CSV output + stub contract.

---

## File Classification

| Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---------------|------|-----------|----------------|---------------|
| `bin/blast_parse.R` (§6/§7 → new per-subtype support roll-up) | R utility (evidence emitter) | transform / batch (BLAST hits → per-subtype CSV) | §7 `summary_tbl` block in the SAME file (lines 253-305) + `scaf_top` (148-156) | exact (same file, same idiom) |
| `bin/summarize.R` (new genotype-level candidate join + NA-fill) | R utility (aggregator) | CRUD / transform (read CSV → `left_join` → `Summary.csv`) | `df_denovo` read+typed-empty block (441-466) and the `final` `left_join` chain (717-735) in the SAME file | exact (PLUMB-02 pattern is the literal spec) |
| `bin/genotype_utils.R::genotype_from_subtype()` | R utility (pure helper) | transform (subtype → genotype) | — (reused unchanged on BOTH join sides; likely NOT modified) | n/a — consume as-is |
| `modules/local/blastparse/main.nf` (declare new support CSV if emitted separately) | Nextflow module (routing) | request-response (process I/O contract) | the existing `*blastparse.csv`/`*blast_out.csv` `emit:` + stub lines in the SAME file (24, 60-61) | exact (sibling declared output) |

---

## Pattern Assignments

### `bin/blast_parse.R` — neutral per-subtype support roll-up (R utility, transform)

**Analog (same file):** §7 `summary_tbl` writer (lines 253-305) and `scaf_top` (lines 148-156).

**What this phase builds:** A per-subtype support table (D-02 subtype grain), one row per subtype, each row carrying the single best contig **by `sc_length`** (D-03) and **that contig's** `pident`, BLAST `length`, and `kmer_cov`, plus a `genotype` key. This is the dominance-neutral replacement for the §7 `major_ref`/`minor_ref` logic — but §7 itself stays (D-04 legacy shim).

**Start object — `scaf_top` (best hit per scaffold) already exists** (lines 148-156). The roll-up groups `scaf_top` by `subtype` and picks the longest contig per group. This is the object `classify_minor_denovo()` reconstructs inline (`denovo_confirm.R` lines 52-56) — reuse `scaf_top` directly instead.

```r
## ── 4. Top BLAST hit per scaffold (all lengths) ----------------------------
scaf_top <- if (nrow(scaf) > 0) {
scaf %>%
  arrange(evalue, desc(bitscore)) %>%      # best hit = lowest e‑value, highest bitscore
  group_by(qseqid) %>% slice(1) %>% ungroup() %>% # take first hit per scaffold
  arrange(desc(bitscore))
} else {
  tibble()
}
```

**Already-extracted columns available on every `scaf`/`scaf_top` row** (lines 82-87) — no new extraction needed, ASUP-01's four metrics map directly:
- `subtype` (from `separate(sseqid, into = c("subtype", NA))`)
- `sc_length` = full contig length (`str_extract(qseqid, "(?<=_length_)[0-9]+")`) — best-contig selection key (D-03), distinct from BLAST `length`
- `kmer_cov` (`str_extract(qseqid, "(?<=_cov_)[0-9.]+")`)
- `pident`, `length` (BLAST % identity and BLAST alignment length — note the specifics: `length` is reported IN ADDITION TO `sc_length`)

**Best-contig-by-length idiom to copy** (the existing §7 major-contig selection, lines 258-261 / 293-296 — copy the `slice_max(sc_length, n = 1) ... distinct() %>% pull()` shape; the `distinct()` step is load-bearing because one contig can have multiple BLAST hits to the same reference and would otherwise duplicate the length):

```r
major_contig <- scaf_top %>%
  slice_max(sc_length, n = 1) %>%               # longest contig for this reference
  select(qseqid, sc_length) %>% distinct() %>%  # Remove duplicates if several hits against the same reference
  pull(qseqid)
```

For the support roll-up this becomes a per-subtype `group_by(subtype) %>% slice_max(sc_length, n = 1, with_ties = FALSE)` carrying `pident`, `length`, `kmer_cov` from that same winning row (D-03: one coherent contig backs all four numbers — NOT independent per-metric maxima).

**Genotype key (D-02) — derive via the canonical helper, NOT inline `str_sub`.** The existing §7 uses `major_geno <- str_sub(major_name, 1, 1)` (line 257) which is NOT 2k1b-aware — do NOT copy that for the support table. Instead derive `genotype = genotype_from_subtype(subtype)`. CAVEAT: `blast_parse.R` does NOT currently source `genotype_utils.R` (it loads only `tidyverse` + `seqinr`, lines 12-15) and the BLASTPARSE module does NOT stage `genotype_utils.R` (see module input block). Either add the source + a `path(genotype_utils)` module input, OR carry subtype grain to `summarize.R` and derive the genotype key there (where the helper IS already sourced). Per D-02/Discretion the genotype collapse happens at join time in `summarize.R` anyway, so emitting **subtype grain + raw subtype** from `blast_parse.R` and deriving the genotype key in `summarize.R` is the lower-friction route — confirm with planner.

**Empty / no-hits guard to mirror** (the §7 `if (nrow(scaf_top) > 0) { ... } else { NA }` shape, lines 254/268-272, and the §6 `if (nrow(scaf_top) > 0)` guard at line 237). The new support writer MUST emit a typed header-only CSV when `nrow(scaf_top) == 0` so `summarize.R`'s `map_dfr` read does not abort on skip-assembly (mirror line 110 `write_csv(tibble(), ...)` and the §6 guard).

**Output write idiom** (line 305):
```r
write_csv(summary_tbl, paste0(prefix, ".blastparse.csv"))
```
The support CSV uses the same `write_csv(<tbl>, paste0(prefix, ".<suffix>.csv"))` form. Discretion (CONTEXT D-04 / line 47): fold into the existing `*_blast_out.csv` (already read into `df_blast_out` and already carries `subtype`/`sc_length`/`kmer_cov`/`pident`/`length` per contig) OR emit a new `*.assembly_support.csv`. The `*_blast_out.csv` is the least-disruptive read path since `summarize.R` already ingests it; a dedicated pre-rolled `*.assembly_support.csv` is cleaner for the join. Planner to pick.

---

### `bin/summarize.R` — genotype-level candidate join + NA-fill (R utility, CRUD/transform)

**Analog (same file):** the `df_denovo` read + typed-empty-tibble block (lines 441-466, the PLUMB-02 pattern) and the `final` `left_join` chain (lines 717-735).

**This is the literal spec.** The new join MUST mirror the PLUMB-02 typed-empty-tibble guard exactly: read the support CSV(s) with `map_dfr(..., read_csv)` when files exist, else build a zero-row tibble that declares **every** support column with its correct type so the `left_join` always emits them and existing sample rows NA-fill (criterion #3: explicit `assembly_support = "none"` status + NA metrics, no row loss).

**Read + typed-empty pattern to copy verbatim in shape** (lines 441-466):
```r
blastparse_files <- list.files(path = path_denovo, pattern = "blastparse.csv$", full.names = TRUE)

if (length(blastparse_files) > 0) {
  df_denovo <- map_dfr(blastparse_files, read_csv) %>%
    rename(
      sampleName                 = sample,
      denovo_major_ref           = major_ref,
      ...
    )
} else {
  # PLUMB-02: ... Declare all columns with their blast_parse.R types so the join
  # always emits them; with zero rows here every existing sample row NA-fills.
  df_denovo <- tibble(
    sampleName                 = character(),
    denovo_major_ref           = character(),
    denovo_major_contig_length = integer(),
    ...
  )
}
```
The `df_blast_out` block right below (lines 476-492) is the SAME pattern keyed differently (basename-derived `sampleName` via `str_remove(basename(.x), "_blast_out.csv$")`) — use whichever read key matches the chosen output shape. The T-03-01 DoS guard ("never abort on skip-assembly; a sample absent → `unconfirmed`/`none`, never refuted") generalizes to the new support columns.

**`left_join` keyed by sampleName — the existing chain** (lines 717-735):
```r
final <- input_samplesheet %>%
  left_join(id_df, join_by(sampleName)) %>%
  ...
  left_join(df_distance_wide, join_by(sampleName)) %>%
  # Samplesheet anchors the left side so a sample with no de novo output keeps its
  # row with NA de novo fields (PLUMB-02).
  left_join(df_denovo, join_by(sampleName))
```

**NEW join shape (this phase's core add).** The candidate set (Phase 6 long-format, `cand_1..cand_n`) is the LEFT side; the per-genotype support roll-up is the RIGHT side. The join key is `(sampleName, genotype)` at the parameterized `denovo_match_level`. Both genotype keys are computed by the SAME helper (Discretion):
- Candidate side: `genotype_from_subtype(candidate_ref)` — note the Phase 6 candidates CSV ALREADY emits `candidate_genotype` and `candidate_subtype` precomputed (see `modules/local/parsefirstmapping/main.nf` line 74 header and `bin/summarize_mapping_to_all_references.R` lines 80-88), so the candidate-side key may already exist — prefer reusing `candidate_genotype` over recomputing.
- Support side: `genotype_from_subtype(subtype)` collapsing subtype-grain rows up to genotype at join time when `denovo_match_level == "genotype"` (D-02). When `denovo_match_level == "subtype"`, join on subtype directly (the finer grain is retained precisely so this stays honorable).

**Match-level switch idiom to copy** — `classify_minor_denovo()` already encodes the exact `match_level` branch the join must reuse (`denovo_confirm.R` lines 62-66):
```r
match_key = if (match_level == "subtype") subtype else genotype_from_subtype(subtype)
```
Apply this on BOTH sides before the join so keys are computed identically.

**Collapse-to-genotype before join (when match_level == genotype, D-03 carried through):** group support rows by `(sampleName, genotype)` and `slice_max(sc_length, n = 1, with_ties = FALSE)` so one best contig per genotype backs the four metrics — same single-best-contig discipline as `blast_parse.R`. This satisfies criterion #4: at default flags + N=2 the join's `cand_2` support must equal today's `denovo_minor_contig_length` (length) on the regression fixtures.

**Where to insert / what NOT to touch (D-04 legacy shim):** ADD the new join into the `final` chain (after the existing `df_denovo` join). Leave `apply_denovo_layer()` (lines 811-815), `minor_denovo_status`, `coinfection_flag`, and the `review_flag` `pmap_chr` build (lines 937-965) UNCHANGED — the legacy minor-coupled path runs in parallel this phase and its `denovo_minor_*` columns appear alongside the new per-candidate support columns in `Summary.csv` (intentional duplication, removed Phase 8/9). The `final %>% select(...)` reorder block (line 968+) must be extended to carry the new support columns through (mirror how the four `denovo_*_subtype*` columns are listed at lines 983-986).

**Reading the candidate CSV (new read needed):** `summarize.R` does NOT yet read `*.candidates.csv` — it currently consumes Phase 6 selection via `parsefirstmapping_df` (lines 153-207, the `*.parsefirstmapping.csv` read). The new long-format candidate read should follow the SAME `list.files(path = path_3, pattern = ...)` + `map_dfr(read_csv)` idiom used for `first_mapping_files` (line 153) / `blastparse_files` (line 441), with a typed-empty fallback. Planner to confirm whether the candidate long-table is published into an existing staged dir.

---

### `bin/genotype_utils.R::genotype_from_subtype()` — pure helper (reused, likely unmodified)

**No change expected** — consume as-is on both join sides. It is a sourced helper (no `commandArgs`), 2k1b-aware (`if_else(subtype == "2k1b", subtype, substr(subtype, 1, 1))`), already `source()`d at `summarize.R` line 10 and staged into `summarize.R`'s workdir. Critical contrast: `blast_parse.R` §7 uses the NON-2k1b-aware `str_sub(major_name, 1, 1)` (line 257) — the new support genotype key must use `genotype_from_subtype()` instead. If the genotype key is derived in `blast_parse.R`, that script must additionally `source("genotype_utils.R")` and the BLASTPARSE module must stage it as a `path` input (see next section); if derived in `summarize.R`, no change here.

```r
genotype_from_subtype <- function(subtype) {
  if_else(subtype == "2k1b", subtype, substr(subtype, 1, 1))
}
```

---

### `modules/local/blastparse/main.nf` — declare new support output (Nextflow routing)

**Analog (same file):** the existing `*blast_out.csv` / `*blastparse.csv` `emit:` declarations (lines 21, 24) and their stub lines (60-61).

**ONLY needed if the support table is a NEW CSV** (not folded into `*_blast_out.csv`). If folded into the existing `*_blast_out.csv`, NO module change is required — that channel already exists (`emit: blast_res`, line 21) and `summarize.R` already reads it.

**Declared-output pattern to copy** (lines 19-26):
```groovy
output:
tuple val(meta), path("*contigs.fa")    , emit: contigs    , optional: true
tuple val(meta), path('*blast_out.csv') , emit: blast_res
tuple val(meta), path("*major.fa")      , emit: major_fasta, optional: true
tuple val(meta), path("*minor.fa")      , emit: minor_fasta, optional: true
tuple val(meta), path("*blastparse.csv"), emit: csv
tuple val(meta), path("*.png")          , emit: png
path "versions.yml"                     , emit: versions
```
A new support CSV adds e.g. `tuple val(meta), path("*assembly_support.csv"), emit: support`. Match the existing alignment/comma style (Prettier `printWidth: 120`).

**Stub contract MUST be updated in lockstep** (the established discipline — every declared non-optional output needs a deterministic stub artifact, lines 60-61). Copy the stub `printf "<header>\n" > ${prefix}.<suffix>.csv` form:
```groovy
printf "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore\n" > ${prefix}.blast_out.csv
printf "id,header\n" > ${prefix}.blastparse.csv
```
The new support stub line must carry the real header columns the planner settles on (so a `-stub-run` of the workflow yields a parseable file). If staging `genotype_utils.R`, also add it to the `input:` block mirroring the PARSEFIRSTMAPPING `path(genotype_utils)` precedent (referenced in 06-CONTEXT line 83).

---

## Shared Patterns

### Genotype key derivation (2k1b-aware) — apply on BOTH join sides
**Source:** `bin/genotype_utils.R` (`genotype_from_subtype`, lines 24-26)
**Apply to:** support-table genotype key AND candidate-side genotype key — computed identically so join keys match (CONTEXT Discretion, D-02). Never use `str_sub(name, 1, 1)`.

### Typed-empty-tibble + `left_join` NA-fill (PLUMB-02 / T-03-01 DoS guard)
**Source:** `bin/summarize.R` `df_denovo` block (lines 441-466) and `df_blast_out` block (476-492)
**Apply to:** the new support-table read and the candidate-join read. Zero rows must declare every column with its type; the samplesheet-anchored `left_join` then NA-fills — no row loss, never abort on skip-assembly.

### Single-best-contig-by-`sc_length` selection (D-03)
**Source:** `bin/blast_parse.R` §7 major-contig idiom (lines 258-261, `slice_max(sc_length, n = 1) %>% select(...) %>% distinct() %>% pull()`)
**Apply to:** picking the one contig that backs all four metrics, both in the `blast_parse.R` roll-up and in the genotype-collapse step in `summarize.R`. The `distinct()` de-dup of multi-hit contigs is load-bearing.

### Match-level branch
**Source:** `bin/denovo_confirm.R` lines 62-66 (`if (match_level == "subtype") subtype else genotype_from_subtype(subtype)`)
**Apply to:** the join-key computation on both sides; `denovo_match_level` param already plumbed into `summarize.R` (arg 7, line 28; passed via `conf/modules_hcv.config` line 328).

### Additive + legacy-column retention (Phase 6 D-06 → this phase D-04)
**Source:** the existing `denovo_minor_*` columns and `apply_denovo_layer()`/`review_flag` chain in `summarize.R` (lines 811-965) stay UNCHANGED
**Apply to:** all four touchpoints — new support columns are ADDED alongside legacy columns; nothing legacy is rewired this phase.

---

## No Analog Found

None. Every construct this phase introduces has an exact in-repo precedent (listed above). The only genuinely new shape is the **genotype-level join key on two long-format tables** — but its two halves (typed-empty `left_join` and the `match_level` key branch) each have a literal analog, so the planner composes rather than invents.

---

## Metadata

**Analog search scope:** `bin/` (R helpers), `modules/local/blastparse/`, `modules/local/parsefirstmapping/`, `conf/modules_hcv.config`
**Files scanned:** `bin/blast_parse.R`, `bin/summarize.R`, `bin/genotype_utils.R`, `bin/denovo_confirm.R`, `bin/summarize_mapping_to_all_references.R`, `modules/local/blastparse/main.nf`, `modules/local/parsefirstmapping/main.nf`, `conf/modules_hcv.config`
**Pattern extraction date:** 2026-06-13
