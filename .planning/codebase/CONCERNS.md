# Concerns & Technical Debt

**Analysis Date:** 2026-06-05

This document captures technical debt, known issues, fragile areas, and reproducibility concerns in `hcvtyper`. Severity reflects impact on correctness/maintainability, not effort.

## Outstanding Design Work (High Priority)

The repo contains an untracked handoff specification — `hcvtyper_handoff_denovo_informed_selection.md` (2026-06-05) — describing unimplemented improvements driven by manuscript benchmarking (Access Microbiology ACMI-D-26-00042). This is the primary planned work, not yet started:

1. **Major-gating not implemented** — a minor strain can be called even when the major strain fails its `minRead`/`minCov` thresholds. Real failure mode: sample ERR1810469 (expected single 3a) — major 3a failed (8.8% breadth, 3.5× depth) yet a 1a minor was still called. Selection logic should gate minor evaluation on the major passing. *Fix is small and self-contained (handoff Change 1).*

2. **De novo / BLAST evidence computed but never consumed** — SPAdes assembly + BLAST run in parallel as analyst-facing QC only (`workflows/hcvtyper.nf` runs `SPADES` → `BLAST_BLASTN` → `BLASTPARSE`, emitting `*.blastparse.csv` / `*_blast_out.csv`), but the result is never fed back into the major/minor call. Benchmarking showed this orthogonal evidence would refute cross-mapping artefacts (e.g. a spurious 1a:4g co-infection where de novo produced only 1a contigs). Integrating it is handoff Change 2 (the core change), behind a configurable flag with an asymmetric fall-back rule (absence of a minor contig only refutes when de novo otherwise assembled the major).

3. **De novo-informed reference selection (stretch)** — using the de novo BLAST hit to choose the targeted-mapping reference instead of read-count alone (handoff Change 3). Highest risk; deferred to its own phase.

## Known Bugs / Correctness Risks

- **`bin/summarize_mapping_to_all_references.R:144`** — `if (length(minor_ref > 0))` is almost certainly a logic error. `minor_ref > 0` is evaluated first (elementwise/coercion), then `length(...)` of the result is always ≥ 1, so the branch effectively always runs. Intent was likely `if (length(minor_ref) > 0)`. Could write a spurious minor FASTA when there is no minor reference.
- **Cross-mapping false minors** — documented in handoff: 1a reads cross-mapping to a divergent 4g reference clear all mapping thresholds. No `minRead`/`minCov` cutoff separates artefacts from genuine low-abundance co-infections (they overlap on every mapping axis), so a stricter threshold is not a valid fix.

## Fragile Areas

- **Hardcoded relative paths in R reporting scripts** — `bin/summarize.R` hardcodes input directories (`path_1 <- "trimmed/"`, `path_2 <- "kraken_classified/"`, `path_10 <- "variation/"`, `path_11 <- "consensus_distance/"`, etc.). These rely on Nextflow staging the right files into the work dir under those exact names; renaming a publish/stage path silently breaks aggregation with no error (uses `list.files` + `warning()`, not `stop()`).
- **`errorStrategy "ignore"` on optional tools** (e.g. `modules/local/tanoti.nf`) — failures are swallowed; downstream channels must tolerate missing outputs. Easy to mask real failures.
- **Channel ordering / join assumptions** — workflow relies on `join` operations keyed on meta maps (e.g. `BLASTPARSE.out.major_fasta.join(...)`); mismatched keys silently drop samples rather than erroring.

## Test Coverage Gaps

- **No unit tests for R or Python helper scripts** — only Black/isort formatting is checked in CI. R scripts (the bulk of `bin/`, including the buggy `summarize_mapping_to_all_references.R`) have zero automated tests → high risk of silent failures.
- **Single end-to-end nf-test** (`tests/default.nf.test`) with snapshot matching is the only behavioural test. The planned de novo-confirmation logic will need targeted tests against the handoff evidence table.
- **Docker-only CI** — no Singularity or Conda profile testing despite both being supported config paths.

## Performance Notes

- **SPAdes de novo assembly** is the heaviest per-sample step; high-coverage samples can be slow. Since its output is currently unused for calling, there is wasted compute until handoff Change 2 lands (after which it becomes load-bearing).
- **Kraken2 database** is large; first-run download/setup is a one-time cost worth documenting for reproducibility.

## Reproducibility Concerns

- **Handoff spec is untracked** (`hcvtyper_handoff_denovo_informed_selection.md` shows as `??` in git status). It should be moved into `docs/` and version-controlled so the planned work and its rationale are preserved.
- **Benchmark validation data lives outside the repo** — truth/results in `/home/jon.brate/HCV_nextflow_pipeline/.../combined_analysis.tsv` and mounted dataset paths (`/mnt/N/...`). Re-benchmarking per the handoff requires these external mounts; not reproducible from the repo alone.
- **Reference database versioning** — reference sets used for mapping/BLAST lack an in-repo changelog tying results to a specific reference version.

## Branch / Working State

- Active branch `dev` (PRs target `master`). Working tree shows recent consensus-distance work (major/minor) committed, plus the untracked handoff doc and a modified `.gitignore`.

---

*Concerns analysis: 2026-06-05*
