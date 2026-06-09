# HCVTyper — De novo-Informed Strain Selection

## What This Is

HCVTyper is a Nextflow DSL2 pipeline that genotypes Hepatitis C virus (HCV) from NGS data and detects major/minor strain co-infections. This milestone closes a known gap: the pipeline already runs de novo assembly + BLAST in parallel, but the result is **output-only QC that is never fed back into the strain call**. This work makes the de novo/BLAST evidence (and a major-gating rule) actually inform which minor strains get reported — fixing two benchmarked failure modes that coverage thresholds alone cannot fix.

For: the HCVTyper maintainers/analysts at FHI and downstream users running clinical/research HCV genotyping. Tied to manuscript revision Access Microbiology ACMI-D-26-00042.

## Core Value

A reported minor strain (co-infection) must be backed by orthogonal de novo/BLAST evidence — never reported on top of a failed major, and never reported when de novo refutes it as a cross-mapping artefact — while genuine low-abundance co-infections are preserved.

## Requirements

### Validated

<!-- Inferred from existing codebase (brownfield) -->

- ✓ QC/trim (fastp|cutadapt) + Kraken2 HCV read selection — existing
- ✓ First-pass mapping (bowtie2|tanoti) of HCV reads vs full reference set — existing
- ✓ Major/minor strain selection from read-count + coverage thresholds (`parsefirstmapping`) — existing
- ✓ De novo assembly (SPAdes) + BLAST vs reference set, parsed to per-sample major/minor genotype (`blastparse` / `blast_parse.R`) — existing, **currently output-only**
- ✓ Targeted re-mapping + iVar consensus + HCV-GLUE genotyping/RAVs — existing
- ✓ Cross-sample contamination detection (all-vs-all BLAST) — existing
- ✓ Result aggregation (SUMMARIZE) + MultiQC report — existing
- ✓ 1a/1b co-infection exception and 2k/1b recombinant suppression in selection logic — existing

### Active

<!-- This milestone. Hypotheses until shipped and validated. -->

- [ ] **Change 1 — Major-gate:** Only evaluate/report a minor strain if the major passes both `minRead` and `minCov`. If the major fails, report first-mapping stats but issue no minor call (and no major genotype call, matching current failed-major behaviour).
- [ ] **Change 2 — De novo confirmation of the candidate minor:** After first-mapping selection and de novo/BLAST, cross-check the candidate minor before reporting a co-infection:
  - **Confirm** if de novo produced a substantial contig whose best BLAST hit is the **same genotype** (not subtype) as the candidate minor.
  - **Refute** (→ report as single infection, surface candidate in QC status only) if de novo assembled a substantial *major*-genotype contig but **no** substantial contig of the candidate minor's genotype.
  - **Fall back** (keep mapping-only call, flag "de novo unconfirmed") if de novo failed overall (no good major contig either) — must NOT suppress genuine low-yield minors.
- [ ] **Configurable + non-breaking:** Gate Change 2 behind a flag (e.g. `--denovo_confirm_minor`), **default ON**, so old behaviour is reproducible by toggling it off.
- [ ] **Reporting status fields:** Add explicit minor-strain status to the genotype report — `confirmed_by_denovo` / `unconfirmed` / `refuted` — so the analyst sees the basis for the call.
- [ ] **Preserve existing exceptions:** 1a/1b co-infection exception and 2k/1b recombinant suppression must continue to work after both changes.
- [ ] **Tests:** Unit/integration tests covering the major-gate, confirm/refute/fall-back branches, and the genotype-match (not subtype) logic.

### Out of Scope

- **Change 3 — De novo-informed reference selection** — deferred to a later milestone. Uses the de novo contig's best BLAST hit to pick the targeted-mapping reference. Highest risk, changes many results, needs its own re-benchmarking pass.
- **Full manuscript re-benchmarking as part of DoD** — re-running sim + real datasets and diffing against `combined_analysis.tsv` is a separate manual analyst step. This project ships code + tests; the analyst validates afterward.
- **Deciding co-infection vs contamination biology** — when a minor is confirmed genuine, whether it represents co-infection or contamination is a biology question, not a pipeline one.
- **Changing `minRead`/`minCov` thresholds as the fix** — benchmarking proved thresholds cannot separate artefact minors from genuine ones; the discriminator must be orthogonal (de novo), not a stricter cutoff.

## Context

- **Brownfield.** Codebase already maps cleanly (`.planning/codebase/`). Core selection logic lives in `modules/local/parsefirstmapping/` (with R helpers `summarize_mapping_to_all_references.R`); de novo evidence is parsed in `modules/local/blastparse/` (`bin/blast_parse.R`) producing per-sample `major_ref`/`minor_ref` genotype calls. Final reporting is `modules/local/summarize/` (`bin/summarize.R`).
- **Architecture note / correction:** the handoff guessed the parse scripts were Python; they are **R** (`blast_parse.R`, `summarize_mapping_to_all_references.R`, `summarize.R`). Verify exact selection/reporting touchpoints during planning.
- **The enabling fact:** `blastparse` already reduces de novo contigs to a major/minor genotype call per sample. Change 2 is largely about *consuming* existing output at decision time, not new computation.
- **Benchmark evidence (why this matters):**
  - Failure A: sim 1a single-infection → first-pass flagged a false 4g minor (53,279 reads, 67.7% breadth, 1,212× depth) clearing all thresholds; de novo produced only 1a contigs, no 4g → refute.
  - Failure B: ERR1810469 (3a) → major 3a failed coverage (8.8% breadth, 3.5× depth) yet a 1a minor was still called → major-gate fixes it.
  - Corollary: most flagged minors in single-genotype samples are REAL (de novo recovered full second-genotype genomes, e.g. ERR1810447 full 9,207 bp 2b). Goal is confirm genuine + refute artefacts, NOT suppress all minors.
- **Datasets (read-only mounts, for the analyst's later re-benchmark):** sim `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCVTyper_sim_data/`; real `/mnt/N/.../HCVTyper_SRA_data/`; truth+results `/home/jon.brate/HCV_nextflow_pipeline/data/thomson2016/combined_analysis/combined_analysis.tsv`.

## Constraints

- **Tech stack**: Nextflow DSL2 + nf-core conventions; HCV-specific logic in local modules with R/Python helper scripts in `bin/`. New logic should follow existing module/script patterns.
- **Compatibility**: Non-breaking — existing behaviour must be reproducible via the `--denovo_confirm_minor` flag (default ON). Existing 1a/1b and 2k/1b exceptions must be preserved.
- **Discriminator design**: Minor confirmation must match at **genotype** level, not subtype (subtype matching is too strict — cf. 3b→3a). Match level should be a parameter.
- **"Substantial contig" calibration**: A minimum contig length (≈≥1,000 bp or ≥X% of genome) plus possibly minimum BLAST identity/length and a k-mer-coverage floor, tuned to confirm partial genuine contigs (e.g. ERR1810453's 2,949 bp 2b) while ignoring ~300 bp / ~1× noise contigs.
- **Timeline**: Manuscript resubmission deadline 2026-08-11 — implementing + releasing before then lets future-work claims be upgraded to features, but the manuscript stands as-is if not.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Scope = Changes 1 + 2 only; defer Change 3 | Phased by risk; 1+2 are lower-risk and fix the two benchmarked failures. Change 3 needs its own validation. | — Pending |
| Refuted minor → report as single infection + QC flag | Don't report an unconfirmed co-infection; surface the candidate in QC status fields for analyst visibility. | — Pending |
| De novo-confirmation flag defaults ON | Make the improved behaviour standard; old behaviour reproducible by toggling off. | — Pending |
| Re-benchmarking NOT in DoD | Ship code + tests; full dataset re-benchmark is a separate manual analyst step. | — Pending |
| Match de novo↔minor at genotype level (parameterized) | Subtype matching too strict; genotype level handles divergent/novel subtypes (3b→3a). | — Pending |
| Asymmetric refute rule | Absence of a minor contig only refutes when de novo otherwise assembled a substantial major contig; protects genuine low-yield co-infections. | — Pending |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-06-05 after initialization*
