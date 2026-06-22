# HCVTyper — Hand-off: de novo-informed strain selection + major-gating

*Created 2026-06-05. Source: HCV manuscript revision (Access Microbiology ACMI-D-26-00042). This document is the bridge between the manuscript/benchmarking workspace (`/home/jon.brate/HCV_nextflow_pipeline`) and the pipeline code repo (github.com/folkehelseinstituttet/hcvtyper). Open it inside the hcvtyper session — it is self-contained and assumes no prior context.*

---

## 1. One-paragraph summary

HCVTyper currently selects the major and minor strains purely from first-pass **read-mapping** counts + coverage thresholds, and runs **de novo assembly + BLAST in parallel as QC output that is never fed back into the call**. Benchmarking for the manuscript revealed two concrete failure modes that the *already-computed* de novo/BLAST evidence would fix. This work: (1) **gate minor-strain evaluation on the major strain passing its thresholds**, and (2) **use the de novo/BLAST result to confirm or refute a candidate minor strain** before it is reported as a co-infection. A stretch goal (3) is to let de novo/BLAST inform the *choice of reference* for targeted mapping. None of this is implemented yet.

---

## 2. Why — the evidence from benchmarking

Two real failure modes were found in the benchmark data:

**Failure mode A — cross-mapping produces a false minor that clears all thresholds.**
In the simulated genotype **1a single-infection** sample, first-pass mapping flagged a **genotype 4g** "minor": **53,279 reads, 67.7 % coverage breadth (≥10×), mean depth 1,212×** — far above the defaults (≥500 reads, ≥30 % breadth). On mapping evidence alone it would be reported as a **1a:4g co-infection**. It is spurious: the reads are 1a reads cross-mapping to a divergent 4g reference. **De novo assembly produced only 1a contigs (one 9,189 bp contig, 100 % identity to 1a) and NO 4g contig.** The de novo evidence already in the pipeline output unambiguously refutes the minor — it just isn't used.

**Failure mode B — a minor is called on top of a failed major.**
Real sample **ERR1810469** (expected single 3a): the major 3a **failed** coverage (8.8 % breadth, depth 3.5×, 248 reads — below thresholds), yet a 1a minor (1,030 reads, 78 % breadth) was still called. A minor should never be evaluated when the major has not passed.

**Key negative result (why thresholds alone can't fix this).**
Across the benchmark, the false/artefactual minor calls overlap genuine co-infection minors on *every* mapping axis — reads 1,030–53,279 vs 589–277,320; breadth 67.7–98.2 % vs 12.5–100 %. No `minRead`/`minCov` setting separates them; raising thresholds removes genuine low-abundance co-infections before it removes artefacts. **The discriminator must be orthogonal (de novo/BLAST), not a stricter cutoff.**

**Crucial corollary — most flagged "minors" in single-genotype samples are REAL.** De novo confirmed genuine second-genotype genomes in several real single-genotype samples (e.g. ERR1810447: full **9,207 bp 2b** genome; ERR1810475: multi-kb **1a** contigs). A full/multi-kb second-genotype genome cannot come from cross-mapping. So the goal is **not** to suppress all minors — it is to **confirm genuine ones and refute artefacts**, using de novo. (Whether the real ones are co-infection vs contamination is left open; that's a biology question, not a pipeline one.)

### Evidence table (mapping call vs de novo)

| Sample | Mapping minor (reads; cov ≥10×) | De novo minor-genotype contig | Correct outcome |
|---|---|---|---|
| sim1_1a_single (sim, 1a) | 4g — 53,279; 67.7 % | **None** (only 1a; 9,189 bp, 100 %) | **Refute → no co-infection** |
| ERR1810469 (3a) | 1a — 1,030; 78.4 % | 1a 4,503 bp + 3a 4,154 bp; **major 3a failed** | **Major-gate → no minor** |
| ERR1810447 (1b) | 2b — 4,199; 90.3 % | **Full 2b genome — 9,207 bp, 91 %** | **Confirm minor** |
| ERR1810475 (2b) | 1a — 11,646; 98.2 % | 1a — 4,734 + 4,429 bp | **Confirm minor** |
| ERR1810453 (1a) | 2b — 3,263; 85.0 % | 2b — partial, 2,949 bp | Confirm (fragmented) |
| sim1 (true 1a:1b) | 1b | Full 1b (9,339 bp) + full 1a (9,076 bp) | Confirm (closely related — works) |
| sim2 (true 2a:3a) | 2a | Full 2a + full 3a (~9.5 kb) | Confirm |
| ERR1810505 (true 1a:3a) | 3a | Full 1a + full 3a (~9.1 kb) | Confirm |
| ERR1810511 (true 1a:2a) | 2a | 1a 8,177 bp + 2a 7,517 bp | Confirm |

---

## 3. Current architecture (VERIFY against the code first)

This is my best understanding from the benchmark outputs; **confirm exact module/script names in the repo before editing.**

```
QC/trim → Kraken2 (HCV read selection)
   → FIRST MAPPING: Bowtie2, all HCV reads vs full reference set
   → parsefirstmapping (script): count reads/coverage per reference;
        • candidate MAJOR = reference with most reads, IF passes minRead & minCov
        • candidate MINOR = highest-coverage reference of a DIFFERENT genotype, IF passes minRead & minCov
        • exceptions: 1a/1b allowed as co-infection; 2k/1b recombinant suppresses gt1+gt2 co-infection
   → TARGETED_MAPPING (subworkflow): re-map reads to selected major (+ minor) reference(s)
   → consensus (iVar) ; RAVs (HCV-GLUE, from BAM) ; reporting
   → de novo (SPAdes) + BLAST vs same reference set  ← runs in parallel, OUTPUT ONLY, not fed back
```

Likely files to find (names approximate): the `parsefirstmapping` process/script (Python?), the subworkflow that emits the major/minor reference choice into `TARGETED_MAPPING`, and the blast-parsing process that produces `*.blastparse.csv` (columns: `sample, major_ref, major_contig_length, minor_ref, minor_contig_length`) and `*_blast_out.csv` (per-contig BLAST hits: `qseqid, sseqid, subtype, pident, length, …, kmer_cov`).

**The enabling fact:** `blastparse` already reduces de novo contigs to a `major_ref` / `minor_ref` genotype call per sample. The integration is largely about *consuming* that existing output at decision time.

---

## 4. Proposed changes (phased by risk)

### Change 1 — Major-gate (small, low-risk, do first)
In the selection logic: **only evaluate/report a minor strain if the major strain passes both `minRead` and `minCov`.** If the major fails, report first-mapping stats but issue no minor call (and no major genotype call, matching current behaviour for failed majors). Removes failure mode B (ERR1810469). Affects only samples where the major fails — expected to change **only ERR1810469** in the benchmark.

### Change 2 — De novo/BLAST confirmation of the candidate minor (the core change)
After first-mapping selection **and** de novo/BLAST, cross-check the candidate minor before reporting a co-infection:

- **Confirm** the minor if de novo produced a contig whose best BLAST hit is the **same genotype** as the candidate minor, and the contig is "substantial" (see thresholds below).
- **Refute / downgrade** the minor if de novo *succeeded for the major* (assembled a substantial major-genotype contig) but produced **no** substantial contig of the candidate minor's genotype. → report as single infection, or report the minor as "unconfirmed — possible cross-mapping" rather than a co-infection.
- **Fall back** (keep mapping-only call, flag "de novo unconfirmed") if de novo failed overall (no good major contig either) — e.g. genuine but low-yield samples. **Do not suppress a minor just because de novo assembled nothing**, or you will lose real low-abundance co-infections (e.g. IVT extreme ratios).

The asymmetry is the important part: **absence of a minor contig only counts as evidence of an artefact when de novo otherwise worked** (assembled the major). The 4g case has a perfect 1a major contig + no 4g contig → strong refute. A low-yield sample with nothing assembled → inconclusive, keep the call but flag.

### Change 3 — De novo-informed reference selection (stretch goal, highest risk)
Use the de novo contig's best BLAST hit to pick the *best* major/minor reference for `TARGETED_MAPPING` (instead of read-count alone). Could improve divergent/novel-subtype handling (e.g. the 3b→3a misassignment). **This changes many results and needs the most re-benchmarking — keep it as a separate later phase, not bundled with 1–2.**

---

## 5. Design details & edge cases

- **Genotype vs subtype matching.** Match the de novo contig to the candidate minor at **genotype** level, not subtype — subtype matching is too strict (cf. the 3b→3a case, where the nearest reference is a different subtype but correct genotype). Make the match level a parameter.
- **"Substantial contig" thresholds.** Need a minimum contig length (e.g. ≥ ~1,000 bp, or ≥ X % of genome) and possibly a minimum BLAST identity/length, to ignore the hundreds of ~300 bp, ~1× k-mer-coverage noise contigs. Tune against the evidence table: must confirm ERR1810453's 2,949 bp partial 2b and the multi-kb 1a in ERR1810475, while not being fooled by the small noise contigs. The genuine minor contigs had k-mer coverage clearly above the ~1× background; consider a k-mer-coverage floor too.
- **Preserve existing logic:** the 1a/1b co-infection exception and the 2k/1b recombinant suppression must still work. Test sim1 (1a:1b) — de novo recovers both full genomes there, so confirmation should pass.
- **Closely related pairs:** verify de novo doesn't merge the two genomes for 1a:1b; benchmark shows it separated them (NODE_1=1b, NODE_2=1a), but confirm the confirmation rule is lenient enough.
- **Make it configurable / non-breaking:** consider a flag (e.g. `--denovo_confirm_minor true/false`, default decided after benchmarking) so the behaviour can be toggled and the old behaviour reproduced. Important for the manuscript: you want to be able to show before/after.
- **Reporting:** add explicit status fields to the genotype report — e.g. minor `confirmed_by_denovo` / `unconfirmed` / `refuted` — so the analyst sees the basis. This matches the manuscript framing (de novo as confirmatory evidence).

---

## 6. Re-benchmarking plan & expected outcomes

Re-run the same datasets and diff against the current `combined_analysis.tsv`. Expected changes **only** in:

| Sample | Now | After changes 1+2 |
|---|---|---|
| sim1_1a_single | reported 1a:4g | **1a single** (4g refuted: no de novo 4g contig) |
| ERR1810469 | reported 3a:1a | **major fails → no call / major-only** (gate) |
| ERR1810447 | 1b:2b | 1b:2b, now **de novo-confirmed** (unchanged call, added status) |
| ERR1810475 | 2b:1a | 2b:1a, **confirmed** (unchanged call) |
| ERR1810453 | 1a:2b | 1a:2b, **confirmed** (unchanged call) |

**Everything else must be unchanged.** Specifically verify NO regression on the true co-infections: sim1, sim2, ERR1810505, ERR1810511, the IVT mixtures (incl. extreme ratios where the minor is genuine but may not assemble — these must NOT be suppressed; confirm the fall-back rule protects them), and the Qiu SRR1762352–55 1b:3a samples.

Datasets (read-only mounts used for the manuscript):
- Simulated: `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCVTyper_sim_data/`
- Real-world: `/mnt/N/.../HCVTyper_SRA_data/`
- Per-sample de novo BLAST tables to validate logic against: `…/blastparse/*.blastparse.csv` and `…/_blast_out.csv`
- Combined analysis (truth + current results): `/home/jon.brate/HCV_nextflow_pipeline/data/thomson2016/combined_analysis/combined_analysis.tsv`

---

## 7. How to run this in the hcvtyper session

1. **Copy this file into the hcvtyper repo** (e.g. `docs/denovo_informed_selection_handoff.md`) or paste its path. The hcvtyper repo is a different directory, so this session's memory will not carry over — this doc is the context.
2. Open the hcvtyper repo as the workspace and start with a discovery step. Suggested first prompt:
   > "Read `docs/denovo_informed_selection_handoff.md`. Then locate the actual code for: (a) the first-mapping parse / major-minor selection logic, (b) the subworkflow feeding `TARGETED_MAPPING`, and (c) the blast-parsing process that emits `*.blastparse.csv`. Map the data flow and confirm or correct the architecture in §3 before proposing edits."
3. **Implement Change 1 (major-gate) first**, on a branch — it's small and self-contained. Add a test.
4. **Then Change 2 (de novo confirmation)**, behind a configurable flag, with the asymmetric fall-back rule in §4. Add the reporting status fields (§5).
5. **Re-benchmark** per §6 and diff against `combined_analysis.tsv`. Confirm only the five samples above change.
6. **Defer Change 3** unless time allows; it needs its own validation pass.
7. When validated and merged (and a release is cut), come back to the manuscript and **upgrade the future-work text to a claimed feature** — see the revision package note below.

---

## 8. Tie-back to the manuscript

- The manuscript currently frames major-gating + de novo-informed selection as **future work** (not implemented) — see `Access_Microbiology/Revision_1/revision_package.md` § V-b(d), and the single-infection Results paragraph in `revised_results_section.md`.
- If you implement, re-benchmark, and release **before resubmission** (deadline 2026-08-11), you can upgrade these to claimed features and show ERR1810469's minor suppressed and the 4g call refuted automatically. Update `revision_package.md` § V-b accordingly and regenerate any affected tables/numbers.
- If you do **not** finish before resubmission, the manuscript stands as-is (de novo = analyst-facing QC; gating = future work), which is already accurate and defensible.

---

## 9. Open questions for the implementer

1. What exactly should a refuted minor be reported as — single infection, or co-infection-candidate flagged "unconfirmed"? (Manuscript leans toward: don't report a co-infection, but surface the unconfirmed candidate in QC.)
2. Contig-length / k-mer-coverage / identity thresholds for "substantial contig" — calibrate against the evidence table.
3. Genotype-level match confirmed as the right granularity (not subtype)?
4. Should Change 1 (major-gate) also suppress the *major* genotype call when the major fails (current behaviour) or is that already the case? Verify.
5. Default on/off for the de novo-confirmation flag.
