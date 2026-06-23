# HCVTyper pipeline — corrections & checks from simulated-data benchmarking

**Date:** 2026-06-23
**Reporter:** Jon Brate (FHI)
**Scope:** Code/behaviour issues in the **HCVTyper pipeline** (repo `folkehelseinstituttet/hcvtyper`). Manuscript-side changes and HCV-GLUE behaviour are tracked separately and are **out of scope here**.

## Context

`pipeline_version = folkehelseinstituttet/hcvtyper dev`, Summary.csv dated 2026-06-23 08:10.
Input: 7 simulated samples (InSilicoSeq, MiSeq model), known truth: 5 single-infections (1a, 1b, 2a, 3a, 4a) + 2 co-infections (sim1 = 1a:1b 7:3, sim2 = 2a:3a 7:3), each with a known NS5A RAV inserted at 20 % frequency.

What works correctly and needs no change: genotype/subtype calling (7/7 correct), de novo confirmation of the major (7/7), de novo refutation of spurious minors (single-infections all resolve to monoinfection), and RAV detection wherever the strain is typed (all 5 inserted RAVs recovered). The items below are about **reporting/column logic**, not the core calls.

---

## 1. Co-infection: role classifier and per-strain statistic columns disagree ⚠️ (bug)

**Where:** role-classification stage vs the legacy `Major_*` / `Minor_*` columns and the GLUE step, in `Summary.csv`.

**Symptom (sim2, true mix 2a 70 % : 3a 30 %):**

| Field | Value | Actually describes |
|---|---|---|
| `Major_role_subtype` | **2a** (dominance 12.054 > 3a 11.652) — correct | 2a (true major) |
| `Major_genotype_mapping` | **3a** | 3a (true minor) |
| `Reads_nodup_mapped_major` | 231,132 | 3a |
| `Major_consensus_similarity_pct` | 100 | 3a |
| `GLUE_subtype` / NS5A | 3a / M28K | 3a |
| `Reads_nodup_mapped_minor` | 507,708 | 2a |

The dominance/role classifier correctly identifies 2a as the dominant strain, but that reassignment is **not propagated** to the per-strain read-count, coverage, consensus and GLUE/resistance columns, which remain keyed to the *mapping*-major/minor (here 3a). The row therefore contradicts itself: `Major_role_subtype = 2a` paired with 3a's read counts and 3a's resistance result. (sim1, where mapping-major and role-major are both 1a, is internally consistent — the bug only shows when the classifier flips the mapping order.)

**Fix options:** (a) re-key the `Major_*`/`Minor_*`, consensus, and GLUE columns to follow the role classifier so every "major" column refers to the dominant strain (recommended); or (b) keep them mapping-keyed but rename so the split is explicit. Either way, **run GLUE/resistance on the dominant (role-major) strain**, not the mapping-major — otherwise a co-infection can report the minor strain's resistance profile under "major".

---

## 2. `major_typable = NO` for well-covered, corroborated single-infections ⚠️ (likely bug)

**Symptom:** sim11asingle (1a), sim22asingle (2a), sim23asingle (3a) report `major_typable = NO` despite ~690k mapped reads, correct GLUE subtype, and `denovo_major_subtype_match = YES`. The 1b and 4a single-infections report `major_typable = YES`.

**Correlate:** the three `NO` samples all have `Major_cov_breadth_min_10 = NA`, while the `YES` samples have it populated (99.93 / 99.99). Typability appears to depend on a `Major_cov_breadth_min_10` value that is missing for these three, rather than on any real coverage deficit.

**Check:** why is `Major_cov_breadth_min_10` (and `_min_5`) NA for these three samples specifically? Trace whether typability keys off that field and treats NA as fail. These are clean, high-coverage single infections and should be typable.

---

## 3. Resistance reported for only one strain in a co-infection ℹ️ (design check)

`Summary.csv` carries a single GLUE/resistance result per sample, run on the mapping-major consensus. For sim2 that is the 3a strain; the co-infecting 2a strain is not separately typed for resistance. Decide whether a confirmed co-infection should produce a resistance profile for **both** strains. (See also #1 — at minimum it should be the dominant strain.)

---

## 4. `minor_typable` inconsistency across de-novo-refuted minors ℹ️ (minor)

For the spurious minors that de novo refutes in the single-infections, `minor_typable` is `UNKNOWN` for some (sim11asingle 1i, sim22asingle 1m, sim23asingle 1n) and `NO` for others (sim11bsingle 2k1b, sim3 "4"). The outcome (`overall_sample_call = monoinfection`) is correct in all cases, but the `UNKNOWN` vs `NO` split looks arbitrary. Confirm the intended semantics and make it consistent.

---

## Summary

| # | Item | Type | Priority |
|---|---|---|---|
| 1 | Role classifier not propagated to stat/consensus/GLUE columns; GLUE on mapping-major not dominant strain | bug | high |
| 2 | `major_typable = NO` driven by missing `Major_cov_breadth_min_10` on clean samples | bug | high |
| 3 | Only one strain typed for resistance in co-infections | design | medium |
| 4 | `minor_typable` UNKNOWN vs NO inconsistent for refuted minors | polish | low |

**Validated as correct (no action):** genotype/subtype calls, de novo confirm/refute, monoinfection vs co-infection classification, and RAV detection for typed strains.
