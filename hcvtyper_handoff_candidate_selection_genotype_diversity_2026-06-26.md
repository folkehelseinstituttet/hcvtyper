# Handoff: genotype-diverse candidate selection — root cause (ERR1810475) + proposed design

**Date:** 2026-06-26
**Author:** Jon Bråte (FHI)
**Audience:** HCVTyper development team
**Scope:** Documents (1) the code-level root cause of the ERR1810475 minor-displacement regression first
flagged in `hcvtyper_handoff_rerun_comparison_2026-06-25.md` §3b, confirmed against the latest run
(`folkehelseinstituttet/hcvtyper 1.3.0-g9e0fe57`, SRA dir Summary.csv 2026-06-26 05:31); (2) that the same
mechanism produces the cosmetic same-genotype "Minor" seen on monoinfections including the simulated set;
(3) a proposed redesign (genotype-diverse candidate set) and a critical evaluation of it; and (4) the
validation checks that must pass before it is adopted. Code references are to the **dev branch of
`/home/jon.brate/hcvtyper`**.

---

## 1. Root cause — confirmed in code

ERR1810475 (Thomson truth = single 2b) carries a genuine, de-novo-confirmed genotype-1a minor
(`denovo_minor_subtype = 1a`, 4,734 bp contig, `denovo_minor_ref = 1a_HQ850279`) that the current pipeline
**does not surface in the mapping minor call**. It is reported as a 2b monoinfection. This is **not** a
coverage problem — it is one of the deepest samples in the set (4.96 M raw reads; 2b major at 219,719
dedup reads, ~3,167× depth, 99.98 % breadth @≥10×; 9,549 bp 2b de novo contig). The cause is the
candidate-selection topology.

**`summary/candidates.csv` for ERR1810475 — only two candidates, both genotype 2:**

| rank | candidate_ref | genotype | reads | role | role_reason |
|---|---|---|---|---|---|
| 1 | `2b_D10988` | 2 | 4.46 M | dominant | dominant |
| 2 | `2_JF735117` | 2 | 47,102 | background | `same_genotype_as_dominant` |

The 1a is **not in the candidate list at all** — it survives only as a de novo contig.

### Why — two deliberate design choices in v1.3.0

1. **Fixed 2-slot, genotype-blind selection.** `n_candidates = 2` (`nextflow.config:83`).
   `bin/summarize_mapping_to_all_references.R:129-143` takes the top reference per **subtype**, ranks
   subtypes by first-mapping read count, and keeps the **top-2 subtypes by reads**. The in-code comment
   (`:125-128`) is explicit: the genotype-validity rules were *removed* from selection and moved into
   Phase-8 classification — *"Selection is now mechanical (read-recruitment ranking only)."* Nothing stops
   both slots being filled by the same genotype.

2. **Validity deferred to Phase-8 demotion.** `bin/classify_roles.R` (`is_valid_minor()` :137; D-12
   :394-406) demotes a same-genotype (non-1a/1b) slot-2 to `role = background`,
   `role_reason = same_genotype_as_dominant`, `minor_typable = NO`. This **protects the call** (475 stays a
   correct monoinfection) but only *labels* the redundant candidate — it does not free the slot for a
   genuine different-genotype minor.

3. **Rescue replaces, never adds.** `bin/rescue_evaluation.R` can swap a *mismatched* candidate for a
   de-novo-supported reference, but cannot add a 3rd slot. For 475 both genotype-2 slots are internally
   concordant (mapping gt 2 ≈ dominant de novo 2b), so no rescue trigger fires. The 4,734 bp 1a is
   structurally **unpromotable**: not a top-2 read recruiter, and no mismatched slot to rescue into.

**Net:** the genuine 1a loses the only non-dominant slot to a redundant genotype-2 reference because
selection is read-ranked, capped at 2, and genotype-blind. This is the intended trade-off (robustly
suppress same-genotype duplicates) biting the one real sample where a true low-recruitment second genotype
exists.

---

## 2. Same mechanism on the simulated set (cosmetic there)

Every clean simulated monoinfection fills slot 2 with a redundant same-genotype/recombinant reference,
all correctly demoted — so the genotype calls are unaffected, but the Summary `Minor` column is populated:

| Sample | Major (slot 1) | Slot 2 | role_reason | call |
|---|---|---|---|---|
| sim3 (4a) | `4a_Y11604` | `4_JF735135` (gt 4) | `same_genotype_as_dominant` | monoinfection ✓ |
| sim11asingle (1a) | `1a_HQ850279` | `1i_KJ439772` (gt 1) | `same_genotype_as_dominant` | monoinfection ✓ |
| sim11bsingle (1b) | `1b_EU781827` | `2k1b_AY587845` | `discordant_identity` | monoinfection ✓ |
| sim22asingle (2a) | `2a_D00944` | `1m_KJ439778` (gt 1) | `refuted_denovo` | monoinfection ✓ |
| sim23asingle (3a) | `3a_D17763` | `1n_KJ439775` (gt 1) | `refuted_denovo` | monoinfection ✓ |
| sim1 / sim2 (true co-inf) | 1a / 2a | 1b / 3a | `corroborated` | co-infection ✓ |

**Key difference:** on the simulations the slot-2 same-genotype candidate is **purely cosmetic** (clean
single genotypes — no hidden minor to displace). At ERR1810475 the same slot-2 same-genotype candidate
**actively crowds out a real, de-novo-confirmed 1a**. Same code path; benign on simulations, harmful on
the one real sample with a genuine low-recruitment second genotype.

> **Manuscript note:** the Summary `Minor`/`Minor_reference` columns list a same-genotype reference for
> most monoinfections (real and simulated). Any table/supplementary built from those columns must filter
> on `minor_typable == YES` (or `role`), not on the raw `Minor` column, or it will show spurious
> same-genotype minors.

---

## 3. Proposed design — genotype-diverse candidate set

**Intent (JB):** the candidate set should represent distinct biological entities. After the rescue step
there should be **up to `n_candidates` references from different genotypes** (with the 1a/1b pair allowed,
and the 2k1b recombinant blocked as today). If only one genotype is identified, **one reference** is the
correct output. First-mapping and de novo/BLAST should both be able to nominate different-genotype
candidates.

This is essentially relocating the existing `is_valid_minor()` predicate (`classify_roles.R:137`,
rules 1-3) **upstream** into / after selection, reversing the D-05 decision that made selection
genotype-blind.

### Recommended refinement — enforce diversity *after rescue*, not at first-mapping

Keep nomination **broad and neutral** (first-mapping AND de novo/BLAST free to nominate, as now), and add
a **single genotype-collapse step after rescue** (= JB's step 4). At that step:

- keep ≤ `n_candidates` references spanning **distinct genotypes** (best-**subtype** reference per
  genotype — do **not** collapse to a generic genotype, or 4a→generic "4" subtype resolution is lost);
- apply the `is_valid_minor` allow/block rules (1a/1b through; 2k1b blocked);
- break "which genotype earns a slot" ties using **de-novo support**, not read count alone — this is what
  rescues ERR1810475 when the genuine minor under-recruits;
- **exempt the divergent-strain case**: if a single genotype's best reference is a poor match (low
  consensus-to-ref identity + de-novo-vs-mapping subtype conflict — the ERR1810450/451 signature), allow a
  second *same-genotype* reference for within-genotype resolution rather than forcing a different genotype
  into the slot.

Enforcing after rescue (rather than filtering at first-mapping) preserves the same-genotype /
different-subtype information that the divergent-strain QC and subtype calling still need, keeps the
auditable "what was nominated vs what was judged" separation D-05 was after, and still delivers the
genotype-diverse end state.

---

## 4. Critical evaluation — risks to control

1. **Safety burden shifts to de-novo refutation, and it fails open on low coverage.** Today a clean
   monoinfection's same-genotype slot-2 *cannot* produce a false co-infection. Under genotype-diverse
   selection, every monoinfection's slot-2 becomes its top different-genotype cross-mapping/contamination
   signal, and the only thing preventing a false co-infection is `denovo_confirm` / `refuted_denovo`. That
   refute (`classify_roles.R:379-389`, D-11) demotes an unsupported different-genotype candidate **only if
   the dominant assembled substantially**. On a **low-coverage** monoinfection where the dominant does not
   assemble, the fallback is `uncorroborated_kept` → **co-infection**. So the model trades a structural
   guarantee for a coverage-dependent one and adds a false-positive surface on exactly the low-yield
   (HiSeq/metagenomic) samples that dominate the Thomson dropouts. ERR1810475 itself is safe (high
   coverage; both strains assemble).

2. **"One reference per genotype" can degrade divergent strains (ERR1810450/451).** Those genotype-4 runs
   currently carry two *same-genotype* references (`4a_DQ418789` + generic `4_JF735135`); competitive
   mapping lets each read pick the better-matching genotype-4 reference, which matters when the strain sits
   12-14 % off the nearest reference. Forcing one reference per genotype removes that within-genotype
   per-read choice and could reduce coverage/consensus on the divergent samples the manuscript showcases.
   Hence the divergent-strain exemption in §3.

3. **Re-coupling selection to validity.** D-05 deliberately decoupled these. Applying the constraint at the
   **post-rescue collapse** (not first-mapping) keeps the decoupling benefit (a full record of what was
   nominated before judgement).

4. **De novo as a primary nominator is a real scope increase.** De novo/BLAST is currently used for
   confirmation/rescue, not primary candidate nomination. Letting it nominate different-genotype candidates
   (JB's intent) is the right direction but is more than a one-line change to `denovo_layer.R` /
   `rescue_evaluation.R`.

---

## 5. Validation gate — must all pass before adoption

Re-run the full single-infection set under the new rule and confirm:

- **(a)** all 50 single-infection runs still call `monoinfection` (no new false co-infections — the
   §4.1 low-coverage failure mode);
- **(b)** ERR1810450/451 keep subtype **4a** and their 86-88 % divergent-reference flag (the §4.2
   divergent-strain regression);
- **(c)** ERR1810475 now reports the **1a minor** in the mapping call (the target fix);
- **(d)** simulated monoinfections (sim3, sim11asingle, sim22asingle, sim23asingle) no longer emit a
   same-genotype `Minor`, and sim1/sim2 co-infections are unchanged.

If (a)-(d) hold, the model is good.

---

## 6. Manuscript dependency

Until this is resolved, **hold** the ERR1810475 row of Table 2 and the line-224 Results text: the
biological claim (genuine 1a, multi-kb contig) is still supported by the latest data, but "produced a
minor 1a **mapping**" no longer matches the pipeline. See
`hcvtyper_handoff_rerun_comparison_2026-06-25.md` for the rest of the Table 2 refresh status (469 fix and
450/451→4a are ready; the §3a coverage-loss side-effect was resolved in the 2026-06-26 run).
