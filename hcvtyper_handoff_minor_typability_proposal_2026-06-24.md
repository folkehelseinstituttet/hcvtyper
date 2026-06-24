# Handoff addendum (2026-06-24): ERR1810469 under v1.3.0 joint mapping

**Date:** 2026-06-24
**Author:** Jon Bråte (FHI)
**Scope:** This is a timestamped addendum to `hcvtyper_handoff_minor_typability_proposal.md`.
It describes **only** how the new pipeline results differ from the old ones for the single sample
**ERR1810469**. It does not revisit any other sample or any of the proposal's general logic.

---

## 1. What changed in the pipeline between the two runs

The earlier results (the ones reflected in the current manuscript tables, v1.1.6–v1.2.0 lineage) mapped
reads **independently** against each selected candidate reference (one targeted mapping per candidate).
A read could therefore map to *both* the major and the minor candidate, so per-candidate read counts
were not mutually exclusive and the summed "mapped reads" across candidates could exceed the number of
reads in the sample.

The new run (`folkehelseinstituttet/hcvtyper 1.3.0`) performs **joint competitive mapping**: the selected
candidate references are concatenated into a single index and the reads are mapped together, so each read
is assigned to only its single best-matching reference (read-partitioning). Counts across candidates are
now mutually exclusive.

This is the intended architecture, and on the controlled data it improves quantification (e.g. the
simulated 7:3 co-infections now recover 69.3 % / 68.7 % major fraction, near-exact). ERR1810469 is the
case where the same change, combined with the role-classification logic, produces a **worse** call than
before.

---

## 2. ERR1810469 — old result vs new result

**Ground truth (Thomson 2016):** single **genotype 3a** infection (one of the 26 single-genotype patients).

| | **Old run** (independent mapping; current manuscript) | **New run** (v1.3.0 joint mapping) |
|---|---|---|
| `overall_sample_call` | 3a major + 1a minor, **flagged for review** | **monoinfection** |
| Major call | **3a** (matches Thomson) | **1a** (1a_HQ850279) |
| Major support | 248 nodup reads, 8.8 % breadth @≥10×, 3.5× depth | 1030 nodup reads, 78.4 % breadth @≥10×, 15.2× depth, 94.2 % identity to ref |
| Minor call | 1a (1030 reads, 78.4 % breadth) — surfaced | 3a **demoted to de novo only** (`denovo_minor_subtype = 3a`, 4154 bp contig); `minor_typable = NO` |
| `review_flag` | raised (dominance could not be confidently assigned) | **NA** |

**Net effect:** a call that was **concordant with the published 3a truth and explicitly flagged** became
an **incorrect 1a monoinfection with no flag**. The 3a strain is now only visible as a de novo contig.

---

## 3. Why it flipped — the evidence points to 3a, the call went to 1a

The unbiased first-mapping (all references, before competition) shows 3a is the dominant genotype by raw
recruitment:

| genotype | first-mapping reads (with-dup, top refs) |
|---|---|
| 3a | ~4960 (3a_D17763 = 2175; +D28917 1416; +X76918 1370) |
| 1a | ~1070 (1a_HQ850279 = 662; +AF009606 141; +M62321 111; …) |

So both the **first-mapping recruitment** (3a ≈ 5× the 1a reads) and the **published truth** say 3a.
After candidate selection + competitive mapping, however, the divergent / duplicate-heavy 3a collapses to
248 unique reads at 8.8 % breadth, while the cleanly-matching 1a retains 1030 reads at 78.4 % breadth — so
the coverage-led dominance score crowns **1a** as major, and (because the minor is now same-genotype-as-major
or sub-gate) no second genotype is reported and no flag is raised.

The 3a's pattern — high first-mapping recruitment but low targeted breadth, short de novo contigs at high
k-mer coverage — is the background/contamination signature already described in the parent proposal (§3).
The 1a is plausibly a well-matched contaminant that outranks the genuine-but-poorly-sequenced 3a on
coverage alone.

---

## 4. Implication

This is exactly the dominance-conflict case the parent proposal's **§3 `indeterminate` role** and **§7.6 Q2
reads-vs-k-mer disagreement trigger** are designed to catch:

- reads / breadth / contig length favour **1a**
- first-mapping recruitment / de novo k-mer coverage favour **3a**

Target behaviour (per parent proposal §7.4): both genotypes reported as present, `overall_sample_call =
co-infection (indeterminate dominance)`, `review_flag` raised — never a silent 1a monoinfection. The v1.3.0
run does not yet apply that trigger to ERR1810469: it forces a coverage-led major and drops the flag. The
regression here is in **role assignment + flagging**, not in the joint mapping itself.
