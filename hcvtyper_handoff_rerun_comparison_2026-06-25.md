# Handoff: re-run comparison (2026-06-25) — ERR1810469 fix + a new same-genotype side-effect

**Date:** 2026-06-25
**Author:** Jon Bråte (FHI)
**Audience:** HCVTyper development team
**Scope:** Compares the latest re-run (`folkehelseinstituttet/hcvtyper 1.3.0`, SRA dir
`…/2026/HCVTyper_SRA_data/`, Summary.csv 2026-06-25 07:08) against the previous run that the current
manuscript tables were built from. Documents (1) the ERR1810469 fix — confirmed working, (2) an
improvement at ERR1810450/451, and (3) a **new side-effect** introduced alongside the fix that the dev
team should review before the manuscript numbers are finalised. **NB:** the result files were still being
copied to the output dir when this was written; Summary.csv is complete, per-position depth files mostly
present. Manuscript was **not** edited pending dev review of the side-effect.

---

## 1. ERR1810469 — FIXED ✓ (target behaviour achieved)

Thomson truth = single 3a.

| | Previous run | **New run** |
|---|---|---|
| `overall_sample_call` | `monoinfection` | **`co-infection (indeterminate dominance)`** |
| Major | **1a** (forced on coverage) | **3a** (matches Thomson) — 190 reads, 2.7× depth, 4.2% breadth@≥10× |
| Minor | 1a (same as major / NA) | **1a** — 1017 reads, 15.1× depth, 77.7% breadth@≥10×, `minor_typable = YES` |
| de novo | dn_min 3a | dn_maj=1a, dn_min=3a (both genotypes assembled) |
| `review_flag` | **NA** (dropped) | **"Dominance ordering uncertain — read-count and k-mer-coverage rankings disagree. Both genotypes reported as present; co-infection vs contamination agnostic. Please review."** |

The reads-vs-k-mer disagreement trigger and the indeterminate-dominance role both fire as designed: the
better-covered 1a no longer steamrolls the genuine (low-coverage) 3a, both genotypes are reported, and the
flag is back. This is exactly the §3/§7.6-Q2 target from the typability proposal.

---

## 2. ERR1810450 / ERR1810451 — subtype now resolved to 4a (improvement), still flagged divergent

| | Previous run | **New run** |
|---|---|---|
| Major call | genotype **"4"** (subtype unresolved) | **"4a"** (subtype resolved) |
| Major reference | `4_JF735135` (generic genotype 4) | `4a_DQ418789` |
| consensus-to-ref identity | 88.2 / 86.4 | 88.6 / 86.3 (essentially unchanged) |
| de novo major subtype | — | `4` (genotype only) → `denovo_major_subtype_match = NO` |
| `review_flag` | (none / divergent) | **"Major subtype conflict between de novo assembly and mapping — possible reference mismatch or highly divergent strain. Please review."** |

Good outcome: the subtype is now correctly called 4a, **and** two independent QC signals still flag the
sample as divergent — the low consensus-to-reference identity (86–88% vs the 91–95% typical of a
well-matched subtype) and the de-novo-vs-mapping subtype conflict. So the divergent-strain detection story
holds; the genotype-4 lineage is simply ~12–14% diverged from the nearest available 4a reference.

**Note (read-splitting, see §3):** these samples now also carry a `4_JF735135` *minor* candidate that
takes ~half the reads (450: major 7,024 vs minor-cand 7,404; 451: 28,630 vs 28,043) — two genotype-4
references competing for the same reads.

---

## 3. NEW SIDE-EFFECT — same-genotype second candidates on monoinfections

The re-run introduces a second candidate **of the same genotype as the major** on a large number of
monoinfection samples (~23 of the single-infection set). The minor is correctly left `minor_typable = NO`
(so the genotype call is unaffected), but the **joint competitive mapping splits reads between the two
same-genotype references**, with two consequences:

### 3a. Coverage loss on clean single infections

Most affected samples are unchanged (the second reference is divergent enough that reads stay with the
major), but where the second reference is close (notably major `1a_EF407457`/`1a_M62321` paired with a
`1c_AY651061` or second-1a minor candidate) the major loses a large fraction of its reads:

| Sample | Major | minor candidate | breadth@≥10× prev → new | Δ |
|---|---|---|---|---|
| ERR1810445 | 1a (EF407457) | 1a/1c | 91.8 → **72.0** | **−19.8** |
| ERR1810484 | 1a (M62321) | 1a/1c | 98.3 → **82.0** | **−16.2** |
| ERR1810479 | 1a (EF407457) | 1a/1c | 25.6 → **9.8** | **−15.8** |
| ERR1810452 | 1a (EF407457) | 1a/1c | 98.7 → 95.9 | −2.8 |

(Most other monoinfections with a same-genotype candidate — e.g. 458/459/466/467/487 with a 1c candidate,
and the 2b set — show **no** coverage change, so the effect is specific to certain close reference pairs.)

### 3b. Displacement of a genuine second-genotype minor

At **ERR1810475** (2b, expected) the previous run reported a genuine **1a** minor (de novo-corroborated
with multi-kb 1a contigs — a real second-genotype signal). In the new run the minor slot is taken by a
**same-genotype `2_JF735117`** candidate (`minor_typable = NO`), and the sample is reported as a 2b
monoinfection. The 1a is still assembled (de novo `dn_min = 1a`), but the genuine second-genotype minor is
**no longer surfaced in the mapping minor call** — displaced by the same-genotype candidate. This is a
sensitivity regression for genuine co-infection/contamination detection.

### Likely cause / suggested check

The candidate-selection step now appears to admit a second reference of the same genotype as the major
into the competitive-mapping set. On a single-genotype sample this (a) splits reads across two
near-identical references (coverage loss) and (b) can out-compete a genuine different-genotype minor for
the minor slot. Suggest restricting the minor candidate to a **different genotype** from the major (as the
co-infection logic intends), or de-duplicating same-genotype references before the joint mapping, so that
a single genotype's reads are not partitioned across two of its own references.

---

## 4. Summary for the manuscript

- ERR1810469 fix is correct and stable → ready to fold into Table 2 (3a major / 1a minor / indeterminate
  dominance / review flag) once the rest is settled.
- ERR1810450/451 → 4a is a clean improvement → ready to fold in (Concordant = Yes; keep the divergent-ref
  QC narrative, which still applies).
- **Hold the full Table 2 numeric refresh** until §3 is resolved: the coverage drops (445/452/479/484) and
  the ERR1810475 minor displacement are side-effects of the same-genotype-candidate behaviour and will
  change again after a fix + re-run. Baking them in now would record transient regressions.
- Simulated tables (1 & 3) unchanged.
