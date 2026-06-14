# Phase 8: Dominance Scoring + Strain-Role Classification - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-13
**Phase:** 08-dominance-scoring-strain-role-classification
**Areas discussed:** Dominance score + dominant pick, Co-infection abundance floor, Corroboration + exceptions, Sample call + legacy retire

---

## Dominance score + dominant pick

### What does "dominant" mean structurally?
| Option | Description | Selected |
|--------|-------------|----------|
| Top score, gated | One dominant = highest score that passes major-gate; none passing → indeterminate | ✓ |
| Top score, ungated | Highest score always dominant; gate only flags quality | |
| Absolute threshold | Dominant requires absolute score cutoff; allows zero/multiple | |

### How to combine the three signals?
| Option | Description | Selected |
|--------|-------------|----------|
| Weighted sum, normalized | Normalize/transform each component, weighted sum, largest weight on evenness | ✓ |
| Evenness-gated, reads tie-break | Evenness as primary multiplier; reads/k-mer only rank even candidates | |
| Relative-to-top shares | Express each component as share relative to top candidate | |

### How to quantify breadth-evenness?
| Option | Description | Selected |
|--------|-------------|----------|
| CV over full reference | CV of per-position depth incl. zeros, mapped to 0–1 (1/(1+CV)) | ✓ |
| Fraction-within-band | Fraction of positions within a band of the median | |
| Gini / evenness index | Gini coefficient over per-position depth | |

### How to expose/calibrate the weights?
| Option | Description | Selected |
|--------|-------------|----------|
| Params, calibrated defaults | Named score_weight_* params, defaults calibrated to evidence table | ✓ |
| Hardcoded constants | Calibrated constants at top of summarize.R | |
| You decide during planning | Leave exposure choice to researcher/planner | |

### How does k-mer cov enter when assembly_support='none'?
| Option | Description | Selected |
|--------|-------------|----------|
| Bonus-only, no penalty | Capped positive boost; absence = no boost, never penalty | ✓ |
| NA → contributes 0 | Missing k-mer cov = 0 in the sum | |
| Exclude from dominant pick | k-mer cov only enters corroboration, not dominance | |

### How to pick dominant on near-equal scores?
| Option | Description | Selected |
|--------|-------------|----------|
| Highest wins, deterministic | Tie-break by reads then ref name; dominant label cosmetic | ✓ |
| Flag co-dominant | Surface balanced/co-dominant note within a margin | |
| You decide | Leave to researcher | |

**User's choice:** As marked above.
**Notes:** Read counts alone are contamination-prone (false 4g had more reads than genuine minors), so evenness must dominate. Per-position depth is already available (`cov$X3` from `depth/*.tsv`). Don't demote a genuine low-yield dominant for failing de novo (fall-back rule).

---

## Co-infection abundance floor

### What should the floor be?
| Option | Description | Selected |
|--------|-------------|----------|
| Reuse minRead/minCov | Same major-gate thresholds; de novo is the discriminator | ✓ |
| Relative-to-dominant | ≥ X% of dominant's reads + minRead/minCov | |
| Dominance-score floor | Candidate score must clear a cutoff | |

### Which mapping's coverage?
| Option | Description | Selected |
|--------|-------------|----------|
| Targeted (second) mapping | Per-candidate targeted coverage; same source as evenness score | ✓ |
| First-pass recruitment | Selection-time reads + first-pass breadth | |
| You decide | Researcher picks per wiring | |

### Same thresholds for dominant gate and co-infection floor?
| Option | Description | Selected |
|--------|-------------|----------|
| Same minRead/minCov | One threshold concept; reproduces v1.0 | ✓ |
| Separate lenient floor | coinf_min_read/cov params, lower than dominant gate | |
| You decide | Researcher decides per evidence table | |

**User's choice:** As marked above.
**Notes:** Thresholds cannot separate artefact from genuine (handoff key negative result) — the floor is a minimal presence gate; de novo corroboration does the discrimination ("guilty until corroborated").

---

## Corroboration + exceptions

### What counts as "has genotype-level assembly support"?
| Option | Description | Selected |
|--------|-------------|----------|
| Reuse denovo_* floors | contig len 1000 AND kmer 2.0 AND identity 90 at genotype level | ✓ |
| New corroboration thresholds | Fresh role-specific cutoffs | |
| You decide | Researcher confirms reproduction | |

### Background vs flagged-co-infection when candidate clears floor but no support?
| Option | Description | Selected |
|--------|-------------|----------|
| Asymmetric (de novo-worked) | Background only if de novo worked for the dominant; else keep flagged uncorroborated | ✓ |
| Strict guilty-until-corroborated | Any candidate lacking own support → background | |
| You decide | Researcher confirms vs evidence table | |

### How to encode different-genotype + 1a/1b allow + 2k1b block?
| Option | Description | Selected |
|--------|-------------|----------|
| Port is_valid_minor verbatim | Lift existing logic as post-scoring special-cases | ✓ |
| Config-driven pair lists | coinfection_allow_pairs / block_pairs params | |
| You decide | Researcher decides port vs generalize | |

### How to represent the role reason?
| Option | Description | Selected |
|--------|-------------|----------|
| Coded field + derived sentence | Controlled-vocab role_reason + human sentence | ✓ |
| Coded field only | Just the controlled-vocab column | |
| Human sentence only | Extend review_flag prose | |

**User's choice:** As marked above.
**Notes:** Asymmetric refute rule preserves genuine low-yield/IVT minors and tempers the literal CLASS-02 wording — verifier must check against the asymmetric rule. `is_valid_minor()` port preserves COMPAT-04 exactly.

---

## Sample call + legacy retire

### How to derive the overall sample call (uncorroborated case)?
| Option | Description | Selected |
|--------|-------------|----------|
| 3 values, flag carries nuance | Exactly monoinfection/co-infection/indeterminate; nuance in role_reason | ✓ |
| Add possible_coinfection | 4th value for uncorroborated-only case | |
| You decide | Researcher finalizes derivation | |

### Retire legacy confirmation path now, or keep additive until Phase 9?
| Option | Description | Selected |
|--------|-------------|----------|
| Retire confirmation logic now | Replace classify_minor_denovo/minor_denovo_status; rewire review_flag onto roles | ✓ |
| Keep additive until Phase 9 | Legacy path runs in parallel through Phase 8 | |
| You decide | Researcher decides per review_flag entanglement | |

### Where to write per-candidate roles/scores/reasons?
| Option | Description | Selected |
|--------|-------------|----------|
| Wide Summary + long candidate detail | Summary.csv stays wide; full per-candidate detail in *.candidates.csv | ✓ |
| Long-format Summary | Summary.csv becomes one-row-per-candidate | |
| You decide | Researcher picks the split | |

**User's choice:** As marked above.
**Notes:** Deliberately departs from the Phase 6/7 additive precedent to avoid contradictory confirmation columns in Summary.csv. Phase 9 retains responsibility for Major_*/Minor_* column aliasing + filename migration + the golden regression suite.

---

## Claude's Discretion

- Exact column names/types for score, role, and role_reason fields.
- Normalization constants/transforms (log base, breadth source, k-mer-cov bonus cap).
- New sourceable R helper vs inline in summarize.R.
- Precise measurement of "de novo worked for the dominant".
- Final role_reason controlled vocabulary and derived-sentence wording.

## Deferred Ideas

- 4th `possible_coinfection` sample-call value — rejected to keep the 3-value contract.
- Co-dominant / balanced flag for near-ties — deferred (deterministic tie-break chosen).
- Config-driven coinfection_allow_pairs/block_pairs — deferred in favor of verbatim port.
- Filename migration + Major_*/Minor_* column aliasing + full CI golden regression — Phase 9.
