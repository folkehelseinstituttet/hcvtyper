# HCVTyper Handoff — De-novo Rescue Swaps Out the Dominant Strain (2679090-HCV)

- **Date:** 2026-07-01 05:59:02
- **Branch:** dev
- **Sample:** `2679090-HCV` (`2679090-HCVS47L001`)
- **Run:** `NGS_SEQ-20260210-01`
- **Author:** jon.brate@fhi.no (analysis via Claude Code)
- **Status:** Bug diagnosed; fix **drafted, NOT yet applied** (apply via `/gsd-quick`)

## Paths

- Original (v1.1.5) results:
  `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/NGS_SEQ-20260210-01/`
- Dev-branch results:
  `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon/NGS_SEQ-20260210-01/`
- Rescue logic: `bin/rescue_evaluation.R`
- Tests: `bin/tests/test_rescue_evaluation.R`
- Params: `nextflow.config` (lines 63-67, 83), `conf/modules_hcv.config` (line 108)

---

## 1. Symptom

The dev-branch run turned a clean genotype-1 (subtype **1a**) sample into a spurious
**1 + 3a co-infection**, and reported a "95.3% mapped reads" figure for 3a that is
actually genotype 1a's number.

| Field | Original (v1.1.5) | Dev branch |
|---|---|---|
| overall call | 1a / 1b (effectively **mono genotype-1**) | **1 + 3a co-infection** |
| Major reference | 1a_HQ850279 | 1_KJ439780 |
| Major subtype | **1a** | **1** (unresolved) |
| Minor | 1b | **3a** |
| review_flag | none | "assignment uncertain — please review" |

### Evidence trail

**First mapping** (`parsefirstmapping/2679090-HCVS47L001.candidates.csv`):
- cand 1 = `1a_HQ850279`, **6,534,184 reads**, cov 100
- cand 2 = `1_KJ439780`, 55,438 reads, cov 91

**After rescue** (`blastparse/2679090-HCVS47L001.rescued.candidates.csv`):
- cand 1 = `3a_X76918`, **6,534,184 reads, cov 100** (← 1a's numbers, carried over),
  `rescued_from=1a_HQ850279`, trigger: `denovo 3a_X76918 contig 4189bp pident=92.357 aln=4187bp kmer_cov=2.72797 (replaced 1a_HQ850279)`
- cand 2 = `1_KJ439780`, 55,438 reads, cov 91

**Assembly support** (`blastparse/2679090-HCVS47L001.assembly_support.csv`):

| subtype | best_ref | contig len | pident | aln | **kmer_cov** |
|---|---|---|---|---|---|
| **1a** (own) | 1a_HQ850279 | 2118 | 92.30 | 2118 | **31647×** |
| 1b | 1b_EU781827 | 1824 | 90.74 | 1804 | 1.26 |
| 3a (alt) | 3a_X76918 | 4189 | 92.36 | 4187 | **2.73×** |

**Final summary** (`summary/candidates.csv`) — stale first-map numbers next to real targeted numbers:

| rank | ref | candidate_reads (first-map, stale) | targeted_reads_nodup (real) | role |
|---|---|---|---|---|
| 1 | 3a_X76918 | 6,534,184 | **3,180** | co-infection (corroborated) |
| 2 | 1_KJ439780 | 55,438 | **97,242** | dominant |

The "95.34%" that `Summary.csv` shows as `percent_mapped_reads_minor_firstmapping`
for 3a is identical to 1a's original `percent_mapped_reads_major_firstmapping`
(`6,534,184 / 6,853,211 = 95.3%`). It followed the label swap, not the biology.

Real 3a signal: `Percent_reads_mapped_of_trimmed_with_dups_minor = 0.36%`,
25,727 with-dup / 3,180 nodup reads, breadth min-10 = 61%, consensus similarity 89.9%.
Genotype 1 captured ~95% of all reads; the true 3a reference captures ~0.36%.

---

## 2. Root cause

### 2a. Why a 6.5M-read / 95%-cov candidate got swapped out

`evaluate_row()` in `bin/rescue_evaluation.R` (lines 190-263) decides rescues using
**only** `assembly_support.csv` (de-novo contig metrics). It never reads
`candidate_reads` / `candidate_cov`. So the mapping dominance of 1a is invisible to
the rescue. With configured floors `rescue_min_length=3000, rescue_min_pident=85,
rescue_min_aln_length=3000, rescue_min_kmer_cov=2.0`:

1. **Own-subtype confirmation FAILS on length** (lines 206-209): 1a's contig is
   2118 bp < 3000, so `floors_ok` returns FALSE — 1a is treated as "unconfirmed"
   and becomes eligible for replacement, despite kmer_cov 31,647× (unmistakably the
   dominant strain).
2. **Best alternative = long 3a** (lines 215-218): 4189 bp.
3. **3a passes all four floors** (line 256): 4189>3000, 92.36>85, 4187≥3000,
   kmer_cov 2.73≥2.0 (clears by 0.7×). Rescue fires; cand 1 rewritten 1a→3a, 1a's
   6.5M reads / 100% cov carried over unchanged (lines 317-322).

The floors reward contig LENGTH and ignore an ~11,000× k-mer-coverage disparity
(31,647 vs 2.73). Contig length is not monotonic with read depth, so a genuine
dominant strain with a fragmented assembly is declared "unconfirmed."

### 2b. Why 1a couldn't come back

After the overwrite, additive nomination (lines 383-433) only surfaces genotypes
**not already present**. Candidate genotypes were then {3, 1} — genotype 1 was
"present" via the generic `1_KJ439780` — so nothing re-nominated the proper 1a
reference. The subtype-resolved 1a was lost; genotype 1 survived only as the coarse
`1_KJ439780`.

### 2c. The "true bug" in candidates.csv numbers

On REPLACE, `candidate_ref` is overwritten but `candidate_reads` / `candidate_cov`
are **not reset** — they still describe the displaced reference's first mapping.
That is why 3a displays 1a's 6.5M / 100% / 95.3%. The nomination path already emits
`NA_real_` for these; the replace path should too.

### 2d. Implication for the final call

- A false **3a co-infection** was invented from a thin cross-mapping contig
  (kmer_cov 2.7×) — exactly the artefact the milestone is meant to reject.
- The dominant strain was **de-resolved** 1a → generic "1" (1_KJ439780), major
  consensus similarity dropping to 83.6%.
- The minor changed **1b → 3a** (different genotype, not just subtype).
- Partial self-correction: role logic uses *real* targeted reads, so 1_KJ439780
  (97,242) was promoted to Major over the mislabelled 3a (3,180), and the
  `review_flag` fired. But the sample is still auto-reported as a co-infection.

---

## 3. Key finding that reshaped the fix

The existing test suite (`bin/tests/test_rescue_evaluation.R`) **Subtest 1
("rescue-fires")** deliberately encodes a cross-genotype REPLACE as *correct*: a
lone `1a` candidate with **no** `1a` assembly support, where de novo says `2b`, is
*supposed* to be remapped `1a → 2b` (first-mapping-got-the-genotype-wrong
correction).

So a blunt "block all cross-genotype replacement" (original option #1) breaks
intended behaviour. The real discriminator between subtest 1 (replace is right) and
2679090-HCV (replace is wrong) is **whether the candidate's own subtype is
assembled**: subtest 1 has no `1a` contig; 2679090 has one (kmer 31,647×), just
fragmented below the length floor. All guards below are scoped around that
discriminator so they are **non-breaking** against the 13 existing subtests
(guards are inert when the extra args are absent or the own subtype is unassembled).

---

## 4. Drafted fix (NOT yet applied)

### 4.0 Stale first-mapping stats — reset on REPLACE

`bin/rescue_evaluation.R`, insert immediately before line 485 (`out_cols <- ...`):

```r
# --- Stale first-mapping-stat guard ------------------------------------------
# A REPLACED candidate's first-mapping read count / coverage described the
# DISPLACED reference, not the reference now shown in candidate_ref. Blank them
# (as the additive-nomination path already does) so no downstream table reports
# the old ref's reads / first-mapping % against the new ref. Targeted mapping
# recomputes the real per-candidate numbers.
out <- out %>%
  mutate(
    candidate_reads = if_else(!is.na(rescued_from), NA_real_, candidate_reads),
    candidate_cov   = if_else(!is.na(rescued_from), NA_real_, candidate_cov)
  )
```

Placed after collapse/re-rank so ranking still uses the real first-map reads.

### 4.1 New parameters + wiring

`nextflow.config`, after line 67:

```groovy
    rescue_kmer_cov_ratio       = 10        // Fix #3: block a REPLACE when the candidate's own-subtype contig k-mer coverage exceeds the rescue target's by >= this factor (cross-mapping-noise guard). 0/blank => disabled.
    rescue_dominant_protect_cov = 90        // Fix #2: a first-mapping candidate at >= this coverage whose own subtype IS assembled (identity+k-mer floors met) is self-confirming and never replaced. Blank => disabled.
```

`conf/modules_hcv.config`, line 108 (append 11th/12th positional args):

```groovy
        ext.args   = { "${params.rescue_min_length} ${params.rescue_min_pident} ${params.rescue_min_aln_length} ${params.rescue_min_kmer_cov} ${params.rescue_1a1b_length} ${params.n_candidates} ${params.rescue_kmer_cov_ratio} ${params.rescue_dominant_protect_cov}" }
```

`bin/rescue_evaluation.R`, after line 83 (`n_candidates_cap` block):

```r
# Optional 11th/12th positional args (Fix #2/#3). Absent => guard disabled, so
# the 9-arg subprocess test and any legacy caller preserve prior behaviour.
rescue_kmer_cov_ratio <- if (length(args) >= 11) suppressWarnings(as.numeric(args[11])) else NA_real_
if (is.na(rescue_kmer_cov_ratio) || rescue_kmer_cov_ratio <= 0) rescue_kmer_cov_ratio <- Inf
dominant_protect_cov  <- if (length(args) >= 12) suppressWarnings(as.numeric(args[12])) else NA_real_
if (is.na(dominant_protect_cov)) dominant_protect_cov <- Inf
```

### 4.2 Fix #2 (primary) — mapping-aware own-subtype confirmation

`bin/rescue_evaluation.R`, replace lines 200-210 inside `evaluate_row()`:

```r
  sample_support <- support %>% filter(sample == sample_id)
  if (nrow(sample_support) == 0) return(no_rescue)

  # Best own-subtype contig depth — reference point for Fix #2 / #3.
  own_rows <- sample_support %>% filter(subtype == cand_sub)
  own_kmer <- if (nrow(own_rows) > 0) max(own_rows$best_contig_kmer_cov, na.rm = TRUE) else NA_real_

  # Own-subtype confirmation. Skipped for 2k1b (D-03 fires regardless).
  if (!isTRUE(cand_sub == "2k1b") && nrow(own_rows) > 0) {
    # Standard: longest own-subtype contig passes all four floors.
    own_longest <- own_rows %>% arrange(desc(best_contig_length)) %>% slice(1)
    if (floors_ok(own_longest, rescue_min_length)) return(no_rescue)
    # Fix #2 — mapping-aware confirmation. A candidate already well covered by
    # first-mapping reads is a real dominant strain; its de-novo contig merely
    # assembling SHORT must not make it eligible for replacement. Confirm if the
    # own subtype clears the QUALITY floors (identity + k-mer depth) even when it
    # fails the LENGTH / ALN floors. (Subtest 1 stays a replace: there the own
    # subtype has NO contig, so this cannot fire.)
    own_best_q <- own_rows %>% arrange(desc(best_contig_kmer_cov)) %>% slice(1)
    own_quality_ok <- isTRUE(
      own_best_q$best_contig_pident[1]   >  rescue_min_pident &
      own_best_q$best_contig_kmer_cov[1] >= rescue_min_kmer_cov
    )
    if (isTRUE(row$candidate_cov >= dominant_protect_cov) && own_quality_ok)
      return(no_rescue)
  }
```

Effect: 1a (cov 100, own contig pident 92.3 / kmer 31,647×) is confirmed → not
replaced. Genotype 3 is then still absent, so the **existing nomination path**
re-surfaces 3a additively (NA reads), and collapse keeps
`{1a_HQ850279 (Major), 3a_X76918 (minor)}`. **1a preserved as Major — core bug gone.**

Note: this + the existing nomination path already achieve original option #1's goal
(keep the dominant, surface the de-novo strain additively) without breaking
subtest 1 — so the standalone cross-genotype block (#1) is dropped.

### 4.3 Fix #3 (backstop) — relative k-mer-coverage guard on REPLACE only

Add inside `evaluate_row()` after the `no_rescue` / `own_kmer` setup:

```r
  # Fix #3 — a REPLACE target whose k-mer depth is dwarfed by the candidate's own
  # assembled subtype is cross-mapping noise sitting on a real strain; refuse.
  # Inert when the own subtype has no contig (own_kmer NA) — preserves subtest 1.
  kmer_ratio_blocks <- function(target_kmer) {
    !is.na(own_kmer) && is.finite(own_kmer) && !is.na(target_kmer) &&
      target_kmer > 0 && (own_kmer / target_kmer) >= rescue_kmer_cov_ratio
  }
```

2k1b path — line 235, add to the condition:

```r
        if (!is.na(rescue_ref) && rescue_ref != orig_ref &&
            !kmer_ratio_blocks(g2_row$best_contig_kmer_cov[1])) {
```

Normal four-floor path — after line 256:

```r
  if (!floors_ok(alt_support, length_floor)) return(no_rescue)
  if (kmer_ratio_blocks(alt_support$best_contig_kmer_cov[1])) return(no_rescue)  # Fix #3
```

For 2679090: own 1a 31,647× vs target 3a 2.73× → ratio ≈ 11,600 ≥ 10 → blocked.
**Applies to the REPLACE path only — never nomination** (a relative-kmer guard on
nomination would break subtest 12 / ERR1810475, which nominates a thin 7.1× 1a
contig under a deep 30,307× 2b dominant).

### 4.4 New unit tests

`bin/tests/test_rescue_evaluation.R` — add `extra_args` to `run_rescue()`:

```r
run_rescue <- function(case, prefix, cands, support, refs, prestage = character(0), extra_args = NULL) {
  ...
  exit <- system2("Rscript",
    c(shQuote(script), shQuote(prefix), shQuote(cand_path), shQuote(support_path),
      shQuote(refs_path),
      TH_LENGTH, TH_PIDENT, TH_ALN, TH_KMER, TH_1A1B_LEN,
      if (is.null(extra_args)) character(0) else extra_args),
    stdout = FALSE, stderr = FALSE)
  ...
}

REFS_ERR <- c(REFS, list("1a_HQ850279" = strrep("G",60),
                         "3a_X76918"   = strrep("T",60),
                         "1_KJ439780"  = strrep("A",60)))

# Subtest 14: dominant strain protected (2679090-HCV shape).
rA <- run_rescue("err_dominant", "ERRDOM",
  cands = mk_cands("ERRDOM",
    mk_cand(1, "1a_HQ850279", "1a", "1", 6534184, 100, "pass"),
    mk_cand(2, "1_KJ439780",  "1",  "1",   55438,  91, "pass")),
  support = mk_support("ERRDOM",
    mk_support_row("1a", "1a_HQ850279", 2118, 92.304, 2118, 31647.19),
    mk_support_row("3a", "3a_X76918",   4189, 92.357, 4187, 2.728)),
  refs = REFS_ERR, extra_args = c("2", "10", "90"))
assert_schema("err_dominant", rA$cands_out)
a1 <- rA$cands_out %>% filter(candidate_ref == "1a_HQ850279")
if (nrow(a1) != 1 || !is.na(a1$rescued_from[1]))
  fail("err-dominant: dominant 1a must survive UNREPLACED (Fix #2)")
if (!any(rA$cands_out$candidate_ref == "3a_X76918"))
  fail("err-dominant: 3a must be surfaced additively via nomination")
if (any(rA$cands_out$candidate_ref == "1_KJ439780"))
  fail("err-dominant: redundant genotype-1 slot must collapse")
if ((rA$cands_out %>% arrange(candidate_rank) %>% slice(1))$candidate_ref[1] != "1a_HQ850279")
  fail("err-dominant: Major (rank 1) must remain 1a_HQ850279")
ok("dominant-protect -> well-covered assembled 1a kept as Major; thin 3a additive (2679090-HCV)")

# Subtest 15: Fix #3 in isolation (below the #2 cov floor).
rB <- run_rescue("kmer_guard", "KMERG",
  cands = mk_cands("KMERG",
    mk_cand(1, "1a_HQ850279", "1a", "1", 500000, 50, "pass")),
  support = mk_support("KMERG",
    mk_support_row("1a", "1a_HQ850279", 2118, 92.3,   2118, 31647.19),
    mk_support_row("3a", "3a_X76918",   4189, 92.357, 4187, 2.728)),
  refs = REFS_ERR, extra_args = c("2", "10", "90"))
b1 <- rB$cands_out %>% filter(candidate_ref == "1a_HQ850279")
if (nrow(b1) != 1 || !is.na(b1$rescued_from[1]))
  fail("kmer-guard: 1a REPLACE must be blocked by the relative k-mer guard (Fix #3)")
ok("kmer-guard -> own-subtype depth >> target depth blocks the REPLACE (Fix #3)")

# Subtest 16: legitimate correction still fires AND stale stats are blanked.
rC <- run_rescue("stale_reads", "STALEREADS",
  cands = mk_cands("STALEREADS",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("STALEREADS",
    mk_support_row("2b", "2b_AY232748", 4120, 96.2, 3840, 4.1)),
  refs = REFS, extra_args = c("2", "10", "90"))
c1 <- rC$cands_out %>% filter(candidate_rank == 1)
if (is.na(c1$rescued_from[1]))
  fail("stale-reads: legitimate 1a->2b correction must still fire with guards enabled")
if (!is.na(c1$candidate_reads[1]) || !is.na(c1$candidate_cov[1]))
  fail("stale-reads: a REPLACED candidate must not carry the displaced ref's reads/cov")
ok("stale-reads -> replace still fires; displaced ref's first-mapping reads/cov blanked")
```

---

## 5. Expected behaviour after the fix

| | Before | After |
|---|---|---|
| Candidate set | {3a_X76918, 1_KJ439780} | {**1a_HQ850279** (Major), 3a_X76918 (minor)} |
| Major subtype | "1" (unresolved) | **1a** (recovered) |
| 3a's reads/cov/first-map% | 6.5M / 100 / 95.3% (stale) | **NA** (never first-mapped) |
| subtests 1-13 | pass | **pass** (guards inert unless args passed / own subtype assembled) |

The spurious 3a *minor* nomination itself remains (same thin-contig shape subtest 12
legitimately keeps) — a downstream role/concordance calibration question, separate
from this bug, and already what triggers `review_flag`.

---

## 6. Next steps

1. Apply via `/gsd-quick` (creates planning artifact, atomic commits, runs nf-test).
2. Run the R contract test:
   `Rscript bin/tests/test_rescue_evaluation.R` → expect `ALL PASS`.
3. Re-run the module nf-test: `modules/local/rescueevaluation/tests/main.nf.test`
   (uses `~/.nf-test` 0.9.3, `--profile docker`; mind host disk — see memory
   "Host disk near-full").
4. Update the nf-core schema for the two new params
   (`nf-core pipelines schema build` in the NEXTFLOW conda env) — deferred/lint.
5. Re-run 2679090-HCV end-to-end and confirm Major=1a, Minor=3a additive with NA
   first-map stats and the review_flag still raised.

## 7. Open questions

- Is `rescue_dominant_protect_cov = 90` the right floor? Tune against the benchmark
  panel (must not resurrect the failure modes the milestone closed).
- Should the spurious 3a *minor* be suppressible via a concordance rule (targeted
  nodup reads vs first-map), independent of the nomination floors?
