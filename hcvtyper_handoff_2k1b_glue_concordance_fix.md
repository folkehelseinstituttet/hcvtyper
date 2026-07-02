# Handoff: Fix 2k/1b samples wrongly failing the Phase-8 concordance gate (empty MultiQC genotype)

**Branch:** `dev` (current tip when this was raised: `a953f61`) · **Date raised:** 2026-07-02 · **Severity:** Release blocker (silently zeroes out every genuine 2k/1b call)

This file is self-contained — it can be worked from on another machine without the originating chat.

---

## 1. Symptom

TEST run: `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/TEST/`
(samples `2143147-HCV`, `2144338-HCV`, `2302582-HCV` — three samples previously identified as
the HCV 2k/1b recombinant).

All three now come out `overall_sample_call = indeterminate` with an **empty `Genotype` column
in the MultiQC report** (`multiqc/multiqc_data/multiqc_summary.txt`) and
`review_flag = "No candidate passed the major-gate — overall sample call indeterminate."`

The 2k/1b call itself is **not lost** — both mapping and de novo still agree on it:

| Sample | `Major_genotype_mapping` | `Major_reference` | `denovo_major_subtype` | `denovo_major_subtype_match` | `GLUE_genotype` / `GLUE_subtype` |
|---|---|---|---|---|---|
| 2143147-HCV | 2k1b | 2k1b_AY587845 | 2k1b | YES | 1 / 1b |
| 2144338-HCV | 2k1b | 2k1b_AY587845 | 2k1b | YES | 1 / 1b |
| 2302582-HCV | 2k1b | 2k1b_AY587845 | 2k1b | YES | 1 / 1b |

`summary/candidates.csv` shows the mechanism directly:
`concordance_status = discordant`, `concordance_reason = discordant_mapping_vs_glue`,
`role = background`, `role_reason = discordant_identity`, `below_floor = TRUE`.

---

## 2. Root cause

**File:** `bin/classify_roles.R`, function `apply_concordance()` (lines 70–124).

This is the Phase-8 "concordance gate" added in commits `6f58ae8` (`feat(concordance): add
apply_concordance()...`) and `6d33abf` (`feat(classify): wire concordance_status gate...`),
2026-06-16. It compares three identity legs at genotype level — mapping, GLUE, de novo — and
marks a candidate `discordant` if any two legs disagree:

```r
map_gt   <- as.character(df$candidate_genotype)        # mapping-derived, "2k1b"-aware
glue_gt  <- as.character(df$candidate_glue_genotype)    # from HCV-GLUE, always "1" or "2" for 2k/1b
...
glue_conflict   <- has_glue   & !is.na(glue_gt)   & map_gt != glue_gt
denovo_conflict <- has_denovo & !is.na(denovo_gt) & map_gt != denovo_gt
```

A `discordant` candidate is then excluded from the eligible pool for the major-gate
(`classify_roles.R:317-341`, `concordance_ok <- ... concordance_status != "discordant"`). With
the sole candidate discordant, **no candidate is eligible**, so `classify_roles()` returns
`overall_sample_call = "indeterminate"` (D-14 fallback ladder) and the legacy `Major`/`Minor`
alias columns stay `NA` (`summarize.R:1655`, `Major = Major_subtype`). MultiQC's `Genotype`
column is literally `Major` (`summarize.R:1738`, `Genotype = Major`), so it renders empty.

**Why this isn't a one-off — it is structural and will hit every 2k/1b sample:**
HCV-GLUE's genotype call comes from phylogenetic clade placement
(`json[["phdrReport"]][["samReferenceResult"]][["genotypingResult"]][["genotypeCladeCategoryResult"]][["shortRenderedName"]]`,
parsed in `bin/GLUE_json_parser.R:104-107`). GLUE's reference tree has **no CRF_02k/1b clade** —
recombination detection across gene regions is not part of its genotyping model. A genuine 2k/1b
recombinant (2k from 5′UTR–NS2, 1b from NS3–3′UTR — the 1b portion is the majority of genome
length) is therefore always placed in genotype **1** (subtype **1b**), sometimes genotype 2,
never `2k1b`. So `map_gt == "2k1b"` vs `glue_gt %in% c("1","2")` will be flagged `discordant`
for **every** correctly-called 2k/1b sample, not just these three.

**This directly violates the CLAUDE.md Compatibility constraint** ("Existing 1a/1b and 2k/1b
exceptions must be preserved"). That constraint was honoured for **co-infection pairing**
(`is_valid_minor()`, `classify_roles.R:126-150` — blocks 2k1b paired with genotype 1/2/2k1b) and
for the **rescue** path (D-03 in `rescue_evaluation.R:256-282`), both of which already have
2k1b-aware special-casing. `apply_concordance()` was added later (2026-06-16) as a *third*
identity-comparison surface and was never given the same exception — that's the gap.

**Not implicated:** the de novo leg is fine as-is. `denovo_gt` is computed via the
already-2k1b-aware `genotype_from_subtype()` (single-sourced in `bin/genotype_utils.R`), so when
a de novo contig genuinely matches a 2k1b reference, `denovo_gt == "2k1b" == map_gt` — no
conflict. Confirmed by the data above: `denovo_major_subtype_match = YES` on all three samples.
Only the GLUE leg needs the exception.

**Also not weakened by this fix:** the D-11 de novo-refutation safety net (the pipeline's
headline "never reported when de novo refutes it as a cross-mapping artefact" capability) is a
completely separate mechanism — it runs on `.own_substantial` (the candidate's own de novo
support), computed independently of `apply_concordance()` (`classify_roles.R:356-386`). A
spurious 2k1b **mapping** artefact with no genuine 2k1b de novo contig behind it will still be
correctly refuted via `refuted_denovo`, regardless of this fix.

---

## 3. The fix (surgical, mirrors the existing 2k1b exception pattern)

**Edit `bin/classify_roles.R`, inside `apply_concordance()`** (around current lines 95–121).

Add a GLUE-only 2k1b exemption before computing `glue_conflict`, and thread a distinct
`concordance_reason` through so the exception is auditable in `candidates.csv` (not silently
indistinguishable from a real 3-leg string match):

```r
  map_gt   <- as.character(df$candidate_genotype)
  glue_gt  <- as.character(df$candidate_glue_genotype)
  denovo_gt <- ifelse(
    has_denovo,
    vapply(df$assembly_support_subtype, function(s) {
      if (is.na(s) || !nzchar(s)) NA_character_ else as.character(genotype_from_subtype(s))
    }, character(1L)),
    NA_character_
  )

  # 2k1b structural exception (CLAUDE.md Constraints): HCV-GLUE's clade-placement
  # tree has no CRF_02k/1b category, so a genuine 2k/1b recombinant is ALWAYS
  # reported by GLUE as genotype 1 or 2 (whichever region/majority-length portion
  # dominates the consensus), never "2k1b". This is expected GLUE behaviour, not
  # evidence of a wrong mapping/de novo call. Mirrors the 2k1b-aware exception
  # already applied to co-infection pairing in is_valid_minor() (below) and to
  # the D-03 rescue rule in rescue_evaluation.R. GLUE leg only — the de novo leg
  # already agrees natively via genotype_from_subtype()'s 2k1b-aware rule.
  glue_2k1b_exempt <- map_gt == "2k1b" & glue_gt %in% c("1", "2")

  glue_conflict   <- has_glue   & !is.na(glue_gt)   & map_gt != glue_gt & !glue_2k1b_exempt
  denovo_conflict <- has_denovo & !is.na(denovo_gt) & map_gt != denovo_gt

  status <- character(nrow(df))
  reason <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    if (glue_conflict[i] || denovo_conflict[i]) {
      status[i] <- "discordant"
      reason[i] <- if (glue_conflict[i] && denovo_conflict[i]) {
        "discordant_all_legs"
      } else if (glue_conflict[i]) {
        "discordant_mapping_vs_glue"
      } else {
        "discordant_mapping_vs_denovo"
      }
    } else if (has_glue[i] && has_denovo[i]) {
      status[i] <- "confirmed"
      reason[i] <- if (glue_2k1b_exempt[i]) "confirmed_2k1b_recombinant" else "all_legs_concordant"
    } else if (has_glue[i] || has_denovo[i]) {
      status[i] <- "unconfirmed"
      reason[i] <- if (has_glue[i]) {
        if (glue_2k1b_exempt[i]) "two_legs_2k1b_recombinant_glue_only" else "two_legs_glue_only"
      } else {
        "two_legs_denovo_only"
      }
    } else {
      status[i] <- "unconfirmed"
      reason[i] <- "no_corroborating_legs"
    }
  }
```

That's the entire diff. No changes needed to `classify_roles.R`'s gate logic
(`concordance_ok`/`eligible`, lines 317-341) — `status[i]` still resolves to `"confirmed"` or
`"unconfirmed"`, both of which already pass the `!= "discordant"` eligibility check.

**Deliberately scoped narrow:** the exemption only fires for `glue_gt %in% c("1","2")`, not any
GLUE genotype — a 2k1b mapping call disagreeing with, say, GLUE genotype 3 or 4 would be
biologically implausible and should still surface as `discordant` for review.

---

## 4. gsd-quick implementation plan

This is a single self-contained bugfix confined to one function in one file plus its unit tests
— a good fit for `/gsd-quick` rather than a full `/gsd-plan-phase` cycle. Suggested invocation
and scope:

```
/gsd-quick fix 2k1b samples failing the concordance gate against GLUE's genotype call —
apply_concordance() in bin/classify_roles.R needs a GLUE-only 2k1b exemption, see
hcvtyper_handoff_2k1b_glue_concordance_fix.md for full root cause + exact diff
```

### Task breakdown (what the gsd-quick executor should do)

1. **Apply the fix** — edit `bin/classify_roles.R::apply_concordance()` exactly as in §3 above
   (add `glue_2k1b_exempt`, gate `glue_conflict` on it, thread the two new reason strings
   through the `status`/`reason` loop). No other file needs a code change.

2. **Extend the unit test** — `bin/tests/test_classify_roles.R`. Test 11 (lines 333-389)
   already unit-tests `apply_concordance()` directly with a `mk_conc()` builder. Add a **Test
   12** immediately after it, same style:
   - `S6`: `mk_conc("S6", "2k1b", "1", "supported", "2k1b")` → expect
     `status == "confirmed"`, `reason == "confirmed_2k1b_recombinant"` (this is the exact shape
     of the three TEST-run samples: mapping 2k1b + GLUE gt1 + de novo 2k1b).
   - `S7`: `mk_conc("S7", "2k1b", "2", "none", NA)` → expect `status == "unconfirmed"`,
     `reason == "two_legs_2k1b_recombinant_glue_only"` (GLUE says gt2 instead of gt1 — still
     exempt; no de novo leg present).
   - `S8`: `mk_conc("S8", "2k1b", "3", "supported", "2k1b")` → expect `status == "discordant"`,
     `reason == "discordant_mapping_vs_glue"` (GLUE genotype 3 is NOT in the exemption list —
     confirms the exemption is narrowly scoped, not "always exempt when map_gt is 2k1b").
   - Regression guard: re-assert Test 11's `S1`/`S2` (non-2k1b subtypes) are unaffected —
     already covered by re-running the existing test, just confirm it still passes unchanged.

3. **Add one end-to-end test** covering the full `apply_concordance()` → `classify_roles()`
   chain (Test 11 only unit-tests the helper in isolation; nothing today proves the gate
   actually admits a 2k1b candidate as dominant). Add as **Test 13**, modelled on Test 6's
   `mk_cand()` + `classify()` pattern (lines ~202-212) but for a *single* 2k1b candidate with a
   GLUE leg attached (extend `mk_cand()` with an optional `glue_gt = NA` parameter, or build the
   row with `mk_cand()` then `mutate(candidate_glue_genotype = "1")`, then run
   `apply_concordance()` before `classify()`):
   - One candidate: `subtype = "2k1b"`, reads/cov well above the 500/2.0 defaults, own
     substantial de novo support, `candidate_glue_genotype = "1"`.
   - Before the fix: `overall_sample_call == "indeterminate"`, `role == "background"`,
     `role_reason == "discordant_identity"`.
   - After the fix: `overall_sample_call == "monoinfection"`, `role == "dominant"`.
   - This is the test that would have caught the TEST-run regression directly — worth keeping
     even though Test 12 covers the `apply_concordance()` unit separately.

4. **Run the test file:**
   ```
   Rscript bin/tests/test_classify_roles.R
   ```
   Expect `ALL PASS` (currently 11 tests; will be 13 after this change).

5. **Check for nf-test snapshot fallout.** `apply_concordance()`/`classify_roles.R` are sourced
   inline by the `SUMMARIZE` process (staged as a `path` input, not their own module), so the
   only nf-test snapshot that could be affected is:
   ```
   modules/local/summarize/tests/main.nf.test.snap
   ```
   Check whether `modules/local/summarize/tests/data/` contains a 2k1b + GLUE-gt1/2 fixture
   scenario (a quick `grep -rl 2k1b modules/local/summarize/tests/data/` — at last check it did
   **not**, only two incidental `2k1b` mentions in `blast/Test_*_blast_out.csv` unrelated to
   this candidate). If the fixture doesn't exercise this path, the snapshot should be
   byte-identical and does **not** need regeneration. If it does, regenerate per the project's
   nf-test workflow (`~/.nf-test` pinned 0.9.3, invoke via PATH prefix + `--profile docker` —
   see memory `nf-test-local-install`) and diff the snapshot change is confined to the expected
   rows before accepting it.

6. **End-to-end verification against the real regression.** Re-run (or ask to re-run) the TEST
   samplesheet (`/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/TEST/samplesheet.csv`,
   3 samples) and confirm in the new `summary/Summary.csv` / `summary/candidates.csv`:
   - `2143147-HCV`, `2144338-HCV`, `2302582-HCV` → `overall_sample_call = monoinfection`,
     `Major_subtype = 2k1b` (or whatever the current subtype-vs-genotype column naming is at
     fix time), `concordance_status = confirmed`, `concordance_reason = confirmed_2k1b_recombinant`.
   - `multiqc/multiqc_data/multiqc_summary.txt` → `Genotype` column populated (`2k1b`) for all
     three, `review_flag` empty/OK instead of the major-gate failure sentence.
   - Spot-check a handful of **non**-2k1b samples from the same or another recent run to confirm
     `concordance_status`/`overall_sample_call` are byte-identical to before the fix (this
     change should be a no-op for every sample where `candidate_genotype != "2k1b"`).

### Files touched (summary)

| File | Lines | Change |
|---|---|---|
| `bin/classify_roles.R` | ~85–121 (inside `apply_concordance()`) | Add `glue_2k1b_exempt`, gate `glue_conflict` on it, add two new `concordance_reason` values |
| `bin/tests/test_classify_roles.R` | append after line 389 | Test 12 (`apply_concordance()` 2k1b cases) + Test 13 (end-to-end `classify_roles()` 2k1b-now-dominant case) |
| `modules/local/summarize/tests/main.nf.test.snap` | — | Regenerate ONLY if the fixture data is confirmed to exercise a 2k1b+GLUE candidate; check first, likely a no-op |

No changes needed to: `bin/genotype_utils.R` (already 2k1b-aware, reused as-is),
`bin/GLUE_json_parser.R` (GLUE's genotype-1/1b call for 2k/1b is correct GLUE behaviour, not a
parsing bug), `bin/rescue_evaluation.R` / the `is_valid_minor()` co-infection-pairing exception
(unaffected — different identity-comparison surface, already correct).

---

## 5. CLAUDE.md update (blocked for agents — do manually)

The general principle this bug exposes ("HCV-GLUE structurally cannot call `2k1b`; every new
identity-comparison surface must carry the same exception the codebase already applies in
`is_valid_minor()` and the D-03 rescue rule") should be captured in CLAUDE.md's Constraints
section so it isn't rediscovered the next time a strain-identity check is added.

**This could not be applied directly in this session** — a runtime "ARS scope guard" hook
hard-blocks any agent from writing `CLAUDE.md` ("part of the enforcement infrastructure and may
not be written by any agent"). The exact bullet to paste in manually (after the existing
"Compatibility" bullet in `### Constraints`) is staged in:

```
CLAUDE_md_addition_2k1b_glue.md
```

at the repo root — copy that bullet into
`CLAUDE.md` by hand, then delete the staging file. If working from a different machine, that
staging file may not exist there yet (it's an untracked file created alongside this handoff) —
the full bullet text is duplicated below for convenience:

> - **GLUE has no 2k1b clade (structural, not a bug)**: HCV-GLUE's genotyping is phylogenetic
>   clade placement (`genotypeCladeCategoryResult`/`subtypeCladeCategoryResult` in the JSON,
>   parsed by `bin/GLUE_json_parser.R`) against a tree with no CRF_02k/1b category. A genuine
>   2k/1b recombinant is therefore reported by GLUE as genotype **1** (or **2**), never `2k1b` —
>   this is expected, not a parser defect. **Any code that compares GLUE's genotype/subtype call
>   against the mapping- or de novo-derived genotype (concordance checks, discordance gates,
>   `is_valid_minor()`-style pairing rules, rescue logic, future strain-identity checks) MUST
>   special-case `map_gt == "2k1b" && glue_gt %in% c("1","2")` as non-conflicting.** Skipping this
>   exception silently fails the major-gate for every true 2k/1b sample
>   (`overall_sample_call = "indeterminate"`, empty `Genotype` in MultiQC) — this exact
>   regression shipped once already (see `hcvtyper_handoff_2k1b_glue_concordance_fix.md` for the
>   incident and fix). The 2k1b-aware genotype string itself is single-sourced in
>   `bin/genotype_utils.R::genotype_from_subtype()` — reuse it rather than re-deriving
>   genotype-from-subtype logic ad hoc.

---

## 6. Context notes

- This is unrelated to the earlier GLUE-parser regression (`hcvtyper_handoff_glue_parser_fix.md`,
  2026-06-15, cand1/cand2 filename migration) — that one caused GLUE columns to be **entirely
  empty**; here GLUE is populated and correct (`GLUE_genotype=1, GLUE_subtype=1b` is the right
  GLUE answer for a 2k/1b sample), the bug is in how a *downstream* gate interprets that correct
  GLUE answer.
- The `hcvtyper_v3.0_vs_v1.2.0_results_comparison.md` validation doc (pre-dates the concordance
  gate, written before 2026-06-16) already flagged "2k/1b recombinant exception preserved
  (`ERR1810443`)" as a checked-off item — but that check only exercised the co-infection-pairing
  exception (`is_valid_minor()`), not this later concordance gate. Worth adding a 2k1b
  **monoinfection** case (like the three TEST samples here) to whatever validation cohort is used
  for the next full re-run, since the pairing exception and the concordance exception are
  independent code paths that can silently diverge again in the future.
