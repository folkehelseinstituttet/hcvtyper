# Handoff — `evidence_state == "refuted"` is structurally unreachable

**Date:** 2026-08-05
**Branch:** `dev` (manifest version 1.3.0); line numbers as of `1acf637`
**Reported by:** manuscript review of the Methods description of the four-state evidence scheme
**Status:** analysis only — no code changed. Needs confirmation against the simulated and Thomson datasets.

---

## Summary

A candidate's `evidence_state` can never be `refuted`, in either supported
`--denovo_match_level` mode, because **contig attribution and contig contradiction
are keyed on the same field**.

A contig is attached to a candidate only when the genotype (or subtype) of its best
BLAST hit *equals* that of the candidate's reference. Refutation requires the
attached contig's genotype to *differ* from that of the candidate's reference. The
two conditions are mutually exclusive by construction, so the refuted branch cannot
be taken.

The practical consequence is not that a wrong call is produced — it is that a
designed safety check silently does nothing. A candidate whose mapping genotype is
contradicted by the sample's only substantial contig is currently reported as
`weak` (no attributed assembly) rather than `refuted`, and the contradicting contig
is invisible to that candidate.

---

## Code path

**1. Attribution is keyed on genotype equality.**
`bin/assembly_support_join.R`

```r
# candidate side (L90-92)
.match_key = as.character(
  if (match_level == "subtype") candidate_subtype else candidate_genotype
)

# support side (L109-112)
.match_key = as.character(
  if (match_level == "subtype") subtype else genotype_from_subtype(subtype)
)

# join (L129-131)
left_join(support_collapsed, by = c("sampleName", ".match_key"))
```

`assembly_support_subtype` is then the winning contig's subtype *within that key
group* (L120). At `match_level = "genotype"` (the default) this guarantees
`genotype_from_subtype(assembly_support_subtype) == candidate_genotype`. At
`match_level = "subtype"` it guarantees the stronger
`assembly_support_subtype == candidate_subtype`. Either way the genotypes are equal.

**2. Contradiction is keyed on genotype inequality of the same column.**
`bin/classify_roles.R:80-95`

```r
own_denovo_conflict <- function(map_genotype, assembly_support, assembly_support_subtype) {
  ...
  denovo_gt <- genotype_from_subtype(asup_subtype)     # L86-92
  conflict  <- has_denovo & !is.na(denovo_gt) & !is.na(map_gt) & map_gt != denovo_gt   # L93
}
```

`map_gt` is `candidate_genotype`; `denovo_gt` is derived from
`assembly_support_subtype`, which step 1 has just guaranteed to be in the same
genotype. `conflict` is therefore always `FALSE`.

**3. The refuted band depends on that predicate.**
`bin/classify_roles.R:577-586`

```r
evidence_state_col <- ifelse(
  !asup_exists, "weak",
  ifelse(denovo_contradicts & quality_fails_state, "refuted",   # unreachable
    ifelse(asup_score >= evidence_hi_cut, "confirmed",
      ifelse(asup_score >= evidence_lo_cut, "probable", "weak"))))
```

A candidate with no attributed contig short-circuits to `weak` at the first
`ifelse`; a candidate with one can never satisfy `denovo_contradicts`.

---

## Empirical confirmation

All `summary/candidates.csv` files under
`/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCV_paper_revisjon/`
(7 result directories, 318 candidate rows):

| Check | Result |
|---|---|
| `evidence_state` distribution | 266 `confirmed`, 31 `weak`, 21 `probable`, **0 `refuted`** |
| Rows where `genotype(assembly_support_subtype) != candidate_genotype` | **0 of 318** |
| `assembly_support` | 293 `supported`, 25 `none` |
| `concordance_status` | 247 `confirmed`, 71 `unconfirmed`, **0 `discordant`** |
| `concordance_reason` | 247 `all_legs_concordant`, 46 `two_legs_denovo_only`, 25 `no_corroborating_legs` |

293 candidates had an attributed contig and not one contradicted its candidate —
consistent with the structural argument rather than with a rare event that merely
failed to occur.

---

## Relationship to the existing CR-02 note

`bin/summarize.R:1682-1693` already records that `any_refuted_denovo` is
"effectively dead on real data", attributing this to the `discordant_identity` hard
gate in `classify_one_sample()` pre-empting the refuted path, and states that
`any_discordant_identity` "is the trigger that ACTUALLY fires in production for a
genuine own-assembly-vs-mapping contradiction".

**That explanation is incomplete, and the second half does not hold on this data.**
`concordance_status` is never `discordant` in any of the 318 rows, so
`discordant_identity` does not fire either. The reason is the same one described
above: `apply_concordance()`'s de novo leg reads the same joined
`assembly_support_subtype` column via the same shared `own_denovo_conflict()`
predicate (`classify_roles.R:126`), so its de novo leg is dead for the same
structural reason.

So the correct statement is not "refuted is pre-empted by discordant_identity" but
"**both** de novo contradiction paths are dead, because the column they test cannot
contradict." The pre-emption ordering is real but moot.

**Not affected:** `apply_concordance()`'s **GLUE** leg is live. It compares
`candidate_glue_genotype` against `candidate_genotype`, neither of which is
constrained by the assembly join. It was populated on 247 of 318 rows here and
simply happened to agree every time — that leg can still fire and should not be
touched.

---

## Where the genuine contradiction case *is* handled

The real scenario — mapping selects genotype X while the sample's assembly says
genotype Y — is handled upstream, and correctly, by `RESCUE_EVALUATION`
(`bin/rescue_evaluation.R`). That module filters the support table by **sample
only** (`filter(sample == sample_id)`, L207) and then explicitly looks for the best
contig of a *different* subtype (L246-249). Because it is not keyed on genotype
equality, it can see the contradiction that `classify_roles()` cannot.

This suggests `refuted` may be a vestige of the design as it stood before rescue and
nomination existed: at that time `classify_roles()` was the only place a
mapping-vs-assembly contradiction could be caught. It is worth deciding whether the
state is still needed at all, rather than assuming it must be repaired.

---

## Downstream code that is consequently dead

- `classify_roles.R:669-671` — the `refuted` → `background` / `refuted_denovo` role branch
- `classify_roles.R:823, 888` — `refuted_denovo` role-reason and review-note strings
- `classify_roles.R:871` — `is_conflict` via `refuted_denovo`
- `summarize.R:1698, 2027` — `any_refuted_denovo` rollup and its `call_confidence = provisional` trigger
- `denovo_layer.R:65` — `minor_denovo_status == "refuted"` (legacy path; check whether still reachable via `denovo_confirm.R:78`, which computes its own separate `"refuted"` and is **not** the same predicate — worth confirming these two have not been conflated)

---

## Suggested tests

The dev team has the simulated and Thomson datasets; these should discriminate
cleanly.

1. **Construct the case the state was designed for.** Take a sample where mapping
   selects genotype X and the only substantial contig is genotype Y, with the Y
   contig failing the substantiality floors (`--denovo_min_contig_length`,
   `--denovo_min_kmer_cov`, `--denovo_min_blast_identity`). Expected under the
   current code: candidate X gets `assembly_support = "none"` → `weak`, and the Y
   contig is never associated with X. Confirm this is what happens.
2. **Check the unit tests are not passing on synthetic input that the join cannot
   produce.** `bin/tests/test_classify_roles.R` exercises the refuted band by
   constructing frames directly. If those frames carry an
   `assembly_support_subtype` whose genotype differs from `candidate_genotype`,
   they encode a state the production join can never emit, and the tests are
   green on an impossible input. Test20 (cited in the CR-02 note as the call-order
   proof) is worth inspecting first.
3. **Re-run the Thomson co-infection set** and confirm `evidence_state` never takes
   the value `refuted` and `concordance_status` never takes `discordant`, matching
   the 318-row result above.

---

## Options (not a recommendation — needs a design decision)

1. **Attribute contigs more broadly than the contradiction test.** Give
   `classify_roles()` access to the sample's full contig set, as
   `RESCUE_EVALUATION` already has, so a contradicting contig can be seen. This
   restores the intended behaviour but re-introduces a sample-scoped input to a
   function whose current design deliberately keeps each candidate's evidence to
   itself (D-11 / EVID-02).
2. **Delegate refutation entirely to rescue/nomination** and remove `refuted`, the
   `refuted_denovo` role reason, and `any_refuted_denovo`. Simplest, and arguably
   honest about where the check now lives. Requires a Summary.csv schema note since
   `Major_evidence_state` / `Minor_evidence_state` are documented as four-valued.
3. **Leave as-is and document.** Acceptable only if the case is genuinely covered
   by rescue in all configurations — including when rescue is disabled or its
   floors are not met.

---

## Manuscript impact

The revision under review (`Access_Microbiology/Revision_1/revised_manuscript_complete.md`)
currently:

- describes the four-state scheme in Methods, "*De novo* minor-strain confirmation";
- lists "a minor was refuted" as a `call_confidence = provisional` trigger;
- states in Results that within the real-world single-infection set "every retained
  minor was supported by substantial assembly evidence and none was refuted" —
  which reads as an empirical finding but is structural.

None of these change a reported result, since no call depends on the unreachable
branch. But the description should be reconciled with whichever option above is
chosen before submission (deadline 2026-08-11).
