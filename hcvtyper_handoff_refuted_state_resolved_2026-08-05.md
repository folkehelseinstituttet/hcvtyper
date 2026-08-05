# Handoff — `refuted` removed (option 2), with two corrections to the original analysis

**Date:** 2026-08-05
**Follows:** `hcvtyper_handoff_refuted_state_unreachable_2026-08-05.md` (`f576b44`)
**Shipped:** `5c8874d` (code + tests), `103253a` (docs)
**Status:** implemented and verified. **Manuscript wording still outstanding — see §6.**

---

## 1. What was decided and done

Option 2: refutation is delegated entirely to `RESCUE_EVALUATION` and the `refuted`
state is removed. `evidence_state` is now three-valued: `confirmed` / `probable` / `weak`.

The reporting handoff's structural argument was **confirmed in full** — attribution and
contradiction are keyed on the same field, so the band could never be taken.

---

## 2. Correction 1 — `discordant_identity` is NOT dead. It was nearly deleted.

The original handoff states:

> The same predicate backs `apply_concordance()`'s de novo leg, so
> `discordant_identity` is dead for the same reason.

**This is wrong, and following it literally would have removed a live safety check.**

`apply_concordance()` sets `discordant` on `glue_conflict || denovo_conflict`. Only the
**de novo** leg was dead. The **GLUE** leg compares `candidate_genotype` against
`candidate_glue_genotype`, neither of which is constrained by the assembly join, and it
is live — as the same handoff correctly says two paragraphs later when discussing the
GLUE leg in isolation.

The two statements cannot both hold. The empirical observation that backed the claim —
0 `discordant` in 318 rows — is explained by GLUE happening to agree on every row, not
by unreachability. `concordance_status == "discordant"`, the `discordant_identity` role
reason and its hard gate are therefore **retained**.

What *was* removed from that path: the de novo conflict leg itself, and the two reason
strings that depended on it (`discordant_all_legs`, `discordant_mapping_vs_denovo`),
which were unreachable. `discordant_mapping_vs_glue` is now the only reason string, and
`discordant` has exactly one source.

**Lesson for the next reader:** "0 occurrences in the data" and "unreachable by
construction" are different claims. This document's §3 argument establishes the second
for `refuted`; it does not establish it for `discordant`, and the row counts alone never
could.

---

## 3. Correction 2 — four unit tests were green on impossible input, not three

The original handoff's suggested test 2 asked whether the unit tests pass on frames the
join cannot produce. They do, and there was one more than it flagged:

| test | fixture | what it asserted |
|---|---|---|
| `test_classify_roles.R` Test17(a) | `assembly_support_subtype = "3a"` on a genotype-1 candidate | contradiction → `refuted` |
| Test19d | same shape | `refuted` → `background` / `refuted_denovo` |
| Test20 | same shape | the CR-02 call-order claim |
| **Test11 `S3`** | `mk_conc("S3", "1a", NA, "supported", "3a")` | mapping-vs-de-novo → `discordant` |

`mk_cand()` and `mk_conc()` both build the *joined* frame directly and never call
`assembly_support_join()`, so the invariant that makes the state impossible was never
exercised. That is why this survived, and it is why the CR-02 note in `summarize.R` —
reasoned from Test20's behaviour — asserted as production fact something that never
happened in production.

**Test20 has been rewritten to compose `join_assembly_support()`** at both
`--denovo_match_level` settings and assert the invariant directly:

- an off-genotype contig is never attributed (`assembly_support == "none"`)
- an attributed contig always shares the candidate's genotype
- the candidate lands `weak` / `no_own_assembly` — the same demotion the deleted arm
  produced, by a truthful reason
- a same-genotype contig *is* attributed and does yield `confirmed` (the join works)

If anyone later widens attribution, that test fails and the design decision gets made
deliberately rather than by accident.

`Test11 S3` is repurposed to cover de-novo-only-and-agreeing, a case the suite was
missing entirely.

---

## 4. Why it was not repaired

Worth recording, because "restore the intended behaviour" is the intuitive option and it
is the wrong one.

The band required `denovo_contradicts & quality_fails_state`, where
`quality_fails_state = !own_substantial`. It fired **only when the contradicting contig
FAILED the substantiality floors** (< 500 bp, < 2.0× k-mer, or < 90% identity). That is
the profile of assembly noise, not of a second strain.

The IVT dilution series measures this directly. Those nine mixtures contain only 1a and
2a, yet carry six off-genotype contigs:

| sample | contig | length | k-mer |
|---|---|---|---|
| IVT5 | `3a_D17763` | 320 bp | 1.22× |
| IVT2 | `3a_D17763` | 267 bp | 0.80× |
| IVT5 | `3a_X76918` | 257 bp | 0.85× |
| IVT7 | `3a_D17763` | 248 bp | 0.79× |
| IVT6 | `3a_D17763` | 226 bp | 1.02× |
| IVT4 | `7b_KX092342` | 109 bp | 652× |

Every one fails `own_substantial`. Under option 1 — widening attribution so a
contradicting contig can be seen — each would have set both legs TRUE and **refuted a
correct dominant candidate** in IVT2, IVT5, IVT6 and IVT7 on the strength of noise.

The case genuinely worth catching, a *substantial* contradicting contig, is owned by
`RESCUE_EVALUATION`: unconditional (no skip parameter), sample-scoped, and explicitly
seeking a different-subtype contig at 3,000 bp / 85% / 3,000 bp aln / 2.0× k-mer. The
two mechanisms cover opposite ends of the evidence range, and `refuted` covered the
wrong one.

A candidate whose only contig is off-genotype still gets `assembly_support = "none"` →
`weak` → `background` / `no_own_assembly`. It is demoted either way; only the label
differs.

---

## 5. Verification

| check | result |
|---|---|
| Structural argument re-derived in code | confirmed; both call sites pass the genotype-equal pair |
| Empirical, 3 further result dirs (27 samples, 40 candidate rows) | 0 refuted, 0 discordant, 0 genotype mismatches — **358 rows** with the original 318 |
| `bin/tests/run_all.sh` | **ALL PASS**, 0 failures |
| `nf-test tests/default.nf.test` | **PASSED, snapshot UNCHANGED** |
| Replay of `summarize.R` over a 10-sample run | **`Summary.csv` and `candidates.csv` BYTE-IDENTICAL** |

The last row is the load-bearing one. The pre-removal baseline contained zero `refuted`,
zero `discordant` and zero `refuted_denovo`, so if any removed path had been reachable
the output would have moved. It did not.

---

## 6. Still outstanding: the manuscript

**Not addressed here — the manuscript is not in this repo** and was not available on the
machine this work was done on. Deadline **2026-08-11**.

Three passages need reconciling, in priority order:

1. **Results — the load-bearing one.** *"every retained minor was supported by
   substantial assembly evidence and none was refuted"* reads as an empirical finding.
   It was **true by construction** and could not have come out otherwise. The other two
   items are merely stale; this one makes a claim about the data that the data never had
   a chance to falsify. If only one thing is fixed, fix this.
2. **Methods, "*De novo* minor-strain confirmation"** — describes a four-state scheme.
   It is now three-valued (`confirmed` / `probable` / `weak`).
3. **The `call_confidence = provisional` trigger list** — remove "a minor was refuted".

Suggested replacement for (1), if the underlying point is worth keeping:

> Candidates were demoted to background when no contig of their own genotype cleared the
> corroboration floors. Contradiction between a candidate's mapping genotype and a
> substantial contig of a different genotype is handled upstream by the de novo rescue
> step rather than at classification.

No reported result changes, since no call ever depended on the removed branch.

---

## 7. Also worth knowing

`denovo_confirm.R`'s `"refuted"` (CONF-02: "the major's subtype is substantially
assembled and the minor's is not") is a **different, live predicate** and was left
untouched, as were `denovo_layer.R` and the tests covering them. The original handoff was
right to ask whether the two had been conflated; in the code they had not — they are
homonyms.

Two documented values were found stale *independently* of this change and corrected in
`103253a`: `uncorroborated_kept`, documented in `README.md` (twice) and
`docs/output_interpretation.md` as a live `role_reason` though Phase 12 retired it; and a
MultiQC conditional-formatting rule matching `review_flag` on the substring `"refuted"`,
which no review-note string has ever contained. Both suggest the README's vocabulary
lists are not checked against the code — a lint asserting that every `role_reason` and
`evidence_state` named in the docs is producible would have caught all three.
