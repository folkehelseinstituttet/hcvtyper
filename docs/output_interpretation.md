# Interpreting HCVTyper results

This guide explains how to read the `Summary.csv` table and the **Results summary**
section of the MultiQC report, and — importantly — how to _confirm a call_ against the
underlying evidence the pipeline publishes.

HCV NGS data is full of edge cases (cross-mapping, fragmented de-novo assemblies,
low-yield co-infections, recombinants). The report is therefore built around a simple
principle: **the pipeline states what the data indicate, and separately states how much
to trust it.** A concrete call is still produced for every sample — but a marginal call
is labelled as such so it reads differently from a clean one at a glance.

---

## The two-axis model: Final call + Confidence

Read these two columns together first. Everything else is supporting detail.

| Column                                 | Question it answers           | Values                                                                                                  |
| -------------------------------------- | ----------------------------- | ------------------------------------------------------------------------------------------------------- |
| `overall_sample_call` (**Final call**) | _What do the data indicate?_  | `monoinfection`, `co-infection`, `co-infection (indeterminate dominance)`, `indeterminate`, `untypable` |
| `call_confidence` (**Confidence**)     | _How much should I trust it?_ | `high`, `provisional`, `review`, `indeterminate`                                                        |

`call_confidence` is not new biology — it is an explicit roll-up of signals that are
otherwise scattered across `review_flag`, `gate_flag`, `rescue_effect`, the
subtype-match columns and the per-candidate roles. Its whole purpose is so you do not
have to reconstruct confidence in your head.

### Final call vocabulary

| Value                                    | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `monoinfection`                          | One dominant strain; any other candidate was demoted to background.                                                                                                                                                                                                                                                                                                                                                                                    |
| `co-infection`                           | Two (or more) strains cleared the abundance floor **and** have de-novo/assembly support.                                                                                                                                                                                                                                                                                                                                                               |
| `co-infection (indeterminate dominance)` | Two strains are present but read-count and k-mer-coverage rankings disagree on which is dominant. Reported as present; do not rely on the major/minor ordering.                                                                                                                                                                                                                                                                                        |
| `indeterminate`                          | No candidate is eligible to be called dominant: coverage exists, but no candidate combines measurable coverage with a mapping identity that agrees with GLUE. Note this is **not** a `--minRead`/`--minCov` test — that floor has been informational since D-07/D-09 and gates nothing here.                                                                                                                                                           |
| `untypable`                              | No usable coverage on any candidate, on **either** mapping — neither the first-pass panel-wide breadth (`candidate_cov`) nor the targeted post-competition breadth (`cand_cov_breadth`) shows anything above zero. A candidate measured on only one of the two axes is still typable; up to and including 2.0.0 a de novo-rescued candidate, whose first-pass numbers are deliberately blank, could reach this state with excellent targeted coverage. |

### Confidence tiers

| Tier            | What it means                                                                                                                                                                                                              | Typical action                                                                     |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `high`          | No conflicting signals fired; the call stands on clean evidence. `Review note` is empty.                                                                                                                                   | Accept.                                                                            |
| `provisional`   | A **soft** caveat: identity not corroborated, a co-infection kept without corroboration, a minor-subtype conflict, or a minor-slot rescue. Usable, but note the caveat.                                                    | Skim the `Review note` and evidence columns.                                       |
| `review`        | A **hard** conflict that should block automatic reporting: the major failed mapping QC, de-novo and mapping disagree on the major subtype, the de-novo rescue overrode the **Major** reference, or dominance is ambiguous. | Manually adjudicate before reporting — see `Review note` and the underlying files. |
| `indeterminate` | No actionable call (indeterminate/untypable).                                                                                                                                                                              | Inspect raw QC; the sample may need re-sequencing.                                 |

The two are kept consistent by construction: a sample with a non-empty `Review note`
(`review_flag`) is **never** `high` — it is demoted to at least `provisional`. So a
`high` call always has an empty Review note. (The converse does not hold: a few
`provisional` samples carry no prose note — e.g. a minor-slot rescue — because
`provisional` is a softer bucket than a full review sentence.)

---

## `review_flag` — the human-readable "why"

`review_flag` (**Review note** in MultiQC, highlighted orange) is empty for clean
samples and otherwise contains one or more full sentences joined by `|`, each ending
in a prompt to review. The triggers:

1. **Co-infection, but major/minor assignment uncertain** — de novo and mapping disagree on which strain is dominant.
2. **Major subtype conflict (monoinfection)** — de-novo assembly and mapping disagree on the major subtype; possible reference mismatch or divergent strain.
3. **Different-genotype de-novo contig under a monoinfection** — possible missed co-infection or contamination.
4. **Genotype call provisional** — dominant identity not corroborated by GLUE or de novo (mapping evidence only).
5. **Co-infection kept without de-novo corroboration** — de novo inconclusive for both strains; the minor may be a genuine low-yield co-infection.
6. **Indeterminate dominance** — read-count and k-mer rankings disagree.
7. **No candidate eligible to be called dominant** — no candidate combines measurable coverage with a mapping identity that agrees with GLUE; call indeterminate.
8. **Major strain failed mapping quality thresholds.**
9. **De-novo rescue overrode the Major reference** — the primary call was reassigned automatically; confirm against `rescue_audit.csv` before reporting.
10. **Co-infection minor looks like background reads** — the minor holds a small share of the reads and/or its contigs are far less PCR-duplicated than the dominant genotype's (see below). The call is unchanged; check `contamination_check/`.

### Trigger 10: a co-infection minor that looks like background

High-depth probe-capture runs make it easy for reads bled in from another sample — index hopping,
demultiplexing bleed, carry-over — to clear absolute floors like `--minRead` and to assemble into
long contigs that corroborate a minor. Two signals, computed within the sample alone, mark such a
minor. They never change the call; they add a sentence such as

> Co-infection candidate 2 (2b_D10988) holds 1.23% of the sample's HCV reads (below 2.00%, the range
> of cross-sample background), and its contigs are 1.80% PCR duplicates against 90.9% for the dominant
> genotype, so its reads probably did not go through this library's PCR — consistent with background reads from
> another sample (index hopping or carry-over); check contamination_check/ before reporting.

**Read fraction** (`Minor_read_fraction`, `--minor_min_read_fraction`). Bleed is roughly a constant
fraction of a run's reads, so a minor's _share_ of the sample's reads is comparable across depths in a
way a read count is not. A low share is not proof: in the IVT dilution series the genuine 1:50 minors
hold 1.95 % and 2.7 % of reads, overlapping the 0.17–1.7 % of the background minors in the Thomson
MiSeq run. So the sentence always says what the duplicate test found — a low-share minor that went
through the library's PCR reads as "a genuine low-level strain, or contamination before library
preparation", not as background.

**Duplicate discordance** (`Minor_dup_rate`, `Major_dup_rate`, `Minor_dup_discordant`). Each sample's
reads are mapped back to its own contigs and PCR duplicates are marked. Within one library every
molecule is amplified and sequenced alike, so a genuine minor carries the library's duplicate rate.
Reads that arrived from another library after its PCR are a sparse random draw from a large pool and
are almost all unique. A minor at a small fraction of the dominant's duplicate rate therefore did not
go through this sample's library preparation. This is the stronger of the two signals: in the IVT
dilution series every genuine minor (1:1 to 1:50) was duplicated like its dominant (96.9–98.2 % vs
97.8–98.7 %), while every background minor in the Thomson run was far below its library (1.8–8.6 % vs
91–96 %). Its limits:

- **One-directional.** Contamination that entered _before_ library PCR (sample handling, extraction)
  is amplified with the sample and looks concordant. `FALSE` means "not shown to be foreign", not
  "genuine".
- **Needs a saturated library.** If the dominant genotype's contigs are barely duplicated
  (`--dup_min_library_rate`), duplicates say nothing either way and the result is `NA`. Capture
  libraries are usually saturated; low-input metagenomic libraries often are not.
- **Needs reads.** Fewer than `--dup_min_fragments` read pairs on the minor's contigs gives `NA`.
- **Genotype level.** Contigs are pooled by the genotype of their best BLAST hit, so a minor of the
  dominant's own genotype (1a/1b) cannot be tested and gives `NA`.

The contamination check (`contamination_check/`) is the natural next look: it shows whether the
minor's contigs are near-identical to a much higher-coverage sample in the same run.

### Every sentence names its candidate and its numbers

A review sentence is only useful if it can be adjudicated without re-running anything. Each one
therefore names **which** candidate it refers to (rank/slot, reference, or subtype) and carries the
**measured value** next to the floor or expectation it missed — for example _"candidate 2 contig matched
2c but mapping says 3a"_, or _"identity 89.0 below 90 floor"_.

The two contig-based triggers (2 and 3) go further and render all four measured metrics of the
conflicting contig, with an interpretation clause that switches on the **aligned fraction**:

> …different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0;
> only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly
> supported — the contig may be largely non-HCV, chimeric, or too divergent to type.

A contig that aligns over its full length keeps the straightforward co-infection wording; genuine minors
align over 99–100% of their contigs. Nothing is suppressed by this — the analyst decides, with the
numbers present. For the major-conflict context the interpretation is deliberately different: a 69 bp
anchor cannot support a subtype call, so it cannot support a _disagreement_ with one either. The honest
reading is that the conflict may be phantom, not that the contig is junk.

### Trigger 3 has a substantiality floor

The "different-genotype contig under a monoinfection" sentence is gated on two conditions:

- the contig must reach `--review_min_offgenotype_contig_length` (default **1000 bp**), and
- the pair must not be a 2k/1b recombinant against a genotype 1 or 2 major.

Without the length floor this trigger fired on 51% of a 140-sample cohort — median triggering contig
606 bp, shortest 142 bp — and accounted for 92% of all `provisional` calls. A flag that fires on half a
cohort is a flag reviewers learn to ignore. The floor is deliberately **separate from and higher than**
`--denovo_min_contig_length`: that floor confirms a minor that mapping already supports, whereas here
the contig is the _only_ evidence. Length is the only leg used — genuine low-yield minors sit _below_
`--denovo_min_kmer_cov`, so a k-mer leg would suppress exactly the samples most worth reviewing.

---

## Per-candidate roles and `role_reason` (analyst glossary)

Each candidate is classified into a `role` with a coded `role_reason`. These appear in
`Summary.csv` (`Major_role_reason` / `Minor_role_reason`) and in the per-sample
`*.candidates.csv` / `*.rescued.candidates.csv`. Plain-language mapping:

| `role`          | `role_reason`                      | Plain language                                                                               |
| --------------- | ---------------------------------- | -------------------------------------------------------------------------------------------- |
| `dominant`      | `dominant`                         | The primary strain (highest dominance score among gate-passing candidates).                  |
| `co-infection`  | `corroborated`                     | A second strain confirmed by de-novo/assembly support.                                       |
| `background`    | `same_genotype_as_dominant`        | Demoted: same genotype as the dominant strain (not counted as a distinct co-infection).      |
| `background`    | `recombinant_2k1b`                 | Demoted by the 2k/1b recombinant exception rule.                                             |
| `background`    | `discordant_identity`              | The candidate's mapping identity and its GLUE identity disagree.                             |
| `indeterminate` | `indeterminate_dominance_conflict` | Part of an indeterminate-dominance pair (reads favour one strain, k-mer coverage the other). |

Background candidates are **surfaced explicitly with their reason**, never silently
dropped — so a demoted strain is always visible in the candidate files.

---

## De-novo rescue columns

The rescue step compares first-mapping against de-novo assembly + contig BLAST and can
reassign a candidate's reference. Because that is an automated override, it is traceable
at three levels of detail:

| Column / file               | Granularity                                        | Use                                                                                                                                                                                                                            |
| --------------------------- | -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `rescue_flag`               | boolean                                            | Did any rescue survive into the call?                                                                                                                                                                                          |
| `rescue_effect`             | `none` / `minor_ref_changed` / `major_ref_changed` | **Where** a surviving rescue landed. `major_ref_changed` always forces `review` confidence and a review note — it reassigned the dominant strain.                                                                              |
| `{sample}.rescue_audit.csv` | full ledger                                        | Authoritative from/to + trigger for **every** rescue, nomination and block decision — **including** rescues that fired then were dropped by genotype-collapse or the candidate cap (which leave no trace in `candidates.csv`). |

When `rescue_effect` is not `none`, open `rescue_audit.csv` to see exactly which
reference was swapped for which, and why.

---

## `Major_evidence` / `Minor_evidence` — the basis for each strain

These two columns condense the basis for each reported strain into one string, e.g.:

```
1a_HQ850279 | 97,242 reads (nodup) | 91% breadth@10x | de novo 1a | 99.2% consensus id | ref from first-mapping
```

Tokens are only shown when the underlying value exists. The final token states whether
the reference came straight from first-mapping or was **reassigned by the de-novo
rescue** — the quickest way to spot a rescued call. `Minor_evidence` is empty for
monoinfections.

---

## How to confirm a call

For any sample tagged `provisional` or `review`, the evidence to adjudicate it is
published:

1. **Read the `Review note`** — it names the specific conflict.
2. **Check `Major_evidence` / `Minor_evidence`** — reads, breadth, de-novo subtype, consensus identity, and whether the reference was rescued.
3. **If `rescue_effect` ≠ `none`** — open `blastparse/{sample}.rescue_audit.csv` for the from/to and trigger.
4. **Coverage plots** — `plots/` per-strain depth/evenness across the genome (a spiky high-read candidate with poor breadth is a cross-mapping artefact).
5. **Candidate table** — `parsefirstmapping/{sample}.candidates.csv` and `blastparse/{sample}.rescued.candidates.csv` list every candidate with its role and reason, including demoted background strains.
6. **GLUE report** — the genotype/subtype and resistance call on the strain consensus.

A `high`-confidence call needs none of this; the point of the confidence axis is to tell
you which samples do.

> **`denovo_minor_ref`, `denovo_minor_contig` and `denovo_minor_contig_length` describe one
> contig.** They previously did not: the reference was selected from the per-contig table
> (restricted to contigs whose _own_ top hit is off-genotype) while the contig name was
> re-derived from the full hit table, so a conserved 5′UTR/core region on an unrelated contig
> could win on bitscore and send an analyst to the wrong sequence. If you are working from
> results produced before this fix, verify the contig name against the length before trusting it.

---

## Tuning the thresholds

Every floor referenced above is a pipeline parameter, documented with its default and rationale
in the [README](../README.md#optional-parameters):

| What you want to change                               | Parameters                                                                                                   |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| How many candidate references are selected and mapped | `--n_candidates`                                                                                             |
| The abundance floor a candidate must clear            | `--minRead`, `--minCov`                                                                                      |
| When a contig counts as corroborating evidence        | `--denovo_min_contig_length`, `--denovo_min_kmer_cov`, `--denovo_min_blast_identity`, `--denovo_match_level` |
| When a contig may reassign a candidate's reference    | the `--rescue_*` family                                                                                      |
| What dominance means (reads vs. evenness)             | the `--score_weight_*` family                                                                                |
| How noisy the off-genotype-contig review flag is      | `--review_min_offgenotype_contig_length`                                                                     |

Raising a floor makes the pipeline quieter and riskier; lowering it makes it noisier and safer.
The defaults are calibrated on a 140-sample routine cohort and are the recommended starting
point — with an expert reviewing the output, losing a real co-infection is the worse error.
