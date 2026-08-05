# Interpreting HCVTyper results

This guide explains how to read the `Summary.csv` table and the **Results summary**
section of the MultiQC report, and — importantly — how to *confirm a call* against the
underlying evidence the pipeline publishes.

HCV NGS data is full of edge cases (cross-mapping, fragmented de-novo assemblies,
low-yield co-infections, recombinants). The report is therefore built around a simple
principle: **the pipeline states what the data indicate, and separately states how much
to trust it.** A concrete call is still produced for every sample — but a marginal call
is labelled as such so it reads differently from a clean one at a glance.

---

## The two-axis model: Final call + Confidence

Read these two columns together first. Everything else is supporting detail.

| Column | Question it answers | Values |
|---|---|---|
| `overall_sample_call` (**Final call**) | *What do the data indicate?* | `monoinfection`, `co-infection`, `co-infection (indeterminate dominance)`, `indeterminate`, `untypable` |
| `call_confidence` (**Confidence**) | *How much should I trust it?* | `high`, `provisional`, `review`, `indeterminate` |

`call_confidence` is not new biology — it is an explicit roll-up of signals that are
otherwise scattered across `review_flag`, `gate_flag`, `rescue_effect`, the
subtype-match columns and the per-candidate roles. Its whole purpose is so you do not
have to reconstruct confidence in your head.

### Final call vocabulary

| Value | Meaning |
|---|---|
| `monoinfection` | One dominant strain; any other candidate was demoted to background. |
| `co-infection` | Two (or more) strains cleared the abundance floor **and** have de-novo/assembly support. |
| `co-infection (indeterminate dominance)` | Two strains are present but read-count and k-mer-coverage rankings disagree on which is dominant. Reported as present; do not rely on the major/minor ordering. |
| `indeterminate` | No candidate passed the major-gate (first-mapping quality thresholds). |
| `untypable` | No usable coverage on any candidate. |

### Confidence tiers

| Tier | What it means | Typical action |
|---|---|---|
| `high` | No conflicting signals fired; the call stands on clean evidence. `Review note` is empty. | Accept. |
| `provisional` | A **soft** caveat: identity not corroborated, a co-infection kept without corroboration, a minor-subtype conflict, or a minor-slot rescue. Usable, but note the caveat. | Skim the `Review note` and evidence columns. |
| `review` | A **hard** conflict that should block automatic reporting: the major failed mapping QC, de-novo and mapping disagree on the major subtype, the de-novo rescue overrode the **Major** reference, or dominance is ambiguous. | Manually adjudicate before reporting — see `Review note` and the underlying files. |
| `indeterminate` | No actionable call (indeterminate/untypable). | Inspect raw QC; the sample may need re-sequencing. |

The two are kept consistent by construction: a sample with a non-empty `Review note`
(`review_flag`) is **never** `high` — it is demoted to at least `provisional`. So a
`high` call always has an empty Review note. (The converse does not hold: a few
`provisional` samples carry no prose note — e.g. a minor-slot rescue — because
`provisional` is a softer bucket than a full review sentence.)

---

## `review_flag` — the human-readable "why"

`review_flag` (**Review note** in MultiQC, highlighted orange) is empty for clean
samples and otherwise contains one or more full sentences joined by ` | `, each ending
in a prompt to review. The triggers:

1. **Co-infection, but major/minor assignment uncertain** — de novo and mapping disagree on which strain is dominant.
2. **Major subtype conflict (monoinfection)** — de-novo assembly and mapping disagree on the major subtype; possible reference mismatch or divergent strain.
3. **Different-genotype de-novo contig under a monoinfection** — possible missed co-infection or contamination.
4. **Genotype call provisional** — dominant identity not corroborated by GLUE or de novo (mapping evidence only).
5. **Co-infection kept without de-novo corroboration** — de novo inconclusive for both strains; the minor may be a genuine low-yield co-infection.
6. **Indeterminate dominance** — read-count and k-mer rankings disagree.
7. **No candidate passed the major-gate** — call indeterminate.
8. **Major strain failed mapping quality thresholds.**
10. **De-novo rescue overrode the Major reference** — the primary call was reassigned automatically; confirm against `rescue_audit.csv` before reporting.

### Every sentence names its candidate and its numbers

A review sentence is only useful if it can be adjudicated without re-running anything. Each one
therefore names **which** candidate it refers to (rank/slot, reference, or subtype) and carries the
**measured value** next to the floor or expectation it missed — for example *"candidate 2 contig matched
2c but mapping says 3a"*, or *"identity 89.0 below 90 floor"*.

The two contig-based triggers (2 and 3) go further and render all four measured metrics of the
conflicting contig, with an interpretation clause that switches on the **aligned fraction**:

> …different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0;
> only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly
> supported — the contig may be largely non-HCV, chimeric, or too divergent to type.

A contig that aligns over its full length keeps the straightforward co-infection wording; genuine minors
align over 99–100% of their contigs. Nothing is suppressed by this — the analyst decides, with the
numbers present. For the major-conflict context the interpretation is deliberately different: a 69 bp
anchor cannot support a subtype call, so it cannot support a *disagreement* with one either. The honest
reading is that the conflict may be phantom, not that the contig is junk.

### Trigger 3 has a substantiality floor

The "different-genotype contig under a monoinfection" sentence is gated on two conditions:

- the contig must reach `--review_min_offgenotype_contig_length` (default **1000 bp**), and
- the pair must not be a 2k/1b recombinant against a genotype 1 or 2 major.

Without the length floor this trigger fired on 51% of a 140-sample cohort — median triggering contig
606 bp, shortest 142 bp — and accounted for 92% of all `provisional` calls. A flag that fires on half a
cohort is a flag reviewers learn to ignore. The floor is deliberately **separate from and higher than**
`--denovo_min_contig_length`: that floor confirms a minor that mapping already supports, whereas here
the contig is the *only* evidence. Length is the only leg used — genuine low-yield minors sit *below*
`--denovo_min_kmer_cov`, so a k-mer leg would suppress exactly the samples most worth reviewing.

---

## Per-candidate roles and `role_reason` (analyst glossary)

Each candidate is classified into a `role` with a coded `role_reason`. These appear in
`Summary.csv` (`Major_role_reason` / `Minor_role_reason`) and in the per-sample
`*.candidates.csv` / `*.rescued.candidates.csv`. Plain-language mapping:

| `role` | `role_reason` | Plain language |
|---|---|---|
| `dominant` | `dominant` | The primary strain (highest dominance score among gate-passing candidates). |
| `co-infection` | `corroborated` | A second strain confirmed by de-novo/assembly support. |
| `background` | `same_genotype_as_dominant` | Demoted: same genotype as the dominant strain (not counted as a distinct co-infection). |
| `background` | `recombinant_2k1b` | Demoted by the 2k/1b recombinant exception rule. |
| `background` | `discordant_identity` | The candidate's mapping identity and its GLUE identity disagree. |
| `indeterminate` | `indeterminate_dominance_conflict` | Part of an indeterminate-dominance pair (reads favour one strain, k-mer coverage the other). |

Background candidates are **surfaced explicitly with their reason**, never silently
dropped — so a demoted strain is always visible in the candidate files.

---

## De-novo rescue columns

The rescue step compares first-mapping against de-novo assembly + contig BLAST and can
reassign a candidate's reference. Because that is an automated override, it is traceable
at three levels of detail:

| Column / file | Granularity | Use |
|---|---|---|
| `rescue_flag` | boolean | Did any rescue survive into the call? |
| `rescue_effect` | `none` / `minor_ref_changed` / `major_ref_changed` | **Where** a surviving rescue landed. `major_ref_changed` always forces `review` confidence and a review note — it reassigned the dominant strain. |
| `{sample}.rescue_audit.csv` | full ledger | Authoritative from/to + trigger for **every** rescue, nomination and block decision — **including** rescues that fired then were dropped by genotype-collapse or the candidate cap (which leave no trace in `candidates.csv`). |

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
> (restricted to contigs whose *own* top hit is off-genotype) while the contig name was
> re-derived from the full hit table, so a conserved 5′UTR/core region on an unrelated contig
> could win on bitscore and send an analyst to the wrong sequence. If you are working from
> results produced before this fix, verify the contig name against the length before trusting it.

---

## Tuning the thresholds

Every floor referenced above is a pipeline parameter, documented with its default and rationale
in the [README](../README.md#optional-parameters):

| What you want to change | Parameters |
|---|---|
| How many candidate references are selected and mapped | `--n_candidates` |
| The abundance floor a candidate must clear | `--minRead`, `--minCov` |
| When a contig counts as corroborating evidence | `--denovo_min_contig_length`, `--denovo_min_kmer_cov`, `--denovo_min_blast_identity`, `--denovo_match_level` |
| When a contig may reassign a candidate's reference | the `--rescue_*` family |
| What dominance means (reads vs. evenness) | the `--score_weight_*` family |
| How noisy the off-genotype-contig review flag is | `--review_min_offgenotype_contig_length` |

Raising a floor makes the pipeline quieter and riskier; lowering it makes it noisier and safer.
The defaults are calibrated on a 140-sample routine cohort and are the recommended starting
point — with an expert reviewing the output, losing a real co-infection is the worse error.
