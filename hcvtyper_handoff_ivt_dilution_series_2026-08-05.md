# Handoff — the IVT dilution series: sensitivity, its reporting limit, and why the default should stay

**Date:** 2026-08-05
**Branch:** `dev` @ `1acf637`
**Purpose:** supply the manuscript with a defensible sensitivity statement, and record why the
observed assembly sensitivity is deliberately *not* being converted into a parameter change.
**Status:** measurement complete; **no code change proposed, and none wanted.**

---

## 1. Summary

The nine IVT mixtures were downloaded from ENA and run end to end. Two results, and they point in
opposite directions:

1. **De novo assembly recovers the minority strain at every ratio, 1:1 through 1:5000** — down to a
   505 bp contig at 0.94× k-mer coverage. Every predicted contig length and k-mer value from the
   260804-dnrank handoff §8 reproduced to the digit.
2. **The pipeline only *reports* the minority down to ~1:500.** Beyond that the sample is reported as
   a clean monoinfection at `high` confidence with an empty `review_flag`, even though the assembly
   contains the second strain.

The gap between those two is `review_min_offgenotype_contig_length` (1,000 bp).

**The recommendation is to leave that default alone.** §5 shows why: in this same controlled material,
the false-positive contig population overlaps the true low-abundance minority population so closely
that no threshold separates them. The IVT series is the cleanest input this pipeline will ever see,
and even here the margin is thin.

---

## 2. What was run

Nine public ENA run accessions from Thomson *et al.* (2016) — `ERR1810511`, `513`, `515`, `517`,
`519`, `521`, `523`, `525`, `527`, labelled `IVT1`–`IVT9` in ENA metadata. *In vitro* RNA transcript
mixtures of **H77 (1a, AF009606)** and **JFH-1 (2a, AB047639)**.

```bash
# 1. fetch (464 MB, a few minutes)
nextflow run nf-core/fetchngs -revision 1.12.0 -profile docker \
  --input ids_nohdr.csv --outdir results --nf_core_pipeline rnaseq

# 2. run at 1acf637, stock parameters
nextflow run hcvtyper -profile docker --input samplesheet.csv --outdir . \
  --references data/blast_db/HCVgenosubtypes_8.5.19_clean.fa \
  --kraken_focused_db data/kraken_db/db_hepacivirus
```

Defaults in force: `denovo_min_contig_length 500`, `denovo_min_kmer_cov 2.0`,
`denovo_min_blast_identity 90`, `review_min_offgenotype_contig_length 1000`.

**Two `fetchngs` input traps.** `--input` rejects a `.txt` extension outright, and it reads a CSV
**header row as an accession** (`Mixture of ids provided via --input: id`). Use a headerless `.csv`.

---

## 3. Assembly sensitivity — the minority is recovered at every ratio

| | ratio (1a:2a) | minority | minority contig | k-mer cov |
|---|---|---|---|---|
| IVT1 | 1:1 | — | 8,177 bp / 7,517 bp | ≈1,550× both |
| IVT2 | 1:5 | 1a | 4,585 bp | 373× |
| IVT3 | 5:1 | 2a | 2,578 bp | 326× |
| IVT4 | 1:50 | 1a | 1,563 bp | 111× |
| IVT5 | 50:1 | 2a | 3,075 bp | 232× |
| IVT6 | 1:500 | 1a | 1,395 bp | 16.7× |
| IVT7 | 500:1 | 2a | 349 bp | 1.4× |
| IVT8 | 1:5000 | 1a | 505 bp | 0.94× |
| IVT9 | 5000:1 | 2a | 920 bp | 1.7× |

Every value matches 260804-dnrank §8 exactly.

**Caveat on the ratio labels.** ENA metadata gives only `IVT1`–`IVT9`; it carries no ratios. The ratio
column above is inferred by matching measured contig lengths and k-mer coverages to §8's table, and is
corroborated by internal structure — the minority direction alternates cleanly from IVT2 (1a, 2a, 1a,
2a, …) and minority k-mer coverage falls monotonically within each direction (1a: 373 → 111 → 16.7 →
0.94). **Confirm the ratio↔IVT mapping against Thomson *et al.* (2016) before publication.** It is
inference from data, not an authoritative label.

---

## 4. Reporting sensitivity — three tiers, one threshold

| ratios | reported call | confidence | `review_flag` | minority visible? |
|---|---|---|---|---|
| 1:1 – 50:1 | **co-infection** | high | empty | **yes** — reported as the Minor strain |
| 1:500 | monoinfection | provisional | off-genotype contig note | **yes** — as a caveat |
| 500:1, 1:5000, 5000:1 | **monoinfection** | **high** | **empty** | **no** |

IVT6 (1:500) surfaces only because its minority contig is 1,395 bp, clearing the 1,000 bp threshold.
IVT9 misses by 80 bp — a 920 bp contig.

So the two limits differ and should be stated separately:

- **Co-infection reporting limit: between 1:50 and 1:500.**
- **Assembly recovery limit: at least 1:5000**, at sub-1× k-mer coverage.

---

## 5. Why the default must NOT be lowered

The obvious inference from §3 — "assembly sees 349 bp, so lower the threshold to 349 bp" — is wrong,
and this dataset refutes it directly.

These mixtures contain **only genotypes 1a and 2a**. Every contig of any other genotype is a false
positive by construction. There are six:

| sample | contig | length | k-mer | pident |
|---|---|---|---|---|
| IVT5 | `3a_D17763` | 320 bp | 1.22× | 96.2% |
| IVT2 | `3a_D17763` | 267 bp | 0.80× | 92.0% |
| IVT5 | `3a_X76918` | 257 bp | 0.85× | 94.2% |
| IVT7 | `3a_D17763` | 248 bp | 0.79× | 94.0% |
| IVT6 | `3a_D17763` | 226 bp | 1.02× | 96.9% |
| IVT4 | `7b_KX092342` | 109 bp | 652× | 100% |

Now place them beside the true minorities that currently go unreported:

| | length | k-mer |
|---|---|---|
| **true** minority, IVT7 (500:1) | 349 bp | 1.38× |
| **true** minority, IVT8 (1:5000) | 505 bp | 0.94× |
| **true** minority, IVT9 (5000:1) | 920 bp | 1.71× |
| **false** 3a, IVT5 | 320 bp | 1.22× |
| **false** 3a, IVT6 | 226 bp | 1.02× |

**The two populations are not separable on either axis.** The lowest true minority (349 bp, 1.38×) and
the highest false positive (320 bp, 1.22×) differ by 29 bp and 0.16× k-mer. A threshold low enough to
report IVT7's genuine 2a minority also admits IVT5's spurious 3a contig. K-mer coverage does not
separate them either — the true minorities sit at 0.94–1.71×, squarely inside the false population's
0.79–1.22× band.

**And this is the easy case.** The IVT samples are synthetic RNA transcripts: no patient background,
no host material, no co-infecting quasispecies, no library contamination, two known strains. Real
clinical material is not this clean, and the same Thomson cohort shows it — see §6. Any margin
measured here is an upper bound on what would survive in routine samples.

Hence: **the 1,000 bp default stays.** What this dataset licenses is a *statement about sensitivity*,
not a change to the operating point.

---

## 6. What real samples look like by comparison

From the 10-sample verification run at `d64d35a` (7 simulated + 3 Thomson), for contrast with the
controlled mixtures above:

- **`ERR1810451`** — a 4a monoinfection carrying a **1,538 bp 2b contig at k-mer 2.3, 91.8% identity**.
  It clears the 1,000 bp threshold and is flagged as *"possible missed co-infection or contamination"*.
  Whether it is a real minor strain or background is exactly the question the threshold exists to
  triage, and at 1,538 bp it is well above anything seen in the IVT false population.
- **`ERR1810469`** — called `co-infection (indeterminate dominance)`: 190 vs 1,017 nodup reads, major
  breadth 4.18% at 10×. Read-count and k-mer rankings disagree. Nothing about this sample is clean.
- **`sim22asingle` / `sim23asingle`** — even *simulated* single-genotype data attracts cross-mapping
  second candidates (`1m_KJ439778`, `1n_KJ439775`) at 6% breadth, retained as rank-2 background.

The IVT series answers "what can the assembler see under ideal conditions". It does not answer "what
can be trusted in a routine sample", and the manuscript should not let the first stand in for the
second.

---

## 7. Bonus: §2's rejected approach is now empirically confirmed

260804-dnrank §2 rejected re-ranking the minor slot by k-mer coverage, reporting that `ERR1810521`'s
flag would be lost because the selected contig would change from **1,395 bp @ 11.0–16.7×** to
**845 bp @ 18.4×**, dropping below the 1,000 bp threshold. That claim was listed as unverifiable in
§10.4 for want of the accession.

`ERR1810521` is IVT6, and its contig inventory contains **both** contigs:

```
1a_AF009606   1395 bp   k-mer 16.74   pident 98.71    <- current selection, flag fires
1a_AF009606    845 bp   k-mer 18.40   pident 99.64    <- k-mer ranking would pick this, flag lost
```

The k-mer-ranked alternative is indeed shorter *and* higher-coverage, exactly as §2 describes, and
845 bp is below 1,000. **§2's mechanism is confirmed on its own data.** Do not revisit that approach.

This also sharpens §5: the very property that makes k-mer ranking attractive — preferring high
coverage — systematically selects shorter contigs, and contig length is the only axis on which true
and false low-abundance signals separate at all.

---

## 8. Suggested manuscript framing

Two claims, kept apart:

> De novo assembly recovers the minority strain at every mixture ratio tested, from 1:1 to 1:5000,
> including a 505 bp contig at 0.94× k-mer coverage.

> Co-infection is *reported* down to approximately 1:500. Below that, the minority contig falls under
> the 1,000 bp off-genotype reporting threshold and the sample is called a monoinfection.

Then the honest reason the gap is not closed:

> The threshold is not set by assembler sensitivity but by specificity. In these controlled mixtures
> the spurious off-genotype contigs (226–320 bp, 0.79–1.22× k-mer) overlap the genuine sub-1:500
> minorities (349–920 bp, 0.94–1.71×) closely enough that no length or coverage cut separates them.
> Because clinical specimens carry more background than *in vitro* transcript mixtures, the threshold
> is deliberately conservative.

---

## 9. Suggested tests

The series is now a usable fixture: inputs are public, and the accession→ratio mapping is in §3.

1. **Assert the reported outcome per ratio, not merely that a contig exists.** A test checking only
   "minority recovered" passes on IVT7–IVT9 while the Summary says `monoinfection` / `high` / no flag.
   Assert the call, the confidence tier and the flag.
2. **IVT6 as the off-genotype-flag fixture.** Assert the flag fires with a ≥1,000 bp contig, and — per
   §7 — that the selected contig is the 1,395 bp one and not the 845 bp one. That single assertion
   locks out the §2 regression permanently.
3. **A specificity assertion.** No genotype-3 or genotype-7 candidate may ever be *reported* for these
   samples, regardless of what the assembler emits. This is the check that would fail first if the
   threshold were lowered.

---

## 10. Reproduction

Everything above is public and reproducible:

- Accessions: `ERR1810511`/`513`/`515`/`517`/`519`/`521`/`523`/`525`/`527`
- Pipeline: `folkehelseinstituttet/hcvtyper` @ `1acf637`, `-profile docker`, stock parameters
- Contig inventories: `blastparse/<sample>_blast_out.csv`, collapsed to one best hit per contig
  (lowest e-value, then highest bitscore) — the same `scaf_top` frame `blast_parse.R` ranks on
- Calls and flags: `summary/Summary.csv`, `summary/candidates.csv`
