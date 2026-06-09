---
status: complete
phase: milestone-v1
source: phases/01-04 SUMMARY.md files + full 93-sample run
started: 2026-06-09T00:00:00Z
updated: 2026-06-09T00:00:00Z
results_path: /mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCVTyper_SRA_data/summary/Summary.csv
---

## Current Test

[testing complete]

## Ground Truth Note

ERR1810469 (sP575531 in Thomson et al. 2016): confirmed 3a single infection.
The pipeline now correctly calls minor_typable=NO via the secondary major-gate.
Pre-milestone it falsely called a 3a/1a co-infection.

## Tests

### 1. sim11asingle refuted (Change 2 — de novo contradicts 4g minor)
expected: minor_denovo_status = refuted, minor_typable = NO
result: pass

### 2. ERR1810469 major-gated (Change 1 — failed major blocks minor)
expected: minor_typable = NO. Targeted major (3a) mapping shows only 248 deduplicated reads and 24.21% coverage — both below minRead (499) and minCov (29). The secondary gate fires and suppresses the minor call. minor_denovo_status = confirmed_by_denovo (de novo confirms 1a reads are real, but 3a is the artefactual "major").
result: pass

### 3. Known true co-infections preserved
expected: ERR1810447 (1b/2b), ERR1810453 (1a/2b), ERR1810475 (2b/1a), ERR1810505 (3a/1a), sim1 (1a/1b) all show minor_typable = YES and minor_denovo_status = confirmed_by_denovo. Major targeted stats well above thresholds (35k–544k reads, 98–100% coverage).
result: pass

### 4. ERR1810507 refuted (de novo short contig)
expected: ERR1810507 shows minor_typable = NO, minor_denovo_status = refuted. Was YES in pre-milestone results. Minor de novo contig is 829 bp — below the 1000 bp substantiality floor — so de novo evidence is insufficient to confirm. The pipeline correctly downgrades. (Note: strong mapping evidence with 12191 reads / 99.32% major cov; this is a threshold-calibration edge case worth reviewing.)
result: issue
reported: "This is theoretically a genuine 1a:3a (or 1a:3b) mixture — lab error or capture differences may explain the 3b→3a switch. The de novo evidence is fairly good: 829 bp contig, k-mer cov ~20, and almost the entire 3a reference covered by contigs in total. User suggests the 1000 bp threshold may be too strict — 800 bp is quite good. Total reference coverage across multiple contigs should also be taken into account, not just single-contig length."
severity: major

### 5. Simulated co-infections and singles
expected: sim1 (1a/1b) and sim2 (2a/3a) both show minor_typable = YES. Simulated single-infection samples (sim11asingle, sim11bsingle, sim22asingle, sim23asingle, sim3) all show either minor_typable = NO/UNKNOWN with no spurious co-infection call.
result: pass

### 6. minor_denovo_status column present and well-formed
expected: Summary.csv contains the minor_denovo_status column. Values across the 93 samples are exclusively from the set {confirmed_by_denovo, refuted, unconfirmed, not_evaluated, NA}. No unexpected strings.
result: pass
notes: Actual counts — 14 confirmed_by_denovo, 2 refuted (sim11asingle + ERR1810507), 77 NA. Initial expected count of 13 was wrong; ERR1810469 has confirmed_by_denovo (de novo confirms 1a reads) while minor_typable=NO is driven by the secondary major-gate, not de novo.

### 7. No new false-positive co-infections in single-infection SRA samples
expected: Samples that had no minor reference candidate in the pre-milestone run (Major_reference present, Minor_reference = NA) still show no minor call in the new results. No sample has minor_typable flipped to YES unexpectedly.
result: pass
notes: YES list = ERR1810447, ERR1810453, ERR1810475, ERR1810503, ERR1810505, ERR1810510-11-13-15-17-19, sim1, sim2. ERR1810469 + ERR1810507 correctly removed; sim1 + sim2 correctly added. No unexpected new YES.

### 8. Row count stable — no dropped samples
expected: Summary.csv has 93 rows (same count as before). All samples from the samplesheet produce a row, including samples that mapped no reads (UNKNOWN major).
result: pass

## Summary

total: 8
passed: 7
issues: 1
pending: 0
skipped: 0
blocked: 0

## Gaps

- truth: "ERR1810507 — genuine 1a:3a co-infection should be confirmed_by_denovo and minor_typable=YES"
  status: failed
  reason: "User reported: genuine mixture (possibly 1a:3b with 3b→3a reclassification); de novo contig is 829 bp (k-mer cov ~20) and cumulative contig coverage spans nearly the full 3a reference. 1000 bp single-contig floor is too strict — 800 bp would pass this case. Cumulative reference coverage across contigs is not currently considered."
  severity: major
  test: 4
  root_cause: ""
  artifacts: []
  missing: []
  debug_session: ""
