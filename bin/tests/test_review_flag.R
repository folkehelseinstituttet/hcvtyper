#!/usr/bin/env Rscript

# test_review_flag.R ------------------------------------------------------
# Phase 13 Plan 01 Task 3 (EVID-06). Unit tests for the two pure review_flag
# helpers in bin/classify_roles.R:
#   - candidate_review_fragment(): per-candidate; NA for a clean candidate, else a
#     "candidate <rank> (<subtype>_<ref>): <reason with the concrete value>" fragment.
#     Fires for weak/probable/discordant AND the D-08 demotions.
#   - sample_review_message(): sample-level; D-10 triggers stay generic (no candidate
#     name), D-11 triggers (dominant_unconfirmed, major_ref_changed) name the candidate,
#     merges the pre-collapsed candidate fragment, NA when nothing fires.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")

frag <- function(role, role_reason, evidence_state, rank = 2L,
                 ref = "2c_JX227949", subtype = "2c",
                 len = NA_real_, pid = NA_real_, kmer = NA_real_,
                 concordance_status = NA_character_) {
  candidate_review_fragment(
    role = role, role_reason = role_reason, evidence_state = evidence_state,
    concordance_status = concordance_status,
    candidate_rank = rank, candidate_ref = ref, candidate_subtype = subtype,
    best_contig_length = len, best_contig_pident = pid, best_contig_kmer_cov = kmer
  )
}

# --- candidate_review_fragment: per-trigger candidate naming + concrete value ----
# Each flaggable trigger fires in isolation and must name the candidate (rank + ref)
# and carry a concrete value (measured contig numbers where a contig exists).
weak  <- frag("background", "weak_own_assembly_below_floor", "weak", len = 200, pid = 50, kmer = 0.5)
prob  <- frag("co-infection", "corroborated", "probable", len = 3000, pid = 80, kmer = 5)
# 260805: the refuted_denovo fragment case was removed with the unreachable
# evidence_state band. discordant_identity remains — it is live via the GLUE leg —
# but its evidence_state is now `weak`, since `refuted` no longer exists.
disc  <- frag("background", "discordant_identity", "weak", len = 700, pid = 84, kmer = 4)

for (nm in c("weak", "prob", "disc")) {
  s <- get(nm)
  if (is.na(s) || !nzchar(s))
    fail(sprintf("candidate_review_fragment[%s] must fire (non-NA fragment)", nm))
  if (!grepl("candidate 2", s, fixed = TRUE) || !grepl("2c_JX227949", s, fixed = TRUE))
    fail(sprintf("candidate_review_fragment[%s] must name the candidate rank + ref, got: %s", nm, s))
}
# Concrete measured values present per state.
if (!grepl("200", weak, fixed = TRUE) || !grepl("50.0", weak, fixed = TRUE))
  fail(sprintf("weak fragment must carry the measured contig length/identity, got: %s", weak))
if (!grepl("probable", prob, fixed = TRUE))
  fail(sprintf("probable fragment must carry the evidence_state, got: %s", prob))
if (!grepl("84.0", disc, fixed = TRUE) && !grepl("700", disc, fixed = TRUE))
  fail(sprintf("discordant fragment must carry a concrete contig value, got: %s", disc))
ok("ReviewFlag-1 (EVID-06/D-05/D-07): weak/probable/discordant each name the candidate + carry a concrete value")

# --- D-08 demotion pair: each produces a NAMED fragment despite good own evidence ---
dem_same <- frag("background", "same_genotype_as_dominant", "confirmed", len = 9000, pid = 94, kmer = 30)
dem_2k1b <- frag("background", "recombinant_2k1b", "confirmed", len = 9000, pid = 95, kmer = 30)
if (is.na(dem_same) || !grepl("candidate 2", dem_same, fixed = TRUE) || !grepl("demoted", dem_same, fixed = TRUE))
  fail(sprintf("same_genotype_as_dominant must produce a named 'demoted' fragment, got: %s", dem_same))
if (!grepl("94.0", dem_same, fixed = TRUE))
  fail(sprintf("same_genotype_as_dominant fragment must carry its own contig evidence value, got: %s", dem_same))
if (is.na(dem_2k1b) || !grepl("candidate 2", dem_2k1b, fixed = TRUE) || !grepl("2k/1b", dem_2k1b, fixed = TRUE))
  fail(sprintf("recombinant_2k1b must produce a named 2k/1b fragment, got: %s", dem_2k1b))
ok("ReviewFlag-2 (D-08): same_genotype_as_dominant + recombinant_2k1b each flag with a named fragment despite good evidence")

# --- Clean candidates produce NO fragment (D-06 clean path) --------------------
if (!is.na(frag("dominant", "dominant", "confirmed", len = 9189, pid = 100, kmer = 40)))
  fail("a clean dominant candidate must yield NA (no fragment)")
if (!is.na(frag("co-infection", "corroborated", "confirmed", len = 9479, pid = 89, kmer = 514)))
  fail("a clean confirmed co-infection must yield NA (no fragment)")
if (!is.na(frag("indeterminate", "indeterminate_dominance_conflict", "confirmed", len = 8000, pid = 92, kmer = 20)))
  fail("indeterminate_dominance_conflict is a sample-level generic message, not a per-candidate fragment (D-10)")
ok("ReviewFlag-3 (D-06/D-10): clean dominant/confirmed-co-infection and the dominance-conflict candidate yield NO fragment")

# --- sample_review_message: D-10 triggers stay GENERIC (no rank/ref token) ------
msg_indet <- sample_review_message(
  overall_sample_call = "indeterminate",
  denovo_major_subtype_match = NA, denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = NA,
  major_subtype = NA, rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref")
msg_gate <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = NA, denovo_minor_subtype_match = NA,
  gate_flag = "low_cov", denovo_minor_subtype = NA, denovo_major_subtype = NA,
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref")
msg_idom <- sample_review_message(
  overall_sample_call = "co-infection (indeterminate dominance)",
  denovo_major_subtype_match = NA, denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = NA,
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref")
for (nm in c("msg_indet", "msg_gate", "msg_idom")) {
  s <- get(nm)
  if (is.na(s) || !nzchar(s)) fail(sprintf("D-10 trigger %s must fire a generic message", nm))
  # Generic prose may contain the word "candidate" ("No candidate passed…"); what it
  # must NOT contain is a candidate IDENTITY token — "candidate <rank>" or a ref token.
  if (grepl("candidate [0-9]", s) || grepl("3a_ref", s, fixed = TRUE))
    fail(sprintf("D-10 trigger %s must stay GENERIC (no candidate rank/ref token), got: %s", nm, s))
}
ok("ReviewFlag-4 (D-10): is_indet / gate_flag / is_indet_dom stay generic — no candidate name injected")

# --- D-11: dominant_unconfirmed + major_ref_changed NAME their candidate --------
msg_unconf <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = NA, denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = NA,
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = TRUE,
  dominant_rank = 1L, dominant_ref = "3a_ref")
if (is.na(msg_unconf) || !grepl("candidate 1", msg_unconf, fixed = TRUE) || !grepl("3a_ref", msg_unconf, fixed = TRUE))
  fail(sprintf("dominant_unconfirmed must NAME the dominant candidate (D-11), got: %s", msg_unconf))
msg_resc <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = NA, denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = NA,
  major_subtype = "3a", rescue_effect = "major_ref_changed", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref")
if (is.na(msg_resc) || !grepl("candidate 1", msg_resc, fixed = TRUE) || !grepl("3a_ref", msg_resc, fixed = TRUE))
  fail(sprintf("major_ref_changed must NAME the reassigned candidate (D-11), got: %s", msg_resc))
ok("ReviewFlag-5 (D-11): dominant_unconfirmed + major_ref_changed name their candidate")

# --- Nothing fires -> NA_character_ (MultiQC NA sentinel) ----------------------
msg_clean <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "YES", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = "3a",
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref", candidate_fragment = NA_character_)
if (!is.na(msg_clean))
  fail(sprintf("a clean sample with no trigger must return NA_character_, got: %s", msg_clean))
ok("ReviewFlag-6 (Pitfall 2): a clean sample returns NA_character_, never a literal 'NA'")

# --- Success criterion 4: monoinfection with ONE non-dominant flagged candidate --
# still surfaces that candidate's fragment via the merged candidate_fragment arg.
cand_frag <- frag("background", "weak_own_assembly_below_floor", "weak",
                  rank = 2L, ref = "4d_DQ418786", subtype = "4d", len = 200, pid = 50, kmer = 0.5)
msg_mono_flag <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "YES", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = "3a",
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_ref", candidate_fragment = cand_frag)
if (is.na(msg_mono_flag) || !grepl("candidate 2", msg_mono_flag, fixed = TRUE) ||
    !grepl("4d_DQ418786", msg_mono_flag, fixed = TRUE))
  fail(sprintf("monoinfection with a non-dominant flagged candidate must surface that fragment (criterion 4), got: %s", msg_mono_flag))
ok("ReviewFlag-7 (criterion 4): a monoinfection sample with one non-dominant flagged candidate still surfaces its named fragment")

# --- Production collapse shape: pmap_chr over candidate rows -> per-sample collapse --
cand_rows <- tibble(
  role = c("dominant", "background", "co-infection"),
  role_reason = c("dominant", "no_own_assembly", "corroborated"),
  evidence_state = c("confirmed", "weak", "probable"),
  concordance_status = NA_character_,
  candidate_rank = c(1L, 2L, 3L),
  candidate_ref = c("3a_ref", "4d_DQ418786", "2c_JX227949"),
  candidate_subtype = c("3a", "4d", "2c"),
  best_contig_length = c(9189, NA, 3000),
  best_contig_pident = c(100, NA, 80),
  best_contig_kmer_cov = c(40, NA, 5)
)
frags <- pmap_chr(cand_rows, candidate_review_fragment)
collapsed <- paste(frags[!is.na(frags)], collapse = " | ")
if (sum(!is.na(frags)) != 2)
  fail(sprintf("expected exactly 2 flagged fragments (weak + probable), got %d", sum(!is.na(frags))))
if (!grepl("candidate 2", collapsed, fixed = TRUE) || !grepl("candidate 3", collapsed, fixed = TRUE))
  fail(sprintf("collapsed fragment must retain both flagged candidates, got: %s", collapsed))
ok("ReviewFlag-8: candidate_review_fragment vectorizes via pmap_chr and collapses NA-filtered per sample")

# --- 260803-ogc: offgenotype_contig_reviewable() -------------------------------
# The gate for the monoinfection different-genotype-contig review sentence. Two
# legs: the 2k1b recombinant pair exclusion (is_valid_minor rule 2) and the
# contig-length floor. Empirical basis: the 2026-08-03 threshold sweep over 5 runs
# / 140 samples, on which the flag fired on 72 samples with no floor applied.

ogc <- offgenotype_contig_reviewable

# Genuine off-genotype contigs still fire.
if (!ogc("1b", "3a")) fail("OGC-1: a 1b contig against a 3a major must remain reviewable")
if (!ogc("3a", "1a")) fail("OGC-1: a 3a contig against a 1a major must remain reviewable")
# Same genotype never fires (the within-genotype-1 artefact class stays suppressed:
# this is why the trigger must NOT delegate wholesale to is_valid_minor(), whose
# rule 1 admits 1a/1b as a valid co-infection pair).
if (ogc("1a", "1b")) fail("OGC-1: a 1a contig against a 1b major must NOT fire (rule-1 trap)")
if (ogc("3a", "3a")) fail("OGC-1: same subtype must not fire")
ok("OGC-1: genotype-difference test keeps genuine off-genotype contigs and rejects within-genotype pairs")

# 2k1b recombinant pairs are excluded against genotype 1 and 2, in both directions.
for (mj in c("1a", "1b", "2c", "2k1b")) {
  if (ogc("2k1b", mj))
    fail(sprintf("OGC-2: a 2k1b contig against a %s major must be excluded (is_valid_minor rule 2)", mj))
  if (ogc(mj, "2k1b"))
    fail(sprintf("OGC-2: a %s contig against a 2k1b major must be excluded (is_valid_minor rule 2)", mj))
}
# ...but 2k1b against a genotype OUTSIDE {1,2,2k1b} is a real conflict and survives
# (observed: Sample2K-3a, a 1116 bp 2k1b contig against a 3a major).
if (!ogc("2k1b", "3a")) fail("OGC-2: a 2k1b contig against a 3a major must still fire")
if (!ogc("3a", "2k1b")) fail("OGC-2: a 3a contig against a 2k1b major must still fire")
ok("OGC-2: 2k1b pairs excluded against genotype {1,2,2k1b} only, both directions")

# The genotype_from_subtype()-alone trap: swapping the old substr(x,1,1) comparison
# for genotype_from_subtype() does NOT suppress 2k1b, because it maps "2k1b" ->
# "2k1b" and so "2k1b" != "1" still differs. Assert the pair test is what does the
# work, and that we did not regress to the naive comparison.
if (as.character(genotype_from_subtype("2k1b")) != "2k1b")
  fail("OGC-3: genotype_from_subtype('2k1b') must return '2k1b' (2k1b-aware rule)")
if (genotype_from_subtype("2k1b") == genotype_from_subtype("1b"))
  fail("OGC-3: precondition — gfs('2k1b') and gfs('1b') must differ, which is why the pair test is required")
if (ogc("2k1b", "1b"))
  fail("OGC-3: gfs alone would let 2k1b-vs-1b through; the explicit pair test must suppress it")
# substr() would NOT have fired on 2k1b-vs-2c ("2" == "2"); gfs does, and the pair
# test must suppress it so behaviour is unchanged for that case.
if (ogc("2k1b", "2c"))
  fail("OGC-3: 2k1b-vs-2c must stay unflagged (substr parity), suppressed by the pair test not by gfs")
ok("OGC-3: the 2k1b suppression comes from the explicit pair test, not from genotype_from_subtype()")

# Contig-length floor.
if (!ogc("1b", "3a", contig_length = 1000, min_length = 1000))
  fail("OGC-4: a contig exactly at the floor must be reviewable (>=, not >)")
if (ogc("1b", "3a", contig_length = 999, min_length = 1000))
  fail("OGC-4: a contig below the floor must be suppressed")
if (!ogc("1b", "3a", contig_length = 142, min_length = 0))
  fail("OGC-4: min_length = 0 must reproduce the pre-260803-ogc behaviour")
# The three must-keep samples: the legacy typable=YES minors this build demotes to
# monoinfection. Any floor that drops one of these is disqualified.
for (L in c(4467, 2787, 2706)) {
  if (!ogc("1b", "3a", contig_length = L, min_length = 1000))
    fail(sprintf("OGC-4: must-keep contig of %d bp must survive the 1000 bp floor", L))
}
ok("OGC-4: contig-length floor is inclusive, suppresses sub-floor contigs, retains all three must-keeps")

# NA handling: an unmeasurable length fails OPEN (keeps the flag); an NA subtype
# never fires.
if (!ogc("1b", "3a", contig_length = NA_real_, min_length = 1000))
  fail("OGC-5: NA contig length must fail OPEN and keep the flag")
if (ogc(NA_character_, "3a", contig_length = 5000, min_length = 1000))
  fail("OGC-5: NA contig subtype must never fire")
if (ogc("1b", NA_character_, contig_length = 5000, min_length = 1000))
  fail("OGC-5: NA major subtype must never fire")
ok("OGC-5: NA contig length fails open; NA subtypes never fire")

# Vectorised over a frame, returning a bare logical with no NAs (it is used inside
# a mutate() across the whole sample frame).
vec <- ogc(
  contig_subtype = c("1b",  "2k1b", "1b",  NA,    "3a"),
  major_subtype  = c("3a",  "1b",   "3a",  "3a",  "3a"),
  contig_length  = c(2706,  9000,   606,   9000,  NA),
  min_length     = 1000
)
if (!is.logical(vec) || anyNA(vec))
  fail("OGC-6: must return a bare logical vector with no NAs")
if (!identical(vec, c(TRUE, FALSE, FALSE, FALSE, FALSE)))
  fail(sprintf("OGC-6: vectorised result wrong, got: %s", paste(vec, collapse = ",")))
ok("OGC-6: vectorises over a frame and returns a bare NA-free logical")

# End-to-end through sample_review_message(): the caller masks the length leg, so a
# sub-floor contig arrives as NA and no sentence is emitted; a substantial one fires.
srm_mono <- function(dv_minor, maj = "3a") sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "YES", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = dv_minor, denovo_major_subtype = maj,
  major_subtype = maj, rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = paste0(maj, "_ref"), candidate_fragment = NA_character_)

msg_fires <- srm_mono("1b")
if (is.na(msg_fires) || !grepl("different-genotype contig (1b)", msg_fires, fixed = TRUE))
  fail(sprintf("OGC-7: a substantial off-genotype contig must still emit the sentence, got: %s", msg_fires))
# The length leg is applied by the caller, which masks a sub-floor contig to NA.
if (!is.na(srm_mono(NA_character_)))
  fail("OGC-7: a length-masked (NA) contig must emit no sentence — the caller's floor is what suppresses")
# The 2k1b leg is applied INSIDE the helper: 2k1b vs a 1b major is suppressed...
if (!is.na(srm_mono("2k1b", maj = "1b")))
  fail("OGC-7: a 2k1b contig against a 1b major must emit no sentence (pair rule)")
# ...while 2k1b vs a 3a major is a real conflict and still fires (Sample2K-3a).
if (is.na(srm_mono("2k1b", maj = "3a")))
  fail("OGC-7: a 2k1b contig against a 3a major must still emit the sentence")
ok("OGC-7: sample_review_message() honours the masked slot, the pair rule, and still fires on a real conflict")

# --- 260803-ogc option C: contig_evidence_note() -------------------------------
# The measured-evidence clause. Real numbers from run 20251212-01 / 20260521-01.

ocn <- contig_evidence_note

# SampleSA-1a: 1620 bp contig, 69 bp aligned (4%). Independently confirmed against nt —
# the contig's only HCV-like region is a ~212 bp tail closest to 1a, the same
# genotype as the major, so there is no second strain and the 6i label is an
# artefact of a short anchor.
note_junk <- ocn(contig_length = 1620, aln_length = 69, pident = 91.30, kmer_cov = 1.02)
for (tok in c("1620 bp contig", "69 bp aligned", "(4%)", "91.3% identity", "k-mer cov 1.0")) {
  if (!grepl(tok, note_junk, fixed = TRUE))
    fail(sprintf("OGC-8: note must carry '%s', got: %s", tok, note_junk))
}
if (!grepl("weakly supported", note_junk, fixed = TRUE))
  fail(sprintf("OGC-8: a 4%%-aligned contig must be described as weakly supported, got: %s", note_junk))
if (grepl("possible missed co-infection", note_junk, fixed = TRUE))
  fail(sprintf("OGC-8: a 4%%-aligned contig must NOT assert a possible co-infection, got: %s", note_junk))
ok("OGC-8: a poorly-aligned contig carries all four metrics and is described as weakly supported")

# SampleFL-1a: 4467 bp contig, 4448 bp aligned (100%) — a genuine co-infection signal.
note_real <- ocn(contig_length = 4467, aln_length = 4448, pident = 91.93, kmer_cov = 1.42)
if (!grepl("(100%)", note_real, fixed = TRUE))
  fail(sprintf("OGC-9: a fully-aligned contig must report 100%%, got: %s", note_real))
if (!grepl("possible missed co-infection or contamination", note_real, fixed = TRUE))
  fail(sprintf("OGC-9: a fully-aligned contig must keep the co-infection wording, got: %s", note_real))
if (grepl("weakly supported", note_real, fixed = TRUE))
  fail(sprintf("OGC-9: a fully-aligned contig must NOT be called weakly supported, got: %s", note_real))
# All three must-keeps align >=99% and must all read as genuine.
for (m in list(c(4467, 4448), c(2787, 2762), c(2706, 2705))) {
  n <- ocn(contig_length = m[1], aln_length = m[2], pident = 92, kmer_cov = 1.5)
  if (grepl("weakly supported", n, fixed = TRUE))
    fail(sprintf("OGC-9: must-keep contig %d/%d must not be called weakly supported", m[2], m[1]))
}
ok("OGC-9: fully-aligned contigs keep the co-infection wording; all three must-keeps read as genuine")

# Partial metrics degrade gracefully; nothing known yields NA so the caller falls back.
if (!is.na(ocn(1620, NA, NA, NA)) && !grepl("1620 bp contig", ocn(1620, NA, NA, NA), fixed = TRUE))
  fail("OGC-10: a note with only a contig length must still render that length")
if (grepl("aligned", ocn(1620, NA, NA, NA), fixed = TRUE))
  fail("OGC-10: an unknown alignment length must not be rendered")
if (!is.na(ocn(NA, NA, NA, NA)))
  fail("OGC-10: with no metrics at all the note must be NA so the caller falls back")
# A zero-length contig must not divide by zero.
if (grepl("%)", ocn(0, 0, 90, 1), fixed = TRUE))
  fail("OGC-10: a zero-length contig must not render a fraction")
ok("OGC-10: partial metrics degrade gracefully, no metrics yields NA, no divide-by-zero")

# End-to-end: the clause reaches review_flag, and the fallback wording is preserved
# when the metrics are unavailable.
srm_note <- function(note) sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "YES", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = "6i", denovo_major_subtype = "1a",
  major_subtype = "1a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "1a_AF009606", candidate_fragment = NA_character_,
  offgeno_note = note)

msg_annotated <- srm_note(note_junk)
if (!grepl("different-genotype contig (6i)", msg_annotated, fixed = TRUE) ||
    !grepl("69 bp aligned (4%)", msg_annotated, fixed = TRUE) ||
    !grepl("Please review.", msg_annotated, fixed = TRUE))
  fail(sprintf("OGC-11: the annotated sentence must name the subtype, carry the metrics and end with the call to action, got: %s", msg_annotated))
msg_bare <- srm_note(NA_character_)
if (!grepl("— possible missed co-infection or contamination. Please review.", msg_bare, fixed = TRUE))
  fail(sprintf("OGC-11: with no metrics the original wording must be preserved verbatim, got: %s", msg_bare))
# Default argument: an OLD-STYLE call omitting offgeno_note entirely must still work.
msg_legacy <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "YES", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = "6i", denovo_major_subtype = "1a",
  major_subtype = "1a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "1a_AF009606", candidate_fragment = NA_character_)
if (!identical(msg_legacy, msg_bare))
  fail("OGC-11: omitting offgeno_note must be identical to passing NA (backward compatibility)")
ok("OGC-11: the clause reaches review_flag, the bare wording is preserved, and the arg is back-compatible")

# --- 260803-ogc follow-up 1: the MAJOR subtype-conflict sentence ---------------
# Same helper, context = "major_conflict". This trigger sets call_confidence =
# "review" on its own, so it is the hardest sample-level signal, and it named two
# subtypes and no numbers. Note the metrics cannot come from Major_best_contig_*:
# those describe the CANDIDATE's genotype group, whereas the contig that caused the
# disagreement belongs to denovo_major_subtype's group.

# A short-anchor conflict must be described as possibly artefactual, NOT as a
# reference mismatch — a 69 bp anchor cannot support a subtype call, so it cannot
# support a disagreement with one either.
mc_weak <- contig_evidence_note(1816, 120, 89.0, 1.3, context = "major_conflict")
if (!grepl("120 bp aligned (7%)", mc_weak, fixed = TRUE))
  fail(sprintf("OGC-12: major-conflict note must carry the aligned fraction, got: %s", mc_weak))
if (!grepl("artefact of a short anchor", mc_weak, fixed = TRUE))
  fail(sprintf("OGC-12: a short-anchor conflict must be flagged as possibly artefactual, got: %s", mc_weak))
if (grepl("largely non-HCV", mc_weak, fixed = TRUE))
  fail(sprintf("OGC-12: major-conflict context must not reuse the off-genotype wording, got: %s", mc_weak))

# A full-length conflict keeps the original interpretation.
mc_real <- contig_evidence_note(4000, 3950, 92.0, 8.0, context = "major_conflict")
if (!grepl("possible reference mismatch or highly divergent strain", mc_real, fixed = TRUE))
  fail(sprintf("OGC-12: a full-length conflict must keep the reference-mismatch wording, got: %s", mc_real))
if (grepl("artefact", mc_real, fixed = TRUE))
  fail(sprintf("OGC-12: a full-length conflict must not be called artefactual, got: %s", mc_real))

# The two contexts must not bleed into each other.
og_weak <- contig_evidence_note(1816, 120, 89.0, 1.3, context = "offgenotype")
if (grepl("artefact of a short anchor", og_weak, fixed = TRUE))
  fail(sprintf("OGC-12: off-genotype context must not use the major-conflict wording, got: %s", og_weak))
if (!grepl("largely non-HCV", og_weak, fixed = TRUE))
  fail(sprintf("OGC-12: off-genotype context must keep its own wording, got: %s", og_weak))
ok("OGC-12: contig_evidence_note honours both contexts and does not cross-contaminate their wording")

# End-to-end through sample_review_message(): the conflict sentence carries the
# clause, falls back verbatim without it, and stays back-compatible when omitted.
srm_conf <- function(note, omit = FALSE) {
  a <- list(overall_sample_call = "monoinfection",
            denovo_major_subtype_match = "NO", denovo_minor_subtype_match = NA,
            gate_flag = "ok", denovo_minor_subtype = NA, denovo_major_subtype = "1a",
            major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
            dominant_rank = 1L, dominant_ref = "3a_D17763",
            candidate_fragment = NA_character_)
  if (!omit) a$majconf_note <- note
  do.call(sample_review_message, a)
}
mc_msg <- srm_conf(mc_weak)
if (!grepl("Major subtype conflict for candidate 1 (3a_D17763)", mc_msg, fixed = TRUE) ||
    !grepl("mapping (3a) vs contig (1a)", mc_msg, fixed = TRUE) ||
    !grepl("120 bp aligned (7%)", mc_msg, fixed = TRUE) ||
    !grepl("Please review.", mc_msg, fixed = TRUE))
  fail(sprintf("OGC-13: the annotated conflict sentence must keep its identity tokens and carry the metrics, got: %s", mc_msg))
mc_bare <- srm_conf(NA_character_)
if (!grepl("vs contig (1a). Possible reference mismatch or highly divergent strain. Please review.",
           mc_bare, fixed = TRUE))
  fail(sprintf("OGC-13: without metrics the original wording must be preserved verbatim, got: %s", mc_bare))
if (!identical(srm_conf(NULL, omit = TRUE), mc_bare))
  fail("OGC-13: omitting majconf_note must equal passing NA (backward compatibility)")
# Both contig sentences can fire on one sample and must each get their own clause.
both <- sample_review_message(
  overall_sample_call = "monoinfection",
  denovo_major_subtype_match = "NO", denovo_minor_subtype_match = NA,
  gate_flag = "ok", denovo_minor_subtype = "6i", denovo_major_subtype = "1a",
  major_subtype = "3a", rescue_effect = "none", dominant_unconfirmed = FALSE,
  dominant_rank = 1L, dominant_ref = "3a_D17763", candidate_fragment = NA_character_,
  offgeno_note = note_junk, majconf_note = mc_weak)
if (!grepl("Major subtype conflict", both, fixed = TRUE) ||
    !grepl("different-genotype contig (6i)", both, fixed = TRUE))
  fail(sprintf("OGC-13: both contig sentences must survive the merge, got: %s", both))
if (!grepl("69 bp aligned (4%)", both, fixed = TRUE) ||
    !grepl("120 bp aligned (7%)", both, fixed = TRUE))
  fail(sprintf("OGC-13: each sentence must carry ITS OWN contig's metrics, not the other's, got: %s", both))
ok("OGC-13: the conflict sentence is annotated, falls back verbatim, and both contig sentences keep their own metrics")

cat("ALL PASS\n")
