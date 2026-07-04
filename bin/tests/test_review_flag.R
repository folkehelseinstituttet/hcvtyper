#!/usr/bin/env Rscript

# test_review_flag.R ------------------------------------------------------
# Phase 13 Plan 01 Task 3 (EVID-06). Unit tests for the two pure review_flag
# helpers in bin/classify_roles.R:
#   - candidate_review_fragment(): per-candidate; NA for a clean candidate, else a
#     "candidate <rank> (<subtype>_<ref>): <reason with the concrete value>" fragment.
#     Fires for weak/probable/refuted/discordant AND the D-08 demotions.
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
refu  <- frag("background", "refuted_denovo", "refuted", len = 600, pid = 85, kmer = 3)
disc  <- frag("background", "discordant_identity", "refuted", len = 700, pid = 84, kmer = 4)

for (nm in c("weak", "prob", "refu", "disc")) {
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
if (!grepl("85.0", refu, fixed = TRUE) && !grepl("600", refu, fixed = TRUE))
  fail(sprintf("refuted fragment must carry a concrete contig value, got: %s", refu))
if (!grepl("84.0", disc, fixed = TRUE) && !grepl("700", disc, fixed = TRUE))
  fail(sprintf("discordant fragment must carry a concrete contig value, got: %s", disc))
ok("ReviewFlag-1 (EVID-06/D-05/D-07): weak/probable/refuted/discordant each name the candidate + carry a concrete value")

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

cat("ALL PASS\n")
