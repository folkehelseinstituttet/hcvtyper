#!/usr/bin/env Rscript

# test_rescue_evaluation.R -----------------------------------------------
# Subprocess-contract test (Phase 10, denovo-subtype-rescue) for the REAL
# de-novo subtype-rescue evaluator in bin/rescue_evaluation.R.
#
# rescue_evaluation.R now uses assembly_support.csv with subtype-level
# matching instead of the legacy blastparse.csv major/minor rank-slot
# assignment. For each candidate it looks up the candidate's own subtype in
# assembly_support: if the own subtype passes floors, the de novo confirms
# the reference and no rescue fires. Otherwise it finds the best-supported
# alternative subtype (longest contig, different subtype) and rescues to
# that subtype's best_ref if it passes the floors.
#
# The rescue contract encoded here:
#
#   (1) rescue-fires : own subtype absent from support AND best alternative
#       subtype passes all four floors -> rescued_from set, rescue_trigger
#       non-NA, confirmation_status forced to "pass", _cand{rank}. FASTA
#       written for the rescue reference.
#   (2) floor-fail   : alternative subtype fails a floor -> rescued_from NA,
#       confirmation_status unchanged.
#   (3) 2k1b-rule    : a 2k1b candidate over a genotype-2 contig meeting the
#       floors triggers rescue to the genotype-2 reference (D-03 special
#       recombinant rule).
#   (4) 1a1b-floor   : the 1a/1b boundary uses the higher rescue_1a1b_length
#       floor (D-04).
#   (5) skip-assembly: empty support -> candidates pass through unchanged
#       with rescued_from / rescue_trigger NA-filled (D-10 DoS guard).
#   (6) forced-pass  : a below_threshold candidate rescued from its slot gets
#       confirmation_status forced to "pass" (D-05/D-06).
#   (7) stale-passthrough: ref-changing rescue removes the stale module-staged
#       pass-through FASTA so no duplicate @SQ line reaches BOWTIE2_BUILD.
#   (8) collapse-guard: rescue blocked when rescue_ref is already the original
#       ref of another candidate slot (ERR1810469 regression).
#   (9) own-subtype-confirmed: no rescue fires when the candidate's own
#       subtype has strong assembly support (floors pass).
#  (10) same-genotype-blocked: rescue to 1a_AF009606 blocked because another
#       candidate slot already holds genotype 1 (1a_EF407457).
#  (11) 1a1b-exception: rescue to 1b_D90208 fires even though cand1 is 1a
#       (1a/1b co-infection pair is an explicitly permitted exception).
#
# D-02 default thresholds: length 3000, pident 85, aln_length 3000,
# kmer_cov 2, 1a1b_length 5000.
#
# Run from any cwd:
#   Rscript bin/tests/test_rescue_evaluation.R
# Exits 0 and prints "ALL PASS" when all assertions hold.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
script   <- file.path(bin_dir, "rescue_evaluation.R")

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")

if (!file.exists(script)) {
  cat("RED: bin/rescue_evaluation.R does not exist yet.\n")
  cat("FAIL: rescue_evaluation.R is not yet implemented\n")
  quit(status = 1)
}

# The 8-column long candidates schema + two rescue audit columns.
CAND_IN_COLS <- c("sample", "candidate_rank", "candidate_ref", "candidate_subtype",
                  "candidate_genotype", "candidate_reads", "candidate_cov",
                  "confirmation_status")
RESCUE_COLS  <- c("rescued_from", "rescue_trigger")

# assembly_support.csv schema — now includes best_ref (the closest reference
# for the best contig of each subtype, used as the rescue target).
SUPPORT_COLS <- c("sample", "subtype", "best_ref", "best_contig_length",
                  "best_contig_pident", "best_contig_aln_length",
                  "best_contig_kmer_cov")

TH_LENGTH   <- 3000
TH_PIDENT   <- 85
TH_ALN      <- 3000
TH_KMER     <- 2
TH_1A1B_LEN <- 5000

# -------------------------------------------------------------------------
# Subprocess harness: stage candidates / support / references fixtures in a
# tempdir, invoke the REAL rescue_evaluation.R, read back the rewritten
# <prefix>.rescued.candidates.csv and list any written FASTA files.
#
#   cands   : tibble with the 8-column long candidate schema.
#   support : tibble with the 7-column assembly_support schema.
#   refs    : named character vector subtype-tagged accession -> dummy sequence.
#
# Returns list(exit, cands_out, fastas, path).
# -------------------------------------------------------------------------
run_rescue <- function(case, prefix, cands, support, refs, prestage = character(0)) {
  wd <- tempfile(paste0("rescue_", case, "_")); dir.create(wd)

  for (f in prestage) writeLines(c(paste0(">stub"), "ACGT"), file.path(wd, f))

  cand_path <- file.path(wd, paste0(prefix, ".candidates.in.csv"))
  write_csv(cands, cand_path)

  support_path <- file.path(wd, paste0(prefix, ".assembly_support.csv"))
  write_csv(support, support_path)

  ref_lines <- unlist(lapply(names(refs), function(h) c(paste0(">", h), refs[[h]])))
  refs_path <- file.path(wd, paste0(prefix, "_references.fa"))
  writeLines(ref_lines, refs_path)

  old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
  exit <- system2(
    "Rscript",
    c(shQuote(script), shQuote(prefix),
      shQuote(cand_path), shQuote(support_path),
      shQuote(refs_path),
      TH_LENGTH, TH_PIDENT, TH_ALN, TH_KMER, TH_1A1B_LEN),
    stdout = FALSE, stderr = FALSE
  )

  out_path  <- file.path(wd, paste0(prefix, ".rescued.candidates.csv"))
  cands_out <- if (file.exists(out_path)) read_csv(out_path, show_col_types = FALSE) else NULL
  fastas    <- list.files(wd, pattern = "\\.fa$", full.names = FALSE)
  list(exit = exit, cands_out = cands_out, fastas = fastas, path = out_path)
}

mk_cand <- function(rank, ref, subtype, genotype, reads, cov, status) {
  tibble(
    candidate_rank      = as.integer(rank),
    candidate_ref       = ref,
    candidate_subtype   = subtype,
    candidate_genotype  = genotype,
    candidate_reads     = reads,
    candidate_cov       = cov,
    confirmation_status = status
  )
}
mk_cands <- function(sample, ...) {
  bind_rows(...) %>% mutate(sample = sample) %>% select(all_of(CAND_IN_COLS))
}

mk_support <- function(sample, ...) {
  bind_rows(...) %>% mutate(sample = sample) %>% select(all_of(SUPPORT_COLS))
}
# best_ref: the database reference closest to the best contig for this subtype.
mk_support_row <- function(subtype, best_ref, length, pident, aln, kmer) {
  tibble(subtype = subtype, best_ref = best_ref,
         best_contig_length = length, best_contig_pident = pident,
         best_contig_aln_length = aln, best_contig_kmer_cov = kmer)
}

REFS <- list(
  "2b_AY232748"  = strrep("A", 60),
  "2k1b_AF177036" = strrep("C", 60),
  "1a_M62321"    = strrep("G", 60),
  "1b_D90208"    = strrep("T", 60),
  "2a_AB047639"  = strrep("A", 60)
)

REFS_DUPGENO <- c(REFS, list(
  "1a_EF407457"  = strrep("G", 60),
  "1a_AF009606"  = strrep("C", 60),
  "1c_AY651061"  = strrep("T", 60)
))

empty_support <- tibble(
  sample = character(0), subtype = character(0), best_ref = character(0),
  best_contig_length = double(0), best_contig_pident = double(0),
  best_contig_aln_length = double(0), best_contig_kmer_cov = double(0)
)

assert_schema <- function(case, cands_out) {
  if (is.null(cands_out))
    fail(paste(case, ": no <prefix>.rescued.candidates.csv written"))
  missing <- setdiff(c(CAND_IN_COLS, RESCUE_COLS), colnames(cands_out))
  if (length(missing) > 0)
    fail(paste(case, ": output missing columns:", paste(missing, collapse = ",")))
}

# =========================================================================
# Subtest 1: rescue-fires
# =========================================================================
# Candidate is 1a but own subtype has no assembly support; the best alternative
# is 2b with a contig clearing all four floors -> rescue to 2b_AY232748.
r1 <- run_rescue(
  "fires", "FIRES",
  cands = mk_cands("FIRES",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("FIRES",
    mk_support_row("2b", "2b_AY232748", 4120, 96.2, 3840, 4.1)),
  refs = REFS
)
assert_schema("fires", r1$cands_out)
row1 <- r1$cands_out %>% filter(candidate_rank == 1)
if (is.na(row1$rescued_from[1]))
  fail("fires: rescued_from must be set when rescue fires")
if (row1$rescued_from[1] != "1a_M62321")
  fail(paste("fires: rescued_from must be the original ref 1a_M62321, got", row1$rescued_from[1]))
if (is.na(row1$rescue_trigger[1]))
  fail("fires: rescue_trigger must be non-NA when rescue fires")
if (row1$confirmation_status[1] != "pass")
  fail(paste("fires: confirmation_status must be forced to 'pass', got", row1$confirmation_status[1]))
if (!any(grepl("_cand1\\.", r1$fastas)))
  fail(paste("fires: a _cand1. FASTA must be written; got:", paste(r1$fastas, collapse = ",")))
ok("rescue-fires -> rescued_from set, rescue_trigger non-NA, status forced pass, _cand1. FASTA written")

# =========================================================================
# Subtest 2: floor-fail (pident below threshold -> no rescue)
# =========================================================================
r2 <- run_rescue(
  "floorfail", "FLOORFAIL",
  cands = mk_cands("FLOORFAIL",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("FLOORFAIL",
    mk_support_row("2b", "2b_AY232748", 4120, 80.0, 3840, 4.1)),  # pident 80 < 85
  refs = REFS
)
assert_schema("floorfail", r2$cands_out)
f2 <- r2$cands_out %>% filter(candidate_rank == 1)
if (!is.na(f2$rescued_from[1]))
  fail(paste("floorfail: rescued_from must be NA when a floor fails, got", f2$rescued_from[1]))
ok("floor-fail -> rescued_from NA (no rescue on sub-floor pident)")

# =========================================================================
# Subtest 3: 2k1b special rule (D-03)
# =========================================================================
r3 <- run_rescue(
  "k2b", "K2B",
  cands = mk_cands("K2B",
    mk_cand(1, "2k1b_AF177036", "2k1b", "2k1b", 60000, 92, "pass")),
  support = mk_support("K2B",
    mk_support_row("2a", "2a_AB047639", 4500, 97.0, 4200, 6.0)),
  refs = REFS
)
assert_schema("k2b", r3$cands_out)
k3 <- r3$cands_out %>% filter(candidate_rank == 1)
if (is.na(k3$rescued_from[1]))
  fail("2k1b-rule: a 2k1b candidate over a genotype-2 contig meeting floors must rescue")
if (k3$rescued_from[1] != "2k1b_AF177036")
  fail(paste("2k1b-rule: rescued_from must be 2k1b_AF177036, got", k3$rescued_from[1]))
if (is.na(k3$rescue_trigger[1]))
  fail("2k1b-rule: rescue_trigger must be non-NA")
ok("2k1b-rule -> genotype-2 contig over a 2k1b candidate triggers rescue (D-03)")

# =========================================================================
# Subtest 4: 1a/1b boundary uses the higher length floor (D-04)
# =========================================================================
r4a <- run_rescue(
  "ab_short", "ABSHORT",
  cands = mk_cands("ABSHORT",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("ABSHORT",
    mk_support_row("1b", "1b_D90208", 3500, 99.0, 3400, 8.0)),  # 3500 < 5000 floor
  refs = REFS
)
assert_schema("ab_short", r4a$cands_out)
a4 <- r4a$cands_out %>% filter(candidate_rank == 1)
if (!is.na(a4$rescued_from[1]))
  fail(paste("1a1b-floor: 3500bp must NOT rescue (below 5000bp floor), got rescued_from",
             a4$rescued_from[1]))

r4b <- run_rescue(
  "ab_long", "ABLONG",
  cands = mk_cands("ABLONG",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("ABLONG",
    mk_support_row("1b", "1b_D90208", 5200, 99.0, 5100, 8.0)),  # 5200 >= 5000 floor
  refs = REFS
)
assert_schema("ab_long", r4b$cands_out)
b4 <- r4b$cands_out %>% filter(candidate_rank == 1)
if (is.na(b4$rescued_from[1]))
  fail("1a1b-floor: 5200bp (>= 5000bp floor) must rescue")
ok("1a1b-floor -> 3500bp no rescue, 5200bp rescues (D-04)")

# =========================================================================
# Subtest 5: skip-assembly (empty support -> pass-through with NA columns)
# =========================================================================
r5 <- run_rescue(
  "skip", "SKIP",
  cands = mk_cands("SKIP",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass"),
    mk_cand(2, "1b_D90208", "1b", "1", 40000, 70, "below_threshold")),
  support = empty_support,
  refs = REFS
)
if (r5$exit != 0)
  fail(paste("skip-assembly: script must exit 0, got exit", r5$exit))
assert_schema("skip", r5$cands_out)
if (nrow(r5$cands_out) != 2)
  fail(paste("skip-assembly: both candidates must pass through, got", nrow(r5$cands_out), "rows"))
if (any(!is.na(r5$cands_out$rescued_from)))
  fail("skip-assembly: rescued_from must be NA for every row")
ok("skip-assembly -> candidates pass through unchanged, rescue columns NA (D-10 DoS guard)")

# =========================================================================
# Subtest 6: forced-pass — below_threshold candidate rescued, status forced
# =========================================================================
# Single below_threshold candidate with no own subtype support. The best
# alternative (2b) passes floors -> rescue fires and status is forced to "pass"
# regardless of the original confirmation_status. The genotype-diverse
# finalization re-ranks the lone survivor to rank 1 (dominant-first), so the
# rescued FASTA basename is _cand1.
r6 <- run_rescue(
  "forced", "FORCED",
  cands = mk_cands("FORCED",
    mk_cand(2, "1a_M62321", "1a", "1", 20000, 40, "below_threshold")),
  support = mk_support("FORCED",
    mk_support_row("2b", "2b_AY232748", 4120, 96.2, 3840, 4.1)),
  refs = REFS
)
assert_schema("forced", r6$cands_out)
if (nrow(r6$cands_out) != 1)
  fail(paste("forced-pass: exactly one candidate must survive, got", nrow(r6$cands_out)))
f6 <- r6$cands_out %>% filter(candidate_rank == 1)
if (nrow(f6) != 1)
  fail("forced-pass: the lone survivor must be re-ranked to rank 1")
if (is.na(f6$rescued_from[1]))
  fail("forced-pass: the below_threshold candidate must be rescued")
if (f6$confirmation_status[1] != "pass")
  fail(paste("forced-pass: confirmation_status must be forced to 'pass', got",
             f6$confirmation_status[1]))
if (!any(grepl("_cand1\\.", r6$fastas)))
  fail(paste("forced-pass: the rescued FASTA basename must match '_cand1.'; got:",
             paste(r6$fastas, collapse = ",")))
ok("forced-pass -> below_threshold row rescued, status forced pass, re-ranked to _cand1. (D-05/D-06)")

# =========================================================================
# Subtest 7: stale pass-through FASTA removed on ref-changing rescue
# =========================================================================
r7 <- run_rescue(
  "stale", "STALE",
  cands = mk_cands("STALE",
    mk_cand(1, "1a_M62321", "1a", "1", 80000, 95, "pass")),
  support = mk_support("STALE",
    mk_support_row("2b", "2b_AY232748", 4120, 96.2, 3840, 4.1)),
  refs = REFS,
  prestage = "STALE.1a_M62321_cand1.fa"
)
assert_schema("stale", r7$cands_out)
s7 <- r7$cands_out %>% filter(candidate_rank == 1)
if (is.na(s7$rescued_from[1]))
  fail("stale-passthrough: rank-1 candidate must be rescued to trigger the dedup check")
cand1_fastas <- grep("_cand1\\.", r7$fastas, value = TRUE)
if (any(grepl("1a_M62321_cand1\\.", cand1_fastas)))
  fail(paste("stale-passthrough: stale pass-through must be removed; got:",
             paste(cand1_fastas, collapse = ",")))
if (!any(grepl("2b_AY232748_cand1\\.", cand1_fastas)))
  fail(paste("stale-passthrough: rescued FASTA must remain; got:",
             paste(cand1_fastas, collapse = ",")))
if (length(cand1_fastas) != 1)
  fail(paste("stale-passthrough: exactly ONE _cand1. FASTA must survive; got:",
             paste(cand1_fastas, collapse = ",")))
ok("stale-passthrough -> ref-changing rescue removes the {orig_ref} FASTA (dup-@SQ regression)")

# =========================================================================
# Subtest 8: collapse guard — rescue blocked when rescue_ref is the original
#            ref of another candidate slot (ERR1810469 regression)
# =========================================================================
REFS_COLLAPSE <- c(REFS, list("3a_D17763" = strrep("T", 60)))

r8 <- run_rescue(
  "collapse", "COLLAPSE",
  cands = mk_cands("COLLAPSE",
    mk_cand(1, "3a_D17763", "3a", "3", 5009, 46, "pass"),
    mk_cand(2, "1a_M62321",  "1a", "1", 1069, 69, "pass")),
  support = mk_support("COLLAPSE",
    mk_support_row("1a", "1a_M62321", 4503, 93.1, 4200, 5.4)),
  refs = REFS_COLLAPSE
)
assert_schema("collapse", r8$cands_out)
c8_cand1 <- r8$cands_out %>% filter(candidate_rank == 1)
c8_cand2 <- r8$cands_out %>% filter(candidate_rank == 2)
if (!is.na(c8_cand1$rescued_from[1]))
  fail(paste("collapse-guard: cand1 must NOT be rescued when rescue_ref is cand2's ref, got rescued_from",
             c8_cand1$rescued_from[1]))
if (c8_cand1$candidate_ref[1] != "3a_D17763")
  fail(paste("collapse-guard: cand1 must remain 3a_D17763, got", c8_cand1$candidate_ref[1]))
if (c8_cand2$candidate_ref[1] != "1a_M62321")
  fail(paste("collapse-guard: cand2 must remain 1a_M62321, got", c8_cand2$candidate_ref[1]))
if (nrow(r8$cands_out) != 2)
  fail(paste("collapse-guard: both candidates must survive, got", nrow(r8$cands_out), "rows"))
ok("collapse-guard -> rescue blocked when rescue_ref already held by another slot (ERR1810469 regression)")

# =========================================================================
# Subtest 9: own-subtype-confirmed — no rescue when candidate's own subtype
#            has strong assembly support (floors pass)
# =========================================================================
# cand1 = 2b with strong 2b support (floors pass). Even though there is also
# 1a support, the own-subtype guard fires first and no rescue is attempted.
r9 <- run_rescue(
  "ownconfirm", "OWNCONFIRM",
  cands = mk_cands("OWNCONFIRM",
    mk_cand(1, "2b_AY232748", "2b", "2", 80000, 95, "pass")),
  support = mk_support("OWNCONFIRM",
    mk_support_row("2b", "2b_AY232748", 4500, 96.0, 4300, 5.0),
    mk_support_row("1a", "1a_M62321",   5000, 94.0, 4800, 3.0)),
  refs = REFS
)
assert_schema("ownconfirm", r9$cands_out)
o9 <- r9$cands_out %>% filter(candidate_rank == 1)
if (!is.na(o9$rescued_from[1]))
  fail(paste("own-subtype-confirmed: must NOT rescue when own subtype passes floors, got rescued_from",
             o9$rescued_from[1]))
if (o9$candidate_ref[1] != "2b_AY232748")
  fail(paste("own-subtype-confirmed: candidate_ref must remain 2b_AY232748, got", o9$candidate_ref[1]))
ok("own-subtype-confirmed -> no rescue when candidate's own subtype has strong assembly support")

# =========================================================================
# Subtest 10: same-genotype collapse — a redundant same-genotype slot is
#             dropped entirely (not kept-but-unrescued) at selection.
# =========================================================================
# cand1 = 1a_EF407457 (pass, genotype 1, dominant by reads). cand2 = 1c_AY651061
# (pass, genotype 1) is the same genotype as the dominant and is NOT the permitted
# 1a/1b cross-pair, so the genotype-collapse step removes it. Only the dominant
# 1a survives — a within-genotype "minor" the pipeline cannot resolve never
# reaches JOINT_MAPPING.
r10 <- run_rescue(
  "dupgeno", "DUPGENO",
  cands = mk_cands("DUPGENO",
    mk_cand(1, "1a_EF407457", "1a", "1", 10000, 88, "pass"),
    mk_cand(2, "1c_AY651061", "1c", "1",   200, 12, "pass")),
  support = mk_support("DUPGENO",
    mk_support_row("1a", "1a_AF009606", 4000, 95.0, 3500, 3.0)),
  refs = REFS_DUPGENO
)
assert_schema("dupgeno", r10$cands_out)
if (nrow(r10$cands_out) != 1)
  fail(paste("same-genotype-collapse: only the dominant genotype-1 ref must survive, got",
             nrow(r10$cands_out), "rows"))
if (r10$cands_out$candidate_ref[1] != "1a_EF407457")
  fail(paste("same-genotype-collapse: survivor must be the dominant 1a_EF407457, got",
             r10$cands_out$candidate_ref[1]))
if (any(grepl("1c_AY651061", r10$fastas)))
  fail(paste("same-genotype-collapse: the collapsed 1c FASTA must be removed; got:",
             paste(r10$fastas, collapse = ",")))
ok("same-genotype-collapse -> redundant same-genotype slot dropped at selection (within-genotype unresolvable)")

# =========================================================================
# Subtest 11: 1a1b-exception (negative control) — rescue to 1b fires even
#             though cand1 is 1a (1a/1b is the permitted co-infection pair).
# =========================================================================
# cand1 = 1a_M62321 (genotype 1, subtype 1a). cand2 = 1c_AY651061 (genotype 1,
# subtype 1c). Assembly support has only 1b passing floors with best_ref=1b_D90208.
# cand1 own-subtype (1a) has no support → rescue to 1b attempted but blocked by
# same-genotype guard (cand2=1c is not 1a/1b). cand2 own-subtype (1c) has no
# support → rescue to 1b_D90208; same-genotype guard checks: rescue_sub=1b,
# other_sub=1a → 1a/1b exception applies → rescue fires.
r11 <- run_rescue(
  "ab_exception", "ABEXCEP",
  cands = mk_cands("ABEXCEP",
    mk_cand(1, "1a_M62321",  "1a", "1", 10000, 88, "pass"),
    mk_cand(2, "1c_AY651061","1c", "1",   200, 12, "pass")),
  support = mk_support("ABEXCEP",
    mk_support_row("1b", "1b_D90208", 5200, 99.0, 5100, 8.0)),
  refs = REFS_DUPGENO
)
assert_schema("ab_exception", r11$cands_out)
e11_cand2 <- r11$cands_out %>% filter(candidate_rank == 2)
if (is.na(e11_cand2$rescued_from[1]))
  fail("1a1b-exception: cand2 must be rescued to 1b_D90208 (1a/1b cross-subtype pair allowed)")
if (e11_cand2$candidate_ref[1] != "1b_D90208")
  fail(paste("1a1b-exception: cand2 candidate_ref must be 1b_D90208, got",
             e11_cand2$candidate_ref[1]))
ok("1a1b-exception -> rescue to 1b fires when other slot is 1a (permitted cross-subtype pair)")

# =========================================================================
# Subtest 12: de-novo nomination (ERR1810475 shape) — a genuine different-
#             genotype strain absent from first-mapping is APPENDED as a new
#             candidate from assembly_support.
# =========================================================================
# Single genotype-2 candidate (2b, well-supported own subtype -> no rescue). The
# assembly support ALSO carries a strong genotype-1 contig (1a, >3kb, passes the
# four floors) that first-mapping never ranked. Genotype 1 is absent from the
# candidate set, so it is nominated and appended at rank 2 with status forced to
# "pass" and a _cand2. FASTA written.
REFS_NOMINATE <- c(REFS, list("1a_HQ850279" = strrep("G", 60)))
r12 <- run_rescue(
  "nominate", "NOMINATE",
  cands = mk_cands("NOMINATE",
    mk_cand(1, "2b_AY232748", "2b", "2", 4460744, 100, "pass")),
  support = mk_support("NOMINATE",
    mk_support_row("2b", "2b_AY232748", 9549, 91.2, 9471, 30307.0),
    mk_support_row("1a", "1a_HQ850279", 4734, 93.1, 4738, 7.1)),
  refs = REFS_NOMINATE
)
assert_schema("nominate", r12$cands_out)
if (nrow(r12$cands_out) != 2)
  fail(paste("nominate: the nominated 1a must be appended (2 candidates), got",
             nrow(r12$cands_out), "rows"))
n12_dom <- r12$cands_out %>% filter(candidate_rank == 1)
n12_min <- r12$cands_out %>% filter(candidate_rank == 2)
if (n12_dom$candidate_ref[1] != "2b_AY232748")
  fail(paste("nominate: dominant must remain 2b_AY232748, got", n12_dom$candidate_ref[1]))
if (n12_min$candidate_ref[1] != "1a_HQ850279")
  fail(paste("nominate: the nominated minor must be 1a_HQ850279, got", n12_min$candidate_ref[1]))
if (n12_min$candidate_genotype[1] != "1")
  fail(paste("nominate: nominated minor genotype must be 1, got", n12_min$candidate_genotype[1]))
if (n12_min$confirmation_status[1] != "pass")
  fail("nominate: nominated minor confirmation_status must be forced to 'pass'")
if (is.na(n12_min$rescue_trigger[1]) || !grepl("denovo-nomination", n12_min$rescue_trigger[1]))
  fail("nominate: nominated minor must carry a 'denovo-nomination' rescue_trigger")
if (!any(grepl("1a_HQ850279_cand2\\.", r12$fastas)))
  fail(paste("nominate: a 1a_HQ850279_cand2. FASTA must be written; got:",
             paste(r12$fastas, collapse = ",")))
ok("denovo-nomination -> different-genotype de-novo strain appended as a mapped candidate (ERR1810475)")

# =========================================================================
# Subtest 13: same-genotype collapse, divergent gt-4 (ERR1810450/451 shape) —
#             two same-genotype references collapse to a single best ref.
# =========================================================================
# cand1 = 4a_DQ418789 (dominant), cand2 = a generic genotype-4 reference (subtype
# "4"). Both are genotype 4 and not a 1a/1b pair, so the collapse keeps only the
# best read recruiter (4a) and drops the generic genotype-4 second reference, so a
# single genotype-4 reference enters JOINT_MAPPING.
REFS_GT4 <- c(REFS, list("4a_DQ418789" = strrep("A", 60), "4_JF735135" = strrep("C", 60)))
r13 <- run_rescue(
  "gt4collapse", "GT4",
  cands = mk_cands("GT4",
    mk_cand(1, "4a_DQ418789", "4a", "4", 80000, 89, "pass"),
    mk_cand(2, "4_JF735135",  "4",  "4", 30000, 77, "pass")),
  support = mk_support("GT4",
    mk_support_row("4", "4_JF735135", 9506, 85.6, 9417, 130.7)),
  refs = REFS_GT4
)
assert_schema("gt4collapse", r13$cands_out)
if (nrow(r13$cands_out) != 1)
  fail(paste("gt4-collapse: only one genotype-4 reference must survive, got",
             nrow(r13$cands_out), "rows"))
if (r13$cands_out$candidate_ref[1] != "4a_DQ418789")
  fail(paste("gt4-collapse: survivor must be the dominant 4a_DQ418789, got",
             r13$cands_out$candidate_ref[1]))
if (any(grepl("4_JF735135_cand", r13$fastas)))
  fail(paste("gt4-collapse: the dropped generic genotype-4 FASTA must be removed; got:",
             paste(r13$fastas, collapse = ",")))
ok("gt4-collapse -> two same-genotype refs collapse to the single best (ERR1810450/451)")

cat("\nALL PASS\n")
