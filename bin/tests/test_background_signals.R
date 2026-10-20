#!/usr/bin/env Rscript

# test_background_signals.R -----------------------------------------------
# Unit tests for bin/background_signals.R (quick task 261020-q6v): two
# sample-level signals that a co-infection minor may be background reads.
#
#   BS-1  candidate_read_fractions(): share of the sample's candidate reads,
#         summing to 1 per sample; NA when a sample has no reads
#   BS-2  genotype_dup_rates() on real data (ERR1810447, Thomson 2016, public
#         SRA): the sample's own 1b genome is ~91 % PCR duplicates, the 2b genome
#         carried over from a run-mate ~2 % — pooled by each contig's top hit
#   BS-3  candidate_dup_signal(): TRUE for that 2b; FALSE when concordant; NA
#         (never FALSE) when the test has no power — unsaturated library, too few
#         fragments, no contigs — or the candidate is the dominant genotype
#   BS-4  background_note(): sentence carries the numbers and the duplicate
#         verdict; a concordant low-share minor (IVT 1:50) never reads as
#         background; NA when nothing fires
#   BS-5  empty inputs give typed-empty frames
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
fx_dir   <- file.path(this_dir, "fixtures", "background")

source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "background_signals.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
ok <- function(msg) cat("PASS:", msg, "\n")
near <- function(a, b, tol = 1e-6) isTRUE(abs(a - b) < tol)

# ── BS-1 ──────────────────────────────────────────────────────────────────────
reads <- tribble(
  ~sampleName, ~candidate_rank, ~reads_withdup,
  "S1",        1L,              457069,
  "S1",        2L,              5677,
  "S2",        1L,              0,
  "S2",        2L,              0
)
fr <- candidate_read_fractions(reads)
s1 <- fr %>% filter(sampleName == "S1")
if (!near(sum(s1$candidate_read_fraction), 1)) fail("BS-1 fractions must sum to 1")
if (!near(s1$candidate_read_fraction[s1$candidate_rank == 2], 5677 / (457069 + 5677)))
  fail("BS-1 minor fraction wrong")
if (!all(is.na(fr$candidate_read_fraction[fr$sampleName == "S2"]))) fail("BS-1 zero-read sample must be NA")
ok(sprintf("BS-1 read fractions (ERR1810447-like minor = %.4f)", s1$candidate_read_fraction[2]))

# ── BS-2 ──────────────────────────────────────────────────────────────────────
s <- "ERR1810447"
blast_out <- read_csv(file.path(fx_dir, paste0(s, "_blast_out.csv")), show_col_types = FALSE,
                      col_types = cols(qseqid = col_character(), sseqid = col_character(),
                                       subtype = col_character(), .default = col_double())) %>%
  mutate(sampleName = s)
dups_dir <- tempfile("bs_"); dir.create(dups_dir)
file.copy(file.path(fx_dir, paste0(s, ".contig_dups.tsv")), dups_dir)
contig_dups <- read_contig_dups(dups_dir)
if (nrow(contig_dups) < 10 || !all(contig_dups$sampleName == s)) fail("BS-2 read_contig_dups did not read the fixture")
gd <- genotype_dup_rates(contig_dups, blast_out)
g1 <- gd %>% filter(genotype == "1"); g2 <- gd %>% filter(genotype == "2")
if (nrow(g1) != 1 || nrow(g2) != 1) fail("BS-2 expected genotype 1 and 2 pools")
if (!(g1$dup_rate > 0.85 && g2$dup_rate < 0.1))
  fail(sprintf("BS-2 expected 1 ~0.91 and 2 ~0.02, got %.3f / %.3f", g1$dup_rate, g2$dup_rate))
ok(sprintf("BS-2 real-data pools: genotype 1 %.3f (%g frags), genotype 2 %.3f (%g frags)",
           g1$dup_rate, g1$dup_fragments, g2$dup_rate, g2$dup_fragments))

# ── BS-3 ──────────────────────────────────────────────────────────────────────
cands <- tibble(sampleName = s, candidate_rank = 1:2, candidate_ref = c("1b_EU781827", "2b_D10988"))
dom   <- tibble(sampleName = s, dominant_ref = "1b_EU781827")
sig <- candidate_dup_signal(cands, dom, gd)
if (!isTRUE(sig$dup_discordant[sig$candidate_rank == 2])) fail("BS-3 the carried-over 2b must be discordant")
if (!is.na(sig$dup_discordant[sig$candidate_rank == 1])) fail("BS-3 the dominant itself must be NA, not FALSE")
ok(sprintf("BS-3 ERR1810447 2b discordant (ratio %.3f)", sig$dup_ratio[sig$candidate_rank == 2]))

synth <- tribble(
  ~sampleName, ~genotype, ~dup_fragments, ~dup_duplicates,
  "C",         "1",       10000,          9000,   # concordant minor
  "C",         "3",       500,            420,
  "U",         "1",       10000,          2500,   # unsaturated library (0.25)
  "U",         "3",       500,            5,
  "F",         "1",       10000,          9000,   # too few minor fragments
  "F",         "3",       12,             0,
  "G",         "1",       10000,          9000    # same genotype minor (1a/1b)
) %>% mutate(dup_rate = dup_duplicates / dup_fragments)
cs <- tibble(sampleName = c("C", "U", "F", "G", "N"), candidate_rank = 2L,
             candidate_ref = c("3a_D17763", "3a_D17763", "3a_D17763", "1a_HQ850279", "3a_D17763"))
ds <- tibble(sampleName = c("C", "U", "F", "G", "N"), dominant_ref = "1b_EU781827")
r <- candidate_dup_signal(cs, ds, synth) %>% arrange(sampleName)
want <- c(C = FALSE, F = NA, G = NA, N = NA, U = NA)
got  <- setNames(r$dup_discordant, r$sampleName)[names(want)]
if (!identical(unname(got), unname(want)))
  fail(paste("BS-3 power/same-genotype handling wrong:", paste(names(got), got, collapse = " ")))
ok("BS-3 concordant -> FALSE; unsaturated / few fragments / same genotype / no contigs -> NA")

# ── BS-4 ──────────────────────────────────────────────────────────────────────
# discordant (+ low share): background wording
n_both <- background_note(2L, "2b_D10988", 0.0123, 0.018, 0.909, TRUE, 0.02)
if (!str_detect(n_both, "Co-infection candidate 2 \\(2b_D10988\\) holds 1.23% of the sample's HCV reads \\(below 2.00%") ||
    !str_detect(n_both, "1.80% PCR duplicates against 90.9% for the dominant genotype, so its reads did not go through this library's PCR") ||
    !str_detect(n_both, "check contamination_check/"))
  fail(paste("BS-4 discordant sentence wrong:", n_both))
# discordant only (share above the floor, e.g. ERR1810469 at 17.5 %)
n_dup <- background_note(2L, "1a_HQ850279", 0.175, 0.027, 0.96, TRUE, 0.02)
if (str_detect(n_dup, "holds") || !str_detect(n_dup, "\\(1a_HQ850279\\): its contigs are 2.70% PCR duplicates") ||
    !str_detect(n_dup, "did not go through")) fail(paste("BS-4 duplicate-only sentence wrong:", n_dup))
# low share, concordant (IVT5: genuine 1:50 at 1.95 %): must NOT read as background
n_conc <- background_note(2L, "2a_AB047639", 0.0195, 0.972, 0.987, FALSE, 0.02)
if (!str_detect(n_conc, "holds 1.95%") || !str_detect(n_conc, "went through this library's PCR") ||
    str_detect(n_conc, "did not go through") || str_detect(n_conc, "consistent with background"))
  fail(paste("BS-4 concordant low-share sentence wrong:", n_conc))
# low share, untested
n_na <- background_note(2L, "3a_D17763", 0.01, NA, 0.3, NA, 0.02)
if (!str_detect(n_na, "the duplicate test was not informative")) fail(paste("BS-4 untested sentence wrong:", n_na))
if (!is.na(background_note(2L, "3a_D17763", 0.3, 0.8, 0.9, FALSE, 0.02))) fail("BS-4 nothing fired must be NA")
if (!is.na(background_note(2L, "3a_D17763", NA, NA, NA, NA, 0.02))) fail("BS-4 all-NA must be NA")
ok("BS-4 sentence variants: discordant, duplicate-only, concordant low share, untested, none")

# ── BS-5 ──────────────────────────────────────────────────────────────────────
if (nrow(candidate_read_fractions(reads[0, ])) != 0) fail("BS-5 empty reads")
if (nrow(read_contig_dups(tempfile())) != 0) fail("BS-5 missing dir")
if (nrow(genotype_dup_rates(contig_dups[0, ], blast_out)) != 0) fail("BS-5 empty dups")
if (nrow(candidate_dup_signal(cands[0, ], dom, gd)) != 0) fail("BS-5 empty candidates")
ok("BS-5 empty inputs are typed-empty")

cat("\nAll background_signals tests passed.\n")
