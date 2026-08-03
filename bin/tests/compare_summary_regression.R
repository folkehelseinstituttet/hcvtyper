#!/usr/bin/env Rscript

# compare_summary_regression.R --------------------------------------------
# Old-vs-new Summary.csv differ, for proving that a change altered ONLY the
# columns it was supposed to alter (quick task 260803-ogc).
#
# Deliberately NOT named test_*.R: bin/tests/run_all.sh globs test_*.R and runs
# each with no arguments, and this script requires two file paths. It is a
# developer/CI tool invoked explicitly, not part of the unit suite.
#
# Usage:
#   Rscript bin/tests/compare_summary_regression.R OLD.csv NEW.csv [allowed_cols]
#
#   allowed_cols  Comma-separated columns permitted to differ. Defaults to the
#                 260803-ogc expectation:
#                   review_flag,call_confidence,Major_genotype,Minor_genotype
#                 Pass "review_flag,call_confidence" to check the flag change
#                 alone (i.e. with the Major_genotype fix reverted or excluded).
#
# Exit status: 0 when every difference is confined to the allowed columns and
# the sample set is unchanged; 1 otherwise. Column ADDITIONS are reported and
# tolerated when the added column is in the allowed set (the Major_genotype fix
# adds those columns to runs that lacked them); any other schema change fails.
#
# Everything is read as character so that a run missing a column, or a column
# whose type readr infers differently between the two files, cannot abort the
# comparison or produce a spurious 1 vs 1.0 difference.
# -------------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compare_summary_regression.R OLD.csv NEW.csv [allowed_cols]\n")
  quit(status = 2)
}
old_path <- args[1]
new_path <- args[2]
allowed <- if (length(args) >= 3 && nzchar(args[3])) {
  trimws(str_split(args[3], ",")[[1]])
} else {
  c("review_flag", "call_confidence", "Major_genotype", "Minor_genotype")
}

for (p in c(old_path, new_path)) if (!file.exists(p)) {
  cat("FAIL: no such file:", p, "\n"); quit(status = 2)
}

read_summary <- function(p) read_csv(p, col_types = cols(.default = col_character()),
                                     progress = FALSE)
old <- read_summary(old_path)
new <- read_summary(new_path)

say  <- function(...) cat(..., "\n", sep = "")
rule <- function(t) say("\n", strrep("-", 70), "\n", t, "\n", strrep("-", 70))
status <- 0
fail <- function(...) { say("FAIL: ", ...); status <<- 1 }

say("old: ", old_path, "  (", nrow(old), " rows x ", ncol(old), " cols)")
say("new: ", new_path, "  (", nrow(new), " rows x ", ncol(new), " cols)")
say("columns permitted to differ: ", paste(allowed, collapse = ", "))

if (!"sampleName" %in% names(old) || !"sampleName" %in% names(new)) {
  cat("FAIL: both files must carry a sampleName column\n"); quit(status = 2)
}

## -- 1. Sample set --------------------------------------------------------
rule("1. Sample set")
only_old <- setdiff(old$sampleName, new$sampleName)
only_new <- setdiff(new$sampleName, old$sampleName)
if (length(only_old)) fail(length(only_old), " sample(s) only in OLD: ",
                           paste(head(only_old, 10), collapse = ", "))
if (length(only_new)) fail(length(only_new), " sample(s) only in NEW: ",
                           paste(head(only_new, 10), collapse = ", "))
if (anyDuplicated(old$sampleName) || anyDuplicated(new$sampleName))
  fail("duplicate sampleName values — comparison would be ambiguous")
if (!length(only_old) && !length(only_new)) say("OK: identical sample set (", nrow(new), ")")

## -- 2. Schema ------------------------------------------------------------
rule("2. Schema")
dropped <- setdiff(names(old), names(new))
added   <- setdiff(names(new), names(old))
if (length(dropped)) fail("column(s) DROPPED: ", paste(dropped, collapse = ", "))
if (length(added)) {
  unexpected <- setdiff(added, allowed)
  if (length(unexpected)) {
    fail("column(s) ADDED that are not in the allowed set: ",
         paste(unexpected, collapse = ", "))
  } else {
    say("OK (allowed): column(s) added: ", paste(added, collapse = ", "))
  }
}
if (!length(dropped) && !length(added)) say("OK: identical column set")

## -- 3. Per-column value differences --------------------------------------
rule("3. Value differences, per column")
shared_cols <- intersect(names(old), names(new))
key_cols    <- intersect(c("sampleName"), shared_cols)

# Align both frames on sampleName so row order cannot masquerade as a difference.
common <- intersect(old$sampleName, new$sampleName)
o <- old %>% filter(sampleName %in% common) %>% arrange(sampleName)
n <- new %>% filter(sampleName %in% common) %>% arrange(sampleName)

# NA == NA is equal; NA vs a value is a difference.
differs <- function(a, b) xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b)

diff_tbl <- map_dfr(setdiff(shared_cols, key_cols), function(cl) {
  d <- differs(o[[cl]], n[[cl]])
  tibble(column = cl, n_differing = sum(d), allowed = cl %in% allowed)
}) %>% filter(n_differing > 0) %>% arrange(desc(n_differing))

if (nrow(diff_tbl) == 0) {
  say("No column differs on any of the ", length(common), " shared samples.")
} else {
  print(as.data.frame(diff_tbl), row.names = FALSE)
  violations <- diff_tbl %>% filter(!allowed)
  if (nrow(violations) > 0) {
    fail(nrow(violations), " column(s) changed that should not have: ",
         paste(violations$column, collapse = ", "))
    rule("Examples of disallowed differences (up to 10 per column)")
    for (cl in violations$column) {
      d <- differs(o[[cl]], n[[cl]])
      ex <- tibble(sampleName = o$sampleName[d],
                   old = o[[cl]][d], new = n[[cl]][d]) %>% head(10)
      say("\n== ", cl, " ==")
      print(as.data.frame(ex), row.names = FALSE)
    }
  }
}

## -- 4. Verdict -----------------------------------------------------------
rule("Verdict")
if (status == 0) {
  say("PASS — every difference is confined to: ", paste(allowed, collapse = ", "))
} else {
  say("FAIL — see above.")
}
quit(status = status)
