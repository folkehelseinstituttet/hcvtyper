#!/usr/bin/env Rscript

# genotype_utils.R --------------------------------------------------------
# Canonical, single-sourced 2k1b-aware genotype-from-subtype rule.
#
# This is a SOURCED helper, not an entrypoint: it defines a function and has
# no commandArgs parser. Consumers source() it from their task workdir (the
# file is staged there as a declared `path` process input by the calling
# Nextflow module), then call genotype_from_subtype() on a subtype vector.
#
# The rule is a verbatim extraction of the historical inline idiom in
# summarize_mapping_to_all_references.R:37 — for the 2k1b recombinant we keep
# the whole subtype name as the genotype; for every other subtype the genotype
# is the first character of the subtype (e.g. "3a" -> "3", "1b" -> "1").
# -------------------------------------------------------------------------

# Defensive: the consumer (summarize_mapping_to_all_references.R) already loads
# tidyverse, which provides if_else(). This guard only fires if sourced into a
# session that has not (Assumption A5).
if (!exists("if_else")) {
  library(tidyverse)
}

genotype_from_subtype <- function(subtype) {
  if_else(subtype == "2k1b", subtype, substr(subtype, 1, 1))
}
