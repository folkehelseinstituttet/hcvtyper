#!/usr/bin/env Rscript
#
# blast_parse.R  —  tidy BLAST‑tab output, basic QC plots,
#                   per‑subtype scaffold FASTAs (≥500 bp),
#                   and an “alignment” bar‑plot of top hits.
#
# Usage: blast_parse.R <prefix> <blast_out> <contigs> <references> <agens>
#        * <references> IS used (read and consumed by write_ref_fasta).
#        * <agens> is retained for CLI compatibility but no longer used.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)   # readr, dplyr, tidyr, ggplot2, purrr
  library(seqinr)      # FASTA I/O
})

## ── 1. Command‑line args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: blast_parse.R <prefix> <blast_out> <contigs> <references> <agens>",
    call. = FALSE
  )
}
prefix     <- args[1]
blast_out  <- args[2]
contigs    <- args[3]
references <- args[4]
# agens      <- args[5]   # not used

if (!file.exists(references) || file.size(references) == 0) {
  stop("Reference FASTA file '", references, "' does not exist or is empty.", call. = FALSE)
}
ref_fa     <- read.fasta(             # DNA FASTA with HCV references
  file    = references,
  seqtype = "DNA"
)

## ── 2. Input files ----------------------------------------------------------
# Contigs FASTA (for sequence export)
if (!file.exists(contigs) || file.size(contigs) == 0) {
  stop("Contigs FASTA file '", contigs, "' does not exist or is empty.", call. = FALSE)
}
contigs_fa <- read.fasta(file = contigs, seqtype = "DNA")

# BLAST outfmt 6 table
# Set up the empty tibble first, in case the blast_out is empty
empty_scaf <- tibble(
  qseqid   = character(0),
  sseqid   = character(0),
  pident   = double(0),
  length   = integer(0),
  mismatch = integer(0),
  gapopen  = integer(0),
  qstart   = integer(0),
  qend     = integer(0),
  sstart   = integer(0),
  send     = integer(0),
  evalue   = double(0),
  bitscore = double(0),
  subtype  = character(0),
  sc_length = double(0),
  kmer_cov  = double(0)
)

# Read BLAST outfmt 6 table
if (!file.exists(blast_out) || file.size(blast_out) == 0) {
  message("BLAST output '", blast_out, "' is missing or empty — continuing with no hits.")
  scaf <- empty_scaf # Set scaf to empty_scaf if blast_out is empty
} else {
  scaf <- tryCatch(
    {
      read_tsv(
        blast_out,
        col_names = FALSE,
        show_col_types = FALSE
      ) %>%
        rename(qseqid  = X1,  sseqid  = X2,  pident   = X3,  length   = X4,
               mismatch = X5, gapopen = X6,  qstart   = X7,  qend     = X8,
               sstart   = X9, send    = X10, evalue   = X11, bitscore = X12) %>%
        # pull subtype from reference header (e.g. 3a_D1776 → subtype = "3a")
        separate(sseqid, into = c("subtype", NA), remove = FALSE) %>%
        # extract scaffold length & kmer coverage from header: NODE_?_length_<len>_cov_<cov>
        mutate(
          sc_length = as.numeric(str_extract(qseqid, "(?<=_length_)[0-9]+")),
          kmer_cov  = as.numeric(str_extract(qseqid, "(?<=_cov_)[0-9.]+"))
        )
    },
    error = function(e) {
      warning("Failed to parse BLAST output '", blast_out, "': ", conditionMessage(e))
      empty_scaf
    }
  )
}

# Write reformatted BLAST output
write_csv(scaf, paste0(prefix, "_blast_out.csv"))

## ── 3. Quick QC plots -------------------------------------------------------
# Make empty plots if no blast_out
if (nrow(scaf) == 0) {
  # placeholder plots
  p_blank <- ggplot() + theme_void() +
    ggtitle(paste0(prefix, " — No BLAST hits found"))
  ggsave(paste0(prefix, ".bitscore_plot.png"), plot = p_blank, dpi = 300, width = 9, height = 4, bg = "white")
  ggsave(paste0(prefix, ".hitlength_plot.png"), plot = p_blank, dpi = 300, width = 9, height = 4, bg = "white")
  ggsave(paste0(prefix, ".alignment_plot.png"), plot = p_blank, dpi = 300, width = 10, height = 4, bg = "white")

  # empty top hits
  write_csv(tibble(), paste0(prefix, "_top_hits.csv"))
} else {
# 3a. Top‑30 bitscores
scaf %>%
  arrange(desc(bitscore)) %>% slice_head(n = 30) %>%
  ggplot(aes(x = reorder(sseqid, -bitscore), y = bitscore)) +
  geom_point() +
  labs(
    title = paste0(prefix, " - top 30 BLAST bitscores"),
    x = "Reference",
    y = "Bitscore"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
  paste0(prefix, ".bitscore_plot.png"),
  dpi = 300, width = 9, height = 4, bg = "white"
)

# 3b. Top‑30 hit lengths
scaf %>%
  arrange(desc(length)) %>% slice_head(n = 30) %>%
  ggplot(aes(x = reorder(sseqid, -length), y = length)) +
  geom_point() +
  labs(
    title = paste0(prefix, " - top 30 BLAST hit lengths"),
    x = "Reference",
    y = "Hit length (bp)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
  paste0(prefix, ".hitlength_plot.png"),
  dpi = 300, width = 9, height = 4, bg = "white"
)
}

## ── 4. Top BLAST hit per scaffold (all lengths) ----------------------------
scaf_top <- if (nrow(scaf) > 0) {
scaf %>%
  arrange(evalue, desc(bitscore)) %>%      # best hit = lowest e‑value, highest bitscore
  group_by(qseqid) %>% slice(1) %>% ungroup() %>% # take first hit per scaffold
  arrange(desc(bitscore)) # Arrange again by bitscore, because the order was unset after the previous line
} else {
  tibble()
}
write_csv(scaf_top, paste0(prefix, "_top_hits.csv"))

## ── 4b. Neutral per‑subtype assembly‑support roll‑up (ASUP‑01, D‑01/D‑02/D‑03) ----
# Dominance‑neutral replacement for the §7 major/minor logic: for every subtype
# seen in the de novo contigs, summarise the SINGLE best contig by sc_length
# (D‑03) and carry THAT contig's four ASUP‑01 metrics — full contig length,
# BLAST % identity, BLAST alignment length, and k‑mer coverage. Raw metrics ONLY:
# no denovo_min_* threshold floor is applied (D‑01 — the substantiality verdict
# is Phase 8). The raw subtype token is carried; genotype derivation is deferred
# to summarize.R (D‑02). §6/§7 below stay UNCHANGED (D‑04 legacy shim).
if (nrow(scaf_top) > 0) {
  support_tbl <- scaf_top %>%
    group_by(subtype) %>%
    # single best contig per subtype, by full contig length (D‑03). distinct() on
    # qseqid/sc_length is load‑bearing: one contig can have several BLAST hits to
    # the same reference and would otherwise duplicate the winning row.
    slice_max(sc_length, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(subtype, sseqid, qseqid, sc_length, pident, length, kmer_cov) %>%
    distinct() %>%
    transmute(
      sample                = prefix,
      subtype,
      best_ref              = sseqid,
      best_contig_length    = sc_length,
      best_contig_pident    = pident,
      best_contig_aln_length = length,
      best_contig_kmer_cov  = kmer_cov
    )
} else {
  # T‑07‑01 DoS guard: zero hits / skip‑assembly → typed header‑only CSV, exit 0,
  # never abort (mirror the §6/§7 empty guards and the line‑110 empty‑write idiom).
  support_tbl <- tibble(
    sample                 = character(0),
    subtype                = character(0),
    best_ref               = character(0),
    best_contig_length     = double(0),
    best_contig_pident     = double(0),
    # WR-03: double (not integer) to match the populated path (`length` from
    # read_tsv) and the join helper's typed-empty support frame, so the same
    # logical column has ONE consistent type everywhere.
    best_contig_aln_length = double(0),
    best_contig_kmer_cov   = double(0)
  )
}
write_csv(support_tbl, paste0(prefix, ".assembly_support.csv"))

## ── 5. Alignment‑style bar plot (100 top hit contigs) ----------------------------
# Create scaffold factor levels sorted by subtype, then by sstart
if (nrow(scaf_top) > 0) {
scaf_ordered <- scaf_top %>%
  slice_head(n = 100) %>% # Only use the 100 top hits
  arrange(subtype, sstart, qseqid) %>%
  mutate(y_pos = row_number())  # numeric y position

# Plot: alignment-like overview of scaffold BLAST hits
p_align <- scaf_ordered %>%
  ggplot(aes(xmin = pmin(sstart, send),
             xmax = pmax(sstart, send),
             ymin = y_pos - 0.4,
             ymax = y_pos + 0.4,
             fill = subtype)) +
  geom_rect() +
  scale_y_continuous(
    breaks = scaf_ordered$y_pos,
    labels = scaf_ordered$qseqid
  ) +
  scale_fill_viridis_d(option = "D") +
  theme_minimal() +
  labs(
    title = paste0(prefix, ": Blast hit regions (sorted by subtype)"),
    x = "Reference position",
    y = "Contig",
    fill = "Subtype"
  )

ggsave(paste0(prefix, ".alignment_plot.png"),
       plot = p_align,
       width = 10,
       bg = "white",
       height = max(4, 0.2 * nrow(scaf_ordered)),  # scale with number of contigs
       dpi = 300)
}

## ── 5b. Contig length vs coverage dot plot ----------------------------
if (nrow(scaf_top) > 0) {
  # Limit to the same 100 top contigs as in the alignment plot
  scaf_dot <- scaf_top %>%
    arrange(desc(bitscore)) %>%
    slice_head(n = 100) %>%  # limit to top 100 contigs
    arrange(desc(kmer_cov), desc(sc_length)) %>%
    mutate(
      contig_key   = paste(subtype, qseqid, sep = "::"),
      contig_order = factor(contig_key, levels = rev(unique(contig_key))) # highest coverage on top within subtype
    )

  p_dot <- ggplot(scaf_dot, aes(x = sc_length, y = contig_order)) +
    geom_point(aes(size = kmer_cov, color = subtype), alpha = 0.8) +
    scale_y_discrete(labels = function(x) sub(".*::", "", x)) +  # show only contig ID on axis
    scale_size_continuous(range = c(2, 10)) +
    scale_color_viridis_d(option = "D") +
    labs(
      title = paste0(prefix, ": Contig length vs coverage by subtype"),
      x = "Contig length (bp)",
      y = "Contigs (within subtype: by coverage, then length)",
      size = "K-mer coverage",
      color = "Subtype"
    ) +
    theme_minimal() +
    theme(
      axis.text.y   = element_text(size = 6),
      strip.text.y  = element_text(face = "bold")
    ) +
    facet_grid(rows = vars(subtype), scales = "free_y", space = "free_y")

  ggsave(
    paste0(prefix, ".contig_dot_plot.png"),
    plot = p_dot,
    width = 9,
    bg = "white",
    height = max(4, 0.25 * n_distinct(scaf_dot$subtype) + 0.18 * nrow(scaf_dot)),  # scale with content
    dpi = 300
  )
}

## ── 6. Contig FASTAs ≥500 bp, grouped by subtype -------------------------
if (nrow(scaf_top) > 0) {
scaf_top_long <- scaf_top %>% filter(sc_length >= 500)

# Write one FASTA per subtype
scaf_top_long %>%
  group_by(subtype) %>%
  group_walk(~{
    subtype_name <- .y$subtype
    seqs <- contigs_fa[.x$qseqid]
    write.fasta(
      sequences = seqs,
      names     = names(seqs),
      file.out  = paste0(prefix, ".", subtype_name, "_contigs.fa")
    )
  })
}
# --- 7. Major / minor reference summary (display-only, consumed by summarize.R) ---
if (nrow(scaf_top) > 0) {
# a) pick closest major and (optionally) minor reference names
# 260805 (§9.4): the major slot now uses the same ONE ROW discipline as the minor
# slot below. scaf_top row 1 IS the row that defines major_name (the frame is sorted
# by descending bitscore), so reading the reference, the contig and the contig length
# off that single row guarantees all three describe the same contig.
#
# Previously the three were derived three different ways:
#   major_name          = scaf_top$sseqid[1]                        -- row 1
#   major_contig        = longest contig in scaf_top, UNFILTERED    -- any contig
#   major_contig_length = longest contig in the FULL scaf table
#                         filtered to sseqid == major_name          -- any contig
# On sim1 that reported major_ref = 1a_HQ850279 (whose own contig is 9,076 bp)
# alongside major_contig_length = 9,339 bp — the length of the 1b contig, which
# merely carries a secondary 78.7%-identity hit against 1a_HQ850279.
major_row    <- scaf_top %>% slice(1)
major_name   <- major_row$sseqid[1]              # best overall hit
major_geno   <- str_sub(major_name, 1, 1)
major_contig <- major_row$qseqid[1]
major_len    <- major_row$sc_length[1]

# The minor selection is captured as ONE ROW, and the reference, the contig name and
# the contig length are all read off that row (260803-ogc). Previously only sseqid was
# pulled here and the contig name was re-derived downstream in summarize.R from the
# full BLAST table as "the best-bitscore hit to this reference" — a DIFFERENT grain.
# scaf_top holds one row per contig (its own top hit), so this filter keeps contigs
# whose OWN top hit is off-genotype; the downstream re-derivation searched all contigs
# unrestricted and could therefore return a contig that this filter had excluded.
#
# SampleSA-1a is the case in the wild: NODE_3 (1620 bp, top hit 6i_DQ835770,
# bitscore 97) wins here, but NODE_2 — a 1a contig whose 5'UTR/core region hits the
# same 6i reference at bitscore 1074 — won the downstream lookup. Summary.csv reported
# denovo_minor_ref = 6i_DQ835770 and denovo_minor_contig_length = 1620 (both NODE_3)
# next to denovo_minor_contig = NODE_2_length_3232. Three fields, two contigs, and an
# analyst sent to the wrong sequence.
#
# That particular re-derivation could only bite the MINOR slot, because it searched
# by BITSCORE and major_name is the globally best hit — its row is necessarily also
# the top-bitscore row for that reference.
#
# 260805 (§9.4): the major slot had the same disease from a different vector. Its
# length was derived by slice_max(sc_length) — by LENGTH, not bitscore — so the
# argument above never protected it, and a longer contig carrying a weak secondary
# hit to major_name won. Fixed above by reading the major slot off ONE ROW too.
minor_row  <- scaf_top %>%
  filter(!str_starts(subtype, major_geno)) %>%   # must be different genotype
  slice_head(n = 1)
minor_name   <- if (nrow(minor_row) == 0) NA_character_ else minor_row$sseqid[1]
minor_contig <- if (nrow(minor_row) == 0) NA_character_ else minor_row$qseqid[1]
minor_len    <- if (nrow(minor_row) == 0) NA_real_      else minor_row$sc_length[1]
} else {
  major_name <- NA_character_
  major_contig <- NA_character_
  major_len <- NA_real_
  minor_name <- NA_character_
  minor_contig <- NA_character_
  minor_len <- NA_real_
}

# b) summary CSV
summary_tbl <- tibble(
  sample       = prefix,
  major_ref    = major_name,
  # Read off major_row (260805, §9.4) — same contig as major_ref, by construction.
  major_contig_length = major_len,
  minor_ref    = minor_name,
  # Both read straight off minor_row, so minor_ref / minor_contig /
  # minor_contig_length always describe ONE contig (260803-ogc).
  #
  # This REPLACES a lookup that re-queried scaf_top for the longest contig sharing
  # minor_ref as its top hit, excluding major_contig. Where exactly one contig has
  # that reference as its top hit — the overwhelming majority — the value is
  # unchanged. Where several do, the reported length is now the contig that actually
  # won the selection (highest bitscore) rather than the longest of the group, and
  # the old form could also yield a zero-length pull when the winner happened to be
  # major_contig.
  minor_contig = minor_contig,
  minor_contig_length = minor_len
)
write_csv(summary_tbl, paste0(prefix, ".blastparse.csv"))

