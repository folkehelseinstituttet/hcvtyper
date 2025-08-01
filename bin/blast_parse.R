#!/usr/bin/env Rscript
#
# blast_parse.R  —  tidy BLAST‑tab output, basic QC plots,
#                   per‑subtype scaffold FASTAs (≥500 bp),
#                   and an “alignment” bar‑plot of top hits.
#
# Usage: blast_parse.R <prefix> <blast_out> <contigs> <references> <agens>
#        * <references> and <agens> are kept for CLI compatibility
#          but no longer used by this script.
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
contigs  <- args[3]
references <- args[4]
# agens      <- args[5]   # not used

ref_fa     <- read.fasta(             # DNA FASTA with HCV references
  file    = references,
  seqtype = "DNA"
)

## ── 2. Input files ----------------------------------------------------------
# Contigs FASTA (for sequence export)
contigs_fa <- read.fasta(file = contigs, seqtype = "DNA")

# BLAST outfmt 6 table
scaf <- read_tsv(
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

# Write reformatted BLAST output
write_csv(scaf, paste0(prefix, "_blast_out.csv"))

## ── 3. Quick QC plots -------------------------------------------------------
# 3a. Top‑30 bitscores
scaf %>%
  arrange(desc(bitscore)) %>% slice_head(n = 30) %>%
  ggplot(aes(x = reorder(sseqid, -bitscore), y = bitscore)) +
  geom_point() +
  labs(
    title = paste0(prefix, " – top 30 BLAST bitscores"),
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
    title = paste0(prefix, " – top 30 BLAST hit lengths"),
    x = "Reference",
    y = "Hit length (bp)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(
  paste0(prefix, ".hitlength_plot.png"),
  dpi = 300, width = 9, height = 4, bg = "white"
)

## ── 4. Top BLAST hit per scaffold (all lengths) ----------------------------
scaf_top <- scaf %>%
  arrange(evalue, desc(bitscore)) %>%      # best hit = lowest e‑value, highest bitscore
  group_by(qseqid) %>% slice(1) %>% ungroup()

write_csv(scaf_top, paste0(prefix, "_top_hits.csv"))

## ── 5. Alignment‑style bar plot (ALL contigs) ----------------------------
# Create scaffold factor levels sorted by subtype, then by sstart
scaf_ordered <- scaf_top %>%
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
    y = "Scaffold",
    fill = "Subtype"
  )

ggsave(paste0(prefix, ".alignment_plot.png"),
       plot = p_align,
       width = 10,
       bg = "white",
       height = max(4, 0.2 * nrow(scaf_ordered)),  # scale with number of contigs
       dpi = 300)

## ── 6. Scaffold FASTAs ≥500 bp, grouped by subtype -------------------------
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

# --- 7. Major / minor reference summary + FASTA export ---------------------

# a) pick closest major and (optionally) minor reference names
major_name <- scaf_top$sseqid[1]                 # best overall hit
major_geno <- str_sub(major_name, 1, 1)

minor_vec  <- scaf_top %>%
  filter(!str_starts(subtype, major_geno)) %>%   # must be different genotype
  slice_head(n = 1) %>%
  pull(sseqid)
minor_name <- if (length(minor_vec) == 0) NA_character_ else minor_vec

# b) FASTA export -----------------------------------------------------------
# Helper that writes the sequence only if it exists
write_ref_fasta <- function(ref_name, tag) {
  if (!is.na(ref_name) && ref_name %in% names(ref_fa)) {
    write.fasta(
      sequences = ref_fa[ref_name],
      names     = ref_name,
      file.out  = paste0(prefix, ".", ref_name, "_", tag, ".fa")
    )
  }
}

write_ref_fasta(major_name, "major")
write_ref_fasta(minor_name, "minor")

# c) summary CSV
summary_tbl <- tibble(
  sample       = prefix,
  major_ref    = major_name,
  major_length = scaf %>% filter(sseqid == major_name) %>%
                   slice_max(sc_length, n = 1) %>% pull(sc_length),
  minor_ref    = minor_name,
  minor_length = if (is.na(minor_name)) NA_integer_ else
                   scaf %>% filter(sseqid == minor_name) %>%
                     slice_max(sc_length, n = 1) %>% pull(sc_length)
)
write_csv(summary_tbl, paste0(prefix, ".blastparse.csv"))

