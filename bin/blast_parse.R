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
  scaf_dot <- scaf_top %>%
    arrange(desc(kmer_cov), desc(sc_length)) %>%
    mutate(contig_order = factor(qseqid, levels = rev(unique(qseqid)))) # Reverse order so highest coverage is on top

  p_dot <- ggplot(scaf_dot, aes(x = sc_length, y = contig_order)) +
    geom_point(aes(size = kmer_cov, color = subtype), alpha = 0.8) +
    scale_size_continuous(range = c(2, 10)) +
    scale_color_viridis_d(option = "D") +
    labs(
      title = paste0(prefix, ": Contig length vs coverage"),
      x = "Contig length (bp)",
      y = "Contigs (ordered by coverage, then length)",
      size = "K-mer coverage",
      color = "Subtype"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 6)
    )

  ggsave(
    paste0(prefix, ".contig_dot_plot.png"),
    plot = p_dot,
    width = 9,
    bg = "white",
    height = max(4, 0.2 * nrow(scaf_dot)),  # scale height with number of contigs
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
# --- 7. Major / minor reference summary + FASTA export ---------------------
if (nrow(scaf_top) > 0) {
# a) pick closest major and (optionally) minor reference names
major_name <- scaf_top$sseqid[1]                 # best overall hit
major_geno <- str_sub(major_name, 1, 1)
major_contig <- scaf_top %>%
  slice_max(sc_length, n = 1) %>%               # longest contig for this reference
  select(qseqid, sc_length) %>% distinct() %>% # Remove duplicates if several hits against the same reference
  pull(qseqid)

minor_vec  <- scaf_top %>%
  filter(!str_starts(subtype, major_geno)) %>%   # must be different genotype
  slice_head(n = 1) %>%
  pull(sseqid)
minor_name <- if (length(minor_vec) == 0) NA_character_ else minor_vec
} else {
  major_name <- NA_character_
  major_contig <- NA_character_
  minor_name <- NA_character_
}

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
  major_contig_length = scaf %>% filter(sseqid == major_name) %>%
                   slice_max(sc_length, n = 1) %>%
                   # If the major contig have multiple blast hits against the same reference, the length will be duplicated
                   select(qseqid, sc_length) %>% distinct() %>% pull(sc_length),
  minor_ref    = minor_name,
  minor_contig_length = if (is.na(minor_name)) NA_integer_ else
                   scaf %>% filter(sseqid == minor_name)  %>%
                    filter(qseqid != major_contig) %>%  # Exclude the major contig if it is also a minor hit
                    slice_max(sc_length, n = 1) %>%
                    # If the minor contig have multiple blast hits against the same reference, the length will be duplicated
                    select(qseqid, sc_length) %>% distinct() %>% pull(sc_length)
)
write_csv(summary_tbl, paste0(prefix, ".blastparse.csv"))

