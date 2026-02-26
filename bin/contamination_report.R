#!/usr/bin/env Rscript
# contamination_report.R
# Parses all-vs-all BLAST output (outfmt 6 + custom columns) from the
# contamination check subworkflow, identifies cross-sample contig pairs,
# and writes a TSV, a heatmap PNG, and a MultiQC-compatible JSON table.
#
# Usage: contamination_report.R <blast_txt> <prefix>

library(tidyverse)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: contamination_report.R <blast_txt> <prefix>")
}

blast_file <- args[1]
prefix     <- args[2]

# Column names matching the BLAST -outfmt used in the subworkflow:
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
col_names <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

blast <- read_tsv(
    blast_file,
    col_names = col_names,
    col_types  = "ccdiiiiiiidid",
    comment    = "#"
)

# Parse sample IDs from headers.
# Header format written by CONTIG_FILTER_RENAME: sampleID|original_contig_header
blast <- blast %>%
    mutate(
        q_sample = sub("\\|.*", "", qseqid),
        s_sample = sub("\\|.*", "", sseqid),
        q_contig = sub("^[^|]*\\|", "", qseqid),
        s_contig = sub("^[^|]*\\|", "", sseqid)
    )

# Keep only cross-sample hits (different samples) and remove self-hits
# (identical query and subject), then deduplicate reciprocal pairs.
cross <- blast %>%
    filter(q_sample != s_sample, qseqid != sseqid) %>%
    mutate(
        snp_distance = mismatch + gapopen,
        # Canonical key: sort both sides alphabetically so A-vs-B == B-vs-A
        pair_key = paste(
            pmin(qseqid, sseqid),
            pmax(qseqid, sseqid),
            sep = "|||"
        )
    ) %>%
    distinct(pair_key, .keep_all = TRUE) %>%
    select(
        sample1          = q_sample,
        contig1          = q_contig,
        sample2          = s_sample,
        contig2          = s_contig,
        pident,
        alignment_length = length,
        mismatches       = mismatch,
        gapopen,
        snp_distance,
        evalue,
        bitscore
    ) %>%
    arrange(snp_distance, desc(pident))

write_tsv(cross, paste0(prefix, ".contamination_pairs.tsv"))

# ── Heatmap ────────────────────────────────────────────────────────────────────
all_samples <- sort(unique(c(blast$q_sample, blast$s_sample)))

if (nrow(cross) > 0) {
    # Count shared contig pairs per sample-pair (symmetric)
    heatmap_counts <- bind_rows(
        cross %>% select(sample1, sample2),
        cross %>% select(sample1 = sample2, sample2 = sample1)
    ) %>%
        count(sample1, sample2, name = "n_pairs")
} else {
    heatmap_counts <- tibble(
        sample1  = character(),
        sample2  = character(),
        n_pairs  = integer()
    )
}

# Full grid so every cell is drawn (including diagonal = sample vs itself)
grid <- expand_grid(sample1 = all_samples, sample2 = all_samples) %>%
    left_join(heatmap_counts, by = c("sample1", "sample2")) %>%
    mutate(
        n_pairs = replace_na(n_pairs, 0L),
        # Diagonal cells (same sample) shown as NA to distinguish from 0
        n_pairs = if_else(sample1 == sample2, NA_integer_, n_pairs)
    )

p <- ggplot(grid, aes(x = sample2, y = sample1, fill = n_pairs)) +
    geom_tile(color = "grey70", linewidth = 0.4) +
    geom_text(
        aes(label = if_else(!is.na(n_pairs) & n_pairs > 0,
                            as.character(n_pairs), "")),
        size = 5
    ) +
    scale_fill_gradient(
        low      = "white",
        high     = "#d73027",
        na.value = "grey90",
        name     = "Shared\ncontigs",
        limits   = c(0, NA)
    ) +
    labs(
        title    = "Cross-sample contig sharing (potential contamination)",
        subtitle = paste0(
            "Contigs \u2265 1000 bp, filtered by BLAST local alignment"
        ),
        x = "Sample", y = "Sample"
    ) +
    theme_bw() +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank()
    )

ggsave(
    paste0(prefix, ".contamination_heatmap.png"),
    plot   = p,
    width  = max(6, length(all_samples) * 0.9 + 2),
    height = max(5, length(all_samples) * 0.9 + 1.5),
    dpi    = 150
)

# ── MultiQC custom content ─────────────────────────────────────────────────────
if (nrow(cross) > 0) {
    mqc_data <- cross %>%
        mutate(
            row_id = paste(sample1, "contig", contig1, "vs", sample2, "contig", contig2)
        ) %>%
        select(row_id, sample1, contig1, sample2, contig2,
               pident, alignment_length, mismatches, gapopen, snp_distance,
               evalue, bitscore)

    # Build list of named lists (one per row) keyed by row_id
    mqc_data_list <- setNames(
        lapply(seq_len(nrow(mqc_data)), function(i) {
            as.list(mqc_data[i, -1])   # drop row_id from values
        }),
        mqc_data$row_id
    )
} else {
    mqc_data_list <- setNames(list(), character(0))
}

mqc <- list(
    id           = "contamination_check",
    section_name = "Cross-sample Contamination Check",
    description  = paste0(
        "Contigs \u2265 1000 bp shared between samples, detected by all-vs-all BLAST. ",
        "Each row is a unique cross-sample contig pair. ",
        "SNP distance = mismatches + gap openings. ",
        "Any entries here should be investigated as potential sample contamination."
    ),
    plot_type    = "table",
    pconfig      = list(
        id    = "contamination_table",
        title = "Cross-sample contig pairs"
    ),
    headers      = list(
        sample1          = list(title = "Sample 1"),
        contig1          = list(title = "Contig 1"),
        sample2          = list(title = "Sample 2"),
        contig2          = list(title = "Contig 2"),
        pident           = list(title = "% Identity",    format = "{:.1f}"),
        alignment_length = list(title = "Aln. length"),
        mismatches       = list(title = "Mismatches"),
        gapopen          = list(title = "Gap openings"),
        snp_distance     = list(title = "SNP distance"),
        evalue           = list(title = "E-value",       format = "{:.2e}"),
        bitscore         = list(title = "Bitscore",      format = "{:.1f}")
    ),
    data         = mqc_data_list
)

write(
    toJSON(mqc, auto_unbox = TRUE, pretty = TRUE),
    paste0(prefix, ".contamination_mqc.json")
)
