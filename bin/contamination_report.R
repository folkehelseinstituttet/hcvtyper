#!/usr/bin/env Rscript
# contamination_report.R
# Parses all-vs-all BLAST output from the contamination check subworkflow,
# identifies cross-sample contig pairs, infers contamination direction from
# SPAdes coverage, optionally flags index-hopping noise using fastp read counts,
# and optionally annotates pairs with GLUE resistance mutation similarity.
#
# Usage: contamination_report.R <blast_txt> <prefix> [hop_rate] [min_dir_ratio] [genome_size] [kmer_size] [min_aln_cov]
#
# min_aln_cov: minimum fraction of the shorter contig that must be covered by
#   the BLAST alignment to retain the hit.  E.g. 0.9 means the alignment must
#   span ≥90% of whichever of the two contigs is shorter.  Hits that are short
#   partial matches are excluded.  Default: 0.9.
#   SPAdes reports k-mer coverage (not read depth) in contig headers.
#   The conversion is: read_depth = kmer_cov / (read_length - kmer_size + 1)
#   For Illumina 150 bp reads the default SPAdes final k is 127, giving a
#   multiplier of 24. This value is used to scale the hop threshold into the
#   same k-mer coverage units before comparison. Default: 127.
#
# Optional inputs auto-detected in working directory:
#   fastp/   *.fastp.json                    - per-sample fastp JSON files
#   glue/    *.major.major.*.nodup.json      - per-sample major GLUE JSON files

library(tidyverse)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: contamination_report.R <blast_txt> <prefix> [hop_rate] [min_dir_ratio] [genome_size]")
}

blast_file    <- args[1]
prefix        <- args[2]
hop_rate      <- if (length(args) >= 3 && nchar(args[3]) > 0) as.numeric(args[3]) else 0.001
min_dir_ratio <- if (length(args) >= 4 && nchar(args[4]) > 0) as.numeric(args[4]) else 10.0
genome_size   <- if (length(args) >= 5 && nchar(args[5]) > 0) as.numeric(args[5]) else 9500
# SPAdes k-mer size: used to convert the read-depth hop threshold into k-mer
# coverage units so that both sides of the flag_index_hopping comparison are
# in the same units. Formula: kmer_cov = read_depth * (read_length - kmer_size + 1)
kmer_size     <- if (length(args) >= 6 && nchar(args[6]) > 0) as.numeric(args[6]) else 127
min_aln_cov   <- if (length(args) >= 7 && nchar(args[7]) > 0) as.numeric(args[7]) else 0.9

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

# Keep only cross-sample hits (different samples) and remove self-hits,
# then deduplicate reciprocal pairs (A-vs-B == B-vs-A).
cross <- blast %>%
    filter(q_sample != s_sample, qseqid != sseqid) %>%
    mutate(
        snp_distance = mismatch + gapopen,
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

# ── Parse SPAdes contig metadata (length + coverage) ─────────────────────────
# SPAdes names contigs: NODE_<n>_length_<bp>_cov_<depth>
if (nrow(cross) > 0) {
    cross <- cross %>%
        mutate(
            contig1_len = as.numeric(sub(".*_length_(\\d+)_cov_.*", "\\1", contig1)),
            contig1_cov = as.numeric(sub(".*_cov_([0-9.]+).*",      "\\1", contig1)),
            contig2_len = as.numeric(sub(".*_length_(\\d+)_cov_.*", "\\1", contig2)),
            contig2_cov = as.numeric(sub(".*_cov_([0-9.]+).*",      "\\1", contig2))
        ) %>%
        mutate(
            # Direction: higher-coverage contig = likely source sample
            source_sample    = if_else(contig1_cov >= contig2_cov, sample1, sample2),
            recipient_sample = if_else(contig1_cov >= contig2_cov, sample2, sample1),
            source_cov       = pmax(contig1_cov, contig2_cov),
            recipient_cov    = pmin(contig1_cov, contig2_cov),
            source_len       = if_else(contig1_cov >= contig2_cov, contig1_len, contig2_len),
            recipient_len    = if_else(contig1_cov >= contig2_cov, contig2_len, contig1_len),
            cov_ratio        = source_cov / pmax(recipient_cov, 0.01)
        ) %>%
        # ── Alignment coverage filter ─────────────────────────────────────────
        # Keep only hits where the alignment spans >= min_aln_cov of the shorter
        # contig. Short partial BLAST hits between unrelated regions are excluded.
        mutate(
            aln_cov_frac = alignment_length / pmin(contig1_len, contig2_len)
        ) %>%
        filter(aln_cov_frac >= min_aln_cov)
    message(sprintf(
        "[filter] alignment coverage >= %.0f%% of shorter contig: %d pairs retained",
        min_aln_cov * 100, nrow(cross)
    ))
}

# ── Optional: index-hopping threshold from fastp JSON files ──────────────────
# Nextflow stages *.fastp.json files into fastp/ subdir when provided.
fastp_files <- list.files("fastp", pattern = "\\.fastp\\.json$", full.names = TRUE)

if (length(fastp_files) > 0) {
    fastp_data <- tibble(path = fastp_files) %>%
        mutate(
            sample_id   = sub("\\.fastp\\.json$", "", basename(path)),
            total_reads = map_dbl(path, function(f) {
                d <- tryCatch(fromJSON(f), error = function(e) NULL)
                if (is.null(d)) NA_real_
                else d[["summary"]][["after_filtering"]][["total_reads"]]
            })
        ) %>%
        filter(!is.na(total_reads))

    n_samples           <- nrow(fastp_data)
    run_reads           <- sum(fastp_data$total_reads)
    # Expected coverage of index-hopped reads reaching a single recipient sample:
    # total run reads * hop_rate = reads that can hop to any well;
    # divide by n_samples for per-well share, convert to coverage with 150bp reads.
    # hop_cov_threshold is in READ DEPTH units:
    #   (total_run_reads * hop_rate * read_length) / (n_samples * genome_size)
    # SPAdes contig coverage is K-MER coverage, not read depth.
    # To compare them directly we scale up to k-mer units:
    #   kmer_hop_threshold = read_depth_threshold * (read_length - kmer_size + 1)
    read_length         <- 150
    kmer_multiplier     <- read_length - kmer_size + 1      # e.g. 24 for k=127
    hop_cov_threshold   <- (run_reads * hop_rate * read_length) / (n_samples * genome_size)
    hop_kmer_threshold  <- hop_cov_threshold * kmer_multiplier
    has_fastp           <- TRUE
    message(sprintf(
        "[fastp] %d samples | %s total reads | index-hop threshold: %.3fx read-depth | %.1fx k-mer (k=%d, multiplier=%d, hop_rate=%.3f%%)",
        n_samples, format(run_reads, big.mark = ","),
        hop_cov_threshold, hop_kmer_threshold, kmer_size, kmer_multiplier, hop_rate * 100
    ))
} else {
    fastp_data          <- tibble(sample_id = character(), total_reads = numeric())
    hop_cov_threshold   <- NA_real_
    hop_kmer_threshold  <- NA_real_
    kmer_multiplier     <- NA_real_
    n_samples           <- NA_integer_
    run_reads           <- NA_real_
    has_fastp           <- FALSE
}

# Flag pairs where recipient contig k-mer coverage is below the index-hopping
# threshold. The threshold is stored in read-depth units (human-readable) but
# the comparison is done in k-mer units to match the SPAdes contig coverage.
# contig1_cov / contig2_cov are SPAdes k-mer coverage values.
if (nrow(cross) > 0) {
    cross <- cross %>%
        mutate(
            index_hop_threshold = if (has_fastp) round(hop_cov_threshold, 4) else NA_real_,
            flag_index_hopping  = if (has_fastp)
                                      (contig1_cov < hop_kmer_threshold |
                                       contig2_cov < hop_kmer_threshold)
                                  else NA
        )
}

# ── Optional: GLUE resistance mutation similarity ─────────────────────────────
# Nextflow stages *.major.major.*.nodup.json files into glue/ subdir when provided.
glue_files <- list.files("glue",
                         pattern = "_major\\.major\\.nodup\\.json$",
                         full.names = TRUE)

# Load a GLUE JSON file, skipping the debug-log preamble that precedes the JSON
load_glue_json <- function(f) {
    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
    start <- which(startsWith(trimws(lines), "{"))[1]
    if (is.na(start)) return(NULL)
    tryCatch(
        fromJSON(paste(lines[start:length(lines)], collapse = "\n"),
                 simplifyVector = FALSE),
        error = function(e) NULL
    )
}

# Return character vector of present, sufficiently-covered RAS mutation names
extract_ras <- function(f) {
    d <- load_glue_json(f)
    if (is.null(d)) return(character(0))
    ras <- tryCatch(
        d[["phdrReport"]][["samReferenceResult"]][["rasScanResults"]],
        error = function(e) list()
    )
    if (length(ras) == 0) return(character(0))
    present <- Filter(
        function(r) isTRUE(r[["present"]]) && isTRUE(r[["sufficientCoverage"]]),
        ras
    )
    if (length(present) == 0) return(character(0))
    vapply(present, function(r) as.character(r[["variationName"]]), character(1))
}

# Filename pattern: {sampleID}.{subtype}_{ref}_major.major.nodup.json
# e.g. 2665137-HCV.1b_EU781827_major.major.nodup.json  ->  2665137-HCV
get_sample_from_glue <- function(f) {
    sub("\\.[^.]+_[^.]+_major\\.major\\.nodup\\.json$", "", basename(f))
}

jaccard_sim <- function(a, b) {
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) == 0 && length(b) == 0) return(NA_real_)
    length(intersect(a, b)) / length(union(a, b))
}

if (length(glue_files) > 0) {
    glue_ras <- list()
    for (f in glue_files) {
        sid <- get_sample_from_glue(f)
        if (!nzchar(sid) || sid == basename(f)) next  # regex didn't match
        ras <- extract_ras(f)
        glue_ras[[sid]] <- union(glue_ras[[sid]], ras)
    }
    has_glue <- length(glue_ras) > 0
    if (has_glue) message(sprintf("[glue] RAS parsed for %d samples", length(glue_ras)))
} else {
    glue_ras <- list()
    has_glue <- FALSE
}

# Add GLUE RAS Jaccard similarity and shared mutation columns to cross-sample pairs
if (nrow(cross) > 0 && has_glue) {
    cross <- cross %>%
        rowwise() %>%
        mutate(
            .r1 = list(if (!is.null(glue_ras[[sample1]])) glue_ras[[sample1]] else character(0)),
            .r2 = list(if (!is.null(glue_ras[[sample2]])) glue_ras[[sample2]] else character(0)),
            shared_mutations   = paste(sort(intersect(.r1, .r2)), collapse = "; "),
            n_shared_mutations = length(intersect(.r1, .r2)),
            jaccard_similarity = jaccard_sim(.r1, .r2)
        ) %>%
        ungroup() %>%
        select(-.r1, -.r2)
}

write_tsv(cross, paste0(prefix, ".contamination_pairs.tsv"))

# ── Directional contamination summary ─────────────────────────────────────────
all_samples <- sort(unique(c(blast$q_sample, blast$s_sample)))

if (nrow(cross) > 0) {
    dir_summary <- cross %>%
        group_by(source_sample, recipient_sample) %>%
        summarise(
            n_contig_pairs           = n(),
            n_above_hop_threshold    = if (has_fastp)
                                           sum(!flag_index_hopping, na.rm = TRUE)
                                       else NA_integer_,
            median_cov_ratio         = median(cov_ratio),
            median_source_cov        = median(source_cov),
            median_recipient_cov     = median(recipient_cov),
            median_jaccard_ras       = if (has_glue)
                                           median(jaccard_similarity, na.rm = TRUE)
                                       else NA_real_,
            shared_mutations_union   = if (has_glue) {
                                           all_mut <- unique(unlist(strsplit(
                                               shared_mutations[nchar(shared_mutations) > 0], "; "
                                           )))
                                           paste(sort(all_mut), collapse = "; ")
                                       } else NA_character_,
            .groups = "drop"
        ) %>%
        arrange(desc(n_contig_pairs), desc(median_cov_ratio))
} else {
    dir_summary <- tibble(
        source_sample = character(), recipient_sample = character(),
        n_contig_pairs = integer(), n_above_hop_threshold = integer(),
        median_cov_ratio = numeric(), median_source_cov = numeric(),
        median_recipient_cov = numeric(), median_jaccard_ras = numeric(),
        shared_mutations_union = character()
    )
}

# Remove columns that are entirely NA (i.e., optional inputs not provided).
# Guard: on a zero-row tibble all(is.na(.)) == TRUE for every column (vacuous
# truth), so we skip the select to preserve the schema.
if (nrow(dir_summary) > 0) {
    dir_summary <- dir_summary %>%
        select(where(~ !all(is.na(.))))
}

write_tsv(dir_summary, paste0(prefix, ".contamination_direction.tsv"))

# ── Directional bubble plot ────────────────────────────────────────────────────
# Only show high-confidence pairs: median_cov_ratio >= min_dir_ratio.
# Low cov-ratio pairs represent shared HCV sequence biology, not contamination.
dir_plot <- dir_summary %>%
    filter(median_cov_ratio >= min_dir_ratio)

if (nrow(dir_plot) > 0) {
    # Build per-bubble label: total pairs (pairs above index-hop threshold)
    if (has_fastp && "n_above_hop_threshold" %in% names(dir_plot)) {
        dir_plot <- dir_plot %>%
            mutate(bubble_label = paste0(n_contig_pairs, " (", n_above_hop_threshold, ")"))
    } else {
        dir_plot <- dir_plot %>%
            mutate(bubble_label = as.character(n_contig_pairs))
    }

    # Colour by GLUE Jaccard if available, otherwise by coverage ratio
    use_glue_colour <- has_glue && "median_jaccard_ras" %in% names(dir_plot)

    p_dir <- ggplot(dir_plot,
                    aes(x = source_sample, y = recipient_sample,
                        size = n_contig_pairs)) +
        {
            if (use_glue_colour)
                aes(color = median_jaccard_ras)
            else
                aes(color = log10(pmax(median_cov_ratio, 1) + 1))
        }

    p_dir <- ggplot(dir_plot,
                    aes(x     = source_sample,
                        y     = recipient_sample,
                        size  = n_contig_pairs,
                        color = if (use_glue_colour) median_jaccard_ras
                                else log10(pmax(median_cov_ratio, 1) + 1))) +
        geom_point(alpha = 0.85) +
        geom_text(aes(label = bubble_label),
                  size = 3.5, color = "black", vjust = -1.1) +
        scale_size_area(max_size = 14, name = "Shared\ncontig pairs") +
        {
            if (use_glue_colour)
                scale_color_gradient(low = "#fee090", high = "#a50026",
                                     na.value = "grey80",
                                     name = "Mutation\nJaccard")
            else
                scale_color_gradient(low = "#fee090", high = "#d73027",
                                     name = expression(log[10](cov~ratio)))
        } +
        labs(
            title    = "Inferred contamination direction (high-confidence pairs)",
            subtitle = paste0(
                sprintf("Only pairs with median coverage ratio \u2265 %.0fx shown. ", min_dir_ratio),
                if (has_fastp)
                    sprintf("Index-hop threshold = %.3fx (%.1f%% hop rate, %s reads).\n",
                            hop_cov_threshold, hop_rate * 100,
                            format(run_reads, big.mark = ","))
                else "\n",
                "X-axis = likely source (higher SPAdes coverage). ",
                "Y-axis = likely recipient (lower coverage). ",
                if (has_fastp) "Label = total (above hop threshold)." else "",
                if (use_glue_colour) "\nColour = Jaccard similarity of GLUE RAS mutations." else ""
            ),
            x = "Likely source sample",
            y = "Likely recipient sample"
        ) +
        theme_bw() +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank()
        )

    n_src <- length(unique(dir_plot$source_sample))
    n_rec <- length(unique(dir_plot$recipient_sample))
    ggsave(
        paste0(prefix, ".contamination_direction.png"),
        plot   = p_dir,
        width  = max(6, n_src * 0.9 + 3),
        height = max(5, n_rec * 0.9 + 3),
        dpi    = 150
    )
}

# ── Heatmap (symmetric, all pairs) ────────────────────────────────────────────
if (nrow(cross) > 0) {
    heatmap_counts <- bind_rows(
        cross %>% select(sample1, sample2),
        cross %>% select(sample1 = sample2, sample2 = sample1)
    ) %>%
        count(sample1, sample2, name = "n_pairs")
} else {
    heatmap_counts <- tibble(sample1 = character(), sample2 = character(), n_pairs = integer())
}

grid <- expand_grid(sample1 = all_samples, sample2 = all_samples) %>%
    left_join(heatmap_counts, by = c("sample1", "sample2")) %>%
    mutate(
        n_pairs = replace_na(n_pairs, 0L),
        n_pairs = if_else(sample1 == sample2, NA_integer_, n_pairs)
    )

p_heatmap <- ggplot(grid, aes(x = sample2, y = sample1, fill = n_pairs)) +
    geom_tile(color = "grey70", linewidth = 0.4) +
    geom_text(
        aes(label = if_else(!is.na(n_pairs) & n_pairs > 0, as.character(n_pairs), "")),
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
        subtitle = "Contigs \u2265 1000 bp, filtered by BLAST local alignment",
        x = "Sample", y = "Sample"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank()
    )

ggsave(
    paste0(prefix, ".contamination_heatmap.png"),
    plot   = p_heatmap,
    width  = max(6, length(all_samples) * 0.9 + 2),
    height = max(5, length(all_samples) * 0.9 + 1.5),
    dpi    = 150
)

# ── MultiQC custom content ─────────────────────────────────────────────────────
mqc_cols <- c("sample1", "contig1", "sample2", "contig2",
              "pident", "alignment_length", "aln_cov_frac", "mismatches", "gapopen", "snp_distance",
              "evalue", "bitscore")
if (has_fastp && "flag_index_hopping" %in% names(cross)) {
    mqc_cols <- c(mqc_cols, "flag_index_hopping", "index_hop_threshold")
}
if (has_glue && "n_shared_mutations" %in% names(cross)) {
    mqc_cols <- c(mqc_cols, "n_shared_mutations", "jaccard_similarity", "shared_mutations")
}

if (nrow(cross) > 0) {
    mqc_data <- cross %>%
        mutate(row_id = paste(sample1, "contig", contig1, "vs", sample2, "contig", contig2)) %>%
        select(row_id, any_of(mqc_cols))

    mqc_data_list <- setNames(
        lapply(seq_len(nrow(mqc_data)), function(i) as.list(mqc_data[i, -1])),
        mqc_data$row_id
    )
} else {
    mqc_data_list <- setNames(list(), character(0))
}

hop_desc <- if (has_fastp) {
    sprintf(" Index-hop cov threshold = %.3fx (rate %.1f%%, %s total run reads).",
            hop_cov_threshold, hop_rate * 100, format(run_reads, big.mark = ","))
} else ""

glue_desc <- if (has_glue) " Jaccard similarity of GLUE RAS mutations shown where available." else ""

mqc_headers <- list(
    sample1          = list(title = "Sample 1"),
    contig1          = list(title = "Contig 1"),
    sample2          = list(title = "Sample 2"),
    contig2          = list(title = "Contig 2"),
    pident           = list(title = "% Identity",    format = "{:.1f}"),
    alignment_length = list(title = "Aln. length"),
    aln_cov_frac     = list(title = "Aln. cov. (frac)", format = "{:.3f}"),
    mismatches       = list(title = "Mismatches"),
    gapopen          = list(title = "Gap openings"),
    snp_distance     = list(title = "SNP distance"),
    evalue           = list(title = "E-value",       format = "{:.2e}"),
    bitscore         = list(title = "Bitscore")
)
if (has_fastp) {
    mqc_headers$flag_index_hopping  <- list(title = "Index-hop flag")
    mqc_headers$index_hop_threshold <- list(title = "Hop threshold (x read-depth)")
}
if (has_glue) {
    mqc_headers$n_shared_mutations <- list(title = "Shared RAS (n)")
    mqc_headers$jaccard_similarity <- list(title = "RAS Jaccard", format = "{:.3f}")
    mqc_headers$shared_mutations   <- list(title = "Shared RAS mutations")
}

mqc <- list(
    id           = "contamination_check",
    section_name = "Cross-sample Contamination Check",
    description  = paste0(
        "Contigs \u2265 1000 bp shared between samples, detected by all-vs-all BLAST. ",
        "Each row is a unique cross-sample contig pair. ",
        "SNP distance = mismatches + gap openings.",
        hop_desc, glue_desc,
        " Any entries should be investigated as potential sample contamination."
    ),
    plot_type = "table",
    pconfig   = list(id = "contamination_table", title = "Cross-sample contig pairs"),
    headers   = mqc_headers,
    data      = mqc_data_list
)

write(
    toJSON(mqc, auto_unbox = TRUE, pretty = TRUE),
    paste0(prefix, ".contamination_mqc.json")
)
