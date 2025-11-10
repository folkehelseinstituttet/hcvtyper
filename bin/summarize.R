#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# Define variables --------------------------------------------------------

samplesheet      <- args[1]
stringency_1     <- args[2]
stringency_2     <- args[3]
pipeline_version <- args[4]
pipeline_name    <- args[5]

script_name_version <- if (!is.na(pipeline_version) && nzchar(trimws(pipeline_version))) {
  paste(pipeline_name, pipeline_version)
} else {
  paste(pipeline_name, "(version unknown)")
}

path_1 <- "trimmed/"
path_2 <- "kraken_classified/"
path_3 <- "parsefirst_mapping/"
path_4 <- "stats_withdup/"
path_5 <- "stats_markdup/"
path_6 <- "depth/"
path_7 <- "blast/"
path_8 <- "glue/"
path_9 <- "id/"
path_10 <- "variation/"



# Trimmed ----------------------------------------------------------------

# List files
trimmed_files <- list.files(path = path_1, pattern = ".log$", full.names = TRUE)

# Empty df
trimmed_df <- tibble(
  sampleName = rep(NA_character_, length(trimmed_files)),
  total_raw_reads = rep(NA_real_, length(trimmed_files)),
  total_trimmed_reads = rep(NA_real_, length(trimmed_files))
)

for (i in seq_along(trimmed_files)) {
  f <- trimmed_files[i]
  # Get the sampleName
  trimmed_df$sampleName[i] <- str_split(basename(f), "\\.")[[1]][1]

  if (!file.exists(f) || file.size(f) == 0) {
    warning(glue::glue("Log file missing or empty: {f}"))
    next
  }

  try(rm(trimmed_stats), silent = TRUE)

  # Read the log file from either cutadapt or fastp
  trimmed_stats <- read_lines(f) %>% tibble(value = .)

  # Check if the string fastp is found in the log file. Then process accordingly
  if (any(str_detect(trimmed_stats$value, "fastp v"))) {
    ## ---- fastp logs ----

    # Raw reads = Read1 before + Read2 before
    raw1_idx <- which(str_detect(trimmed_stats$value, "Read1 before filtering:")) # Finds the line number
    raw2_idx <- which(str_detect(trimmed_stats$value, "Read2 before filtering:"))

    raw1_val <- if (length(raw1_idx) > 0) {
      as.numeric(str_extract(trimmed_stats$value[raw1_idx + 1], "\\d+")) # Pull out the number on the line after
    } else NA
    raw2_val <- if (length(raw2_idx) > 0) {
      as.numeric(str_extract(trimmed_stats$value[raw2_idx + 1], "\\d+"))
    } else NA

    if (!is.na(raw1_val) && !is.na(raw2_val)) {
      trimmed_df$total_raw_reads[i] <- raw1_val + raw2_val
    }

    # Trimmed reads = "reads passed filter:"
    passed <- trimmed_stats %>%
      filter(str_detect(value, "^reads passed filter:"))
    if (nrow(passed) > 0) {
      tmp <- as.numeric(str_extract(passed$value, "\\d+"))
      trimmed_df$total_trimmed_reads[i] <- tmp
    }

  } else {
    ## ---- Cutadapt logs ----
    raw <- filter(trimmed_stats, str_detect(value, "Total read pairs processed:"))
    if (nrow(raw) > 0) {
      tmp <- str_split(raw$value, "\\s+")[[1]] # Split on white space and get the list content
      tmp <- as.numeric(str_remove_all(tmp[length(tmp)], ",")) # extract the last element which contains the pair number, remove commas and create a numeric
      tmp <- tmp * 2 # Double to get reads
      trimmed_df$total_raw_reads[i] <- tmp
    }

    trimmed <- filter(trimmed_stats, str_detect(value, "Pairs written \\(passing filters\\)"))
    if (nrow(trimmed) > 0) {
      tmp <- str_split(trimmed$value, "\\s+")[[1]]
      tmp <- as.numeric(str_remove_all(tmp[length(tmp)-1], ","))
      tmp <- tmp * 2
      trimmed_df$total_trimmed_reads[i] <- tmp
    }
  }
}

trimmed_df <- as_tibble(trimmed_df)

# Kraken ------------------------------------------------------------------
# List files
kraken_files <- list.files(path = path_2, pattern = "kraken2.report.txt$", full.names = TRUE)

# Empty df
kraken_df <- as.data.frame(matrix(nrow = length(kraken_files), ncol = 2))
colnames(kraken_df) <- c("sampleName", "total_classified_reads")

for (i in 1:length(kraken_files)) {
  try(rm(kraken_stats))
  # Get sample name
  kraken_df$sampleName[i] <- str_split(basename(kraken_files[i]), "\\.")[[1]][1]

  # Get total of trimmed sequences put in to the mapping. Adding try() if no root sequences
  try(kraken_stats <- read_tsv(kraken_files[i], col_names = FALSE) %>% filter(X6 == "root") %>% pull(X2))
  if (exists("kraken_stats") & length(kraken_stats) > 0) {
    kraken_df$total_classified_reads[i] <- kraken_stats*2 # Kraken reports read pairs (fragments)
  }
}
kraken_df <- as_tibble(kraken_df)

# Total mapped reads to all references, with duplicates -------------------
first_mapping_files <- list.files(path = path_3, pattern = "parsefirstmapping.csv$", full.names = TRUE)

# Empty df
parsefirstmapping_df <- tibble(
  sampleName = rep(NA_character_, length(first_mapping_files)),
  total_mapped_reads = rep(NA_real_, length(first_mapping_files)),
  major_mapped_reads = rep(NA_real_, length(first_mapping_files)),
  minor_mapped_reads = rep(NA_real_, length(first_mapping_files))
)

# If the length of parsefirstmapping_files is non-zero
if (length(first_mapping_files) > 0) {
  for (i in 1:length(first_mapping_files)) {
    try(rm(sample_parsefirstmapping))
    # Get sample name
    parsefirstmapping_df$sampleName[i] <- str_split(basename(first_mapping_files[i]), "\\.")[[1]][1]

    # Read the mapping stats
    sample_parsefirstmapping <- read_csv(first_mapping_files[i])

    # Get number of mapped reads before duplicate removal
    parsefirstmapping_df$total_mapped_reads[i] <- sample_parsefirstmapping %>% pull(total_mapped_reads)

    # Get the number of mapped reads against all major references belonging to the major subtype
    parsefirstmapping_df$major_mapped_reads[i] <- sample_parsefirstmapping %>% pull(major_reads)

    # Get the number of mapped reads against all minor references belonging to the minor subtype
    parsefirstmapping_df$minor_mapped_reads[i] <- sample_parsefirstmapping %>% pull(minor_reads)
  }
}

parsefirstmapping_df <- as_tibble(parsefirstmapping_df) %>%
  # Calculate the median of total mapped reads accross all samples, and then for each sample the fraction of mapped reads compared to the median
  # In addition calculate the fraction of mapped reads against major and minor references of the total mapped reads
  mutate(
    median_mapped = median(total_mapped_reads),
    fraction_mapped_reads_vs_median = total_mapped_reads / median_mapped,
    percent_mapped_reads_major_firstmapping = round(major_mapped_reads / total_mapped_reads * 100, digits = 2),
    percent_mapped_reads_minor_firstmapping = round(minor_mapped_reads / total_mapped_reads * 100, digits = 2)
  ) %>%
  select(sampleName, total_mapped_reads, fraction_mapped_reads_vs_median, percent_mapped_reads_major_firstmapping, percent_mapped_reads_minor_firstmapping)

# Second mapping, reads mapped with duplicates ----------------------------
# List files
stats_files <- list.files(path = path_4, pattern = "\\withdup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "first_major_minor", "trimmed_reads_withdups_mapped")

for (i in 1:length(stats_files)) {
  try(rm(map_stats))
  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][2]

  # Get major or minor
  tmp_df$first_major_minor[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][3]

  # Read the mapping stats
  map_stats <- read_tsv(stats_files[i], col_names = FALSE, comment = "#")

  # Get number of mapped reads before duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)

  mapped_reads <- as.numeric(mapped_reads)
  tmp_df$trimmed_reads_withdups_mapped[i] <- mapped_reads
}
tmp_df <- as_tibble(tmp_df)

# Add number of raw and trimmed reads - needed for calculation of percentages
tmp_df <- left_join(tmp_df, trimmed_df, by = "sampleName")

# Add number of classified reads from Kraken2 -  - needed for calculation of percentages
tmp_df <- left_join(tmp_df, kraken_df, by = "sampleName")

df_with_dups <- tmp_df %>%
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_withdup_mapped_major = case_when(first_major_minor == "major" ~ trimmed_reads_withdups_mapped)) %>%
  mutate(Reads_withdup_mapped_minor = case_when(first_major_minor == "minor" ~ trimmed_reads_withdups_mapped)) %>%
  # Don't include number of reads mapped in the first mapping. Info must be taken from another process if we should include
  #mutate(Reads_withdup_mapped_first_mapping = case_when(first_major_minor == "first_mapping" ~ trimmed_reads_withdups_mapped)) %>%
  select(-trimmed_reads_withdups_mapped) %>%
  # Create columns for the major and minor references
  mutate(Major_reference = case_when(first_major_minor == "major" ~ reference)) %>%
  mutate(Minor_reference = case_when(first_major_minor == "minor" ~ reference)) %>%
  mutate(Major_reference = str_remove(Major_reference, "_major"),
         Minor_reference = str_remove(Minor_reference, "_minor")) %>%
  select(-reference) %>%
  # Calculate percent of the trimmed reads mapped
  #mutate(total_trimmed_reads_with_dups = as.integer(total_trimmed_reads_with_dups),
  #       Reads_withdup_mapped_major = as.integer(Reads_withdup_mapped_major),
  #       Reads_withdup_mapped_minor = as.integer(Reads_withdup_mapped_minor)) %>%
         #Reads_withdup_mapped_first_mapping = as.integer(Reads_withdup_mapped_first_mapping)) %>%
  mutate(Percent_reads_mapped_of_trimmed_with_dups_major = Reads_withdup_mapped_major / total_trimmed_reads * 100,
         Percent_reads_mapped_of_trimmed_with_dups_minor = Reads_withdup_mapped_minor / total_trimmed_reads * 100) %>%
         #Percent_reads_mapped_with_dups_first_mapping = Reads_withdup_mapped_first_mapping / total_trimmed_reads_with_dups * 100) %>%
  # Create one row per sample
  select(-first_major_minor) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Reads mapped no duplicates ----------------------------------------------
# List files
stats_files <- list.files(path = path_5, pattern = "nodup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "first_major_minor", "trimmed_reads_nodups_mapped")

for (i in 1:length(stats_files)) {
  try(rm(mapped_reads))

  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][2]

  # Get major or minor
  tmp_df$first_major_minor[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][3]

  # Read the mapping stats
  map_stats <- read_tsv(stats_files[i], col_names = FALSE, comment = "#")

  # Get number of mapped reads after duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)

  mapped_reads <- as.numeric(mapped_reads)
  tmp_df$trimmed_reads_nodups_mapped[i] <- mapped_reads

}
tmp_df <- as_tibble(tmp_df)

df_nodups <- tmp_df %>%
  # Create columns for major and minor
  separate(reference, into = c("genotype", NA), sep = "_", remove = F) %>%
  mutate(Major_genotype_mapping = case_when(first_major_minor == "major" ~ genotype)) %>%
  mutate(Minor_genotype_mapping = case_when(first_major_minor == "minor" ~ genotype)) %>%
  select(-genotype) %>%
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_nodup_mapped_major = case_when(first_major_minor == "major" ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_minor = case_when(first_major_minor == "minor" ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_first_mapping = case_when(first_major_minor == "first_mapping" ~ trimmed_reads_nodups_mapped)) %>%
  select(-trimmed_reads_nodups_mapped) %>%
  #mutate(Percent_mapped_major = case_when(first_major_minor == "major" ~ Percent_trimmed_reads_mapped)) %>%
  #mutate(Percent_mapped_minor = case_when(first_major_minor == "minor" ~ Percent_trimmed_reads_mapped)) %>%
  # Create columns for the major and minor references
  mutate(Major_reference = case_when(first_major_minor == "major" ~ reference)) %>%
  mutate(Minor_reference = case_when(first_major_minor == "minor" ~ reference)) %>%
  mutate(Major_reference = str_remove(Major_reference, "_major"),
         Minor_reference = str_remove(Minor_reference, "_minor")) %>%
  # Create one row per sample
  select(-reference, -first_major_minor) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Combine mapped reads data
df_mapped_reads <- full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Remove columns for total_raw_reads, total_trimmed_reads and total_classified_reads.
  # These will be added later to ensure info is kept for samples that were filtered out before the second mapping
  select(-total_raw_reads, -total_trimmed_reads, -total_classified_reads)

# Coverage ----------------------------------------------------------------

# Add both breadth (in percent) and depth (average depth)
# All this is without duplicates

# List files
cov_files <- list.files(path = path_6, pattern = "tsv$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(cov_files), ncol = 7))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_1", "cov_breadth_min_5", "cov_breadth_min_10", "first_major_minor", "avg_depth")

for (i in 1:length(cov_files)) {
  try(rm(cov))

  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][1]

  # Get reference name
  tmp_df$reference[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][2]

  # Get major or minor
  tmp_df$first_major_minor[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][3]

  # Read the depth per position
  cov <- read_tsv(cov_files[i], col_names = FALSE)

  # Reference length
  ref_length <- nrow(cov)

  # Average depth
  tmp_df$avg_depth[i] <- mean(cov$X3)

  # Nr. of positions with coverage >=1, >= 5 and > 9
  # If ref_length is zero it means that no reads were mapped. Set coverage to zero.
  # Coverage may also be zero if there are reads mapped, but never more than 5 per position
  if (ref_length > 0) {
    pos_1 <- nrow(
      cov %>%
        filter(X3 >= 1)
    )

    pos_5 <- nrow(
      cov %>%
        filter(X3 >= 5)
    )

    pos_10 <- nrow(
      cov %>%
        filter(X3 > 9)
    )
    # Coverage breadth
    breadth_1 <- round(pos_1 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_1[i] <- breadth_1

    breadth_5 <- round(pos_5 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_5[i] <- breadth_5

    breadth_10 <- round(pos_10 / ref_length * 100, digits = 2)
    tmp_df$cov_breadth_min_10[i] <- breadth_10
  } else if (ref_length == 0) {
    tmp_df$cov_breadth_min_1[i] <- 0
    tmp_df$cov_breadth_min_5[i] <- 0
    tmp_df$cov_breadth_min_10[i] <- 0
  }
}

# Create column for subtype and Sample_ref
tmp_df <- as_tibble(tmp_df)

df_coverage <- tmp_df %>%
  # Don't need first mapping data
  filter(reference != "first_mapping") %>%
  # Create columns for major and minor coverage
  mutate(Major_cov_breadth_min_1 = case_when(first_major_minor == "major" ~ cov_breadth_min_1)) %>%
  mutate(Minor_cov_breadth_min_1 = case_when(first_major_minor == "minor" ~ cov_breadth_min_1)) %>%
  mutate(Major_cov_breadth_min_5 = case_when(first_major_minor == "major" ~ cov_breadth_min_5)) %>%
  mutate(Minor_cov_breadth_min_5 = case_when(first_major_minor == "minor" ~ cov_breadth_min_5)) %>%
  mutate(Major_cov_breadth_min_10 = case_when(first_major_minor == "major" ~ cov_breadth_min_10)) %>%
  mutate(Minor_cov_breadth_min_10 = case_when(first_major_minor == "minor" ~ cov_breadth_min_10)) %>%
  # Create columns for major and minor average depth
  mutate(Major_avg_depth = case_when(first_major_minor == "major" ~ avg_depth)) %>%
  mutate(Minor_avg_depth = case_when(first_major_minor == "minor" ~ avg_depth)) %>%
  # Create columns for the major and minor references
  mutate(Major_reference = case_when(first_major_minor == "major" ~ reference)) %>%
  mutate(Minor_reference = case_when(first_major_minor == "minor" ~ reference)) %>%
  mutate(Major_reference = str_remove(Major_reference, "_major"),
         Minor_reference = str_remove(Minor_reference, "_minor")) %>%
  # Create one row per sample
  select(-reference, -first_major_minor, -cov_breadth_min_1, -cov_breadth_min_5, -cov_breadth_min_10, -avg_depth) %>%
  group_by(sampleName) %>%
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)


# Length of contigs -------------------------------------------------------

# Print the length of the longest scaffold matching the given reference
blast_files <- list.files(path = path_7, pattern = "txt$", full.names = TRUE)

# Start with an empty tibble having correct column types
df_contigs <- tibble(
  scaffold_length = numeric(),
  reference       = character(),
  sampleName      = character()
)

if (length(blast_files) == 0) {
  message("No BLAST files found in: ", path_7)
} else {
  for (bf in blast_files) {
    sampleName <- str_split(basename(bf), "\\.")[[1]][1]

    # If the file is missing or zero-length, record a row with NA values
    if (!file.exists(bf) || file.size(bf) == 0) {
      df_contigs <- bind_rows(
        df_contigs,
        tibble(scaffold_length = NA_real_, reference = NA_character_, sampleName = sampleName)
      )
      next
    }

    # Try to read the file; read everything as character to avoid parsing errors
    dat <- tryCatch(
      read_tsv(bf, col_names = FALSE, col_types = cols(.default = col_character()), progress = FALSE),
      error = function(e) {
        warning("Failed to read '", bf, "': ", conditionMessage(e))
        NULL
      }
    )

    # If read failed or produced no rows, add NA row and continue
    if (is.null(dat) || nrow(dat) == 0) {
      df_contigs <- bind_rows(
        df_contigs,
        tibble(scaffold_length = NA_real_, reference = NA_character_, sampleName = sampleName)
      )
      next
    }

    # Ensure expected columns exist (X1 = query header, X2 = subject header)
    if (!"X1" %in% names(dat)) dat$X1 <- NA_character_
    if (!"X2" %in% names(dat)) dat$X2 <- NA_character_

    # Extract genotype (text before first underscore) and scaffold_length using regex
    processed <- dat %>%
      mutate(
        genotype = str_extract(X2, "^[^_]+"),
        scaffold_length = as.numeric(str_extract(X1, "(?<=_length_)[0-9]+"))
      ) %>%
      # prefer rows with largest scaffold_length per genotype (NA scaffold_length sorts last)
      arrange(desc(scaffold_length)) %>%
      group_by(genotype) %>%
      # Select the row with the longest scaffold lengths for each genotype/blast hit
      slice_head(n = 1) %>%
      ungroup() %>%
      # remove duplicate queries (same X1) keeping the first
      # Sometimes the same contigs has two or more hits
      distinct(X1, .keep_all = TRUE) %>%
      select(scaffold_length, reference = X2) %>%
      mutate(sampleName = sampleName)

    # If processing resulted in zero rows, add NA row; otherwise append results
    if (nrow(processed) == 0) {
      df_contigs <- bind_rows(
        df_contigs,
        tibble(scaffold_length = NA_real_, reference = NA_character_, sampleName = sampleName)
      )
    } else {
      df_contigs <- bind_rows(df_contigs, processed)
    }
  }
}

# GLUE --------------------------------------------------------------------

glue_file <- list.files(path = path_8, pattern = "GLUE_collected_report_major.tsv$", full.names = TRUE)
glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character()))

# Collect also the minor GLUE report
glue_file_minor <- list.files(path = path_8, pattern = "GLUE_collected_report_minor.tsv$", full.names = TRUE)
glue_report_minor <- read_tsv(glue_file_minor, col_types = cols(GLUE_subtype = col_character()))

# Extract the GLUE genotypes and subtypes for major and minor and compare them

if (nrow(glue_report) > 0) {
  major_gt <- glue_report %>%
    select(Sample, GLUE_genotype, GLUE_subtype) %>%
    rename(Major_genotype = GLUE_genotype,
           Major_subtype = GLUE_subtype)
}

if (nrow(glue_report_minor) > 0) {
  minor_gt <- glue_report_minor %>%
    select(Sample, GLUE_genotype, GLUE_subtype) %>%
    rename(Minor_genotype = GLUE_genotype,
           Minor_subtype = GLUE_subtype)
}

# Check if the same genotype has been called for major and minor. If Yes, then minor is not typable
if (exists("major_gt") & exists("minor_gt")) {
  gt_check <- major_gt %>%
    left_join(minor_gt, by = "Sample") %>%
    mutate(
      identical_geno = case_when(
        Major_genotype == Minor_genotype ~ "YES",
        is.na(Minor_genotype) ~ NA,
        .default = "NO"
      ),
      identical_subgeno = case_when(
        Major_subtype == Minor_subtype ~ "YES",
        is.na(Minor_subtype) ~ NA,
        .default = "NO"
      )
    )
}

# Sequencer ID ------------------------------------------------------------
id_files <- list.files(path = path_9, pattern = "sequencerID.tsv$", full.names = TRUE)

if (length(id_files > 0)) {
    # Empty df
    id_df <- as.data.frame(matrix(nrow = length(id_files), ncol = 2))
    colnames(id_df) <- c("sampleName", "sequencer_id")
}

for (i in 1:length(id_files)) {
  try(rm(id))

  # Get sample name
  id_df$sampleName[i] <- str_split(basename(id_files[i]), "\\.")[[1]][1]

  # Get the sequencer id
  id <- read_tsv(id_files[i], col_names = FALSE)

  # Check if file is empty (i.e., there were no FASTQ reads for this sample)
  if (nrow(id) == 0 || ncol(id) == 0) {
    # Set sequencer ID to NA
    id_df$sequencer_id[i] <- NA
  } else if (str_detect(id$X1, "^@SRR")) {
    # Extract the first field if header start with '@SRR'
    id_df$sequencer_id[i] <- id %>% pull(X1)
  } else {
    id_df$sequencer_id[i] <- id %>%
      pull(X1) %>%
      # Extract string up to the first ":".
      # The "?" means a "lazy", or non-greedy, match to get the shortest string that satisfies the criteria.
      # This is useful because there are several ":"
      str_extract("^.*?:") %>%
      # Remove the leading "@" and the trailing ":"
      str_remove_all("^@|:$")
  }
}

id_df <- as_tibble(id_df)

# Variation ---------------------------------------------------------------
# Read both the major and minor variation plots
variation_plot_files <- list.files(path = path_10, pattern = ".*variation_plot.*\\.png$", full.names = TRUE)

# Create R-code that will gather all the variation plots, then plot them as a grid with four columns and as many rows as needed.
# Make separate grids for files containing the string "major" and "minor"
if (length(variation_plot_files) > 0) {
  # Create a grid of plots for major and minor
  major_plots <- variation_plot_files[grepl("major", variation_plot_files)]
  minor_plots <- variation_plot_files[grepl("minor", variation_plot_files)]

  # Create a grid of plots for major
  if (length(major_plots) > 0) {
    # Read images and convert to grobs
    image_list <- lapply(major_plots, function(file) {
      img <- png::readPNG(file)
      grid::rasterGrob(img, interpolate = TRUE)
    })

    # Create a list of ggplot objects containing the images
    plot_list <- lapply(image_list, function(g) {
      ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void()
    })

    # Arrange the plots in a grid with 4 columns
    grid_arrange_major <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

    # Save the arranged grid to a PNG file with a white background
    ggsave(filename = "Variation_plot_major.png", plot = grid_arrange_major, width = 12, height = 8, dpi = 300, bg = "white")
  }

  # Create a grid of plots for minor
  if (length(minor_plots) > 0) {    # Read images and convert to grobs
    image_list <- lapply(minor_plots, function(file) {
      img <- png::readPNG(file)
      grid::rasterGrob(img, interpolate = TRUE)
    })

    # Create a list of ggplot objects containing the images
    plot_list <- lapply(image_list, function(g) {
      ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void()
    })

    # Arrange the plots in a grid with 4 columns
    grid_arrange_minor <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

    # Save the arranged grid to a PNG file with a white background
    ggsave(filename = "Variation_plot_minor.png", plot = grid_arrange_minor, width = 12, height = 8, dpi = 300, bg = "white")
  }
}

# Join dataframes ---------------------------------------------------------

# Start with the original input samplesheeet and extract only the sample names. This is to ensure that all samples are included in the final summary, even if they have no data.
input_samplesheet <- read_csv(samplesheet) %>%
  select("sampleName" = sample)

final <- input_samplesheet %>%
  # Add sequencer id
  left_join(id_df, join_by(sampleName)) %>%
  # Add number of raw and trimmed reads
  left_join(trimmed_df, by = "sampleName") %>%
  # Add number of classified reads from Kraken2
  left_join(kraken_df, by = "sampleName") %>%
  # Add total mapped reads from first mapping
  left_join(parsefirstmapping_df, join_by(sampleName)) %>%
  # Add mapped reads stats
  left_join(df_mapped_reads, join_by(sampleName)) %>%
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference))

if (nrow(glue_report) > 0) {
  final <- final %>%
    # Add glue result. Only Majority currently
    left_join(glue_report, by = c("sampleName" = "Sample"))
}

# Add script name and version
final <- final %>%
  add_column("pipeline_version" = script_name_version)


# Decide if a sample is "typable" or not
final <- final %>%
  mutate(major_typable = case_when(
    Major_cov_breadth_min_1 >= 10 & Major_avg_depth >= 2 ~ "YES",
    .default = "NO"
  )) %>%
  mutate(minor_typable = case_when(
    Minor_cov_breadth_min_1 >= 10 & Minor_avg_depth >= 2 ~ "YES",
    .default = "NO"
  ))

# If minor genotype is the same as major, then not typable. But only possible if there are minor glue reports available
# But allow for the co-infection of 1a and 1b even though these belong to the same genotype
if (nrow(glue_report) > 0 & exists("gt_check")) {
  final <- final %>%
    left_join(gt_check, by = c("sampleName" = "Sample")) %>%
    mutate(minor_typable = case_when(
      identical_geno == "NO" ~ "YES",                              # Different genotypes, so minor is typable
      identical_geno == "YES" & identical_subgeno == "NO" &
        ((Major_subtype == "1a" & Minor_subtype == "1b") |
         (Major_subtype == "1b" & Minor_subtype == "1a")) ~ "YES", # If the genotype is the same and subtypes are different, but must be 1a and 1b combination. Then allow typable Minor
      identical_geno == "YES" & identical_subgeno == "NO" ~ "NO",  # If the genotype is the same and subtypes are different, but not 1a and 1b combination. Then not typable Minor
      identical_geno == "YES" & identical_subgeno == "YES" ~ "NO", # Same genotype & same subtype → not typable
      is.na(identical_geno) ~ "UNKNOWN"
    ))
}

  # Add scaffold length info - for the moment not included
  # left_join(df_contigs, join_by(sampleName)) %>%
  # mutate(test = case_when(Majority_reference == reference ~ "OK",
  #                         Minority_reference == reference ~ "OK")) %>%
  # filter(test == "OK") %>%

# If the GLUE report is missing, and GLUE columns with NAs
if (!"GLUE_genotype" %in% colnames(final)) {
  final <- final %>%
    add_column("Reference" = NA_character_,
               "GLUE_genotype" = NA_character_,
               "GLUE_subtype" = NA_character_,
               "glecaprevir" = NA_character_,
               "glecaprevir_mut" = NA_character_,
               "glecaprevir_mut_short" = NA_character_,
               "grazoprevir" = NA_character_,
               "grazoprevir_mut" = NA_character_,
               "grazoprevir_mut_short" = NA_character_,
               "paritaprevir" = NA_character_,
               "paritaprevir_mut" = NA_character_,
               "paritaprevir_mut_short" = NA_character_,
               "voxilaprevir" = NA_character_,
               "voxilaprevir_mut" = NA_character_,
               "voxilaprevir_mut_short" = NA_character_,
               "NS34A" = NA_character_,
               "NS34A_short" = NA_character_,
               "daclasvir" = NA_character_,
               "daclasvir_mut" = NA_character_,
               "daclasvir_mut_short" = NA_character_,
               "elbasvir" = NA_character_,
               "elbasvir_mut" = NA_character_,
               "elbasvir_mut_short" = NA_character_,
               "ledipasvir" = NA_character_,
               "ledipasvir_mut" = NA_character_,
               "ledipasvir_mut_short" = NA_character_,
               "ombitasvir" = NA_character_,
               "ombitasvir_mut" = NA_character_,
               "ombitasvir_mut_short" = NA_character_,
               "pibrentasvir" = NA_character_,
               "pibrentasvir_mut" = NA_character_,
               "pibrentasvir_mut_short" = NA_character_,
               "velpatasvir" = NA_character_,
               "velpatasvir_mut" = NA_character_,
               "velpatasvir_mut_short" = NA_character_,
               "NS5A" = NA_character_,
               "NS5A_short" = NA_character_,
               "dasabuvir" = NA_character_,
               "dasabuvir_mut" = NA_character_,
               "dasabuvir_mut_short" = NA_character_,
               "sofosbuvir" = NA_character_,
               "sofosbuvir_mut" = NA_character_,
               "sofosbuvir_mut_short" = NA_character_,
               "NS5B" = NA_character_,
               "NS5B_short" = NA_character_,
               "HCV project version" = NA_character_,
               "GLUE engine version" = NA_character_,
               "PHE drug resistance extension version" = NA_character_,
               )
}

# Reorder columns
final <- final %>%
  select(sampleName,
         total_raw_reads,
         total_trimmed_reads,
         total_classified_reads,
         total_mapped_reads,
         fraction_mapped_reads_vs_median,
         Major_genotype_mapping,
         Major_reference,
         Minor_genotype_mapping,
         Minor_reference,
         major_typable,
         minor_typable,
         Reads_withdup_mapped_major,
         Reads_nodup_mapped_major,
         Percent_reads_mapped_of_trimmed_with_dups_major,
         Major_cov_breadth_min_5,
         Major_cov_breadth_min_10,
         percent_mapped_reads_major_firstmapping,
         Reads_withdup_mapped_minor,
         Reads_nodup_mapped_minor,
         Percent_reads_mapped_of_trimmed_with_dups_minor,
         Minor_cov_breadth_min_5,
         Minor_cov_breadth_min_10,
         percent_mapped_reads_minor_firstmapping,
         everything()) %>%
  distinct() %>% # Remove any duplicated rows from the different joins
  # If there are no minor genotype reports, then the identical_geno and identical_subgeno columns will not exist. Therefore use any_of in case they are not there
  select(-any_of(c("Major_minor", "identical_geno", "identical_subgeno")))

# Write file
write_csv(final, file = "Summary.csv")

# Write file for MultiQC
# Add MultiQC info lines
header <- c("# id: 'summary'",
            "# section_name: 'Summary'",
            "# description: 'These statistics are generated from the process SUMMARIZE and the R script summarize.R",
            "# format: 'csv'")

# Convert final data to data frame
tt <- as.data.frame(final)

# Set up file name for writing to (NB, can't use capital S i summary for MultiQC to pick it up)
file <- "summary_mqc.csv"

# Add MultiQC header to file
#write_lines(header, file)

# Add the column names to file
tt %>% colnames() %>% paste0(collapse = ",") %>% write_lines(file, append = TRUE)

# Write the data to file
write_csv(tt, file, append = TRUE) # colnames will not be included

