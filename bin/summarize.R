#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# Define variables --------------------------------------------------------

# NB! Script name is currently hard coded. Needs to be changed
script_name_version <- "viralseq v1.0.1"
stringency_1 <- args[1]
stringency_2 <- args[2]

path_1 <- "cutadapt/"
path_2 <- "kraken_classified/"
path_3 <- "stats_withdup/"
path_4 <- "stats_markdup/"
path_5 <- "depth/"
path_6 <- "blast/"
path_7 <- "glue/"
path_8 <- "id/"



# Cutadapt ----------------------------------------------------------------

# List files
cutadapt_files <- list.files(path = path_1, pattern = ".log$", full.names = TRUE)

# Empty df
cutadapt_df <- as.data.frame(matrix(nrow = length(cutadapt_files), ncol = 3))
colnames(cutadapt_df) <- c("sampleName", "total_raw_reads", "total_trimmed_reads")

for (i in 1:length(cutadapt_files)) {
  try(rm(cutadapt_stats))

  # Get sample name
  cutadapt_df$sampleName[i] <- str_split(basename(cutadapt_files[i]), "\\.")[[1]][1]

  # Get number of raw and trimmed read pairs
  try(cutadapt_stats <- as_tibble(read_lines(cutadapt_files[i])))

  raw <- filter(cutadapt_stats, str_detect(value, "Total read pairs processed:"))

  if (nrow(raw) > 0) {
    tmp <- str_split(raw, "\\s+")[[1]] # Split on white space and get the list content
    tmp <- as.numeric(str_remove_all(tmp[length(tmp)], ",")) # extract the last element which contains the pair number, remove commas and create a numeric
    tmp <- tmp*2 # Double to get reads
    cutadapt_df$total_raw_reads[i] <- tmp
  }

  trimmed <- filter(cutadapt_stats, str_detect(value, "Pairs written \\(passing filters\\)"))

  if (nrow(trimmed) > 0) {
    tmp <- str_split(trimmed, "\\s+")[[1]] # Split on white space and get the list content
    tmp <- as.numeric(str_remove_all(tmp[length(tmp)-1], ",")) # extract the second last element which contains the pair number, remove commas and create a numeric
    tmp <- tmp*2 # Double to get reads
    cutadapt_df$total_trimmed_reads[i] <- tmp
  }
}
cutadapt_df <- as_tibble(cutadapt_df)

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

# Reads mapped with duplicates --------------------------------------------
# List files
stats_files <- list.files(path = path_3, pattern = "\\withdup.stats$", full.names = TRUE)

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

# Add number of raw and trimmed reads
tmp_df <- left_join(tmp_df, cutadapt_df, by = "sampleName")

# Add number of classified reads from Kraken2
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
stats_files <- list.files(path = path_4, pattern = "nodup.stats$", full.names = TRUE)

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


# Coverage ----------------------------------------------------------------

# Add both breadth (in percent) and depth (average depth)
# All this is without duplicates

# List files
cov_files <- list.files(path = path_5, pattern = "tsv$", full.names = TRUE)

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


# Length of scaffolds -----------------------------------------------------

# Print the length of the longest scaffold matching the given reference
blast_files <- list.files(path = path_6, pattern = "txt$", full.names = TRUE)

if (length(blast_files > 0)) {

# Empty df
df_scaffolds <- tribble(
  ~"scaffold_length", ~"reference", ~"sampleName",
  )

for (i in 1:length(blast_files)) {

  # Read the length of scaffolds
  blast_out <- read_tsv(blast_files[i], col_names = FALSE) %>%
    # Separate the genotype from the subject header
    separate(X2, into = c("genotype", NA), remove = FALSE) %>%
    # Get scaffold length info
    separate(X1, c(NA, NA, NA, "scaffold_length", NA, NA), sep = "_", remove = FALSE) %>%
    mutate(scaffold_length = as.numeric(scaffold_length)) %>%
    # Select the row with the longest scaffold lengths for each genotype/blast hit
    group_by(genotype) %>%
    slice_max(order_by = scaffold_length, n = 1) %>%
    ungroup() %>%
    # Sometimes the same scaffolds has two or more hits
    distinct(X1, .keep_all = TRUE) %>%
    select(scaffold_length,
           "reference" = X2) %>%
    # Legg til en kolonne med Sample Name
    add_column("sampleName" = str_split(basename(blast_files[i]), "\\.")[[1]][1])

  # Add to df_scaffolds
  df_scaffolds <- bind_rows(df_scaffolds, blast_out)

}
}
# GLUE --------------------------------------------------------------------

glue_file <- list.files(path = path_7, pattern = "GLUE_collected_report.tsv$", full.names = TRUE)
glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character()))

if (nrow(glue_report) > 0) {
  glue_report <- glue_report %>%
    # Only keep the majority reports for the summary
    filter(Major_minor == "major")
}

# Sequencer ID ------------------------------------------------------------
id_files <- list.files(path = path_8, pattern = "sequencerID.tsv$", full.names = TRUE)

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

  # Extract the first field if header does not start with '@SRR'
  if (str_detect(id$X1, "^@SRR")) {
    id_df$sequencer_id[i] <- id %>% pull(X1)
  } else {
    id_df$sequencer_id[i] <- id %>%
      # Extract string up to the first ":".
      # The "?" means a "lazy", or non-greedy, match to get the shortest string that satisfies the criteria.
      # This is useful because there are several ":"
      str_extract("^.*?:") %>%
      # Remove the leading "@" and the trailing ":"
      str_remove_all("^@|:$")
  }
}

id_df <- as_tibble(id_df)

# Join dataframes ---------------------------------------------------------

final <-
  # Combine mapped reads data
  full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Add sequencer id
  left_join(id_df, join_by(sampleName))

if (nrow(glue_report) > 0) {
  final <- final %>%
    # Add glue result. Only Majority currently
    left_join(glue_report, by = c("sampleName" = "Sample"))
}

# Add percentage of total classified reads belonging to major and minor genotype. Duplicates included. As a proxy for abundance.
final <- final %>%
  mutate(abundance_major = ( Reads_withdup_mapped_major / total_classified_reads) * 100 ) %>%
  mutate(abundance_minor = ( Reads_withdup_mapped_minor / total_classified_reads) * 100 )

# Add script name and version
final <- final %>%
  add_column("script_name_stringency" = script_name_version)


# Decide if a sample is "typbar" or not

final <- final %>%
  mutate(major_typbar = case_when(
    Major_cov_breadth_min_1 >= 10 & Major_avg_depth >= 2 ~ "YES",
    .default = "NO"
  )) %>%
  mutate(minor_typbar = case_when(
    Minor_cov_breadth_min_1 >= 10 & Minor_avg_depth >= 2 ~ "YES",
    .default = "NO"
  ))

  # Add scaffold length info - for the moment not included
  # left_join(df_scaffolds, join_by(sampleName)) %>%
  # mutate(test = case_when(Majority_reference == reference ~ "OK",
  #                         Minority_reference == reference ~ "OK")) %>%
  # filter(test == "OK") %>%

# Reorder columns
final <- final %>%
  select(sampleName,
         total_raw_reads,
         total_trimmed_reads,
         total_classified_reads,
         Major_genotype_mapping,
         Major_reference,
         Minor_genotype_mapping,
         Minor_reference,
         Reads_withdup_mapped_major,
         Reads_nodup_mapped_major,
         Percent_reads_mapped_of_trimmed_with_dups_major,
         Major_cov_breadth_min_5,
         Major_cov_breadth_min_10,
         abundance_major,
         Reads_withdup_mapped_minor,
         Reads_nodup_mapped_minor,
         Percent_reads_mapped_of_trimmed_with_dups_minor,
         Minor_cov_breadth_min_5,
         Minor_cov_breadth_min_10,
         abundance_minor,
         everything())

# Write file
write_csv(final, file = "Genotype_mapping_summary_long.csv")

# Write file for MultiQC
# Add MultiQC info lines
header <- c("# id: 'summary'",
            "# section_name: 'Summary'",
            "# description: 'These statistics are generated from the process SUMMARIZE and the R script summarize.R",
            "# format: 'csv'")

# Convert final data to data frame
tt <- as.data.frame(final)

# Set up file name for writing to
file <- "Genotype_mapping_summary_long_mqc.csv"

# Add MultiQC header to file
write_lines(header, file)

# Add the column names to file
tt %>% colnames() %>% paste0(collapse = ",") %>% write_lines(file, append = TRUE)

# Write the data to file
write_csv(tt, file, append = TRUE) # colnames will not be included


# Create LW import --------------------------------------------------------

lw_import <- final %>%
    # Change "." to ","
    #mutate(
    #    Percent_reads_mapped_of_trimmed_with_dups_major = str_replace(Percent_reads_mapped_of_trimmed_with_dups_major, "\\.", ","),
    #    Major_cov_breadth_min_1 = str_replace(Major_cov_breadth_min_1, "\\.", ","),
    #    Major_avg_depth = str_replace(Major_avg_depth, "\\.", ","),
    #    Major_cov_breadth_min_5 = str_replace(Major_cov_breadth_min_5, "\\.", ","),
    #    Major_cov_breadth_min_10 = str_replace(Major_cov_breadth_min_10, "\\.", ","),
    #    abundance_minor = str_replace(abundance_minor, "\\.", ","),
    #    Minor_cov_breadth_min_5 = str_replace(Minor_cov_breadth_min_5, "\\.", ","),
    #    Minor_avg_depth = str_replace(Minor_avg_depth, "\\.", ","),
    #    Minor_cov_breadth_min_5 = str_replace(Minor_cov_breadth_min_5, "\\.", ","),
    #    Minor_cov_breadth_min_10 = str_replace(Minor_cov_breadth_min_10, "\\.", ",")
  #) %>%
  select("Sample" = sampleName,
         "Percent mapped reads of trimmed:" = Percent_reads_mapped_of_trimmed_with_dups_major,
         "Majority genotype:" = Major_genotype_mapping,
         "Number of mapped reads:" = Reads_withdup_mapped_major,
         "Percent covered:" = Major_cov_breadth_min_1,
         "Number of mapped reads without duplicates:" = Reads_nodup_mapped_major,
         "Average depth without duplicates:" = Major_avg_depth,
         "Percent covered above depth=5 without duplicates:" = Major_cov_breadth_min_5,
         "Percent covered above depth=9 without duplicates:" = Major_cov_breadth_min_10,
         "Most abundant minority genotype:" = Minor_genotype_mapping,
         "Percent most abundant minority genotype:" = abundance_minor,
         "Number of mapped reads minor:" = Reads_withdup_mapped_minor,
         "Percent covered minor:" = Minor_cov_breadth_min_5,
         "Number of mapped reads minor without duplicates:" = Reads_nodup_mapped_minor,
         "Average depth minor without duplicates:" = Minor_avg_depth,
         "Percent covered above depth=5 minor without duplicates:" = Minor_cov_breadth_min_5,
         "Percent covered above depth=9 minor without duplicates:" = Minor_cov_breadth_min_10,
         "Script name and stringency:" = script_name_stringency,
         "Total number of reads before trim:" = total_raw_reads,
         "Total number of reads after trim:" = total_trimmed_reads,
         "Majority quality:" = major_typbar,
         "Minor quality:" = minor_typbar,
         everything()
         ) %>%
  select(-total_classified_reads,
         -Major_reference,
         -Minor_reference,
         -abundance_major,
         -Percent_reads_mapped_of_trimmed_with_dups_minor,
         -Reads_nodup_mapped_first_mapping,
         -Minor_cov_breadth_min_1
         ) %>%
  # Convert "YES" to "Typbar" and "NO" to "Ikke typbar" in the "Majority quality:" and "Minor quality:" columns
  mutate(`Majority quality:` = case_when(`Majority quality:` == "YES" ~ "Typbar",
                                         `Majority quality:` == "NO" ~ "Ikke typbar")) %>%
  mutate(`Minor quality:` = case_when(`Minor quality:` == "YES" ~ "Typbar",
                                      `Minor quality:` == "NO" ~ "Ikke typbar"))

# Remove column "Major_minor" if exists
# This column is not present if GLUE is dropped
if ("^Major_minor$" %in% colnames(lw_import)) {
  lw_import <- lw_import %>% select(-Major_minor)
}

# Write file
write_csv(lw_import, file = "Genotype_mapping_summary_long_LW_import.csv")


