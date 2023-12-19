#!/usr/bin/env Rscript

library(tidyverse)

# Number of mapped reads --------------------------------------------------
path_1 <- "kraken_classified/"
path_2 <- "stats_withdup/"
path_3 <- "stats_markdup/"
path_4 <- "depth/"
path_5 <- "blast/"
path_6 <- "glue/"


# Total trimmed reads with duplicates -------------------------------------
# List files
kraken_files <- list.files(path = path_1, pattern = "kraken2.report.txt$", full.names = TRUE)

# Empty df
kraken_df <- as.data.frame(matrix(nrow = length(kraken_files), ncol = 2))
colnames(kraken_df) <- c("sampleName", "total_trimmed_reads_with_dups")

for (i in 1:length(kraken_files)) {
  try(rm(kraken_stats))
  # Get sample name
  kraken_df$sampleName[i] <- str_split(basename(kraken_files[i]), "\\.")[[1]][1]
  
  # Get total of trimmed sequences put in to the mapping. Adding try() if no root sequences
  try(kraken_stats <- read_tsv(kraken_files[i], col_names = FALSE) %>% filter(X6 == "root") %>% pull(X2))
  if (exists("kraken_stats") & length(kraken_stats) > 0) {
    kraken_df$total_trimmed_reads_with_dups[i] <- kraken_stats 
  }
}
kraken_df <- as_tibble(kraken_df)

# Reads mapped with duplicates --------------------------------------------
# List files
stats_files <- list.files(path = path_2, pattern = "\\withdup.stats$", full.names = TRUE)

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
  
  # Get number of mapped reads after duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)
  tmp_df$trimmed_reads_withdups_mapped[i] <- mapped_reads
}
tmp_df <- as_tibble(tmp_df)

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
  mutate(total_trimmed_reads_with_dups = as.integer(total_trimmed_reads_with_dups),
         Reads_withdup_mapped_major = as.integer(Reads_withdup_mapped_major),
         Reads_withdup_mapped_minor = as.integer(Reads_withdup_mapped_minor)) %>% 
         #Reads_withdup_mapped_first_mapping = as.integer(Reads_withdup_mapped_first_mapping)) %>% 
  mutate(Percent_reads_mapped_with_dups_major = Reads_withdup_mapped_major / total_trimmed_reads_with_dups * 100,
         Percent_reads_mapped_with_dups_minor = Reads_withdup_mapped_minor / total_trimmed_reads_with_dups * 100) %>% 
         #Percent_reads_mapped_with_dups_first_mapping = Reads_withdup_mapped_first_mapping / total_trimmed_reads_with_dups * 100) %>% 
  # Create one row per sample
  select(-first_major_minor) %>% 
  group_by(sampleName) %>% 
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Reads mapped no duplicates ----------------------------------------------
# List files
stats_files <- list.files(path = path_3, pattern = "nodup.stats$", full.names = TRUE)

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

# List files
cov_files <- list.files(path = path_4, pattern = "tsv$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(cov_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_5", "first_major_minor")

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
  
  # Nr. of positions with coverage >= 5
  pos <- nrow(
    cov %>% 
      filter(X3 >= 5)
  )
  
  # Coverage breadth
  breadth <- round(pos / ref_length * 100, digits = 2)
  tmp_df$cov_breadth_min_5[i] <- breadth
  
}

# Create column for subtype and Sample_ref
tmp_df <- as_tibble(tmp_df)

df_coverage <- tmp_df %>%
  # Don't need first mapping data
  filter(reference != "first_mapping") %>% 
  # Create columns for major and minor coverage
  mutate(Major_cov_breadth_min_5 = case_when(first_major_minor == "major" ~ cov_breadth_min_5)) %>% 
  mutate(Minor_cov_breadth_min_5 = case_when(first_major_minor == "minor" ~ cov_breadth_min_5)) %>% 
  # Create columns for the major and minor references
  mutate(Major_reference = case_when(first_major_minor == "major" ~ reference)) %>% 
  mutate(Minor_reference = case_when(first_major_minor == "minor" ~ reference)) %>% 
  mutate(Major_reference = str_remove(Major_reference, "_major"),
         Minor_reference = str_remove(Minor_reference, "_minor")) %>% 
  # Create one row per sample
  select(-reference, -first_major_minor, -cov_breadth_min_5) %>% 
  group_by(sampleName) %>% 
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1) 


# Length of scaffolds -----------------------------------------------------

# Print the length of the longest scaffold matching the given reference
blast_files <- list.files(path = path_5, pattern = "txt$", full.names = TRUE)

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

# GLUE --------------------------------------------------------------------

glue_file <- list.files(path = path_6, pattern = "GLUE_collected_report.tsv$", full.names = TRUE)
glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character())) %>% 
  # Only keep the majority reports for the summary
  filter(Major_minor == "major")

# glue_reports <- list.files(path = path_5, pattern = "GLUE_report.tsv$", full.names = TRUE) %>% 
#   # Keep the file names as the names of the list elements
#   set_names() %>% 
#   map(read_tsv, col_types = cols(GLUE_subtype = col_character())) %>% 
#   # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
#   # The column name will be "sampleName"
#   bind_rows(.id = "sampleName") %>% 
#   # Clean up sampleName
#   mutate(sampleName = str_remove(sampleName, "json//")) %>% # "json//
#   # Create a new column that keeps the sample name, reference for mapping and major/minor
#   mutate("Sample_ref" = str_remove(sampleName, "_GLUE_report\\.tsv")) %>% 
#   # Keep the reference in a separate column
#   separate(sampleName, into = c("sampleName", "reference", "major_minor"), sep = "\\.") %>% 
#   mutate(major_minor = str_remove(major_minor, "_GLUE_report")) 

# Join dataframes ---------------------------------------------------------

final <- 
  # Combine mapped reads data
  full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference)) %>% 
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference)) %>% 
  # Add scaffold length info - for the moment not included
  # left_join(df_scaffolds, join_by(sampleName)) %>% 
  # mutate(test = case_when(Majority_reference == reference ~ "OK",
  #                         Minority_reference == reference ~ "OK")) %>% 
  # filter(test == "OK") %>%
  # Add glue result. Only Majority currently
  left_join(glue_report, by = c("sampleName" = "Sample"))




# # Join on sampleName, reference, major_minor
# tmp_df <- left_join(tmp_df, glue_reports)
# 
# # Join all objects on 
# # Join with stats object if sampleName and references are the same.
# df <- tmp_df %>%
#   dplyr::left_join(stats, by = join_by(sampleName, Majority_reference, Minority_reference))

# Join the blast/scaffold lengths only on sampleName and reference

# Write file
write_csv(final, file = "Genotype_mapping_summary_long.csv")
