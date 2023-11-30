#!/usr/bin/env Rscript

library(tidyverse)

# Number of mapped reads --------------------------------------------------
path_1 <- "stats_withdup/"
path_2 <- "stats_markdup/"
path_3 <- "depth/"
path_4 <- "blast/"

# Reads mapped with duplicates --------------------------------------------
# List files
stats_files <- list.files(path = path_1, pattern = "\\withdup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 5))
colnames(tmp_df) <- c("sampleName", "reference", "first_major_minor", "total_trimmed_reads_with_dups", "trimmed_reads_withdups_mapped")

for (i in 1:length(stats_files)) {
  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][1]
  
  # Get reference name
  tmp_df$reference[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][2]
  
  # Get major or minor
  tmp_df$first_major_minor[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][3]
  
  # Read the mapping stats
  map_stats <- read_tsv(stats_files[i], col_names = FALSE, comment = "#") 
  
  # Get total of trimmed sequences put in to the mapping
  raw_reads <- map_stats %>% filter(X2 == "raw total sequences:") %>% pull(X3)
  tmp_df$total_trimmed_reads_with_dups[i] <- raw_reads
  
  # Get number of mapped reads after duplicate removal
  mapped_reads <- map_stats %>% filter(X2 == "reads mapped:") %>% pull(X3)
  tmp_df$trimmed_reads_withdups_mapped[i] <- mapped_reads
}
tmp_df <- as_tibble(tmp_df)

df_with_dups <- tmp_df %>% 
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_withdup_mapped_majority = case_when(first_major_minor == "majority" ~ trimmed_reads_withdups_mapped)) %>%
  mutate(Reads_withdup_mapped_minority = case_when(first_major_minor == "minority" ~ trimmed_reads_withdups_mapped)) %>%
  mutate(Reads_withdup_mapped_first_mapping = case_when(first_major_minor == "first_mapping" ~ trimmed_reads_withdups_mapped)) %>% 
  select(-trimmed_reads_withdups_mapped) %>% 
  # Create columns for the major and minor references
  mutate(Majority_reference = case_when(first_major_minor == "majority" ~ reference)) %>% 
  mutate(Minority_reference = case_when(first_major_minor == "minority" ~ reference)) %>% 
  mutate(Majority_reference = str_remove(Majority_reference, "_major"),
         Minority_reference = str_remove(Minority_reference, "_minor")) %>% 
  select(-reference) %>% 
  # Calculate percent of the trimmed reads mapped
  mutate(total_trimmed_reads_with_dups = as.integer(total_trimmed_reads_with_dups),
         Reads_withdup_mapped_majority = as.integer(Reads_withdup_mapped_majority),
         Reads_withdup_mapped_minority = as.integer(Reads_withdup_mapped_minority),
         Reads_withdup_mapped_first_mapping = as.integer(Reads_withdup_mapped_first_mapping)) %>% 
  mutate(Percent_reads_mapped_with_dups_majority = Reads_withdup_mapped_majority / total_trimmed_reads_with_dups * 100,
         Percent_reads_mapped_with_dups_minority = Reads_withdup_mapped_minority / total_trimmed_reads_with_dups * 100,
         Percent_reads_mapped_with_dups_first_mapping = Reads_withdup_mapped_first_mapping / total_trimmed_reads_with_dups * 100) %>% 
  # Create one row per sample
  select(-first_major_minor) %>% 
  group_by(sampleName) %>% 
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# Reads mapped no duplicates ----------------------------------------------
# List files
stats_files <- list.files(path = path_2, pattern = "\\markdup.stats$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(stats_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "first_major_minor", "trimmed_reads_nodups_mapped")

for (i in 1:length(stats_files)) {
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
  mutate(Majority_genotype_mapping = case_when(first_major_minor == "majority" ~ genotype)) %>%
  mutate(Minority_genotype_mapping = case_when(first_major_minor == "minority" ~ genotype)) %>% 
  select(-genotype) %>% 
  # Create columns for reads mapped to major and minor genotype
  mutate(Reads_nodup_mapped_majority = case_when(first_major_minor == "majority" ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_minority = case_when(first_major_minor == "minority" ~ trimmed_reads_nodups_mapped)) %>%
  mutate(Reads_nodup_mapped_first_mapping = case_when(first_major_minor == "first_mapping" ~ trimmed_reads_nodups_mapped)) %>% 
  select(-trimmed_reads_nodups_mapped) %>% 
  #mutate(Percent_mapped_majority = case_when(first_major_minor == "majority" ~ Percent_trimmed_reads_mapped)) %>% 
  #mutate(Percent_mapped_minority = case_when(first_major_minor == "minority" ~ Percent_trimmed_reads_mapped)) %>%
  # Create columns for the major and minor references
  mutate(Majority_reference = case_when(first_major_minor == "majority" ~ reference)) %>% 
  mutate(Minority_reference = case_when(first_major_minor == "minority" ~ reference)) %>% 
  mutate(Majority_reference = str_remove(Majority_reference, "_major"),
         Minority_reference = str_remove(Minority_reference, "_minor")) %>% 
  # Create one row per sample
  select(-reference, -first_major_minor) %>% 
  group_by(sampleName) %>% 
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1)

# stats <- list.files(path = path_1, pattern = "\\.stats$", full.names = TRUE) %>% 
#   # Keep the file names as the names of the list elements
#   set_names() %>% 
#   map(read_tsv, col_names = FALSE, comment = "#") %>% 
#   # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
#   # The column name will be "sampleName"
#   bind_rows(.id = "sampleName") %>% 
#   # Clean up sampleName
#   mutate(sampleName = basename(sampleName)) %>%
#   separate(sampleName, into = c("sampleName", "reference", "first_major_minor"), sep = "\\.") %>% 
#   # Extract relevant mapping info and get it on one row
#   filter(X2 == "raw total sequences:" | X2 == "reads mapped:") %>% # Keep total reads in and reads mapped
#   mutate(X2 = str_remove(X2, ":")) %>% 
#   pivot_wider(names_from = X2, values_from = X3) %>% 
#   # Create a new column that keeps the sample name, reference for mapping and major/minor together
#   unite("Sample_ref", c("sampleName", "reference", "first_major_minor"), sep = "_", remove = FALSE) %>% 
#   # Select relevant columns and rename
#   select(sampleName,
#          reference,
#          first_major_minor,
#          Sample_ref,
#          "Reads_mapped" = `reads mapped`,
#          "Total_trimmed_reads" = `raw total sequences`) %>% 
#   mutate(Reads_mapped = as.numeric(Reads_mapped),
#          Total_trimmed_reads = as.numeric(Total_trimmed_reads)) %>% 
#   # Calculate percent of trimmed reads mapped
#   mutate("Percent_trimmed_reads_mapped" = Reads_mapped / Total_trimmed_reads * 100) %>% 
#   # Create columns for major and minor
#   separate(reference, into = c("genotype", NA), sep = "_", remove = F) %>% 
#   mutate(Majority_genotype_mapping = case_when(first_major_minor == "majority" ~ genotype)) %>%
#   mutate(Minority_genotype_mapping = case_when(first_major_minor == "minority" ~ genotype)) %>% 
#   select(-genotype) %>% 
#   # Create columns for reads mapped to major and minor genotype
#   mutate(Reads_mapped_majority = case_when(first_major_minor == "majority" ~ Reads_mapped)) %>%
#   mutate(Reads_mapped_minority = case_when(first_major_minor == "minority" ~ Reads_mapped)) %>%
#   mutate(Percent_mapped_majority = case_when(first_major_minor == "majority" ~ Percent_trimmed_reads_mapped)) %>% 
#   mutate(Percent_mapped_minority = case_when(first_major_minor == "minority" ~ Percent_trimmed_reads_mapped)) %>%
#   # Create columns for the major and minor references
#   mutate(Majority_reference = case_when(first_major_minor == "majority" ~ reference)) %>% 
#   mutate(Minority_reference = case_when(first_major_minor == "minority" ~ reference)) %>% 
#   mutate(Majority_reference = str_remove(Majority_reference, "_major"),
#          Minority_reference = str_remove(Minority_reference, "_minor")) %>% 
#   # Create one row per sample
#   select(-Sample_ref, -reference, -first_major_minor, -Reads_mapped, -Percent_trimmed_reads_mapped) %>% 
#   group_by(sampleName) %>% 
#   # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
#   fill(everything(), .direction = "downup") %>%
#   slice(1)

  
# Coverage ----------------------------------------------------------------

# List files
cov_files <- list.files(path = path_3, pattern = "txt\\.gz$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(cov_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_5", "first_major_minor")

for (i in 1:length(cov_files)) {
  
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
  mutate(Majority_cov_breadth_min_5 = case_when(first_major_minor == "majority" ~ cov_breadth_min_5)) %>% 
  mutate(Minority_cov_breadth_min_5 = case_when(first_major_minor == "minority" ~ cov_breadth_min_5)) %>% 
  # Create columns for the major and minor references
  mutate(Majority_reference = case_when(first_major_minor == "majority" ~ reference)) %>% 
  mutate(Minority_reference = case_when(first_major_minor == "minority" ~ reference)) %>% 
  mutate(Majority_reference = str_remove(Majority_reference, "_major"),
         Minority_reference = str_remove(Minority_reference, "_minor")) %>% 
  # Create one row per sample
  select(-reference, -first_major_minor, -cov_breadth_min_5) %>% 
  group_by(sampleName) %>% 
  # Fill missing values per group (i.e. sampleName. Direction "downup" fill values from both rows)
  fill(everything(), .direction = "downup") %>%
  slice(1) 


# Length of scaffolds -----------------------------------------------------

# Print the length of the longest scaffold matching the given reference
blast_files <- list.files(path = path_4, pattern = "txt$", full.names = TRUE)

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

glue_file <- list.files(pattern = "GLUE_collected_report.tsv$", full.names = TRUE)
glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character()))

# glue_reports <- list.files(path = path_4, pattern = "GLUE_report.tsv$", full.names = TRUE) %>% 
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
  full_join(df_with_dups, df_nodups, join_by(sampleName, Majority_reference, Minority_reference)) %>% 
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Majority_reference, Minority_reference))  
  # Add scaffold length info - for the moment not included
  #left_join(df_scaffolds, join_by(sampleName)) %>% 
  #mutate(test = case_when(Majority_reference == reference ~ "OK",
  #                        Minority_reference == reference ~ "OK")) %>% View()
  #filter(test == "OK") %>%
  # Add glue resulst - not currently. Will add for Majority and Minority separately
  #left_join(glue_report, by = c("sampleName" = "Sample")) %>% View()




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
