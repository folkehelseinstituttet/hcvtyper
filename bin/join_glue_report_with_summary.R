#!/usr/bin/env Rscript

library(tidyverse)

glue_file <- list.files(pattern = "^GLUE_collected_report.tsv$", full.names = TRUE, recursive=TRUE)

glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character()))

if (nrow(glue_report) > 0) {
  glue_report <- glue_report %>%
    # Only keep the majority reports for the summary
    filter(Major_minor == "major")
}

# Join dataframes ---------------------------------------------------------

final <-
  # Combine mapped reads data
  full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference)) %>%
  # Add coverage
  left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference))

if (nrow(glue_report) > 0) {
  final <- final %>%
    # Add glue result. Only Majority currently
    left_join(glue_report, by = c("sampleName" = "Sample"))
}

  # Add scaffold length info - for the moment not included
  # left_join(df_scaffolds, join_by(sampleName)) %>%
  # mutate(test = case_when(Majority_reference == reference ~ "OK",
  #                         Minority_reference == reference ~ "OK")) %>%
  # filter(test == "OK") %>%

# Reorder columns
final <- final %>%
  select(sampleName,
         total_raw_read_pairs,
         total_trimmed_read_pairs,
         total_classified_read_pairs,
         Major_genotype_mapping,
         Major_reference,
         Minor_genotype_mapping,
         Minor_reference,
         Read_pairs_withdup_mapped_major,
         Read_pairs_nodup_mapped_major,
         Percent_read_pairs_mapped_of_trimmed_with_dups_major,
         Major_cov_breadth_min_5,
         Read_pairs_withdup_mapped_minor,
         Read_pairs_nodup_mapped_minor,
         Percent_read_pairs_mapped_of_trimmed_with_dups_minor,
         Minor_cov_breadth_min_5,
         everything()) #%>%
  #select(-Reference, -Major_minor)

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
