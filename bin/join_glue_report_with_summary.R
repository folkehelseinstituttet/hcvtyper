#!/usr/bin/env Rscript

library(tidyverse)

summary_file <- list.files(path = "/summarize", pattern = "^Genotype_mapping_summary_long_LW_import.csv$", full.names = TRUE, recursive = TRUE)

summary <- read_csv(summary_file)

glue_file <- list.files(path = "/hcvglue/", pattern = "^GLUE_collected_report.tsv$", full.names = TRUE, recursive = TRUE)

glue_report <- read_tsv(glue_file, col_types = cols(GLUE_subtype = col_character()))

if (nrow(glue_report) > 0) {
  glue_report <- glue_report %>%
    # Only keep the majority reports for the summary
    filter(Major_minor == "major")
}

# Join dataframes ---------------------------------------------------------
if (nrow(glue_report) > 0) {
  final <- summary %>%
    # Add glue result. Only Majority currently
    left_join(glue_report, by = join_by(Sample)) %>%
    select(-Major_minor)
}

# Convert "." to "," as decimal separators
final <- final %>%
    mutate(across(where(is.numeric), ~ format(., decimal.mark = ",", scientific = FALSE)))

# Write file
write_tsv(final, file = "Genotype_mapping_summary_long_LW_import_with_glue.tsv")
