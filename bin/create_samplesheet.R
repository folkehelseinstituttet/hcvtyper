#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript create_samplesheet.R <path_to_fastq_folders> <samplesheet_name> <agent>", call.=FALSE)
}

folder  <- args[1] # "/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run820_Virus/Run820/"
outfile <- args[2] # "2023.01.19-HCV_Run829.csv"
agens   <- args[3] # "HCV"

# Get the fastq files
fastq <- list.files(folder,
           recursive = TRUE,
           full.names = TRUE,
           pattern = "gz$")

# Filter on agens
if (agens == "HCV") {
    fastq <- fastq[str_detect(fastq, agens)]
} else if (agens == "ROV" | agens == "RoV") {
    agens <- "R[Oo]V" # Set agens to allow both RoV and ROV
    fastq <- fastq[str_detect(fastq, agens)]
} else {
    print("Agens is not correct")
}


# Extract sample ID = first string before first underscore
extract_sample_id <- function(x) str_extract(basename(x), "^[^_]+")

# Split R1 and R2
R1 <- fastq[str_detect(fastq, "_R1")]
R2 <- fastq[str_detect(fastq, "_R2")]

# Create tibbles with sample ID
df_R1 <- tibble(
  sample = map_chr(R1, extract_sample_id),
  fastq_1 = R1
)

df_R2 <- tibble(
  sample = map_chr(R2, extract_sample_id),
  fastq_2 = R2
)

# Join R1 and R2 on sample ID. This ensures R1 and R2 are correctly paired
df <- left_join(df_R1, df_R2, by = "sample")

# Check for missing pairs
if (any(is.na(df$fastq_1)) || any(is.na(df$fastq_2))) {
  stop("Some samples are missing R1 or R2 files", call. = FALSE)
}

# Check that the R1 and R2 files are correctly paired
tmp <- df %>%
  mutate(tmpR1 = map_chr(fastq_1, ~ gsub("_.*", "", basename(.))),
         tmpR2 = map_chr(fastq_2, ~ gsub("_.*", "", basename(.)))) %>%
  select(tmpR1, tmpR2)

if (identical(tmp$tmpR1, tmp$tmpR2)) {
  print("All files are correctly paired")
} else {
  print("R1 and R2 files not correctly paired")
}

write_csv(df, outfile)
