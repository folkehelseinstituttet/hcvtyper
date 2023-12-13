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
fastq <- fastq[str_detect(fastq, agens)]

R1 <- sort(fastq[grep("R1", fastq)])
R2 <- sort(fastq[grep("R2", fastq)])

df <- as_tibble(cbind(R1, R2))

# Check that the R1 and R2 files are correctly paired
tmp <- df %>%
  mutate(tmpR1 = gsub("_.*", "", basename(R1)),
         tmpR2 = gsub("_.*", "", basename(R2))) %>%
  select(tmpR1, tmpR2)

if (identical(tmp$tmpR1, tmp$tmpR2)) {
  df <- df %>%
    mutate(sample_id = gsub("_.*", "", basename(R1))) %>%
    select("sample" = sample_id,
           "fastq_1" = R1,
           "fastq_2" = R2)
} else {
  print("R1 and R2 files not correctly paired")
}

write_csv(df, outfile)
