#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# Define variables --------------------------------------------------------

# NB! Script name is currently hard coded. Needs to be changed
script_name <- "viralseq"

path_1 <- "id/"
path_2 <- "parse_phylogeny/"
path_3 <- "alignment_metrics/"
path_4 <- "depth/"
path_5 <- "stats_withdup/"
path_6 <- "stats_markdup/"
path_7 <- "cutadapt/"
path_8 <- "kraken_classified/"



# Sequencer ID ------------------------------------------------------------
id_files <- list.files(path = path_1, pattern = "sequencerID.tsv$", full.names = TRUE)

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


# Parse phylogeny, alignment and mapping statistics -----------------------

## Parse phylogeny

# Summarize metrics from parsing the phylogenetic tree including all contigs for a given gene and corresponding references.

# List files
parse_phylogeny_files <- list.files(path = path_2, pattern = ".csv$", full.names = TRUE)

# Empty df
parse_phylogeny_df <- tribble(
  ~"sampleName", ~"gene", ~"contig_name", ~"target_clade", ~"closest_sequence", ~"ratio",
)

for (i in 1:length(parse_phylogeny_files)) {
  try(rm(tmp))
  tmp <- read_csv(parse_phylogeny_files[i]) %>%
    add_column("sampleName" = str_split(basename(parse_phylogeny_files[i]), "\\.")[[1]][1]) %>%
    add_column("gene" = str_split(basename(parse_phylogeny_files[i]), "\\.")[[1]][3]) %>%
    select(sampleName,
           gene,
           "contig_name" = Sequence,
           "target_clade" = `Target clade prefix`,
           "closest_sequence" = `Closest sequence`,
           "ratio" = `Closest Count / Total Count`) %>%
    # Remove everything after "." in the contig name
    separate(contig_name, into = c("contig_name", NA), sep = "\\.")

  # Add to df
  parse_phylogeny_df <- bind_rows(parse_phylogeny_df, tmp)
}

# Join with sequencer ID
joined_df <- full_join(id_df, parse_phylogeny_df, by = join_by(sampleName))

## Alignment metrics

# Summarize metrics from the pairwise alignment between the de novo contig and the closest reference in the phylogeny

# List files
alignment_metrics_files <- list.files(path = path_3, pattern = ".csv$", full.names = TRUE)

# Empty df
alignment_metrics_df <- tribble(
  ~"sampleName", ~"gene", ~"contig_name", ~"closest_sequence", ~"percent_similarity", ~"aligned_length", ~"total_length",
)

for (i in 1:length(alignment_metrics_files)) {
  try(rm(tmp))
  tmp <- read_csv(alignment_metrics_files[i]) %>%
    add_column("sampleName" = str_split(basename(alignment_metrics_files[i]), "\\.")[[1]][1]) %>%
    add_column("gene" = str_split(basename(alignment_metrics_files[i]), "\\.")[[1]][3]) %>%
    # Clean up closest sequence name. Keep everything before the first ":"
    mutate(`Closest sequence` = str_split(`Closest sequence`, ":")[[1]][1]) %>%
    select(sampleName,
          gene,
          "contig_name" = Sequence,
          "closest_sequence" = `Closest sequence`,
          "percent_similarity" = `Percent Similarity`,
          "aligned_length" = `Aligned Length`,
          "total_length" = `Total Length`) %>%
    # Remove everything after "." in the contig name
    separate(contig_name, into = c("contig_name", NA), sep = "\\.")



  # Add to df
  alignment_metrics_df <- bind_rows(alignment_metrics_df, tmp)

}

# Join parse phylogeny and alignment metrics on the both sampleName, gene, contig name and closest sequence.
# This is to ensure that phylogeny results and alignment stats are generated on the same results.
# use full_join to keep all observations
joined_df <- full_join(parse_phylogeny_df, alignment_metrics_df, by = join_by(sampleName, gene, contig_name, closest_sequence))


## Mapping statistics

# Summarize statistics from the bowtie2 mapping

# Breadth of coverage without duplicates:
# Input files has the following format: <sampleName>.<reference>.markdup.csv
# And the content is:
# sampleName,reference,cov_breadth_min_1,cov_breadth_min_5,cov_breadth_min_10,avg_depth
# ROV2,NODE_19_length_792_cov_86.636090,100,100,100,,222.98175787728027
# The problem here is that I need to know also the contig name in the info above.

# List files
depth_files <- list.files(path = path_4, pattern = "markdup.csv$", full.names = TRUE)

# Empty df
depth_df <- read_csv(depth_files) %>%
    rename(contig_name = reference) %>%
    # Remove everything after "." in the contig name
    separate(contig_name, into = c("contig_name", NA), sep = "\\.")

# Join
joined_df <- full_join(joined_df, depth_df, by = join_by(sampleName, contig_name))

## Mapped reads with duplicates

# List files
stats_withdup_files <- list.files(path = path_5, pattern = ".*withdup.*\\.csv$", full.names = TRUE)

stats_withdup_df <- read_csv(stats_withdup_files) %>%
  select(sampleName,
         #gene,
         "contig_name" = contig,
         "trimmed_reads_withdups_mapped" = mapped_reads)

# Join
#joined_df <- full_join(joined_df, stats_withdup_df, by = join_by(sampleName, gene, contig_name))
joined_df <- full_join(joined_df, stats_withdup_df, by = join_by(sampleName, contig_name))


## Mapped reads without duplicates

# List files
stats_markdup_files <- list.files(path = path_6, pattern = ".*markdup.*\\.csv$", full.names = TRUE)

stats_markdup_df <- read_csv(stats_markdup_files) %>%
  select(sampleName,
         #gene,
         "contig_name" = contig,
         "trimmed_reads_nodups_mapped" = mapped_reads)

# Join
#joined_df <- full_join(joined_df, stats_markdup_df, by = join_by(sampleName, gene, contig_name))
joined_df <- full_join(joined_df, stats_markdup_df, by = join_by(sampleName, contig_name))

# Cutadapt ----------------------------------------------------------------

# List files
cutadapt_files <- list.files(path = path_7, pattern = ".log$", full.names = TRUE)

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

# Join
joined_df <- full_join(joined_df, cutadapt_df, by = join_by(sampleName))

# Kraken ------------------------------------------------------------------
# List files
kraken_files <- list.files(path = path_8, pattern = "kraken2.report.txt$", full.names = TRUE)

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

# Join
joined_df <- full_join(joined_df, kraken_df, by = join_by(sampleName))


# # Summarize data per sample. For the moment keep one sample per row.
# parse_phylogeny_df <- parse_phylogeny_df %>%
#   # Check how many genes/segments present per sample
#   group_by(sampleName, gene) %>%
#   summarise(
#     gene_count = n(), # Count the number of the same gene per sample
#     target_clade = paste(target_clade, collapse = ";"), # Make a new column with all target clades per gene
#     closest_sequence = paste(closest_sequence, collapse = ";"), # Make a new column with all closest sequences per gene
#     ratio = paste(ratio, collapse = ";") # Make a new column with all ratios per gene
#   ) %>%
#   ungroup()


# Reorder columns
final <- joined_df %>%
  select(sampleName,
         total_raw_reads,
         total_trimmed_reads,
         total_classified_reads,
         gene,
         contig_name,
         "clade" = target_clade,
         ratio,
         "closest_reference" = closest_sequence,
         percent_similarity,
         aligned_length,
         total_length,
         cov_breadth_min_1,
         cov_breadth_min_5,
         cov_breadth_min_10,
         "average_depth" = avg_depth,
         trimmed_reads_withdups_mapped,
         trimmed_reads_nodups_mapped)

# Write file
write_csv(final, file = "Genotype_mapping_summary_long.csv")

# Write file for MultiQC
# Add MultiQC info lines
header <- c("# id: 'summary'",
            "# section_name: 'Summary'",
            "# description: 'These statistics are generated from the process SUMMARIZE_ROV and the R script summarize_rov.R",
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


# # Create LW import --------------------------------------------------------

# # NEED TO ADD:
# #  - "Majority quality:" og "Minor quality" (typbar/ikke typbar)

# lw_import <- final %>%
#   select("Sample" = sampleName,
#          "Percent mapped reads of trimmed:" = Percent_reads_mapped_of_trimmed_with_dups_major,
#          "Majority genotype:" = Major_genotype_mapping,
#          "Number of mapped reads:" = Reads_withdup_mapped_major,
#          "Percent covered:" = Major_cov_breadth_min_1,
#          "Number of mapped reads without duplicates:" = Reads_nodup_mapped_major,
#          "Average depth without duplicates:" = Major_avg_depth,
#          "Percent covered above depth=5 without duplicates:" = Major_cov_breadth_min_5,
#          "Percent covered above depth=9 without duplicates:" = Major_cov_breadth_min_10,
#          "Most abundant minority genotype" = Minor_genotype_mapping,
#          "Percent most abundant minority genotype:" = abundance_minor,
#          "Number of mapped reads minor:" = Reads_withdup_mapped_minor,
#          "Percent covered minor:" = Minor_cov_breadth_min_5,
#          "Number of mapped reads minor without duplicates:" = Reads_nodup_mapped_minor,
#          "Average depth minor without duplicates:" = Minor_avg_depth,
#          "Percent covered above depth=5 minor without duplicates:" = Minor_cov_breadth_min_5,
#          "Percent covered above depth=9 minor without duplicates:" = Minor_cov_breadth_min_10,
#          "Script name and stringency:" = script_name_stringency,
#          "Total number of reads before trim:" = total_raw_reads,
#          "Total number of reads after trim:" = total_trimmed_reads,
#          "Majority quality:" = major_typbar,
#          "Minor quality:" = minor_typbar,
#          everything()
#   ) %>%
#   select(-`Most abundant minority genotype`,
#          -total_classified_reads,
#          -Major_reference,
#          -Minor_reference,
#          -abundance_major,
#          -Percent_reads_mapped_of_trimmed_with_dups_minor,
#          -Reads_nodup_mapped_first_mapping,
#          -Minor_cov_breadth_min_1
#   )

# # Remove column "Major_minor" if exists
# # This column is not present if GLUE is dropped
# if ("^Major_minor$" %in% colnames(lw_import)) {
#   lw_import <- lw_import %>% select(-Major_minor)
# }

# Write file
#write_csv(lw_import, file = "Genotype_mapping_summary_long_LW_import.csv")


