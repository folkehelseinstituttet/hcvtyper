#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_bam_variation.R <bam file> <prefix>", call.=FALSE)
}

library(Rsamtools)
library(tidyverse)

# Define the bam file
bam_file <- args[1]
prefix <- args[2]

# Read the sample name from the input file. Use this as a control that the right bam file is inputted
sampleName <- unlist(str_split(basename(bam_file), pattern = "\\."))[1]

# Check that prefix and sampleName are identical. If not, stop the script and print an error message
if (prefix != sampleName) {
  stop("Error: Prefix and sample name do not match. Please check the input files.", call.=FALSE)
}

# Extract if it's a major or minor variant from the bam file name
major_minor <- unlist(str_split(basename(bam_file), pattern = "\\."))[3]

# Set variation threshold for highlighting (modify as needed)
variation_threshold <- 0.15

# Set minimum coverage for considering a position
min_coverage <- 10

# Open the BAM file
bam <- BamFile(bam_file)

# Read the BAM header
header <- scanBamHeader(bam_file)

# Extract and view the reference lengths
ref_length <- header[[1]]$targets
ref_name <- names(header[[1]]$targets)

# Set up pileup parameters: here we distinguish nucleotides and ignore strand information for simplicity.
p_param <- PileupParam(max_depth = 1000000,
                       min_nucleotide_depth = 0,
                       distinguish_nucleotides = TRUE,
                       distinguish_strands = FALSE,
                       ignore_query_Ns = TRUE)

# Generate the pileup data
# The result is a data frame with columns such as: seqnames, pos, nucleotide, and count.
pileup_result <- pileup(bam, pileupParam = p_param)


# Calculate total counts and percentage per position
pileup_summary <- pileup_result %>%
  group_by(seqnames, pos) %>% # Group by reference and position
  mutate(total = sum(count)) %>%               # Total count at that position
  mutate(percentage = (count / total) * 100) %>% # Percentage of each nucleotide
  ungroup()

# Create a complete grid of positions and nucleotides (A, C, G, T)
all_positions <- expand.grid(
  seqnames = ref_name,
  pos = 1:ref_length,
  nucleotide = c("A", "C", "G", "T"),
  stringsAsFactors = FALSE
)

# Merge the complete grid with the pileup result
complete_pileup <- all_positions %>%
  left_join(pileup_result, by = c("seqnames", "pos", "nucleotide")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# Compute variation per position
variation_data <- complete_pileup %>%
  group_by(seqnames, pos) %>%
  summarise(
    total_count = sum(count),
    max_nucleotide_count = max(count, na.rm = TRUE),
    variation = ifelse(total_count > 0, 1 - (max_nucleotide_count / total_count), 0)
  ) %>%
  ungroup()

# Identify positions with zero coverage
zero_coverage <- variation_data %>% filter(total_count == 0)

# Identify positions with total count of 9 or less
low_coverage <- variation_data %>% filter(total_count > 0 & total_count < min_coverage)

# Set variation to zero for low coverage positions
variation_data <- variation_data %>%
  mutate(variation = ifelse(total_count > 0 & total_count <= 9, 0, variation))

# Identify positions with variation above the threshold
high_variation <- variation_data %>% filter(variation > variation_threshold)

# Plot variation per position
p <- ggplot(variation_data, aes(x = pos, y = variation)) +
  geom_segment(aes(xend = pos, yend = 0), color = "black") +  # Black bars for variation
  geom_point(data = zero_coverage, aes(x = pos, y = 0), color = "blue", size = 2) +  # Blue dots for zero coverage
  geom_point(data = low_coverage, aes(x = pos, y = 0), color = "grey", size = 2) +  # Grey dots for low coverage
  geom_point(data = high_variation, aes(x = pos, y = variation), color = "red", size = 2) +  # Red dots for high variation
  labs(title = paste0(sampleName, ".", major_minor,".", ref_name, ". Variation threshold ", variation_threshold, ". Min. coverage ", min_coverage),
       x = "Position",
       y = "Nucleotide variation") +
  ylim(0, 1.0) +  # Hardcode the y-axis limit to 1.0
  theme_minimal()

# Save the plot as a PNG file
output_filename <- paste0(sampleName,".variation_plot_", ref_name, "_", major_minor, ".png")
ggsave(output_filename, plot = p, width = 10, height = 5, dpi = 300, bg = "white")

# Write the raw variation data to a tsv file
output_data_filename <- paste0(sampleName,".bam_variation_", ref_name, "_", major_minor, ".tsv")
write_tsv(variation_data, output_data_filename)

# Write the complete pileup data to a tsv file
output_pileup_filename <- paste0(sampleName,".bam_pileup_", ref_name, "_", major_minor, ".tsv")
write_tsv(complete_pileup, output_pileup_filename)
