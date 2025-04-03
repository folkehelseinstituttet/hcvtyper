#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: plot_bam_variation.R <bam file>", call.=FALSE)
}

library(Rsamtools)
library(tidyverse)

# Set noise threshold for highlighting (modify as needed)
noise_threshold <- 0.15  

bam_file <- args[1]

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

# Compute noise per position
noise_data <- complete_pileup %>%
  group_by(seqnames, pos) %>%
  summarise(
    total_count = sum(count),
    max_nucleotide_count = max(count, na.rm = TRUE),
    noise = ifelse(total_count > 0, 1 - (max_nucleotide_count / total_count), 0)
  ) %>%
  ungroup()

# Identify positions with zero coverage
zero_coverage <- noise_data %>% filter(total_count == 0)

# Identify positions with noise above the threshold
high_noise <- noise_data %>% filter(noise > noise_threshold)

# Plot noise per position
p <- ggplot(noise_data, aes(x = pos, y = noise)) +
  geom_segment(aes(xend = pos, yend = 0), color = "black") +  # Black bars for variation
  geom_point(data = zero_coverage, aes(x = pos, y = 0), color = "blue", size = 2) +  # Blue dots for zero coverage
  geom_point(data = high_noise, aes(x = pos, y = noise), color = "red", size = 2) +  # Red dots for high noise
  labs(title = paste("Variation per Position in", ref_name),
       x = "Position",
       y = "Noise") +
  theme_minimal()

# Save the plot as a PNG file
output_filename <- paste0("variation_plot_", ref_name, "_major.png")
ggsave(output_filename, plot = p, width = 10, height = 5, dpi = 300, bg = "white")
