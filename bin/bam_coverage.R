#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: bam_coverage.R <depth file> <prefix>", call.=FALSE)
}

depth <- args[1]
prefix <- args[2]
# Get the sample name, reference, major_minor from the depth filename.
# Use this to control against the prefix and info in bam file.
sampleName <- unlist(str_split(basename(depth), pattern = "\\."))[1]
major_minor <- unlist(str_split(basename(depth), pattern = "\\."))[3]
ref_from_filename <- unlist(str_split(basename(depth), pattern = "\\."))[2] %>%
    str_remove(stringr::fixed(paste0("_", major_minor)))

df <- read_tsv(depth, col_names = FALSE)

# Check that prefix and sampleName are identical. If not, stop the script and print an error message
if (prefix != sampleName) {
  stop("Error: Prefix and sample name do not match. Please check the input files.", call.=FALSE)
}

# If there are no mapped reads df is empty
if (nrow(df) > 0) {
  df <- df %>%
    # Rename columns
    rename("Reference" = X1,
           "Position" = X2,
           "Coverage" = X3)

# Get the mapped reference from the tsv-file produced by samtools depth
reference <- df %>% distinct(Reference) %>% pull(Reference)

# Check that ref_from_filename and reference are identical. If not, stop the script and print an error message
if (ref_from_filename != reference) {
  stop("Error: Reference from filename and reference from depth file do not match. Please check the input files.", call.=FALSE)
}

# Define the coverage cutoff
coverage_cutoff <- 10

# Add column for highlighting
df <- df %>%
  mutate(below_cutoff = ifelse(Coverage < coverage_cutoff, Coverage, NA_real_))

# Start base plot
plot <- ggplot(df, aes(x = Position)) +
  # Add ribbon only if some values are below cutoff
  {if (any(df$Coverage < coverage_cutoff, na.rm = TRUE))
    geom_ribbon(aes(ymin = 0, ymax = below_cutoff), fill = "red", alpha = 0.2)
  } +

  # Main coverage line
  geom_line(aes(y = Coverage), color = "steelblue", linewidth = 0.6) +

  # Horizontal line at cutoff
  geom_hline(yintercept = coverage_cutoff, color = "darkgreen", linetype = "dotted") +

  # Label cutoff line
  annotate("text", x = max(df$Position, na.rm = TRUE) * 0.01,
         y = -max(df$Coverage, na.rm = TRUE) * 0.02,  # 2% below the baseline
         label = paste("Coverage cutoff =", coverage_cutoff),
         hjust = 0, color = "darkgreen", size = 3.5) +

  # Labels
  labs(
    x = "Genome Position",
    y = "Read Coverage",
    title = paste0(sampleName, ".", major_minor, ".", reference)
  ) +

  # Theme tweaks
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save the plot
ggsave(plot = plot,
       file = paste0(sampleName, ".", major_minor, ".", reference, ".png"),
         device = "png",
         dpi = 300,
         bg = "white")
}
