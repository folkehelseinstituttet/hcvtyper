#!/usr/bin/env python3

from Bio import SeqIO
import csv
import sys

# Define sample name
sample_name = sys.argv[1]

# Define gene name
gene_name = sys.argv[2]

# Define file names
aligned_fasta_file = sys.argv[3]

# Extract the gene name from the input fasta and csv files. Compare this to the gene_name variable as a sanity check
gene_name_2 = aligned_fasta_file.split(".")[1]

def calculate_percent_similarity_and_update_csv(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    if len(sequences) < 2:
        print("Error: The FASTA file must contain at least two sequences.")
        return

    # Get the names of the two sequences calculate_pairwise_alignment_metrics
    ref_name    = sequences[0].name
    # Replace "|" and "/" with "_" in the full_header and keep everything before the first colon
    sanitized_ref_name = ref_name.replace("|", "_").replace("/", "_").split(':')[0]
    contig_name = sequences[1].name
    # Sanitize the sequence name to keep only the first two elements. E.g. "NODE_18"
    sanitized_contig_name = "_".join(contig_name.split("_")[:2])

    # Extract the two sequences
    seq1 = str(sequences[0].seq).upper()
    seq2 = str(sequences[1].seq).upper()

    # Find the start and end of the aligned regions, excluding leading and trailing gaps
    start = next((i for i in range(min(len(seq1), len(seq2))) if seq1[i] != '-' or seq2[i] != '-'), None)
    end = next((i for i in range(min(len(seq1), len(seq2))-1, -1, -1) if seq1[i] != '-' or seq2[i] != '-'), None)

    if start is None or end is None:
        print("No alignment found.")
        return

    # Slicing the sequences to the aligned region
    seq1 = seq1[start:end+1]
    seq2 = seq2[start:end+1]

    # Calculate the percent similarity and aligned length excluding internal leading and trailing gaps
    aligned_length = sum(1 for nt1, nt2 in zip(seq1, seq2) if nt1 != '-' and nt2 != '-')
    differences = sum(nt1 != nt2 for nt1, nt2 in zip(seq1, seq2) if nt1 != '-' and nt2 != '-')
    percent_similarity = ((1 - (differences / aligned_length)) * 100) if aligned_length else 0

    # Output to CSV
    csv_file = f"{sample_name}.alignment_metrics.{gene_name}.{sanitized_contig_name}.{sanitized_ref_name}.csv"

    with open(csv_file, "w", newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(["Sequence", "Closest sequence", "Percent Similarity", "Aligned Length", "Total Length"])
        # Aappend the new data
        writer.writerow([contig_name, ref_name, round(percent_similarity, 2), aligned_length, len(seq1)])

    print(f"Percentage similarity saved to {csv_file}")

# Run the function
if gene_name == gene_name_2:
    calculate_percent_similarity_and_update_csv(aligned_fasta_file)
