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
csv_ratio          = sys.argv[4]

# Extract the gene name from the input fasta and csv files. Compare this to the gene_name variable as a sanity check
gene_name_2 = aligned_fasta_file.split(".")[1]
gene_name_3 = csv_ratio.split("_")[-1].split("." )[0]

def calculate_percent_similarity_and_update_csv(fasta_file, prefix_csv):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    if len(sequences) < 2:
        print("Error: The FASTA file must contain at least two sequences.")
        return

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

    # Load existing data from prefix_ratio.csv
    with open(prefix_csv, 'r') as file:
        reader = csv.reader(file)
        existing_data = list(reader)

    # Output to CSV
    csv_file = '{}.genotyping_result.{}.csv'.format(sample_name, gene_name)

    with open(csv_file, "w", newline='') as file:
        writer = csv.writer(file)
        # Write the existing headers and append new headers
        writer.writerow(existing_data[0] + ["Percentage Similarity", "Aligned Length", "Total Length"])
        # Write the existing data and append the new data
        writer.writerow(existing_data[1] + [round(percent_similarity, 2), aligned_length, len(seq1)])

    print(f"Percentage similarity saved to {csv_file}")

# Run the function
if gene_name == gene_name_2 == gene_name_3:
    calculate_percent_similarity_and_update_csv(aligned_fasta_file, csv_ratio)
