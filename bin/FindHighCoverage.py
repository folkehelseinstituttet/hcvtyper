#!/usr/bin/env python3


import sys
from Bio import SeqIO

# Define input files
sample_id = sys.argv[1]
fasta     = sys.argv[2]

# Define gene names to search for
gene_names = ["VP7", "VP4", "VP6", "VP1", "VP2", "VP3", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5"]

def extract_highest_coverage(input_file):
    highest_cov = 0
    best_record = None

    # Read sequences and find the one with the highest coverage
    for record in SeqIO.parse(input_file, "fasta"):
        # Extract coverage value from the header
        header_parts = record.description.split('_')
        try:
            cov_index = header_parts.index('cov') + 1
            coverage = float(header_parts[cov_index])
            # Update if this record has higher coverage
            if coverage > highest_cov:
                highest_cov = coverage
                best_record = record
        except (ValueError, IndexError):
            continue  # Skip if the header is not in the expected format or coverage is missing

    if best_record:
        # Calculate sequence length
        seq_length = len(best_record.seq)
        # Create new header with sequence length
        new_header = f">{best_record.description} NewLength{seq_length}"
        best_record.description = best_record.id = "target"  # Reset description and id

        # Write the best record to the output file
        output_file = f"{sample_id}_{gene}_highest_cov.fasta"
        with open(output_file, "w") as output_handle:
            SeqIO.write(best_record, output_handle, "fasta")

        # Write new header to the header file
        output_header_file = f"{sample_id}_{gene}_highest_cov_header.txt"
        with open(output_header_file, "w") as header_handle:
            header_handle.write(new_header + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python script.py <gene>")
        sys.exit(1)
    for gene_name in gene_names:
        extract_highest_coverage(fasta)


