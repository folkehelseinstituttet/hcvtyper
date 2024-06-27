#!/usr/bin/env python3

from Bio import SeqIO
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py <gene_name>")
    sys.exit(1)

gene_name = sys.argv[1]

def extract_and_combine_sequences(header_file, references_file, output_type_file, temp_output_file):
    """
    Extracts sequences from reference and type files based on modified headers and combines them into a single output file.
    """
    # Step 1: Read the header from the CSV file
    try:
        with open(header_file, 'r') as file:
            full_header = file.readline().strip().split(':')[0]  # Consider part before ':' only
    except FileNotFoundError:
        print(f"Error: The file {header_file} was not found.")
        return

    # Step 2: Read references.fasta and extract the matching sequence
    matching_sequence = None
    try:
        for record in SeqIO.parse(references_file, "fasta"):
            # Compare only the part of the header before ':'
            if record.description.split(':')[0] == full_header:
                matching_sequence = record
                break
    except FileNotFoundError:
        print(f"Error: The file {references_file} was not found.")
        return

    if not matching_sequence:
        print("No matching sequence found in the reference file.")
        return

    # Step 3: Write the matching sequence to a new file
    try:
        with open(temp_output_file, 'w') as output_handle:
            SeqIO.write([matching_sequence], output_handle, "fasta")
    except Exception as e:
        print(f"Error writing to {temp_output_file}: {str(e)}")
        return

    # Step 4: Append the sequence from output_x_1.fasta
    try:
        type_sequence = next(SeqIO.parse(output_type_file, "fasta"), None)
        if type_sequence:
            with open(temp_output_file, 'a') as output_handle:
                SeqIO.write([type_sequence], output_handle, "fasta")
        else:
            print("No sequence found in the type output file.")
    except FileNotFoundError:
        print(f"Error: The file {output_type_file} was not found.")

# Define file paths
csv_file = 'temp_prefix_header_{}.csv'.format(gene_name)
references_file = '/home/torstein/RoV/RotavAr/v0/6/References_{}.fasta'.format(gene_name)
output_type_file = 'output_{}_1.fasta'.format(gene_name)
temp_output_file = 'temp_percentagecalc_{}.fasta'.format(gene_name)

# Call the function with correct parameters
extract_and_combine_sequences(csv_file, references_file, output_type_file, temp_output_file)
