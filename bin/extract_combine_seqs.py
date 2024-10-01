#!/usr/bin/env python3

from Bio import SeqIO
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py <gene_name>")
    sys.exit(1)

# Define sample name
sample_name       = sys.argv[1]

# Define gene name
gene_name         = sys.argv[2]

# Define files
references_file   = sys.argv[3]
gff_extract_fasta = sys.argv[4]
csv_file          = sys.argv[5]

#csv_file         = '{}.parse_phylo.{}.csv'.format(sample_name, gene_name)
#references_file = '/home/torstein/RoV/RotavAr/v0/6/References_{}.fasta'.format(gene_name)
#output_type_file = 'output_{}_1.fasta'.format(gene_name)
#output_file  = '{}.temp_closest_ref.{}.fasta'.format(sample_name, gene_name)

# Extract the gene name from the input csv file. Compare this to the gene_name variable as a sanity check
gene_name_2 = csv_file.split(".")[-2]

def extract_and_combine_sequences(csv_file, references_file, gff_extract_fasta):
    """
    Extracts sequences from reference and type files based on modified headers and combines them into a single output file.
    """
    # Step 1: Read the contents of the csv file from FindGT6_2.py
    try:
        with open(csv_file, 'r') as file:
            lines = file.readlines()[1:]  # Skip the header
        # Debugging
        print(f"Read {len(lines)} lines from {csv_file}")
    except FileNotFoundError:
        print(f"Error: The file {csv_file} was not found.")
        return

    # Step 2: Process each line in the csv file, find the closest sequence and extract it from the reference fasta file
    for line in lines:
        columns = line.strip().split(',')
        if len(columns) < 3:
            print(f"Invalid line format: {line}")
            continue

        # Extract the third column that holds the name of the nearest reference sequence to the de novo sequence
        full_header = columns[2]

        # Replace "|" and "/" with "_" in the full_header
        sanitized_header = full_header.replace("|", "_").replace("/", "_")

        # Extract the sequence name from the first column
        sequence_name = columns[0]

        # Sanitize the sequence name to keep only the first two elements. E.g. "NODE_18"
        sanitized_sequence_name = "_".join(sequence_name.split("_")[:2])

        # Debugging
        print(f"Processing full_header: {full_header}, sequence_name: {sequence_name}")

        # Find the corresponding sequence in the references file
        reference_sequence = None
        try:
            for record in SeqIO.parse(references_file, "fasta"):
                # Compare only the part of the header before ':'
                if record.description.split(':')[0] == full_header:
                    reference_sequence = record
                    break
            # Debugging
            if reference_sequence:
                print(f"Found reference sequence for {full_header}")
            else:
                print(f"No matching reference sequence found for {full_header}")
        except FileNotFoundError:
            print(f"Error: The file {references_file} was not found.")
            return

        if not reference_sequence:
            print("No matching sequence found in the reference file.")
            return

        # Find the corresponding sequence in the gff_extract_fasta file
        contig_sequence = None
        try:
            for record in SeqIO.parse(gff_extract_fasta, "fasta"):
                if record.id == sequence_name:
                    contig_sequence = record
                    break
            # Debugging
            if contig_sequence:
                print(f"Found contig sequence for {sequence_name}")
            else:
                print(f"No matching contig sequence found for {sequence_name}")
        except FileNotFoundError:
            print(f"Error: The file {gff_extract_fasta} was not found.")
            return

        if not contig_sequence:
            print(f"No matching sequence found in the gff_extract_fasta file for {sequence_name}.")
            continue

        # Step 3: Write the matching sequence to a new file
        combined_output_file = f"{sample_name}.{gene_name}.{sanitized_sequence_name}.{sanitized_header}.fasta"
        # Debugging
        print(f"Writing combined sequences to {combined_output_file}")
        try:
            with open(combined_output_file, 'w') as output_handle:
                SeqIO.write([reference_sequence, contig_sequence], output_handle, "fasta")
            # Debugging
            print(f"Successfully wrote to {combined_output_file}")
        except Exception as e:
            print(f"Error writing to {combined_output_file}: {str(e)}")
            return

# Check for gene name mismatch and exit with an error message if there is a mismatch
if gene_name != gene_name_2:
    print(f"Error: Gene name mismatch: {gene_name} != {gene_name_2}")
    sys.exit(1)

# Call the function with correct parameters
if gene_name == gene_name_2:
    extract_and_combine_sequences(csv_file, references_file, gff_extract_fasta)
else:
    print(f"Gene name mismatch: {gene_name} != {gene_name_2}")
