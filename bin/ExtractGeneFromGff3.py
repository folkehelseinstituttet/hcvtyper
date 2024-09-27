#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

# Define input files
sample_id = sys.argv[1]
gff3      = sys.argv[2]
fasta     = sys.argv[3]

def extract_sequences(name, thresholds):
    try:
        # Set the thresholds for the current gene, or use default values if not specified
        lower_threshold, upper_threshold = thresholds.get(name, (100, 10000))  # Default thresholds

        with open(gff3) as gff_file, open(fasta) as fasta_file:
            sequences = SeqIO.parse(fasta_file, "fasta")
            fasta_dict = {seq.id: seq for seq in sequences}

            seq_diffs = []
            for line in gff_file:
                columns = line.strip().split("\t")
                if len(columns) < 9:
                    continue

                if "Name=" + name in columns[8]:
                    seqid = columns[0]
                    start_pos = int(columns[3]) - 1  # Convert to 0-based index
                    end_pos = int(columns[4])
                    strand = columns[6]

                    if seqid in fasta_dict:
                        full_sequence = fasta_dict[seqid]
                        gene_sequence = full_sequence.seq[start_pos:end_pos]  # Slice to get gene sequence

                        if strand == "-":
                            gene_sequence = gene_sequence.reverse_complement()

                        new_seq_record = SeqRecord(gene_sequence,
                                                   id=full_sequence.id,
                                                   description=full_sequence.description.replace(seqid + " ", "", 1))

                        seq_diffs.append(new_seq_record)
                    else:
                        print(f"Sequence ID {seqid} not found in FASTA file")

            # Filter sequences based on the lower and upper threshold and sort by length
            seq_diffs = [seq for seq in seq_diffs if lower_threshold <= len(seq.seq) <= upper_threshold]
            seq_diffs.sort(key=lambda x: len(x.seq), reverse=True)

            # Save all sequences in a single file if there are any sequences that meet the criteria
            # Modify output file name
            if seq_diffs:
                output_file = f"{sample_id}.{name}.fasta"
                with open(output_file, "w") as out_file:
                    SeqIO.write(seq_diffs, out_file, "fasta")

    except FileNotFoundError:
        print("One of the files was not found. Please check the file paths.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Create a dictionary with gene names and length thresholds with lower and upper bounds
        thresholds = {
            "VP7": (490, 1300),
            "VP4": (1100, 2500),
            "VP6": (590, 1500),
            "VP1": (3000, 3500),
            "VP2": (2300, 2800),
            "VP3": (2200, 2700),
            "NSP1": (1200, 1600),
            "NSP2": (800, 1200),
            "NSP3": (800, 1200),
            "NSP4": (500, 900),
            "NSP5": (500, 750)
            # Add more gene names and their thresholds as needed
        }
        # Loop through all genes and extract sequences
        for name in thresholds.keys():
            extract_sequences(name, thresholds)
    else:
        print("Please provide the gene name parameter.")
