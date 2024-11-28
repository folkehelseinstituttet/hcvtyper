#!/usr/bin/env python3

from Bio import SeqIO
import sys

def extract_sequence(input_fasta, sequence_name):
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == sequence_name:
                output_fasta = f"{sequence_name}.fasta"
                with open(output_fasta, "w") as outfile:
                    SeqIO.write(record, outfile, "fasta")
                print(f"Sequence {sequence_name} has been written to {output_fasta}")
                return
        print(f"Sequence {sequence_name} not found in {input_fasta}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: prepare_bowtie2_build.py <input_fasta> <sequence_name>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    sequence_name = sys.argv[2]

    extract_sequence(input_fasta, sequence_name)
