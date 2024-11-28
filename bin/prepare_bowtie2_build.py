#!/usr/bin/env python3

from Bio import SeqIO
import sys

def extract_node_sequences(input_fasta, prefix, gene_name):
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            # The de novo assembled contig always starts with "NODE"
            if record.id.startswith("NODE"):
                # Sanitize the sequence name to keep only the first two elements
                # sanitized_contig_name = "_".join(record.id.split("_")[:2])
                output_fasta = f"{prefix}.{gene_name}.{record.id}.fasta"
                with open(output_fasta, "w") as outfile:
                    SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: prepare_bowtie2_build.py <input_fasta> <prefix> <gene_name>")
        sys.exit(1)

    fasta     = sys.argv[1]
    prefix    = sys.argv[2]
    gene_name = sys.argv[3]

    extract_node_sequences(fasta, prefix, gene_name)
