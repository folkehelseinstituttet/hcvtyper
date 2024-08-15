#!/usr/bin/env python3

from Bio import Phylo
import csv
import sys

# Define sample name
sample_name = sys.argv[1]

# Define input
tree_file = sys.argv[2]
gene_name = sys.argv[3]

# Extract the gene name from the input file name. Compare this to the gene_name variable as a sanity check
gene_name_2 = tree_file.split(".")[1]

# Define the search pattern to identify the sample branch in the tree
target_name = 'NODE'  # Part of the name of your target sequence. Should be unique among all sequences

def list_and_sort_by_auto_prefix_save_csv(tree_file, target_name):
    tree = Phylo.read(tree_file, "newick")

    # Find all target clades
    target_clades = [clade for clade in tree.find_clades() if clade.name and target_name in clade.name]

    if not target_clades:
        print(f"Target '{target_name}' not found.")
        return

    # Prepare the csv file for writing
    csv_file_path = '{}.parse_phylo.{}.csv'.format(sample_name, gene_name)

    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write header
        writer.writerow(["Sequence", "Target clade prefix", "Closest sequence", "Closest Count / Total Count"])

        for target_clade in target_clades:
            # Calculate distances from all clades to the target
            distances = []
            for clade in tree.find_clades():
                if clade.name and clade != target_clade:
                    distance = tree.distance(target_clade, clade)
                    distances.append((clade.name, distance))

            # Sort by distance
            sorted_distances = sorted(distances, key=lambda x: x[1])

            if sorted_distances:
                # Process the closest header to stop at the first underscore
                closest_header = sorted_distances[0][0]
                underscore_index = closest_header.find('_')
                if underscore_index != -1:
                    closest_header = closest_header[:underscore_index]  # Truncate at the first underscore

            # Filter by this dynamically found prefix
            prefix = closest_header.split('|')[0] if '|' in closest_header else closest_header
            filtered_names = [name for name, _ in sorted_distances if name.startswith(prefix)]

            # Count total sequences with the dynamically determined prefix
            total_with_prefix = sum(1 for clade in tree.find_clades() if clade.name and clade.name.startswith(prefix))

            # Calculate the ratio of the closest count divided by the total count
            ratio = len(filtered_names) / total_with_prefix if total_with_prefix else 0

            # Append the sequence, ptarget clade name, closest sequence, and ratio to the CSV file
            writer.writerow([target_clade.name, prefix, closest_header, ratio])

    print(f"Data saved to {csv_file_path}")

if __name__ == "__main__":
    if len(sys.argv) > 1 and gene_name == gene_name_2:
        list_and_sort_by_auto_prefix_save_csv(tree_file, target_name)
    else:
        print("Please provide the gene name parameter.")
