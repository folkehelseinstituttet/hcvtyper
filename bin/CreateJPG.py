#!/usr/bin/env python3

import os
from ete3 import Tree, TreeStyle, NodeStyle, AttrFace
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Required for headless environments

sample_name = sys.argv[1]
gene_name   = sys.argv[2]
treefile    = sys.argv[3]

# Extract the gene name from the input fasta and csv files. Compare this to the gene_name variable as a sanity check
gene_name_2 = treefile.split(".")[1]

def display_tree(tree_file):
    # Check if the tree file exists
    if not os.path.exists(tree_file):
        print(f"Error: Tree file {tree_file} not found.")
        return

    try:
        # Load the tree
        tree = Tree(tree_file)
        print("Tree file loaded successfully.")

        # Root the tree at a node with prefix 'OutgroupRVC'
        root_node = None
        for node in tree.traverse():
            if node.name.startswith('OutgroupRVC'):
                root_node = node
                break

        if root_node:
            tree.set_outgroup(root_node)
            print(f"Rooted the tree at node: {root_node.name}")
        else:
            print("No node with prefix 'OutgroupRVC' found. Tree will not be re-rooted.")

        # Set up tree style
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.show_branch_length = True
        ts.show_branch_support = True

        # Traverse each node, applying custom styles as necessary
        for node in tree.traverse():
            if node.is_leaf():
                if node.name == "target":
                    # Define a special face for the target node: red and bold
                    face = AttrFace("name", fgcolor="red", fstyle='bold')
                else:
                    # Regular style for other nodes
                    face = AttrFace("name")

                # Add the face to the node
                node.add_face(face, column=0, position="branch-right")

        # Render tree in JPG format
        tree.render(tree_jpg, w=1600, units="px", dpi=300, tree_style=ts)
        print(f"Tree image saved as {tree_jpg}")

    except Exception as e:
        print(f"Failed to render tree: {e}")

if __name__ == '__main__':
    if gene_name == gene_name_2:
        tree_jpg = f'{sample_name}.temp2_tree_visualization.{gene_name}.jpg'
        display_tree(treefile)
