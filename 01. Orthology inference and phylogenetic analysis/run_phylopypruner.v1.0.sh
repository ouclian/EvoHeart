#!/bin/bash
# ==================================================
# Script name: run_phylopypruner.v1.0.sh
# Author:      Wentao Han(hanwentao@ouc.edu.cn)
# Date created: 2024-03-05
# Last modified: 2024-03-05
# Version:     1.0
# ==================================================

# Set the working directory containing Orthogroup sequences
WORK_DIR="./"
PHYLOPYPRUNER_PATH="/mnt/inspurfs/home/hanwentao/miniconda3/bin/phylopypruner"

# Change to the working directory
cd "$WORK_DIR" || { echo "Directory not found: $WORK_DIR"; exit 1; }

# Step 1: Perform multiple sequence alignment on all FASTA files (*.fa)
echo "Starting multiple sequence alignment..."
for i in *.fa; do 
    # Check if there are any .fa files
    if [[ -f $i ]]; then
        mafft --maxiterate 1000 --localpair "$i" > "${i%%.fa}.fas"
        echo "Aligned $i to ${i%%.fa}.fas"
    else
        echo "No FASTA files found for alignment."
    fi
done

# Step 2: Build a phylogenetic tree for each aligned sequence file (*.fas) using FastTree
echo "Constructing phylogenetic trees..."
for i in *.fas; do 
    # Check if there are any .fas files
    if [[ -f $i ]]; then
        FastTree "$i" > "${i%%.fas}.nwk"
        echo "Constructed tree for $i and saved as ${i%%.fas}.nwk"
    else
        echo "No aligned FASTA files found for tree construction."
    fi
done

# Step 3: Run phylopypruner to prune the phylogenetic trees based on specified parameters
echo "Pruning phylogenetic trees in all OG directories..."
for OG_DIR in /mnt/inspurfs/home/hanwentao/heart/01_species_tree/OrthoFinder/*; do
    if [[ -d $OG_DIR ]]; then
        echo "Processing directory: $OG_DIR"
        "$PHYLOPYPRUNER_PATH" --dir "$OG_DIR" --min-len 100 --prune MI
    else
        echo "No orthogroup directories found."
    fi
done

echo "Pipeline completed successfully."
