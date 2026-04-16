#!/bin/bash
# ==================================================
# Script name: run_orthofinder.v1.0.sh
# Author:      Wentao Han(hanwentao@ouc.edu.cn)
# Date created: 2024-03-05
# Last modified: 2024-03-05
# Version:     1.0
# ==================================================

# Step 1: Infer orthologs using OrthoFinder

orthofinder -f pep -t 40 -S blast -M msa

# Step 2: Construct a phylogenetic tree using RAxML

raxmlHPC-SSE3 -s aligned_sequences.fasta -m PROTGAMMAILGF -o Nve -n species_tree_root_20240305 -T 48 -N 1000 -p 20240305 -f a -x 20240305