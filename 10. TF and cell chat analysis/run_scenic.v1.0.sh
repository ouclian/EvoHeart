#!/bin/bash
# ==================================================
# Script name: run_scenic.v1.0.sh
# Author:      Qun Liu(liuqun@genomics.com)
# Date created: 2025-04-14
# Last modified: 2025-04-14
# Version:     1.0
# ==================================================


source activate scenic
in_loom=$1
outpre=${in_loom%%.*}


pyscenic grn \
--num_workers 20 \
--output ${outpre}".grn.tsv" \
--method grnboost2 \
${in_loom} \
allTFs_hg38.txt

pyscenic ctx \
${outpre}".grn.tsv" \
hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ${in_loom} \
--mode "dask_multiprocessing" \
--output ${outpre}".reg.csv" \
--num_workers 20 \
--mask_dropouts

pyscenic aucell \
${in_loom} \
${outpre}".reg.csv" \
--output ${outpre}".SCENIC.loom" \
--num_workers 20

