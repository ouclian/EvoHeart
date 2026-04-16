#!/bin/bash
# ==================================================
# Script name: run_saturn.v1.0.sh
# Author:      Qun Liu(liuqun@genomics.com)
# Date created: 2025-04-14
# Last modified: 2025-04-14
# Version:     1.0
# ==================================================

source activate saturn
python3 train-saturn.py --in_data=run_heart.csv \
                              --in_label_col=celltype \
			      --ref_label_col=celltype \
                              --num_macrogenes=2000  \
			      --hv_genes=8000  \
                              --centroids_init_path=run_heart.pkl \
                              --work_dir=. \
                              --device cpu \
