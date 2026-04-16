# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# ==================================================
# Script name: run_samap.v1.0.py
# Author:      Qun Liu(liuqun@genomics.com)
# Date created: 2025-04-14
# Last modified: 2025-04-14
# Version:     1.0
# ==================================================

import argparse as arg
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                                sankey_plot, chord_plot, CellTypeTriangles, 
                                ParalogSubstitutions, FunctionalEnrichment,
                                convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pickle

parser = arg.ArgumentParser(formatter_class=arg.RawDescriptionHelpFormatter, description = '''This programe is used to perform SAMap.''')
parser.add_argument('-specieslist', action='store', required = "True", help='the first species name of blast')
parser.add_argument('-fnlist', action='store', required = "True",help='the h5ad of 1st species')
parser.add_argument('-celltypelist', action='store', required = "True",help='Celltype1.')
parser.add_argument('-NUMITERS', type=int, default=3, action='store', help='the number of iters')
parser.add_argument('-cpu', type=int, default=10, help='Number of cpu used.')
parser.add_argument('-fmap', required = "True",  help='the blast output directory')
parser.add_argument('-outprefix',required = "True",  help='the prefix of output.')
#parser.add_argument('-mappingtable', action='store', help='Mappingtable path.')
#parser.add_argument('-enrichedGene', action='store', help='Enriched genes table')
args = parser.parse_args()

specieslist = args.specieslist.split(",")
fnlist = args.fnlist.split(",")
celltypelist = args.celltypelist.split(",")

filenames = dict(zip(specieslist, fnlist))
keys = dict(zip(specieslist,celltypelist))

#filenames = {args.species1:args.fn1, args.species2:args.fn2, args.species3:args.fn3}
#keys={args.species1:args.celltype1, args.species2:args.celltype2, args.species3:args.celltype3}
sm = SAMAP(
    filenames,
    f_maps = args.fmap,
    save_processed=False, #if False, do not save the processed results to `*_pr.h5ad`
    keys=keys
   )
#neigh_from_keys={args.species1:True, args.species2:True}
#sm.run(NUMITERS=args.NUMITERS, neigh_from_keys=neigh_from_keys, ncpus=args.cpu)
sm.run(NUMITERS=args.NUMITERS, ncpus=args.cpu)

with open('all.sm', 'wb') as f:
    pickle.dump(sm, f)

D,MappingTable=get_mapping_scores(sm, keys, n_top=0)
MappingTable.to_csv(args.outprefix + "mappingtable.xls", sep='\t')

plt.figure(figsize=(100,150))
sc.pl.umap(sm.smap.samap.adata, color='species')
plt.savefig('umap_species.png', bbox_inches='tight', dpi=600)

# plot one specie umap
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.sams[args.species1].adata, color=args.celltype1)
#plt.savefig('umap_'+ args.species1+ '.png', bbox_inches='tight', dpi=600)
#
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.sams[args.species2].adata, color=args.celltype2)
#plt.savefig('umap_'+ args.species2+ '.png', bbox_inches='tight', dpi=600)

# plot one specie umap with samap matrix
#sm.sams[args.species1].adata.obsm['X_umap']=sm.sams[args.species1].adata.obsm['X_umap_samap']
#sm.sams[args.species2].adata.obsm['X_umap']=sm.sams[args.species2].adata.obsm['X_umap_samap']
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.sams[args.species1].adata, color=args.celltype1)
#plt.savefig('umap_'+ args.species1+ '_samap.png', bbox_inches='tight', dpi=600)
#
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.sams[args.species2].adata, color=args.celltype2)
#plt.savefig('umap_'+ args.species2+ '_samap.png', bbox_inches='tight', dpi=600)

# gene-pair finder
gpf=GenePairFinder(sm,keys=keys)
gene_pairs=gpf.find_all(align_thr=0.1, thr=0.05)
gene_pairs.to_csv(args.outprefix + "enrichedGene.xls", sep='\t', header=True, index=None)

#sankey_plot
plt.figure(figsize=(100,150))
sankey_plot(MappingTable, align_thr=0.05, species_order = [args.species1,args.species2])
plt.savefig(args.outprefix + '_sankey.png', bbox_inches='tight', dpi=600)

## plot umap for all
#celltype=[args.species1+ '_'+ i for i in sm.smap.samap.adata.obs[args.species1+ '_'+ args.celltype1].to_list() if i !='unassigned'] + [args.species2+ '_'+ i for i in sm.smap.samap.adata.obs[args.species2+ '_'+ args.celltype2].to_list() if i !='unassigned']
#sm.smap.samap.adata.obs['CellType']=celltype
#plt.figure(figsize=(100,150))
#sc.pl.umap(sm.smap.samap.adata, color='CellType')
#plt.savefig('umap_'+ args.species1+ '_'+ args.species2+ '_celltype.png', bbox_inches='tight', dpi=600)


