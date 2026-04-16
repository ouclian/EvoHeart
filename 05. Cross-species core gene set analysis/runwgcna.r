library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

args = commandArgs(T)

seurat_obj <- readRDS(args[1])
seurat_obj <- NormalizeData(seurat_obj) %>% FindVariableFeatures(nfeatures=3000) %>% ScaleData() %>% RunPCA(verbose=FALSE)

#meta = read.csv(args[2],row.names=1)
#seurat_obj@meta.data = meta
celltype = read.csv(args[2])
seurat_obj$celltype = celltype$celltype[match(seurat_obj$RNA_snn_res.0.5,celltype$cluster)]
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'hdWGCNA'
)

# set up dummy variable so we can run MetacellsByGroups for all clusters together
seurat_obj$all_cells <- 'all'

# run hdWGCNA metacell aggregation
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("orig.ident","celltype"),
  ident.group = "celltype",
  k = 25,
  target_metacells=1000,
  max_shared=5
)

seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)

pdf(paste0("metacell_celltype.pdf"))
DimPlotMetacells(seurat_obj, group.by='celltype') + umap_theme() + ggtitle("Cell Type")
dev.off()
pdf(paste0("metacell_orig.ident.pdf"))
DimPlotMetacells(seurat_obj, group.by='orig.ident') + umap_theme() + ggtitle("Sample")
dev.off()

# setup expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  use_metacells=TRUE,
  assay = 'RNA',
  slot = 'data'
)


# test soft power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
plot_list <- PlotSoftPowers(seurat_obj)
pdf(paste0("softpower_threshold.pdf"))
wrap_plots(plot_list, ncol=2)
dev.off()


# compute the co-expression network
seurat_obj <- ConstructNetwork(seurat_obj)

# compute module eigengenes and eigengene-based connectivity
seurat_obj <- ModuleEigengenes(seurat_obj,group.by.vars="orig.ident")
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
hMEs <- GetMEs(seurat_obj)
saveRDS(seurat_obj, file=paste0('hdWGCNA_object.rds'))

seurat_obj <- ModuleConnectivity(seurat_obj)
pdf(paste0("KME.genes.pdf"))
PlotKMEs(seurat_obj, ncol=20)
dev.off()
# rename modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = 'hd-M',
  wgcna_name='hdWGCNA'
)

# plot module eigengenes
plot_list <- ModuleFeaturePlot(seurat_obj, order=TRUE, raster=TRUE, alpha=1, restrict_range=FALSE,features='hMEs')

# assemble plots
pdf(paste0("hub_gene.pdf"))
wrap_plots(plot_list, ncol=6)
dev.off()
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module);
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
write.table(modules,"module_genes.csv",sep="\t",quote=F)

hub_df <- GetHubGenes(seurat_obj, n_hubs = 500)
write.table(hub_df,"hub_genes.head500.csv",sep="\t",quote=F)
