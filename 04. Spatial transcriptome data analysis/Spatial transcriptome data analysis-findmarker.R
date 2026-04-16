# ==============================================================================
# Environment Setup
# ==============================================================================
rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

library(Seurat)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(clustree)
library(SeuratWrappers)

# Setup Conda Environment
use_condaenv("./SAMap", required = TRUE)
py_config()
print(py_config())
reticulate::py_last_error()

# Load Gene ID Dictionary
id_dict <- read.table("./geneid_dict.csv")
id_dict$Tsc <- gsub("gene-", "", id_dict$Tsc)

# List input files
files <- list.files("./spc_cellbin/00spc", pattern = 'h5ad$', full.names = TRUE)

# ==============================================================================
# Main Processing Loop (per file/species)
# ==============================================================================
for (t in seq_along(files)) {
  print(paste("Processing file", t, ":", files[t]))
  
  filename <- files[t]
  spcise <- substr(basename(filename), 1, 3)
  
  # Load and convert AnnData to Seurat
  adata <- sceasy::convertFormat(filename, from = "anndata", to = "seurat")
  
  # Clean up gene names (replace underscores with dashes for Seurat compatibility)
  rownames(adata@assays$RNA@data) <- gsub("_", "-", rownames(adata@assays$RNA@data))
  rownames(adata@assays$RNA@meta.features) <- gsub("_", "-", rownames(adata@assays$RNA@meta.features))
  
  # Standard Seurat Workflow
  adata <- NormalizeData(adata, verbose = FALSE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc.filt@var.genes, npcs = 30, verbose = FALSE)
  
  adata <- FindNeighbors(adata, reduction = "pca", dims = 1:30)
  

  res_range <- seq(0.1, 1.5, 0.1)
  adata <- FindClusters(adata, resolution = res_range)
  

  p_tree <- clustree(adata, prefix = "RNA_snn_res.") +
    theme(legend.position = "bottom") + 
    scale_color_brewer(palette = "Set1") +
    scale_edge_color_continuous(low = "grey80", high = "red")
  
  out_base <- file.path("./01_spc_cellbin_deg", spcise)
  dir.create(file.path(out_base, "01Dendrogram"), recursive = TRUE, showWarnings = FALSE)
  
  print("Saving Dendrogram...")
  ggsave(file.path(out_base, "01Dendrogram.png"), p_tree)
  
  # 2. Load EggNOG annotations (Moved outside inner loop for efficiency) ------
  eggnog_path <- file.path("./Eggnog", spcise, "eggnog-result.emapper.annotations")
  egg_anno <- read.table(eggnog_path, sep = "\t", header = TRUE, comment.char = "#")
  
  egg_cols <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", 
                "max_annot_lvl", "COG_category", "Description", "Preferred_name", 
                "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", 
                "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", 
                "BiGG_Reaction", "PFAMs")
  colnames(egg_anno) <- egg_cols
  
  # Clean and handle duplicate Preferred names
  egg_anno <- egg_anno[, c("query", "Preferred_name")]
  egg_anno$Preferred_name <- make.unique(as.character(egg_anno$Preferred_name), sep = "_")
  
  # 3. Process each resolution ------------------------------------------------
  for (i in res_range) {
    res_col <- paste0("RNA_snn_res.", i)
    Idents(adata) <- res_col
    
    # Calculate markers
    cluster_markers <- SeuratWrappers::RunPrestoAll(
      adata, test.use = "wilcox", only.pos = TRUE, 
      logfc.threshold = 0.25, min.pct = 0.25
    )
    
    # Save Spatial Plot
    spatial_dir <- file.path(out_base, "02Spatial Graph")
    dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
    
    p_spatial <- ggplot(adata@meta.data, aes(x, y, color = !!sym(res_col))) +
      geom_point(size = 0.2)
    
    ggsave(file.path(spatial_dir, paste0(spcise, "_", res_col, "_空间.png")), p_spatial)
    
    # Save UMAP Plot
    umap_dir <- file.path(out_base, "03umap")
    dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
    
    p_umap <- DimPlot(adata, label = TRUE)
    ggsave(file.path(umap_dir, paste0(spcise, "_", res_col, "_umap.png")), p_umap)
    
    # Annotation mapping: Gene ID Symbol
    idx_id <- match(cluster_markers$gene, id_dict[[spcise]])
    cluster_markers$symbol <- id_dict$genesymbol[idx_id]
    
    # Annotation mapping: EggNOG Symbol
    query_genes <- cluster_markers$gene
    if (spcise == 'Tsc') {
      query_genes <- paste0('gene-', query_genes)
    }
    
    idx_egg <- match(query_genes, egg_anno$query)
    cluster_markers$egg_symbol <- egg_anno$Preferred_name[idx_egg]
    
    # Save Marker CSV
    print(paste("Saving marker file for resolution", i))
    write.csv(cluster_markers, file.path(out_base, paste0(spcise, "_", res_col, "_marker.csv")))
    
    # Highlight individual clusters in spatial plots
    for (q in unique(adata@meta.data[[res_col]])) {
      highlight_cluster <- as.character(q)
      
      # Prepare metadata for highlighting
      adata@meta.data$highlight <- ifelse(adata@meta.data[[res_col]] == highlight_cluster, highlight_cluster, "Other")
      adata@meta.data$highlight <- factor(adata@meta.data$highlight, levels = c("Other", highlight_cluster))
      
      p_highlight <- ggplot(adata@meta.data, aes(x, y, color = highlight)) +
        geom_point(size = 0.2) +
        scale_color_manual(values = setNames(c("gray", "red"), c("Other", highlight_cluster))) +
        labs(color = "Cluster") +
        theme_minimal()
      
      # Save individual cluster plot
      cluster_plot_dir <- file.path(spatial_dir, spcise, as.character(i))
      dir.create(cluster_plot_dir, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(cluster_plot_dir, paste0(highlight_cluster, ".png")), p_highlight)
    }
  }
  
  # 4. Save Final AnnData -----------------------------------------------------
  print('Saving processed H5AD...') 
  out_h5ad <- file.path("./spc_cellbin/02h5ad", paste0(basename(filename), "_spc_cellbin.h5ad"))
  dir.create(dirname(out_h5ad), recursive = TRUE, showWarnings = FALSE)
  
  sceasy::convertFormat(
    adata, from = "seurat", to = "anndata", 
    drop_single_values = FALSE, outFile = out_h5ad
  )
}