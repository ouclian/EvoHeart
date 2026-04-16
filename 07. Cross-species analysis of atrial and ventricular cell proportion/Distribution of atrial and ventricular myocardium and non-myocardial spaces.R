rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

library(Seurat)
library(sceasy)
library(reticulate)
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)
library(plotly)
args <- commandArgs(T)

setwd('/.')

# List all h5ad files in the input directory
files <- list.files("./01h5ad/", pattern = 'h5ad$', full.names = T)

# ==============================================================================
# Main Processing Loop (per file/species)
# ==============================================================================
for (q in 1:length(files)) {
  filename <- files[q]
  spcise <- substr(basename(filename), 1, 3)
  
  print("reading data")
  print(filename)
  adata <- sceasy::convertFormat(filename, from = "anndata", to = "seurat")

  # Remove red blood cells from both atria and ventricle
  adata <- subset(adata, !(Chamber_celltypes %in% c("Atria_Red blood cell", "Ventricle_Red blood cell")))
  table(adata$Chamber_celltypes)
  
  # Extract metadata from the seurat object
  metadata <- adata@meta.data
  metadata$merge <- ifelse(
    (metadata$Chamber_celltypes == "Atria_Cardiomyocyte"), 
    "Atria_Cardiomyocyte",
    ifelse(
      (metadata$Chamber_celltypes == "Ventricle_Cardiomyocyte"), 
      "Ventricle_Cardiomyocyte", 
      paste0(metadata$Chamber_separation_all, "_", "others")
    )
  )

  if (q == 4) {
    metadata$merge <- ifelse(
      (metadata$Chamber_celltypes == "Heart tube_Cardiomyocyte"),
      "Heart tube_Cardiomyocyte",
      "other"
    )
  }
   
  # Define color mapping for different cell types
  merge_colors <- c(
    "Atria_Cardiomyocyte" = "#005B99",          
    "Ventricle_Cardiomyocyte" = "#d252c4",         
    "Atria_others" = "#FFA500",         
    "Ventricle_others" = "#14f5d6",         
    "Other_others" = "grey",    
    "others" = "#14f5d6",     
    "Heart tube_Cardiomyocyte" = "#d252c4"    
  )
  
  # Determine point size based on file index
  size <- ifelse(q == 6, 1.5, ifelse(q == 1, 0.3, 0.01))

  # Create spatial distribution plot with no legend
  p1 <- ggplot()
  for (merge_value in unique(metadata$merge)) {
    p1 <- p1 + geom_point(
      data = subset(metadata, merge == merge_value),
      aes(x = x, y = y),
      color = merge_colors[merge_value],   
      size = size
    )
  }  
  p1 <- p1 + theme_void() +
    theme(legend.position = "right") 
  

  dir.create('./01Spatial_distribution/', showWarnings = F)
  dir.create('./01Spatial_distribution/01Legend-free/', showWarnings = F)
  
  ggsave(
    paste0("./01Spatial_distribution/01Legend-free/", spcise, "_Spatial_pattern_diagram.png"), 
    p1, 
    width = 10, 
    height = 10
  )
  ggsave(
    paste0("./01Spatial_distribution/01Legend-free/", spcise, "_Spatial_pattern_diagram.pdf"), 
    device = "pdf", 
    p1, 
    width = 10, 
    height = 10
  )

  # Create spatial distribution plot with legend
  p2 <- ggplot()  

  # Add points for each merge value with color mapping
  for (merge_value in unique(metadata$merge)) {
    p2 <- p2 + geom_point(
      data = subset(metadata, merge == merge_value),
      aes(x = x, y = y, color = merge), 
      size = size
    )
  }
   
  # Apply color scale and theme settings to the plot with legend
  p2 <- p2 + scale_color_manual(values = merge_colors) +
    labs(color = "celltype") + 
    theme_minimal()
  dir.create('./01Spatial_distribution/02Legend/', showWarnings = F)

  ggsave(
    paste0("./01Spatial_distribution/02Legend/", spcise, "_Spatial_pattern_diagram.png"), 
    p2
  )
  ggsave(
    paste0("./01Spatial_distribution/02Legend/", spcise, "_Spatial_pattern_diagram.pdf"), 
    device = "pdf", 
    p2
  )
    
  # Convert plot to interactive HTML version
  p1_plotly <- ggplotly(p1)
    dir.create('./01Spatial_distribution/03html/', showWarnings = F)
  htmlwidgets::saveWidget(
    p1_plotly, 
    paste0("./01Spatial_distribution/03html/", spcise, "_Spatial_pattern_diagram.html")
  )
}