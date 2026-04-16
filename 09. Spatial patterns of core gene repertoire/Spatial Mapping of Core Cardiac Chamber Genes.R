rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

library(ggplot2)
library(Seurat)


base_path <- "./"
setwd(base_path)

# Parameters
bin <- 100
size <- 0.4
size_merge <- 0.4

# ==============================================================================
# 2. File and Data Preparation
# ==============================================================================

# Select RDS files based on bin size
if (bin == 50) {
  rds_dir <- "./01bin50"
} else if (bin == 100) {
  rds_dir <- "./01bin100"
}

# Maintain specific file order from original script
all_rds_files <- list.files(rds_dir, pattern = 'rds$', full.names = TRUE)
file_idx <- c(2, 1, 11, 4, 8, 6, 5, 10, 3, 14, 13, 12, 9, 7)
files <- all_rds_files[file_idx]

# Load Core Gene Lists
atrgene <- read.csv(file.path(base_path, "Atrium_core_genes.csv"))
vtrgene <- read.csv(file.path(base_path, "ventricle_core_genes.csv"))

# Clean Gene Names
atrgene$Mmu <- gsub('gene-', '', atrgene$Mmu)
atrgene$Hsa <- toupper(atrgene$X)
vtrgene$Mmu <- gsub('gene-', '', vtrgene$Mmu)
vtrgene$Hsa <- toupper(vtrgene$X)

# Define Species Vector
specisse_list <- c("Afu", "Afa", "Pye", "Cin", "Lre", "Dre", "Cpl", 
                   "Pan", "Ame", "Xla", "Tsc", "Tgu", "Mmu", "Hsa")

# Define Color Palette
merge_colors <- c("atr" = "#2e84c6", "vtr" = "#b42d45", "atr_vtr" = "#d8d4cf", "others" = "#d8d4cf")

# ==============================================================================
# 3. Main Processing Loop
# ==============================================================================

for (i in seq_along(files)) {
  message("--- Processing file: ", files[i], " ---")
  
  # Load Seurat object
  sc <- readRDS(files[i])
  specise <- tools::toTitleCase(specisse_list[i])
  message("Species: ", specise)
  
  # 3.1 Scoring ----------------------------------------------------------------
  
  # Atrial Scoring
  atr_features <- atrgene[[specise]]
  sc <- AddModuleScore(object = sc, features = list(atr_features), name = "atr")
  
  captured_atr <- sum(atr_features %in% rownames(sc))
  message("Atrial genes: total ", length(atr_features), ", captured ", captured_atr)
  
  # Ventricular Scoring
  vtr_features <- vtrgene[[specise]]
  sc <- AddModuleScore(object = sc, features = list(vtr_features), name = "vtr")
  
  captured_vtr <- sum(vtr_features %in% rownames(sc))
  message("Ventricular genes: total ", length(vtr_features), ", captured ", captured_vtr)
  
  # 3.2 Thresholding-------------------------------------------------
  
  # Atrial threshold
  atr_threshold <- quantile(sc@meta.data$atr1, 0.9)
  sc$chambers_atr <- ifelse(sc$atr1 > atr_threshold, "atr", "others")
  
  # Ventricular threshold
  vtr_threshold <- quantile(sc$vtr1, 0.9)
  sc$chambers_vtr <- ifelse(sc$vtr1 > vtr_threshold, "vtr", "others")
  
  # Define Merge column
  sc$merge <- ifelse(sc$chambers_atr != "others" & sc$chambers_vtr != "others", "atr_vtr",
                ifelse(sc$chambers_atr != "others", "atr",
                  ifelse(sc$chambers_vtr != "others", "vtr", "others")))
  
  metadata <- sc@meta.data
  
  # 3.3 Plotting ---------------------------------------------------------------
  
  plot_and_save <- function(sub_metadata, col_name, folder_name, plot_size, suffix) {
    # Base background layer
    p <- ggplot(sub_metadata, aes(x = x, y = y)) +
      geom_tile(fill = "#d8d4cf") +
      theme_void() +
      coord_fixed(ratio = 1) +
      theme(legend.position = "none")
    
    # Add colored points layer by layer
    unique_vals <- unique(sub_metadata[[col_name]])
    for (val in unique_vals) {
      # Specific logic for i=4 (Cin) preserved from original script
      if (i == 4 && col_name %in% c("chambers_vtr", "merge") && val == "others") next
      
      p <- p + geom_point(
        data = sub_metadata[sub_metadata[[col_name]] == val, ],
        aes(x = x, y = y),
        color = merge_colors[val],
        size = plot_size
      )
    }
    
    # Directory and Save
    out_dir <- file.path(base_path, "01Image", folder_name, as.character(plot_size))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, paste0(specise, suffix, ".pdf")), plot = p, device = "pdf")
  }
  
  # Execute Plots
  plot_and_save(metadata, "chambers_atr", "01Atrium", size)
  plot_and_save(metadata, "chambers_vtr", "02ventricle",size)
  plot_and_save(metadata, "merge", ", size_merge)
    
  message("Done with ", specise, "\n")
  rm(sc, metadata, sc_atr, sc_vtr)
  gc()
}