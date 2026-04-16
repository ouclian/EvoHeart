rm(list = ls())
options(stringsAsFactors = FALSE)
gc()
setwd('.')

# --- Load Libraries ---
library(ggplot2)
library(Seurat)
library(dplyr)

# ==============================================================================
# Function: plot_c_scores
# Description: Calculates signature scores using AddModuleScore and plots them spatially
# ==============================================================================
plot_c_scores <- function(i, c1 = TRUE, c2 = FALSE, c3 = FALSE, c4 = TRUE, 
                          output = "./02ouput_result/",
                          prefix = "c1_c4_bind",
                          bin = 50,
                          addblast = TRUE,
                          savepdf = FALSE,
                          separation = TRUE,
                          size = 0.5) {
  
  # --- 1. File Path Configuration ---
  if (bin == 50) {
    rds_dir <- './00data/01bin50RDS/00rds/'
  } else {
    rds_dir <- './00data/02bin100RDS/02.bin100rds/'
  }
  
  dict_dir <- ifelse(addblast, "./00data/c1-c4-add-dict", "./00data/c1-c4-dict")
  
  # List files and reorder to maintain specific species sequence
  files <- list.files(rds_dir, pattern = 'rds$', full.names = TRUE)
  files <- files[c(2, 1, 11, 4, 8, 6, 5, 10, 3, 14, 13, 12, 9, 7)]
  
  dirs <- list.dirs(dict_dir, full.names = TRUE)
  
  # Logging current file
  print(paste0("Reading RDS: ", files[i - 1]))
  
  # --- 2. Data Loading ---
  clas_files <- list.files(dirs[i], full.names = TRUE)
  sc <- readRDS(files[i - 1])
  
  # Selected signatures based on function arguments
  selected_cs <- c("C1", "C2", "C3", "C4")[c(c1, c2, c3, c4)]
  selected_files <- list()
 
  # Read gene lists for each selected signature
  for (sig in selected_cs) {
    # Match filename based on substring pattern "C1", "C2", etc.
    match_idx <- grep(sig, substr(basename(clas_files), 8, 9))
    if (length(match_idx) > 0) {
      selected_files[[sig]] <- read.csv(clas_files[match_idx])
      print(paste0("Reading table: ", clas_files[match_idx]))
    }
  }
  
  # Extract species name from filename for labeling
  c1_idx <- grep("C1", substr(basename(clas_files), 8, 9))
  species <- substr(basename(clas_files[c1_idx]), 4, 6)

  # --- 3. Gene List Pre-processing ---
  # Standardize gene names and filter genes present in the Seurat object
  process_gene_list <- function(sig_name) {
    if (!is.null(selected_files[[sig_name]]) && nrow(selected_files[[sig_name]]) != 0) {
      if (i == 12) {
        selected_files[[sig_name]][, 5] <<- gsub("gene-", "", selected_files[[sig_name]][, 5])
      }
      genes <- selected_files[[sig_name]][selected_files[[sig_name]][, 5] %in% rownames(sc), 5]
      return(genes)
    }
    return(NULL)
  }

  # --- 4. Module Scoring ---
  if (separation) {
    # Calculate scores for each signature independently
    for (sig in selected_cs) {
      genes <- process_gene_list(sig)
      if (length(genes) >= 2) {
        print(paste0("Begin scoring: ", sig))
        sc <- AddModuleScore(object = sc, features = list(genes), name = tolower(sig))
        
        # Binary classification based on 90th percentile threshold
        score_col <- paste0(tolower(sig), "1")
        threshold <- quantile(sc@meta.data[[score_col]], 0.9)
        sc@meta.data[[paste0("chambers_", tolower(sig))]] <- ifelse(
          sc@meta.data[[score_col]] > threshold, tolower(sig), "others"
        )
      } else {
        print(paste0("Signature ", sig, ": Insufficient features"))
      }
    }
  } else {
    # Merge all signature gene lists into one combined score
    all_genes <- unique(unlist(lapply(selected_cs, process_gene_list)))
    print("Begin merged scoring")
    merged_name <- paste0(tolower(selected_cs), collapse = "_")
    sc <- AddModuleScore(object = sc, features = list(all_genes), name = merged_name)
    
    score_col <- paste0(merged_name, "1")
    threshold <- quantile(sc@meta.data[[score_col]], 0.9)
    sc@meta.data[[paste0("chambers_", merged_name)]] <- ifelse(
      sc@meta.data[[score_col]] > threshold, merged_name, "others"
    )
  }

  # --- 5. Result Labeling & Merging ---
  # Create output directory structure
  addb_suffix <- ifelse(addblast, "addblast", "no_addblast")
  final_output_path <- file.path(output, paste0("bin", bin, "_", addb_suffix), prefix)
  dir.create(final_output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Handle coordinate mapping for specific species (e.g., Ciona)
  if (i == 5) {
    sc$x <- sc$array_col
    sc$y <- sc$array_row
  }
  
  metadata <- sc@meta.data
  
  # Determine labels for cells based on which signatures they passed the threshold for
  if (separation) {
    condition_matrix <- sapply(selected_cs, function(sig) metadata[[paste0("chambers_", tolower(sig))]] != "others")
    merge_labels <- apply(condition_matrix, 1, function(row) {
      if (!any(row)) return("others")
      paste0(selected_cs[row], collapse = "_")
    })
  } else {
    merge_labels <- metadata[[paste0("chambers_", merged_name)]]
  }
  metadata$merge <- merge_labels
  
  # --- 6. Visualization ---
  # Define standard color palette for signature combinations
  color_palette <- c(
    "C1" = "red", "C2" = "green", "C3" = "purple", "C4" = "blue", 
    "C1_C2" = "orange", "C1_C3" = "pink", "C1_C4" = "yellow", 
    "C2_C3" = "cyan", "C2_C4" = "magenta", "C3_C4" = "brown", 
    "C1_C2_C3" = "gray", "C1_C2_C4" = "lightblue", "C1_C3_C4" = "lightgreen", 
    "C2_C3_C4" = "lightpink", "C1_C2_C3_C4" = "black", "others" = "#D3D3D3"
  )
  
  # Base plot with background (others)
  p1 <- ggplot(metadata, aes(x = x, y = y)) +
    geom_tile(fill = "grey") +
    theme_void() +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none")

  # Overlay points for each signature combination
  unique_labels <- setdiff(unique(metadata$merge), "others")
  for (label in unique_labels) {
    p1 <- p1 + geom_point(
      data = subset(metadata, merge == label),
      aes(x = x, y = y),
      color = color_palette[label],
      size = size
    )
  }
  
  # Special point-based base plot for index 5
  if (i == 5) {
    p1 <- ggplot(metadata, aes(x = x, y = y, color = merge)) +
      geom_point(size = size) +
      scale_color_manual(values = color_palette) +
      theme_void() +
      coord_fixed(ratio = 1) +
      theme(legend.position = "none")
  }

  # --- 7. Save Output ---
  save_dir <- file.path(final_output_path, as.character(size))
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  file_base <- file.path(save_dir, paste0(species, "_", paste0(selected_cs, collapse = "_")))
  
  print(paste0("Saving: ", file_base, ".png"))
  ggsave(paste0(file_base, ".png"), plot = p1, width = 10, height = 10, dpi = 200)
  
  if (savepdf) {
    ggsave(paste0(file_base, ".pdf"), plot = p1, width = 10, height = 10)
  }

  # Cleanup memory
  rm(sc, metadata, p1)
  gc()
}

# ==============================================================================
# Execution Loop: Run Multiple Parameter Combinations
# ==============================================================================

# Signatures Glossary: 
# c1 = CoreDT, c2 = CorePT, c3 = CoreBA2, c4 = CoreBA1

sizes_to_test <- c(2.5, 0.5, 0.3, 0.1, 2)
bins_to_test <- c(50, 100)

for (sz in sizes_to_test) {
  for (bn in bins_to_test) {
    
    # 1. Single Signature Analysis
    # C2 only (index 2-4)
    for (idx in 2:4) {
      plot_c_scores(i = idx, c1 = F, c2 = T, c3 = F, c4 = F, prefix = "c2_CorePT", bin = bn, size = sz)
    }
    # C1 only (index 5-15)
    for (idx in 5:15) {
      plot_c_scores(i = idx, c1 = T, c2 = F, c3 = F, c4 = F, prefix = "c1_CoreDT", bin = bn, size = sz)
    }
    # C3 only (index 2-15)
    for (idx in 2:15) {
      plot_c_scores(i = idx, c1 = F, c2 = F, c3 = T, c4 = F, prefix = "c3_CoreBA2", bin = bn, size = sz)
    }
    # C4 only (index 2-15)
    for (idx in 2:15) {
      plot_c_scores(i = idx, c1 = F, c2 = F, c3 = F, c4 = T, prefix = "c4_CoreBA1", bin = bn, size = sz)
    }
    
    # 2. Combined Signature Analysis (Binding)
    # C1 + C4
    for (idx in 5:15) {
      plot_c_scores(i = idx, c1 = T, c2 = F, c3 = F, c4 = T, prefix = "CoreDT_CoreBA1_bind", bin = bn, size = sz, savepdf = T)
    }
    # C2 + C4
    for (idx in 2:4) {
      plot_c_scores(i = idx, c1 = F, c2 = T, c3 = F, c4 = T, prefix = "CorePT_CoreBA1_bind", bin = bn, size = sz, savepdf = T)
    }
  }
}