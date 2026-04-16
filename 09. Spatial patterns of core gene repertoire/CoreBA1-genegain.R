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
# ==============================================================================
plot_c_scores <- function(i, c1 = FALSE, c2 = FALSE, c3 = FALSE, c4 = TRUE, 
                          genegain = TRUE, lend = TRUE,
                          output = "./09genegain_C4/00result/",
                          prefix = "gain_c4_bind",
                          bin = 100,
                          addblast = TRUE,
                          savepdf = FALSE,
                          size = 0.5,
                          separation = TRUE) {
  
  # --- 1. Path and File Selection ---
  # Define base directories based on parameters
  rds_base <- if (bin == 50) "./00data/01bin50RDS/00rds/" else "./00data/02bin100RDS/02.bin100rds/"
  dict_base <- if (addblast) "./00data/c1-c4-add-dict" else "./00data/c1-c4-dict"
  
  files <- list.files(rds_base, pattern = 'rds$', full.names = TRUE)
  
  # Apply specific species ordering as per original logic
  if (bin == 50 || addblast) {
    files <- files[c(2, 1, 11, 4, 8, 6, 5, 10, 3, 14, 13, 12, 9, 7)]
  } else if (bin == 100 && !addblast) {
    files <- files[c(2, 1, 11, 4, 8, 6, 5, 10, 3, 14, 13, 12, 9, 7)]
  }
  
  dirs <- list.dirs(dict_base, full.names = TRUE)
  if (genegain) dirs_gain <- list.dirs("./07gain-gene", full.names = TRUE)

  print(paste0("Processing RDS: ", files[i - 1]))
  
  # --- 2. Load Data ---
  clas_files <- list.files(dirs[i], full.names = TRUE)
  sc <- readRDS(files[i - 1])
  
  selected_cs <- c("C1", "C2", "C3", "C4")[c(c1, c2, c3, c4)]
  selected_files <- list()
  
  # Load signature gene lists (Column 5)
  for (sig in selected_cs) {
    match_idx <- grep(sig, substr(basename(clas_files), 8, 9))
    if (length(match_idx) > 0) {
      temp_df <- read.csv(clas_files[match_idx])
      gene_names <- temp_df[, 5]
      # Species-specific gene name cleaning
      if (i == 12) gene_names <- gsub("gene-", "", gene_names)
      # Filter genes present in the dataset
      selected_files[[sig]] <- gene_names[gene_names %in% rownames(sc)]
    }
  }
  
  # Load Gene Gain list (Column 6)
  if (genegain) {
    gain_file <- list.files(dirs_gain[i], full.names = TRUE)
    if (length(gain_file) > 0) {
      temp_gain <- read.csv(gain_file)
      gain_genes <- temp_gain[, 6]
      # Species-specific prefix/suffix handling
      if (i == 5) gain_genes <- paste0("KY21:", gain_genes)
      if (i == 12) gain_genes <- gsub("gene-", "", gain_genes)
      selected_files[["gain"]] <- gain_genes[gain_genes %in% rownames(sc)]
    }
  }

  species <- substr(basename(clas_files[grep("C4", substr(basename(clas_files), 8, 9))]), 4, 6)

  # --- 3. Module Scoring ---
  all_targets <- if (genegain) c(selected_cs, "gain") else selected_cs
  
  if (separation) {
    # Score each signature separately
    for (sig in all_targets) {
      if (length(selected_files[[sig]]) >= 2) {
        print(paste0("Scoring signature: ", sig))
        sc <- AddModuleScore(object = sc, features = list(selected_files[[sig]]), name = tolower(sig))
        
        # Calculate threshold (90th percentile)
        score_col <- paste0(tolower(sig), "1")
        threshold <- quantile(sc@meta.data[[score_col]], 0.9)
        sc@meta.data[[paste0("chambers_", tolower(sig))]] <- ifelse(
          sc@meta.data[[score_col]] > threshold, tolower(sig), "others"
        )
      }
    }
  } else {
    # Merge all genes into a single score
    pooled_genes <- unique(unlist(selected_files))
    merged_name <- paste0(tolower(selected_cs), collapse = "_")
    sc <- AddModuleScore(object = sc, features = list(pooled_genes), name = merged_name)
    
    score_col <- paste0(merged_name, "1")
    threshold <- quantile(sc@meta.data[[score_col]], 0.9)
    sc@meta.data[[paste0("chambers_", merged_name)]] <- ifelse(
      sc@meta.data[[score_col]] > threshold, merged_name, "others"
    )
  }

  # --- 4. Result Path Setup ---
  addb_flag <- ifelse(addblast, "addblast", "no_addblast")
  len_flag <- ifelse(lend, "lend", "no_glend")
  final_output <- file.path(output, paste0("bin", bin, "_", addb_flag), prefix, len_flag)
  dir.create(final_output, showWarnings = FALSE, recursive = TRUE)

  # Coordinate mapping for specific species (Ciona)
  if (i == 5) {
    sc$x <- sc$array_col
    sc$y <- sc$array_row
  }
  
  metadata <- sc@meta.data
  
  # Combine individual scores into one label per cell
  condition_cols <- if (separation) all_targets else merged_name
  merge_conditions <- sapply(condition_cols, function(x) metadata[[paste0("chambers_", tolower(x))]] != "others")
  
  metadata$merge <- apply(merge_conditions, 1, function(row) {
    if (!any(row)) return("others")
    paste0(all_targets[row], collapse = "_")
  })
  
  # --- 5. Visualization Setup ---
  # Define complex color palette for all possible combinations
  merge_colors <- c(
    "C1" = "#FF4500", "C2" = "#32CD32", "C3" = "#800080", "C4" = "blue", "gain" = "red",
    "C1_C2" = "#FF6347", "C1_C3" = "#FF69B4", "C1_C4" = "#FFD700", "C1_gain" = "#87CEEB",
    "C2_C3" = "#20B2AA", "C2_C4" = "#00FA9A", "C2_gain" = "#7FFFD4", "C3_C4" = "#A52A2A",
    "C3_gain" = "#BA55D3", "C4_gain" = "#D3D3D3", "others" = "#D3D3D3"
  )

  # Plotting logic
  if (lend) {
    # Standard plot with legend
    p1 <- ggplot(metadata, aes(x = x, y = y, color = merge)) +
      geom_point(size = size) +
      scale_color_manual(values = merge_colors) +
      theme_void() + coord_fixed(ratio = 1) +
      guides(color = guide_legend(title = "Chambers", override.aes = list(size = 2)))
  } else {
    # Manual layered plot without legend
    p1 <- ggplot(metadata, aes(x = x, y = y)) +
      geom_tile(fill = "grey90") + # Background
      theme_void() + coord_fixed(ratio = 1) +
      theme(legend.position = "none")
    
    # Layer each group
    unique_labels <- setdiff(unique(metadata$merge), "others")
    for (lbl in unique_labels) {
      p1 <- p1 + geom_point(data = subset(metadata, merge == lbl), 
                            aes(x = x, y = y), color = merge_colors[lbl], size = size)
    }
  }

  # --- 6. Save Output ---
  if (savepdf) {
    pdf_dir <- file.path(final_output, "pdf", paste0("size_", size))
    dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
    pdf_path <- file.path(pdf_dir, paste0(species, "_", paste0(all_targets, collapse = "_"), "_gain.pdf"))
    print(paste0("Saving PDF: ", pdf_path))
    ggsave(pdf_path, plot = p1, width = 10, height = 10)
  }
  
  # Memory cleanup
  rm(sc, metadata, p1)
  gc()
  print("Done.")
}

# ==============================================================================
# Execution Section: Loop through parameters
# ==============================================================================

# Target: C4 and Gain genes
for (sz in c(1, 1.5, 2, 2.5, 3, 5)) {
  for (bn in c(50, 100)) {
    for (idx in 2:14) {
      plot_c_scores(i = idx, 
                    c1 = FALSE, c2 = FALSE, c3 = FALSE, c4 = TRUE, 
                    genegain = TRUE, 
                    output = "./09genegain_C4/00result/",
                    prefix = "C4_gain_bind",
                    bin = bn,
                    addblast = FALSE,
                    savepdf = TRUE,
                    size = sz,
                    separation = TRUE,
                    lend = FALSE)
    }
  }
}