rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working directory
wd <- "."
setwd(wd)

# --- Load Libraries ---
library(readr)
library(reticulate)
library(dplyr)

# --- Data Loading and Dictionary Setup ---
# Load gene ID mapping dictionary
id_map <- read.table("./05/geneid_dict.csv", header = TRUE, sep = ",")
id_map$Tsc <- gsub("gene-", "", id_map$Tsc)

# Load the merged count matrix (Large CSV)
raw_counts <- read_csv("./A_Vmerge.csv")
data_matrix <- raw_counts[, -1] 

# --- Python Integration: Aggregate rows by Geneid ---
# Sum up counts for identical Geneids using Pandas for efficiency
pd <- import("pandas")
py_data <- r_to_py(data_matrix)
py_data_agg <- py_data$groupby("Geneid")$sum()
data_aggregated <- py_to_r(py_data_agg)

# Add Geneid column back from rownames for subsetting
data_aggregated$Geneid <- rownames(data_aggregated)


# Function to generate formatted txt files for individual samples
# format: genesymbol, Geneid, Readcount, geneexonlength, RPKM
save_sample_txt <- function(sample_id, output_dir) {
  # Prepare a 5-column dataframe as required by downstream tools
  # Note: Readcount, length, and RPKM are set to the same count value in this logic
  sample_df <- data_aggregated[, c("Geneid", "Geneid", sample_id, sample_id, sample_id)]
  colnames(sample_df) <- c("genesymbol", "Geneid", "Readcount", "geneexonlength", "RPKM")
  
  # Round numerical columns
  sample_df[, 3:5] <- round(sample_df[, 3:5])
  
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Write to tab-separated file without headers
  write.table(sample_df, file.path(output_dir, paste0(sample_id, ".txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

# ==============================================================================
# Core Function: getdiff
# Purpose: Orchestrate sample file generation and execute AWK filtering for DEGs
# ==============================================================================

getdiff <- function(alist, vlist) {
  cat(paste0("\nProcessing Atrium list (size: ", length(alist), ") and Ventricle list (size: ", length(vlist), ")\n"))
  
  # 1. Export individual sample text files
  target_dir <- "./12.single_sample_deg/01a_v/"
  for (sample in c(alist, vlist)) {
    save_sample_txt(sample, target_dir)
  }
  
  # 2. Define directory paths for DEG results
  deg_input_dir <- "./12.single_sample_deg/02Single_sample_deg/00deg/"
  atrium_out_dir <- "./12.single_sample_deg/02Single_sample_deg/01Atrium/"
  ventricle_out_dir <- "./12.single_sample_deg/02Single_sample_deg/02ventricle/"
  
  dir.create(atrium_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(ventricle_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct the comparison filename prefix used by the external tool
  comp_name <- paste0(alist[1], "num", length(alist), "_VS_", vlist[1], "num", length(vlist))
  input_file <- paste0(deg_input_dir, comp_name, ".txt")
  
  # 3. Filter Atrium-high genes (log2FC > 1.5) using AWK
  cat("Screening for differentially expressed genes in Atria...\n")
  awk_atria <- sprintf(
    "awk '{if($3 > 1.5 && $4 == 1) print $1\",\"$2\",\"$3\",\"$4}' %s | sort -k3,3nr > %s",
    input_file, paste0(atrium_out_dir, comp_name, "_A.csv")
  )
  system(awk_atria)
  
  # 4. Filter Ventricle-high genes (log2FC < -1.5) using AWK
  cat("Screening for differentially expressed genes in Ventricles...\n")
  awk_ventricle <- sprintf(
    "awk '{if($3 < -1.5 && $4 == 1) print $1\",\"$2\",\"$3\",\"$4}' %s | sort -k3,3nr > %s",
    input_file, paste0(ventricle_out_dir, comp_name, "_V.csv")
  )
  system(awk_ventricle)
  
  # 5. Post-process CSVs: Add headers and re-save
  process_csv <- function(file_path) {
    if (file.exists(file_path) && file.info(file_path)$size > 0) {
      tmp_data <- read.table(file_path, quote = "#", sep = ",")
      colnames(tmp_data) <- c("Row.names", "Geneid_Ref", "log2FoldChange", "Flag")
      write.csv(tmp_data, file_path, row.names = FALSE)
    } else {
      warning(paste("File not found or empty:", file_path))
    }
  }
  
  process_csv(paste0(atrium_out_dir, comp_name, "_A.csv"))
  process_csv(paste0(ventricle_out_dir, comp_name, "_V.csv"))
}

# ==============================================================================
# Execution Section: Species Analysis
# ==============================================================================

# Species: 12_Cpl
getdiff(alist = "12_Cpl_A", vlist = "12_Cpl_V")

# Species: 13_Ler
getdiff(alist = "13_Ler_A", vlist = "13_Ler_V")

# Species: 15_Xtr
getdiff(alist = "15_Xtr_A", vlist = c("15_Xtr_V", "15_Xtr_V.1", "15_Xtr_V.2"))

# Species: 17_Gja
getdiff(alist = "17_Gja_A", vlist = c("17_Gja_V", "17_Gja_V.1", "17_Gja_V.2"))

cat("\n--- Processing Finished ---\n")