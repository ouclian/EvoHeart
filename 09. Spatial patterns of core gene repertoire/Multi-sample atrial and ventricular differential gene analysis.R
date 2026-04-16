rm(list = ls())
options(stringsAsFactors = FALSE)

# Set working directory
wd <- "."
setwd(wd)

# --- Load Libraries ---
library(readr)
library(DESeq2)
library(ggfortify)
library(cluster)
library(reticulate)
library(dplyr)

# --- Data Loading and Pre-processing ---

# Load gene ID dictionary
id <- read.table("./geneid_dict.csv", header = TRUE, sep = ",")
id$Tsc <- gsub("gene-", "", id$Tsc)

# Load main count matrix (large CSV)
tem <- read_csv("./00data/A_Vmerge.csv")
data_raw <- tem[, -1] # Remove the first index column

# --- Python Integration for Data Aggregation ---
# Use Python's pandas to sum rows with the same Geneid efficiently
pd <- import("pandas")
py_data <- r_to_py(data_raw)
py_data_agg <- py_data$groupby("Geneid")$sum()
data_aggregated <- py_to_r(py_data_agg)

# Check aggregated data structure
head(data_aggregated)
colnames(data_aggregated)

# ==============================================================================
# Function: diffAV
# Purpose: Perform DESeq2 analysis, save results
# Arguments: alist (atrium samples), vlist (ventricle samples)
# ==============================================================================
diffAV <- function(alist, vlist) {
  
  # 1. Prepare count matrix and metadata
  # Subset data for specific species samples and round to integers for DESeq2
  test_counts <- round(data_aggregated[, c(alist, vlist)])
  
  # Define experimental conditions
  condition <- factor(c(rep("A", length(alist)), rep("V", length(vlist))))
  coldata <- data.frame(row.names = colnames(test_counts), condition)
  
  # 2. Initialize DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = test_counts, 
                                 colData = coldata, 
                                 design = ~ condition)
  
  # Set 'Ventricle' as the reference level
  dds$condition <- relevel(dds$condition, ref = "V")
  
  # 3. Data Filtering and Transformation
  # Keep genes with at least 10 counts in at least one sample
  keep <- rowSums(counts(dds) >= 10) >= 1
  dds <- dds[keep, ] 
  
  # Variance Stabilizing Transformation (for PCA/visualization)
  vst_data <- vst(dds, blind = TRUE)  
  
  # 4. Run Differential Expression Analysis
  dds <- DESeq(dds) 
  res <- results(dds)
  
  # 5. Result Merging and Filtering
  # Merge stats with VST normalized counts
  res_df <- as.data.frame(res)
  vst_matrix <- as.data.frame(assay(vst_data))
  resdata_merged <- merge(res_df, vst_matrix, by = "row.names", sort = FALSE)
  
  # Filter by significance (padj < 0.05) and remove NAs
  resdata_sig <- resdata_merged %>%
    filter(padj < 0.05) %>%
    filter(!is.na(Row.names)) %>%
    arrange(desc(log2FoldChange))
  
  # 6. Define Up-regulated (Atrium) and Down-regulated (Ventricle) genes
  # Fold change threshold: |log2FC| > 1.5
  adata <- resdata_sig %>% filter(log2FoldChange > （1.5）)
  vdata <- resdata_sig %>% filter(log2FoldChange < （-1.5）)
  
  # 7. Output Result Files
  dir.create("./00data/01atrium", showWarnings = FALSE, recursive = TRUE)
  dir.create("./00data/02ventricle", showWarnings = FALSE, recursive = TRUE)
  
  write.csv(adata, paste0("./00data/01atrium/", alist[1], ".csv"), row.names = FALSE)
  write.csv(vdata, paste0("./00data/02ventricle/", vlist[1], ".csv"), row.names = FALSE)
}

# ==============================================================================
# Execution Section: Process Multiple Species
# ==============================================================================

# Species: 01_Afu
diffAV(alist = c("01_Afu_A", "01_Afu_A.1", "01_Afu_A.2", "01_Afu_A.3", "01_Afu_A.4"),
       vlist = c("01_Afu_V", "01_Afu_V.1", "01_Afu_V.2", "01_Afu_V.3", "01_Afu_V.4", "01_Afu_V.5"))

# Species: 03_Tcr
diffAV(alist = c("03_Tcr_A", "03_Tcr_A.1", "03_Tcr_A.2"),
       vlist = c("03_Tcr_V", "03_Tcr_V.1", "03_Tcr_V.2"))

# Species: 04_Cgi
diffAV(alist = c("04_Cgi_A", "04_Cgi_A.1", "04_Cgi_A.2"),
       vlist = c("04_Cgi_V", "04_Cgi_V.1", "04_Cgi_V.2"))

# Species: 05_Pye
diffAV(alist = c("05_Pye_A", "05_Pye_A.1", "05_Pye_A.2", "05_Pye_A.3", "05_Pye_A.4", "05_Pye_A.5"),
       vlist = c("05_Pye_V", "05_Pye_V.1", "05_Pye_V.2", "05_Pye_V.3"))

# Species: 10_Lre
diffAV(alist = c("10_Lre_A", "10_Lre_A.1", "10_Lre_A.2"),
       vlist = c("10_Lre_V", "10_Lre_V.1", "10_Lre_V.2"))

# Species: 11_Dre
diffAV(alist = c("11_Dre_A", "11_Dre_A.1", "11_Dre_A.2"),
       vlist = c("11_Dre_V", "11_Dre_V.1", "11_Dre_V.2"))

# Species: 14_Ame
diffAV(alist = c("14_Ame_A", "14_Ame_A.1", "14_Ame_A.2"),
       vlist = c("14_Ame_V", "14_Ame_V.1", "14_Ame_V.2", "14_Ame_V.3", "14_Ame_V.4", "14_Ame_V.5"))

# Species: 16_Xla
diffAV(alist = c("16_Xla_A", "16_Xla_A.1", "16_Xla_A.2", "16_Xla_A.3", "16_Xla_A.4", "16_Xla_A.5", "16_Xla_A.6"),
       vlist = c("16_Xla_V", "16_Xla_V.1", "16_Xla_V.2", "16_Xla_V.3", "16_Xla_V.4"))

# Species: 19_Gga
diffAV(alist = c("19_Gga_A", "19_Gga_A.1", "19_Gga_A.2"),
       vlist = c("19_Gga_V", "19_Gga_V.1", "19_Gga_V.2"))

# Species: 20_Tgu
diffAV(alist = c("20_Tgu_A", "20_Tgu_A.1", "20_Tgu_A.2"),
       vlist = c("20_Tgu_V", "20_Tgu_V.1", "20_Tgu_V.2"))

# Species: 21_Rno
diffAV(alist = c("21_Rno_A", "21_Rno_A.1", "21_Rno_A.2"),
       vlist = c("21_Rno_V", "21_Rno_V.1", "21_Rno_V.2"))

# Species: 22_Mmu
diffAV(alist = c("22_Mmu_A", "22_Mmu_A.1"),
       vlist = c("22_Mmu_V", "22_Mmu_V.1", "22_Mmu_V.2"))

# Species: 23_Tbe
diffAV(alist = c("23_Tbe_A", "23_Tbe_A.1", "23_Tbe_A.2"),
       vlist = c("23_Tbe_V", "23_Tbe_V.1", "23_Tbe_V.2"))