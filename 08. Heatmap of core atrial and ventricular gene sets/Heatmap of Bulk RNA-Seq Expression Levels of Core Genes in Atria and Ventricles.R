rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Set working directory
wd <- "."
setwd(wd)

# --- Load Libraries ---
library(pheatmap)
library(reticulate)
library(ComplexHeatmap)
library(circlize)

# ==============================================================================
# Part 1: Data Aggregation via Python (reticulate)
# ==============================================================================
raw_data_path <- "./bulk_a_v_merge.csv"
data <- read.csv(raw_data_path)

# Standardize column names (remove 'X' prefix from numeric strings)
colnames(data) <- gsub("X", "", colnames(data))
col_order <- order(colnames(data))
data <- data[, c(1, col_order[1:(length(col_order) - 1)])]
data[is.na(data)] <- 0

# --- Python Processing via Pandas ---
pd <- import("pandas")
py_df <- r_to_py(data)
# Sum up log2FC values grouped by Row.names
py_agg <- py_df$groupby("Row.names")$sum()
data_agg <- py_to_r(py_agg)

mergedata <- data_agg
mergedata$Row.names <- rownames(data_agg)
# Reorganize columns to put Row.names first
mergedata <- mergedata[, c("Row.names", setdiff(colnames(mergedata), "Row.names"))]
colnames(mergedata) <- gsub("X", "", colnames(mergedata))

# Set rownames for downstream filtering
rownames(mergedata) <- mergedata$Row.names
data_to_flip <- mergedata

# ==============================================================================
# Part 2: Adjust Sign Polarity 
# ==============================================================================
# Ensure Atria remains positive and Ventricle becomes negative
for (i in seq(2, 32, 2)) {
  data_to_flip[, i] <- ifelse(mergedata[, i] != 0, mergedata[, i], 0)
  data_to_flip[, i + 1] <- ifelse(mergedata[, i + 1] != 0, -mergedata[, i + 1], 0)
}

data_processed <- data_to_flip

# ==============================================================================
# Part 3: Filtering Core Genes
# ==============================================================================
# Load core gene list
core_gene_path <- "./core_genes.csv"
core_list <- read.csv(core_gene_path)

# Filter data to keep only core genes
data_filtered <- data_processed[data_processed$Row.names %in% core_list$X, ]

# Match row order to the core list order
match_idx <- match(trimws(core_list$X), trimws(data_filtered$Row.names))
data_filtered <- data_filtered[match_idx, ]
data_filtered <- data_filtered[!duplicated(data_filtered$Row.names) & !is.na(data_filtered$Row.names), ]

rownames(data_filtered) <- data_filtered$Row.names
plot_data <- data_filtered[, -1] 

# Select specific species for single-cell cross-referencing
selected_species <- c("Afu_A", "Pye_A", "Lre_A", "Dre_A", "Cpl_A", "Ame_A", "Xla_A", "Tgu_A", "Mmu_A",
                      "Afu_V", "Pye_V", "Lre_V", "Dre_V", "Cpl_V", "Ame_V", "Xla_V", "Tgu_V", "Mmu_V")

plot_data <- plot_data[, selected_species]
plot_data <- plot_data[apply(plot_data, 1, function(x) sd(x) != 0), ]
df_scaled <- t(scale(t(plot_data)))

# ==============================================================================
# Part 4: Heatmap Setup (Colors & Annotations)
# ==============================================================================

# Define Color mapping
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# --- Row Annotations (Left Side) ---
# Group genes into Atrial vs Ventricular categories
left_group <- c(rep("Atrial core genes", 75), rep("Ventricular core genes", 121))
left_color_map <- c("Atrial core genes" = "#94ce91", "Ventricular core genes" = "#f9a025")

# --- Top Annotations (Column Side) ---
# Identify Chamber type based on column suffix
top_group <- ifelse(grepl("_A$", colnames(df_scaled)), "Atrial", "Ventricular")
top_color_map <- c("Atrial" = "#66c2a5", "Ventricular" = "#8da0cb")

# Define genes to be labeled specifically
genes_to_label <- c("ALDH1A2", "C1QB", "FLRT2", "PITX2", "COL1A2", "EFEMP1", 
                    "OSR1", "IQGAP1", "ARPC3", "MYO18A", "LDB3", "PRDM16", 
                    "CACNA1C", "ANK2", "SLC8A1", "FHOD3", "CDH2", "ESRRG", 
                    "IRX2", "SORBS2")

# ==============================================================================
# Part 5: Final Heatmap Generation & Export
# ==============================================================================

# Create main heatmap object
main_ht <- Heatmap(df_scaled,
              name = "Exp",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = TRUE,
              column_names_rot = 45,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 8),
              
              # Left annotation (Gene Category)
              left_annotation = rowAnnotation(
                "Gene type" = left_group,
                col = list("Gene type" = left_color_map),
                show_annotation_name = FALSE
              ),
              
              # Top annotation (Chamber)
              top_annotation = HeatmapAnnotation(
                "Chamber" = top_group,
                col = list("Chamber" = top_color_map),
                show_annotation_name = FALSE
              ),
              rect_gp = gpar(col = "white", lwd = 0)
)

# Define gene labels for the right side
marker_anno <- rowAnnotation(
  link = anno_mark(
    at = which(rownames(df_scaled) %in% genes_to_label),
    labels = rownames(df_scaled)[rownames(df_scaled) %in% genes_to_label],
    labels_gp = gpar(fontsize = 10)
  )
)

# Export to PDF
pdf("./bulk_Heatmap.pdf", height = 10, width = 13)
draw(main_ht + marker_anno, merge_legend = TRUE)
dev.off()

message("Success: Heatmap exported to ./bulk_Heatmap.pdf")