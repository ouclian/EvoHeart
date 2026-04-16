rm(list = ls())
options(stringsAsFactors = FALSE)
gc()

# --- Load Libraries ---
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(patchwork)
library(reticulate)
library(ComplexHeatmap)
library(circlize)
library(paletteer)

# --- 1. Data Ingestion & Pre-processing ---
# Load core genes and merged data
core <- read.csv("core_genes.csv")
mergedata <- read.csv("./00data/py_merged_data.csv")[, -1]
mergedata[is.na(mergedata)] <- 0

# --- 2. Python Integration for Data Aggregation ---
# Use pandas to calculate the mean expression grouped by Gene_id
pd <- import("pandas")
py_mydata <- r_to_py(mergedata)
py_mydata_agg <- py_mydata$groupby("Gene_id")$mean()
data_aggregated <- py_to_r(py_mydata_agg)

# Clean column names and structure
mergedata <- data_aggregated
mergedata$Geneid <- rownames(data_aggregated)
# Move Geneid to the first column
mergedata <- mergedata[, c("Geneid", setdiff(colnames(mergedata), "Geneid"))]

# Standardize column prefixes (remove leading 'X')
colnames(mergedata) <- gsub("^X0", "0", colnames(mergedata))
colnames(mergedata) <- gsub("^X1", "1", colnames(mergedata))

# Sort columns alphabetically and filter by core genes
mergedata <- mergedata[, c("Geneid", sort(setdiff(colnames(mergedata), "Geneid")))]
mergedata <- mergedata[mergedata$Geneid %in% core[, "X"], ]
mydata <- column_to_rownames(mergedata, "Geneid")

# --- 3. Cell Type Filtering & Renaming ---
# Select specific cell types using regex
target_cells <- c("_Cardiomyocyte", "_Endothelial.cell", "_Fibroblast", "_Epicardium")
regex_pattern <- paste(target_cells, collapse = "|")
mydata <- mydata[, grep(regex_pattern, colnames(mydata))]

# Clean suffixes and standardize names
colnames(mydata) <- gsub("\\d+_", "", colnames(mydata))
colnames(mydata) <- gsub("_lfc", "", colnames(mydata))
colnames(mydata) <- gsub("Endothelial.cell", "Endothelial cell", colnames(mydata))

# --- 4. Species Name Mapping (Abbreviation to Full Name) ---
species_mapping <- c(
  "Dme" = "Fruit_fly", "Pye" = "Scallop", "Afa" = "Octopus", "Afu" = "Snail",
  "Bja" = "Amphioxus", "Cin" = "Ascidian", "Lre" = "Lamprey", "Cpl" = "Shark",
  "Dre" = "Zebrafish", "Pan" = "Lungfish", "Ame" = "Axolotl", "Xla" = "Xenopus",
  "Tsc" = "Turtle", "Tgu" = "Zebrafinch", "Mmu" = "Mouse", "Hsa" = "Human"
)

for (abbr in names(species_mapping)) {
  colnames(mydata) <- gsub(paste0("^", abbr, "_"), 
                           paste0(species_mapping[abbr], "_"), 
                           colnames(mydata))
}

# Re-order by cell type groups for better visualization
mydata <- mydata[, c(grep("_Cardiomyocyte", colnames(mydata)),
                     grep("_Endothelial", colnames(mydata)),
                     grep("_Fibroblast", colnames(mydata)),
                     grep("_Epicardium", colnames(mydata)))]

# --- 5. Scaling and Cleaning ---
index <- match(core$X, rownames(mydata))
mydata <- mydata[index, ]
mydata <- mydata[!is.na(mydata[, 1]), ] # Remove NA rows
mydata_cleaned <- mydata[apply(mydata, 1, function(x) sd(x) != 0), ]

# Z-score normalization (Scale by row)
data_matrix <- as.matrix(mydata_cleaned)
df_scaled <- t(scale(t(data_matrix)))

# --- 6. Heatmap Annotation Setup ---

# Row Annotations: Gene categories (Atrial vs Ventricular core)
row_split_idx <- 58 
left_group <- c(rep("Atrial core genes", row_split_idx), 
                rep("Ventricular core genes", nrow(df_scaled) - row_split_idx))
left_color <- c("Atrial core genes" = "#94ce91", "Ventricular core genes" = "#f9a025")

# Column Annotations: Cell type groups
top_group <- case_when(
  grepl("_Cardiomyocyte$", colnames(df_scaled)) ~ "Cardiomyocyte",
  grepl("_Endothelial cell$", colnames(df_scaled)) ~ "Endothelial cell",
  grepl("_Fibroblast$", colnames(df_scaled)) ~ "Fibroblast",
  TRUE ~ "Epicardium"
)
top_color <- c("Cardiomyocyte" = "#b2d18e", "Endothelial cell" = "#de823c", 
               "Fibroblast" = "#8bc3af", "Epicardium" = "#60418a")

# Define color gradient for expression values
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# --- 7. Building the Heatmap ---

# Identify specific genes to label on the right side
genes_to_label <- c("ALDH1A2", "C1QB", "FLRT2", "PITX2", "COL1A2", "EFEMP1", 
                    "OSR1", "IQGAP1", "ARPC3", "MYO18A", "LDB3", "PRDM16", 
                    "CACNA1C", "ANK2", "SLC8A1", "FHOD3", "CDH2", "ESRRG", 
                    "IRX2", "SORBS2")

ht <- Heatmap(df_scaled, 
              name = "Exp",
              col = col_fun,
              
              # Clustering and Labels
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_row_names = FALSE, 
              show_column_names = TRUE,
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 8),
              
              # Row Annotation (Left)
              left_annotation = rowAnnotation(
                "Gene type" = left_group,
                col = list("Gene type" = left_color),
                show_annotation_name = FALSE
              ),
              
              # Column Annotation (Top)
              top_annotation = HeatmapAnnotation(
                "Celltype" = top_group,
                col = list("Celltype" = top_color),
                show_annotation_name = FALSE
              ),
              
              # Legend Parameters
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                title_position = "topcenter"
              ),
              rect_gp = gpar(col = NA) # No cell borders
)

# Add right-side gene labels with lines (anno_mark)
row_anno_link <- rowAnnotation(
  link = anno_mark(
    at = which(rownames(df_scaled) %in% genes_to_label),
    labels = rownames(df_scaled)[rownames(df_scaled) %in% genes_to_label],
    labels_gp = gpar(fontsize = 10)
  )
)

# --- 8. Final Drawing & Export ---
final_plot <- ht + row_anno_link

# Export to PDF
pdf("./01pdf/Heatmap_Refined.pdf", height = 10, width = 13)
draw(final_plot, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right", 
     merge_legend = TRUE,
     padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()

message("Heatmap successfully exported to ./01pdf/Heatmap_Refined.pdf")