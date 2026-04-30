library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyverse)
setwd('./')
gene_matrix_sc<-read.csv('.csv')

species_order <- c("Hsa", "Mmu", "Tgu", "Tsc", "Xla", "Ame", "Pan", "Dre", "Cpl", "Lre", "Cin", "Bja", "Pye", "Afu", "Afa", "Dme")

cell_types <- colnames(gene_matrix_sc)[-1] %>% 
  str_extract("_[^_]+$") %>% 
  str_remove("_") %>% 
  unique()


new_col_order <- map_dfr(cell_types, function(ct) {
 
  ct_cols <- colnames(gene_matrix_sc)[str_detect(colnames(gene_matrix_sc), paste0("_", ct, "$"))]
  

  tibble(
    col_name = ct_cols,
    species = str_extract(ct_cols, "^[A-Za-z]{2,3}"),
    cell_type = ct
  ) %>%
    mutate(species_rank = match(species, species_order)) %>%
    arrange(species_rank) %>%
    select(col_name)
}) %>% pull(col_name)


gene_matrix_sc_sorted <- gene_matrix_sc %>%
  select(gene, all_of(new_col_order))

get_celltype <- function(gene) {
  str_extract(gene, "_[^_]+$") %>% str_remove("_")
}

#df<-gene_matrix_sc_sorted[,c(1,126:129,2:125)]
#split_column <- function(value) {
#  if (is.na(value)) {
#    return(c(NA, NA))
#  }
#  split_result <- strsplit(value, "_", fixed = TRUE)[[1]]
#  if (length(split_result) == 1) {
#    return(c(split_result[1], NA))
#  }
#  return(split_result)
#}

#df$symbol <- sapply(df$row_id, function(x) split_column(x)[1])
#df$celltype <- sapply(df$row_id, function(x) split_column(x)[2])

#filtered_df <- df %>%
#  group_by(celltype) %>%
#  mutate(row_sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
#  slice_max(order_by = row_sum, n = 20, with_ties = FALSE)
  
filtered_df<-gene_matrix_sc_sorted
filtered_df <- data.frame(lapply(filtered_df, function(x) ifelse(is.na(x), 0, x)))
filtered_df <- filtered_df[!duplicated(filtered_df$gene), ]
rownames(filtered_df)<-filtered_df$gene

mat <- filtered_df[,c(2:153)]
original_rownames <- rownames(mat)
mat <- as.data.frame(apply(mat, 2, as.numeric))
rownames(mat) <- original_rownames

is_all_zero <- function(x) {
  all(x == 0, na.rm = TRUE)
}
mat_clean <- mat[, !sapply(mat, is_all_zero)]
mat_clean[mat_clean > 5] <- 5
mat_clean <- mat_clean[rowSums(mat_clean, na.rm = TRUE) != 0, ]
mat_clean<-mat_clean[,c(6:49, 62:80, 1:5, 129:143, 119:128, 81:102, 111:118,50:61,103:110, 144:150)]
write.csv(mat_clean,'new-sc-celltypemarker.csv')
colors <- colorRampPalette(c('#DBDBDCFF','firebrick3'))(1000)
p <- pheatmap(mat_clean, color = colors,scale = 'none', cluster_cols = FALSE,cluster_row = FALSE,border_color = "T")
ggsave('new-sc-celltypemarker.pdf',p,wi=20,he=60,limitsize = FALSE)


########################################################################################################
gene_matrix_st<-read.csv('core_genes-ST-ALL.csv')

species_order <- c("Hsa", "Mmu", "Tgu", "Tsc", "Xla", "Ame", "Pan", "Dre", "Cpl", "Lre", "Cin", "Bja", "Pye", "Afu", "Afa", "Dme")

# 1. 
cell_types <- colnames(gene_matrix_st)[-1] %>% 
  str_extract("_[^_]+$") %>% 
  str_remove("_") %>% 
  unique()

# 2. 
new_col_order <- map_dfr(cell_types, function(ct) {
  ct_cols <- colnames(gene_matrix_st)[str_detect(colnames(gene_matrix_st), paste0("_", ct, "$"))]
  
  tibble(
    col_name = ct_cols,
    species = str_extract(ct_cols, "^[A-Za-z]{2,3}"),
    cell_type = ct
  ) %>%
    mutate(species_rank = match(species, species_order)) %>%
    arrange(species_rank) %>%
    select(col_name)
}) %>% pull(col_name)

gene_matrix_st_sorted <- gene_matrix_st %>%
  select(gene, all_of(new_col_order))

dt<-gene_matrix_st_sorted
dt <- data.frame(lapply(dt, function(x) ifelse(is.na(x), 0, x)))
dt <- dt[!duplicated(dt$gene), ]
rownames(dt)<-dt$gene
matt <- dt[,c(2:108)]

is_all_zero <- function(x) {
  all(x == 0, na.rm = TRUE)
}
matt_clean <- matt[, !sapply(matt, is_all_zero)]
matt_clean[matt_clean > 2] <- 2
matt_clean <- matt_clean[rowSums(matt_clean, na.rm = TRUE) != 0, ]
matt_clean<-matt_clean[,c(1:13, 26:49,56,14:25,50:55,76:81, 87:94, 60:70, 82:86, 95:103, 71:75, 104:106)]
matt_clean<-matt_clean[,c(1:101,103)]
matt_clean<-matt_clean[,c(1:70,87,71:86,88:102)]

write.csv(matt_clean,'new-st-celltypemarker.csv')
p1 <- pheatmap(matt_clean, color = colors, scale = 'none',cluster_cols = FALSE,cluster_row = FALSE,border_color = "T")
ggsave('new-st-celltypemarker.pdf',p1,wi=15,he=50,limitsize=F)

sc<-read.csv('new-sc-celltypemarker-select.csv',row.names=1)
colors <- colorRampPalette(c('#DBDBDCFF','firebrick3'))(1000)
p <- pheatmap(sc, color = colors,scale = 'none', cluster_cols = FALSE,cluster_row = FALSE,border_color = "T")
ggsave('sc-celltypemarker-select.pdf',p,wi=20,he=60,limitsize = FALSE)

st<-read.csv('new-st-celltypemarker-selected.csv',row.names=1)
p1 <- pheatmap(st, color = colors,scale = 'none', cluster_cols = FALSE,cluster_row = FALSE,border_color = "T")
ggsave('st-celltypemarker-select.pdf',p1,wi=20,he=60,limitsize = FALSE)
