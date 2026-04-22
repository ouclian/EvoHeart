library(Seurat)
library(ggplot2)
library(dplyr)
library(ggforce)
library(tidyverse)
library(pheatmap)
library(data.table)

#####添加人的marker
setwd("./")
data<-read.csv('Human_top100marker.csv')
df<- data[,c('genename.SYMBOL','celltype','avg_log2FC')]
data_filtered <- df[!is.na(df$genename.SYMBOL), ]

dt <- as.data.table(data_filtered)
dt_wide <- dcast(
  dt,
  genename.SYMBOL ~ celltype,
  value.var = "avg_log2FC",
  fun.aggregate = max, 
  fill = 0
)
rownames(dt_wide)<-dt_wide$genename.SYMBOL
df <- dt_wide[,c(4,3)]
colnames(df)<-c('Hsa_Ventricular.Cardiomyocytes..CACNA1C.','Hsa_Ventricular.Cardiomyocytes')
all<-read.csv('01.ALL_marker_spciecs_sort.csv',row.names = 1)
all$genename<-rownames(all)
df$genename<-rownames(df)
merged_df <- merge(all, df, by = "genename", all = TRUE)

all_cols <- names(merged_df)
new_order <- c(all_cols[1:32], all_cols[52], all_cols[33:51])
merged_df <- merged_df[, new_order]

###重新画图
all_data <- merged_df
rownames(all_data) <- all_data$genename
all_data <- all_data[,-1]
all_data[is.na(all_data)] <- 0

non_zero_count_all <- apply(all_data != 0, 1, sum)
non_zero_percentage <- non_zero_count_all / ncol(all_data)
all_cm <- all_data[non_zero_percentage >= 0.2, ]

non_zero_A<- apply(all_data[, 11:19] != 0, 1, all)
non_zero_ratio_A <- rowSums(all_data[, 11:19] != 0) / 9
filtered_rows_A <- all_data[non_zero_ratio_A >= 0.3, ]

non_zero_V<- apply(all_data[, 20:32] != 0, 1, all)
non_zero_ratio_V <- rowSums(all_data[, 20:32] != 0) / 13
filtered_rows_V <- all_data[non_zero_ratio_V >= 0.3, ]

non_zero_CA<- apply(all_data[, 33:51] != 0, 1, all)
non_zero_ratio_CA <- rowSums(all_data[, 33:51] != 0) / 19
filtered_rows_CA <- all_data[non_zero_ratio_CA >= 0.3, ]

non_zero_CA_A<- apply(all_data[, 36:42] != 0, 1, all)
non_zero_ratio_CA_A <- rowSums(all_data[, 36:42] != 0) / 7
filtered_rows_CA_A <- all_data[non_zero_ratio_CA_A >= 0.3, ]

non_zero_CA_V<- apply(all_data[, 43:51] != 0, 1, all)
non_zero_ratio_CA_V <- rowSums(all_data[, 43:51] != 0) / 9
filtered_rows_CA_V <- all_data[non_zero_ratio_CA_V >= 0.3, ]

all.list<-rownames(all_cm)
a.list<-rownames(filtered_rows_A)
v.list<-rownames(filtered_rows_V)
ca.list<-rownames(filtered_rows_CA)
caa.list<-rownames(filtered_rows_CA_A)
cav.list<-rownames(filtered_rows_CA_V)

# 使用循环删除gene_list1中出现在gene_list2中的基因
for (gene in all.list) {
  a.list <- a.list[a.list != gene]
}

for (gene in all.list) {
  v.list <- v.list[v.list != gene]
}

for (gene in ca.list) {
  caa.list <- caa.list[caa.list != gene]
}

for (gene in ca.list) {
  cav.list <- cav.list[cav.list != gene]
}

for (gene in ca.list) {
  all.list <- all.list[!all.list == gene]
}

for (gene in a.list) {
  if (gene %in% v.list) {
    all.list <- c(all.list, gene)
    a.list <- a.list[!a.list == gene]
    v.list <- v.list[!v.list == gene]
  }
}

gene_lists <- list(all.list,a.list,v.list,ca.list,caa.list,cav.list)
all_genes <- unlist(gene_lists)

matching_rows <- rownames(all_data) %in% all_genes
all_gene_data <- all_data[matching_rows, ]

all_gene_data$annotation <- NA

matching_rows2 <- rownames(all_gene_data) %in% all.list
all_gene_data$annotation[matching_rows2] <- "ALL"

matching_rows3 <- rownames(all_gene_data) %in% a.list
all_gene_data$annotation[matching_rows3] <- "CM(A)"
matching_rows4 <- rownames(all_gene_data) %in% v.list
all_gene_data$annotation[matching_rows4] <- "CM(V)"
matching_rows5 <- rownames(all_gene_data) %in% ca.list
all_gene_data$annotation[matching_rows5] <- "CM(CACNA1C)"
matching_rows6 <- rownames(all_gene_data) %in% caa.list
all_gene_data$annotation[matching_rows6] <- "CM(CACNA1C)_A"
matching_rows7 <- rownames(all_gene_data) %in% cav.list
all_gene_data$annotation[matching_rows7] <- "CM(CACNA1C)_V"

all_gene_data$annotation <- factor(all_gene_data$annotation, levels = c("ALL", "CM(A)",'CM(V)','CM(CACNA1C)','CM(CACNA1C)_A','CM(CACNA1C)_V'))
sorted_index <- order(all_gene_data$annotation)
sorted_data <- all_gene_data[sorted_index, ]

annotation_row <- as.data.frame(all_gene_data$annotation)
rownames(annotation_row) <- rownames(all_gene_data)
colnames(annotation_row)<-c('Geneset')
annotation_col <- data.frame(
  celltype = rep(c("Cardiomyocytes", "Atrial Cardiomyocytes",'Ventricular Cardiomyocytes','Cardiomyocytes (CACNA1C)','Atrial Cardiomyocytes (CACNA1C)','Ventricular Cardiomyocytes (CACNA1C)'), times = c(10,9,13,3,7,9)))
rownames(annotation_col) <- colnames(all_gene_data[,1:51])

col_labels <- substr(colnames(sorted_data[,1:51]), 1, 3)
anno_colors <- list(
  celltype = c('Atrial Cardiomyocytes'='#FB9A99FF',
             'Ventricular Cardiomyocytes'='#57C3F3',
             'Cardiomyocytes'='#B2DF8AFF',
             'Cardiomyocytes (CACNA1C)'='#016D06',
             'Atrial Cardiomyocytes (CACNA1C)'='#8B0000',
             'Ventricular Cardiomyocytes (CACNA1C)'='#12126E'),
  Geneset = c("ALL" = "#BD6263", "CM(A)" = "#A9D179","CM(V)"="#84CAC0", "CM(CACNA1C)" = "#F5AE6B","CM(CACNA1C)_A"="#BCB8D3","CM(CACNA1C)_V"="#4387B5"))

colors <- colorRampPalette(c("#E6E6E6FF",'red','darkred'))(1000)


p<-pheatmap(sorted_data[,1:51], annotation_row = annotation_row,annotation_col=annotation_col,scale='none',cluster_cols=F,cluster_rows=F,color = colors, labels_col = col_labels,gaps_col = c(10,19,32,35,42),gaps_row = c(5,16,66,108,116),annotation_colors = anno_colors)
ggsave('cm_subtype.pdf',p,wi=14,he=18)

cagene<- sorted_data[67:136,]
write.csv(rownames(cagene),'cagene.csv')
