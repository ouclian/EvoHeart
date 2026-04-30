library(Seurat)
library(ggplot2)
library(dplyr)
library(ggforce)
library(tidyr)

setwd('./')
data <- readRDS("./.rds")

st_color_celltype = c(
'CM1(const.)' = '#B2DF8AFF',
'CM3(a const.)' = '#FB9A99FF',
'CM5(v const.)' = '#57C3F3',
'CM2(trans.)' = '#016D06',
'CM4(a trans.)' = '#8B0000',
'CM6(v trans.)' = '#12126E')
p <- DimPlot(data, group.by = 'labels3', cols = st_color_celltype, reduction='tsne',raster=FALSE)
ggsave('saturn.pdf', p, width = 20, height = 16, dpi = 600, bg = "transparent")


mycolors <- c('#A6CEE3FF','#1F78B4FF', '#B2DF8AFF', '#33A02CFF' ,'#FB9A99FF' ,'#E31A1CFF' ,'#FDBF6FFF' ,'#FF7F00FF','#CAB2D6FF' ,'#6A3D9AFF' ,'#FFFF99FF' ,'#B15928FF','#DA70D6FF','#A33C79',  '#57C3F3', '#4E4646','#110D8E', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175','#FFA500FF', '#800080FF', '#008000FF', '#FFFF00FF', '#00FFFFFF',  '#FF00FFFF', '#8B4513FF', '#A9A9A9FF', '#DEB887FF', '#5F9EA0FF',  '#7FFF00FF', '#D2691EFF', '#6495EDFF', '#FFF8DCFF', '#DC143CFF',  '#00FFFFFF', '#00008BFF', '#008B8BFF', '#B8860BFF', '#A0522DFF',  '#4B0082FF', '#556B2FFF', '#FF8C00FF', '#9932CCFF', '#8B008BFF',  '#E9967AFF', '#8FBC8FFF', '#4169E1FF', '#FFA07AFF', '#DAA520FF',  '#808000FF', '#BDB76BFF', '#800000FF', '#483D8BFF', '#2F4F4FFF',  '#228B22FF', '#FAF0E6FF', '#FFD700FF', '#ADD8E6FF', '#F08080FF',  '#E6E6FAFF', '#708090FF', '#FFF0F5FF', '#F0E68CFF', '#EEEEEEFF',  '#7B68EEFF', '#F5FFFAFF', '#FFE4E1FF', '#00FF00FF', '#98FB98FF',  '#FFEBCDFF', '#D3D3D3FF', '#F8F9FAFF', '#FFB6C1FF', '#DA70D6FF')

data$spe_celltype <- paste0(data$species, '_', data$labels3)
p1 <- DimPlot(data, group.by = 'species', cols =mycolors , label=TRUE, reduction='tsne',raster=FALSE)
ggsave('saturn_species.pdf', p1, width = 20, height = 16, dpi = 600, bg = "transparent")

subset_data <- subset(
  x = data,
  subset = labels3 %in% c("CM1(const.)", "CM2(trans.)", "CM3(a const.)", "CM4(a trans.)", "CM5(v const.)", "CM6(v trans.)")
)

avg_exp <- AverageExpression(
  subset_data,
  assays = "RNA",
  #features = variable_genes,
  group.by = "spe_celltype",
  slot = "data"  
)
exp_matrix <- avg_exp$RNA
exp_matrix <- as.matrix(avg_exp$RNA)
cor_matrix <- cor(exp_matrix,use = "complete.obs", method = "pearson")

write.csv(cor_matrix,
          file = "cor_matrix.csv",
          row.names = TRUE)
data<- cor_matrix         
data <- read.csv('cor_matrix.csv',row.names=1)
name <- read.csv('rename.txt',sep='\t')

name_map <- setNames(name$allname, name$name)
replace_abbreviation <- function(x, map) {
  # 提取前3个字母作为可能的缩写
  abbr <- substr(x, 1, 3)
  # 如果缩写在映射表中，则替换前3个字母，否则保持不变
  if (abbr %in% names(map)) {
    paste0(map[abbr], substr(x, 4, nchar(x)))  # 拼接全称+剩余文字
  } else {
    x  # 不匹配时保持原样
  }
}

rownames(data) <- sapply(rownames(data), replace_abbreviation, map = name_map)
colnames(data) <- sapply(colnames(data), replace_abbreviation, map = name_map)

write.csv(data,
          file = "cor_matrix_rename.csv",
          row.names = TRUE)
library(pheatmap)
p1 <- pheatmap(data,
         color = colorRampPalette(c('#0a34f1', 'white', '#ec3527'))(100),
         breaks = seq(-1, 1, length.out = 100),  
         clustering_method = "complete",
         clustering_distance_rows = as.dist(1-data),
         clustering_distance_cols = as.dist(1-data),
         display_numbers = FALSE,  
         border_color = NA,  
         treeheight_row = 50,
         treeheight_col = 50,
         angle_col = 90,
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Correlation between Cell Types across Species")

ggsave("celltype_correlation_heatmap.pdf", p1,width = 18, height = 18,dpi=300)
