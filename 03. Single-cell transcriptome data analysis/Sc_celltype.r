library(Seurat)
library(ggplot2)
setwd('F:/02.待发表数据/01.Evolution_BGI/00.AAA最终版/02.返修/01.sc-rds/')

sc_color_celltype = c(
'Adipocyte' = '#fff8dc',
'Atrial Cardiomyocyte' = '#e85e59',
'Ventricular Cardiomyocyte' = '#e85e59',
'Cardiomyocyte' = '#e85e59',

'Lymphatic Endothelial cell' = '#91d0be',
'Endothelial cell' = '#91d0be',
'Epicardium' = '#062ad7',
'Epithelial cell' = '#8c564b',
'Fibroblast' =  '#00BBFFFF',
'Immune cell' = '#6778ae',
'Lymphoid'= '#6778ae',
'Myeloid' = '#6778ae',
'Neural cell' = '#2DDBA7',
'Proliferating cell' = '#008000',
'Red blood cell' = '#FE1D05FF',
'Smooth muscle cell' = '#f4be1d',
'Valve' = '#ffff9a')

rds_files <- list.files(pattern = "\\.rds$")

for (rds_file in rds_files) {
  data <- readRDS(rds_file)
  celltype_table <- prop.table(table(data$celltype))
  filename<-basename(rds_file)
  write.csv(celltype_table, paste0(filename,'_celltype_cellratio.csv'))

  p1 <- DimPlot(data,group.by='celltype',cols=sc_color_celltype)
  ggsave(paste0(filename,'_sc.pdf'), p1, width = 10, height = 10, dpi = 300,limitsize = FALSE)  
}
