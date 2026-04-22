library(Seurat)
library(ggplot2)
library(gghighlight)
setwd('F:/02.待发表数据/01.Evolution_BGI/00.AAA最终版/02.返修/03.st-bin50-rds/')

st_color_celltype = c(
    'Epicardium' = '#062ad7',
    'Atrial Epicardium' = '#062ad7',
    'Ventricular Epicardium' = '#062ad7',
    'Epicardium & Fibroblasts' = '#062ad7',

    'Ventricular Cardiomyocytes' = '#F59E9BFF',
    'Ventricular Cardiomyocyte' = '#F59E9BFF',
    'Right_Ventricular_Cardiomyocyte' = '#F59E9BFF',

    'Atrial Cardiomyocytes' = '#3181be',
    'Atrial Cardiomyocyte' = '#3181be',
    'Atrial Cardiomyocytes & Fibroblasts' = '#3181be',

    'Cardiomyocytes' = '#e85e59',
    'Cardiomyocyte' = '#e85e59',
    'Cardiomycyte' = '#e85e59',

    'Red blood cells' = '#FE1D05FF',
    'Red blood cell' = '#FE1D05FF',

    'Endocardium' = '#91d0be',
    'Endothelial cells' = '#91d0be',
    'Endothelial cell ' = '#91d0be',
    'Endothelial cell' = '#91d0be',
    'Lymphatic Endothelial cell' = '#91d0be',

    'Valve' = '#ffff9a',

    'Fibroblasts' =  '#00BBFFFF',
    'Fibroblast' = '#00BBFFFF',
    'Fibroblasts & Cardiomyocytes' =  '#00BBFFFF',

    'Smooth muscle cells' = '#f4be1d',
    'Smooth muscle cell' = '#f4be1d',

    'Immune cells'= '#6778ae',
    'Immune cell' = '#6778ae',
    'Myeloid' = '#6778ae',
    'Myeloid ' = '#6778ae',
    'Lymphoid' = '#6778ae',

    'Proliferating cells' = '#008000',
    'Proliferating cell' = '#008000',

    'Epithelial cells' = '#8c564b',
    'Epithelial cell' = '#8c564b',
    'Ciliated Epithelial cell' = '#8c564b',

    'Endothelial_Neural cell'= '#2DDBA7',
    'Endothelial_neuron cell' = '#2DDBA7',
    'Neural cell' = '#2DDBA7',
    'Neurons' = '#2DDBA7',

    'Adipocyte' = '#fff8dc',
    'Platelet' = '#EBCBF8FF',
    'Perivascular cells' = '#D0D78EFF',
    'Mast cell' = '#E19A80FF',
    'NA' = '#B9B9B8FF',
    'others' = '#B9B9B8FF'
)

rds_files <- list.files(pattern = "\\.rds$")

for (rds_file in rds_files) {
  data <- readRDS(rds_file)
  celltype_table <- prop.table(table(data$celltype))
  filename<-basename(rds_file)
  write.csv(celltype_table, paste0(filename,'_celltype_cellratio.csv'))

  sample_meta.data<-data@meta.data
  p1 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 1) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0(filename,'_spatial_1.pdf'), p1, width = 20, height = 20, dpi = 300,limitsize = FALSE)  
  
  p2 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 1.5) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0(filename,'_spatial_1.5.pdf'), p2, width = 20, height = 20, dpi = 300,limitsize = FALSE)  

  p3 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 2) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0(filename,'_spatial_2.pdf'), p3, width = 20, height = 20, dpi = 300,limitsize = FALSE)  
}



rds_files <- list.files(pattern = "\\.rds$")

for (rds_file in rds_files) {
  data <- readRDS(rds_file)
  filename<-basename(rds_file)
  sample_meta.data <- data@meta.data
  celltype <- unique(data$celltype)

  for (i in celltype) {
    sizes <- c(1, 1.5, 2)
    for (s in sizes) {
      p <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) +
        geom_point(shape = 19, size = s) +
        theme_void() +
        coord_fixed(ratio = 1) +
        theme(legend.position = 'right') +
        scale_color_manual(values = st_color_celltype) +
        gghighlight(celltype == i)

      ggsave(paste0(filename, '_', i, '_', s, '.pdf'),
             p, width = 20, height = 20, dpi = 300, limitsize = FALSE)
    }
  }
}



st_color_celltype = c(
    'Ventricular Cardiomyocytes' = '#F59E9BFF',
    'Ventricular Cardiomyocyte' = '#F59E9BFF',
    'Right_Ventricular_Cardiomyocyte' = '#F59E9BFF',

    'Atrial Cardiomyocytes' = '#3181be',
    'Atrial Cardiomyocyte' = '#3181be',
    'Atrial Cardiomyocytes & Fibroblasts' = '#3181be',

    'Cardiomyocytes' = '#e85e59',
    'Cardiomyocyte' = '#e85e59',
    'Cardiomycyte' = '#e85e59'
)

rds_files <- list.files(pattern = "\\.rds$")

for (rds_file in rds_files) {
  data <- readRDS(rds_file)
  celltype_table <- prop.table(table(data$celltype))
  filename<-basename(rds_file)
  write.csv(celltype_table, paste0('./01.Cardio/',filename,'_celltype_cellratio.csv'))

  sample_meta.data<-data@meta.data
  p1 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 1) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0('./01.Cardio/',filename,'_spatial_1.pdf'), p1, width = 20, height = 20, dpi = 300,limitsize = FALSE)  
  
  p2 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 1.5) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0('./01.Cardio/',filename,'_spatial_1.5.pdf'), p2, width = 20, height = 20, dpi = 300,limitsize = FALSE)  

  p3 <- ggplot(sample_meta.data, aes(x = x, y = y, color = celltype)) + geom_point(shape = 19, size = 2) + theme_void() + coord_fixed(ratio = 1) + theme(legend.position = 'right') + scale_color_manual(values = st_color_celltype)

  ggsave(paste0('./01.Cardio/',filename,'_spatial_2.pdf'), p3, width = 20, height = 20, dpi = 300,limitsize = FALSE)  
}




