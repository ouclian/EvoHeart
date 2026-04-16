library(dplyr)
library(factoextra)
library(ggplot2)

################# step1. load data and PCA  
## Other cell types also use the same data preprocessing pipeline as CM, so only the PCA of CM is shown here.
cm = read.csv('.csv') 

rownames(cm) = cm$symbol
cm = cm[-c(1,2)]
data_t <- t(cm)
data_t[is.na(data_t)]=0
data_t = data_t[!grepl('Hsa',rownames(data_t)),]

zero_var_cols <- apply(data_t, 2, var) == 0
clean_data <- data_t[, !zero_var_cols]

expression_matrix <- scale(clean_data)
pca <- prcomp(expression_matrix,center=F,scale=F)
species_colors <- c(
  'Hsa' = '#75001a',
  'Mmu' = '#e1372f',
  'Tgu' = '#ee6a4d',
  'Tsc' = '#e16c00',
  'Xla' = '#eeb343',
  'Ame' = '#f6dc7b',
  'Pan' = '#a3a000',
  'Dre' = '#307331',
  'Cpl' = '#5fae5f',
  'Lre' = '#d1e7c8',
  'Cin' = '#9b8b85',
  'Bja' = '#b8afa7',
  'Pye' = '#a8cce3',
  'Afu' = '#4077b9',
  'Afa' = '#2c59a1',
  'Dme' = '#8482be'
)

################# step2. species group and plot
df = data.frame(pca$x)[c(1,2,3)]
df$species = gsub('\\..*','',rownames(df))
df$group <- ifelse(
  df$species %in% c("Cin", "Bja"),           
  "Protochordata",
  ifelse(                                    
    df$species %in% c( "Tsc", "Tgu", "Mmu"),
    "Amniota",
    "Anamniotes"                                 
  )
)

centroids <- df %>%
  group_by(species) %>%
  summarise(
    PC1_mean = mean(PC1),
    PC2_mean = mean(PC2),
    group = first(group)  
  )

p = ggplot(df, aes(x = PC1, y = PC2, color = species, fill = group)) +
  stat_ellipse(
    geom = "path",
    linetype = "dashed",
    color = "gray50",  
    aes(color = NULL)       
  ) +
  geom_point(size = 0.5, alpha = 0.8) + 
  geom_point(
    data = centroids,
    aes(x = PC1_mean, y = PC2_mean),
    shape = 18,        
    size = 3,
    color = 'black'
  ) +
  geom_jitter(width = 5, height = 3, alpha = 0.5, size = 0.5) +
  geom_text(
    data = centroids,
    aes(x = PC1_mean, y = PC2_mean, label = species),
    color = "black",
    vjust = -0.5,        
    fontface = "bold"
  ) +
  scale_color_manual(values = species_colors) +
  
  labs(title = "PCA with Confidence Ellipses by Group") +
  theme_minimal()
ggsave('CM_pca.pdf',p)

