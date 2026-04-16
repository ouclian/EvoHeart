library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(pheatmap)

################### Step1: Read data and build orthogroup-gene correspondence ###################
og = read.table('./orthogroup_proid2geneid.txt',sep = '\t',quote = "")
symbols = read.table('./06.filtered_genesymbol_toupper2.txt',quote = "",sep = '\t')
match_gene_to_og <- function(gene_list, og, species) {
  og <- as.data.table(og)
  expanded <- og[, .(gene = trimws(unlist(strsplit(get(species), ",")))), by = Orthogroup]
  expanded <- expanded[gene != "" & !is.na(gene)]
  mapping <- setNames(expanded$Orthogroup, expanded$gene)
  return(mapping[gene_list])
}
files = list.files('findmarkers/Cardiomyocyte/',full.names = T)

all = data.frame()
for(i in files[-10]){
  data = read.csv(i)
  species = paste0(substr(strsplit(i,'/|_')[[1]][3],1,1),substr(strsplit(i,'/|_')[[1]][4],1,2))
  new_table <- data %>%
    mutate(
      expr = case_when(
        pct.1 > 0.25 & avg_log2FC >= 1 ~ 1,
        pct.1 > 0.1 & avg_log2FC < 1 ~ 0.5,
        TRUE ~ 0  
        # Assign expression level grades (expr value) to each gene, divided into three grades:
        # 1. expr = 1: Highly expressed gene (pct.1 >= 0.25 and avg_log2FC >= 1)
        # 2. expr = 0.5: Moderately expressed gene (pct.1 >= 0.1 and avg_log2FC < 1)
        # 3. expr = 0: Low/non-expressed gene (all other cases)
      ),
      species = species
    ) 
  new_table$orthogroup = match_gene_to_og(new_table$gene, og, species)  
  new_table$symbols = symbols$genesymbol[match(new_table$gene,symbols[[species]])]
  new_table = new_table[c('gene','symbols','species','orthogroup','expr')]
all = rbind(all,new_table)  
}
write.csv(all,'coreMarkers/CM_all_pct.csv')

########### Step2. Find orthogroups where symbols have expr=1 in more than half of the species

all_species_list <- unique(all$species)
protostome_species <- c("Dme", "Afa", "Lfu", "Pye")  
deuterostome_species <- c("Bja","Cro",  "Lre", "Cpl", 
                          "Dre", "Pan", "Ame", "Xla",
                          "Tsc", "Tgu", "Mmu", "Hsa")  

# Alternative: create a species classification dataframe
species_class <- data.frame(
  species = all_species_list,
  group = c(rep("protostome", 4), rep("deuterostome", 12))  # Adjust if your order is different
)

# Merge with original data
all_with_group <- all %>%
  left_join(species_class, by = "species")

# Calculate statistics for each symbol and orthogroup
symbol_og_stats <- all_with_group %>%
  group_by(symbols, orthogroup) %>%
  summarise(
    # For protostomes: expr=1 in >=2 species
    protostome_expr1 = sum(group == "protostome" & expr == 1, na.rm = TRUE),
    protostome_total = sum(group == "protostome", na.rm = TRUE),
    
    # For deuterostomes: expr=1 in >=6 species
    deuterostome_expr1 = sum(group == "deuterostome" & expr == 1, na.rm = TRUE),
    deuterostome_total = sum(group == "deuterostome", na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  # Filter: (protostome_expr1 >= 2) OR (deuterostome_expr1 >= 6)
  filter(protostome_expr1 >= 2 | deuterostome_expr1 >= 6) %>%
  # Calculate additional statistics
  mutate(
    protostome_ratio = ifelse(protostome_total > 0, protostome_expr1/protostome_total, 0),
    deuterostome_ratio = ifelse(deuterostome_total > 0, deuterostome_expr1/deuterostome_total, 0)
  )

# Filter the original data to only include these symbol-og pairs
filtered_data <- all %>%
  filter(paste(symbols, orthogroup, sep = "_") %in% 
           paste(symbol_og_stats$symbols, symbol_og_stats$orthogroup, sep = "_"))

filtered_data$og_symbol <- paste(filtered_data$orthogroup, filtered_data$symbols, sep = "|")

all_species <- unique(all$species)
all_og_symbols <- unique(filtered_data$og_symbol)

########### Step3. Add WGCNA module information to the symbols (kME>0.3)
hub_files <- list.files("hub_genes/", pattern = "\\.csv$", full.names = TRUE)
hub_all <- do.call(rbind, lapply(hub_files, function(file) {
  data <- read.csv(file)
  data$species <- gsub("\\.csv$", "", basename(file))
  return(data)
}))

hub_all <- hub_all %>%
  left_join(species_class, by = "species")

high_kmE <- hub_all %>%
  filter(kmE_value > 0.3)

gene_counts <- high_kmE %>%
  group_by(gene) %>%
  summarise(
    protostome_count = sum(group == "protostome", na.rm = TRUE),
    deuterostome_count = sum(group == "deuterostome", na.rm = TRUE),
    .groups = 'drop'
  )


selected_genes <- gene_counts %>%
  filter(protostome_count >= 2 & deuterostome_count >= 6) %>%
  pull(gene)



result_df <- data.frame(
  species = all_species,
  stringsAsFactors = FALSE
)

for(og_symbol in all_og_symbols) {
  parts <- strsplit(og_symbol, "\\|")[[1]]
  orthogroup_val <- parts[1]
  symbol_val <- parts[2]
  
  species_in_og <- unique(og[[orthogroup_val]])
  temp_values <- rep(NA, length(all_species))
  
  for(j in seq_along(all_species)) {
    current_species <- all_species[j]
    
    if(current_species %in% names(og)) {
      species_genes <- og[og$Orthogroup == orthogroup_val, current_species]
      
      if(!is.na(species_genes) && species_genes != "") {
        gene_list <- unlist(strsplit(as.character(species_genes), ","))
        gene_list <- trimws(gene_list)
        
        matched_row <- filtered_data %>%
          filter(species == current_species, 
                 orthogroup == orthogroup_val, 
                 symbols == symbol_val)
        
        if(nrow(matched_row) > 0) {
          temp_values[j] <- matched_row$gene[1]
        } else {
          og_genes <- og[og$Orthogroup == orthogroup_val, current_species]
          if(!is.na(og_genes) && og_genes != "") {
            # Take the first gene from the orthogroup
            all_genes <- unlist(strsplit(as.character(og_genes), ","))
            temp_values[j] <- trimws(all_genes[1])
          } else {
            temp_values[j] <- NA
          }
        }
      } else {
        temp_values[j] <- NA
      }
    } else {
      temp_values[j] <- NA
    }
  }
  
  # Add to result dataframe
  col_name <- paste0(orthogroup_val, "_", symbol_val)
  result_df[[col_name]] <- temp_values
}

# Save the result
write.csv(result_df, "coreMarkers/CM_core_gene_set.csv", row.names = FALSE)

############ Step4. Visualize with a heatmap

rownames(result_df) = result_df$Symbol
cm_matrix <- as.matrix(result_df[3:18])


cm_matrix_numeric <- apply(cm_matrix, 2, as.numeric)
rownames(cm_matrix_numeric) <- result_df$Symbol  


pheatmap(cm_matrix_numeric,cluster_rows = F,cluster_cols = F,na_col = "white",
         color = c('#c7e0e5','#dfaa80','#db9597'),
         width = 10,height = 20,fontsize_row = 8,
         filename = 'CM_core_genes_heatmap.pdf')
