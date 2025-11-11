#mi pregunta es, existe algun gen hub que sea un hub diferenciador, es decir, solo presente como hub en la enfermedad de Alzheimer
#y que este se conserve a lo largo de las redes en las diferentes regiones anatomicas?

#Get packages
pacman::p_load("tidyverse",
               "clusterProfiler", 
               "gridExtra", 
               "biomaRt",
               "cowplot")

# 1. Cargar archivos
data_dir <- "~/Desktop/local_work/results_topos"
files <- list.files(data_dir, pattern = "_nodes_summary.csv$", full.names = TRUE)

# 2. Leer archivos y agregar metadatos
node_data <- map_df(files, function(file) {
  df <- fread(file)
  df$file <- basename(file)
  df$Region <- str_remove(df$file, "_counts_.*")
  df$Phenotype <- ifelse(str_detect(df$file, "_AD_"), "AD", "Control")
  return(df)
})

# 3. Definir hubs: top 5% degree por red
node_data <- node_data %>%
  group_by(file) %>%
  mutate(
    hub_degree = degree >= quantile(degree, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

# 4. Tabla: gene × region × phenotype
hubs_summary <- node_data %>%
  filter(hub_degree) %>%
  group_by(node, Region, Phenotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(hub = TRUE)

# 5. Pivot: gene × Region × Phenotype
hubs_pivot <- hubs_summary %>%
  pivot_wider(names_from = Phenotype, values_from = hub, values_fill = FALSE)

# 6. Filtrar genes AD que no son hubs en ningún control
exclusive_AD_hubs <- hubs_pivot %>%
  filter(AD == TRUE & Control == FALSE) %>%
  group_by(node) %>%
  summarise(n_regions_AD = n(), regions = paste(Region, collapse = ",")) %>%
  arrange(desc(n_regions_AD))

# Resultado final
print(exclusive_AD_hubs)

# 8. Crear heatmap gene × region
heatmap_matrix <- hubs_pivot %>%
  filter(AD == TRUE & Control == FALSE) %>%
  select(node, Region) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Region, values_from = value, values_fill = 0) %>%
  column_to_rownames("node") %>%
  as.matrix()

pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "AD-Exclusive Hub Genes (Any Metric)",
         fontsize_row = 7,
         color = colorRampPalette(c("gray90", "firebrick3"))(100))

# Agregar columna 'Sum' con número de regiones donde el gen es hub en AD
heatmap_matrix.c <-cbind(heatmap_matrix, Sum = rowSums(heatmap_matrix)) %>% 
 as.data.frame()

table(heatmap_matrix.c$Sum)

wide_hubs <-  heatmap_matrix.c %>% filter(Sum >= 3)

rownames(wide_hubs)
w<-enrichGO(gene = rownames(wide_hubs),
         OrgDb = "org.Hs.eg.db", 
         keyType = 'ENSEMBL',
         readable = TRUE,
         #universe = universe, 
         ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
         pvalueCutoff = 0.05, 
         qvalueCutoff = 0.10)

# 9. Exportar tabla larga para análisis funcional
hub_genes_long <- node_data %>%
  filter(hub_total) %>%
  distinct(node, Region, Phenotype)

fwrite(hub_genes_long, "all_hub_genes_long.csv")
