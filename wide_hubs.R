#mi pregunta es, existe algun gen hub que sea un hub diferenciador, es decir, solo presente como hub en la enfermedad de Alzheimer
#y que este se conserve a lo largo de las redes en las diferentes regiones anatomicas?
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("mygene")
library(mygene)

#Get packages
pacman::p_load("tidyverse",
               "clusterProfiler", 
               "gridExtra", 
               "biomaRt",
               "cowplot", 
               "vroom", 
               "pheatmap", "org.Hs.eg.db",
               "AnnotationDbi")

#Functions


#Get files
input_dir <- "/STORAGE/csbig/networks_final/fomo_networks/results_topos_infomap/"
# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}
files <- list.files(input_dir, pattern = "_nodes_summary.csv$", full.names = TRUE)
#Get data
node_data <- map_df(files, function(file) {
  df <- vroom::vroom(file)
  df$file <- basename(file)
  #df$Region <- str_remove(df$file, "_counts_.*")
  df$Region <- sub("_counts_.*", "", df$file)
  df$Region <- gsub("^(Mayo_|ROSMAP_)", "", df$Region)
  df$Phenotype <- ifelse(str_detect(df$file, "_AD_"), "AD", "Control")
  df$network <- paste0(df$Region, "_", df$Phenotype) 
  return(df)
})

#PREGUNTA 1:¿Existen genes que son hubs exclusivos en AD y que se conservan a lo largo de diferentes regiones cerebrales?

#Define hubs: top 10% degree and pagerank
node_data <- node_data %>%
  group_by(file) %>%
  mutate(
    hub_degree = degree >= quantile(degree, 0.90, na.rm = TRUE)
  ) %>%
  mutate(
    hub_pagerank = pagerank >= quantile(pagerank, 0.90, na.rm = TRUE)
  ) %>%
  ungroup()

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = node_data$node,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

symbols.df <- tibble(
  node   = names(symbols),
  symbol = unname(symbols)
)

node_data <- node_data %>%
  left_join(symbols.df, by = "node") %>%
  mutate(symbol = coalesce(symbol, node))

#Build table gene × region × phenotype
hubs_degree.s <- node_data %>%
  filter(hub_degree & hub_pagerank) %>%
  #filter(hub_pagerank) %>%
  group_by(symbol, Region, Phenotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(hub = TRUE)

#Pivot gene × Region × Phenotype
hubs_degree_pivot <- hubs_degree.s %>%
  pivot_wider(names_from = Phenotype,
              values_from = hub, 
              values_fill = FALSE)

#Filter genes that are hubs in any AD net but not hubs in any control
exclusive_AD_hubs_degree <- hubs_degree_pivot %>%
  filter(AD == TRUE & Control == FALSE) %>%
  group_by(symbol) %>%
  summarise(n_regions_AD = n(), regions = paste(Region, collapse = ",")) %>%
  arrange(desc(n_regions_AD))

exclusive_AD_hubs_degree_filter <- exclusive_AD_hubs_degree %>% 
  filter(n_regions_AD >=2)

hubs_degree_pivot_filter <- hubs_degree_pivot %>%
  filter(symbol %in% exclusive_AD_hubs_degree_filter$symbol)

#Heatmap gene × region
heatmap_matrix <- hubs_degree_pivot_filter %>%
  filter(AD == TRUE & Control == FALSE) %>%
  dplyr::select(symbol, Region) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Region, 
              values_from = value, 
              values_fill = 0) %>%
  column_to_rownames("symbol") %>%
  as.matrix()

#Plot heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         legend = FALSE,
         main = NA,
         fontsize_row = 10,
         color = colorRampPalette(c("gray90", "royalblue4"))(100), 
         file = "heatmap_hubs.jpeg")

#PREGUNTA 2: ¿Qué funciones moleculares están asociadas a estos genes?

#Enrichment of genes
wide_hubs <-unique(hubs_degree_pivot_filter$symbol)

w_enrich <-enrichGO(gene = wide_hubs,
         OrgDb = "org.Hs.eg.db", 
         keyType = 'SYMBOL',
         readable = TRUE,
         #universe = universe, 
         ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
         pvalueCutoff = 0.05, 
         qvalueCutoff = 0.10)

as.data.frame(w_enrich)

#dotplot(w_enrich_edox, showCategory=30)
  
cnetplot <- cnetplot(w_enrich_edox, node_label_size = NULL)
cnetplot

#PREGUNTA 3: ¿Cuáles genes son los hubs más persistentes en AD?

exclusive_AD_hubs_degree_filter %>%
  top_n(20, n_regions_AD) %>%
  ggplot(aes(x = reorder(symbol, n_regions_AD), y = n_regions_AD)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() +
  labs(title = "Top AD-Exclusive Hubs by Region Recurrence",
       x = "Gene", y = "Number of AD Regions") +
  theme_minimal()

#PREGUNTA 4: ¿Estos genes están correlacionados entre regiones?
#Ver si los hubs se repiten de forma coordinada entre regiones

# Jaccard similarities between AD hub sets per region
library(reshape2)
hub_sets <- hubs_degree_pivot_filter %>%
  filter(AD == TRUE) %>%
  group_by(Region) %>%
  summarise(genes = list(unique(node))) %>%
  deframe()

# Compute Jaccard Index
jaccard_index <- function(a, b) length(intersect(a, b)) / length(union(a, b))

region_names <- names(hub_sets)
jac_mat <- matrix(0, nrow = length(region_names), ncol = length(region_names),
                  dimnames = list(region_names, region_names))

for (i in region_names) {
  for (j in region_names) {
    jac_mat[i, j] <- jaccard_index(hub_sets[[i]], hub_sets[[j]])
  }
}

pheatmap(jac_mat, main = "Jaccard Similarity of Hub Genes Across AD Regions",
         color = colorRampPalette(c("white", "darkred"))(100))

#PREGUNTA 5: ¿Dónde aparecen estos genes en la red?
#¿Estos genes hubs se agrupan en módulos específicos? ¿Están dispersos?

wide_hub_genes <- exclusive_AD_hubs_degree_filter$node

# Solo redes AD
wide_nodes_AD <- node_data %>%
  filter(Phenotype == "AD") %>%
  filter(node %in% wide_hub_genes)

modules_distribution <- wide_nodes_AD %>%
  group_by(network, membership_infomap) %>%
  summarise(n_genes = n(),
            genes = paste(node, collapse = ","),
            .groups = "drop")

# Asegúrate de que tiene la columna 'membership_infomap'
stopifnot("membership_infomap" %in% colnames(nodes_AD))

module_summary <- modules_distribution %>%
  group_by(network) %>%
  summarise(
    total_hubs = sum(n_genes),
    n_modules_with_hubs = n(),
    max_module_size = max(n_genes),
    top_module_prop = max(n_genes) / sum(n_genes)
  )

genes_modules <- nodes_AD %>%
  dplyr::select(node, network, membership_infomap)

genes_modules_matrix <- genes_modules %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = network, values_from = value, values_fill = 0) %>%
  column_to_rownames("node") %>%
  as.matrix()

pheatmap(genes_modules_matrix,
         main = "Presencia de genes hub AD-exclusivos en redes",
         color = colorRampPalette(c("white", "red3"))(100))

#Grand panel


