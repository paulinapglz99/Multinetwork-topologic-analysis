#My question is, is there any gene hub that is a differentiating hub,
#that is, only present as a hub in Alzheimer's disease
#and that is conserved throughout the networks in different anatomical regions?

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#Get packages
pacman::p_load("tidyverse",
               "clusterProfiler", 
               "gridExtra", 
               "biomaRt",
               "cowplot", 
               "vroom", 
               "pheatmap", 
               "org.Hs.eg.db",
               "AnnotationDbi",
               "optparse")

#Functions
#Read_network function
read_network <- function(path, type = opt$type) {
  message("DEBUG - received type: ", paste0(type, collapse = ","))
  type <- match.arg(type)
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    
    if (type == "edgelist") {
      g <- if (ncol(df) >= 3) igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
      else igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      
    } else if (type == "adjacency") {
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      
    } else if (type == "auto") {
      looks_like_adj <- all(sapply(df[,-1, drop=FALSE], is.numeric))
      if (looks_like_adj) {
        rownames(df) <- make.unique(as.character(df[, 1]))
        mat <- as.matrix(df[, -1, drop=FALSE])
        storage.mode(mat) <- "numeric"
        g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      } else {
        g <- if (ncol(df) >= 3) igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
        else igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
    }
    
  } else {
    stop("Unsupported file format: ", ext)
  }
  
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

#Extract metadata
extract_info <- function(filename) {
  base <- basename(filename)
  if (grepl("_counts_", base)) {
    region <- sub("_counts_.*", "", base)
    phenotype <- ifelse(grepl("_AD_", base), "AD", "Control")
    return(data.frame(filename = filename, region = region, phenotype = phenotype))
  } else {
    return(NULL)
  }
}

#Get files
input_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/"
#"/STORAGE/csbig/networks_final/fomo_networks/results_topos_infomap/"

#Argument parsing
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Directory with network files", metavar = "path"),
  make_option(c("-d", "--out_dir"), type = "character", default = "jaccard_output", help = "Output directory [default: %default]", metavar = "path"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$", help = "File name pattern to match [default: %default]", metavar = "regex"),
  make_option(c("-t", "--type"), type = "character", default = "auto", help = "Network file type: edgelist, adjacency, auto [default: %default]", metavar = "type")
)

opt <- parse_args(OptionParser(option_list = option_list))
opt$out_dir <- "~/Desktop/local_work/fomo_networks/results_wide_hubs"
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
opt$input_dir <- "~/Desktop/local_work/fomo_networks/fomo_networks/"
opt$pattern <- "\\.tsv$"

setwd(opt$out_dir)

#Define function to convert from ENSMBL to SYMBOL
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
  filter(n_regions_AD >=2) %>%
  arrange(desc(n_regions_AD))

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

# Reordenar filas del heatmap según gene_order
heatmap_matrix <- heatmap_matrix[exclusive_AD_hubs_degree_filter$symbol, ]

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

# Asegura que los genes estén ordenados igual que en el heatmap
exclusive_AD_hubs_degree_filter$symbol <- factor(
  exclusive_AD_hubs_degree_filter$symbol, 
  levels = rev(exclusive_AD_hubs_degree_filter$symbol)
)

persistent <- exclusive_AD_hubs_degree_filter %>%
  #top_n(20, n_regions_AD) %>%
  ggplot(aes(x = symbol, y = n_regions_AD)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  geom_text(aes(label = regions), 
            hjust = 1.2, 
            size = 4,
            color="white" 
              ) +  # Mueve etiqueta a la derecha
  coord_flip() +
  labs(title = "Top AD-Exclusive Hubs by Region Recurrence",
       x = " ", y = "Number of Regions") +
  theme_minimal()

ggsave("persistent_hubs.jpeg", 
       persistent)

#PREGUNTA 4: ¿Estos genes están correlacionados entre regiones?
#Ver si los hubs se repiten de forma coordinada entre regiones
# Calcular la matriz de correlación entre regiones

# Lista de archivos de redes AD
ad_net_files <- list.files(opt$input_dir, pattern = "_AD_.*\\.tsv$", full.names = TRUE)

# Extraer metadatos de cada archivo
meta <- map_df(ad_net_files, extract_info)
meta$region <- gsub("^(Mayo_|ROSMAP_)", "", meta$region)
meta$Network <- paste0(meta$region, "_", meta$phenotype)

# Leer redes
ad_networks <- map(ad_net_files, ~ read_network(.x, type = "auto"))
names(ad_networks) <- meta$region

# Genes persistentes
wide_hub_genes <- exclusive_AD_hubs_degree_filter$symbol

# Filtrar redes por genes persistentes
ad_networks_persistent <- map(ad_networks, function(g) {
  igraph::induced_subgraph(g, vids = igraph::V(g)[name %in% wide_hub_genes])
})


#PREGUNTA 5: ¿Dónde aparecen estos genes en la red?
#¿Estos genes hubs se agrupan en módulos específicos? ¿Están dispersos?

wide_hub_genes <- exclusive_AD_hubs_degree_filter$symbol



