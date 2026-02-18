#!/usr/bin/env Rscript

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
               "optparse", 
               "cowplot", 
               "ggplotify")

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
#input_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/"
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

################### ################### ################### ################### 

#QUESTION 1: Are there genes that are exclusive hubs in AD and that are conserved across different brain regions?

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

hubs_degree_pivot_filter_AD <- hubs_degree_pivot %>%
  filter(symbol %in% exclusive_AD_hubs_degree_filter$symbol)

#Heatmap gene × region
heatmap_matrix <- hubs_degree_pivot_filter_AD %>%
  filter(AD == TRUE & Control == FALSE) %>%
  dplyr::select(symbol, Region) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Region, 
              values_from = value, 
              values_fill = 0) %>%
  column_to_rownames("symbol") %>%
  as.matrix()

#Reorder rows of the heatmap according to gene_order
heatmap_matrix <- heatmap_matrix[exclusive_AD_hubs_degree_filter$symbol, ]

#Plot heatmap
ad.ht <- pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         legend = FALSE,
         main = NA,
         fontsize_row = 10,
         color = colorRampPalette(c("gray95", "firebrick"))(100), 
         file = "heatmap_hubs.jpeg")

ad.ht <-ggplotify::as.ggplot(ad.ht)
  
################### AHORA CON LOS CONTROLES ##########################

exclusive_Control_hubs_degree <- hubs_degree_pivot %>%
  filter(Control == TRUE & AD == FALSE) %>%
  group_by(symbol) %>%
  summarise(
    n_regions_Control = n(),
    regions = paste(Region, collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(n_regions_Control))

exclusive_Control_hubs_degree_filter <- exclusive_Control_hubs_degree %>%
  filter(n_regions_Control >= 2) %>%
  arrange(desc(n_regions_Control))

hubs_degree_pivot_control_filter <- hubs_degree_pivot %>%
  filter(symbol %in% exclusive_Control_hubs_degree_filter$symbol)

heatmap_matrix_control <- hubs_degree_pivot_control_filter %>%
  filter(Control == TRUE & AD == FALSE) %>%
  dplyr::select(symbol, Region) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Region,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("symbol") %>%
  as.matrix()
heatmap_matrix_control <- heatmap_matrix_control[
  exclusive_Control_hubs_degree_filter$symbol,
]

#Vis

ctrl.ht <-pheatmap(
  heatmap_matrix_control,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  legend = FALSE,
  main = NA,
  fontsize_row =4,
  color = colorRampPalette(c("gray95", "cornflowerblue"))(100),
  file = "heatmap_hubs_Control.jpeg"
)

ctrl.ht <-ggplotify::as.ggplot(ctrl.ht)

#Vis both

hts <- cowplot::plot_grid(
  ad.ht, ctrl.ht,
  labels = c("A", "B"),
  label_size = 10,
  ncol = 2
)

#Vis

hts

#Save plots

ggsave(filename = "hub_heatmaps.jpeg",
       plot = hts,
       width = 4,
       height = 6)

################### ################### ################### ################### 

#QUESTION 2: What molecular functions are associated with these genes?

#Enrichment of genes
wide_hubs_AD <- unique(hubs_degree_pivot_filter$symbol)

w_enrich_AD <- enrichGO(gene = wide_hubs_AD,
         OrgDb = "org.Hs.eg.db",
         keyType = 'SYMBOL',
         readable = TRUE,
         #universe = universe,
         ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.10)

as.data.frame(w_enrich_AD)

#dotplot(w_enrich_edox, showCategory=30)

cnetplot <- cnetplot(w_enrich_AD, node_label_size = NULL)
cnetplot

################### AHORA CON LOS CONTROLES ##########################

#Enrichment of genes
wide_hubs_ctr <-unique(hubs_degree_pivot_control_filter$symbol)

w_enrich_ctr <-enrichGO(gene = wide_hubs_ctr,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'SYMBOL',
                    readable = TRUE,
                    #universe = universe,
                    ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.10)

as.data.frame(w_enrich)

#dotplot(w_enrich_edox, showCategory=30)

cnetplot <- cnetplot(w_enrich_ctr, node_label_size = NULL)
cnetplot

#QUESTION 4: Are these genes correlated between regions?
#See if hubs repeat in a coordinated manner between regions.
#Calculate the correlation matrix between regions.
# 
# # Lista de archivos de redes AD
# ad_net_files <- list.files(opt$input_dir, pattern = "_AD_.*\\.tsv$", full.names = TRUE)
# 
# # Extraer metadatos de cada archivo
# meta <- map_df(ad_net_files, extract_info)
# meta$region <- gsub("^(Mayo_|ROSMAP_)", "", meta$region)
# meta$Network <- paste0(meta$region, "_", meta$phenotype)
# 
# # Leer redes
# ad_networks <- map(ad_net_files, ~ read_network(.x, type = "auto"))
# names(ad_networks) <- meta$region
# 
# # Genes persistentes
# wide_hub_genes <- exclusive_AD_hubs_degree_filter$symbol
# 
# # Filtrar redes por genes persistentes
# ad_networks_persistent <- map(ad_networks, function(g) {
#   igraph::induced_subgraph(g, vids = igraph::V(g)[name %in% wide_hub_genes])
# })

################### ################### ################### ################### 

#QUESTION 5: Where do these genes appear in the network?
#Are these hub genes grouped into specific modules? Are they scattered?
#----------------------------
# Helper function
#----------------------------
get_membership_comparison <- function(gene_vector, hub_label) {
  
  ad_df <- node_data %>%
    filter(
      Phenotype == "AD",
      symbol %in% gene_vector
    ) %>%
    dplyr::distinct(symbol, Region, membership) %>%
    dplyr::rename(membership_AD = membership)
  
  ctr_df <- node_data %>%
    filter(
      Phenotype == "Control",
      symbol %in% gene_vector
    ) %>%
    dplyr::distinct(symbol, Region, membership) %>%
    dplyr::rename(membership_Control = membership)
  
  full_join(ad_df, ctr_df,
            by = c("symbol", "Region")) %>%
    dplyr::mutate(hub_type = hub_label) %>%
    arrange(symbol, Region)
}

#----------------------------
# Apply to both hub sets
#----------------------------

membership_all <- bind_rows(
  get_membership_comparison(wide_hubs_AD,  "AD_exclusive"),
  get_membership_comparison(wide_hubs_ctr, "Control_exclusive")
)

membership_all

#Visualization
sankey_data <- membership_all %>%
  filter(
    !is.na(membership_Control),
    !is.na(membership_AD)
  ) %>%
  mutate(
    membership_Control = as.factor(membership_Control),
    membership_AD = as.factor(membership_AD)
  )

#Plot sankey with facet-wrap
igv_palette <- c(
  "#1f78b4",  # blue
  "#33a02c",  # green
  "#e31a1c",  # red
  "#ff7f00",  # orange
  "#6a3d9a",  # purple
  "#b15928",  # brown
  "#a6cee3",  # light blue
  "#b2df8a",  # light green
  "#fb9a99",  # pink
  "#fdbf6f",  # light orange
  "#cab2d6",  # lavender
  "#ffff99"   # yellow
)

sankey_data$membership_Control <- as.factor(sankey_data$membership_Control)
sankey_data$membership_AD <- as.factor(sankey_data$membership_AD)

sankey.p <- ggplot(sankey_data,
       aes(axis1 = membership_Control,
           axis2 = membership_AD)) +
  
  geom_alluvium(aes(fill = membership_Control),
                width = 0.25,
                alpha = 0.9,
                knot.pos = 0.35) +
  
  geom_stratum(width = 0.35,
               fill = "white",
               color = "black",
               size = 0.3) +
  
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3,
            fontface = "bold") +
  
  scale_fill_manual(values = igv_palette) +
  
  scale_x_discrete(limits = c("Control", "AD"),
                   expand = c(.08, .08)) +
  
  facet_wrap(~ Region, nrow = 1) +
  
  labs(
    title = "Module Transitions of Exclusive Hubs",
    subtitle = "Control → AD modular reassignment",
    x = NULL,
    y = "Hub count",
    fill = "Control Module"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )
#Vis

sankey.p

################### ################### ################### ################### 
#que funciones tienen

enrich_files <- list.files("~/Desktop/local_work/fomo_networks/results_topos_louvain/results_comm/",
                           pattern = "_enrichment\\.csv$",
                           full.names = TRUE)

enrich_data <- lapply(enrich_files, vroom::vroom)

enrich_df <- map2_dfr(
  enrich_files,
  enrich_data,
  ~ .y %>%
    mutate(
      file = basename(.x),
      Region = str_extract(file, "(?<=Mayo_)[A-Z]+"),
      Phenotype = ifelse(str_detect(file, "_AD_"), "AD", "Control")
    )
)

enrich_top <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Region, Phenotype, CommunityID) %>%
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Region, Phenotype, CommunityID, Description)

membership_all_clean <- membership_all %>%
  dplyr::rename(
    CommunityID_AD = membership_AD,
    CommunityID_Control = membership_Control
  )

membership_annotated <- membership_all_clean %>%
  
  # Función del módulo en AD
  left_join(
    enrich_top %>%
      filter(Phenotype == "AD") %>%
      dplyr::rename(
        CommunityID_AD = CommunityID,
        Function_AD = Description
      ),
    by = c("Region", "CommunityID_AD")
  ) %>%
  
  # Función del módulo en Control
  left_join(
    enrich_top %>%
      filter(Phenotype == "Control") %>%
      dplyr::rename(
        CommunityID_Control = CommunityID,
        Function_Control = Description
      ),
    by = c("Region", "CommunityID_Control")
  )

#END