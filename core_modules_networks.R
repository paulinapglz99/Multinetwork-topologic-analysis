#!/usr/bin/env Rscript

#core_modules.R

#Packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  "data.table", 
  "tidyverse",
  "tools",
  "optparse",
  "igraph"
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directory with files *_nodes_summary.csv", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_core_modules",
              help = "Output directory", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex to select archives", metavar = "character"), 
  make_option(c("-r", "--exclusive_modules_file"), type = "character", default = "./modules_exclusivos.csv",
              help = "File with exclusive modules", metavar = "character"),
  make_option(c("-n", "--network_dir"), type = "character", default = "./networks",
              help = "Network directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

opt$input_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain"
opt$output_dir <- "~/Desktop/local_work/fomo_networks/results_core_modules_networks"
opt$exclusive_modules_file <- "~/Desktop/local_work/fomo_networks/results_core_modules/exclusive_modules_AD_CTRL.csv"
opt$network_dir <- "~/Desktop/local_work/fomo_networks/translated_graphs/"

#Aux functions
extract_info <- function(filename) {
  base <- basename(filename)
  region <- sub("_counts_.*", "", base)
  phenotype <- ifelse(grepl("_AD_", base), "AD", "Control")
  return(data.frame(filename = filename, region = region, phenotype = phenotype))
}

load_modules <- function(file) {
  df <- fread(file)
  stopifnot(all(c("node", "membership") %in% names(df)))
  region <- sub("_nodes_summary.*", "", basename(file))
  phenotype <- ifelse(grepl("_AD_", file), "AD", "Control")
  df$filename <-basename(file)
  df <- df %>%
    mutate(
      Region = str_split(filename, "_", simplify = TRUE)[, 2],
      Phenotype = str_extract(filename, "(?<=_counts_)[^_]+"),
      unique_module_id = paste(Region, Phenotype, membership, sep = "_")
    )
  
  #df$unique_module_id <- paste(region, phenotype, df$membership, sep = "_")
  df
}

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

annotate_graph_with_modules <- function(graph, module_df) {
  V(graph)$name <- as.character(V(graph)$name)
  node_data <- module_df[match(V(graph)$name, module_df$node), ]
  V(graph)$module <- node_data$membership
  V(graph)$unique_module_id <- node_data$unique_module_id
  V(graph)$exclusive_status <- node_data$exclusive_status
  return(graph)
}

export_filtered_graphs <- function(graph, status, output_dir, tag) {
  filtered_nodes <- V(graph)[exclusive_status == status]
  if (length(filtered_nodes) == 0) return(NULL)
  subgraph <- induced_subgraph(graph, filtered_nodes)
  fname <- file.path(output_dir, paste0("filtered_", tag, "_", status, ".graphml"))
  write_graph(subgraph, fname, format = "graphml")
}

#Create directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

#Get module info
module_files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
modules_df <- rbindlist(lapply(module_files, load_modules))

#Get file of exclusive modules
exclusive_modules <- vroom::vroom(opt$exclusive_modules_file)
exclusive_modules <- exclusive_modules %>%
  mutate(Module = str_remove(Module, "^[^_]+_")) %>% 
  rename(is_exclusive = Phenotype) %>% 
  dplyr::select(-Region)

modules_df <- modules_df %>% left_join(exclusive_modules, by = c("unique_module_id"= "Module"))

#Get network names
network_names <- basename(network_files) %>%
  str_replace(".graphml$", "") %>%
  str_extract_all("[^_]+") %>%
  map_chr(~ paste(.x[2], .x[4], sep = "_"))

#Get networks
network_files <- list.files(opt$network_dir, pattern = "\\.graphml", full.names = TRUE)
networks <- lapply(network_files, read_network)
names(networks) <- network_names

graph <- networks[[1]]

#annotate_graph_with_modules <- function(graph, modules_df) {
  V(graph)$name <- as.character(V(graph)$name)
  node_data <- modules_df[match(V(graph)$name, modules_df$node), ]
  #node_data <- modules_df %>% slice(match(V(graph)$name, node))
  V(graph)$module <- node_data$membership
  V(graph)$unique_module_id <- node_data$unique_module_id
  V(graph)$exclusive_status <- node_data$exclusive_status
  return(graph)
#}

for (netfile in network_files) {
  cat("Procesando red:", netfile, "\n")
  graph <- read_network(netfile)
  graph <- annotate_graph_with_modules(graph, modules_df)
  
  # # Nombre base para guardar
  # tag <- tools::file_path_sans_ext(basename(netfile))
  # 
  # # Guardar subgrafos filtrados
  # export_filtered_graphs(graph, "AD_exclusive", opt$output_dir, tag)
  # export_filtered_graphs(graph, "Control_exclusive", opt$output_dir, tag)
}
