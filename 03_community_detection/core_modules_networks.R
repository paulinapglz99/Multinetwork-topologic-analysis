#!/usr/bin/env Rscript

#core_modules_networks.R

#Aim: Find connected modules specific to the phenotypes
#Strategy:
#1. Label the nodes to their module
#2. Label nodes by their classification
#3. Filter nodes that have the label "AD_exclusive" or "Control_exclusive"
#4. Export those subgraphs

#Packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  "data.table", 
  "tidyverse",
  "tools",
  "optparse",
  "igraph")

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
              help = "Network directory", metavar = "character"),
  make_option(c("-t", "--type"), type = "character", default = "auto",
              help = "Network format: edgelist, adjacency, or auto", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

opt$input_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain"
opt$output_dir <- "~/Desktop/local_work/fomo_networks/results_core_modules_networks"
opt$exclusive_modules_file <- "~/Desktop/local_work/fomo_networks/results_core_modules/module_classification.csv"
opt$network_dir <- "~/Desktop/local_work/fomo_networks/translated_graphs/"
setwd(opt$output_dir)

#Aux functions
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

#Annotation of the nodes with their module and module classification
annotate_graph_with_modules <- function(graph, modules_df, region, phenotype) {
  #Filter modules that correspond to the region and phenotype
  filtered_modules <- modules_df %>%
    filter(Region == region, Phenotype == phenotype) %>%
    distinct(node, .keep_all = TRUE)  # Elimina duplicados si existieran
  
  #Connect network nodes with their module attributes
  node_df <- tibble(name = V(graph)$name) %>%
    left_join(filtered_modules, by = c("name" = "node")) %>%
    mutate(
      membership = replace_na(membership, -1),
      unique_module_id = replace_na(unique_module_id, "none"),
      exclusive_status = replace_na(exclusive_status, "Not_exclusive")
    )
  
  V(graph)$module <- node_df$membership
  V(graph)$unique_module_id <- node_df$unique_module_id
  V(graph)$exclusive_status <- node_df$exclusive_status
  
  return(graph)
}

#This validates if the graph was well annotated
validate_annotation <- function(g, modules_subset,
                                node_col = "node",    #change if your column has a different name
                                module_col = "membership") {
  #Build a dataframe with nodes and nodes assigned in the graph
  df <- data.frame(node = V(g)$name,
    assigned = V(g)$module)
  #Then merge them with the original modules
  merged <- left_join(df,
    modules_subset,
    by = c("node" = node_col))
  #And compare assigned module vs expected module
  return(table(merged$assigned == merged[[module_col]]))
}

filter_graphs <- function(graph, status, output_dir, tag) {
  #Check if attribute exists
  if (!"exclusive_status" %in% vertex_attr_names(graph)) {
    warning("The graph does not have the “exclusive_status” attribute. It cannot be filtered.")
    return(NULL)
  }
  
  #Filter nodes that meet the status
  filtered_nodes <- V(graph)[V(graph)$exclusive_status == status]
  #If there's no nodes with this statis, exit
  if (length(filtered_nodes) == 0) {
    message(paste("There are no nodes with status:", status, "in the network", tag))
    return(NULL)
  }
  #Make subgraph
  subgraph <- induced_subgraph(graph, filtered_nodes)
  fname <- file.path(output_dir, paste0("filtered_", tag, "_", status, ".graphml"))
  write_graph(subgraph, fname, format = "graphml")
  message(paste("✓ Subgrafo filtrado guardado:", fname))
}

#Checker
checadore <- function(graph, name) {
  message("Checking: ", name)
  
  n_nodes <- vcount(graph)
  n_edges <- ecount(graph)
  
  status_counts <- table(V(graph)$exclusive_status)
  module_counts <- length(unique(V(graph)$module))
  missing_status <- sum(is.na(V(graph)$exclusive_status) | V(graph)$exclusive_status == "Not_exclusive")
  
  cat("Nodes:", n_nodes, "\n")
  cat("Edges:", n_edges, "\n")
  cat("Unique modules:", module_counts, "\n")
  cat("• Count by 'exclusive_status':\n")
  print(status_counts)
  cat("Non-exclusive nodes", missing_status, "\n")
  cat("───────────────\n")
}

#Create directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

#Get module info
module_files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
modules_df <- rbindlist(lapply(module_files, load_modules))

#Assign exclusive_status to each node
exclusive_modules <- vroom::vroom(opt$exclusive_modules_file)

exclusive_modules <- exclusive_modules %>%
  mutate(Module = str_extract(Module, "[^_]+$"),
    unique_module_id = paste(Region, Phenotype, Module, sep = "_"))

#Assign exclusive_status to each node in modules_df
modules_df <- modules_df %>%
  left_join(
    exclusive_modules %>% dplyr::select(unique_module_id, Classification),
    by = "unique_module_id") %>%
  mutate(exclusive_status = coalesce(Classification, "Not_exclusive")) #Si el módulo tiene una clasificación (por ejemplo "AD_exclusive"), esa se usa. 
  #Si no tiene clasificación (es decir, no estaba en exclusive_modules), se asigna "Not_exclusive".

#Get networks
network_files <- list.files(opt$network_dir, pattern = "\\.graphml", full.names = TRUE)
networks <- lapply(network_files, read_network)
#Get network names
network_names <- basename(network_files) %>%
  str_replace(".graphml$", "") %>%
  str_extract_all("[^_]+") %>%
  map_chr(~ paste(.x[2], .x[4], sep = "_"))
names(networks) <- network_names

#Record each network with module information and unique status
for (i in seq_along(networks)) {
  g <- networks[[i]]
  name <- names(networks)[i]  # Ejemplo: "CRB_AD"
  #Get region and phenotype
  parts <- unlist(str_split(name, "_"))
  region <- parts[1]
  phenotype <- parts[2]
  
  message("Annotating graph: ", name, " (", region, ", ", phenotype, ")")
  
  g_annotated <- annotate_graph_with_modules(g, modules_df, region, phenotype)
  
  networks[[i]] <- g_annotated
}

#Check if attribute is now in the networks
vertex_attr_names(networks[[1]])

#Check if nets are well annotated
for (i in seq_along(networks)) {
  g <- networks[[i]]
  name <- names(networks)[i]
  checadore(g, name)
}

#Finally, filter and save subgraphs for DD_exclusive y Control_exclusive
for (i in seq_along(networks)) {
  g <- networks[[i]]
  name <- names(networks)[i]
  
  message("Filtering and saving subgraphs for: ", name)
  
  filter_graphs(g, status = "AD_exclusive", output_dir = opt$output_dir, tag = name)
  filter_graphs(g, status = "Control_exclusive", output_dir = opt$output_dir, tag = name)
}

#END