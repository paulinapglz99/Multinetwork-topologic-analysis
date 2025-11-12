#!/usr/bin/env Rscript
#Script to compute Jaccard similarity between networks' global connectivity

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "data.table",
  "optparse",
  "tools",
  "igraph",
  "pheatmap"
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL, help = "Input directory with network files", metavar = "path"),
  make_option(c("-o", "--out_file"), type = "character", default = "jaccard_heatmap.pdf", help = "Output heatmap PDF [default: %default]", metavar = "path"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$", help = "File pattern to match [default: %default]", metavar = "regex"),
  make_option(c("-t", "--type"), type = "character", default = "auto", help = "Network file type: edgelist, adjacency, auto [default: %default]", metavar = "type")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$input_dir <- ""

if (is.null(opt$input_dir)) stop("You must specify --input_dir")

#Helper functions
read_network <- function(path, type = opt$type) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    
    if (type == "edgelist") {
      g <- igraph::graph_from_data_frame(df[, 1:min(3, ncol(df))], directed = FALSE)
      
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
        g <- igraph::graph_from_data_frame(df[, 1:min(3, ncol(df))], directed = FALSE)
      }
    }
  } else {
    stop("Unsupported file format: ", ext)
  }
  
  g <- igraph::simplify(igraph::as_undirected(g), remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

#Load all networks

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("Need at least 2 files to compute similarity.")

networks <- lapply(files, read_network)
names(networks) <- basename(files)

#Compute Jaccard index on connected nodes
get_connected_nodes <- function(g) {
  nodes <- V(g)[degree(g) > 0]$name
  return(nodes)
}

connected_sets <- lapply(networks, get_connected_nodes)

n <- length(connected_sets)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(connected_sets)

for (i in 1:n) {
  for (j in i:n) {
    set_i <- connected_sets[[i]]
    set_j <- connected_sets[[j]]
    intersect_size <- length(intersect(set_i, set_j))
    union_size <- length(union(set_i, set_j))
    jaccard <- if (union_size > 0) intersect_size / union_size else 0
    jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- jaccard
  }
}

#Plot heatmap
pheatmap(jaccard_matrix, 
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Jaccard Similarity of Network Connectivity",
         filename = opt$out_file)

message("Saved heatmap to: ", opt$out_file)