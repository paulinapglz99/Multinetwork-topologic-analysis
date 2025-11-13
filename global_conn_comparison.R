#!/usr/bin/env Rscript
# Compute Jaccard similarity between networks (based on edge sets)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

pacman::p_load(
  "data.table", "optparse", "tools",
  "igraph", "pheatmap"
)

# -------------------------------
# Argument parsing
# -------------------------------
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Directory with network files", metavar = "path"),
  make_option(c("-d", "--out_dir"), type = "character", default = "jaccard_output", help = "Output directory [default: %default]", metavar = "path"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$", help = "File name pattern to match [default: %default]", metavar = "regex"),
  make_option(c("-t", "--type"), type = "character", default = "auto", help = "Network file type: edgelist, adjacency, auto [default: %default]", metavar = "type")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
opt$out_dir <- "/STORAGE/csbig/networks_final/fomo_networks/results_conn_comparisons"
opt$input_dir <- "/STORAGE/csbig/networks_final/fomo_networks/"
opt$pattern <- "\\.tsv$"

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

#Jaccard function
jaccard_edges <- function(g1, g2) {
  inter <- igraph::intersection(g1, g2, byname = TRUE)
  uni <- igraph::union(g1, g2, byname = TRUE)
  return(length(E(inter)) / length(E(uni)))
}

#Load networks
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("At least 2 network files are required.")

networks <- lapply(files, read_network)
names(networks) <- basename(files)

#Compute Jaccard similarity matrix
n <- length(networks)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(networks)

for (i in 1:n) {
  for (j in i:n) {
    sim <- jaccard_edges(networks[[i]], networks[[j]])
    jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- sim
  }
}

#Save results
csv_out <- file.path(opt$out_dir, "jaccard_edge_matrix.csv")
pdf_out <- file.path(opt$out_dir, "jaccard_edge_heatmap.pdf")
fwrite(as.data.table(jaccard_matrix, keep.rownames = TRUE), csv_out)

#Plot
pheatmap(jaccard_matrix,
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Edge Jaccard Similarity",
         filename = pdf_out, 
         color = hcl.colors(50, "BluYl")
         )

message("Jaccard matrix saved to: ", csv_out)
message("Heatmap saved to: ", pdf_out)

#END