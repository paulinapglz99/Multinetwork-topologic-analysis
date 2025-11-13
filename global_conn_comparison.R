#!/usr/bin/env Rscript
# Compute Jaccard similarity between edge sets of networks

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

opt$input_dir <- "/STORAGE/csbig/networks_final/networks_filtered/"
opt$out_dir <- "/STORAGE/csbig/networks_final/networks_filtered/results_global_comparisons"
opt$pattern <- "\\.tsv$"

#read_network function

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

# -------------------------------
# Read and process networks
# -------------------------------
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("At least 2 network files are required.")

networks <- lapply(files, read_network)
names(networks) <- basename(files)

# -------------------------------
# Extract sets of unique edges
# -------------------------------
get_edge_set <- function(g) {
  edgelist <- igraph::as_edgelist(g)
  apply(edgelist, 1, function(e) paste(sort(e), collapse = "|"))
}

edge_sets <- lapply(networks, get_edge_set)

# -------------------------------
# Compute Jaccard similarity matrix
# -------------------------------
n <- length(edge_sets)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(edge_sets)

for (i in 1:n) {
  for (j in i:n) {
    a <- edge_sets[[i]]
    b <- edge_sets[[j]]
    jaccard <- if (length(union(a, b)) > 0) {
      length(intersect(a, b)) / length(union(a, b))
    } else {
      NA
    }
    jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- jaccard
  }
}

# -------------------------------
# Save matrix and heatmap
# -------------------------------
csv_out <- file.path(opt$out_dir, "jaccard_edge_matrix.csv")
pdf_out <- file.path(opt$out_dir, "jaccard_edge_heatmap.pdf")

fwrite(as.data.table(jaccard_matrix, keep.rownames = TRUE), csv_out)
pheatmap(jaccard_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, main = "Edge Jaccard Similarity",
         filename = pdf_out)

message("✅ Jaccard matrix saved to: ", csv_out)
message("✅ Heatmap saved to: ", pdf_out)

