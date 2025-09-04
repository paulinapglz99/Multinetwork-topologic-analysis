#!/usr/bin/env Rscript

#compare_networks.R

if (!requireNamespace("pacman", quietly = F)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = F)) install.packages("optparse", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "igraph",
  "data.table",
  "future.apply",
  "ggplot2",
  "tidyverse",
  "jsonlite",
  "stringr",
  "optparse", 
  "tools"
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages loaded correctly.: ",
       paste(names(ok)[!ok], collapse = ", "))
}

#Define option list for inputs

option_list <- list(
  optparse::make_option(c("-i","--input_dir"), type="character", help="Carpeta con redes"),
  optparse::make_option(c("-p","--pattern"), type="character", default=".*\\.(csv|tsv|txt|graphml)$", help="File regex (default: .*\\.(csv|tsv|txt|graphml)$)"),
  optparse::make_option(c("-o","--out_dir"), type="character", default="results_jaccard", help="Output directory"),
  optparse::make_option(c("-w","--workers"), type="integer", default=4, help="Number of parallel workers"),
  optparse::make_option(c("-e","--element"), type="character", default="edges", help="Element to compare: “edges” or “nodes” (default: edges)"),
  optparse::make_option(c("--make_html"), action="store_true", default=FALSE, help="Create HTML report"),
  optparse::make_option(c("--seed"), type="integer", default=42, help="Seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (is.null(opt$input_dir)) stop("Please add --input_dir")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

#Function to read graphs
read_network <- function(path, type = opt$type) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    if (type == "edgelist") {
      if (ncol(df) >= 3) {
        g <- igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
      } else {
        g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
    } else if (type == "adjacency") {
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
    } else if (type == "auto") {
      looks_like_adj <- all(sapply(df[,-1, drop=FALSE], is.numeric))
      if (looks_like_adj) {
        rownames(df) <- make.unique(as.character(df[, 1]))  
        df <- df[, -1, drop = FALSE]
        mat <- as.matrix(df)
        storage.mode(mat) <- "numeric"
        g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      } else {
        if (ncol(df) >= 3) {
          g <- igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
        } else {
          g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
        }
      }
    }
  } else {
    stop("Unsupported file format: ", ext)
  }
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"),
                        remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

#Helper function for Jaccard index calculation

jaccard_simplex <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)   #both empty-> max sim
  if (length(a) == 0 || length(b) == 0) return(0)   #one is empty-> 0
  length(intersect(a, b)) / length(union(a, b))
}

#Build sets
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos 2 archivos para comparar.")
net_names <- basename(files)

#Compare modularity between graphs --- ---
#, applying 
#variation of information "vi"
#normalized mutual information "nmi"
#split-join distance "split-join distance"
#Rand index "Rand index"
#adjusted Rand index "adjusted Rand index"

possible_algos <- c("vi", "nmi", "split.join", "rand", "adjusted.rand")

comparison_methods <- sapply(X = possible_algos, FUN = function(i){
  igraph::compare(comm1 = graphAD_plus_modules,
                  comm2 = graphnoAD_plus_modules,
                  method = i
  )
})

comparison_methods








