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

#Helper function: read a network from CSV/TSV/GraphML
read_network <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "graphml") {
    #GraphML file
    g <- read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- fread(path, sep = sep, header = FALSE, data.table = FALSE)
    
    # Check if first row contains non-numeric column names (adjacency matrix)
    if (!all(sapply(df[1, ], is.numeric))) {
      #First row = column names
      colnames(df) <- as.character(df[1, ])
      df <- df[-1, , drop = FALSE]
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      rownames(mat) <- colnames(mat)
      #Create weighted undirected graph
      g <- graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
    } else if (ncol(df) >= 2) {
      #Treat as edgelist: 2 or 3 columns
      if (ncol(df) >= 3) {
        g <- graph_from_data_frame(df[, 1:3], directed = FALSE)
      } else {
        g <- graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
    } else {
      stop("Cannot interpret the file: ", path)
    }
    
  } else {
    stop("Unsupported file format: ", ext)
  }
  
  #Simplify graph (remove loops and multiple edges)
  g <- simplify(igraph::as_undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

#Helper function for Jaccard index calculation

# Jaccard simple (conjuntos a y b)
jaccard_simplex <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)   # ambos vacíos -> similitud máxima
  if (length(a) == 0 || length(b) == 0) return(0)   # uno vacío -> 0
  length(intersect(a, b)) / length(union(a, b))
}

#Build sets
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos 2 archivos para comparar.")
net_names <- basename(files)

############################

# Serializar aristas no dirigidas a pares canonizados "u||v"
edge_set <- function(g) {
  e <- igraph::as_edgelist(g, names = TRUE)
  if (nrow(e) == 0) return(character(0))
  u <- ifelse(e[,1] <= e[,2], e[,1], e[,2])
  v <- ifelse(e[,1] <= e[,2], e[,2], e[,1])
  paste0(u, "||", v)
}

# Conjunto de nodos (usa nombres si existen, si no, IDs)
node_set <- function(g) {
  nm <- igraph::V(g)$name
  if (is.null(nm)) nm <- as.character(seq_len(igraph::vcount(g)))
  unique(nm)
}

# ========= CARGAR ARCHIVOS Y EXTRAER CONJUNTOS =========
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("Se requieren al menos 2 redes para comparar.")
labels <- basename(files)

message("Leyendo ", length(files), " redes…")
# Leer y convertir a conjuntos en serie para controlar RAM
sets <- vector("list", length(files))
for (i in seq_along(files)) {
  g <- read_network(files[i])
  sets[[i]] <- if (tolower(opt$element) == "nodes") node_set(g) else edge_set(g)
  rm(g); gc(verbose = FALSE)
}
names(sets) <- labels

# ========= MATRIZ DE JACCARD EN PARALELO =========
future::plan(future::multisession, workers = opt$workers)

n <- length(sets)
# computamos sólo triángulo superior en paralelo (por filas)
rows_dt <- future.apply::future_lapply(seq_len(n), function(i) {
  si <- sets[[i]]
  sims <- vector("list", n - i + 1)
  idx <- 1
  for (j in i:n) {
    sj <- sets[[j]]
    sims[[idx]] <- data.table::data.table(i = i, j = j, sim = jaccard_simplex(si, sj))
    idx <- idx + 1
  }
  data.table::rbindlist(sims)
}, future.seed = TRUE)

pairs_dt <- data.table::rbindlist(rows_dt)
# construir matriz simétrica
sim_mat <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
sim_mat[cbind(pairs_dt$i, pairs_dt$j)] <- pairs_dt$sim
sim_mat[cbind(pairs_dt$j, pairs_dt$i)] <- pairs_dt$sim
diag(sim_mat) <- 1

# distancia = 1 - similitud
dist_mat <- 1 - sim_mat

# ========= ORDENACIÓN POR CLUSTERING (para el heatmap) =========
# Si hay NAs residuales (no debería), los ponemos en 1 (máx distancia)
dist_mat[is.na(dist_mat)] <- 1
# dist() requiere objeto "dist"; nos aseguramos de simetría y diagonal cero
d <- as.dist(dist_mat)
hc <- hclust(d, method = "average")
ord <- hc$order
sim_mat_ord <- sim_mat[ord, ord, drop = FALSE]
dist_mat_ord <- dist_mat[ord, ord, drop = FALSE]
labels_ord <- rownames(sim_mat_ord)

# ========= SALIDAS =========
# 1) CSVs
data.table::fwrite(as.data.table(sim_mat, keep.rownames = "network"), file.path(opt$out_dir, "jaccard_similarity_matrix.csv"))
data.table::fwrite(as.data.table(dist_mat, keep.rownames = "network"), file.path(opt$out_dir, "jaccard_distance_matrix.csv"))

# tabla larga de pares
long_dt <- data.table::as.data.table(as.table(sim_mat))
data.table::setnames(long_dt, c("Var1","Var2","N"), c("net_i","net_j","jaccard"))
long_dt <- long_dt[net_i < net_j][order(-jaccard)]
data.table::fwrite(long_dt, file.path(opt$out_dir, "jaccard_pairs_long.csv"))

# 2) HEATMAP (PNG con ggplot2)
heat_dt <- data.table::as.data.table(as.table(sim_mat_ord))
data.table::setnames(heat_dt, c("Var1","Var2","N"), c("net_i","net_j","jaccard"))
heat_dt[, net_i := factor(net_i, levels = labels_ord)]
heat_dt[, net_j := factor(net_j, levels = labels_ord)]

p <- ggplot2::ggplot(heat_dt, ggplot2::aes(net_i, net_j, fill = jaccard)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "Jaccard") +
  ggplot2::coord_fixed() +
  ggplot2::labs(title = sprintf("Matriz de similitud Jaccard (%s)", tolower(opt$element)),
                x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                 panel.grid = ggplot2::element_blank())

ggplot2::ggsave(filename = file.path(opt$out_dir, "jaccard_heatmap.png"),
                plot = p, width = max(6, min(16, n * 0.25)), height = max(6, min(16, n * 0.25)), dpi = 200)

# 3) REPORTE HTML opcional
if (isTRUE(opt$make_html)) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(rmarkdown))
  rmd <- file.path(opt$out_dir, "jaccard_report.Rmd")
  cat(
    "---
title: \"Comparación de redes por índice de Jaccard\"
output: html_document
params:
  out_dir: \"", opt$out_dir, "\"
  element: \"", tolower(opt$element), "\"
---

```{r setup, echo=FALSE, message=FALSE}
library(data.table); library(ggplot2)
out_dir <- params$out_dir
element <- params$element
sim <- fread(file.path(out_dir, 'jaccard_similarity_matrix.csv'))
dist <- fread(file.path(out_dir, 'jaccard_distance_matrix.csv'))
pairs <- fread(file.path(out_dir, 'jaccard_pairs_long.csv'))













