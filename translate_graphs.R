#!/usr/bin/env Rscript

# ================================================================
# Load libraries
# ================================================================
pacman::p_load(
  optparse,
  igraph,
  dplyr,
  data.table,
  tools,
  org.Hs.eg.db,
  AnnotationDbi
)

# ================================================================
# Convert ENSG → SYMBOL using org.Hs.eg.db (offline)
# ================================================================
convert_ens_to_symbol <- function(ensembl_ids) {
  
  # Remove version suffixes .1, .12, etc.
  clean_ids <- gsub("\\..*$", "", ensembl_ids)
  
  # Convert ENSG → SYMBOL
  symbols <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db,
    keys = clean_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  # Replace missing symbols with original ENSG
  symbols[is.na(symbols)] <- clean_ids[is.na(symbols)]
  
  return(symbols)
}

# ================================================================
# Translate vertex names
# ================================================================
translate_vertex_names <- function(graph) {
  
  original_names <- V(graph)$name
  translated <- convert_ens_to_symbol(original_names)
  
  V(graph)$name_trad <- translated
  
  return(graph)
}

# ================================================================
# read_network() — supports graphml, csv, tsv, txt
# ================================================================
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

# ================================================================
# PARSER — run multiple networks from command line
# ================================================================
option_list <- list(
  
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directorio con redes [default: directorio actual]", metavar = "character"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = "./translated_networks",
              help = "Directorio de salida [default: %default]", metavar = "character"),
  
  make_option(c("-p", "--pattern"), type = "character",
              default = "\\.(csv|tsv|txt|graphml)$",
              help = "Regex para seleccionar redes [default: %default]", metavar = "character"),
  
  make_option(c("-t", "--type"), type = "character", default = "auto",
              help = "Tipo de red: auto, edgelist, adjacency [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

opt$input_dir <- "/STORAGE/csbig/networks_final/fomo_networks"
opt$output_dir <- "/STORAGE/csbig/networks_final/fomo_networks/translated_graphs"
opt$pattern <- "\\.tsv$"
opt$type <- "edgelist"

# Output folder
if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# Find network files
files <- list.files(opt$input_dir,
                    pattern = opt$pattern,
                    full.names = TRUE)

message("Archivos encontrados:")
print(files)

if (length(files) == 0)
  stop("❌ No se encontraron redes con el patrón especificado.")

#PROCESSING MULTIPLE NETWORKS
for (f in files) {
  
  # Read
  g <- read_network(f, type = opt$type)
  
  # Translate
  g <- translate_vertex_names(g)
  
  # Output name
  out_name <- file.path(
    opt$output_dir,
    paste0(file_path_sans_ext(basename(f)), "_translated.graphml")
  )
  
  write_graph(g, out_name, format = "graphml")
  
  message("✔️ Guardado: ", out_name)
}

message("Translating finished")