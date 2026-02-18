#!/usr/bin/env Rscript

# community_enrichment.R

############################
# Package management
############################

if (!requireNamespace("pacman", quietly = FALSE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

if (!requireNamespace("optparse", quietly = FALSE))
  install.packages("optparse", repos = "https://cloud.r-project.org")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("clusterProfiler", quietly = FALSE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = FALSE))
  BiocManager::install("org.Hs.eg.db")

ok <- pacman::p_load(
  igraph,
  data.table,
  future.apply,
  ggplot2,
  tidyverse,
  jsonlite,
  stringr,
  optparse,
  tools,
  clusterProfiler,
  org.Hs.eg.db
)

if (!all(ok)) {
  stop("Some packages failed to load: ",
       paste(names(ok)[!ok], collapse = ", "))
}

message("All packages loaded correctly.")

############################
# Command line options
############################

option_list <- list(
  optparse::make_option(
    c("-i","--index_file"),
    type="character",
    help="Index file (csv/tsv with file,node_table)"
  ),
  optparse::make_option(
    c("-o","--out_dir"),
    type="character",
    default="results_enrichment",
    help="Output directory"
  ),
  optparse::make_option(
    c("-w","--workers"),
    type="integer",
    default=4,
    help="Number of parallel workers"
  ),
  optparse::make_option(
    c("--seed"),
    type="integer",
    default=42,
    help="Seed"
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list)
)

if (is.null(opt$index_file))
  stop("--index_file missing")

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(opt$seed)

############################
# Ontologies to evaluate
############################

ontologies <- c("BP", "MF")

############################
# Enrichment function
############################

enricher_fun <- function(nodes, network_universe, ontology) {
  enrichGO(
    gene = nodes,
    OrgDb = org.Hs.eg.db,
    universe = network_universe,
    keyType = "ENSEMBL",
    readable = TRUE,
    ont = ontology,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10
  )
}

############################
# Network processing function
############################

process_network <- function(row) {
  
  fname <- row[["file"]]
  summary_path <- row[["node_table"]]
  
  message("Processing: ", fname)
  
  out_file <- file.path(
    opt$out_dir,
    paste0(tools::file_path_sans_ext(fname), "_enrichment.csv")
  )
  
  summary_df <- fread(summary_path)
  
  if (nrow(summary_df) == 0) {
    message("Skipping ", fname, " (empty node table)")
    fwrite(data.table(message="Empty input"), out_file)
    return(out_file)
  }
  
  if (!"membership" %in% colnames(summary_df))
    stop('"node_table" must contain column: membership')
  
  network_universe <- unique(summary_df$node)
  nodes_by_comm <- split(summary_df$node, summary_df$membership)
  
  enriched_list <- lapply(names(nodes_by_comm), function(community_id) {
    
    nodes <- nodes_by_comm[[community_id]]
    
    # robustness: skip very small communities
    if (length(nodes) < 5) return(NULL)
    
    ont_list <- lapply(ontologies, function(ont) {
      
      res <- tryCatch(
        enricher_fun(nodes, network_universe, ont),
        error = function(e) NULL
      )
      
      if (is.null(res)) return(NULL)
      
      df <- as.data.frame(res)
      
      if (nrow(df) == 0) return(NULL)
      
      df$CommunityID <- community_id
      df$Ontology <- ont
      
      return(df)
    })
    
    data.table::rbindlist(ont_list, fill = TRUE)
  })
  
  final_summary <- data.table::rbindlist(enriched_list, fill = TRUE)
  
  if (nrow(final_summary) > 0) {
    fwrite(final_summary, out_file)
  } else {
    message("No significant enrichment found for: ", fname)
    fwrite(data.table(message="No significant enrichment"), out_file)
  }
  
  return(out_file)
}

############################
# Parallel execution
############################

index_df <- fread(opt$index_file)

plan(multisession, workers = opt$workers)

results <- future_lapply(
  split(index_df, seq_len(nrow(index_df))),
  function(row) process_network(row),
  future.seed = TRUE
)

message("Success. Results stored in: ", opt$out_dir)

#END