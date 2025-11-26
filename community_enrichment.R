#!/usr/bin/env Rscript

#community_enrichment.R

if (!requireNamespace("pacman", quietly = F)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = F)) install.packages("optparse", repos = "https://cloud.r-project.org")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = F)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = F)) BiocManager::install("org.Hs.eg.db")

ok <- pacman::p_load(
  "igraph",
  "data.table",
  "future.apply",
  "ggplot2",
  "tidyverse",
  "jsonlite",
  "stringr",
  "optparse", 
  "tools",
  "clusterProfiler",
  "org.Hs.eg.db"
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages loaded correctly.: ",
       paste(names(ok)[!ok], collapse = ", "))
}

#Define option list for inputs

option_list <- list(
  optparse::make_option(c("-i","--index_file"), type="character", help="Index file (csv/tsv with file,node_table, output of network_topology.R)"),
  #optparse::make_option(c("-u","--universe"), type="character", help="File containing gene universe (txt)"),
  optparse::make_option(c("-o","--out_dir"), type="character", default="results_enrichment", help="Output directory"),
  optparse::make_option(c("-w","--workers"), type="integer", default=4, help="Number of parallel workers"),
  optparse::make_option(c("--seed"), type="integer", default=42, help="Seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (is.null(opt$index_file)) stop("--index_file missing")

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
#universe <- scan(opt$universe, what = character())

#Function to read graphs

enricher_fun <- function(nodes, network_universe) {
  enrichGO(gene = nodes,
           OrgDb = org.Hs.eg.db,
           universe = network_universe,
           keyType = 'ENSEMBL',
           readable = TRUE,
           ont = "BP",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.10)
}

#Function to replace nulls in enrichment lists
replace_null <- function(x) if (is.null(x)) new("enrichResult") else x

#Function to process network
process_network <- function(row) {
  # row: c(file, node_table)
  fname <- row[["file"]]
  summary_path <- row[["node_table"]]
  
  message("Processing: ", fname)
  
  #Define output path
  out_file <- file.path(opt$out_dir, paste0(tools::file_path_sans_ext(fname), "_enrichment.csv"))
  
  #Read the path summary
  summary_df <- fread(summary_path)
  
  #Add a check to ensure that the file is not empty.
  if (nrow(summary_df) == 0) {
    message("Warning: Skipping ", fname, " because its node table is empty.")
    return(paste0(out_file, " (skipped, empty input)"))
  }
  
  if (!"membership" %in% colnames(summary_df)) {
    stop('"node_table" does not have the membership column"')
  }
  
  #Read network universe
  network_universe <- unique(summary_df$node)
  
  #Split the nodes per community to enrich
  nodes_by_comm <- split(summary_df$node, summary_df$membership)
  
  #Enrich each comm
  enriched_list <- lapply(names(nodes_by_comm), function(community_id) {
    nodes <- nodes_by_comm[[community_id]]
    res <- tryCatch(enricher_fun(nodes, network_universe), error = function(e) NULL)
    res <- replace_null(res)
    
#Convert to dataframes
    if (nrow(as.data.frame(res)) > 0) {
      df <- as.data.frame(res)
      df$CommunityID <- community_id
      return(df)
    } else {
      return(NULL)
    }
  })
  
#Combine dataframes:)
final_summary <- data.table::rbindlist(enriched_list)
  
#Save combined result (if not empty)
  if (nrow(final_summary) > 0) {
    fwrite(final_summary, file = out_file)
  } else {
    message("No significant enrichment found for any community in: ", fname)
    #Write an empty file to register the processed network even if it had no sig enrichment
    fwrite(data.table(message = "No significant enrichment"), file = out_file)
  }
  
  return(out_file)
}

#Run in parallel

index_df <- fread(opt$index_file)
plan(multisession, workers = opt$workers)

#Run results
results <- future_lapply(split(index_df, seq_len(nrow(index_df))), function(row) {
  process_network(row)
}, future.seed = TRUE)

message("Success. Results in : ", opt$out_dir)

#END