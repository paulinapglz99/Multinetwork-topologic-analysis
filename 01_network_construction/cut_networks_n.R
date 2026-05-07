#!/usr/bin/env Rscript
# Script to filter edgelists per top-N highest absolute MI values

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "data.table",
  "future.apply",
  "optparse",
  "tools"
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages failed to load: ",
       paste(names(ok)[!ok], collapse = ", "))
}

# Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL, help = "Folder containing edgelist files", metavar = "path"),
  make_option(c("-o", "--out_dir"), type = "character", default = "results", help = "Output directory [default: %default]", metavar = "path"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$", help = "File name pattern to match [default: %default]", metavar = "regex"),
  make_option(c("-n", "--top_n"), type = "integer", default = 100000, help = "Number of top edges by absolute MI to keep [default: %default]", metavar = "int"),
  make_option(c("-w", "--workers"), type = "integer", default = 2, help = "Number of parallel workers [default: %default]", metavar = "int"),
  make_option(c("--seed"), type = "integer", default = 42, help = "Random seed for reproducibility [default: %default]", metavar = "int")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Specify --input_dir")

if (!dir.exists(opt$input_dir)) stop("Input directory not found: ", opt$input_dir)

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

#Read files
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) == 0)
  stop("No files matching pattern '", opt$pattern, "' found in ", opt$input_dir)

message("Found ", length(files), " files to process.")

#Function to filter top N edges by absolute MI
filter_top_n_mi <- function(file_path, top_n = 100000, output_dir) {
  message("Processing: ", basename(file_path))
  
  df <- fread(file_path, col.names = c("gene1", "gene2", "MI"))
  
  df <- df[order(-abs(MI))]
  df_top <- head(df, top_n)
  
  output_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(file_path)), "_topN", top_n, ".tsv")
  )
  
  fwrite(df_top, output_file, sep = "\t")
  message("Saved: ", basename(output_file))
  
  return(output_file)
}

#Parallelize over files
plan(multisession, workers = opt$workers)

results <- future_lapply(
  files,
  filter_top_n_mi,
  top_n = opt$top_n,
  output_dir = opt$out_dir
)

message("\nFiltered networks saved in: ", opt$out_dir)

#END