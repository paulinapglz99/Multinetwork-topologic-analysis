#!/usr/bin/env Rscript

#Script: module_expression_summary.R
#Goal: For each network, compute per-module average expression and eigengene
#Usage: Rscript module_expression_summary.R --input_dir ./expression_matrices --module_topos ./module_info --out_dir ./module_outputs

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")

pacman::p_load("tidyverse", "vroom", "optparse", "tools")

#Define CLI options
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", help="Directory with expression matrices"),
  make_option(c("-m", "--module_topos"), type="character", help="Directory with *_nodes_summary.csv files"),
  make_option(c("-o", "--out_dir"), type="character", default="module_outputs", help="Output directory"))

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

expr_files <- list.files(opt$input_dir, pattern = "\\.tsv$", full.names = TRUE)
module_files <- list.files(opt$module_topos, pattern = "_nodes_summary\\.csv$", full.names = TRUE)

if (length(expr_files) == 0 || length(module_files) == 0) {
  stop("No input files found. Check input directory and patterns.")
}

#Functions 

process_network <- function(expr_file, module_file) {
  message("Processing: ", basename(expr_file))
  
  expr <- vroom(expr_file, delim = "\t", col_types = cols()) %>%
    column_to_rownames(var = 1) %>%
    as.matrix()
  
  if (ncol(expr) < nrow(expr)) {
    expr <- t(expr)
  }
  
  modules <- vroom(module_file, delim = ",", col_types = cols()) %>%
    rename(gene = 1) %>%
    filter(gene %in% rownames(expr))
  
  profiles <- modules %>%
    group_by(membership) %>%
    summarise(
      genes = list(gene),
      .groups = "drop"
    ) %>%
    mutate(
      data = map(genes, ~ expr[.x, , drop = FALSE]),
      mean_profile = map(data, ~ colMeans(.x, na.rm = TRUE)),
      eigengene = map(data, function(mat) {
        if (nrow(mat) < 2) return(rep(NA, ncol(mat)))
        pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
        as.numeric(pca$x[, 1])
      }),
      sample = list(colnames(expr))
    ) %>%
    select(membership, sample, mean_profile, eigengene) %>%
    unnest(cols = c(sample, mean_profile, eigengene)) %>%
    rename(module = membership, mean_expr = mean_profile, eigengene = eigengene)
  
  profiles$network <- file_path_sans_ext(basename(expr_file))
  return(profiles)
}

#Loop through matched files
results <- purrr::map_dfr(module_files, function(mod_file) {
  base <- str_replace(basename(mod_file), "_topN200000_nodes_summary.csv", "")
  expr_file <- expr_files[str_detect(expr_files, fixed(base))]
  
  if (length(expr_file) == 1) {
    process_network(expr_file, mod_file)
  } else {
    warning("Could not match expression file for: ", mod_file)
    return(NULL)
  }
})

# Save results
output_file <- file.path(opt$out_dir, "module_expression_profiles.tsv")
vroom::vroom_write(results, output_file, delim = "\t")
message("✅ Results written to: ", output_file)
