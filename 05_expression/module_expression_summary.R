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

# expr_files <- list.files("~/Desktop/local_work/matrices_originales",
#                          pattern = "\\.tsv$", full.names = TRUE)
# module_files <- list.files("~/Desktop/local_work/fomo_networks/results_topos_louvain",
#                            pattern = "_nodes_summary\\.csv$", full.names = TRUE)

if (length(expr_files) == 0 || length(module_files) == 0) {
  stop("No input files found. Check input directory and patterns.")
}

#Load expression matrix
load_expression_matrix <- function(expr_file) {
  vroom(expr_file, delim = "\t", col_types = cols()) %>%
    column_to_rownames(var = "gene") %>%
    as.matrix() %>%
    { if (ncol(.) < nrow(.)) t(.) else . }
}

#Load module annotations
load_module_annotations <- function(module_file, expr_genes) {
  vroom(module_file, delim = ",", col_types = cols()) %>%
    rename(gene = 1) %>%
    filter(gene %in% expr_genes)
}

#Extract module profiles
extract_module_profiles <- function(expr_mat, module_df) {
  module_df %>%
    group_by(membership) %>%
    summarise(genes = list(gene), .groups = "drop") %>%
    mutate(
      expr_data = map(genes, ~ expr_mat[.x, , drop = FALSE]),
      mean_profile = map(expr_data, ~ colMeans(.x, na.rm = TRUE)),
      eigengene = map(expr_data, function(mat) {
        if (nrow(mat) < 2) return(rep(NA, ncol(mat)))
        prcomp(t(mat), center = TRUE, scale. = TRUE)$x[, 1]
      }),
      sample = list(colnames(expr_mat))
    ) %>%
    select(membership, sample, mean_profile, eigengene) %>%
    unnest(cols = c(sample, mean_profile, eigengene)) %>%
    rename(module = membership, mean_expr = mean_profile)
}

#Take the genes per module and calculate the average expression profile and the eigengene (first PCA component)
extract_module_profiles <- function(expr_mat, module_df) {
  module_df %>%
    group_by(membership) %>%
    summarise(genes = list(gene), .groups = "drop") %>%
    mutate(
      expr_data = map(genes, ~ expr_mat[.x, , drop = FALSE]),
      mean_profile = map(expr_data, ~ colMeans(.x, na.rm = TRUE)),
      eigengene = map(expr_data, function(mat) {
        if (nrow(mat) < 2) return(rep(NA, ncol(mat)))
        prcomp(t(mat), center = TRUE, scale. = TRUE)$x[, 1]
      }),
      sample = list(colnames(expr_mat))
    ) %>%
    select(membership, sample, mean_profile, eigengene) %>%
    unnest(cols = c(sample, mean_profile, eigengene)) %>%
    rename(module = membership, mean_expr = mean_profile)
}

#It takes a pair (expression file and module file), loads and processes them, returning a data.frame.

process_network_pair <- function(expr_file, module_file) {
  message("Processing: ", basename(expr_file))
  
  expr_mat <- load_expression_matrix(expr_file)
  modules <- load_module_annotations(module_file, rownames(expr_mat))
  
  profiles <- extract_module_profiles(expr_mat, modules)
  profiles$network <- file_path_sans_ext(basename(expr_file))
  profiles
}

#It matches names and applies process_network_pair() to all of them
match_files_and_process_all <- function(expr_files, module_files) {
  purrr::map_dfr(module_files, function(mod_file) {
    #base <- str_replace(basename(mod_file), "_topN200000_nodes_summary.csv", "")
    base <- str_replace(basename(mod_file), "_top.*_nodes_summary\\.csv$", "")
    expr_file <- expr_files[str_detect(expr_files, fixed(base))]
    
    if (length(expr_file) == 1) {
      process_network_pair(expr_file, mod_file)
    } else {
      warning("Could not match expression file for: ", mod_file)
      return(NULL)
    }
  })
}

#

results <- match_files_and_process_all(expr_files, module_files)

vroom::vroom_write(results, file.path(opt$out_dir, "module_expression_profiles.tsv"), delim = "\t")

# Save results
output_file <- file.path(opt$out_dir, "module_expression_profiles.tsv")
vroom::vroom_write(results, output_file, delim = "\t")
message("Results written to: ", output_file)



#Probando  las funciones

#las cosas

expr_file <- expr_files[7]
module_file <- module_files[7]

#Las cosas

load_expression_matrix <- function(expr_file) {
  vroom(expr_file) %>%
    column_to_rownames(var = "gene") %>%
    as.matrix()
}

expr_mat <- load_expression_matrix(expr_file)

#Load module annotations
load_module_annotations <- function(module_file, expr_genes) {
  vroom(module_file, delim = ",", col_types = cols()) %>%
    rename(gene = 1) %>%
    filter(gene %in% expr_genes)
}

modules <- load_module_annotations(module_file, rownames(expr_mat))

#Take the genes per module and calculate the average expression profile and the eigengene (first PCA component)
extract_module_profiles <- function(expr_mat, modules) {
  modules %>%
    group_by(membership) %>%
    summarise(genes = list(gene), .groups = "drop") %>%
    mutate(
      expr_data = map(genes, ~ expr_mat[.x, , drop = FALSE]),
      mean_profile = map(expr_data, ~ colMeans(.x, na.rm = TRUE)),
      eigengene = map(expr_data, function(mat) {
        if (nrow(mat) < 2) return(rep(NA, ncol(mat)))
        prcomp(t(mat), center = TRUE, scale. = TRUE)$x[, 1]
      }),
      sample = list(colnames(expr_mat))
    ) %>%
    select(membership, sample, mean_profile, eigengene) %>%
    unnest(cols = c(sample, mean_profile, eigengene)) %>%
    rename(module = membership, mean_expr = mean_profile)
}

profiles <- extract_module_profiles(expr_mat, modules)
