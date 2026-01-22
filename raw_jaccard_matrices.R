#!/usr/bin/env Rscript

#core_modules.R

#Packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  data.table, 
  tidyverse,
  tools,
  optparse,
  ggplot2,
  igraph
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directorio con archivos *_nodes_summary.csv", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard",
              help = "Directorio de salida [default: %default]", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex para seleccionar archivos [default: %default]", metavar = "character"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#Inputs
#opt$input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/"
#opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Functions

#A regex to get brain area and phenotype
extract_info <- function(filename) {
  base <- basename(filename)
  if (grepl("_counts_", base)) {
    region <- sub("_counts_.*", "", base)
    phenotype <- ifelse(grepl("_AD_", base), "AD", "Control")
    return(data.frame(filename = filename, region = region, phenotype = phenotype))
  } else {
    return(NULL)
  }
}

jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_matrix <- function(modA, modB) {
  mat <- matrix(0, nrow = length(modA), ncol = length(modB),
                dimnames = list(names(modA), names(modB)))
  for (i in seq_along(modA)) {
    for (j in seq_along(modB)) {
      mat[i, j] <- jaccard_index(modA[[i]], modB[[j]])
    }
  }
  return(mat)
}

load_modules <- function(file) {
  df <- fread(file)
  if (!all(c("node", "membership_infomap") %in% names(df))) {
    stop("File ", file, " must contain columns 'node' and 'membership_infomap'")
  }
  split(df$node, df$membership_infomap)
}

#Build pairs
build_pairs <- function(meta) {
  ad_files <- meta %>% filter(phenotype == "AD")
  control_files <- meta %>% filter(phenotype == "Control")
  
  pairs_ad_ad <- combn(ad_files$filename, 2, simplify = FALSE)
  pairs_control_control <- combn(control_files$filename, 2, simplify = FALSE)
  
  common_regions <- intersect(ad_files$region, control_files$region)
  pairs_ad_control <- lapply(common_regions, function(region) {
    ad_file <- ad_files %>% filter(region == !!region) %>% pull(filename)
    control_file <- control_files %>% filter(region == !!region) %>% pull(filename)
    list(ad_file, control_file)
  })
  
  return(list(
    AD_vs_AD = pairs_ad_ad,
    Control_vs_Control = pairs_control_control,
    AD_vs_Control = pairs_ad_control
  ))
}

#Compare networks
compare_and_export <- function(pairs, tag) {
  for (p in pairs) {
    fileA <- p[[1]]
    fileB <- p[[2]]
    nameA <- tools::file_path_sans_ext(basename(fileA))
    nameB <- tools::file_path_sans_ext(basename(fileB))
    name_out <- paste(tag, nameA, "vs", nameB, sep = "_")
    message("Comparando: ", name_out)
    
    netA <- load_modules(fileA)
    netB <- load_modules(fileB)
    
    jacmat <- jaccard_matrix(netA, netB)
    
    write.csv(jacmat, file = file.path(output_dir, paste0(name_out, "_jaccard_matrix.csv")))
  }
}

#Read files --- ---

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos dos archivos *_nodes_summary.csv.")

#Get metadata --- ---
meta <- do.call(rbind, lapply(files, extract_info))

print(meta)

#Get pairs!
all_pairs <- build_pairs(meta)
compare_and_export(all_pairs$AD_vs_AD, "ADvsAD")
compare_and_export(all_pairs$Control_vs_Control, "ControlvsControl")
compare_and_export(all_pairs$AD_vs_Control, "ADvsControl")

message("Comparisons DONE!. Files are now in: ", output_dir)