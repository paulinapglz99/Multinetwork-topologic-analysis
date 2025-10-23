#!/usr/bin/env Rscript

#Which modules (gene communities)
#are maintained across AD networks and are not present in control networks? And vice versa.

#AD core modules and control core modules 

#Load packages
pacman::p_load(data.table, tidyverse, optparse, ggplot2, tools, reshape2)

#Parser
option_list <- list(
  make_option(c("-m", "--metadata"), type = "character", help = "CSV file with metadata: file_name, region, condition, file_path", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard", help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#Inputs
#opt$metadata <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/networks_nodes_index.csv"
metadata_file <- opt$metadata
#opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules"
out_dir <- opt$output_dir

message("Reading metadata from: ", metadata_file)
message("Output to: ", out_dir)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#Functions

jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_matrix <- function(redA, redB) {
  mat <- matrix(0, nrow = length(redA), ncol = length(redB),
                dimnames = list(names(redA), names(redB)))
  for (i in seq_along(redA)) {
    for (j in seq_along(redB)) {
      mat[i, j] <- jaccard_index(redA[[i]], redB[[j]])
    }
  }
  return(mat)
}

compare_networks <- function(list_nets, pairs) {
  results <- list()
  for (p in pairs) {
    netA_name <- p[1]
    netB_name <- p[2]
    netA <- list_nets[[netA_name]]
    netB <- list_nets[[netB_name]]
    name <- paste(netA_name, netB_name, sep = "_vs_")
    message("Comparing ", name)
    results[[name]] <- jaccard_matrix(netA, netB)
  }
  return(results)
}

#Read metadata and modules

metadata <- vroom(metadata_file)
stopifnot(all(c("file_name", "region", "condition", "file_path") %in% colnames(metadata)))

modules_nets <- list()
for (i in 1:nrow(metadata)) {
  net_name <- metadata$file_name[i]
  df <- fread(metadata$file_path[i])
  stopifnot(all(c("node", "membership_infomap") %in% colnames(df)))
  modules_nets[[net_name]] <- split(df$node, df$membership_infomap)
}

# -------------------------
# Generate pairs
# -------------------------

# 1. AD vs AD
ad_pairs <- metadata[condition == "AD"]
ad_vs_ad <- unlist(lapply(split(ad_pairs, ad_pairs$region), function(df) {
  combn(df$file_name, 2, simplify = FALSE)
}), recursive = FALSE)

# 2. Control vs Control
ctrl_pairs <- metadata[condition == "control"]
ctrl_vs_ctrl <- unlist(lapply(split(ctrl_pairs, ctrl_pairs$region), function(df) {
  combn(df$file_name, 2, simplify = FALSE)
}), recursive = FALSE)

# 3. AD vs Control (one-to-one per region)
ad_vs_ctrl <- list()
for (reg in unique(metadata$region)) {
  ad_file <- metadata[region == reg & condition == "AD", file_name]
  ctrl_file <- metadata[region == reg & condition == "control", file_name]
  ad_vs_ctrl[[reg]] <- lapply(seq_along(ad_file), function(i) c(ad_file[i], ctrl_file[i]))
}
ad_vs_ctrl <- unlist(ad_vs_ctrl, recursive = FALSE)

# Combine all comparisons
all_pairs <- c(ad_vs_ad, ctrl_vs_ctrl, ad_vs_ctrl)

# -------------------------
# Run comparisons and export
# -------------------------

results <- compare_networks(modules_nets, all_pairs)

for (name in names(results)) {
  mat <- results[[name]]
  out_csv <- file.path(out_dir, paste0(name, "_jaccard_matrix.csv"))
  fwrite(as.data.frame(as.table(mat)), out_csv)
  
  # Optional heatmap
  p <- ggplot(melt(as.data.frame(as.table(mat))), aes(Var1, Var2, fill = Freq)) +
    geom_tile() + scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = paste("Jaccard Matrix:", name), x = "Modules A", y = "Modules B")
  ggsave(file.path(out_dir, paste0(name, "_heatmap.png")), plot = p, width = 8, height = 6)
}

message("✅ Comparisons complete. Results saved to: ", out_dir)
