#!/usr/bin/env Rscript
# =============================================================================
# plots_global.R
# Generate global network topology plots: PCA and metric differences heatmap.
# Also computes the AD vs control Jaccard edge overlap per region (console).
# Requires: networks_summary.csv and raw network .tsv files
# =============================================================================
set.seed(10)

# -----------------------------------------------------------------------------
# 0. Dependencies
# -----------------------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

pacman::p_load(
  igraph, data.table,
  ggplot2, tidyverse,
  tools, purrr,
  ggpubr, pheatmap,
  vroom, ggrepel, cowplot,
  ggdendro, patchwork
)

# -----------------------------------------------------------------------------
# 1. Paths (repo-relative; run from the repository root)
# -----------------------------------------------------------------------------

source("config/paths.R")

SUMMARY_CSV  <- file.path(INPUTS_DIR, "topology", "networks_summary.csv")
NETWORKS_DIR <- file.path(INPUTS_DIR, "networks")
NET_PATTERN  <- "\\.tsv$"
OUT_DIR      <- PLOTS_DIR

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

# -----------------------------------------------------------------------------
# 2. Load & tidy global metrics
# -----------------------------------------------------------------------------

globals <- vroom::vroom(SUMMARY_CSV)

globals <- globals %>%
  mutate(
    Region    = sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", network)),
    Phenotype = sub("^.*_counts_([^_]+).*", "\\1", network),
    Phenotype = factor(Phenotype, levels = c("control", "AD")),
    Network   = paste0(Region, " ", Phenotype)
  ) %>%
  select(where(~ !all(is.na(.x))))   # drop all-NA columns

# -----------------------------------------------------------------------------
# 3. Define metrics of interest
# -----------------------------------------------------------------------------

metric_cols <- c(
  "n_nodes", "n_edges", "avg_path_len", "diameter",
  "global_density", "size_giant_component", "frac_giant_component",
  "n_components", "clustering_local_mean", "clustering_global",
  "assortativity", "deg_mean", "deg_median", "deg_sd",
  "kcore_max", "Q_modularity", "perc_targeted_50",
  "n_communities", "largest_community_size"
)

metric_labels <- c(
  n_nodes                = "Number of nodes",
  avg_path_len           = "Average path length",
  diameter               = "Diameter",
  global_density         = "Global density",
  size_giant_component   = "Giant component size",
  frac_giant_component   = "Fraction of giant component",
  n_components           = "Number of components",
  clustering_local_mean  = "Mean local clustering",
  clustering_global      = "Global clustering",
  assortativity          = "Assortativity",
  deg_mean               = "Mean degree",
  deg_median             = "Median degree",
  deg_sd                 = "SD degree",
  kcore_max              = "Maximum k-core",
  Q_modularity           = "Modularity (Q)",
  perc_targeted_50       = "Targeted attack (50%)",
  n_communities          = "Number of communities",
  largest_community_size = "Largest community size"
)

# -----------------------------------------------------------------------------
# 4. Compute signed log-transformed AD – control differences per region
#    f(x) = sign(x) * log10(|x| + 1)
# -----------------------------------------------------------------------------

diffs_mat_log <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Phenotype, values_from = where(is.numeric), names_sep = "_") %>%
  mutate(across(
    ends_with("_AD"),
    ~ . - get(sub("_AD$", "_control", cur_column())),
    .names = "diff_{.col}"
  )) %>%
  select(Region, starts_with("diff_")) %>%
  column_to_rownames("Region") %>%
  as.matrix() %>%
  { sign(.) * log10(abs(.) + 1) }

# Clean column names for display
colnames(diffs_mat_log) <- colnames(diffs_mat_log) %>%
  gsub("diff_|_AD", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

# -----------------------------------------------------------------------------
# 5. Euclidean distance matrix on standardized metrics
# -----------------------------------------------------------------------------

globals_scaled <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  unite("Network", Region, Phenotype, sep = " ") %>%
  column_to_rownames("Network") %>%
  select(where(~ var(.x, na.rm = TRUE) > 0 & !all(is.na(.x)))) %>%
  scale()

dist_matrix <- dist(globals_scaled, method = "euclidean") %>% as.matrix()

# Intra- vs inter-region distance summary (printed to console)
dist_df <- as.data.frame(dist_matrix) %>%
  rownames_to_column("Net1") %>%
  pivot_longer(-Net1, names_to = "Net2", values_to = "Distance") %>%
  filter(Net1 != Net2) %>%
  mutate(
    pair_type = if_else(
      str_extract(Net1, "^[^ ]+") == str_extract(Net2, "^[^ ]+"),
      "intra_region", "inter_region"
    )
  )

dist_summary <- dist_df %>%
  group_by(pair_type) %>%
  summarise(across(Distance, list(mean = mean, sd = sd, min = min, max = max, median = median)),
            .groups = "drop")

ratio_mean <- with(
  dist_summary,
  Distance_mean[pair_type == "intra_region"] / Distance_mean[pair_type == "inter_region"]
)

message(sprintf("\nIntra/inter-region mean distance ratio: %.3f", ratio_mean))
message(sprintf("Intra-region pairs are on average %.2fx %s than inter-region pairs.\n",
                ratio_mean,
                if_else(ratio_mean < 1,
                        "more similar (anatomical identity dominates)",
                        "more dissimilar")))

# -----------------------------------------------------------------------------
# 6. PCA of global network metrics
# -----------------------------------------------------------------------------

globals_metrics <- globals %>%
  select(all_of(metric_cols)) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0 & !all(is.na(.x))))

pca_res <- prcomp(globals_metrics, center = TRUE, scale. = TRUE)

var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

pca_df <- as.data.frame(pca_res$x) %>%
  mutate(Region = globals$Region, Phenotype = globals$Phenotype)

pca.p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, label = Region)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
  labs(
    title = " ",
    x     = sprintf("PC1 (%.1f%% var)", var_exp[1]),
    y     = sprintf("PC2 (%.1f%% var)", var_exp[2])
  ) +
  theme_cowplot() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

#Vis

pca.p

# Variable contributions (saved separately, not part of the composite panel)
surviving_metrics <- colnames(globals_metrics)
valid_labels      <- metric_labels[names(metric_labels) %in% surviving_metrics]

contrib <- sweep(pca_res$rotation^2, 2, colSums(pca_res$rotation^2), FUN = "/") * 100

contrib.p <- as.data.frame(contrib[, 1:2]) %>%
  rownames_to_column("variable") %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "contribution") %>%
  mutate(variable = factor(variable,
                           levels  = surviving_metrics,
                           labels  = valid_labels[surviving_metrics])) %>%
  group_by(PC) %>%
  arrange(desc(contribution), .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  ungroup() %>%
  ggplot(aes(x = variable, y = contribution, fill = PC)) +
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_manual(values = c("PC1" = "#1446A0", "PC2" = "#7d1538")) +
  labs(x = "Variable", y = "Contribution (%)",
       title = "Variable contributions to PC1 and PC2") +
  theme_cowplot()

# -----------------------------------------------------------------------------
# 7. Jaccard edge-set overlap between AD and control networks (per region)
# -----------------------------------------------------------------------------

# Helper: read a .tsv network file into an igraph object
read_network <- function(path) {
  df <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  
  # Auto-detect edgelist vs adjacency
  is_adj <- all(sapply(df[, -1, drop = FALSE], is.numeric))
  if (is_adj) {
    rownames(df) <- make.unique(as.character(df[, 1]))
    mat <- as.matrix(df[, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
  } else {
    g <- if (ncol(df) >= 3) igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
    else                igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
  }
  
  igraph::simplify(igraph::as_undirected(g, mode = "collapse"),
                   remove.multiple = TRUE, remove.loops = TRUE)
}

# Helper: edge set as sorted character pairs
get_edge_set <- function(g) {
  apply(igraph::as_edgelist(g), 1, function(e) paste(sort(e), collapse = "|"))
}

# Helper: Jaccard between two edge sets
jaccard_edges <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

# Load all networks
net_files <- list.files(NETWORKS_DIR, pattern = NET_PATTERN, full.names = TRUE)
stopifnot("At least 2 network files required." = length(net_files) >= 2)

# Extract region and phenotype from filenames
net_meta <- tibble(path = net_files) %>%
  mutate(
    base      = basename(path),
    region    = gsub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", base)),
    phenotype = if_else(grepl("_AD_", base), "AD", "Control"),
    Network   = paste0(region, "_", phenotype)
  ) %>%
  filter(!is.na(region))

networks   <- setNames(lapply(net_meta$path, read_network), net_meta$Network)
edge_sets  <- lapply(networks, get_edge_set)

# Jaccard similarity between all pairs within the same region (AD vs Control)
jaccard_region <- net_meta %>%
  inner_join(net_meta, by = "region", suffix = c("1", "2")) %>%
  filter(phenotype1 != phenotype2,
         phenotype1 == "Control") %>%           # one row per region
  mutate(
    jaccard = map2_dbl(Network1, Network2,
                       ~ jaccard_edges(edge_sets[[.x]], edge_sets[[.y]]))
  ) %>%
  select(region, jaccard) %>%
  distinct() %>%
  arrange(jaccard) %>%
  mutate(region = factor(region, levels = region))

print(jaccard_region)

# -----------------------------------------------------------------------------
# 8. Heatmap: signed log-transformed differences (AD vs control)
# -----------------------------------------------------------------------------

# Symmetric breaks so that zero difference maps to the neutral midpoint
diff_lim <- max(abs(diffs_mat_log), na.rm = TRUE)

heatmap_global <- grid::grid.grabExpr(
  pheatmap(
    diffs_mat_log,
    color        = colorRampPalette(c("firebrick4", "gray93", "cornflowerblue"))(100),
    breaks       = seq(-diff_lim, diff_lim, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale        = "none",
    border_color = "white",
    fontsize     = 20,
    main         = " "
  )
)

# -----------------------------------------------------------------------------
# 9. Assemble composite figure
# -----------------------------------------------------------------------------

# PCA + differences heatmap  (panels A & B)
final_plot <- plot_grid(
  pca.p,
  heatmap_global,
  labels     = c("A", "B"),
  label_size = 14,
  ncol       = 2
)

#Vis
final_plot

# -----------------------------------------------------------------------------
# 10. Save outputs
# -----------------------------------------------------------------------------

ggsave("FigureS1_contributions_pca.jpeg", contrib.p, width = 15, height = 13, dpi = 300)
ggsave("FigureS1_contributions_pca.pdf",  contrib.p, width = 15, height = 13, dpi = 300)

ggsave("pca_metrics.jpeg",         pca.p,       width =  5, height =  5, dpi = 300)
ggsave("pca_metrics.pdf",          pca.p,       width =  5, height =  5, dpi = 300)

ggsave("heatmap_global_diffs.jpeg", heatmap_global, width = 8, height = 7, dpi = 300)
ggsave("heatmap_global_diffs.pdf",  heatmap_global, width = 8, height = 7, dpi = 300)

ggsave("Figure1_final_global_figure.jpeg", final_plot, width = 16, height = 7, dpi = 300)
ggsave("Figure1_final_global_figure.pdf",  final_plot, width = 16, height = 7, dpi = 300)

# =============================================================================
# END
# =============================================================================