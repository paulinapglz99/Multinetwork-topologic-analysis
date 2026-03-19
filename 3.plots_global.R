#!/usr/bin/env Rscript
# =============================================================================
# plots_global.R
# Generate topological network metric plots and heatmaps
# Requires: output of topological_analysis.R (networks_summary.csv)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Dependencies
# -----------------------------------------------------------------------------

if (!requireNamespace("pacman",   quietly = TRUE)) install.packages("pacman",   repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")

pacman::p_load(
  igraph, data.table,
  ggplot2, tidyverse,
  optparse, tools, purrr,
  ggpubr, broom, pheatmap,
  vroom, ggrepel, cowplot,
  ggdendro, patchwork
)

# -----------------------------------------------------------------------------
# 1. Load & tidy data
# -----------------------------------------------------------------------------

globals <- vroom::vroom("~/Desktop/local_work/fomo_networks/results_topos_louvain/networks_summary.csv")

# Parse region and phenotype from network IDs
globals$Region    <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))
globals$Network   <- paste0(globals$Region, " ", globals$Phenotype)

# Drop all-NA columns
globals <- globals %>% select(where(~ !all(is.na(.x))))

# Fix column names (must happen BEFORE defining metric_cols)
colnames(globals)
# [1] "file"                   "n_nodes"                "n_edges"                "avg_path_len"           "diameter"              
# [6] "global_density"         "size_giant_component"   "frac_giant_component"   "n_components"           "clustering_local_mean" 
# [11] "clustering_global"      "assortativity"          "deg_mean"               "deg_median"             "deg_sd"                
# [16] "kcore_max"              "Q_modularity"           "perc_targeted_50"       "perc_random_50"         "n_communities"         
# [21] "largest_community_size" "largest_community_id"   "network"                "Region"                 "Phenotype"             
# [26] "Network"    

# globals <- globals %>%
#   dplyr::rename(
#     global_clustering = clustering_global,
#     degree_mean = deg_mean,
#     degree_median = deg_median,
#     degree_sd = deg_sd
#   )

# -----------------------------------------------------------------------------
# 2. Define metrics
# -----------------------------------------------------------------------------

metric_cols <- c(
  "n_nodes", "n_edges", "avg_path_len", "diameter",
  "global_density", "size_giant_component", "frac_giant_component",
  "n_components", "clustering_local_mean", "clustering_global",
  "assortativity", "deg_mean", "deg_median", "deg_sd",
  "kcore_max", "Q_modularity", "perc_targeted_50",
  "n_communities", "largest_community_size"
)

length(metric_cols)

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
# 3. Compute AD – control differences per region
# -----------------------------------------------------------------------------

means <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop")

means_wide <- means %>%
  pivot_wider(
    names_from  = Phenotype,
    values_from = where(is.numeric),
    names_sep   = "_"
  )

diffs <- means_wide %>%
  mutate(across(
    ends_with("_AD"),
    ~ . - get(sub("_AD$", "_control", cur_column())),
    .names = "diff_{.col}"
  )) %>%
  select(Region, starts_with("diff_"))

# Log-transform signed differences for visualisation
diffs <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()

diffs_mat_log <- sign(diffs) * log10(abs(diffs) + 1)

# Clean up column names for display
colnames(diffs_mat_log) <- colnames(diffs_mat_log) %>%
  gsub("diff_|_AD", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

# Hierarchical clustering of regions by their difference profiles
dist_mat  <- dist(diffs_mat_log)
hc_diff   <- hclust(dist_mat, method = "ward.D2")
region_order <- hc_diff$labels[hc_diff$order]

# -----------------------------------------------------------------------------
# 4. Heatmap A – global metric differences (AD vs control)
# -----------------------------------------------------------------------------

heatmap_global <- grid::grid.grabExpr({
  pheatmap(
    diffs_mat_log,
    color        = colorRampPalette(c("firebrick4", "gray93", "cornflowerblue"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale        = "none",
    border_color = "white",
    fontsize     = 20,
    main         = " "
  )
})

# -----------------------------------------------------------------------------
# 5. Correlation heatmap between network metrics
# -----------------------------------------------------------------------------

cor_mat <- cor(globals %>% select(all_of(metric_cols)), use = "pairwise.complete.obs")

cor_mat.p <- pheatmap(
  cor_mat,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main  = "Correlation between Network Metrics"
)
cor_mat.p

# -----------------------------------------------------------------------------
# 6. PCA of network metrics
# -----------------------------------------------------------------------------

# Remove zero-variance columns before PCA
globals_metrics <- globals %>%
  select(all_of(metric_cols)) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0 & !all(is.na(.x))))

pca_res <- prcomp(globals_metrics, center = TRUE, scale. = TRUE)

# Scores
pca_df <- as.data.frame(pca_res$x) %>%
  mutate(
    Region    = globals$Region,
    Phenotype = globals$Phenotype
  )

# Variance explained
var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pc1_var <- round(var_exp[1], 1)
pc2_var <- round(var_exp[2], 1)

pca.p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, label = Region)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
  labs(
    title = " ",
    x     = paste0("PC1 (", pc1_var, "% var)"),
    y     = paste0("PC2 (", pc2_var, "% var)")
  ) +
  theme_cowplot() +
  theme(
    legend.position  = "top",
    plot.title       = element_text(hjust = 0.5)
  )

#Vis
pca.p

# Variable contributions to PC1 and PC2
contrib <- sweep(pca_res$rotation^2, 2, colSums(pca_res$rotation^2), FUN = "/") * 100

# Use only metrics that survived the zero-variance filter
surviving_metrics <- colnames(globals_metrics)
valid_labels      <- metric_labels[names(metric_labels) %in% surviving_metrics]

contrib.df <- as.data.frame(contrib[, 1:2]) %>%
  rownames_to_column("variable") %>%
  pivot_longer(
    cols      = starts_with("PC"),
    names_to  = "PC",
    values_to = "contribution"
  ) %>%
  mutate(variable = factor(
    variable,
    levels = surviving_metrics,
    labels = valid_labels[surviving_metrics]
  )) %>%
  group_by(PC) %>%
  arrange(desc(contribution), .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = unique(variable))) %>%
  ungroup()

contrib.p <- ggplot(contrib.df, aes(x = variable, y = contribution, fill = PC)) +
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_manual(values = c("PC1" = "#1446A0", "PC2" = "#7d1538")) +
  labs(
    x     = "Variable",
    y     = "Contribution (%)",
    title = "Variable contributions to PC1 and PC2"
  ) +
  theme_cowplot()

contrib.p

# # -----------------------------------------------------------------------------
# 7. Euclidean distance matrix
# # -----------------------------------------------------------------------------

globals_summary <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  unite("Network", Region, Phenotype, sep = " ") %>%
  column_to_rownames("Network") %>%
  select(where(~ var(.x, na.rm = TRUE) > 0 & !all(is.na(.x)))) %>%  # remove zero-var cols
  scale()

dist_matrix <- dist(globals_summary, method = "euclidean") %>% as.matrix()

# Hierarchical clustering of the distance matrix
hc_dist <- hclust(dist(dist_matrix), method = "ward.D2")
ord     <- hc_dist$order
k       <- 4
clusters      <- cutree(hc_dist, k = k)
cluster_order <- clusters[ord]

dist_reordered <- dist_matrix[ord, ord]

# Dendrogram
dendro <- ggdendrogram(hc_dist, rotate = TRUE, theme_dendro = FALSE) +
  theme_cowplot() +
  labs(title = "", x = "", y = "")

# Heatmap of distance matrix
breaks <- which(diff(cluster_order) != 0)

dist_long <- as.data.frame(dist_reordered) %>%
  rownames_to_column("Network1") %>%
  pivot_longer(-Network1, names_to = "Network2", values_to = "Distance") %>%
  mutate(
    Network1 = factor(Network1, levels = rownames(dist_reordered)),
    Network2 = factor(Network2, levels = colnames(dist_reordered))
  ) %>%
  filter(as.integer(Network1) <= as.integer(Network2))

heatmap_dist <- ggplot(dist_long, aes(Network1, Network2, fill = Distance)) +
  geom_tile() +
  scale_fill_viridis_c(
    guide = guide_colorbar(
      position      = "top",
      barheight     = unit(10,  "pt"),
      barwidth      = unit(140, "pt"),
      title.position = "top"
    )
  ) +
  coord_fixed() +
  geom_hline(yintercept = breaks + 0.5, color = "black", linewidth = 0.7) +
  geom_vline(xintercept = breaks + 0.5, color = "black", linewidth = 0.7) +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(angle = 90, hjust = 1, size = 15),
    axis.text.y   = element_text(size = 15),
    legend.position     = "top",
    legend.title        = element_text(size = 15),
    legend.text         = element_text(size = 15)
  ) +
  labs(title = " ", x = "", y = "")

#vis

heatmap_dist

# -----------------------------------------------------------------------------
# 7b. Intra- vs inter-region distance summary
# -----------------------------------------------------------------------------

dist_df <- as.data.frame(as.matrix(dist_matrix)) %>%
  rownames_to_column("Net1") %>%
  pivot_longer(-Net1, names_to = "Net2", values_to = "Distance") %>%
  filter(Net1 != Net2) %>%
  mutate(
    Region1   = str_extract(Net1, "^[^ ]+"),
    Region2   = str_extract(Net2, "^[^ ]+"),
    pair_type = if_else(Region1 == Region2, "intra_region", "inter_region")
  )

dist_summary <- dist_df %>%
  group_by(pair_type) %>%
  summarise(
    mean = mean(Distance),
    min  = min(Distance),
    max  = max(Distance),
    sd   = sd(Distance),
    median = median(Distance),
    .groups = "drop"
  )

print(dist_summary)

ratio_mean <- dist_summary$mean[dist_summary$pair_type == "intra_region"] /
  dist_summary$mean[dist_summary$pair_type == "inter_region"]

cat(sprintf("\nIntra/inter-region mean distance ratio: %.3f\n", ratio_mean))
cat(sprintf("Interpretation: intra-region pairs are on average %.2fx %s than inter-region pairs\n",
            ratio_mean,
            if_else(ratio_mean < 1, "more similar (anatomical identity dominates)", "more dissimilar")))

# -----------------------------------------------------------------------------
# 8. Per-metric line plots
# -----------------------------------------------------------------------------

# make_metric_plot <- function(df, metric, ylab) {
#   ggplot(df, aes(
#     x      = as.factor(Region),
#     y      = .data[[metric]],
#     color  = Phenotype,
#     group  = Phenotype
#   )) +
#     geom_line(linewidth = 0.8, alpha = 0.8) +
#     geom_point(size = 2.5) +
#     labs(title = " ", x = "", y = ylab) +
#     scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
#     scale_y_continuous(limits = c(0, max(df[[metric]], na.rm = TRUE))) +
#     theme_light() +
#     theme(
#       legend.position = "none",
#       panel.spacing   = unit(0.2, "lines"),
#       plot.margin     = margin(1, 1, 1, 1, "mm"),
#       axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1)
#     )
# }
# 
# plot_specs <- data.frame(
#   var  = metric_cols,
#   ylab = colnames(diffs_mat_log),
#   stringsAsFactors = FALSE
# )
# 
# metric_plots <- map2(
#   plot_specs$var,
#   plot_specs$ylab,
#   ~ make_metric_plot(globals, .x, .y)
# )
# 
# metric_grid <- cowplot::plot_grid(plotlist = metric_plots, ncol = 2)
# 
# #Vis
# 
# metric_grid

# -----------------------------------------------------------------------------
# 9. Assemble final figures
# -----------------------------------------------------------------------------

# Figure 1: distance heatmap + metrics PCA
top_panel <- plot_grid(
  heatmap_dist,
  pca.p,
  labels     = c("A", "B"),
  label_size = 14,
  ncol       = 2,
  rel_widths = c(1, 1)
)

# Figure 2: differences heatmap + global connectivity (jaccard_global_connect must exist)
# middle_panel <- plot_grid(
#   heatmap_global,
#   jaccard_global_connect,
#   labels     = c("C", "D"),
#   label_size = 14,
#   ncol       = 2
# )

# Final composite
# final_plot <- plot_grid(top_panel, middle_panel, ncol = 1)

# -----------------------------------------------------------------------------
# 10. Save outputs
# -----------------------------------------------------------------------------

ggsave("contributions_pca.jpeg",     contrib.p,   width = 15, height = 13, dpi = 300)
ggsave("contributions_pca.pdf",      contrib.p,   width = 15, height = 13, dpi = 300)
ggsave("pca_metrics.jpeg",           pca.p,       width =  5, height =  5, dpi = 300)
ggsave("pca_metrics.pdf",            pca.p,       width =  5, height =  5, dpi = 300)
ggsave("metric_plots_grid.jpeg",     metric_grid, width = 12, height = 18, dpi = 300)
ggsave("metric_plots_grid.pdf",      metric_grid, width = 12, height = 15, dpi = 300)
ggsave("final_global_figure.jpeg",   top_panel,   width = 16, height =  9, dpi = 300)
ggsave("final_global_figure.pdf",    top_panel,   width = 16, height =  9, dpi = 300)

# =============================================================================
# END
# =============================================================================