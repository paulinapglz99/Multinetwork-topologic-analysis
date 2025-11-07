#!/usr/bin/env Rscript
# ============================================================
# Script: plots_global.R
# Purpose: Generate topological network metric plots and heatmaps
# Author: Alejandra Paulina Pérez González
# Repository: https://github.com/paulinapglz99/Multinetwork-topologic-analysis
# ============================================================

# ------------------ #
#   Load packages    #
# ------------------ #
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")

pacman::p_load(
  igraph, data.table, ggplot2, tidyverse,
  optparse, tools, purrr, ggpubr, broom, pheatmap, vroom
)

# ------------------ #
#   Load data        #
# ------------------ #
globals <- vroom::vroom("~/Desktop/local_work/results_topos/networks_summary.csv")

# Parse region and phenotype names from network IDs
globals$Region <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))

# Fix column name for global clustering
colnames(globals)[11] <- "global_clustering"
colnames(globals)[13] <-"degree_mean"
colnames(globals)[14] <- "degree_median"
colnames(globals)[15] <- "degree_sd"

colnames(globals)

# ------------------ #
#   Define metrics   #
# ------------------ #
metric_cols <- colnames(globals)[2:21]
# c(
#   "n_nodes", "n_edges", "n_components", "n_communities",
#   "degree_mean", "avg_path_length", "global_clustering",
#   "mean_clustering_local", "diameter", "global_density",
#   "size_giant_component", "frac_giant_component",
#   "Q_modularity", "perc_targeted_50"
# )

# ---------------------------- #
#   Compute means per group    #
# ---------------------------- #
means <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop")

means_wide <- means %>%
  pivot_wider(
    names_from = Phenotype,
    values_from = where(is.numeric),
    names_sep = "_"
  )

# Compute differences (AD - control)
diffs <- means_wide %>%
  mutate(across(ends_with("_AD"),
                ~ . - get(sub("_AD$", "_control", cur_column())),
                .names = "diff_{.col}")) %>%
  select(Region, starts_with("diff_"))

# ---------------------------- #
#   Hierarchical clustering    #
# ---------------------------- #
diffs_mat <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()

diffs_mat_log <- sign(diffs_mat) * log10(abs(diffs_mat) + 1)
dist_mat <- dist(diffs_mat_log)
hc <- hclust(dist_mat, method = "ward.D2")
region_order <- hc$labels[hc$order]

# ---------------------------- #
#   Prepare data for ggplot    #
# ---------------------------- #
diffs_long <- diffs %>%
  pivot_longer(-Region, names_to = "Metric", values_to = "Difference") %>%
  mutate(
    Difference_log = sign(Difference) * log10(abs(Difference) + 1),
    Metric_name = Metric %>%
      str_remove_all("diff_|_AD") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    Region = factor(Region, levels = region_order)
  )

# ---------------------------- #
#   Heatmap (pheatmap)         #
# ---------------------------- #
colnames(diffs_mat_log) <- colnames(diffs_mat_log) %>%
  gsub("diff_|_AD", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

pheatmap(
  diffs_mat_log,
  color = colorRampPalette(c("cornflowerblue", "gray93", "red4"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "white",
  fontsize = 10,
  main = "Differences in Network Metrics (AD - Control, log10)"
)

# ---------------------------- #
#   Function to plot metrics   #
# ---------------------------- #
make_metric_plot <- function(df, metric, title, ylab) {
  ggplot(df, aes(
    x = as.factor(Region),
    y = .data[[metric]],
    color = Phenotype
  )) +
    geom_segment(aes(
      xend = as.factor(Region),
      y = 0, yend = .data[[metric]]
    ), size = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~Phenotype, scales = "free_y", ncol = 2) +
    labs(title = title, x = "", y = ylab) +
    scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
    scale_y_continuous(limits = c(0, max(df[[metric]]))) +
    theme_pubclean() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.2, "lines"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

# ---------------------------- #
#   Define plot specifications #
# ---------------------------- #
plot_specs <- data.frame(
  var = metric_cols,
  ylab = c(
    "Vertices", "Edges", "Components", "Communities", "Mean degree",
    "Average path length", "Global clustering", "Local clustering",
    "Diameter", "Density", "Size", "Fraction", "Q", "Percolation threshold"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------- #
#   Generate lollipop plots    #
# ---------------------------- #
plots <- map2(plot_specs$var, seq_along(plot_specs$var), function(varname, idx) {
  make_metric_plot(globals, varname, plot_specs$ylab[idx], plot_specs$ylab[idx])
})

x <- ggarrange(plotlist = plots, ncol = 2, nrow = 7, align = "hv")

ggsave("global_metrics_lollipop.pdf", x, width = 12, height = 15, units = "in", dpi = 300)
ggsave("global_metrics_lollipop.jpeg", x, width = 12, height = 15, units = "in", dpi = 300)

# ---------------------------- #
#   Statistical comparisons    #
# ---------------------------- #
results <- metric_cols %>%
  map_df(function(metric) {
    df <- globals %>%
      select(Region, Phenotype, all_of(metric)) %>%
      filter(!is.na(.data[[metric]]))
    test <- t.test(formula(paste(metric, "~ Phenotype")), data = df)
    tidy(test) %>% mutate(metric = metric)
  })

results %>%
  select(metric, estimate1, estimate2, statistic, p.value) %>%
  arrange(p.value)

# Save stats results
write.csv(results, "network_metric_stats.csv", row.names = FALSE)
