#!/usr/bin/env Rscript

#plots_global.R
#Usage: https://github.com/paulinapglz99/Multinetwork-topologic-analysis

if (!requireNamespace("pacman", quietly = FALSE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = FALSE)) install.packages("optparse", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "igraph",
  "data.table",
  #"future.apply",
  "ggplot2",
  "tidyverse",
  #"jsonlite",
  #"stringr",
  "optparse", 
  "tools",
  "purrr",
  "ggpubr",
  "broom", 
  "pheatmap"
)

#Get data

globals <- vroom::vroom("~/Desktop/local_work/results_topos/networks_summary.csv")
globals$Region <-sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))
names(globals$clustering_global)
colnames(globals)[11] <- "global_clustering"

metric_cols <- c("n_nodes", "n_edges", "n_components", 
                 "n_communities", "degree_mean","avg_path_length",
                 "global_clustering", "mean_clustering_local",
                 "diameter", "global_density", 
                 "size_giant_component","frac_giant_component",
                 "Q_modularity", "perc_targeted_50")

#Get means per phenotype and region
means <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(metric_cols, mean, na.rm = TRUE), .groups = "drop")

#Pivot and get columns per phenotype
means_wide <- means %>%
  pivot_wider(
    names_from = Phenotype,
    values_from = where(is.numeric),
    names_sep = "_"
  )

#Get differences for all numeric cols
diffs <- means_wide %>%
  mutate(across(ends_with("_AD"), 
                ~ . - get(sub("_AD$", "_control", cur_column())),
                .names = "diff_{.col}")) %>%
  select(Region, starts_with("diff_"))

# === 4. Matriz numérica para clustering ===
diffs_mat <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()

# log-transformación simétrica
diffs_mat_log <- sign(diffs_mat) * log10(abs(diffs_mat) + 1)

# === 5. Clustering jerárquico ===
dist_mat <- dist(diffs_mat_log)                 # distancia euclídea por región
hc <- hclust(dist_mat, method = "ward.D2")      # clustering jerárquico

# === 6. Reordenar las regiones según el dendrograma ===
region_order <- hc$labels[hc$order]

diffs_long <- diffs %>%
  pivot_longer(-Region, names_to = "Metric", values_to = "Difference") %>%
  mutate(
    Difference_log = sign(Difference) * log10(abs(Difference) + 1),
    Metric_name = Metric %>% str_remove_all("diff_|_AD") %>% 
      str_replace_all("_", " ") %>%
      str_to_title()
  ) %>%
  mutate(Region = factor(Region, levels = region_order))
# 
# #Order and reorder data to plot heatmap
# diffs_long <- diffs %>%
#   pivot_longer(-Region, names_to = "Metric", values_to = "Difference") %>%
#   mutate(Metric = gsub("diff_|_AD", "", Metric)) %>%
#   mutate(
#     Difference_log = sign(Difference) * log10(abs(Difference) + 1)
#   ) %>%
#   mutate(Metric_name = Metric %>%
#            str_replace_all("_", " ") %>%   
#            str_to_title()                  
#   )

#Plot heatmap
# global_heat.p <- ggplot(diffs_long, aes(x = Metric_name, y = Region, fill = Difference_log)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "cornflowerblue", mid = "gray93", high = "red4", midpoint = 0) +
#   theme_minimal(base_size = 12) +
#   labs(
#     title = "Differences in network metrics",
#     fill = "Δ (AD - Control)",
#     x = " ", y = " "
#   ) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# global_heat.p
# 
# ggsave("global_heat.jpeg", 
#        global_heat.p, 
#        width = 7, 
#        height = 10, 
#        units = "in", 
#        dpi= 300)

#Plot ordered by clustering
# global_heat.p <- ggplot(diffs_long, aes(x = Metric_name, y = Region, fill = Difference_log)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "cornflowerblue", mid = "gray93", high = "red4", midpoint = 0) +
#   theme_minimal(base_size = 12) +
#   labs(
#     title = "Differences in Network Metrics (AD - Control)",
#     fill = expression(Delta~"(log10)"),
#     x = " ", y = " "
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     panel.grid = element_blank(),
#     plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
#   )
# 
# global_heat.p

#Function to plot
# make_metric_plot <- function(df, metric, title, ylab, color) {
#   ggplot(df, aes(x = as.factor(network), y = .data[[metric]], group = 1)) +
#     geom_line() +
#     geom_point(color = color, size = 2.5) +
#     labs(title = title, x = "Network", y = ylab) +
#     scale_y_continuous(limits = c(min(df[[metric]]), max(df[[metric]]))) +
#     theme_pubclean() +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
#     )
# }

# Crear matriz log-transformada
diffs_mat <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()
diffs_mat_log <- sign(diffs_mat) * log10(abs(diffs_mat) + 1)

# Limpiar nombres de métricas
colnames(diffs_mat_log) <- gsub("diff_|_AD", "", colnames(diffs_mat_log))
colnames(diffs_mat_log) <- gsub("_", " ", colnames(diffs_mat_log))
colnames(diffs_mat_log) <- tools::toTitleCase(colnames(diffs_mat_log))

# Heatmap con clustering en filas y columnas
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

#Plot metrics in other way
make_metric_plot <- function(df, metric, title, ylab) {
  ggplot(df, aes(
    x = as.factor(Region),
    y = .data[[metric]],
    group = Phenotype,
    color = Phenotype
  )) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~Phenotype, scales = "free_y") +
    labs(title = title, x = " ", y = ylab) +
    # Colores fijos por fenotipo
    scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
    scale_y_continuous(limits = c(min(df[[metric]]), max(df[[metric]]))) +
    theme_pubclean() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#Specifications to plot metrics
plot_specs <- data.frame(
  var = c(
    "n_nodes", "n_edges", 
    "n_components", "n_communities",
    "deg_mean", "avg_path_len",
    "global_clustering", "clustering_local_mean",
    "diameter", "global_density",
    "size_giant_component", "frac_giant_component",
    "Q_modularity", "perc_targeted_50"
  ),
  # title = c(
  #   "Vertices", "Edges", "Components",
  #   "Communities", "Mean degree of nodes", "Average path length",
  #   "Global clustering", "Local clustering coefficient", "Diameter", "Global density",
  #   "Size of the giant component", "Fraction of giant component", "Q modularity", "Percolation threshold"),
  ylab = c(
    "Vertices", "Edges", "Components",
    "Communities", "Mean degree", "Average path length",
    "Global clustering", "Local clustering", "Diameter", "Density",
    "Size", "Fraction", "Q", "Percolation threshold"
  ),
  color = c(
    "red4", "cornflowerblue",
    "mediumpurple", "slateblue", "maroon", "navyblue", "olivedrab",
    "lightpink4", "red4", "tomato4", "lightsalmon", "firebrick", "seagreen", "yellowgreen"
  ),
  stringsAsFactors = FALSE
)

#Plot all metrics
plots <- map2(plot_specs$var, seq_along(plot_specs$var), function(varname, idx) {
  make_metric_plot(globals, varname, plot_specs$title[idx], plot_specs$ylab[idx])
})

x <- ggarrange(plotlist = plots, ncol = 2, nrow = 7)
x

ggsave("global_metrics.pdf", 
       x, 
       width = 12, 
       height = 15, 
       units = "in", 
       dpi= 300)

####OPTION 2 ####

make_metric_plot <- function(df, metric, title, ylab) {
  ggplot(df, aes(
    x = as.factor(Region),
    y = .data[[metric]],
    color = Phenotype
  )) +
    geom_segment(aes(
      x = as.factor(Region),
      xend = as.factor(Region),
      y = 0,                
      yend = .data[[metric]]  
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

plots <- map2(plot_specs$var, seq_along(plot_specs$var), function(varname, idx) {
  make_metric_plot(globals, varname, plot_specs$title[idx], plot_specs$ylab[idx])
})

x <- ggarrange(plotlist = plots, ncol = 2, nrow = 7, align = "hv")
x

#in pdf
ggsave("global_metrics_lollipop.pdf", 
       x, 
       width = 12, 
       height = 15, 
       units = "in", 
       dpi= 300)

#in jpeg
ggsave("global_metrics_lollipop.jpeg", 
       x, 
       width = 12, 
       height = 15, 
       units = "in", 
       dpi= 300)

#some statistics
results <- metric_cols %>%
  map_df(function(metric) {
    df <- globals %>%
      select(Region, Phenotype, all_of(metric)) %>%
      filter(!is.na(.data[[metric]]))
    
    # t-test por métrica
    test <- t.test(formula(paste(metric, "~ Phenotype")), data = df)
    
    tidy(test) %>%
      mutate(metric = metric)
  })
results

#Save stats results


results %>%
  select(metric, estimate1, estimate2, statistic, p.value) %>%
  arrange(p.value)
