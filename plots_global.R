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
  "broom"
)

#Get data

globals <- vroom::vroom("~/Desktop/local_work/results_topos/networks_summary.csv")
globals$Region <-sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))

#Function to plot
make_metric_plot <- function(df, metric, title, ylab, color) {
  ggplot(df, aes(x = as.factor(network), y = .data[[metric]], group = 1)) +
    geom_line() +
    geom_point(color = color, size = 2.5) +
    labs(title = title, x = "Network", y = ylab) +
    scale_y_continuous(limits = c(min(df[[metric]]), max(df[[metric]]))) +
    theme_pubclean() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

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
      "clustering_global", "clustering_local_mean",
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

#some statistics

metric_cols <- c("n_nodes", "n_edges", "n_components", "n_communities", "deg_mean",
                 "avg_path_len", "clustering_global", "clustering_local_mean",
                 "diameter", "global_density", "size_giant_component",
                 "frac_giant_component", "Q_modularity", "perc_targeted_50")


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
