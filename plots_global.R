#!/usr/bin/env Rscript
#plots_global.R
#Generate topological network metric plots and heatmaps

#Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")

pacman::p_load(
  igraph, data.table, ggplot2, tidyverse,
  optparse, tools, purrr, ggpubr, broom, pheatmap, vroom, ggrepel, cowplot, ggrepel
)

#Get data
globals <- vroom::vroom("~/Desktop/local_work/results_topos/networks_summary.csv")

#Parse region and phenotype names from network IDs
globals$Region <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))
globals <- globals %>%  select(where(~ !all(is.na(.x))))
#globals <-globals %>% filter(!perc_random_50)
#Fix column name for global clustering
colnames(globals)[11] <- "global_clustering"
colnames(globals)[13] <-"degree_mean"
colnames(globals)[14] <- "degree_median"
colnames(globals)[15] <- "degree_sd"

colnames(globals)

#Define metrics
#metric_cols <- colnames(globals)[2:20]
metric_cols <-  c("n_nodes"    ,            "n_edges"       ,         "avg_path_len"        ,   "diameter"          ,     "global_density"        ,
"size_giant_component"  , "frac_giant_component",   "n_components"       ,    "clustering_local_mean" , "global_clustering"     ,
"assortativity"  ,        "degree_mean"       ,     "degree_median"   ,           "kcore_max" ,            
 "Q_modularity"    ,       "perc_targeted_50"    ,   "n_communities"     ,     "largest_community_size"
 )

#Compute means per group
means <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop")

means_wide <- means %>%
  pivot_wider(
    names_from = Phenotype,
    values_from = where(is.numeric),
    names_sep = "_")

#Compute differences (AD - control)
diffs <- means_wide %>%
  mutate(across(ends_with("_AD"),
                ~ . - get(sub("_AD$", "_control", cur_column())),
                .names = "diff_{.col}")) %>%
  select(Region, starts_with("diff_"))

#Hierarchical clustering

diffs_mat <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()

diffs_mat_log <- sign(diffs_mat) * log10(abs(diffs_mat) + 1)
dist_mat <- dist(diffs_mat_log)
hc <- hclust(dist_mat, method = "ward.D2")
region_order <- hc$labels[hc$order]

#Prepare data for ggplot
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

#Plot Heatmap (pheatmap)
colnames(diffs_mat_log) <- colnames(diffs_mat_log) %>%
  gsub("diff_|_AD", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

heatmap_global<- grid::grid.grabExpr({
  pheatmap(
  diffs_mat_log,
  color = colorRampPalette(c("cornflowerblue", "gray93", "red4"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "white",
  fontsize = 10,
  main = " "
)
})

#Vis
heatmap_global

#Save plot
ggsave("heatmap_global-a.jpeg",
       plot = heatmap_global,
       device = "jpeg",
       width = 15, 
       height = 10,
       units = "in",
       dpi = 300
       )

############PLOT CORRELATIONS##############
#See if certain metrics are co-varying

cor_mat <- cor(globals %>% select(all_of(metric_cols)), use = "pairwise.complete.obs")

pheatmap(cor_mat,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Correlation between Network Metrics")

############PCA of regions by metrics################

globals_scaled <- globals %>%
  select(all_of(metric_cols)) %>%
  scale()

#PCA
pca_res <- prcomp(globals_scaled, center = TRUE, scale. = TRUE)

#Score df
pca_df <- as.data.frame(pca_res$x) %>%
  mutate(Region = globals$Region,
         Phenotype = globals$Phenotype)

#variance
var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pc1 <- round(var_exp[1], 1)
pc2 <- round(var_exp[2], 1)

#Plot PCA
safe_ellipse <- function(...) {
  tryCatch(stat_ellipse(...), error = function(e) NULL)
}
pca <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                          color = Phenotype,
                          label = Region)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
  labs(
    title = " ",
    x = paste0("PC1 (", pc1, "% var)"),
    y = paste0("PC2 (", pc2, "% var)")
  ) +
  theme_cowplot() +
  #stat_ellipse(aes(group = Region), level = 0.95, geom = "polygon", alpha = 0.2) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

#Vis
pca

#Save plot
ggsave("pca-globals.jpeg",
       plot = pca,
       device = "jpeg",
       width = 5, 
       height = 5,
       units = "in",
       dpi = 300
)

############PLOT ALL METRICS ##############

#Function to plot metrics
make_metric_plot <- function(df, metric, title, ylab) {
  ggplot(df, aes(
    x = as.factor(Region),
    y = .data[[metric]],
    color = Phenotype,
    group = Phenotype 
  )) +
    # geom_segment(aes(
    #   xend = as.factor(Region),
    #   y = 0, yend = .data[[metric]]
    # ), size = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2.5) +
    #facet_wrap(~Phenotype, scales = "free_y", ncol = 2) +
    labs(title = " ", x = "", y = ylab) +
    scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
    scale_y_continuous(limits = c(0, max(df[[metric]]))) +
    theme_light() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.2, "lines"),
      plot.margin = margin(1, 1, 1, 1, "mm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#Define plot specifications
plot_specs <- data.frame(
  var = metric_cols,
  ylab = colnames(diffs_mat_log),
  stringsAsFactors = FALSE
)

#Generate lollipop plots
plots <- map2(plot_specs$var, seq_along(plot_specs$var), function(varname, idx) {
  make_metric_plot(globals, varname, plot_specs$ylab[idx], plot_specs$ylab[idx])
})

#Grid plots
x <-cowplot::plot_grid(plotlist = plots, ncol = 2)
#Vis
x

#Save plots
ggsave("2-global_metrics_lollipop.pdf", x, width = 12, height = 15, units = "in", dpi = 300)
ggsave("3-global_metrics_lollipop.jpeg", x, width = 12, height = 18, units = "in", dpi = 300)

#Final grand plot
right_panel <- cowplot::plot_grid(
  heatmap_global,
  pca,
  ncol = 1,
  labels = c(" ", "C"),
  rel_heights = c(1.2, 1)
)
right_panel

grand <- cowplot::plot_grid(
  x,              
  right_panel,
  labels = c("A", "B"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(1.5, 1)  # el lollipop más ancho
)

#Vis
grand

#Save plot
ggsave("1-grand_global_metrics-pres.jpeg", grand, 
       width = 18, 
       height = 17,
       units = "in",
       dpi = 300)

#END