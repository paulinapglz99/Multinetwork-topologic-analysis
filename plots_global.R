#!/usr/bin/env Rscript
#plots_global.R
#Generate topological network metric plots and heatmaps
#This script needs the output of topological_analysis.R

#Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")

pacman::p_load(
  igraph, data.table, 
  ggplot2, tidyverse,
  optparse, tools, purrr,
  ggpubr, broom, pheatmap, 
  vroom, ggrepel, cowplot,
  ggdendro,patchwork
)

#Get data
globals <- vroom::vroom("~/Desktop/local_work/fomo_networks/results_topos_louvain/networks_summary.csv")

#Parse region and phenotype names from network IDs
globals$Region <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", globals$network))
globals$Phenotype <- sub("^.*_counts_([^_]+).*", "\\1", globals$network)
globals$Phenotype <- factor(globals$Phenotype, levels = c("control", "AD"))
globals$Network <- paste0(globals$Region, " ", globals$Phenotype) 
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
metric_cols <-  c("n_nodes"    ,     "avg_path_len"        ,   "diameter"          ,     "global_density"        ,
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

#Hierarchical clustering of differences

diffs_mat <- diffs %>%
  column_to_rownames("Region") %>%
  as.matrix()

diffs_mat_log <- sign(diffs_mat) * log10(abs(diffs_mat) + 1)
dist_mat <- dist(diffs_mat_log)
hc_diff <- hclust(dist_mat, method = "ward.D2")
region_order <- hc_diff$labels[hc_diff$order]

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

#Plot pheatmap
colnames(diffs_mat_log) <- colnames(diffs_mat_log) %>%
  gsub("diff_|_AD", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

heatmap_global<- grid::grid.grabExpr({
  pheatmap(
    diffs_mat_log,
    color = colorRampPalette(c("firebrick4", "gray93", "cornflowerblue"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    border_color = "white",
    fontsize = 20,
    main = " "
  )
})

#Vis
heatmap_global

#Save plot
ggsave("heatmap_global_diffs.jpeg",
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
# 
# globals_metrics <- globals %>% select(all_of(metric_cols))
# 
# #PCA
# pca_res <- prcomp(globals_metrics, center = TRUE, scale. = TRUE)
# 
# #Score df
# pca_df <- as.data.frame(pca_res$x) %>%
#   mutate(Region = globals$Region,
#          Phenotype = globals$Phenotype)
# 
# #variance
# var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
# pc1 <- round(var_exp[1], 1)
# pc2 <- round(var_exp[2], 1)
# 
# #Plot PCA
# safe_ellipse <- function(...) {
#   tryCatch(stat_ellipse(...), error = function(e) NULL)
# }
# pca <- ggplot(pca_df, aes(x = PC1, y = PC2, 
#                           color = Phenotype,
#                           label = Region)) +
#   geom_point(size = 4) +
#   geom_text_repel(size = 3.5, show.legend = FALSE) +
#   scale_color_manual(values = c("control" = "cornflowerblue", "AD" = "red4")) +
#   labs(
#     title = " ",
#     x = paste0("PC1 (", pc1, "% var)"),
#     y = paste0("PC2 (", pc2, "% var)")
#   ) +
#   theme_cowplot() +
#   #stat_ellipse(aes(group = Region), level = 0.95, geom = "polygon", alpha = 0.2) +
#   theme(
#     legend.position = "top",
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# #Vis
# pca
# 
# #Save plot
# ggsave("pca-globals-louvain.jpeg",
#        plot = pca,
#        device = "jpeg",
#        width = 5, 
#        height = 5,
#        units = "in",
#        dpi = 300
# )

############PCA of euclidean distance matrix ################

#Summary of global features
globals_summary <- globals %>%
  group_by(Region, Phenotype) %>%
  summarise(across(all_of(metric_cols),
                   mean, na.rm = TRUE),
            .groups = "drop") %>%
  unite("Network", Region, Phenotype, sep = " ")

#Euclidean distance matrix
dist_matrix <- globals_summary %>%
  column_to_rownames("Network") %>%
  dist(method = "euclidean") %>%
  as.matrix()

#Distance matrix clustering
hc_dist <- hclust(dist(dist_matrix), method = "ward.D2")
ord <- hc_dist$order

#How many clusters?
k <- 4
clusters <- cutree(hc_dist, k = k)
# calcular cambios de cluster
cluster_order <- clusters[ord]
#Distance matrix ordered
dist_reordered <- dist_matrix[ord, ord]

#Dendrogram of distance matrix
dendro <- ggdendrogram(hc_dist, rotate = TRUE, theme_dendro = FALSE) +
  theme_cowplot() +
  labs(title = "", x = "", y = "")

#Vis
dendro

#Hetamap of matrix distance 

#Get separation lines
breaks <- which(diff(cluster_order) != 0)

#Distance matrix in long format
dist_long <- as.data.frame(dist_reordered) %>%
  rownames_to_column("Network1") %>%
  pivot_longer(-Network1, names_to = "Network2", values_to = "Distance") %>%
  mutate(
    Network1 = factor(Network1, levels = rownames(dist_reordered)),
    Network2 = factor(Network2, levels = colnames(dist_reordered))
  ) %>%
  filter(as.integer(Network1) <= as.integer(Network2))

#Plot heatmap

heatmap_dist <- ggplot(dist_long, aes(Network1, Network2, fill = Distance)) +
  geom_tile() +
  scale_fill_viridis_c(
    guide = guide_colorbar(
      position = "top",
      barheight = unit(4, "pt"),    
      barwidth  = unit(120, "pt"),  
      title.position = "top"
    )
  ) +
  coord_fixed() +
  geom_hline(yintercept = breaks + 0.5, color = "black", size = 0.7) +
  geom_vline(xintercept = breaks + 0.5, color = "black", size = 0.7) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7)
  ) +
  labs(title = " ", x = "", y = "")

#Vis
heatmap_dist

#Convert the matrix distance into a pca object
pca_dist <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)

#Get scores and plot
pca_df <- as.data.frame(pca_dist$x) %>%
  rownames_to_column("RegionPhenotype") %>%
  separate(RegionPhenotype, into = c("Region", "Phenotype"), sep = " ")
var_exp <- (pca_dist$sdev^2 / sum(pca_dist$sdev^2)) * 100
pc1 <- round(var_exp[1], 1)
pc2 <- round(var_exp[2], 1)

#Plot
pca_dist.p <-ggplot(pca_df, aes(x = PC1, y = PC2, color = Phenotype, label = Region)) +
  geom_point(size = 3.8, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.3, point.padding = 0.2,
                  show.legend = FALSE) +
  scale_color_manual(
    values = c("control" = "cornflowerblue", "AD" = "red4"),
    name = "Group"
  ) +
  labs(
    title = "",
    x = paste0("PC1 (", pc1, "% var)"),
    y = paste0("PC2 (", pc2, "% var)")
  ) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    plot.title   = element_text(hjust = 0.5),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 9)
  )

##### FINAL GRID ####

library(cowplot)

# Panel A + B lado a lado (heatmap + dendrograma)
top_panel <- plot_grid(
  heatmap_dist,
  pca_dist.p,
  labels = c("A", "B"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(2, 1)   # A más grande que B
)

# Panel C (miniheatmap)
middle_panel <- plot_grid(
  dendro,
  labels = "C",
  label_size = 14,
  ncol = 1
)

# Panel D (PCA)
bottom_panel <- plot_grid(
  miniheat,
  labels = "D",
  label_size = 14,
  ncol = 1
)

# Figura final: (A|B) / C / D
final_plot <- plot_grid(
  top_panel,
  middle_panel,
  bottom_panel,
  ncol = 1,
  rel_heights = c(1.8, 1,   1.6)   # Ajusta como quieras
)

# Mostrar
final_plot

# Guardar
ggsave(
  "final_global_figure.jpeg",
  plot = final_plot,
  width = 14,
  height = 16,
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
# ggsave("2-global_metrics_lollipop.pdf", x, width = 12, height = 15, units = "in", dpi = 300)
# ggsave("3-global_metrics_lollipop.jpeg", x, width = 12, height = 18, units = "in", dpi = 300)

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
  rel_widths = c(1.5, 1)
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