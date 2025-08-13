#!/usr/bin/env Rscript

# analyze_networks.R
#Packages needed: igraph, data.table, future.apply, ggplot2, rmarkdown (opcional), tidyverse
#Usage:
#Rscript analyze_networks.R --input_dir data/edges --pattern "*.tsv" --out_dir results --workers 4 --per_node TRUE --make_html TRUE

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  igraph,
  data.table,
  future.apply,
  ggplot2,
  jsonlite,
  stringr,
  optparse
)

#Define option list for inputs

option_list <- list(
  make_option(c("-i","--input_dir"), type="character", default=NULL, help="Folder with edge list files"),
  make_option(c("-p","--pattern"), type="character", default=".*\\.(txt|tsv|csv)$", help="Pattern regex for files"),
  make_option(c("-o","--out_dir"), type="character", default="results", help="Output directory"),
  make_option(c("-w","--workers"), type="integer", default=2, help="Number of parallel workers"),
  make_option(c("--per_node"), action="store_true", default=FALSE, help="Save metrics per node (CSV per network)"),
  make_option(c("--make_html"), action="store_true", default=FALSE, help="Generate aggregate HTML report (uses rmarkdown)"),
  make_option(c("--percol_steps"), type="integer", default=51, help="Number of steps in percolation simulation (e.g. 51 -> 0.2%, 4%, ...)"),
  make_option(c("--seed"), type="integer", default=42, help="Seed for reproducibility")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Please add --input_dir")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files)==0) stop("No files matching the pattern were found in input_dir.")

#Helper function: reads edge list (flexible to 2 or 3 columns: source, target[,weight])
read_edgelist <- function(path) {
  dt <- tryCatch(fread(path, header = FALSE, data.table = TRUE), error = function(e) NULL)
  if (is.null(dt)) stop(paste("Can't read", path))
  setnames(dt, names(dt), paste0("V", seq_len(ncol(dt))))
  if (ncol(dt) < 2) stop("The edge list must have at least 2 columns: source target")
  if (ncol(dt) >= 3) {
    g <- graph_from_data_frame(dt[, .(V1, V2, weight = V3)], directed = FALSE)
  } else {
    g <- graph_from_data_frame(dt[, .(V1, V2)], directed = FALSE)
  }
  return(g)
}

# Percolation threshold estimator (fast approximate)
# Strategy: for a series of removal fractions f in [0, fmax], remove fraction of nodes (random or targeted by degree)
# mode: ‘random’ or ‘targeted’. Defines how nodes are removed:
#   
# ‘random’: nodes removed at random.
# 
# ‘targeted’: nodes removed in descending order of degree (the most connected first).
# and compute giant component size fraction. Returns first f where GCC < threshold_frac (or NA).
# steps: number of node fraction values to test between 0 and 1 (default 51).
# 
# threshold_frac: fraction of the original graph size below which we consider the giant component to have collapsed (default 0.5, i.e., less than 50%).
# 
# trials: for "random" mode only: how many times to repeat the random removal to stabilize the result (default 3).

percolation_threshold <- function(g, mode = c("random","targeted"), steps = 51, threshold_frac = 0.5, trials = 3) {
  mode <- match.arg(mode)
  n <- vcount(g)
  if (n == 0) return(NA_real_)
  fracs <- seq(0, 1, length.out = steps)
  # For speed, limit to fmax 0.6 typically
  fracs <- fracs[fracs <= 0.6]
  detect_f <- function(g) {
    if (mode == "random") {
      for (f in fracs) {
        k <- floor(f * n)
        if (k == 0) next
        rem <- sample(V(g), k)
        sg <- delete_vertices(g, rem)
        gcc <- max(components(sg)$csize) / max(1, vcount(sg))
        if (is.na(gcc)) gcc <- 0
        if (gcc < threshold_frac) return(f)
      }
    } else if (mode == "targeted") {
      # remove highest-degree nodes cumulatively (fast: compute degrees once, remove top k)
      degs <- degree(g)
      ord <- order(degs, decreasing = TRUE)
      for (f in fracs) {
        k <- floor(f * n)
        if (k == 0) next
        rem <- V(g)[ord[seq_len(k)]]
        sg <- delete_vertices(g, rem)
        gcc <- if (vcount(sg)>0) max(components(sg)$csize) / vcount(g) else 0
        if (gcc < threshold_frac) return(f)
      }
    }
    return(NA_real_)
  }
  if (mode == "random") {
    res <- replicate(trials, detect_f(g))
    # choose median of trials
    return(median(res, na.rm = TRUE))
  } else {
    return(detect_f(g))
  }
}

#Function to analyze one file -> list with summary and optional node table path
analyze_one <- function(path) {
  nm <- basename(path)
  cat("Processing:", nm, "\n")
  g <- tryCatch(read_edgelist(path), error = function(e) { message(e); return(NULL) })
  if (is.null(g)) return(NULL)
  g <- simplify(as.undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  n <- vcount(g)
  m <- ecount(g)
  comp <- components(g)
  ncomp <- comp$no
  gcc_size <- max(comp$csize)
  gcc_frac <- gcc_size / max(1, n)
  giant_vids <- V(g)[which(comp$membership == which.max(comp$csize))]
  g_giant <- induced_subgraph(g, giant_vids)
  apl <- if (vcount(g_giant) > 1) mean_distance(g_giant, directed = FALSE, unconnected = FALSE) else NA
  diam <- if (vcount(g_giant) > 1) diameter(g_giant, directed = FALSE, weights = NA) else 0
  dens <- edge_density(g, loops = FALSE)
  assort <- tryCatch(assortativity_degree(g, directed = FALSE), error = function(e) NA_real_)
  clustering_local_mean <- transitivity(g, type = "localaverage", isolates = "zero")
  clustering_global <- transitivity(g, type = "global")
  degs <- degree(g)
  deg_mean <- mean(degs)
  deg_med <- median(degs)
  deg_sd <- sd(degs)
  kcore <- coreness(g)
  kcore_max <- max(kcore)
  pr <- page_rank(g)$vector
  pr_top <- sort(pr, decreasing = TRUE)[1:min(5, length(pr))]
  
  #Detecting communities with Infomap (error handling)
  comm_infomap <- tryCatch(
    cluster_infomap(g),
    error = function(e) make_clusters(g, membership = rep(1, n))
  )
  Q_mod <- modularity(comm_infomap)
  
  #Calculate largest community
  comm_sizes <- sizes(comm_infomap)
  largest_comm_size <- max(comm_sizes)
  largest_comm_id <- which.max(comm_sizes)
  
  #Percolation thresholds
  perc_targeted <- percolation_threshold(g, mode = "targeted", steps = opt$percol_steps, threshold_frac = 0.5, trials = 1)
  perc_random <- percolation_threshold(g, mode = "random", steps = opt$percol_steps, threshold_frac = 0.5, trials = 3)
  
  summary_row <- data.table(
    file = nm,
    n_nodes = n,
    n_edges = m,
    avg_path_len = apl,
    diameter = diam,
    global_density = dens,
    size_giant_component = gcc_size,
    frac_giant_component = gcc_frac,
    n_components = ncomp,
    clustering_local_mean = clustering_local_mean,
    clustering_global = clustering_global,
    assortativity = assort,
    deg_mean = deg_mean,
    deg_median = deg_med,
    deg_sd = deg_sd,
    kcore_max = kcore_max,
    Q_modularity = Q_mod,
    perc_targeted_50 = perc_targeted,
    perc_random_50 = perc_random,
    n_communities = length(comm_sizes),        #total communities detected
    largest_community_size = largest_comm_size, #largest community size
    largest_community_id = largest_comm_id     #ID of the largest community
  )
  
  node_table_path <- NA_character_
  if (opt$per_node) {
    node_dt <- data.table(
      node = V(g)$name,
      degree = degs,
      pagerank = pr,
      kcore = kcore,
      membership_infomap = membership(comm_infomap)
    )
    nm_base <- tools::file_path_sans_ext(nm)
    node_table_path <- file.path(opt$out_dir, paste0(nm_base, "_nodes_summary.csv"))
    fwrite(node_dt, node_table_path)
  }
  rm(g, g_giant, comp); gc(verbose = FALSE)
  return(list(summary = summary_row, node_table = node_table_path))
}

#Parallelize over files
plan(multisession, workers = opt$workers)
results <- future_lapply(files, FUN = analyze_one, future.seed = TRUE)

#Consolidate summaries
summaries <- rbindlist(lapply(results, function(x) if (!is.null(x)) x$summary else NULL), fill = TRUE)
summaries$network <- sub("^network_", "", tools::file_path_sans_ext(summaries$file))
fwrite(summaries, file.path(opt$out_dir, "networks_summary.csv"))

#Save metadata about per-node outputs
node_paths <- rbindlist(lapply(results, function(x) {
  if (!is.null(x) && !is.na(x$node_table)) data.table(file = x$summary$file, node_table = x$node_table) else NULL
}), fill = TRUE)
if (nrow(node_paths) > 0) fwrite(node_paths, file.path(opt$out_dir, "networks_nodes_index.csv"))

cat("Analysis completed. Summary saved in:", file.path(opt$out_dir, "networks_summary.csv"), "\n")

#Optionally render an HTML report (simple) if requested
if (opt$make_html) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("rmarkdown not installed; install to generate HTML")
  } else {
    rmd <- file.path(opt$out_dir, "net_report.Rmd")
    cat("---
title: "Network analysis report"
output: html_document
---

```{r, echo=FALSE}
library(data.table); library(ggplot2); library(ggpubr)
summary <- fread('/datos/home/paulinapg/Multinetwork-classic-analysis/results/networks_summary.csv')
print(head(summary))
  
```
```{r, echo=FALSE}
#Plot number of vertices 

v_plot <- summary %>%
  ggplot(aes(x = as.factor(network), y = n_nodes, group = 1)) +
    geom_line() +
    geom_point(color = "cornflowerblue", size = 2.5)+
  labs(title = "Number of vertices (genes)",
       x = "Networks",
       y = "Number of vertices") +
    scale_y_continuous(limits = c(max(summary$n_nodes)/4, max(summary$n_nodes))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  v_plot
```

```{r, echo=FALSE}
#Plot number of edges 

E_plot <- summary %>%
  ggplot(aes(x = as.factor(network), y = n_edges, group = 1)) +
    geom_line() +
    geom_point(color = "mediumpurple", size = 2.5)+
  labs(title = "Number of edges",
       x = "Networks",
       y = "Number of edges") +
    scale_y_continuous(limits = c(max(summary$n_edges)/4, max(summary$n_edges))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
E_plot

```

```{r, echo=FALSE}
#Plot number of components

components_plot <- summary %>%
  ggplot(aes(x = as.factor(network), y = n_components, group = 1)) +
    geom_line() +
    geom_point(color = "slateblue", size = 2.5)+
  labs(title = "Number of components",
       x = "Networks",
       y = "Number of components") +
    scale_y_continuous(limits = c(max(summary$n_components)/4, max(summary$n_components))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
components_plot
```

```{r, echo=FALSE}
#Plot number of communities with INFOMAP

components_INFOMAP_plot <- summary %>%
  ggplot(aes(x = as.factor(network), y = n_communities, group = 1)) +
    geom_line() +
    geom_point(color = "navyblue", size = 2.5)+
  labs(title = "Number of communities (INFOMAP)",
       x = "Network",
       y = "Number of communities") +
    scale_y_continuous(limits = c(max(summary$n_communities)/4, max(summary$n_communities))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
components_INFOMAP_plot
```

```{r, echo=FALSE}
#Plot average path length

av_path.p <- summary %>%
  ggplot(aes(x = as.factor(network), y = avg_path_len, group = 1)) +
    geom_line() +
    geom_point(color = "olivedrab", size = 2.5)+
  labs(title = "Average path length",
       x = "Network",
       y = "Average path length") +
    scale_y_continuous(limits = c(max(summary$avg_path_len)/4, max(summary$avg_path_len))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
av_path.p
```

```{r, echo=FALSE}
#Plot global clustering coefficient

global_clus.p <- summary %>%
  ggplot(aes(x = as.factor(network), y = clustering_global, group = 1)) +
    geom_line() +
    geom_point(color = "red4", size = 2.5)+
  labs(title = "Global clustering",
       x = "Network",
       y = "Global clustering") +
    scale_y_continuous(limits = c(max(summary$clustering_global)/4, max(summary$clustering_global))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
global_clus.p
```

```{r, echo=FALSE}
#Plot local clustering coefficient

local_clus.p <- summary %>%
  ggplot(aes(x = as.factor(network), y = clustering_local_mean, group = 1)) +
    geom_line() +
    geom_point(color = "lightpink4", size = 2.5)+
  labs(title = "Global clustering coefficient",
       x = "Network",
       y = "Global clustering") +
    scale_y_continuous(limits = c(max(summary$clustering_local_mean)/4, max(summary$clustering_local_mean))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
local_clus.p
```

```{r, echo=FALSE}
#Plot Diameter

diameter.p <-  summary %>%
  ggplot(aes(x = as.factor(network), y = diameter, group = 1)) +
    geom_line() +
    geom_point(color = "red4", size = 2.5)+
  labs(title = "Diameter",
       x = "Network",
       y = "Diameter") +
    scale_y_continuous(limits = c(max(summary$diameter)/4, max(summary$diameter))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
diameter.p
```

```{r, echo=FALSE}
#Plot Global density

density.p <-  summary %>%
  ggplot(aes(x = as.factor(network), y = global_density, group = 1)) +
    geom_line() +
    geom_point(color = "lightsalmon", size = 2.5)+
  labs(title = "Global density",
       x = "Network",
       y = "Density") +
    scale_y_continuous(limits = c(max(summary$global_density)/4, max(summary$global_density))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
density.p
```

```{r, echo=FALSE}
#Plot size of giant component

size_giant.p <- summary %>%
  ggplot(aes(x = as.factor(network), y = size_giant_component, group = 1)) +
    geom_line() +
    geom_point(color = "firebrick", size = 2.5)+
  labs(title = "Size of the giant component",
       x = "Network",
       y = "Size") +
    scale_y_continuous(limits = c(max(summary$size_giant_component)/4, max(summary$size_giant_component))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
size_giant.p
```

```{r, echo=FALSE}
#Plot percentage of genes in the max size comp

Q_modularity.p <-summary %>%
  ggplot(aes(x = as.factor(network), y = Q_modularity, group = 1)) +
    geom_line() +
    geom_point(color = "seagreen", size = 2.5)+
  labs(title = "Q modularity",
       x = "Network",
       y = "Q") +
    scale_y_continuous(limits = c(max(summary$Q_modularity)/4, max(summary$Q_modularity))) +
  theme_pubclean() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
Q_modularity.p
```

", file = rmd, append = TRUE)
    rmarkdown::render(rmd, output_file = file.path(opt$out_dir, "net_report.html"), quiet = TRUE)
    cat("HTML report generated in:", file.path(opt$out_dir, "net_report.html"), "\n")
  }
}

q(status = 0)

#END