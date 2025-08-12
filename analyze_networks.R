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

# Function to analyze one file -> list with summary and optional node table path
analyze_one <- function(path) {
  nm <- basename(path)
  cat("Processing:", nm, "\n")
  g <- tryCatch(read_edgelist(path), error = function(e) { message(e); return(NULL) })
  if (is.null(g)) return(NULL)
  # Ensure simple undirected (remove loops/multiedges)
  g <- simplify(as.undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  n <- vcount(g)
  m <- ecount(g)
  comp <- components(g)
  ncomp <- comp$no
  gcc_size <- max(comp$csize)
  gcc_frac <- gcc_size / max(1, n)
  # If giant component has >0 nodes, restrict distances to it for path metrics
  giant_vids <- V(g)[which(comp$membership == which.max(comp$csize))]
  g_giant <- induced_subgraph(g, giant_vids)
  # Average path length & diameter: only on giant component (to avoid Inf)
  apl <- if (vcount(g_giant) > 1) mean_distance(g_giant, directed = FALSE, unconnected = FALSE) else NA
  diam <- if (vcount(g_giant) > 1) diameter(g_giant, directed = FALSE, weights = NA) else 0
  dens <- edge_density(g, loops = FALSE)
  assort <- tryCatch(assortativity_degree(g, directed = FALSE), error = function(e) NA_real_)
  # clustering: local mean and global transitivity (triangle-based)
  clustering_local_mean <- transitivity(g, type = "localaverage", isolates = "zero")  # average local
  clustering_global <- transitivity(g, type = "global")
  # degree distribution summary
  degs <- degree(g)
  deg_mean <- mean(degs)
  deg_med <- median(degs)
  deg_sd <- sd(degs)
  # K-core index (coreness)
  kcore <- coreness(g)
  kcore_max <- max(kcore)
  # PageRank
  pr <- page_rank(g)$vector
  pr_top <- sort(pr, decreasing = TRUE)[1:min(5, length(pr))]
  # Q modularity using Infomap partition
  comm_infomap <- tryCatch(cluster_infomap(g), error = function(e) make_clusters(g, membership = rep(1, vcount(g))))
  Q_mod <- modularity(comm_infomap)
  ncomponents <- ncomp
  # Percolation thresholds (targeted + random)
  perc_targeted <- percolation_threshold(g, mode = "targeted", steps = opt$percol_steps, threshold_frac = 0.5, trials = 1)
  perc_random <- percolation_threshold(g, mode = "random", steps = opt$percol_steps, threshold_frac = 0.5, trials = 3)
  # Compose summary
  summary_row <- data.table(
    file = nm,
    n_nodes = n,
    n_edges = m,
    avg_path_len = apl,
    diameter = diam,
    global_density = dens,
    size_giant_component = gcc_size,
    frac_giant_component = gcc_frac,
    n_components = ncomponents,
    clustering_local_mean = clustering_local_mean,
    clustering_global = clustering_global,
    assortativity = assort,
    deg_mean = deg_mean,
    deg_median = deg_med,
    deg_sd = deg_sd,
    kcore_max = kcore_max,
    Q_modularity = Q_mod,
    perc_targeted_50 = perc_targeted,
    perc_random_50 = perc_random
  )
  # Optionally save per-node metrics (degree, pagerank, kcore, cluster membership)
  node_table_path <- NA_character_
  if (opt$per_node) {
    node_dt <- data.table(
      node = V(g)$name,
      degree = degs,
      pagerank = pr,
      kcore = kcore,
      membership_infomap = membership(comm_infomap)
    )
    node_table_path <- file.path(opt$out_dir, paste0(nm, "_nodes.csv"))
    fwrite(node_dt, node_table_path)
  }
  # free memory
  rm(g, g_giant, comp); gc(verbose = FALSE)
  return(list(summary = summary_row, node_table = node_table_path))
}

# Paraleliza sobre archivos
plan(multisession, workers = opt$workers)
results <- future_lapply(files, FUN = analyze_one, future.seed = TRUE)

# Consolidate summaries
summaries <- rbindlist(lapply(results, function(x) if (!is.null(x)) x$summary else NULL), fill = TRUE)
fwrite(summaries, file.path(opt$out_dir, "networks_summary.csv"))

# Save metadata about per-node outputs
node_paths <- rbindlist(lapply(results, function(x) {
  if (!is.null(x) && !is.na(x$node_table)) data.table(file = x$summary$file, node_table = x$node_table) else NULL
}), fill = TRUE)
if (nrow(node_paths) > 0) fwrite(node_paths, file.path(opt$out_dir, "networks_nodes_index.csv"))

cat("Análisis finalizado. Resumen guardado en:", file.path(opt$out_dir, "networks_summary.csv"), "\n")

# Optionally render an HTML report (simple) if requested
if (opt$make_html) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("rmarkdown not installed; install to generate HTML")
  } else {
    rmd <- file.path(opt$out_dir, "net_report.Rmd")
    cat("---
title: \"Network analysis report\"
output: html_document
---\n\n", file = rmd)
    cat("```{r, echo=FALSE}\nlibrary(data.table); library(ggplot2)\nsummary <- fread('", file.path(opt$out_dir, "networks_summary.csv"), "')\nprint(head(summary))\n# ejemplo: histograma de n_nodes\np <- ggplot(summary, aes(x = n_nodes)) + geom_histogram(bins=40) + theme_minimal()\nprint(p)\n```\n", file = rmd, append = TRUE)
    rmarkdown::render(rmd, output_file = file.path(opt$out_dir, "net_report.html"), quiet = TRUE)
    cat("HTML report generated in:", file.path(opt$out_dir, "net_report.html"), "\n")
  }
}

q(status = 0)
