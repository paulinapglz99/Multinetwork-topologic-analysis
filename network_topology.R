#!/usr/bin/env Rscript

#network_topology.R
#Usage:
#Rscript analyze_networks.R --input_dir data/edges --pattern "*.tsv" --out_dir results --workers 4 --per_node TRUE --make_html TRUE

if (!requireNamespace("pacman", quietly = FALSE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = FALSE)) install.packages("optparse", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "igraph",
  "data.table",
  "future.apply",
  "ggplot2",
  "tidyverse",
  "jsonlite",
  "stringr",
  "optparse", 
  "tools"
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages loaded correctly.: ",
       paste(names(ok)[!ok], collapse = ", "))
}

#Define option list for inputs

option_list <- list(
  make_option(c("-i","--input_dir"), type="character", default=NULL, help="Folder with files"),
  make_option(c("-p","--pattern"), type="character", default=".*\\.(txt|tsv|csv|graphml)$", help="Pattern regex for files"),
  make_option(c("-o","--out_dir"), type="character", default="results", help="Output directory"),
  make_option(c("-w","--workers"), type="integer", default=2, help="Number of parallel workers"),
  make_option(c("--per_node"), action="store_true", default=FALSE, help="Save metrics per node (CSV per network)"),
  make_option(c("--make_html"), action="store_true", default=FALSE, help="Generate aggregate HTML report (uses rmarkdown)"),
  make_option(c("--percol_steps"), type="integer", default=51, help="Number of steps in percolation simulation (e.g. 51 -> 0.2%, 4%, ...)"),
  make_option(c("--seed"), type="integer", default=42, help="Seed for reproducibility"),
  make_option(c("--type"), type="character", default="auto", help="Input format: 'auto' (default), 'edgelist' or 'adjacency'. Graphml is detected automatically. ")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Please add --input_dir")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files)==0) stop("No files matching the pattern were found in input_dir.")

#Helper function: read a network from CSV/TSV/GraphML
read_network <- function(path, type = opt$type) {
  message("DEBUG - type recibido: ", paste0(type, collapse = ","))
  
  type <- match.arg(type)  #From the opt parser
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    
    #EDGELIST
    if (type == "edgelist") {
      if (ncol(df) >= 3) {
        g <- igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
      } else {
        g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
      
      #ADJACENCY MATRIX
    } else if (type == "adjacency") {
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", 
                                               weighted = TRUE)
      
      #AUTO-DETECT
    } else if (type == "auto") {
      #If ALL columns (except the first) are numeric, I treat it as adjacency
      looks_like_adj <- all(sapply(df[,-1, drop=FALSE], is.numeric))
      
      if (looks_like_adj) {
        rownames(df) <- make.unique(as.character(df[, 1]))  
        df <- df[, -1, drop = FALSE]
        mat <- as.matrix(df)
        storage.mode(mat) <- "numeric"
        g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected",
                                                 weighted = TRUE)
        
      } else {
        #edge list: only the first 2 as nodes and optionally the 3rd as weight
        if (ncol(df) >= 3) {
          g <- igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
        } else {
          g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
        }
      }
    }
    
  } else {
    stop("Unsupported file format: ", ext)
  }
  
  #Final treatment, no bucles, loops
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"),
                        remove.multiple = TRUE, remove.loops = TRUE)
  
  return(g)
}
#sheesh

#Percolation threshold estimator (fast approximate)
#Strategy: for a series of removal fractions f in [0, fmax], remove fraction of nodes (random or targeted by degree)
#mode: ‘random’ or ‘targeted’. Defines how nodes are removed:
#   
#‘random’: nodes removed at random.
# 
#‘targeted’: nodes removed in descending order of degree (the most connected first).
#and compute giant component size fraction. Returns first f where GCC < threshold_frac (or NA).
#steps: number of node fraction values to test between 0 and 1 (default 51).
# 
#threshold_frac: fraction of the original graph size below which we consider the giant component to have collapsed (default 0.5, i.e., less than 50%).
# 
#trials: for "random" mode only: how many times to repeat the random removal to stabilize the result (default 3).

percolation_threshold <- function(g, mode = c("random","targeted"), steps = 51, threshold_frac = 0.5, trials = 3) {
  mode <- match.arg(mode)
  n <- vcount(g)
  if (n == 0) return(NA_real_)
  fracs <- seq(0, 1, length.out = steps)
  #For speed, limit to fmax 0.6 typically
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
      #remove highest-degree nodes cumulatively (fast: compute degrees once, remove top k)
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
    #choose median of trials
    return(median(res, na.rm = TRUE))
  } else {
    return(detect_f(g))
  }
}

#Function to analyze one file -> list with summary and optional node table path
analyze_one <- function(path) {
  nm <- basename(path)
  cat("Processing:", nm, "\n")
  g <- tryCatch(read_network(path), error = function(e) { message(e); return(NULL) })
  if (is.null(g)) return(NULL)
  g <- simplify(igraph::as_undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
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
  pr_norm <- (pr- min(pr)) / (max(pr) - min(pr))
  pr_top <- sort(pr, decreasing = TRUE)[1:min(5, length(pr))]  #Adjusting the degree distribution to estimate the gamma exponent
  
  #Degree distribution dataframe
    degs.df <- data.frame(degree = degs) %>%
    dplyr::count(degree, name = "freq") %>%
    dplyr::mutate(prob = freq / sum(freq)) %>%
    dplyr::filter(degree > 1, prob > 0)
  
  deg_dist <- data.table(
    network = nm,
    degree = degs.df$degree,
    Prob = degs.df$prob
  )
  
  if (nrow(degs.df) > 1) {
    fit <- lm(log(prob) ~ log(degree), data = degs.df)
    gamma_val <- -coef(fit)[["log(degree)"]]
  } else {
    gamma_val <- NA_real_
  }
  
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
    n_communities = length(comm_sizes),        
    largest_community_size = largest_comm_size, 
    largest_community_id = largest_comm_id,
    gamma_exponent = gamma_val  #gamma value
  )
  
  node_table_path <- NA_character_
  if (opt$per_node) {
    node_dt <- data.table(
      node = V(g)$name,
      degree = degs,
      pagerank = pr,
      pagerank_norm = pr_norm,
      kcore = kcore,
      membership_infomap = membership(comm_infomap)
    )
    nm_base <- tools::file_path_sans_ext(nm)
    node_table_path <- file.path(opt$out_dir, paste0(nm_base, "_nodes_summary.csv"))
    fwrite(node_dt, node_table_path)
  }
  
  #print(deg_dist)
  
  rm(g, g_giant, comp); gc(verbose = FALSE)
  return(list(summary = summary_row,
              node_table = node_table_path,
              degree_dist = deg_dist))
}


#Parallelize over files
plan(multisession, workers = opt$workers)
results <- future_lapply(files, FUN = analyze_one, future.seed = TRUE)

#Consolidate summaries
summaries <- rbindlist(lapply(results, function(x) if (!is.null(x)) x$summary else NULL), fill = TRUE)
summaries$network <- sub("^network_", "", tools::file_path_sans_ext(summaries$file))
fwrite(summaries, file.path(opt$out_dir, "networks_summary.csv"))

#Consolidate degree distributions
degree_dists <- rbindlist(lapply(results, function(x) if (!is.null(x)) x$degree_dist else NULL), fill = TRUE)
fwrite(degree_dists, file.path(opt$out_dir, "networks_degree_dist.csv"))


#Save metadata about per-node outputs
node_paths <- rbindlist(lapply(results, function(x) {
  if (!is.null(x) && !is.na(x$node_table)) data.table(file = x$summary$file, node_table = x$node_table) else NULL
}), fill = TRUE)
if (nrow(node_paths) > 0) fwrite(node_paths, file.path(opt$out_dir, "networks_nodes_index.csv"))

cat("Analysis completed. Summary saved in:", file.path(opt$out_dir, "networks_summary.csv"), "\n")

#Render an HTML report if requested

cat("Generating report \n")

if (opt$make_html) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown")
  }
  pacman::p_load(ggpubr, purrr, rmarkdown)
  
  rmd <- file.path(opt$out_dir, "network_report.Rmd")
  
  #R Markdown head
  cat(
    "---\n",
    "title: \"Network Analysis Report\"\n",
    "output: html_document\n",
    "---\n\n",
    file = rmd,  sep = ""
  )
  
  #Load data
  cat(
    "```{r, echo=FALSE}\n",
    "if (!requireNamespace('pacman', quietly = TRUE)) install.packages('pacman')\n",
    "pacman::p_load(data.table, ggplot2, ggpubr, purrr)\n",
    "summary <- fread('", file.path("networks_summary.csv"), "')\n",
    "print(summary)\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  #Plot function
  cat(
    "```{r, echo=FALSE}\n",
    "make_metric_plot <- function(df, metric, title, ylab, color) {\n",
    "  ggplot(df, aes(x = as.factor(network), y = .data[[metric]], group = 1)) +\n",
    "    geom_line() +\n",
    "    geom_point(color = color, size = 2.5) +\n",
    "    labs(title = title, x = 'Network', y = ylab) +\n",
    "    scale_y_continuous(limits = c(min(df[[metric]]), max(df[[metric]]))) +\n",
    "    theme_pubclean() +\n",
    "    theme(legend.position = 'none', axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))\n",
    "}\n",
    "```\n\n",
    file = rmd, append = TRUE,  sep = ""
  )
  
  #Degree distribution section
 
  #Load degree distributions
  cat(
    "```{r, echo=FALSE}\n",
    "deg_dists <- fread('", file.path("networks_degree_dist.csv"), "')\n",
    "deg_dists <- merge(deg_dists, summary[, .(file, gamma_exponent)], by.x = 'network', by.y = 'file', all.x = TRUE)\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  # Degree distribution plot function
  cat(
    "```{r, echo=FALSE, fig.width=6, fig.height=5}\n",
    "plot_degree <- function(df, net_name) {\n",
    "  gamma_val <- unique(df$gamma_exponent)\n",
    "  ggplot(df, aes(x = degree, y = Prob)) +\n",
    "    geom_point(color = 'blue', size = 2, alpha = 0.8) +\n",
    "    scale_x_log10() + scale_y_log10() +\n",
    "    labs(x = expression(log(k)), y = expression(log(p(k))),\n",
    "         subtitle = paste0(net_name, '\\nDegree distribution')) +\n",
    "    geom_smooth(method = 'lm', se = FALSE, color = 'black', linetype = 'dashed') +\n",
    "    annotate('text', x = min(df$degree), y = max(df$Prob), hjust = 0, vjust = 1,\n",
    "             label = paste0('γ = ', round(gamma_val, 2)), size = 4) +\n",
    "    theme_minimal(base_size = 14)\n",
    "}\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  cat(
    "```{r, echo=FALSE, results='asis'}\n",
    "nets <- unique(deg_dists$network)\n",
    "for (net in nets) {\n",
    "  df <- deg_dists[network == net]\n",
    "  print(plot_degree(df, net))\n",
    "}\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  cat(
    "```{r, echo=FALSE, fig.width=7, fig.height=6}\n",
    "ggplot(deg_dists, aes(x = degree, y = Prob, color = network)) +\n",
    "  geom_point(size = 2, alpha = 0.7) +\n",
    "  geom_smooth(method = 'lm', se = FALSE, linetype = 'dashed', size = 1) +\n",
    "  scale_x_log10() + scale_y_log10() +\n",
    "  labs(x = expression(log(k)), y = expression(log(p(k))),\n",
    "       title = 'All networks degree distribution') +\n",
    "  theme_minimal(base_size = 14)\n",
    "```\n\n",
    file = rmd, sep = '', append = TRUE
  )
  
    # Define plot specifications
  cat(
  "```{r, echo=FALSE}\n",
  "plot_specs <- data.frame(\n",
  "  var = c(\n",
  "    'n_nodes', 'n_edges', 'n_components', 'n_communities', 'deg_mean', 'avg_path_len',\n",
  "    'clustering_global', 'clustering_local_mean', 'diameter', 'global_density',\n",
  "    'size_giant_component', 'frac_giant_component', 'Q_modularity', 'perc_targeted_50'\n",
  "  ),\n",
  "  title = c(\n",
  "    'Number of vertices (genes)', 'Number of edges', 'Number of components',\n",
  "    'Number of communities (INFOMAP)', 'Mean degree of nodes', 'Average path length',\n",
  "    'Global clustering', 'Local clustering coefficient', 'Diameter', 'Global density',\n",
  "    'Size of the giant component', 'Fraction of giant component', 'Q modularity', 'Percolation threshold'\n",
  "  ),\n",
  "  ylab = c(\n",
  "    'Number of vertices', 'Number of edges', 'Number of components',\n",
  "    'Number of communities', 'Mean degree', 'Average path length',\n",
  "    'Global clustering', 'Local clustering', 'Diameter', 'Density',\n",
  "    'Size', 'Fraction', 'Q', 'Percolation threshold'\n",
  "  ),\n",
  "  color = c(\n",
  "    'cornflowerblue', 'mediumpurple', 'slateblue', 'maroon', 'navyblue', 'olivedrab',\n",
  "    'red4', 'lightpink4', 'red4', 'tomato4', 'lightsalmon', 'firebrick', 'seagreen', 'yellowgreen'\n",
  "  ),\n",
  "  stringsAsFactors = FALSE\n",
  ")\n",
  "```\n\n",
  file = rmd, append = TRUE, sep = ""
  )
  
  #Chunk for iterating plots
  cat(
    "```{r, echo=FALSE, fig.height=4, fig.width=7}\n",
    "walk2(plot_specs$var, seq_along(plot_specs$var), function(varname, idx) {\n",
    "  p <- make_metric_plot(summary, varname, plot_specs$title[idx], plot_specs$ylab[idx], plot_specs$color[idx])\n",
    "  print(p)\n",
    "})\n",
    "```\n\n",
    file = rmd, append = TRUE,  sep = ""
  )
  
#Renderize HTML
cat("Renderizing HTML \n")

rmarkdown::render(rmd, output_file = file.path("net_report.html"), quiet = TRUE)
cat("HTML report saved in:", file.path(opt$out_dir, "net_report.html"), "\n")
}

q(status = 0)

#END