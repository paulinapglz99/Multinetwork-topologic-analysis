#!/usr/bin/env Rscript

#null_distributions.R
#Usage: https://github.com/paulinapglz99/Multinetwork-topologic-analysis
time<- Sys.time()
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
  make_option(c("--seed"), type="integer", default=42, help="Seed for reproducibility"),
  make_option(c("--type"), type="character", default="auto", help="Input format: 'auto', 'edgelist' or 'adjacency'. Graphml is detected automatically.")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Please add --input_dir")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files)==0) stop("No files matching the pattern were found in input_dir.")

#Function to read graphs
read_network <- function(path, type = opt$type) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    if (type == "edgelist") {
      if (ncol(df) >= 3) {
        g <- igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
      } else {
        g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
    } else if (type == "adjacency") {
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
    } else if (type == "auto") {
      looks_like_adj <- all(sapply(df[,-1, drop=FALSE], is.numeric))
      if (looks_like_adj) {
        rownames(df) <- make.unique(as.character(df[, 1]))  
        df <- df[, -1, drop = FALSE]
        mat <- as.matrix(df)
        storage.mode(mat) <- "numeric"
        g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      } else {
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
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"),
                        remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

#Calculate network metrics

calculate_metrics <- function(graph) {
  list(
    clustering = transitivity(graph, type = "global"),
    assortativity = assortativity_degree(graph),
    avg_path_length = mean_distance(graph, directed = FALSE),
    modularity = modularity(cluster_infomap(graph))
  )
}

#Randomize networks
randomize_graph <- function(graph, num_randomizations, niter = 1000) {
  replicate(num_randomizations, rewire(graph, with = keeping_degseq(niter = niter)), simplify = FALSE)
}
#For randomizations that need to maintain distribution
randomize_graph_nodis <- function(graph, num_randomizations, niter = 1000) {
  replicate(num_randomizations, rewire(graph, with = each_edge(0.8, loops = FALSE, multiple = TRUE)), simplify = FALSE)
}

#Calculate gamma
calculate_gamma <- function(graph) {
  degree_distribution <- degree(graph)
  degree_freq <- as.data.frame(table(degree_distribution))
  colnames(degree_freq) <- c("degree", "Freq")
  degree_freq$degree <- as.numeric(as.character(degree_freq$degree))
  degree_freq$Prob <- degree_freq$Freq / sum(degree_freq$Freq)
  degree_freq <- degree_freq[degree_freq$degree > 0, ]
  degree_freq$log_degree <- log(degree_freq$degree)
  degree_freq$log_Prob <- log(degree_freq$Prob)
  fit <- lm(log_Prob ~ log_degree, data = degree_freq)
  -coef(fit)["log_degree"]
}

#Function to analyze one network
null_dist_net <- function(path, num_randomizations = 500) {
  g <- read_network(path)
  metrics_obs <- calculate_metrics(g)
  gamma_obs <- calculate_gamma(g)
  
  random_graphs <- randomize_graph(g, num_randomizations)
  metrics_random <- lapply(random_graphs, calculate_metrics)
  
  clustering_random <- sapply(metrics_random, function(x) x$clustering)
  assort_random    <- sapply(metrics_random, function(x) x$assortativity)
  modu_random      <- sapply(metrics_random, function(x) x$modularity)
  
  random_graphs_nodis <- randomize_graph_nodis(g, num_randomizations)
  gamma_random <- sapply(random_graphs_nodis, calculate_gamma)
  
  results <- tibble::tibble(
    network = basename(path),
    clustering_obs = metrics_obs$clustering,
    clustering_ci_low = quantile(clustering_random, 0.025),
    clustering_ci_high = quantile(clustering_random, 0.975),
    clustering_outside_CI = !(metrics_obs$clustering %between% quantile(clustering_random, c(0.025, 0.975))),
    
    assort_obs = metrics_obs$assortativity,
    assort_ci_low = quantile(assort_random, 0.025),
    assort_ci_high = quantile(assort_random, 0.975),
    assort_outside_CI = !(metrics_obs$assortativity %between% quantile(assort_random, c(0.025, 0.975))),
    
    modularity_obs = metrics_obs$modularity,
    modularity_ci_low = quantile(modu_random, 0.025),
    modularity_ci_high = quantile(modu_random, 0.975),
    modularity_outside_CI = !(metrics_obs$modularity %between% quantile(modu_random, c(0.025, 0.975))),
    
    gamma_obs = gamma_obs,
    gamma_ci_low = quantile(gamma_random, 0.025),
    gamma_ci_high = quantile(gamma_random, 0.975),
    gamma_outside_CI = !(gamma_obs %between% quantile(gamma_random, c(0.025, 0.975)))
  )
  return(results)
}

#Execute

plan(multisession, workers = opt$workers)
summaries <- future_lapply(files, null_dist_net, future.seed = opt$seed)

#One tibble
summaries <- dplyr::bind_rows(summaries)

message("#Saving table")
output <- file.path(opt$out_dir, "null_model_report.csv")
summaries$network <- tools::file_path_sans_ext(summaries$network)
readr::write_csv(summaries, output)
message("Report is in: ", output)

#Render an HTML report if requested

cat("Generating report \n")
if (opt$make_html) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown")
  }
  pacman::p_load(ggpubr, purrr, rmarkdown)
  
  rmd <- file.path(opt$out_dir, "null_dis_network_report.Rmd")
  
  #R Markdown head
  cat(
    "---\n",
    "title: \"Network Null Distribution Analysis Report\"\n",
    "output: html_document\n",
    "---\n\n",
    file = rmd,  sep = ""
  )
  
  #Load data
  cat(
    "```{r, echo=FALSE}\n",
    "if (!requireNamespace('pacman', quietly = TRUE)) install.packages('pacman')\n",
    "pacman::p_load(data.table, ggplot2, ggpubr, purrr, dplyr, tidyr)\n",
    "summary_null <- fread('", file.path("null_model_report.csv"), "')\n",
    "knitr::kable(summary_null, caption = 'Network Null Distribution Analysis Report')\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  #Facet plot
  cat(
    "```{r facet_plot, echo=FALSE, fig.height=6, fig.width=10}\n",
    "summary_null.l <- summary_null %>%\n",
    "  select(network,\n",
    "         clustering_obs, clustering_ci_low, clustering_ci_high,\n",
    "         assort_obs, assort_ci_low, assort_ci_high,\n",
    "         modularity_obs, modularity_ci_low, modularity_ci_high,\n",
    "         gamma_obs, gamma_ci_low, gamma_ci_high) %>%\n",
    "  pivot_longer(-network,\n",
    "               names_to = c('metric', '.value'),\n",
    "               names_pattern = '(.*)_(obs|ci_low|ci_high)')\n\n",
    "ggplot(summary_null.l, aes(x = network, y = obs)) +\n",
    "  geom_point(color = 'red', size = 2) +\n",
    "  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +\n",
    "  coord_flip() +\n",
    "  facet_wrap(~metric, scales = 'free_y') +\n",
    "  theme_minimal() +\n",
    "  labs(y = 'Observed value with null CI',\n",
    "       x = 'Network',\n",
    "       title = 'Null model analysis across metrics')\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  #Heatmap
  cat(
    "```{r heatmap_plot, echo=FALSE, fig.height=5, fig.width=7}\n",
    "summary_sig <- summary_null %>%\n",
    "  select(network, ends_with('outside_CI')) %>%\n",
    "  pivot_longer(-network, names_to = 'metric', values_to = 'outside')\n\n",
    "ggplot(summary_sig, aes(x = metric, y = network, fill = outside)) +\n",
    "  geom_tile(color = 'white') +\n",
    "  scale_fill_manual(values = c('FALSE' = 'grey80', 'TRUE' = 'red')) +\n",
    "  labs(title = 'Metrics outside null CI', subtitle = 'Tiles in red show metrics outside the random distribution',\n",
    "       x = 'Metric', y = 'Network') +\n",
    "  theme_minimal()\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  #Renderize HTML
  cat("Renderizing HTML \n")
  rmarkdown::render(rmd, output_file = file.path("null_dis_network_report.html"), quiet = TRUE)
  cat("HTML report saved in:", file.path(opt$out_dir, "null_dis_network_report.html"), "\n")
}

difftime(Sys.time(), time, units = "mins")
#END