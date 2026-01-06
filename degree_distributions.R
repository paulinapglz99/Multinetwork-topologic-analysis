#!/usr/bin/env Rscript

#degree_distributions.R
#Usage:
#Rscript degree_distributions.R --input_dir data/edges --pattern "*.tsv" --out_dir results

#Packages
if (!requireNamespace("pacman", quietly = F)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("optparse", quietly = F)) install.packages("optparse", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  "igraph",
  "data.table",
  "future.apply",
  "ggplot2",
  "tidyverse",
  "jsonlite",
  "stringr",
  "optparse", 
  "tools",
  "fitdistrplus",
  "VGAM",
  "poweRlaw"
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
  make_option(c("--seed"), type="integer", default=42, help="Seed for reproducibility"),
  make_option(c("--type"), type="character", default="auto", help="Input format: 'auto' (default), 'edgelist' or 'adjacency'. Graphml is detected automatically. "),
  make_option(c("--make_html"), action="store_true", default=FALSE, help="Generate aggregate HTML report (uses rmarkdown)")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Please add --input_dir")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files)==0) stop("No files matching the pattern were found in input_dir.")

#Helper function: read a network from CSV/TSV/GraphML
read_network <- function(path, type = opt$type) {
  message("DEBUG - type: ", paste0(type, collapse = ","))
  
  #type <- match.arg(type)  #From the opt parser
  type <- match.arg(type, choices = c("auto","edgelist","adjacency"))
  
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

#Function to fit distributions

fit_distributions <- function(degree_freq) {
  results <- list()
  
  #Power-Law
  pl_fit <- displ$new(degree_freq$degree)
  est_pl <- estimate_xmin(pl_fit)
  pl_fit$setXmin(est_pl)
  p_value_powl <- bootstrap_p(pl_fit, no_of_sims = 200, threads = 2) # reducimos simulaciones para velocidad
  
  results$PowerLaw <- list(p_value = p_value_powl$p)
  
  #Poisson
  poisson_fit <- fitdist(degree_freq$degree, "pois")
  lambda <- poisson_fit$estimate
  ks_test_pois <- ks.test(degree_freq$degree, "ppois", lambda)
  results$Poisson <- list(p_value = ks_test_pois$p.value)
  
  #Normal
  gaussian_fit <- fitdist(degree_freq$degree, "norm")
  mean_est <- gaussian_fit$estimate["mean"]
  sd_est <- gaussian_fit$estimate["sd"]
  ks_test_norm <- ks.test(degree_freq$degree, "pnorm", mean_est, sd_est)
  results$Gaussian <- list(p_value = ks_test_norm$p.value)
  
  #Pareto
  pareto_fit <- vglm(degree_freq$degree ~ 1, paretoff, trace = FALSE)
  shape <- Coef(pareto_fit)["shape"]
  scale <- Coef(pareto_fit)["scale"]
  ks_test_pareto <- tryCatch(
    ks.test(degree_freq$degree, "ppareto", scale = scale, shape = shape),
    error = function(e) list(p.value = NA)
  )
  results$Pareto <- list(p_value = ks_test_pareto$p.value)
  
  #Exponential
  exp_fit <- fitdist(degree_freq$degree, "exp", method="mme")
  rate_exp <- exp_fit$estimate["rate"]
  ks_test_exp <- ks.test(degree_freq$degree, "pexp", rate = rate_exp)
  results$Exponential <- list(p_value = ks_test_exp$p.value)
  
  # Negative binomial
  negbin_fit <- fitdist(degree_freq$degree, "nbinom")
  size_nb <- negbin_fit$estimate["size"]
  prob_nb <- negbin_fit$estimate["prob"]
  ks_test_nb <- ks.test(degree_freq$degree, "pnbinom", size = size_nb, prob = prob_nb)
  results$NegativeBinomial <- list(p_value = ks_test_nb$p.value)
  
  return(results)
}


#Function to analyze degree distributions
pretty_network_name <- function(path) {
  #only basename
  nm <- basename(path)
  #with no extension
  nm <- tools::file_path_sans_ext(nm)
  #and take - and _ our
  nm <- gsub("[_-]+", " ", nm)
  #colap multiple spaces
  nm <- gsub("\\s+", " ", nm)
  #remove spaces at the beginning and end
  nm <- trimws(nm)
  #Title Case
  nm <- tools::toTitleCase(nm)
  
  nm
}

analyze_degrees <- function(path) {
  nm <- basename(path)
  nm <- pretty_network_name(path)   # nombre "bonito"
  
  g <- tryCatch(read_network(path), error = function(e) { message(e); return(NULL) })
  if (is.null(g)) return(NULL)
  
  degs <- degree(g)
  df <- data.frame(degree = degs) %>%
    dplyr::count(degree, name = "freq") %>%
    dplyr::mutate(prob = freq / sum(freq)) %>%
    dplyr::filter(degree > 1, prob > 0)
  
  #Fit distributions
  fits <- fit_distributions(df)
  
  #Select best dist by p-value
  pvals <- sapply(fits, function(x) x$p_value)
  best_dist <- names(which.max(pvals))
  best_pval <- max(pvals, na.rm = TRUE)
  
  #Calculate gamma val
  gamma_val <- if (nrow(df) > 1) {
    fit <- lm(log(prob) ~ log(degree), data = df)
    -coef(fit)[["log(degree)"]]
  } else {
    NA_real_
  }
  
  list(
    deg_table = data.table(
      network = nm,
      degree = df$degree,
      Prob = df$prob,
      gamma_exponent = gamma_val
    ),
    summary = data.table(
      network = nm,
      best_distribution = best_dist,
      best_p_value = best_pval
    )
  )
}

#Parallelize over files

future::plan(future::multisession, workers = opt$workers)

all_deg <- future_lapply(files, analyze_degrees, future.seed = opt$seed)

deg_tables <- rbindlist(lapply(all_deg, `[[`, "deg_table"), fill = TRUE)
summaries <- rbindlist(lapply(all_deg, `[[`, "summary"), fill = TRUE)

fwrite(deg_tables, file.path(opt$out_dir, "networks_degree_dist.csv"))
fwrite(summaries, file.path(opt$out_dir, "best_fit_summary.csv"))

cat("Degree distributions saved in:", file.path(opt$out_dir, "networks_degree_dist.csv"), "\n")
cat("Best fit summary saved in:", file.path(opt$out_dir, "best_fit_summary.csv"), "\n")

#Render an HTML report if requested

cat("Generating report \n")

if (opt$make_html) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown")
  }
  pacman::p_load(ggpubr, purrr, rmarkdown)
  
  rmd <- file.path(opt$out_dir, "network_report.Rmd")
  
  #R Markdown head
  # cat(
  #   "---\n",
  #   "title: \"Network Degree Distribution Analysis Report\"\n",
  #   "output: html_document\n",
  #   "---\n\n",
  #   file = rmd,  sep = ""
  # )
  
  cat(
    "---\n",
    "title: \"Network Degree Distribution Analysis Report\"\n",
    "output: html_document\n",
    "---\n\n",
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(\n",
    "  fig.path = '',\n",
    "  dev = 'pdf',\n",
    "  dpi = 300\n",
    ")\n",
    "```\n\n",
    file = rmd, sep = ""
  )
  
  #Load data
  cat(
    "```{r, echo=FALSE}\n",
    "if (!requireNamespace('pacman', quietly = TRUE)) install.packages('pacman')\n",
    "pacman::p_load(data.table, ggplot2, ggpubr, purrr)\n",
    "summary <- fread('", file.path("networks_degree_dist.csv"), "')\n",
    "print(summary)\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  #Degree distribution section

  #Load degree distributions
  cat(
    "```{r, echo=FALSE}\n",
    "deg_dists <- fread('", file.path("networks_degree_dist.csv"), "')\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )

  #Degree distribution plot function
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
    "```{r per_network_degree_plots, echo=FALSE, results='asis'}\n",
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
  
  # Resumen de distribuciones
  cat(
    "## Best Distribution per Network\n\n",
    "```{r, echo=FALSE}\n",
    "summary_best <- fread('", file.path("best_fit_summary.csv"), "')\n",
    "knitr::kable(summary_best, caption = 'Best fit distribution for each network (based on KS test p-value)')\n",
    "```\n\n",
    file = rmd, sep = "", append = TRUE
  )
  
  
  #Renderize HTML
  cat("Renderizing HTML \n")
  
  #rmarkdown::render(rmd, output_file = file.path("degree_report.html"), quiet = TRUE)
  rmarkdown::render(
    rmd,
    output_file = "degree_report.html",
    output_dir = opt$out_dir,
    quiet = TRUE
  )
  
  cat("HTML report saved in:", file.path(opt$out_dir, "degree_report.html"), "\n")
}

#END