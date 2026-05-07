#!/usr/bin/env Rscript
# =============================================================================
# edge_threshold_sensitivity.R
# Sensitivity analysis of network topological properties across multiple
# edge-count thresholds (top-N by MI score).
#
# Purpose: Demonstrate that key findings (regional clustering hierarchy,
# contraction-densification pattern, CRB stability) are robust to the
# choice of edge threshold.
#
# Usage:
#   Rscript edge_threshold_sensitivity.R \
#     --input_dir /path/to/raw_edgelists \
#     --out_dir results/sensitivity \
#     --workers 8
#
# Input filenames expected:
#   <cohort>_<region>_counts_<phenotype>.tsv
#   e.g. Mayo_CRB_counts_AD.tsv, ROSMAP_DLPFC_counts_control.tsv
#
# Output: Cut networks, degree distributions, topological summaries, HTML report
# =============================================================================

# --- 0. Packages -------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

pacman::p_load(
  data.table,
  igraph,
  future.apply,
  ggplot2,
  dplyr,
  tidyr,
  stringr,
  purrr,
  scales,
  ggpubr,
  ggsci,
  optparse,
  tools,
  fitdistrplus,
  VGAM,
  poweRlaw,
  rmarkdown,
  knitr,
  ggrepel
)

# --- 1. CLI options ----------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory containing raw ARACNe edgelist TSVs"),
  make_option(c("-o", "--out_dir"), type = "character", default = "results/sensitivity",
              help = "Output directory [default: %default]"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$",
              help = "Regex pattern to match input files [default: %default]"),
  make_option(c("-w", "--workers"), type = "integer", default = 4,
              help = "Parallel workers [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default: %default]"),
  make_option(c("--thresholds"), type = "character",
              default = "10000,50000,100000,200000,300000,500000",
              help = "Comma-separated edge thresholds [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (is.null(opt$input_dir)) stop("Specify --input_dir")
if (!dir.exists(opt$input_dir)) stop("Input directory not found: ", opt$input_dir)

THRESHOLDS <- as.integer(str_split_1(opt$thresholds, ","))
message("Edge thresholds: ", paste(THRESHOLDS, collapse = ", "))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$out_dir, "networks"), showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "degree_dists"), showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "topology"), showWarnings = FALSE)

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) == 0) stop("No files matching pattern found in ", opt$input_dir)
message("Found ", length(files), " edgelist files.")

# --- 2. Helper: parse network metadata from filename -------------------------
#
# Works on BOTH original and cut filenames:
#   Mayo_CRB_counts_AD.tsv          -> region=CRB, phenotype=AD
#   Mayo_CRB_counts_AD_top10000.tsv -> region=CRB, phenotype=AD
#
# Strategy: strip the _topN suffix first, then parse the clean name.

parse_network_meta <- function(filename) {
  nm <- tools::file_path_sans_ext(basename(filename))
  
  # Strip _topNNNNN suffix if present (from cut files)
  clean_nm <- str_replace(nm, "_top\\d+$", "")
  
  parts <- str_split_1(clean_nm, "_")
  
  # Phenotype is the last token after stripping _topN
  phenotype <- parts[length(parts)]
  
  # Region: match known codes
  regions <- c("CRB", "DLPFC", "HCN", "PCC", "TC")
  region <- intersect(toupper(parts), regions)
  region <- if (length(region) > 0) region[1] else "UNKNOWN"
  
  # Cohort is the first token
  cohort <- parts[1]
  
  list(
    name      = clean_nm,
    full_name = nm,
    region    = region,
    phenotype = phenotype,
    cohort    = cohort
  )
}

# --- 3. Step 1: Cut networks at each threshold -------------------------------

message("\n=== STEP 1: Cutting networks at multiple thresholds ===")

cut_one_network <- function(file_path, threshold) {
  meta <- parse_network_meta(file_path)
  out_file <- file.path(
    opt$out_dir, "networks",
    paste0(meta$name, "_top", threshold, ".tsv")
  )
  
  if (file.exists(out_file)) {
    message("  [skip] ", basename(out_file), " already exists")
    return(out_file)
  }
  
  df <- fread(file_path, header = FALSE, col.names = c("gene1", "gene2", "MI"))
  df <- df[order(-abs(MI))][seq_len(min(threshold, nrow(df)))]
  fwrite(df, out_file, sep = "\t", col.names = FALSE)
  message("  [done] ", basename(out_file), " (", nrow(df), " edges)")
  return(out_file)
}

task_grid <- expand.grid(
  file = files,
  threshold = THRESHOLDS,
  stringsAsFactors = FALSE
)

plan(multisession, workers = opt$workers)

cut_results <- future_mapply(
  cut_one_network,
  file_path = task_grid$file,
  threshold = task_grid$threshold,
  SIMPLIFY = FALSE,
  future.seed = opt$seed
)

message("Network cutting complete. ", length(cut_results), " files generated.")

# --- 4. Step 2: Compute topology ---------------------------------------------

message("\n=== STEP 2: Computing topological metrics ===")

compute_topology <- function(net_path) {
  
  meta <- parse_network_meta(net_path)
  thr <- as.integer(str_extract(basename(net_path), "(?<=_top)\\d+"))
  
  df <- fread(net_path, header = FALSE, col.names = c("gene1", "gene2", "MI"))
  g <- igraph::graph_from_data_frame(df[, .(gene1, gene2)], directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  n <- igraph::vcount(g)
  m <- igraph::ecount(g)
  comp <- igraph::components(g)
  gcc_size <- max(comp$csize)
  
  giant_ids <- which(comp$membership == which.max(comp$csize))
  g_giant <- igraph::induced_subgraph(g, giant_ids)
  
  apl <- if (igraph::vcount(g_giant) > 1) igraph::mean_distance(g_giant, directed = FALSE) else NA_real_
  diam <- if (igraph::vcount(g_giant) > 1) igraph::diameter(g_giant, directed = FALSE, weights = NA) else 0L
  
  degs <- igraph::degree(g)
  kcore <- max(igraph::coreness(g))
  
  comm <- tryCatch(
    igraph::cluster_louvain(g),
    error = function(e) igraph::make_clusters(g, membership = rep(1L, n))
  )
  comm_sizes <- igraph::sizes(comm)
  
  data.table(
    network     = meta$name,
    region      = meta$region,
    phenotype   = meta$phenotype,
    threshold   = thr,
    n_nodes     = n,
    n_edges     = m,
    density     = igraph::edge_density(g),
    gcc_size    = gcc_size,
    gcc_frac    = gcc_size / max(1L, n),
    n_components = comp$no,
    avg_path_len = apl,
    diameter    = diam,
    clustering_global = igraph::transitivity(g, type = "global"),
    clustering_local  = igraph::transitivity(g, type = "localaverage", isolates = "zero"),
    deg_mean    = mean(degs),
    deg_median  = median(degs),
    deg_sd      = sd(degs),
    kcore_max   = kcore,
    modularity  = igraph::modularity(comm),
    n_communities = length(comm_sizes),
    largest_community = max(comm_sizes)
  )
}

cut_files <- list.files(
  file.path(opt$out_dir, "networks"),
  pattern = "\\.tsv$", full.names = TRUE
)

topo_results <- future_lapply(cut_files, compute_topology, future.seed = opt$seed)
topo_dt <- rbindlist(topo_results, fill = TRUE)

fwrite(topo_dt, file.path(opt$out_dir, "topology", "sensitivity_topology.csv"))
message("Topological metrics saved: sensitivity_topology.csv")

# --- 5. Step 3: Degree distributions -----------------------------------------

message("\n=== STEP 3: Computing degree distributions ===")

compute_degree_dist <- function(net_path) {
  
  meta <- parse_network_meta(net_path)
  thr <- as.integer(str_extract(basename(net_path), "(?<=_top)\\d+"))
  
  df <- fread(net_path, header = FALSE, col.names = c("gene1", "gene2", "MI"))
  g <- igraph::graph_from_data_frame(df[, .(gene1, gene2)], directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  degs <- igraph::degree(g)
  deg_table <- data.table(degree = degs) %>%
    count(degree, name = "freq") %>%
    mutate(prob = freq / sum(freq)) %>%
    filter(degree > 1, prob > 0)
  
  gamma_val <- if (nrow(deg_table) > 1) {
    fit <- lm(log(prob) ~ log(degree), data = deg_table)
    -coef(fit)[["log(degree)"]]
  } else {
    NA_real_
  }
  
  best_dist <- tryCatch({
    pvals <- c()
    pf <- fitdist(deg_table$degree, "pois")
    pvals["Poisson"] <- ks.test(deg_table$degree, "ppois", pf$estimate)$p.value
    gf <- fitdist(deg_table$degree, "norm")
    pvals["Gaussian"] <- ks.test(deg_table$degree, "pnorm",
                                 gf$estimate["mean"], gf$estimate["sd"])$p.value
    ef <- fitdist(deg_table$degree, "exp", method = "mme")
    pvals["Exponential"] <- ks.test(deg_table$degree, "pexp",
                                    rate = ef$estimate["rate"])$p.value
    nf <- fitdist(deg_table$degree, "nbinom")
    pvals["NegBinomial"] <- ks.test(deg_table$degree, "pnbinom",
                                    size = nf$estimate["size"],
                                    prob = nf$estimate["prob"])$p.value
    names(which.max(pvals))
  }, error = function(e) NA_character_)
  
  data.table(
    network   = meta$name,
    region    = meta$region,
    phenotype = meta$phenotype,
    threshold = thr,
    gamma     = gamma_val,
    best_fit  = best_dist
  )
}

deg_results <- future_lapply(cut_files, compute_degree_dist, future.seed = opt$seed)
deg_dt <- rbindlist(deg_results, fill = TRUE)

fwrite(deg_dt, file.path(opt$out_dir, "degree_dists", "sensitivity_degree_summary.csv"))
message("Degree distribution summary saved.")

# --- 6. Step 4: Generate Rmarkdown report ------------------------------------

message("\n=== STEP 4: Generating sensitivity report ===")

rmd_path <- file.path(opt$out_dir, "sensitivity_report.Rmd")

rmd_content <- '
---
title: "Edge Threshold Sensitivity Analysis"
subtitle: "Robustness of network topological properties across pruning thresholds"
author: "Automated pipeline"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    code_folding: hide
    fig_width: 10
    fig_height: 7
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE, warning = FALSE, message = FALSE,
  fig.width = 10, fig.height = 7, dpi = 150
)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, ggplot2, dplyr, tidyr, ggsci, ggpubr,
               scales, knitr, ggrepel, RColorBrewer)
```

```{r load_data}
topo <- fread("topology/sensitivity_topology.csv")
deg  <- fread("degree_dists/sensitivity_degree_summary.csv")

# Create labels
topo <- topo %>%
  mutate(
    net_label = paste(region, phenotype, sep = " "),
    threshold_k = threshold / 1000
  )
deg <- deg %>%
  mutate(
    net_label = paste(region, phenotype, sep = " "),
    threshold_k = threshold / 1000
  )

# Custom 10-color palette (5 regions x 2 phenotypes)
net_levels <- sort(unique(topo$net_label))
n_nets <- length(net_levels)
net_palette <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(n_nets),
  net_levels
)
```

# 1. Overview

This report evaluates the sensitivity of global network topological properties to
the number of edges retained during network pruning. Networks were pruned at
thresholds of **`r paste(sort(unique(topo$threshold_k)), collapse = ", ")`** thousand
edges (top-N by mutual information score). If the key findings reported in the
manuscript are robust, the **relative ordering** of networks across metrics should
remain stable across thresholds, even as absolute values change.

```{r summary_table}
kable(
  topo %>%
    group_by(region, phenotype) %>%
    summarise(thresholds_tested = n(), .groups = "drop") %>%
    arrange(region, phenotype),
  caption = "Networks and thresholds analysed"
)
```

# 2. Global metrics across thresholds

```{r helper_plot}
# Reusable plot function for line plots across thresholds
plot_metric <- function(df, y_var, y_lab, title_text, palette) {
  ggplot(df, aes(x = threshold_k, y = .data[[y_var]],
                 color = net_label, group = net_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = palette) +
    labs(x = "Edge threshold (thousands)", y = y_lab,
         title = title_text, color = "Network") +
    theme_pubclean(base_size = 13) +
    theme(legend.position = "right")
}
```

## 2.1 Number of nodes

```{r nodes_plot}
plot_metric(topo, "n_nodes", "Number of nodes",
            "Node count across edge thresholds", net_palette)
```

## 2.2 Network density

```{r density_plot}
plot_metric(topo, "density", "Density",
            "Network density across edge thresholds", net_palette)
```

## 2.3 Mean degree

```{r mean_degree}
plot_metric(topo, "deg_mean", "Mean degree",
            "Mean degree across edge thresholds", net_palette)
```

## 2.4 Giant component fraction

```{r gcc_frac}
plot_metric(topo, "gcc_frac", "Giant component fraction",
            "Giant component fraction across edge thresholds", net_palette) +
  scale_y_continuous(labels = percent_format())
```

## 2.5 Average path length

```{r apl}
plot_metric(topo, "avg_path_len", "Average path length",
            "Average path length across edge thresholds", net_palette)
```

## 2.6 Diameter

```{r diameter_plot}
plot_metric(topo, "diameter", "Diameter",
            "Diameter across edge thresholds", net_palette)
```

## 2.7 Clustering coefficient (global)

```{r clustering}
plot_metric(topo, "clustering_global", "Global clustering coefficient",
            "Global clustering across edge thresholds", net_palette)
```

## 2.8 Modularity (Louvain Q)

```{r modularity}
plot_metric(topo, "modularity", "Modularity (Q)",
            "Louvain modularity across edge thresholds", net_palette)
```

## 2.9 Number of communities

```{r n_communities}
plot_metric(topo, "n_communities", "Number of communities",
            "Community count across edge thresholds", net_palette)
```

## 2.10 Max k-core

```{r kcore}
plot_metric(topo, "kcore_max", "Maximum k-core",
            "Maximum k-core across edge thresholds", net_palette)
```

# 3. Key robustness tests

## 3.1 Regional clustering hierarchy stability (PCA)

For each threshold, we compute a PCA of all 17 topological metrics and plot
networks in PC1-PC2 space. If the regional grouping persists across panels, the
hierarchy is not an artifact of the 200K cutoff.

```{r pca_facet, fig.height=14, fig.width=14}
# Define region palette (5 regions)
region_palette <- setNames(
  pal_aaas()(5),
  c("CRB", "DLPFC", "HCN", "PCC", "TC")
)

pca_list <- topo %>%
  group_by(threshold_k) %>%
  group_split() %>%
  purrr::map_dfr(function(df) {
    metric_cols <- c("n_nodes", "n_edges", "density", "gcc_size", "gcc_frac",
                     "n_components", "avg_path_len", "diameter",
                     "clustering_global", "clustering_local",
                     "deg_mean", "deg_median", "deg_sd", "kcore_max",
                     "modularity", "n_communities", "largest_community")
    mat <- df %>% dplyr::select(all_of(metric_cols)) %>% as.matrix()
    keep <- apply(mat, 2, var, na.rm = TRUE) > 0
    mat <- mat[, keep, drop = FALSE]
    mat[is.na(mat)] <- 0
    pca <- prcomp(mat, center = TRUE, scale. = TRUE)
    var_exp <- summary(pca)$importance[2, 1:2] * 100
    data.frame(
      net_label   = df$net_label,
      region      = df$region,
      phenotype   = df$phenotype,
      threshold_k = df$threshold_k[1],
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      var_PC1 = var_exp[1],
      var_PC2 = var_exp[2]
    )
  })

ggplot(pca_list, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = net_label), size = 2.8,
                           max.overlaps = 20, show.legend = FALSE) +
  facet_wrap(~ paste0(threshold_k, "K edges"), scales = "free", ncol = 3) +
  scale_color_manual(values = region_palette) +
  labs(title = "PCA of global network metrics across edge thresholds",
       subtitle = "Regional grouping should remain stable across panels",
       color = "Region") +
  theme_pubclean(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))
```

## 3.2 AD vs Control difference heatmap stability

For each threshold, we compute the signed log-transformed Z-score of AD minus
control metric means per region. Consistent colour patterns across thresholds
confirm robustness of the contraction-densification finding.

```{r diff_heatmaps, fig.height=14, fig.width=12}
metric_cols <- c("n_nodes", "n_edges", "density", "gcc_size", "gcc_frac",
                 "n_components", "avg_path_len", "diameter",
                 "clustering_global", "clustering_local",
                 "deg_mean", "deg_median", "kcore_max",
                 "modularity", "n_communities", "largest_community")

# Compute difference matrices for all thresholds at once
diff_all <- topo %>%
  group_by(threshold_k, region, phenotype) %>%
  summarise(across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(threshold_k, region),
    names_from = phenotype,
    values_from = all_of(metric_cols)
  )

# Compute raw differences (AD - control) for each metric
diff_long <- purrr::map_dfr(metric_cols, function(mc) {
  ad_col <- paste0(mc, "_AD")
  ct_col <- paste0(mc, "_control")
  if (!all(c(ad_col, ct_col) %in% names(diff_all))) return(NULL)
  diff_all %>%
    transmute(
      threshold_k = threshold_k,
      region = region,
      metric = mc,
      raw_diff = .data[[ad_col]] - .data[[ct_col]]
    )
})

# Z-score within each threshold-metric combination, then signed log transform
diff_long <- diff_long %>%
  group_by(threshold_k, metric) %>%
  mutate(
    z_diff = as.numeric(scale(raw_diff)),
    z_diff = ifelse(is.na(z_diff), 0, z_diff),
    log_diff = sign(z_diff) * log10(abs(z_diff) + 1)
  ) %>%
  ungroup()

ggplot(diff_long, aes(x = metric, y = region, fill = log_diff)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ paste0(threshold_k, "K edges"), ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "Signed\nlog(|z|+1)"
  ) +
  labs(title = "AD - Control topological differences across edge thresholds",
       subtitle = "Consistent patterns across panels indicate robustness",
       x = "Metric", y = "Region") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank()
  )
```

## 3.3 Rank correlation of metrics across thresholds

For each metric, we rank all 10 networks at each threshold, then compute the
Spearman correlation of that ranking against the ranking at the reference threshold
(200K). If $\\rho > 0.9$ across all thresholds, the relative ordering of networks is
preserved and the 200K cutoff is not driving the results.

```{r rank_corr, fig.height=7, fig.width=11}
ref_thr <- 200

rank_corrs <- topo %>%
  dplyr::select(net_label, threshold_k, all_of(metric_cols)) %>%
  pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  group_by(metric, threshold_k) %>%
  mutate(rank_val = rank(value, na.last = "keep")) %>%
  ungroup()

ref_ranks <- rank_corrs %>%
  filter(threshold_k == ref_thr) %>%
  dplyr::select(net_label, metric, ref_rank = rank_val)

corr_summary <- rank_corrs %>%
  inner_join(ref_ranks, by = c("net_label", "metric")) %>%
  group_by(metric, threshold_k) %>%
  summarise(
    spearman_rho = cor(rank_val, ref_rank, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

ggplot(corr_summary, aes(x = threshold_k, y = spearman_rho,
                          color = metric, group = metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = ref_thr, linetype = "dotted", color = "gray70") +
  annotate("text", x = ref_thr + 10, y = 0.05, label = "200K\n(reference)",
           hjust = 0, size = 3, color = "gray50") +
  scale_y_continuous(limits = c(-0.2, 1.05)) +
  scale_x_continuous(breaks = sort(unique(corr_summary$threshold_k))) +
  labs(x = "Edge threshold (K)",
       y = expression("Spearman " * rho * " vs 200K reference"),
       title = "Rank stability of network metrics across edge thresholds",
       subtitle = "Each line = one metric. Dashed line = rho 0.9. At 200K, rho = 1 by definition.",
       color = "Metric") +
  theme_pubclean(base_size = 12) +
  theme(legend.position = "right", legend.text = element_text(size = 7))
```

# 4. Degree distribution summary

```{r deg_table}
kable(
  deg %>%
    dplyr::select(network, region, phenotype, threshold_k, gamma, best_fit) %>%
    arrange(region, phenotype, threshold_k),
  caption = "Degree distribution: gamma exponent and best-fit model across thresholds",
  digits = 3
)
```

```{r gamma_plot, fig.height=6}
ggplot(deg, aes(x = threshold_k, y = gamma, color = net_label, group = net_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = net_palette) +
  labs(x = "Edge threshold (thousands)",
       y = expression(gamma * " exponent"),
       title = "Degree distribution gamma exponent across thresholds",
       color = "Network") +
  theme_pubclean(base_size = 13) +
  theme(legend.position = "right")
```

```{r bestfit_table}
bestfit_wide <- deg %>%
  dplyr::select(region, phenotype, threshold_k, best_fit) %>%
  pivot_wider(names_from = threshold_k, values_from = best_fit,
              names_prefix = "thr_") %>%
  arrange(region, phenotype)

kable(bestfit_wide,
      caption = "Best-fit degree distribution model per network across thresholds")
```

# 5. Interpretation guide

**What each section answers for the reviewer:**

1. **Section 2 (line plots):** If lines run roughly parallel without crossing,
   the relative differences between networks are preserved regardless of the
   edge threshold. Vertical spread at each threshold shows the magnitude of
   inter-network differences; what matters is that the *ordering* is stable.

2. **Section 3.1 (PCA facets):** If regional clusters (same colour) stay
   together across panels, the finding that networks group by region rather
   than disease status is robust.

3. **Section 3.2 (heatmaps):** If the colour pattern (which regions show
   increases vs decreases in which metrics) is consistent across panels, the
   contraction-densification pattern and CRB expansion pattern are not artifacts.

4. **Section 3.3 (rank correlations):** Spearman rho = 1 at 200K by
   definition (reference). If rho stays above 0.9 at other thresholds,
   the 200K cutoff produces the same relative ordering as any other cutoff.
   Metrics that drop below 0.9 at extreme thresholds (10K) identify the
   least robust topological features.

5. **Section 4 (degree distributions):** If all networks are Gaussian at all
   thresholds, it means the degree distributions are unimodal and symmetric
   rather than heavy-tailed. This is expected for pruned MI networks because
   retaining only the top-N edges compresses the upper tail of the degree
   distribution --- the pruning itself imposes a ceiling on maximum degree.
   The key question is whether the *gamma exponent* is stable, not whether
   the best-fit model changes.

```{r session_info}
sessionInfo()
```
'

writeLines(rmd_content, rmd_path)
message("Rmarkdown template written to: ", rmd_path)

message("Rendering HTML report...")
tryCatch({
  rmarkdown::render(
    rmd_path,
    output_file = "sensitivity_report.html",
    output_dir = opt$out_dir,
    quiet = FALSE
  )
  message("HTML report saved: ", file.path(opt$out_dir, "sensitivity_report.html"))
}, error = function(e) {
  message("Report rendering failed: ", e$message)
  message("You can render manually with:")
  message("  rmarkdown::render('", rmd_path, "')")
})

message("\n=== Pipeline complete ===")
message("Output directory: ", opt$out_dir)