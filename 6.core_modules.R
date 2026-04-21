#!/usr/bin/env Rscript

# core_modules.R
# Mesoscopic analysis: module classification, conservation, and functional enrichment
# across AD and Control co-expression networks in five brain regions.

# ==============================================================================
# 0. SETUP
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  "data.table", "vroom", "tidyverse", "tools", "optparse",
  "ggplot2", "igraph", "cowplot",
  "clusterProfiler", "org.Hs.eg.db",
  "ggsci",        # for scale_fill_aaas() / scale_color_aaas()
  "ggpubr"        # for stat_compare_means() annotation
)

# --- Paths (hardcoded for development) ---
input_dir  <- "~/Desktop/local_work/fomo_networks/results_topos_louvain"
output_dir <- "~/Desktop/local_work/fomo_networks/results_core_modules"
enrich_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/results_comm/"
pattern    <- "_nodes_summary\\.csv$"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# --- Thresholds ---
upp_thres    <- 0.6
low_thres    <- 0.2
min_genes    <- 3

# --- Region order ---
ordered_regions <- c("HCN", "DLPFC", "PCC", "CRB", "TC")

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

# Extract region and phenotype from filename
extract_metadata <- function(filename) {
  base <- basename(filename)
  if (!grepl("_counts_", base)) return(NULL)
  data.frame(
    filename  = filename,
    region    = sub("_counts_.*", "", base) |> gsub(pattern = "^(Mayo_|ROSMAP_)", replacement = ""),
    phenotype = ifelse(grepl("_AD_", base), "AD", "Control")
  )
}

# Load node/membership table and add a unique module ID
load_modules <- function(file) {
  df <- vroom::vroom(file, show_col_types = FALSE)
  stopifnot(all(c("node", "membership") %in% names(df)))
  
  region    <- sub("_counts_.*", "", basename(file)) |> gsub(pattern = "^(Mayo_|ROSMAP_)", replacement = "")
  phenotype <- ifelse(grepl("_AD_", file), "AD", "Control")
  
  df$unique_module_id <- paste(region, phenotype, df$membership, sep = "_")
  split(df$node, df$unique_module_id)
}

# Jaccard similarity between two gene sets
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0  || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Build a long table of module sizes (one row per module) from the flat modules list
build_module_size_table <- function(modules_list, meta_df) {
  tibble(
    Module  = names(modules_list),
    N_genes = lengths(modules_list)
  ) %>%
    dplyr::mutate(
      Phenotype = dplyr::case_when(
        grepl("_AD_",      Module) ~ "AD",
        grepl("_Control_", Module) ~ "Control"
      ),
      Region = stringr::str_extract(Module, paste(c("HCN", "PCC", "TC", "DLPFC", "CRB"), collapse = "|"))
    ) %>%
    dplyr::filter(!is.na(Region), !is.na(Phenotype))
}

# Shapiro-Wilk normality test per group (Region x Phenotype)
test_normality <- function(size_df) {
  size_df %>%
    dplyr::group_by(Region, Phenotype) %>%
    dplyr::summarise(
      n       = dplyr::n(),
      sw_stat = ifelse(n >= 3 & n <= 5000, shapiro.test(N_genes)$statistic, NA_real_),
      sw_pval = ifelse(n >= 3 & n <= 5000, shapiro.test(N_genes)$p.value,   NA_real_),
      normal  = ifelse(!is.na(sw_pval), sw_pval > 0.05, NA),
      .groups = "drop"
    )
}

# Comparison 1: AD vs Control within each region
# Chooses t-test (both normal) or Wilcoxon (otherwise) per region
test_ad_vs_control <- function(size_df, normality_df) {
  purrr::map_dfr(unique(size_df$Region), function(reg) {
    norm_reg    <- dplyr::filter(normality_df, Region == reg)
    both_normal <- all(norm_reg$normal, na.rm = TRUE)
    
    ad_vals   <- dplyr::filter(size_df, Region == reg, Phenotype == "AD")$N_genes
    ctrl_vals <- dplyr::filter(size_df, Region == reg, Phenotype == "Control")$N_genes
    
    result <- if (both_normal) t.test(ad_vals, ctrl_vals) else wilcox.test(ad_vals, ctrl_vals)
    
    tibble(
      Region    = reg,
      test      = ifelse(both_normal, "t-test", "Wilcoxon"),
      statistic = as.numeric(result$statistic),
      p_value   = result$p.value,
      p_label   = dplyr::case_when(
        result$p.value < 0.001 ~ "***",
        result$p.value < 0.01  ~ "**",
        result$p.value < 0.05  ~ "*",
        TRUE                   ~ "ns"
      )
    )
  })
}

# Comparison 2: across regions within each phenotype
# Kruskal-Wallis global test + Dunn post-hoc with BH correction
test_across_regions <- function(size_df) {
  purrr::map_dfr(c("AD", "Control"), function(pheno) {
    sub <- dplyr::filter(size_df, Phenotype == pheno)
    
    kw <- kruskal.test(N_genes ~ Region, data = sub)
    
    # Dunn post-hoc (requires FSA or rstatix; use rstatix here)
    dunn <- rstatix::dunn_test(sub, N_genes ~ Region, p.adjust.method = "BH") %>%
      dplyr::mutate(
        Phenotype = pheno,
        p_label   = dplyr::case_when(
          p.adj < 0.001 ~ "***",
          p.adj < 0.01  ~ "**",
          p.adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
    
    # Attach KW global result as attribute summary row
    message("Kruskal-Wallis [", pheno, "]: chi2 = ", round(kw$statistic, 3),
            ", df = ", kw$parameter, ", p = ", signif(kw$p.value, 3))
    
    dunn
  })
}

# NMI between AD and Control partitions for a single region
compute_nmi <- function(region_name, meta_df, min_genes = 3) {
  sub <- dplyr::filter(meta_df, region == region_name)
  
  if (!all(c("AD", "Control") %in% sub$phenotype)) {
    message("Skipping ", region_name, ": missing one phenotype.")
    return(NULL)
  }
  
  ad   <- data.table::fread(sub$filename[sub$phenotype == "AD"],   select = c("node", "membership"))
  ctrl <- data.table::fread(sub$filename[sub$phenotype == "Control"], select = c("node", "membership"))
  
  # Filter small modules
  ad   <- dplyr::filter(ad,   membership %in% names(which(table(membership) >= min_genes)))
  ctrl <- dplyr::filter(ctrl, membership %in% names(which(table(membership) >= min_genes)))
  
  merged <- dplyr::inner_join(ad, ctrl, by = "node", suffix = c("_AD", "_CTRL"))
  
  if (nrow(merged) == 0) {
    return(tibble(Region = region_name, n_genes_common = 0L, NMI = NA_real_))
  }
  
  tibble(
    Region         = region_name,
    n_genes_common = nrow(merged),
    NMI            = igraph::compare(merged$membership_AD, merged$membership_CTRL, method = "nmi")
  )
}

# GO:BP enrichment with simplification (ENSEMBL input)
run_go_enrichment <- function(genes, ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2, simplify_cutoff = 0.7) {
  if (length(genes) < 3) return(NULL)
  
  ego <- clusterProfiler::enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  clusterProfiler::simplify(ego, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min)
}

# Build a bipartite igraph (GO terms <-> genes) from an enrichResult and export GraphML
build_cnet_graph <- function(enrich_res, pair_id, output_dir = NULL) {
  if (is.null(enrich_res) || nrow(enrich_res@result) == 0) return(NULL)
  
  res <- as.data.frame(enrich_res)
  
  edges <- res |>
    dplyr::select(ID, geneID) |>
    tidyr::separate_rows(geneID, sep = "/") |>
    dplyr::rename(from = ID, to = geneID) |>
    dplyr::mutate(Pair_ID = pair_id)
  
  nodes_go <- res |>
    dplyr::select(ID, Description, pvalue, p.adjust, qvalue, GeneRatio, Count) |>
    dplyr::distinct() |>
    dplyr::transmute(
      name = ID, label = Description, type = "GO",
      Pair_ID = pair_id,
      pvalue, p_adjust = p.adjust, qvalue, GeneRatio, Count
    )
  
  nodes_gene <- edges |>
    dplyr::select(to) |>
    dplyr::distinct() |>
    dplyr::transmute(
      name = to, label = to, type = "Gene",
      Pair_ID = pair_id,
      pvalue = NA_real_, p_adjust = NA_real_, qvalue = NA_real_,
      GeneRatio = NA_character_, Count = NA_integer_
    )
  
  nodes <- dplyr::bind_rows(nodes_go, nodes_gene)
  g     <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    igraph::write_graph(g, file.path(output_dir, paste0(pair_id, ".graphml")), format = "graphml")
  }
  g
}

# Identify modules whose maximum Jaccard against the opposite phenotype is below upp_thres
find_exclusive_modules <- function(jaccards_tb, phenotype_col, module_col) {
  jaccards_tb |>
    dplyr::group_by(Region, .data[[module_col]]) |>
    dplyr::summarise(has_similar = any(Classification == "Similar"), .groups = "drop") |>
    dplyr::filter(!has_similar) |>
    dplyr::select(Region, Module = .data[[module_col]]) |>
    dplyr::mutate(Phenotype = phenotype_col)
}

# Cross-regional conservation: find exclusive module pairs from different regions
# with Jaccard >= threshold
find_conserved_cross_regional <- function(exclusive_df, modules_list, threshold = 0.6) {
  purrr::map_dfr(unique(exclusive_df$Phenotype), function(pheno) {
    sub <- dplyr::filter(exclusive_df, Phenotype == pheno)
    pairs <- combn(nrow(sub), 2, simplify = FALSE)
    
    purrr::map_dfr(pairs, function(idx) {
      i <- idx[1]; j <- idx[2]
      if (sub$Region[i] == sub$Region[j]) return(NULL)
      
      jval <- jaccard_index(modules_list[[sub$Module[i]]], modules_list[[sub$Module[j]]])
      if (jval < threshold) return(NULL)
      
      data.frame(
        Module_1 = sub$Module[i], Region_1 = sub$Region[i],
        Module_2 = sub$Module[j], Region_2 = sub$Region[j],
        Phenotype = pheno, Jaccard_Index = jval
      )
    })
  })
}

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

files       <- list.files(input_dir, pattern = pattern, full.names = TRUE)
meta        <- purrr::map_dfr(files, extract_metadata)
module_list <- purrr::map(files, load_modules)

# Keep only modules with >= min_genes genes
module_list <- purrr::map(module_list, ~ .x[lengths(.x) >= min_genes])
modules     <- unlist(module_list, recursive = FALSE)   # flat named list

regions <- unique(meta$region)

# Load and annotate community enrichments
enrich_files <- list.files(enrich_dir, pattern = "_enrichment\\.csv$", full.names = TRUE)

enrich_all <- purrr::map_dfr(enrich_files, function(f) {
  fname <- basename(f)
  data.table::fread(f) |>
    dplyr::mutate(
      Region    = sub("_counts_.*", "", fname) |> gsub(pattern = "^(Mayo_|ROSMAP_)", replacement = ""),
      Phenotype = ifelse(grepl("_AD_", fname), "AD", "Control"),
      CommunityID = as.character(CommunityID)
    )
}) |>
  dplyr::mutate(Module = paste(Region, Phenotype, CommunityID, sep = "_"))

# ==============================================================================
# 3. PART 0 — NMI: GLOBAL PARTITION SIMILARITY PER REGION
# ==============================================================================

nmi_results <- purrr::map_dfr(regions, compute_nmi, meta_df = meta, min_genes = min_genes)
print(nmi_results)

# vroom::vroom_write(nmi_results, file.path(output_dir, "NMI_AD_vs_Control_by_region.csv"))

# ==============================================================================
# 4. PART 1 — JACCARD CLASSIFICATION OF MODULE PAIRS
# ==============================================================================

# Count modules per network
module_counts <- tibble(
  Region    = meta$region,
  Phenotype = meta$phenotype,
  N_modules = lengths(module_list),
  Network   = paste(Region, Phenotype, sep = "_")
) |>
  dplyr::group_by(Region) |>
  dplyr::mutate(Total_region_modules = sum(N_modules)) |>
  dplyr::ungroup()

module_gene_counts <- tibble(
  Module  = names(modules),
  N_genes = lengths(modules)
)

# Compute all pairwise Jaccards within each region
jaccards_tb <- purrr::map_dfr(regions, function(reg) {
  ctrl_mods <- modules[grep(paste0(reg, "_Control"), names(modules))]
  ad_mods   <- modules[grep(paste0(reg, "_AD"),      names(modules))]
  
  purrr::map_dfr(names(ctrl_mods), function(cn) {
    purrr::map_dfr(names(ad_mods), function(an) {
      jval <- jaccard_index(ctrl_mods[[cn]], ad_mods[[an]])
      tibble(
        Region         = reg,
        Module_Control = cn,
        Module_AD      = an,
        Jaccard_Index  = jval,
        Classification = dplyr::case_when(
          jval >= upp_thres ~ "Similar",
          jval <= low_thres ~ "Dissimilar",
          TRUE              ~ "Intermediate"
        )
      )
    })
  })
}) |>
  dplyr::mutate(Region = factor(Region, levels = ordered_regions))

print(table(jaccards_tb$Classification))
# vroom::vroom_write(jaccards_tb, file.path(output_dir, "jaccards_all_regions.csv"))

# Classify each module as Similar / Intermediate / Exclusive
classify_modules <- function(jaccards, module_col, label) {
  jaccards |>
    dplyr::group_by(Region, .data[[module_col]]) |>
    dplyr::summarise(
      has_similar      = any(Classification == "Similar"),
      has_intermediate = any(Classification == "Intermediate"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Classification = dplyr::case_when(
        has_similar      ~ "Similar",
        has_intermediate ~ "Intermediate",
        TRUE             ~ paste0(label, "_exclusive")
      )
    ) |>
    dplyr::rename(Module = .data[[module_col]]) |>
    dplyr::mutate(Phenotype = label)
}

classified_all <- dplyr::bind_rows(
  classify_modules(jaccards_tb, "Module_AD",      "AD"),
  classify_modules(jaccards_tb, "Module_Control", "Control")
)

# Summary proportions
# module_counts already carries Region + Phenotype, so join only by Network
# and select only the columns we need to avoid Region.x/y, Phenotype.x/y duplicates
summary_all <- classified_all |>
  dplyr::group_by(Region, Phenotype, Classification) |>
  dplyr::summarise(n_modules = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(Network = paste(Region, Phenotype, sep = "_")) |>
  dplyr::left_join(
    dplyr::select(module_counts, Network, N_modules, Total_region_modules),
    by = "Network"
  ) |>
  dplyr::mutate(
    local_proportion  = round(n_modules / N_modules * 100, 1),
    global_proportion = round(n_modules / Total_region_modules * 100, 1)
  ) %>% 
  dplyr::mutate(Region = factor(Region, levels = ordered_regions))

# vroom::vroom_write(summary_all, file.path(output_dir, "modules_classification_summary.csv"))

# Exclusive modules (Jaccard == 0 with all opposite-phenotype modules)
exclusive_AD   <- find_exclusive_modules(jaccards_tb, "AD",      "Module_AD")
exclusive_CTRL <- find_exclusive_modules(jaccards_tb, "Control", "Module_Control")

exclusive_modules <- dplyr::bind_rows(exclusive_AD, exclusive_CTRL) |>
  dplyr::left_join(module_gene_counts, by = "Module") |>
  dplyr::arrange(Region, Phenotype, desc(N_genes))

# vroom::vroom_write(exclusive_modules, file.path(output_dir, "exclusive_modules_AD_CTRL.csv"))

# ==============================================================================
# 5. PART 2 — CROSS-REGIONAL CONSERVATION OF EXCLUSIVE MODULES
# ==============================================================================

conserved_pairs <- find_conserved_cross_regional(exclusive_modules, modules, threshold = upp_thres) |>
  dplyr::mutate(Pair_ID = paste(Module_1, Module_2, sep = "__"))

print(conserved_pairs)
# vroom::vroom_write(conserved_pairs, file.path(output_dir, "conserved_exclusive_modules_across_regions.csv"))

# Gene set comparisons for each conserved pair
pair_gene_sets <- conserved_pairs |>
  dplyr::rowwise() |>
  dplyr::mutate(
    Genes_1            = list(modules[[Module_1]]),
    Genes_2            = list(modules[[Module_2]]),
    Genes_intersection = list(intersect(Genes_1, Genes_2)),
    Genes_unique_1     = list(setdiff(Genes_1, Genes_2)),
    Genes_unique_2     = list(setdiff(Genes_2, Genes_1)),
    Genes_union        = list(union(Genes_1, Genes_2)),
    N_genes_1          = length(Genes_1),
    N_genes_2          = length(Genes_2),
    N_intersection     = length(Genes_intersection),
    N_unique_1         = length(Genes_unique_1),
    N_unique_2         = length(Genes_unique_2),
    N_union            = length(Genes_union)
  ) |>
  dplyr::ungroup()

# GO enrichment on each intersection (only pairs with sufficient overlap)
pair_intersections <- dplyr::filter(pair_gene_sets, N_intersection >= min_genes)

go_by_pair <- purrr::set_names(
  purrr::map(pair_intersections$Genes_intersection, run_go_enrichment),
  pair_intersections$Pair_ID
)

# Build and export cnet GraphML files
cnet_networks <- purrr::imap(go_by_pair, ~ build_cnet_graph(.x, .y, output_dir = "graphml_cnet_networks"))

# ==============================================================================
# 6. PLOTS
# ==============================================================================

# --- Color palette (AAAS) and custom classification colors ---
aaas_cols <- ggsci::pal_aaas()(5)

fill_colors <- c(
  "Control_exclusive" = aaas_cols[1],
  "Intermediate"      = aaas_cols[3],
  "Similar"           = aaas_cols[4],
  "AD_exclusive"      = aaas_cols[2]
)

# 6A. Module size probability density per region — ridgeline plot
#
# Each ridge is one network (Region x Phenotype). The x-axis shows module size
# (N genes) and the y-axis encodes the probability density, so curves integrate
# to 1 regardless of how many modules each network has.
#
# Statistical tests:
#   - AD vs Control within each region: Shapiro-Wilk → t-test or Wilcoxon
#   - Across regions within each phenotype: Kruskal-Wallis + Dunn post-hoc (BH)

pacman::p_load("ggridges", "rstatix")   # load here to keep setup section clean

module_sizes <- build_module_size_table(modules, meta) %>%
  dplyr::mutate(
    Region    = factor(Region,    levels = ordered_regions),
    Phenotype = factor(Phenotype, levels = c("Control", "AD"))
  )

# --- Run tests ---
normality_results  <- test_normality(module_sizes)
stat_ad_vs_ctrl    <- test_ad_vs_control(module_sizes, normality_results)
stat_across_regions <- test_across_regions(module_sizes)

print(normality_results)
print(stat_ad_vs_ctrl)
print(stat_across_regions)

# vroom::vroom_write(normality_results,   file.path(output_dir, "module_size_normality.csv"))
# vroom::vroom_write(stat_ad_vs_ctrl,     file.path(output_dir, "module_size_ADvsCtrl.csv"))
# vroom::vroom_write(stat_across_regions, file.path(output_dir, "module_size_across_regions.csv"))

# --- Ridgeline plot ---
# One panel per Phenotype; ridges ordered by ordered_regions (top to bottom)
# Median line added per ridge for reference

# Both phenotypes overlap on the same ridge per region.
# alpha < 1 lets the two density curves show through each other.
# The p-value label (AD vs Control) sits to the right of each Region label.

# Attach significance labels to the size table for annotation
module_sizes_annot <- module_sizes %>%
  dplyr::left_join(
    dplyr::select(stat_ad_vs_ctrl, Region, p_label),
    by = "Region"
  ) %>%
  dplyr::mutate(
    # Build a combined y-axis label: region name + significance star
    Region_label = paste0(Region, "  ", p_label),
    Region_label = factor(Region_label,
                          levels = unique(Region_label[order(match(Region, ordered_regions))]))
  )

p_module_ridges <- ggplot(module_sizes_annot,
                          aes(x = N_genes, y = Region_label, fill = Phenotype, color = Phenotype)) +
  ggridges::geom_density_ridges(
    scale          = 1.4,          # overlap between ridges
    quantile_lines = TRUE,
    quantiles      = 2,            # median line
    alpha          = 0.5,
    linewidth      = 0.5,
    rel_min_height = 0.005         # trim very low-density tails
  ) +
  scale_fill_aaas(alpha = 0.6) +
  scale_color_aaas() +
  scale_x_continuous(
    trans  = "log10",
    labels = scales::label_comma()
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(face = "bold")
  ) +
  labs(
    x = "Module size (N genes, log scale)",
    y = ""
  )

# 6B. NMI heatmap
p_nmi <- nmi_results |>
  dplyr::arrange(NMI) |>
  dplyr::mutate(Region = factor(Region, levels = ordered_regions)) |>
  ggplot(aes(x = "AD vs Control", y = Region, fill = NMI)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(NMI, 3)), size = 4) +
  scale_fill_gradient(low = aaas_cols[1], high = "#8EAACC") +
  labs(title = "NMI", x = "", y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid    = element_blank(),
    legend.position = "none"
  )

# 6C. Mean Jaccard heatmap
p_mean_jaccard <- jaccards_tb |>
  dplyr::group_by(Region) |>
  dplyr::summarise(mean_jaccard = mean(Jaccard_Index, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(mean_jaccard) |>
  dplyr::mutate(Region = factor(Region, levels = ordered_regions)) |>
  ggplot(aes(x = "", y = Region, fill = mean_jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_jaccard, 3)), size = 4) +
  scale_fill_gradient(low = aaas_cols[1], high = aaas_cols[2]) +
  labs(title = "Jaccard mean", x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

# 6D. Jaccard index histograms (global + per region)
plot_jaccard_hist <- function(data, title_label, binwidth = 0.05) {
  ggplot(data, aes(x = Jaccard_Index)) +
    geom_histogram(binwidth = binwidth, fill = aaas_cols[1], color = "black", alpha = 0.8) +
    geom_vline(xintercept = low_thres, linetype = "dashed", color = aaas_cols[2], linewidth = 1) +
    geom_vline(xintercept = upp_thres, linetype = "dashed", color = aaas_cols[3], linewidth = 1) +
    scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    labs(title = title_label, x = "Jaccard Index", y = "Frequency (log10)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_jaccard_global  <- plot_jaccard_hist(jaccards_tb, "Global")
p_jaccard_regions <- purrr::map(
  levels(jaccards_tb$Region),
  ~ plot_jaccard_hist(dplyr::filter(jaccards_tb, Region == .x), .x, binwidth = 0.04)
)

p_jaccard_panel <- cowplot::plot_grid(
  p_jaccard_global,
  p_jaccard_regions[[1]], p_jaccard_regions[[2]],
  p_jaccard_regions[[3]], p_jaccard_regions[[4]], p_jaccard_regions[[5]],
  ncol = 3
)

# 6E. Module proportion barplots (local and global)
p_local <- ggplot(summary_all, aes(x = Region, y = local_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = fill_colors,
                    labels = c(
                      "AD_exclusive"      = "AD exclusive",
                      "Control_exclusive" = "Control exclusive",
                      "Intermediate"      = "Intermediate",
                      "Similar"           = "Similar"
                    )) +
  labs(title =  element_blank(), y = "%", x = "") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

p_global <- ggplot(summary_all, aes(x = Region, y = global_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = fill_colors,
                    labels = c("AD exclusive", "Control exclusive", "Intermediate", "Similar")) +
  labs(title = "Global proportion", y = "%", x = "") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

p_proportions <- cowplot::plot_grid(p_local, p_global, labels = c("A", "B"))

# 6F. Frequent genes in GO enrichments for conserved pairs
p_freq_genes <- purrr::imap_dfr(go_by_pair, function(ego, pair_name) {
  if (is.null(ego)) return(NULL)
  genes <- as.data.frame(ego)$geneID |>
    strsplit("/") |>
    unlist() |>
    table()
  tibble(pair = pair_name, gene = names(genes), count = as.numeric(genes))
}) |>
  dplyr::group_by(pair) |>
  dplyr::arrange(desc(count)) |>
  dplyr::ungroup() |>
  ggplot(aes(x = reorder(gene, count), y = count, fill = pair)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~pair, scales = "free_y", ncol = 1) +
  scale_fill_aaas() +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = element_blank(), y = "Frequency")

# ==============================================================================
# 7. PANEL ASSEMBLY AND EXPORT
# ==============================================================================

panel_counts <- cowplot::plot_grid(
  p_module_ridges, p_nmi, p_local,
  labels     = c("A", "B", "C"),
  rel_widths = c(2, 0.8, 2.5),
  label_size = 14,
  ncol       = 3
)
# Dunn post-hoc bubble plot (across regions, per phenotype)
# Each bubble is one region pair; size = -log10(p.adj); color = significance.
# Only pairs with p.adj < 0.05 get a filled bubble so non-significant pairs
# remain visible but clearly muted.
p_dunn_bubble <- stat_across_regions %>%
  dplyr::mutate(
    group1    = factor(group1, levels = ordered_regions),
    group2    = factor(group2, levels = rev(ordered_regions)),
    log_p     = -log10(p.adj + 1e-10),
    sig       = p_label != "ns"
  ) %>%
  ggplot(aes(x = group1, y = group2)) +
  geom_point(aes(size = log_p, fill = sig),
             shape = 21, color = "grey40", stroke = 0.4, alpha = 0.85) +
  geom_text(aes(label = p_label), size = 2.8, fontface = "bold", color = "white") +
  facet_wrap(~Phenotype, ncol = 2) +
  scale_fill_manual(
    values = c("TRUE" = aaas_cols[2], "FALSE" = "grey85"),
    labels = c("TRUE" = "p.adj < 0.05", "FALSE" = "ns"),
    name   = ""
  ) +
  scale_size_continuous(
    range = c(4, 14),
    name  = expression(-log[10](p[adj]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "bottom"
  ) +
  labs(x = "", y = "")


# 6G. cnet network plots for conserved module pairs
#
# Two bipartite graphs (GO terms <-> genes), one per conserved pair:
#   Pair A — PCC_AD_25 / DLPFC_AD_14  : 3 genes, 34 GO:BP terms (synaptic)
#   Pair B — PCC_Control_2 / DLPFC_Control_3 : 119 genes, 199 GO:BP terms (immune)
#
# go_by_pair already holds simplify()-filtered enrichResult objects.
#
# Visual encoding (consistent across both panels):
#   GO nodes   — shape  = square (22)
#                size   = rescaled Count  (n genes in that term)
#                fill   = -log10(p.adjust), gradient low→high
#   Gene nodes — shape  = circle (21)
#                size   = degree within this cnet (n GO terms it connects to)
#                fill   = grey (all genes; no hub distinction — hub data lives

#   Edges      — width  = 1 / n_terms_for_gene  (genes connecting many terms
#                         get thinner edges so the hubs stay readable)
#
# Layout: force-directed "fr" for the AD pair (small, clean star structure);
#         "stress" for the Control pair (better for dense immune graphs).

pacman::p_load("ggraph", "tidygraph", "ggnewscale")

# Build an annotated tidygraph from one cnet igraph
prepare_cnet_tg <- function(g, top_n_go = 20) {
  tg <- tidygraph::as_tbl_graph(g)
  
  # Keep only the top_n_go GO terms by p.adjust, then drop gene nodes
  # that become isolated after the GO subset.
  go_ids_keep <- tg %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble() %>%
    dplyr::filter(type == "GO", !is.na(p_adjust)) %>%
    dplyr::slice_min(p_adjust, n = top_n_go, with_ties = FALSE) %>%
    dplyr::pull(name)
  
  tg <- tg %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(type == "Gene" | name %in% go_ids_keep) %>%
    # Remove genes that lost all edges after the GO filter
    tidygraph::activate(edges) %>%
    tidygraph::filter(tidygraph::edge_is_incident(
      which(tidygraph::activate(tg, nodes) %>%
              tibble::as_tibble() %>%
              dplyr::pull(type) == "GO")
    )) %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(!tidygraph::node_is_isolated())
  
  # Compute degree while nodes are active, store as a node attribute,
  # then reference it from the edges context via .N() — the correct
  # tidygraph idiom for accessing node data inside activate(edges).
  tg %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      is_go      = (type == "GO"),
      fill_go    = ifelse(is_go, -log10(p_adjust + 1e-10), NA_real_),
      degree     = tidygraph::centrality_degree(),
      node_size  = ifelse(
        is_go,
        scales::rescale(tidygraph::replace_na(Count, 1L), to = c(5, 16)),
        scales::rescale(degree, to = c(3, 9))
      )
    ) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(
      # .N() gives access to the node data frame while edges are active
      gene_deg   = ifelse(
        .N()$type[from] == "Gene", .N()$degree[from], .N()$degree[to]
      ),
      edge_width = scales::rescale(1 / (gene_deg + 1), to = c(0.3, 1.5))
    )
}

# Layout helper: GO terms on an inner circle, genes on an outer circle.
# This gives the biologically intuitive "star" reading: functional terms
# in the center, genes radiating outward. Returns a data.frame with x, y
# columns in the same node order as the tidygraph object.
compute_concentric_layout <- function(tg, inner_r = 1, outer_r = 3) {
  nodes <- tg %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  
  go_idx   <- which(nodes$type == "GO")
  gene_idx <- which(nodes$type == "Gene")
  
  n_go   <- length(go_idx)
  n_gene <- length(gene_idx)
  
  # Evenly space angles, starting from the top (pi/2) going counterclockwise
  go_angles   <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n_go   + 1)[seq_len(n_go)]
  gene_angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n_gene + 1)[seq_len(n_gene)]
  
  coords <- data.frame(x = numeric(nrow(nodes)), y = numeric(nrow(nodes)))
  coords$x[go_idx]   <- inner_r * cos(go_angles)
  coords$y[go_idx]   <- inner_r * sin(go_angles)
  coords$x[gene_idx] <- outer_r * cos(gene_angles)
  coords$y[gene_idx] <- outer_r * sin(gene_angles)
  
  coords
}

# Core ggraph plot — shared aesthetics for both panels
# pair_title   : main title, e.g. "PCC AD 25 -- DLPFC AD 14"
# pair_subtitle: pathway label shown below title, e.g. "Synaptic transmission"
# inner_r / outer_r: radii of the GO and gene circles respectively
plot_cnet <- function(tg, pair_title, pair_subtitle,
                      go_color = aaas_cols[3],
                      seed     = 42,
                      inner_r  = 1,
                      outer_r  = 3) {
  set.seed(seed)
  
  go_low  <- scales::alpha(go_color, 0.25)
  go_high <- go_color
  
  layout_coords <- compute_concentric_layout(tg, inner_r = inner_r, outer_r = outer_r)
  
  ggraph::ggraph(tg, layout = "manual",
                 x = layout_coords$x, y = layout_coords$y) +
    
    # Edges: width encodes gene connectivity (thinner = more connected gene)
    ggraph::geom_edge_link(
      aes(width = edge_width),
      color = "grey75",
      alpha = 0.5
    ) +
    ggraph::scale_edge_width_identity() +
    
    # Gene nodes (drawn first — sit underneath GO squares)
    ggraph::geom_node_point(
      data   = function(x) dplyr::filter(x, !is_go),
      aes(size = node_size),
      fill   = "grey80",
      shape  = 21, color = "grey25", stroke = 0.35, alpha = 0.9
    ) +
    
    # Gene labels — bold, repelled outward from the circle center
    ggraph::geom_node_label(
      data          = function(x) dplyr::filter(x, !is_go),
      aes(label     = name),
      size          = 2.4,
      fill          = "grey97",
      color         = "grey10",
      fontface      = "bold",
      label.size    = 0.10,
      label.padding = unit(0.12, "lines"),
      repel         = TRUE,
      max.overlaps  = 40,
      seed          = seed
    ) +
    
    # GO term nodes on top of edges
    ggnewscale::new_scale_fill() +
    ggraph::geom_node_point(
      data   = function(x) dplyr::filter(x, is_go),
      aes(size = node_size, fill = fill_go),
      shape  = 22, color = "grey15", stroke = 0.35
    ) +
    scale_fill_gradient(
      low   = go_low,
      high  = go_high,
      name  = expression(-log[10](p[adj])),
      guide = guide_colorbar(order = 1, barwidth = 0.8, barheight = 5)
    ) +
    
    # GO term labels — plain, white background, repelled inward
    ggraph::geom_node_label(
      data          = function(x) dplyr::filter(x, is_go),
      aes(label     = label),
      size          = 2.5,
      fill          = "white",
      color         = "grey10",
      fontface      = "plain",
      label.size    = 0.12,
      label.padding = unit(0.15, "lines"),
      repel         = TRUE,
      max.overlaps  = 30,
      seed          = seed
    ) +
    
    scale_size_identity() +
    labs(title = pair_title, subtitle = pair_subtitle) +
    theme_void(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold",  size = 12, hjust = 0.5,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(face = "plain", size = 10, hjust = 0.5,
                                   color = "grey40", margin = margin(b = 8)),
      legend.position = "right",
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 8)
    )
}

# --- Build the two panels ---

# Pair A: AD synaptic pair — 3 genes, compact inner/outer radii
tg_ad <- prepare_cnet_tg(cnet_networks[["PCC_AD_25__DLPFC_AD_14"]], top_n_go = 15)
p_cnet_ad <- plot_cnet(
  tg_ad,
  pair_title    = "PCC AD 25 -- DLPFC AD 14",
  pair_subtitle = "Synaptic transmission",
  go_color      = aaas_cols[2],
  inner_r       = 1,
  outer_r       = 2.5
)

# Pair B: Control immune pair — many genes, larger outer radius to avoid crowding
tg_ctrl <- prepare_cnet_tg(cnet_networks[["PCC_Control_2__DLPFC_Control_3"]], top_n_go = 13)
p_cnet_ctrl <- plot_cnet(
  tg_ctrl,
  pair_title    = "PCC Control 2 -- DLPFC Control 3",
  pair_subtitle = "Microglial immune response",
  go_color      = aaas_cols[1],
  inner_r       = 1,
  outer_r       = 4     # wider outer ring for the larger gene set
)
p_cnet_ctrl

# Assemble panel A (AD) left, B (Control) right
panel_cnet <- cowplot::plot_grid(
  p_cnet_ad, p_cnet_ctrl,
  labels     = c("A", "B"),
  label_size = 14,
  ncol       = 2,
  rel_widths = c(1, 1.6)   # Control panel wider — more nodes
)

# Save panels
save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(output_dir, paste0(filename, ".pdf")),  plot = plot, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".jpeg")), plot = plot, width = width, height = height, dpi = 300)
}

save_plot(panel_counts, "panel_module_counts",width = 14, height = 14)
save_plot(p_dunn_bubble, "module_size_dunn",        width = 10, height = 5)
save_plot(p_jaccard_panel, "jaccard_histogram",      width = 18, height = 10)
save_plot(p_freq_genes,    "freq_genes_enrichment",  width = 4,  height = 10)
if (!is.null(panel_cnet)) save_plot(panel_cnet, "cnet_conserved_pairs", width = 14, height = 8)

message("Done! All outputs saved to: ", output_dir)
