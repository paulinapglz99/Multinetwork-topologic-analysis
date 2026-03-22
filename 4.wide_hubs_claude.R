#!/usr/bin/env Rscript
# =============================================================================
# wide_hubs_analysis.R
# Identify phenotype-exclusive hub genes conserved across anatomical regions,
# characterize their functional enrichment, and map their modular localization.
#
# Questions addressed:
#   Q1. Are there genes that are exclusive hubs in AD (or Control) and conserved
#       across multiple brain regions?
#   Q2. What biological processes are associated with those hub genes?
#   Q3. Where do these hubs land in the community structure (module membership)?
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Dependencies
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  tidyverse, vroom,
  pheatmap, cowplot, ggplotify,
  clusterProfiler, org.Hs.eg.db, AnnotationDbi,
  optparse
)

# -----------------------------------------------------------------------------
# 1. Paths  (optparse kept for future CLI use; values are hardcoded below)
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", metavar = "path",
              help = "Directory with *_nodes_summary.csv files"),
  make_option(c("-e", "--enrich_dir"), type = "character", metavar = "path",
              help = "Directory with *_enrichment.csv files"),
  make_option(c("-d", "--out_dir"),   type = "character", metavar = "path",
              default = "results_wide_hubs", help = "Output directory [%default]"),
  make_option(c("-n", "--min_regions"), type = "integer", default = 2L,
              metavar = "int", help = "Minimum regions for a hub to be 'wide' [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Hardcoded values (comment out when using CLI) ---
opt$input_dir  <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/"
opt$enrich_dir <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/results_comm/"
opt$out_dir    <- "~/Desktop/local_work/fomo_networks/results_wide_hubs"
opt$min_regions <- 2L
# ------------------------------------------------------

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$out_dir)

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

# Parse region and phenotype from a filename that contains "_counts_"
parse_filename_meta <- function(filenames) {
  tibble(file = filenames) %>%
    mutate(
      Region    = gsub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", basename(file))),
      Phenotype = if_else(str_detect(file, "_AD_"), "AD", "Control"),
      network   = paste0(Region, "_", Phenotype)
    )
}

# Define dual-criterion hubs (top 10% degree AND top 10% pagerank) per network
flag_hubs <- function(df) {
  df %>%
    group_by(file) %>%
    mutate(
      hub_degree   = degree   >= quantile(degree,   0.90, na.rm = TRUE),
      hub_pagerank = pagerank >= quantile(pagerank, 0.90, na.rm = TRUE),
      is_hub       = hub_degree & hub_pagerank
    ) %>%
    ungroup()
}

# Find exclusive wide hubs for one phenotype
# A hub is exclusive if it appears as a hub in `target` but NEVER in `other`.
# A hub is "wide" if it is exclusive in >= min_regions regions.
find_exclusive_hubs <- function(pivot_df, target, other, min_regions) {
  pivot_df %>%
    filter(.data[[target]] == TRUE, .data[[other]] == FALSE) %>%
    group_by(symbol) %>%
    summarise(
      n_regions = n(),
      regions   = paste(Region, collapse = ", "),
      .groups   = "drop"
    ) %>%
    filter(n_regions >= min_regions) %>%
    arrange(desc(n_regions))
}

# Build a binary gene × region presence matrix (for pheatmap)
build_hub_matrix <- function(pivot_df, wide_hubs, target_col) {
  pivot_df %>%
    filter(symbol %in% wide_hubs$symbol, .data[[target_col]] == TRUE) %>%
    dplyr::select(symbol, Region) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = Region, values_from = value, values_fill = 0L) %>%
    column_to_rownames("symbol") %>%
    as.matrix() %>%
    .[wide_hubs$symbol, , drop = FALSE]   # preserve ranking order
}

# Run GO:BP enrichment for a gene set (SYMBOL keys)
run_enrichment <- function(gene_symbols, label) {
  result <- enrichGO(
    gene          = gene_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  )
  message(sprintf("[%s] GO:BP enrichment: %d terms at FDR < 0.10", label, nrow(as.data.frame(result))))
  result
}

# Plot a pheatmap and return it as a ggplot object
hub_heatmap <- function(mat, fill_high, fontsize_row = 10) {
  ht <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    legend       = FALSE,
    main         = NA,
    fontsize_row = fontsize_row,
    color        = colorRampPalette(c("gray95", fill_high))(100),
    silent       = TRUE
  )
  ggplotify::as.ggplot(ht)
}

# -----------------------------------------------------------------------------
# 3. Load node-level data from all networks
# -----------------------------------------------------------------------------

node_files <- list.files(opt$input_dir, pattern = "_nodes_summary\\.csv$", full.names = TRUE)

node_data <- map_df(node_files, ~ vroom::vroom(.x, show_col_types = FALSE) %>%
                      mutate(file = basename(.x))) %>%
  bind_cols(parse_filename_meta(.$file) %>% dplyr::select(-file))

# -----------------------------------------------------------------------------
# 4. Define hubs and map Ensembl IDs to gene symbols
# -----------------------------------------------------------------------------

node_data <- flag_hubs(node_data)

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys    = node_data$node,
  keytype = "ENSEMBL",
  column  = "SYMBOL",
  multiVals = "first"
)

node_data <- node_data %>%
  left_join(tibble(node = names(symbols), symbol = unname(symbols)), by = "node") %>%
  mutate(symbol = coalesce(symbol, node))   # fall back to Ensembl ID if no symbol

# -----------------------------------------------------------------------------
# 5. Q1: Identify exclusive wide hub genes
# -----------------------------------------------------------------------------

# One row per (symbol, Region, Phenotype) combination
hubs_long <- node_data %>%
  filter(is_hub) %>%
  distinct(symbol, Region, Phenotype) %>%
  mutate(hub = TRUE)

# gene × Region × Phenotype pivot
hubs_pivot <- hubs_long %>%
  pivot_wider(names_from = Phenotype, values_from = hub, values_fill = FALSE)

# Exclusive hubs in each phenotype
wide_AD  <- find_exclusive_hubs(hubs_pivot, "AD",      "Control", opt$min_regions)
wide_CTR <- find_exclusive_hubs(hubs_pivot, "Control", "AD",      opt$min_regions)

message(sprintf("Wide AD-exclusive hubs  : %d genes", nrow(wide_AD)))
message(sprintf("Wide CTR-exclusive hubs : %d genes", nrow(wide_CTR)))

# Binary presence matrices
mat_AD  <- build_hub_matrix(hubs_pivot, wide_AD,  "AD") %>% as.data.frame()
mat_AD$Present <- rowSums(mat_AD)
mat_AD <- mat_AD %>% arrange(desc(Present))
mat_CTR <- build_hub_matrix(hubs_pivot, wide_CTR, "Control") %>% as.data.frame()
mat_CTR$Present <- rowSums(mat_CTR)
mat_CTR <- mat_CTR %>% arrange(desc(Present))

# -----------------------------------------------------------------------------
# 6. Plot heatmaps for exclusive hubs
# -----------------------------------------------------------------------------

ht_ad  <- hub_heatmap(mat_AD[1:5],  fill_high = "firebrick",     fontsize_row = 10)
ht_ctr <- hub_heatmap(mat_CTR[1:5], fill_high = "cornflowerblue", fontsize_row = 4)

hts <- plot_grid(ht_ad, ht_ctr, labels = c("A", "B"), label_size = 10, ncol = 2)

ggsave("hub_heatmaps.jpeg", plot = hts, width = 4, height = 6)

# -----------------------------------------------------------------------------
# 7. Q2: Functional enrichment of exclusive hub genes
# -----------------------------------------------------------------------------

enrich_AD  <- run_enrichment(wide_AD$symbol,  label = "AD-exclusive hubs")
enrich_CTR <- run_enrichment(wide_CTR$symbol, label = "CTR-exclusive hubs")

# cnetplots (only if enrichment returned results)
if (nrow(as.data.frame(enrich_AD)) > 0) {
  cnet_AD <- cnetplot(enrich_AD, node_label_size = NULL)
  ggsave("cnetplot_AD_hubs.jpeg", plot = cnet_AD, width = 10, height = 10)
}

if (nrow(as.data.frame(enrich_CTR)) > 0) {
  cnet_CTR <- cnetplot(enrich_CTR, node_label_size = NULL)
  ggsave("cnetplot_CTR_hubs.jpeg", plot = cnet_CTR, width = 10, height = 10)
}

# -----------------------------------------------------------------------------
# 8. Q3: Module membership of exclusive hub genes
# -----------------------------------------------------------------------------

# Retrieve membership for a set of genes in one phenotype
get_membership <- function(genes, phenotype_label) {
  node_data %>%
    filter(Phenotype == phenotype_label, symbol %in% genes) %>%
    distinct(symbol, Region, membership) %>%
    dplyr::rename(!!paste0("membership_", phenotype_label) := membership)
}

membership_all <- bind_rows(
  # AD-exclusive hubs: track where they land in both AD and Control networks
  full_join(
    get_membership(wide_AD$symbol, "AD"),
    get_membership(wide_AD$symbol, "Control"),
    by = c("symbol", "Region")
  ) %>% mutate(hub_type = "AD_exclusive"),
  
  # Control-exclusive hubs
  full_join(
    get_membership(wide_CTR$symbol, "AD"),
    get_membership(wide_CTR$symbol, "Control"),
    by = c("symbol", "Region")
  ) %>% mutate(hub_type = "Control_exclusive")
)

# -----------------------------------------------------------------------------
# 9. Annotate modules with their top enriched biological process
# -----------------------------------------------------------------------------

enrich_files <- list.files(opt$enrich_dir, pattern = "_enrichment\\.csv$", full.names = TRUE)

enrich_df <- map_df(enrich_files, function(f) {
  vroom::vroom(f, show_col_types = FALSE) %>%
    mutate(
      file      = basename(f),
      Region    = str_extract(file, "(?<=Mayo_)[A-Z]+"),
      Phenotype = if_else(str_detect(file, "_AD_"), "AD", "Control")
    )
})

# Top term per module (lowest adjusted p-value)
module_top_term <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Region, Phenotype, CommunityID) %>%
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Region, Phenotype, CommunityID, top_function = Description)

# Annotate membership table
membership_annotated <- membership_all %>%
  dplyr::rename(CommunityID_AD = membership_AD, CommunityID_Control = membership_Control) %>%
  left_join(
    module_top_term %>% filter(Phenotype == "AD") %>%
      dplyr::rename(CommunityID_AD = CommunityID, Function_AD = top_function) %>%
      dplyr::select(-Phenotype),
    by = c("Region", "CommunityID_AD")
  ) %>%
  left_join(
    module_top_term %>% filter(Phenotype == "Control") %>%
      dplyr::rename(CommunityID_Control = CommunityID, Function_Control = top_function) %>%
      dplyr::select(-Phenotype),
    by = c("Region", "CommunityID_Control")
  )

write_csv(membership_annotated, "membership_annotated.csv")

membership_annotated_AD <- membership_annotated %>% filter(symbol %in% rownames(mat_AD))
write_csv(membership_annotated_AD, "membership_annotated_AD.csv")


membership_annotated_CTR <- membership_annotated %>% filter(symbol %in% rownames(mat_CTR))
write_csv(membership_annotated_CTR, "membership_annotated_CTR.csv")

# =============================================================================
# END
# =============================================================================