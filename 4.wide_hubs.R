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
#
# CHANGES FROM PREVIOUS VERSION (2 targeted fixes):
#   FIX 1 — Step 5: hubs_long and hubs_pivot now use `node` (Ensembl ID) as
#            primary key instead of `symbol`. find_exclusive_hubs groups by
#            `node` + `symbol`. This prevents symbol collisions (e.g. KATNB1
#            mapping to two different Ensembl IDs) from incorrectly placing a
#            gene in both the AD-exclusive and CTR-exclusive lists.
#   FIX 2 — Step 5: After computing wide_AD and wide_CTR, any symbol appearing
#            in both lists is flagged, logged, and removed from both before
#            downstream use.
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
opt$limma_dir <- "~/Desktop/local_work/fomo_networks/DEGS"
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

# [FIX 1] Find exclusive wide hubs for one phenotype.
# Now groups by node (Ensembl ID) + symbol to avoid symbol collision artifacts.
# A hub is exclusive if it appears as a hub in `target` but NEVER in `other`.
# A hub is "wide" if it is exclusive in >= min_regions regions.
find_exclusive_hubs <- function(pivot_df, target, other, min_regions) {
  pivot_df %>%
    filter(.data[[target]] == TRUE, .data[[other]] == FALSE) %>%
    group_by(node, symbol) %>%                        # FIX 1: group by node, not symbol
    summarise(
      n_regions = n(),
      regions   = paste(Region, collapse = ", "),
      .groups   = "drop"
    ) %>%
    filter(n_regions >= min_regions) %>%
    arrange(desc(n_regions))
}

# Build a binary gene x region presence matrix (for pheatmap)
build_hub_matrix <- function(pivot_df, wide_hubs, target_col) {
  pivot_df %>%
    filter(node %in% wide_hubs$node, .data[[target_col]] == TRUE) %>%  # FIX 1: filter by node
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

# [FIX 1] One row per (node, symbol, Region, Phenotype) — node is the primary key
hubs_long <- node_data %>%
  filter(is_hub) %>%
  distinct(node, symbol, Region, Phenotype) %>%       # FIX 1: keep node
  mutate(hub = TRUE)

# [FIX 1] Pivot on node + symbol + Region, not symbol alone
hubs_pivot <- hubs_long %>%
  pivot_wider(
    id_cols     = c(node, symbol, Region),            # FIX 1: node as unique key
    names_from  = Phenotype,
    values_from = hub,
    values_fill = FALSE
  )

# Exclusive hubs in each phenotype
wide_AD  <- find_exclusive_hubs(hubs_pivot, "AD",      "Control", opt$min_regions)
wide_CTR <- find_exclusive_hubs(hubs_pivot, "Control", "AD",      opt$min_regions)

# [FIX 2] Remove any symbol appearing in both lists — these are genuinely
# ambiguous genes (hub in AD in some regions, hub in Control in others)
# and cannot be assigned to either exclusive category
ambiguous_symbols <- intersect(wide_AD$symbol, wide_CTR$symbol)

if (length(ambiguous_symbols) > 0) {
  message(sprintf(
    "WARNING: %d symbol(s) appear in both AD and CTR exclusive lists and will be removed: %s",
    length(ambiguous_symbols),
    paste(ambiguous_symbols, collapse = ", ")
  ))
  wide_AD  <- wide_AD  %>% filter(!symbol %in% ambiguous_symbols)
  wide_CTR <- wide_CTR %>% filter(!symbol %in% ambiguous_symbols)
}

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

membership_annotated_AD  <- membership_annotated %>% filter(symbol %in% rownames(mat_AD))
write_csv(membership_annotated_AD, "membership_annotated_AD.csv")

membership_annotated_CTR <- membership_annotated %>% filter(symbol %in% rownames(mat_CTR))
write_csv(membership_annotated_CTR, "membership_annotated_CTR.csv")

# -----------------------------------------------------------------------------
# 10. Q4: Are AD-exclusive wide hubs also differentially expressed?
#
# For each AD-exclusive wide hub gene, we ask whether it is also a DEG
# (FDR <= 0.05 and |log2FC| >= 0.5) in the regions where it acts as a hub.
# This distinguishes transcriptional rewiring (hub status without expression
# change) from abundance-driven centrality (hub status co-occurring with DEG).
# -----------------------------------------------------------------------------

# Load limma results for all regions
# Files must follow the naming pattern: *_{REGION}_limma_*.tsv
limma_files <- list.files(opt$limma_dir, pattern = "limma.*allGenes\\.tsv$", full.names = TRUE)

deg_thresholds <- list(fdr = 0.05, logfc = 0.5)

limma_df <- map_df(limma_files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      Region = str_extract(basename(f), "CRB|TC|DLPFC|HCN|PCC"),
      is_DEG = adj.P.Val <= deg_thresholds$fdr & abs(logFC) >= deg_thresholds$logfc
    ) %>%
    dplyr::rename(node = gene)
})

# For each AD-exclusive hub, retrieve its DEG status in the regions
# where it was identified as a hub
hub_deg_status <- wide_AD %>%
  # Expand one row per (node, region where it is a hub)
  separate_rows(regions, sep = ", ") %>%
  dplyr::rename(Region = regions) %>%
  left_join(
    limma_df %>% dplyr::select(node, Region, logFC, adj.P.Val, is_DEG),
    by = c("node", "Region")
  ) %>%
  mutate(
    deg_direction = case_when(
      is_DEG & logFC > 0 ~ "up",
      is_DEG & logFC < 0 ~ "down",
      !is_DEG            ~ "not DEG",
      TRUE               ~ NA_character_
    )
  )

write_csv(hub_deg_status, "AD_hubs_deg_status.csv")

# Summary: for each hub gene, how many of its hub-regions are also DEG regions?
hub_deg_summary <- hub_deg_status %>%
  group_by(symbol, n_regions) %>%
  summarise(
    n_hub_regions     = n(),
    n_deg_regions     = sum(is_DEG, na.rm = TRUE),
    n_up_regions      = sum(deg_direction == "up",      na.rm = TRUE),
    n_down_regions    = sum(deg_direction == "down",    na.rm = TRUE),
    n_notdeg_regions  = sum(deg_direction == "not DEG", na.rm = TRUE),
    deg_regions       = paste(Region[is_DEG],  collapse = ", "),
    notdeg_regions    = paste(Region[!is_DEG], collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_deg_regions), desc(n_hub_regions))

write_csv(hub_deg_summary, "AD_hubs_deg_summary.csv")

message(sprintf(
  "AD-exclusive hubs with DEG overlap in >= 1 hub region: %d / %d",
  sum(hub_deg_summary$n_deg_regions >= 1),
  nrow(hub_deg_summary)
))

# Plot: heatmap of hub status × DEG direction per region
# Rows = AD-exclusive hub genes, Columns = regions
# Fill = DEG direction (up / down / not DEG / not hub)

region_order <- c("CRB", "TC", "DLPFC", "HCN", "PCC")

hub_deg_plot_df <- hub_deg_status %>%
  mutate(
    Region = factor(Region, levels = region_order),
    deg_direction = factor(deg_direction, levels = c("up", "down", "not DEG"))
  )

# Order genes by total hub regions (descending), same as mat_AD
gene_order <- rownames(mat_AD)

hub_deg_tile <- hub_deg_plot_df %>%
  mutate(symbol = factor(symbol, levels = rev(gene_order))) %>%
  ggplot(aes(x = Region, y = symbol, fill = deg_direction)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c("up" = "#B2182B", "down" = "#2166AC", "not DEG" = "gray85"),
    na.value = "gray95",
    name = "Expression\nchange in AD"
  ) +
  labs(x = NULL, y = NULL,
       title = "AD-exclusive hub genes: DEG status per hub region") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

#Vis
hub_deg_tile


ggsave("AD_hubs_deg_heatmap.jpeg", plot = hub_deg_tile,
       width = 5, height = 6, dpi = 300)

limma_files <- list.files(opt$limma_dir, pattern = "limma.*allGenes\\.tsv$", full.names = TRUE)

# =============================================================================
# END
# =============================================================================