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
#   Q4. Are AD-exclusive wide hubs also differentially expressed?
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
# 1. Paths
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_dir"),   type = "character", metavar = "path",
              help = "Directory with *_nodes_summary.csv files"),
  make_option(c("-e", "--enrich_dir"),  type = "character", metavar = "path",
              help = "Directory with *_enrichment.csv files"),
  make_option(c("-d", "--out_dir"),     type = "character", metavar = "path",
              default = "results_wide_hubs", help = "Output directory [%default]"),
  make_option(c("-n", "--min_regions"), type = "integer",   default = 2L,
              metavar = "int", help = "Minimum regions for a hub to be 'wide' [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Hardcoded values (comment out when using CLI) ---
opt$input_dir   <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/"
opt$enrich_dir  <- "~/Desktop/local_work/fomo_networks/results_topos_louvain/results_comm/"
opt$out_dir     <- "~/Desktop/local_work/fomo_networks/results_wide_hubs"
opt$limma_dir   <- "~/Desktop/local_work/fomo_networks/DEGS"
opt$min_regions <- 2L
# ------------------------------------------------------

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$out_dir)

region_order <- c("HCN", "DLPFC", "PCC", "CRB", "TC")

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

# Find exclusive wide hubs for one phenotype.
# Groups by node (Ensembl ID) + symbol to avoid symbol collision artifacts.
find_exclusive_hubs <- function(pivot_df, target, other, min_regions) {
  pivot_df %>%
    filter(.data[[target]] == TRUE, .data[[other]] == FALSE) %>%
    group_by(node, symbol) %>%
    summarise(
      n_regions = n(),
      regions   = paste(Region, collapse = ", "),
      .groups   = "drop"
    ) %>%
    filter(n_regions >= min_regions) %>%
    arrange(desc(n_regions))
}

# Build a binary gene x region presence matrix
build_hub_matrix <- function(pivot_df, wide_hubs, target_col) {
  pivot_df %>%
    filter(node %in% wide_hubs$node, .data[[target_col]] == TRUE) %>%
    dplyr::select(symbol, Region) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = Region, values_from = value, values_fill = 0L) %>%
    column_to_rownames("symbol") %>%
    as.matrix() %>%
    .[wide_hubs$symbol, , drop = FALSE]
}

# Build a logFC gene x region matrix
build_logfc_matrix <- function(pivot_df, wide_hubs, target_col, limma_df) {
  pivot_df %>%
    filter(node %in% wide_hubs$node, .data[[target_col]] == TRUE) %>%
    dplyr::select(node, symbol, Region) %>%
    left_join(
      limma_df %>% dplyr::select(node, Region, logFC),
      by = c("node", "Region")
    ) %>%
    dplyr::select(symbol, Region, logFC) %>%
    pivot_wider(names_from = Region, values_from = logFC, values_fill = NA_real_) %>%
    column_to_rownames("symbol") %>%
    as.matrix() %>%
    .[wide_hubs$symbol, , drop = FALSE]
}

# Build a long data frame combining logFC and hub binary status
build_heatmap_df <- function(mat_logfc, mat_binary, col_order) {
  cols_ordered <- intersect(col_order, colnames(mat_logfc))
  
  df_logfc <- mat_logfc[, cols_ordered, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("symbol") %>%
    pivot_longer(-symbol, names_to = "Region", values_to = "logFC")
  
  df_hub <- mat_binary[, cols_ordered, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("symbol") %>%
    pivot_longer(-symbol, names_to = "Region", values_to = "is_hub")
  
  left_join(df_logfc, df_hub, by = c("symbol", "Region")) %>%
    mutate(
      Region    = factor(Region, levels = cols_ordered),
      cell_type = case_when(
        is_hub == 1 & !is.na(logFC) ~ "hub",
        is_hub == 0                 ~ "not_hub",
        is.na(logFC)                ~ "absent"
      )
    )
}

# ggplot heatmap with logFC fill + border on hub cells
# Genes ordered from most to fewest hub regions (top → bottom)
plot_heatmap_gg <- function(mat_logfc, mat_binary, wide_hubs,
                            title_label,
                            col_order    = region_order,
                            fontsize_row = 3,
                            border_color = "black",
                            border_lwd   = 0.6) {
  
  cols_ordered <- intersect(col_order, colnames(mat_logfc))
  
  # En plot_heatmap_gg, reemplaza gene_levels por:
  gene_levels <- tibble(
    symbol    = rownames(mat_logfc),
    n_regions = rowSums(!is.na(mat_logfc))  # cuenta regiones con logFC (= hub regions)
  ) %>%
    arrange(n_regions, symbol) %>%   # ascendente + rev() → más regiones arriba
    pull(symbol) %>%
    rev()
  
  df <- build_heatmap_df(mat_logfc, mat_binary, col_order) %>%
    mutate(symbol = factor(symbol, levels = gene_levels))
  
  max_fc <- max(abs(df$logFC), na.rm = TRUE)
  
  ggplot(df, aes(x = Region, y = symbol)) +
    geom_tile(aes(fill = logFC), color = "white", linewidth = 0.2) +
    geom_tile(
      data      = filter(df, cell_type == "hub"),
      fill      = NA,
      color     = border_color,
      linewidth = border_lwd
    ) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-max_fc, max_fc),
      na.value = "#D9D9D9",
      name     = "log2FC"
    ) +
    labs(title = title_label, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y     = element_text(size = fontsize_row),
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

# Run GO:BP enrichment (SYMBOL keys)
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

# Retrieve module membership for a gene set in one phenotype
get_membership <- function(genes, phenotype_label, node_data) {
  node_data %>%
    filter(Phenotype == phenotype_label, symbol %in% genes) %>%
    distinct(symbol, Region, membership) %>%
    dplyr::rename(!!paste0("membership_", phenotype_label) := membership)
}

# -----------------------------------------------------------------------------
# 3. Load node-level data
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
  keys      = node_data$node,
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)

node_data <- node_data %>%
  left_join(tibble(node = names(symbols), symbol = unname(symbols)), by = "node") %>%
  mutate(symbol = coalesce(symbol, node))

# -----------------------------------------------------------------------------
# 5. Q1: Identify exclusive wide hub genes
# -----------------------------------------------------------------------------

hubs_long <- node_data %>%
  filter(is_hub) %>%
  distinct(node, symbol, Region, Phenotype) %>%
  mutate(hub = TRUE)

hubs_pivot <- hubs_long %>%
  pivot_wider(
    id_cols     = c(node, symbol, Region),
    names_from  = Phenotype,
    values_from = hub,
    values_fill = FALSE
  )

wide_AD  <- find_exclusive_hubs(hubs_pivot, "AD",      "Control", opt$min_regions)
wide_CTR <- find_exclusive_hubs(hubs_pivot, "Control", "AD",      opt$min_regions)

# Remove ambiguous symbols (hub in AD in some regions, Control in others)
ambiguous_symbols <- intersect(wide_AD$symbol, wide_CTR$symbol)

if (length(ambiguous_symbols) > 0) {
  message(sprintf(
    "WARNING: %d symbol(s) appear in both AD and CTR exclusive lists and will be removed: %s",
    length(ambiguous_symbols), paste(ambiguous_symbols, collapse = ", ")
  ))
  wide_AD  <- wide_AD  %>% filter(!symbol %in% ambiguous_symbols)
  wide_CTR <- wide_CTR %>% filter(!symbol %in% ambiguous_symbols)
}

message(sprintf("Wide AD-exclusive hubs  : %d genes", nrow(wide_AD)))
message(sprintf("Wide CTR-exclusive hubs : %d genes", nrow(wide_CTR)))

# Binary presence matrices (sorted by n_regions descending — matches gene_levels in plots)
mat_AD <- build_hub_matrix(hubs_pivot, wide_AD, "AD") %>%
  as.data.frame() %>%
  mutate(Present = rowSums(.)) %>%
  arrange(desc(Present))

mat_CTR <- build_hub_matrix(hubs_pivot, wide_CTR, "Control") %>%
  as.data.frame() %>%
  mutate(Present = rowSums(.)) %>%
  arrange(desc(Present))

# Binary matrices without "Present" column (used in heatmap functions)
mat_AD_bin  <- as.data.frame(mat_AD[,  colnames(mat_AD)  != "Present"]) %>% 
  mutate(n_regions = rowSums(.)) %>% 
  arrange(desc(n_regions)) %>% 
  subset(select = -c(n_regions)) %>% 
  as.matrix()

mat_CTR_bin <- as.data.frame(mat_CTR[, colnames(mat_CTR) != "Present"]) %>% 
  mutate(n_regions = rowSums(.)) %>% 
  arrange(desc(n_regions)) %>% 
  subset(select = -c(n_regions)) %>% 
  as.matrix()

# -----------------------------------------------------------------------------
# 6. Q2: Functional enrichment of exclusive hub genes
# -----------------------------------------------------------------------------

enrich_AD  <- run_enrichment(wide_AD$symbol,  label = "AD-exclusive hubs")
enrich_CTR <- run_enrichment(wide_CTR$symbol, label = "CTR-exclusive hubs")

if (nrow(as.data.frame(enrich_AD)) > 0) {
  cnet_AD <- cnetplot(enrich_AD, node_label_size = NULL)
  ggsave("cnetplot_AD_hubs.jpeg", plot = cnet_AD, width = 10, height = 10)
}

if (nrow(as.data.frame(enrich_CTR)) > 0) {
  cnet_CTR <- cnetplot(enrich_CTR, node_label_size = NULL)
  ggsave("cnetplot_CTR_hubs.jpeg", plot = cnet_CTR, width = 10, height = 10)
}

# -----------------------------------------------------------------------------
# 7. Q3: Module membership of exclusive hub genes
# -----------------------------------------------------------------------------

membership_all <- bind_rows(
  full_join(
    get_membership(wide_AD$symbol,  "AD",      node_data),
    get_membership(wide_AD$symbol,  "Control", node_data),
    by = c("symbol", "Region")
  ) %>% mutate(hub_type = "AD_exclusive"),
  
  full_join(
    get_membership(wide_CTR$symbol, "AD",      node_data),
    get_membership(wide_CTR$symbol, "Control", node_data),
    by = c("symbol", "Region")
  ) %>% mutate(hub_type = "Control_exclusive")
)

# Annotate modules with top enriched biological process
enrich_files <- list.files(opt$enrich_dir, pattern = "_enrichment\\.csv$", full.names = TRUE)

enrich_df <- map_df(enrich_files, function(f) {
  vroom::vroom(f, show_col_types = FALSE) %>%
    mutate(
      file      = basename(f),
      Region    = str_extract(file, "(?<=Mayo_)[A-Z]+"),
      Phenotype = if_else(str_detect(file, "_AD_"), "AD", "Control")
    )
})

module_top_term <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Region, Phenotype, CommunityID) %>%
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Region, Phenotype, CommunityID, top_function = Description)

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

write_csv(membership_annotated,                                              "membership_annotated.csv")
write_csv(membership_annotated %>% filter(symbol %in% rownames(mat_AD)),    "membership_annotated_AD.csv")
write_csv(membership_annotated %>% filter(symbol %in% rownames(mat_CTR)),   "membership_annotated_CTR.csv")

# -----------------------------------------------------------------------------
# 8. Q4: Are AD-exclusive wide hubs also differentially expressed?
# -----------------------------------------------------------------------------

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

hub_deg_status <- wide_AD %>%
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

hub_deg_summary <- hub_deg_status %>%
  group_by(symbol, n_regions) %>%
  summarise(
    n_hub_regions    = n(),
    n_deg_regions    = sum(is_DEG,  na.rm = TRUE),
    n_up_regions     = sum(deg_direction == "up",      na.rm = TRUE),
    n_down_regions   = sum(deg_direction == "down",    na.rm = TRUE),
    n_notdeg_regions = sum(deg_direction == "not DEG", na.rm = TRUE),
    deg_regions      = paste(Region[is_DEG],  collapse = ", "),
    notdeg_regions   = paste(Region[!is_DEG], collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_deg_regions), desc(n_hub_regions))

write_csv(hub_deg_summary, "AD_hubs_deg_summary.csv")

message(sprintf(
  "AD-exclusive hubs with DEG overlap in >= 1 hub region: %d / %d",
  sum(hub_deg_summary$n_deg_regions >= 1), nrow(hub_deg_summary)
))

# DEG tile heatmap
gene_order_ad <- rownames(mat_AD)   # already sorted desc by n_regions

hub_deg_tile <- hub_deg_status %>%
  mutate(
    Region        = factor(Region, levels = region_order),
    deg_direction = factor(deg_direction, levels = c("up", "down", "not DEG")),
    symbol        = factor(symbol, levels = rev(gene_order_ad))
  ) %>%
  ggplot(aes(x = Region, y = symbol, fill = deg_direction)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values   = c("up" = "#B2182B", "down" = "#2166AC", "not DEG" = "gray85"),
    na.value = "gray95",
    name     = "Expression\nchange in AD"
  ) +
  labs(x = NULL, y = NULL,
       title = "AD-exclusive hub genes: DEG status per hub region") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave("AD_hubs_deg_heatmap.jpeg", plot = hub_deg_tile, width = 5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 9. LogFC heatmaps with hub borders + single shared legend
# -----------------------------------------------------------------------------
# Justo antes de llamar plot_heatmap_gg, re-ordena explícitamente
wide_AD  <- wide_AD  %>% arrange(desc(n_regions))
wide_CTR <- wide_CTR %>% arrange(desc(n_regions))

mat_logfc_AD  <- build_logfc_matrix(hubs_pivot, wide_AD,  "AD",      limma_df)
mat_logfc_CTR <- build_logfc_matrix(hubs_pivot, wide_CTR, "Control", limma_df)

# Align binary matrices rows to logfc matrices
mat_AD_bin  <- mat_AD_bin[intersect(rownames(mat_logfc_AD),  rownames(mat_AD_bin)),  , drop = FALSE]
mat_CTR_bin <- mat_CTR_bin[intersect(rownames(mat_logfc_CTR), rownames(mat_CTR_bin)), , drop = FALSE]

# Build panels
ht_ad_fc  <- plot_heatmap_gg(
  mat_logfc_AD,  mat_AD_bin,  wide_AD,
  title_label  = "AD-exclusive hubs",
  fontsize_row = 10
)

ht_ctr_fc <- plot_heatmap_gg(
  mat_logfc_CTR, mat_CTR_bin, wide_CTR,
  title_label  = "Control-exclusive hubs",
  fontsize_row = 3
)

# Extract legend from left panel; render both panels without legend
legend_only     <- cowplot::get_legend(ht_ad_fc + theme(legend.position = "right"))
ht_ad_fc_noleg  <- ht_ad_fc  + theme(legend.position = "none")
ht_ctr_fc_noleg <- ht_ctr_fc + theme(legend.position = "none",
                                     axis.text.y = element_text(size = 7))

hts_fc <- cowplot::plot_grid(
  cowplot::plot_grid(
    ht_ad_fc_noleg, ht_ctr_fc_noleg,
    labels     = c("A", "B"),
    label_size = 10,
    ncol       = 2,
    rel_widths = c(1, 2)
  ),
  legend_only,
  ncol       = 2,
  rel_widths = c(1, 0.08)
)

ggsave("hub_heatmaps_logfc.jpeg", plot = hts_fc, width = 12, height = 8, dpi = 300)

message("Done! All outputs saved to: ", opt$out_dir)

# =============================================================================
# END
# =============================================================================