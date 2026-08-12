#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(vroom)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(biomaRt)
  library(grid)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG (repo-relative; run from the repository root) ===================== #
# NOTE: this is an intermediate component consumed by dge_gsea_panel.R (Figure 6),
# not a final manuscript figure itself, so it stays under outputs/, not plots/.
source("config/paths.R")

in_dir  <- file.path(OUTPUTS_DIR, "dge_limma")
out_dir <- file.path(OUTPUTS_DIR, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 0.5
TOP_N_PER_SIDE <- 7

P_FLOOR <- 1e-300

out_png <- file.path(out_dir, "volcano_by_region_AD_vs_control.png")
out_pdf <- file.path(out_dir, "volcano_by_region_AD_vs_control.pdf")
out_rds <- file.path(out_dir, "volcano_by_region_AD_vs_control.rds")

# ===================== HELPERS ===================== #
read_limma_tsv <- function(path) {
  df <- read_tsv(path, col_types = cols(), progress = FALSE)
  needed <- c("gene", "logFC", "P.Value", "adj.P.Val")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(path), ": ", paste(miss, collapse = ", "))

  region_raw <- basename(path) %>%
    str_replace("_limma_AD_vs_Control_allGenes\\.tsv$", "")

  df %>% mutate(region_raw = region_raw)
}

map_ensembl_to_hgnc <- function(ens_ids) {
  ens_ids <- unique(sub("\\..*$", "", ens_ids))

  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ens_ids,
    mart       = mart
  ) %>%
    mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", NA_character_, hgnc_symbol)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
}

# ===================== 1) LOAD ===================== #
files <- list.files(in_dir, pattern = "_limma_AD_vs_Control_allGenes\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No limma TSVs found in: ", in_dir)

message("Found ", length(files), " region files:")
print(basename(files))

all_df <- map_dfr(files, read_limma_tsv)

region_map <- tibble(
  region_raw = c("Mayo_CRB", "Mayo_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC"),
  region     = c("CRB", "TC", "DLPFC", "HCN", "PCC")
)

all_df <- all_df %>%
  left_join(region_map, by = "region_raw") %>%
  mutate(
    region = ifelse(is.na(region), region_raw, region),
    region = factor(region, levels = c("PCC", "HCN", "DLPFC", "TC", "CRB"))
  )

message("\n[CHECK] Genes per region:")
print(all_df %>% count(region, name = "n_genes") %>% arrange(desc(n_genes)))

# ===================== 2) AXES + SIGNIFICANCE ===================== #
all_df <- all_df %>%
  mutate(
    gene_clean = sub("\\..*$", "", gene),

    x_sig      = -log10(pmax(adj.P.Val, P_FLOOR)),
    y_fc       = as.numeric(logFC),

    pass_fdr   = adj.P.Val <= FDR_CUT,
    pass_lfc   = abs(logFC) >= LFC_CUT,
    pass_both  = pass_fdr & pass_lfc,
    direction = case_when(
      pass_both & logFC > 0 ~ "Overexpressed",
      pass_both & logFC < 0 ~ "Underexpressed",
      TRUE ~ "Not significant"
    )
  )

message("\n[CHECK] Significant counts per region (FDR & |logFC|):")
print(all_df %>%
        group_by(region) %>%
        summarise(
          n_genes = n(),
          n_BOTH  = sum(pass_both),
          .groups = "drop"
        ) %>% arrange(desc(n_BOTH)))

# ===================== 3) ANNOTATION (HGNC) ===================== #
map_cache <- file.path(in_dir, "ensembl_to_hgnc_cache.rds")

if (file.exists(map_cache)) {
  message("\nLoading cached HGNC mapping: ", map_cache)
  map_tbl <- readRDS(map_cache)
} else {
  message("\nQuerying biomaRt for HGNC symbols...")
  map_tbl <- map_ensembl_to_hgnc(all_df$gene_clean)
  saveRDS(map_tbl, map_cache)
  message("Saved mapping cache: ", map_cache)
}

all_df <- all_df %>%
  left_join(map_tbl, by = c("gene_clean" = "ensembl_gene_id")) %>%
  mutate(
    has_hgnc = !is.na(hgnc_symbol),
    label    = hgnc_symbol
  )

all_df_write <- all_df %>%
  filter(pass_both) %>%
  dplyr::select(gene, region, logFC, t, adj.P.Val, B, direction, label)
vroom::vroom_write(all_df_write, file.path(in_dir, "DGE_significant_by_region.tsv"))

# ===================== 4) LABELS: TOP N UP + DOWN PER REGION ===================== #
label_df <- all_df %>%
  filter(pass_both, has_hgnc) %>%
  group_by(region) %>%
  group_modify(~{
    up <- .x %>%
      filter(logFC > 0) %>%
      arrange(adj.P.Val, desc(logFC)) %>%
      slice_head(n = TOP_N_PER_SIDE)

    down <- .x %>%
      filter(logFC < 0) %>%
      arrange(adj.P.Val, logFC) %>%
      slice_head(n = TOP_N_PER_SIDE)

    bind_rows(up, down)
  }) %>%
  ungroup()

message("\n[CHECK] Labels per region:")
print(label_df %>% count(region, name = "n_labels") %>% arrange(desc(n_labels)))

label_over  <- label_df %>% filter(direction == "Overexpressed")
label_under <- label_df %>% filter(direction == "Underexpressed")

# ===================== 5) PLOT ===================== #
COL_OVER  <- "#D7301F"
COL_UNDER <- "#2C7FB8"
COL_NS    <- "grey75"

y_cut_line <- -log10(pmax(FDR_CUT, P_FLOOR))

p <- ggplot(all_df, aes(x = y_fc, y = x_sig)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.35, color = "grey35") +
  geom_hline(yintercept = y_cut_line, linetype = "dashed", linewidth = 0.35, color = "grey35") +

  geom_point(
    data = all_df %>% filter(direction == "Not significant"),
    aes(color = direction),
    size = 0.7, alpha = 0.22
  ) +
  geom_point(
    data = all_df %>% filter(direction != "Not significant"),
    aes(color = direction),
    size = 0.85, alpha = 0.50
  ) +

  # anchor points for labeled genes
  geom_point(
    data = label_df,
    aes(color = direction),
    size = 1.6, alpha = 0.90
  ) +

  # overexpressed labels pushed right
  ggrepel::geom_text_repel(
    data = label_over,
    aes(label = label),
    color = "grey10",
    segment.color = COL_OVER,
    show.legend = FALSE,
    size = 4.2,
    direction = "y",
    nudge_x = 0.4,
    force = 5,
    force_pull = 0,
    box.padding = 0.9,
    point.padding = 0.7,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.75,
    segment.size = 0.35,
    max.time = 10,
    max.iter = 200000,
    seed = 2092
  ) +

  # underexpressed labels pushed left
  ggrepel::geom_text_repel(
    data = label_under,
    aes(label = label),
    color = "grey10",
    segment.color = COL_UNDER,
    show.legend = FALSE,
    size = 4.2,
    direction = "y",
    nudge_x = -0.4,
    force = 5,
    force_pull = 0,
    box.padding = 0.9,
    point.padding = 0.7,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.75,
    segment.size = 0.35,
    max.time = 10,
    max.iter = 200000,
    seed = 2
  ) +

  facet_grid(region ~ .) +

  scale_color_manual(
    values = c(
      "Overexpressed"   = COL_OVER,
      "Not significant" = COL_NS,
      "Underexpressed"  = COL_UNDER
    ),
    breaks = c("Overexpressed", "Not significant", "Underexpressed")
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +

  labs(
    x = expression(log[2]~FC),
    y = expression(-log[10](FDR)),
    color = "Expression"
  ) +

  coord_cartesian(clip = "off") +

  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text.y = element_text(face = "bold", size = 18),

    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),

    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold"),

    panel.grid.minor = element_blank(),

    # extra room for labels placed outside the panel
    plot.margin = margin(10, 30, 10, 30)
  )

saveRDS(p, out_rds)

# ===================== 6) SAVE ===================== #
ggsave(out_png, p, width = 10.5, height = 15, dpi = 600)
ggsave(out_pdf, p, width = 10.5, height = 15)
