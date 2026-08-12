#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(biomaRt)
  library(grid)
  library(scales)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG (repo-relative; run from the repository root) ===================== #
source("config/paths.R")

in_base <- file.path(OUTPUTS_DIR, "centrality_comparison")
out_dir <- PLOTS_DIR  # this script's ggsave output is the final Figure 5
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

centralities <- c("degree")

# Region labels used in the scatter panels
networks_scatter <- c("PCC", "DLPFC", "HCN", "CRB", "TC")

# Region names as stored in the FULL/COMMON tables
networks_bar <- c("PCC", "DLPFC", "HCN", "MAYO_CRB", "MAYO_TC")

TOP_PROP <- 0.003  # top 0.3% by |delta percentile rank|

CENT_LABELS <- c(degree = "Degree")
SIGMA <- 1e-8
plog  <- scales::pseudo_log_trans(base = 10, sigma = SIGMA)

safe_breaks <- function(n = 2) {
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(numeric(0))
    rng <- range(x)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) return(rng[1])
    scales::pretty_breaks(n = n)(x)
  }
}

fmt_si <- scales::label_number(
  accuracy = 1,
  scale_cut = scales::cut_si("")
)

cent_pretty <- function(cent_name) {
  out <- unname(CENT_LABELS[cent_name])
  if (is.na(out) || length(out) != 1) return(cent_name)
  out
}

# ===================== HELPERS (scatter inputs) ===================== #
read_common_tbl <- function(net, cent) {
  file_net <- dplyr::recode(net, CRB = "MAYO_CRB", TC = "MAYO_TC", .default = net)
  std_net  <- net

  f <- file.path(in_base, cent, paste0(file_net, "_", cent, "_COMMON.tsv"))
  if (!file.exists(f)) return(NULL)

  df <- vroom::vroom(f, show_col_types = FALSE)

  need <- c("gene", "AD", "CTL", "perc_AD", "perc_CTL", "delta_perc", "HigherNetwork")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse = ", "))

  df %>%
    mutate(network = std_net, centrality = cent) %>%
    dplyr::select(network, centrality, everything())
}

read_spearman_tbl <- function(net) {
  file_net <- dplyr::recode(net, CRB = "MAYO_CRB", TC = "MAYO_TC", .default = net)
  std_net  <- net

  f <- file.path(in_base, "spearman", paste0(file_net, "_spearman.tsv"))
  if (!file.exists(f)) return(NULL)

  df <- vroom::vroom(f, show_col_types = FALSE)

  need <- c("network", "centrality", "n_common_genes", "spearman_rho", "spearman_pvalue")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse = ", "))

  df %>%
    mutate(
      network    = std_net,
      centrality = as.character(centrality)
    ) %>%
    dplyr::select(network, centrality, n_common_genes, spearman_rho, spearman_pvalue)
}

map_ensembl_to_hgnc <- function(ens_ids) {
  ens_ids <- unique(ens_ids)
  ens_ids <- sub("\\..*$", "", ens_ids)

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

# ===================== LOAD COMMON TABLES (degree) ===================== #
all_common <- purrr::map_dfr(
  networks_scatter,
  \(net) purrr::map_dfr(centralities, \(cent) read_common_tbl(net, cent))
)
if (nrow(all_common) == 0) stop("No *_COMMON.tsv files were loaded. Check paths and filenames.")

# ===================== BIOMART ANNOTATION ===================== #
cat("Annotating ENSG -> HGNC with biomaRt...\n")
ann <- map_ensembl_to_hgnc(all_common$gene)

all_common <- all_common %>%
  mutate(gene_key = sub("\\..*$", "", gene)) %>%
  left_join(ann, by = c("gene_key" = "ensembl_gene_id")) %>%
  mutate(gene_label = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, NA_character_))

# ===================== TOP GENES PER PANEL (by |delta percentile|) ===================== #
top_tbl <- all_common %>%
  group_by(network, centrality) %>%
  mutate(abs_dperc = abs(delta_perc)) %>%
  slice_max(order_by = abs_dperc, prop = TOP_PROP, with_ties = FALSE) %>%
  ungroup()

# ===================== LOAD SPEARMAN RESULTS ===================== #
spearman_all <- purrr::map_dfr(networks_scatter, read_spearman_tbl)
if (nrow(spearman_all) == 0) stop("No spearman/<NET>_spearman.tsv files were loaded. Check paths and filenames.")

# ===================== RHO LABEL POSITION (robust under transforms) ===================== #
rho_pos <- all_common %>%
  group_by(network, centrality) %>%
  summarise(
    x_raw_min = suppressWarnings(min(AD,  na.rm = TRUE)),
    x_raw_max = suppressWarnings(max(AD,  na.rm = TRUE)),
    y_raw_min = suppressWarnings(min(CTL, na.rm = TRUE)),
    y_raw_max = suppressWarnings(max(CTL, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(spearman_all, by = c("network", "centrality")) %>%
  mutate(
    x_raw_min = ifelse(is.finite(x_raw_min), x_raw_min, 0),
    x_raw_max = ifelse(is.finite(x_raw_max), x_raw_max, x_raw_min),
    y_raw_min = ifelse(is.finite(y_raw_min), y_raw_min, 0),
    y_raw_max = ifelse(is.finite(y_raw_max), y_raw_max, y_raw_min),

    x_raw_max = ifelse(x_raw_max == x_raw_min, x_raw_max + 1e-12, x_raw_max),
    y_raw_max = ifelse(y_raw_max == y_raw_min, y_raw_max + 1e-12, y_raw_max),

    rho_abs   = abs(spearman_rho),
    rho_label = sprintf("rho = %.3f", spearman_rho),

    x_tr_min = plog$transform(x_raw_min),
    x_tr_max = plog$transform(x_raw_max),
    y_tr_min = plog$transform(y_raw_min),
    y_tr_max = plog$transform(y_raw_max),

    x_tr = x_tr_min + 0.92 * (x_tr_max - x_tr_min),
    y_tr = y_tr_min + 0.97 * (y_tr_max - y_tr_min),

    x = plog$inverse(x_tr),
    y = plog$inverse(y_tr)
  )

# ===================== BARPLOT (degree) -> first panel ===================== #
make_barplot_degree <- function() {

  NET_LABELS_BAR <- c(
    PCC      = "PCC",
    DLPFC    = "DLPFC",
    HCN      = "HCN",
    MAYO_CRB = "CRB",
    MAYO_TC  = "TC"
  )

  FILL_COLORS <- c(
    "Unique Control" = "#377EB8",
    "Shared"         = "#CFCFCF",
    "Unique AD"      = "#E41A1C"
  )

  BAR_ALPHA <- 0.78
  LABEL_MIN_PROP <- 0.035

  read_full_common <- function(net, cent) {
    f_full   <- file.path(in_base, cent, paste0(net, "_", cent, "_FULL.tsv"))
    f_common <- file.path(in_base, cent, paste0(net, "_", cent, "_COMMON.tsv"))
    if (!file.exists(f_full))   stop("Missing: ", f_full)
    if (!file.exists(f_common)) stop("Missing: ", f_common)

    full   <- readr::read_tsv(f_full,   show_col_types = FALSE)
    common <- readr::read_tsv(f_common, show_col_types = FALSE)

    n_full   <- nrow(full)
    n_common <- nrow(common)

    n_unique_ad  <- sum(!is.na(full$AD)  &  is.na(full$CTL))
    n_unique_ctl <- sum( is.na(full$AD)  & !is.na(full$CTL))

    if ((n_unique_ad + n_unique_ctl + n_common) != n_full) {
      warning("Counts don't add up for ", net)
    }

    tibble(
      network = net,
      n_full = n_full,
      n_common = n_common,
      n_unique_ctl = n_unique_ctl,
      n_unique_ad = n_unique_ad
    )
  }

  summary_tbl <- bind_rows(lapply(networks_bar, read_full_common, cent = "degree")) %>%
    mutate(network_label = unname(NET_LABELS_BAR[network]))

  plot_df <- summary_tbl %>%
    dplyr::select(network_label, n_common, n_unique_ctl, n_unique_ad) %>%
    tidyr::pivot_longer(
      cols = c(n_unique_ctl, n_common, n_unique_ad),
      names_to = "category",
      values_to = "n_genes"
    ) %>%
    mutate(
      category = dplyr::recode(
        category,
        n_unique_ctl = "Unique Control",
        n_common     = "Shared",
        n_unique_ad  = "Unique AD"
      ),
      category = factor(category, levels = c("Unique Control", "Shared", "Unique AD")),
      network_label = factor(network_label, levels = unname(NET_LABELS_BAR[networks_bar]))
    ) %>%
    group_by(network_label) %>%
    mutate(
      prop = n_genes / sum(n_genes),
      label_pct = ifelse(prop >= LABEL_MIN_PROP, scales::percent(prop, accuracy = 0.1), "")
    ) %>%
    ungroup()

  ggplot(plot_df, aes(x = network_label, y = prop, fill = category)) +
    geom_col(
      width = 0.72,
      color = "white",
      linewidth = 0.6,
      alpha = BAR_ALPHA
    ) +
    geom_text(
      aes(label = label_pct),
      position = position_stack(vjust = 0.5),
      size = 4,
      color = "black"
    ) +
    scale_fill_manual(values = FILL_COLORS, drop = FALSE) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = "Gene proportion",
      fill = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

p_bar_degree <- make_barplot_degree()

# ===================== SCATTER (per region) ===================== #
make_scatter_one_network <- function(cent_name, net_name) {

  df_cent  <- all_common %>% filter(centrality == cent_name, network == net_name)
  top_cent <- top_tbl     %>% filter(centrality == cent_name, network == net_name)
  rho_cent <- rho_pos     %>% filter(centrality == cent_name, network == net_name)

  top_cent_lab <- top_cent %>% filter(!is.na(gene_label) & gene_label != "")
  cent_title <- cent_pretty(cent_name)

  ggplot(df_cent, aes(x = AD, y = CTL)) +
    geom_point(aes(color = HigherNetwork), alpha = 0.22, size = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_point(data = top_cent, aes(color = HigherNetwork), alpha = 0.95, size = 2.8) +
    ggrepel::geom_label_repel(
      data = top_cent_lab,
      aes(label = gene_label),
      inherit.aes = TRUE,
      fill = "white",
      label.size = 0.25,
      label.r = unit(0.15, "lines"),
      label.padding = unit(0.22, "lines"),
      color = "black",
      fontface = "plain",
      size = 3.3,
      box.padding = 0.75,
      point.padding = 0.55,
      segment.alpha = 0.55,
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 3.0,
      force_pull = 0.20,
      seed = 1
    ) +
    geom_label(
      data = rho_cent,
      aes(x = x, y = y, label = rho_label, fill = rho_abs),
      inherit.aes = FALSE,
      label.size = 0.25,
      label.r = unit(0.20, "lines"),
      label.padding = unit(0.28, "lines"),
      alpha = 0.98,
      size = 3.9,
      color = "black",
      fontface = "bold"
    ) +
    scale_color_manual(values = c("AD" = "#E41A1C", "Control" = "#377EB8")) +
    scale_fill_gradient(
      limits = c(0, 1),
      low  = "#C9A227",
      high = "#FFF7CC",
      guide = "none"
    ) +
    labs(
      title = net_name,
      x = paste0(cent_title, " (AD network)"),
      y = paste0(cent_title, " (Control network)"),
      color = paste0("Higher\n", cent_title)
    ) +
    scale_x_continuous(trans = plog, breaks = safe_breaks(n = 2), labels = fmt_si) +
    scale_y_continuous(trans = plog, breaks = safe_breaks(n = 2), labels = fmt_si) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.text.y = element_text(margin = margin(r = 4))
    )
}

# ===================== GRID (barplot + 5 scatters) ===================== #
make_grid_degree <- function() {
  cent_name <- "degree"
  cent_title <- cent_pretty(cent_name)

  scat_list <- lapply(networks_scatter, \(nn) make_scatter_one_network(cent_name, nn) +
                        theme(legend.position = "none"))

  p_bar_panel <- p_bar_degree + labs(title = "Proportion of genes by network")

  wrap_plots(
    c(list(p_bar_panel), scat_list),
    ncol = 3
  ) +
    plot_annotation(
      title = paste0(cent_title, " centrality: AD vs control")
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
}

# ===================== RUN + SAVE ===================== #
p <- make_grid_degree()

ggsave(
  filename = file.path(out_dir, "Figure5_degree_AD_vs_control_panels.pdf"),
  plot = p, width = 16, height = 8
)
ggsave(
  filename = file.path(out_dir, "Figure5_degree_AD_vs_control_panels.png"),
  plot = p, width = 16, height = 8, dpi = 300
)

cat("\nDONE\n")
