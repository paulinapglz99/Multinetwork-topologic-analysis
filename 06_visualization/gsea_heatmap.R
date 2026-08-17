#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

# ===================== CONFIG (repo-relative; run from the repository root) ===================== #
# NOTE: intermediate component consumed by dge_gsea_panel.R (Figure 6), stays under outputs/.
source("config/paths.R")

in_tsv  <- file.path(OUTPUTS_DIR, "gsea_reactome", "GSEA_top_terms_by_region.tsv")
out_dir <- file.path(OUTPUTS_DIR, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

REG_ORDER <- c("CRB", "TC", "DLPFC", "HCN", "PCC")

# ===================== LOAD ===================== #
top_terms <- read_tsv(in_tsv, col_types = cols(), progress = FALSE)

plot_df <- top_terms %>%
  group_by(term_simple) %>%
  mutate(absNES_max = max(abs(NES), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    region      = factor(region, levels = REG_ORDER),
    sig         = -log10(p.adjust),
    term_simple = reorder(term_simple, absNES_max),

    star = dplyr::case_when(
      p.adjust <= 1e-3 ~ "***",
      p.adjust <= 1e-2 ~ "**",
      p.adjust <= 5e-2 ~ "*",
      TRUE             ~ ""
    ),

    # continuous alpha to soften weaker terms without hiding NES
    alpha_sig = scales::rescale(sig, to = c(0.25, 1.0), from = range(sig, finite = TRUE))
  )

# ===================== PLOT ===================== #
p <- ggplot(plot_df, aes(x = 1, y = term_simple)) +
  geom_tile(
    aes(fill = NES, alpha = alpha_sig),
    color = "grey92", linewidth = 0.4
  ) +
  geom_text(
    aes(label = star),
    fontface = "bold", size = 4, color = "grey10"
  ) +
  facet_grid(. ~ region, switch = "x") +
  scale_x_continuous(
    NULL,
    breaks = NULL,
    limits = c(0.5, 1.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_gradient2(
    name = "NES",
    low  = "#2C7FB8",
    mid  = "#F7F7F7",
    high = "#D7301F",
    midpoint = 0
  ) +
  scale_alpha(guide = "none") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),

    axis.text.y  = element_text(size = 10, face = "bold"),

    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold"),

    panel.spacing.x  = unit(0.15, "lines"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 14, face = "bold")
  )

p <- p + theme(plot.margin = margin(t = 8, r = 25, b = 8, l = 8))

# ===================== SAVE ===================== #
saveRDS(p, file.path(out_dir, "gsea_heatmap.rds"))

ggsave(
  file.path(out_dir, "gsea_heatmap.png"),
  p, width = 16, height = 9, units = "in", dpi = 600, limitsize = FALSE)

ggsave(
  file.path(out_dir, "gsea_heatmap.pdf"),
  p, width = 16, height = 9, units = "in",
  limitsize = FALSE, device = cairo_pdf)
