#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(clusterProfiler)
  library(msigdbr)
})

# ===================== CONFIG ===================== #
in_dir  <- "results/dge_limma"
out_dir <- "results/gsea_reactome"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pattern_in <- "_limma_AD_vs_Control_allGenes\\.tsv$"

PADJ_CUT   <- 0.05
TOP_N_SIDE <- 4  # top terms per direction (up/down) per region

# ===================== GENE SETS ===================== #
TERM2GENE_REACT <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") |>
  transmute(term = gs_name, gene = db_ensembl_gene) |>
  filter(!is.na(gene), gene != "") |>
  distinct()

# ===================== GSEA PER REGION ===================== #
run_one <- function(path) {
  region <- str_replace(basename(path), pattern_in, "")
  df <- read_tsv(path, col_types = cols(), progress = FALSE) |>
    transmute(gene = sub("\\..*$", "", gene), stat = as.numeric(t)) |>
    filter(!is.na(gene), gene != "", is.finite(stat)) |>
    group_by(gene) |>
    summarise(stat = stat[which.max(abs(stat))], .groups = "drop")

  geneList <- df$stat; names(geneList) <- df$gene
  geneList <- sort(geneList, decreasing = TRUE)

  as_tibble(as.data.frame(
    suppressWarnings(GSEA(
      geneList     = geneList,
      TERM2GENE    = TERM2GENE_REACT,
      pvalueCutoff = 1,
      minGSSize    = 10,
      maxGSSize    = 500,
      eps          = 0,
      verbose      = FALSE
    ))
  )) |>
    mutate(region = region, .before = 1)
}

files <- list.files(in_dir, pattern = pattern_in, full.names = TRUE)
if (length(files) == 0) stop("No limma TSVs found in: ", in_dir)

gsea_full <- map_dfr(files, run_one)

write_tsv(gsea_full, file.path(out_dir, "GSEA_reactome_full_terms.tsv"))

# ===================== FILTER + TOP TERMS ===================== #
gsea_sig <- gsea_full %>%
  mutate(
    region = str_remove(region, "^(Mayo_|ROSMAP_)"),
    term_simple = str_replace_all(str_remove(Description, "^REACTOME_"), "_", " ")
  ) %>%
  filter(p.adjust < PADJ_CUT)

write_tsv(gsea_sig, file.path(out_dir, "GSEA_reactome_significant_terms.tsv"))

top_terms_by_region <- gsea_sig %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
  group_by(region, direction) %>%
  arrange(region, direction, desc(abs(NES))) %>%
  slice_head(n = TOP_N_SIDE) %>%
  ungroup()

write_tsv(
  top_terms_by_region,
  file.path(out_dir, "GSEA_top_terms_by_region.tsv")
)

message("Significant terms: ", nrow(gsea_sig),
        " | unique top terms: ", length(unique(top_terms_by_region$term_simple)))
