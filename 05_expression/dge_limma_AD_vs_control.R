#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(limma)
  library(vroom)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG (repo-relative; run from the repository root) ===================== #
source("config/paths.R")

in_dir  <- file.path(INPUTS_DIR, "expression_matrices")
out_dir <- file.path(OUTPUTS_DIR, "dge_limma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- readRDS(file.path(OUTPUTS_DIR, "metadata_sex.rds"))

# Reporting thresholds (used only for summary; full tables are always saved)
PADJ_CUT   <- 0.05
LFC_CUT    <- 0.5
LOG_OFFSET <- 0.5  # avoids log2(0) = -Inf

# ===================== HELPERS ===================== #

# Read matrix: first column must be 'gene', remaining columns are samples
read_expr <- function(path) {
  df <- read_tsv(path, col_types = cols(), progress = FALSE)
  if (colnames(df)[1] != "gene") stop("Expected first column 'gene' in: ", path)
  mat <- df %>% dplyr::select(-gene) %>% as.matrix()
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- df$gene
  mat
}

print_qc <- function(mat, label) {
  x <- as.numeric(mat)
  message(sprintf(
    "%s | genes=%d samples=%d | dupGenes=%d | NA=%s Inf=%s | min=%.3g med=%.3g max=%.3g",
    label, nrow(mat), ncol(mat),
    sum(duplicated(rownames(mat))),
    any(is.na(x)), any(is.infinite(x)),
    min(x, na.rm = TRUE), median(x, na.rm = TRUE), max(x, na.rm = TRUE)
  ))
}

summarize_de <- function(tt, contrast_name, padj = PADJ_CUT, lfc = LFC_CUT) {
  sig <- tt %>% filter(adj.P.Val <= padj, abs(logFC) >= lfc)
  tibble(
    Contrast      = contrast_name,
    Significant   = nrow(sig),
    Upregulated   = sum(sig$logFC > 0),
    Downregulated = sum(sig$logFC < 0)
  )
}

# ===================== CORE ===================== #

run_limma_pair <- function(file_ad, file_ctl, prefix) {

  message("\n==============================")
  message("Contrast: ", prefix)
  message("  AD file:      ", basename(file_ad))
  message("  Control file: ", basename(file_ctl))

  # 1) Load preprocessed matrices (already normalized + corrected)
  mat_ad  <- read_expr(file_ad)
  mat_ctl <- read_expr(file_ctl)

  # 2) Input QC
  print_qc(mat_ad,  "QC AD")
  print_qc(mat_ctl, "QC Control")

  # 3) Harmonize genes across paired files
  common_genes <- intersect(rownames(mat_ad), rownames(mat_ctl))
  message(sprintf("Gene overlap: common=%d | AD=%d | Control=%d",
                  length(common_genes), nrow(mat_ad), nrow(mat_ctl)))
  if (length(common_genes) == 0) stop("No overlapping genes for: ", prefix)

  mat_ad  <- mat_ad[common_genes, , drop = FALSE]
  mat_ctl <- mat_ctl[common_genes, , drop = FALSE]

  # 4) Combine expression and build design (two-group comparison + sex)
  expr <- cbind(mat_ad, mat_ctl)
  sample_ids <- colnames(expr)

  meta_sub <- metadata %>%
    filter(specimenID %in% sample_ids)

  meta_sub <- meta_sub[match(sample_ids, meta_sub$specimenID), ]
  if (any(is.na(meta_sub$specimenID))) {
    stop("Mismatch between expression matrix and metadata (specimenID)")
  }
  sex <- factor(meta_sub$sex, levels = c("female", "male"))
  if (any(is.na(sex))) {
    stop("Missing sex annotation in metadata")
  }
  group <- factor(c(rep("AD", ncol(mat_ad)), rep("Control", ncol(mat_ctl))),
                  levels = c("Control", "AD"))
  message(sprintf("Samples: AD=%d | Control=%d", sum(group == "AD"), sum(group == "Control")))

  design <- model.matrix(~ 0 + group + sex)
  colnames(design) <- gsub("^group", "", colnames(design))
  message("Design columns: ", paste(colnames(design), collapse = ", "))

  # 5) Log2 transform for limma-trend (stabilizes variance; offset avoids -Inf)
  log_expr <- log2(expr + LOG_OFFSET)
  message(sprintf("Log QC: NA=%s Inf=%s | min=%.3g med=%.3g max=%.3g",
                  any(is.na(log_expr)), any(is.infinite(log_expr)),
                  min(log_expr), median(log_expr), max(log_expr)))

  # 6) Fit limma model + contrast (AD vs Control) with robust eBayes
  fit  <- lmFit(log_expr, design)
  contr <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  # 7) Extract all genes (no filtering here), then save
  tt <- topTable(fit2, coef = "AD_vs_Control", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- tt %>% relocate(gene, .before = everything())

  out_tsv <- file.path(out_dir, paste0(prefix, "_limma_AD_vs_Control_allGenes.tsv"))
  vroom_write(tt, out_tsv, delim = "\t")
  message("Saved: ", out_tsv)

  # 8) DEG count summary (based on PADJ_CUT and LFC_CUT)
  sum_tbl <- summarize_de(tt, prefix)
  print(sum_tbl)

  list(prefix = prefix, out_tsv = out_tsv, summary = sum_tbl)
}

# ===================== DISCOVER FILE PAIRS ===================== #

ad_files  <- list.files(in_dir, pattern = "_counts_AD\\.tsv$", full.names = TRUE)
ctl_files <- list.files(in_dir, pattern = "_counts_control\\.tsv$", full.names = TRUE)
if (length(ad_files) == 0 || length(ctl_files) == 0) stop("No AD/control TSVs found in: ", in_dir)

key <- function(p) {
  bn <- basename(p)
  bn <- str_replace(bn, "_counts_AD\\.tsv$", "")
  bn <- str_replace(bn, "_counts_control\\.tsv$", "")
  bn
}

pairs <- inner_join(
  tibble(ad = ad_files,  key = map_chr(ad_files, key)),
  tibble(ct = ctl_files, key = map_chr(ctl_files, key)),
  by = "key"
) %>% arrange(key)

message("Pairs detected: ", nrow(pairs))
print(pairs %>% dplyr::select(key, ad, ct))

# ===================== RUN ALL ===================== #

results <- pmap(
  list(pairs$ad, pairs$ct, pairs$key),
  ~ run_limma_pair(file_ad = ..1, file_ctl = ..2, prefix = ..3)
)

summary_df <- bind_rows(lapply(results, `[[`, "summary")) %>%
  arrange(desc(Significant))

out_summary <- file.path(out_dir, "limma_DE_summary_AD_vs_Control.tsv")
vroom_write(summary_df, out_summary, delim = "\t")

message("\nDONE")
message("Summary saved: ", out_summary)
print(summary_df)
