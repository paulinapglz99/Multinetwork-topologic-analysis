#!/usr/bin/env Rscript

# compare_networks_jaccard.R
# Calcula matrices de similitud de Jaccard entre módulos
# Empareja automáticamente AD vs Control y genera heatmaps

# Usage:
# Rscript compare_networks_jaccard.R \
#    --input_dir /ruta/a/archivos \
#    --output_dir results_jaccard

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  data.table, tidyverse, tools, optparse, ggplot2
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages did not load correctly: ",
       paste(names(ok)[!ok], collapse = ", "))
}

# Parser
option_list <- list(make_option(c("-i", "--input_dir"), type = "character", default = getwd(), help = "Directory where the *_nodes_summary.csv are", metavar = "character"),
                    make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard", help = "Output directory [default: %default]", metavar = "character"),
                    make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$", help = "Regex to select files [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input_dir
out_dir   <- opt$output_dir
pattern   <- opt$pattern

message("Input dir: ", input_dir)
message("Output dir: ", out_dir)
message("Pattern: ", pattern)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# --- FUNCTIONS --------------------------------------------------------------

# Jaccard index
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Jaccard matrix between two lists of modules
jaccard_matrix <- function(redA, redB) {
  mat <- matrix(0, nrow = length(redA), ncol = length(redB),
                dimnames = list(names(redA), names(redB)))
  for (i in seq_along(redA)) {
    for (j in seq_along(redB)) {
      mat[i, j] <- jaccard_index(redA[[i]], redB[[j]])
    }
  }
  return(mat)
}

# Make readable names for comparisons (robusto: varios sufijos/variantes)
make_pretty_name <- function(x) {
  x %>%
    str_replace_all("(?i)_top10pct_nodes_summary", "") %>%
    str_replace_all("(?i)_nodes_summary", "") %>%
    str_replace_all("(?i)_counts_", " ") %>%
    str_replace_all("(?i)_counts", " ") %>%
    str_replace_all("_vs_", " vs ") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

# Compare networks given a named list and list of pairs (each pair is c("A","B"))
compare_networks <- function(list_nets, pairs) {
  results <- list()
  for (p in pairs) {
    netA_name <- p[1]
    netB_name <- p[2]
    if (!(netA_name %in% names(list_nets))) stop("Missing network: ", netA_name)
    if (!(netB_name %in% names(list_nets))) stop("Missing network: ", netB_name)
    netA <- list_nets[[netA_name]]
    netB <- list_nets[[netB_name]]
    name <- paste(netA_name, netB_name, sep = "_vs_")
    message("Comparing ", name)
    results[[name]] <- jaccard_matrix(netA, netB)
  }
  return(results)
}

# --- READ FILES / BUILD MODULE LISTS ----------------------------------------

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2) stop("You need at least 2 files *_nodes_summary.csv to compare.")

modules_nets <- list()

for (f in files) {
  net_name <- tools::file_path_sans_ext(basename(f))
  df <- data.table::fread(f)
  if (!all(c("node", "membership_infomap") %in% colnames(df))) {
    stop("File ", f, " must contain columns 'node' and 'membership_infomap'")
  }
  modules_nets[[net_name]] <- split(df$node, df$membership_infomap)
}

# --- AUTO-DETECT PAIRS (AD vs CONTROL) -------------------------------------

names_df <- tibble(name = names(modules_nets)) %>%
  mutate(
    # prefix: everything up to "_counts" (lazy). If no "_counts", take part before _AD_/_control_ as fallback
    prefix = coalesce(
      str_extract(name, "^.*?_counts"),
      str_remove(name, "(?i)(_AD_|_ALZ|_ALZHEIMER|_CONTROL_|_CTRL_).*")
    ) %>% str_replace_all("(^_|_$)", "") ,
    condition = case_when(
      str_detect(name, regex("(_AD_|\\bAD\\b|_ALZ_|ALZHEIMER)", ignore_case = TRUE)) ~ "AD",
      str_detect(name, regex("(_CONTROL_|\\bCONTROL\\b|_CTRL_)", ignore_case = TRUE)) ~ "control",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(condition)

pairs_auto <- names_df %>%
  group_by(prefix) %>%
  filter(any(condition == "AD") & any(condition == "control")) %>%
  summarise(
    ad = name[condition == "AD"][1],
    ctrl = name[condition == "control"][1],
    .groups = "drop"
  ) %>%
  mutate(pair = map2(ad, ctrl, ~c(.x, .y))) %>%
  pull(pair)

if (length(pairs_auto) == 0) stop("No pairs detected automatically. Revisa los nombres de archivos y patrones AD/control.")

message("Pairs detected (AD vs control):")
print(pairs_auto)

# --- RUN COMPARISONS --------------------------------------------------------

results <- compare_networks(modules_nets, pairs_auto)

# Save jaccard matrix results
for (n in names(results)) {
  out_file <- file.path(out_dir, paste0("Jaccard_", n, ".csv"))
  write.csv(results[[n]], out_file, row.names = TRUE)
  message("Saved: ", out_file)
}

# --- SUMMARY DE PERFECT MATCHES (Jaccard == 1) -------------------------------

perfect_summary <- tibble()

for (n in names(results)) {
  mat <- results[[n]]
  perfect_count <- sum(mat == 1, na.rm = TRUE)
  total_pairs <- length(mat)                 # m * n
  perc <- round(100 * perfect_count / total_pairs, 2)
  pretty_name <- make_pretty_name(n)
  perfect_summary <- bind_rows(
    perfect_summary,
    tibble(
      comparison = pretty_name,
      perfect_pairs = perfect_count,
      total_pairs = total_pairs,
      percentage = perc
    )
  )
}

out_csv <- file.path(out_dir, "tabla_similitud_perfecta.csv")
write.csv(perfect_summary, out_csv, row.names = FALSE)
message("Saved summary: ", out_csv)
print(perfect_summary)

# --- PLOTS (HEATMAPS) ------------------------------------------------------

message("Plotting heatmaps")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Create pretty names vector aligned with names(results)
pretty_names <- sapply(names(results), make_pretty_name, USE.NAMES = FALSE)

# Loop to plot heatmaps
for (i in seq_along(results)) {
  n <- names(results)[i]
  title_txt <- pretty_names[i]
  
  mat <- results[[n]]
  df_plot <- as.data.frame(as.table(mat))
  colnames(df_plot) <- c("Network_1", "Network_2", "Jaccard")
  
  p <- ggplot(df_plot, aes(Network_1, Network_2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.5) + # if ggplot2 < 3.4, use size = 0.5
    scale_fill_gradient(low = "white", high = "navyblue") +
    labs(
      title = title_txt,
      subtitle = "Jaccard similarity index between modules",
      fill = "Jaccard",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      legend.position = "right",
      legend.key.height = unit(1, "cm")
    )
  
  out_png <- file.path(plot_dir, paste0("Heatmap_", n, ".png"))
  ggsave(out_png, plot = p, width = 15, height = 10, dpi = 300)
  message("Saved plot: ", out_png)
}

message("Done.")
# END