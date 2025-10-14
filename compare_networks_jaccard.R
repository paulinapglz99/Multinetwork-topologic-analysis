#!/usr/bin/env Rscript

#compare_networks_jaccard.R
#alcula matrices de similitud de Jaccard entre módulos
# Empareja automáticamente AD vs Control y genera heatmaps

#Usage:
#Rscript compare_networks_modules_auto.R \
#    --input_dir /ruta/a/archivos \
#    --output_dir results_jaccard

#Packages

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

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directory where the *_nodes_summary.csv, output from the network_topology.R script", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard",
              help = "Output directory [default: %default]", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex to select archives [default: %default]", metavar = "character")
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

#Functions --- ---

#Jaccard function
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#Jaccard matrix
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

#Helper function to create legible names of net pairs
make_pretty_name <- function(x) {
  x %>% 
    str_replace_all("_counts_", " ") %>% 
    str_replace_all("_top10pct_nodes_summary", "")  %>% 
    str_replace_all("_vs_", " vs ")  %>% 
    str_replace_all("_", " ")  %>% 
    str_replace_all("\\s+", " ") %>% 
    str_trim()
}

#Compare networks with jaccard

compare_networks <- function(list_nets, pairs) {
  results <- list()
  for (p in pairs) {
    netA_name <- p[1]
    netB_name <- p[2]
    netA <- list_nets[[netA_name]]
    netB <- list_nets[[netB_name]]
    name <- paste(netA_name, netB_name, sep = "_vs_")
    message("Comparing ", name)
    results[[name]] <- jaccard_matrix(netA, netB)
  }
  return(results)
}

#Get data from modules --- ---

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2)
  stop("You need at least 2 files *_nodes_summary.csv to compare.")

modules_nets <- list()

for (f in files) {
  net_name <- tools::file_path_sans_ext(basename(f))
  df <- data.table::fread(f)
  if (!all(c("node", "membership_infomap") %in% colnames(df))) {
    stop("File ", f, " must contain columns 'node' and 'membership_infomap'")
  }
  modules_nets[[net_name]] <- split(df$node, df$membership_infomap)
}

#Check pairs --- ---

names_df <- tibble(name = names(modules_nets)) %>%
  mutate(
    prefix = str_extract(name, "^[A-Za-z0-9_]+_counts"),
    condition = case_when(
      str_detect(name, "_AD_") ~ "AD",
      str_detect(name, "_control_") ~ "control",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(condition)

pairs_auto <- names_df %>%
  group_by(prefix) %>%
  filter(all(c("AD", "control") %in% condition)) %>%
  summarise(
    ad = name[condition == "AD"][1],
    ctrl = name[condition == "control"][1],
    .groups = "drop"
  ) %>%
  mutate(pair = map2(ad, ctrl, ~c(.x, .y))) %>%
  pull(pair)

if (length(pairs_auto) == 0)  stop("No pairs detected automatically")

print(pairs_auto)

#Run results

results <- compare_networks(modules_nets, pairs_auto)

#Save jaccard matrix results

for (n in names(results)) {
  out_file <- file.path(out_dir, paste0("Jaccard_", n, ".csv"))
  write.csv(results[[n]], out_file, row.names = TRUE)
  message("Saved: ", out_file)
}

#Check similitudes - Perfect (Jaccard = 1)

perfect_summary <- tibble()

for (n in names(results)) {
  mat <- results[[n]]
  
  #How many pairs of modules have Jaccard = 1
  perfect_count <- sum(mat == 1, na.rm = TRUE)
  #Total number of possible pairs
  total_pairs <- length(mat)
  #Percentage of perfect matches
  perc <- round(100 * perfect_count / total_pairs, 2)
  #Nice name for the comparison
  pretty_name <- make_pretty_name(n)
  #Make it a df
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

#Save summary of perfect pairs
out_csv <- file.path(out_dir, "tabla_similitud_perfecta.csv")
write.csv(perfect_summary, out_csv, row.names = FALSE)
print(perfect_summary)

#Plot heatmaps
message("Plotting heatmaps")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

pretty_names <- make_pretty_name(n)

#Loop to plot heatmaps
for (i in seq_along(results)) {
  n <- names(results)[i]
  title_txt <- pretty_names[i]
  
  mat <- results[[n]]
  df_plot <- as.data.frame(as.table(mat))
  colnames(df_plot) <- c("Network_1", "Network_2", "Jaccard")
  
  p <- ggplot(df_plot, aes(Network_1, Network_2, fill = Jaccard)) +
    geom_tile(color = "gray90") +
    scale_fill_gradient(low = "white", high = "navyblue") +
    theme_minimal(base_size = 12) +
    labs(
      title = title_txt,
      subtitle = "Jaccard similarity index between modules",
      fill = " ", x = " ", y = " "
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 13)
    )
  
  out_png <- file.path(plot_dir, paste0("Heatmap_", n, ".png"))
  ggsave(out_png, plot = p, width = 8, height = 6, dpi = 300)
  message("Saved plot: ", out_png)
}

#END