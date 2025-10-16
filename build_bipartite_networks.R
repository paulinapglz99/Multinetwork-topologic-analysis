#!/usr/bin/env Rscript

# build_bipartite_networks.R
# Construye redes bipartitas a partir de matrices de Jaccard generadas por compare_networks_jaccard.R
# Cada red bipartita representa comunidades (módulos) AD y Control conectadas por su similitud Jaccard.

# Uso:
# Rscript build_bipartite_networks.R \
#   --input_dir results_jaccard \
#   --output_dir results_bipartite \
#   --threshold 0.2

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  tidyverse, igraph, ggraph, data.table, optparse
)

if (all(ok)) {
  message("All packages loaded correctly.")
} else {
  stop("Some packages did not load correctly.")
}

# --- PARSER -----------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directory with Jaccard_*.csv matrices", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_bipartite",
              help = "Output directory for bipartite graphs", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.2,
              help = "Minimum Jaccard similarity to include an edge [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input_dir
input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_jaccard_1/"
output_dir <- opt$output_dir
output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/bipartite_nets"
threshold <- opt$threshold
threshold <-0.25

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Input dir: ", input_dir)
message("Output dir: ", output_dir)
message("Threshold: ", threshold)

# ---------------------------------------------------------------------------
# Función auxiliar para nombres bonitos
# ---------------------------------------------------------------------------
make_pretty_name <- function(x) {
  x %>%
    str_replace_all("(?i)Jaccard_", "") %>%
    str_replace_all("_vs_", " vs ") %>%
    str_replace_all("_", " ") %>%
    str_trim()
}

# ---------------------------------------------------------------------------
# Lectura de matrices Jaccard y construcción de redes bipartitas
# ---------------------------------------------------------------------------
files <- list.files(input_dir, pattern = "^Jaccard_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No se encontraron matrices Jaccard_*.csv en ", input_dir)

for (f in files) {
  message("\nProcesando: ", basename(f))
  
  # Leer la matriz (primera columna = nombres de fila)
  df <- vroom::vroom(f)
  first_col <- colnames(df)[1]
  
  mat <- df %>%
    rename(Row = !!first_col) %>%
    column_to_rownames("Row") %>%
    as.matrix()
  
  # Convertir a formato largo
  df_edges <- as.data.frame(as.table(mat))
  colnames(df_edges) <- c("AD_module", "Control_module", "Jaccard")
  
  df_edges <- df_edges %>%
    filter(Jaccard >= threshold)
  
  if (nrow(df_edges) == 0) {
    message("⚠️  No edges above threshold for ", basename(f))
    next
  }
  
  # -------------------------------------------------------------------------
  # Crear grafo bipartito
  # -------------------------------------------------------------------------
  nodes_AD <- unique(df_edges$AD_module)
  nodes_Control <- unique(df_edges$Control_module)
  
  nodes <- tibble(
    name = c(nodes_AD, nodes_Control),
    type = c(rep("AD", length(nodes_AD)), rep("Control", length(nodes_Control)))
  )
  
  g <- graph_from_data_frame(df_edges, vertices = nodes, directed = FALSE)
  
  # -------------------------------------------------------------------------
  # Guardar el grafo
  # -------------------------------------------------------------------------
  region_name <- tools::file_path_sans_ext(basename(f))
  out_graph <- file.path(output_dir, paste0(region_name, "_graph.rds"))
  saveRDS(g, out_graph)
  message("✅ Saved graph object: ", out_graph)
  
  # -------------------------------------------------------------------------
  # Graficar red bipartita
  # -------------------------------------------------------------------------
  p <- ggraph(g, layout = "bipartite") +
    geom_edge_link(aes(width = Jaccard, alpha = Jaccard), color = "gray40") +
    geom_node_point(aes(color = type, shape = type), size = 5) +
    geom_node_text(aes(label = name, color = type), repel = TRUE, size = 3) +
    scale_color_manual(values = c("AD" = "#D55E00", "Control" = "#0072B2")) +
    scale_edge_width_continuous(range = c(0.3, 2)) +
    scale_edge_alpha_continuous(range = c(0.4, 0.9)) +
    labs(
      title = paste("Bipartite network -", make_pretty_name(basename(f))),
      subtitle = paste("Edges: Jaccard >", threshold),
      edge_width = "Jaccard"
    ) +
    theme_void(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  out_png <- file.path(output_dir, paste0(region_name, "_bipartite.png"))
  ggsave(out_png, p, width = 10, height = 8, dpi = 300)
  message("✅ Saved plot: ", out_png)
}
