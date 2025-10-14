# ================================================
# compare_networks_modules_Rstudio.R
# Calcula matrices de similitud de Jaccard entre
# módulos definidos en archivos *_nodes_summary.csv
# ================================================

# -------------------------------
# 1. Cargar paquetes
# -------------------------------
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

ok <- pacman::p_load(
  data.table, tidyverse, tools
)

if (all(ok)) {
  message("✅ All packages loaded correctly.")
} else {
  stop("❌ Some packages did not load correctly: ",
       paste(names(ok)[!ok], collapse = ", "))
}

# -------------------------------
# 2. CONFIGURACIÓN MANUAL
# -------------------------------
# 📂 Directorio donde están tus archivos *_nodes_summary.csv
input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/"

# 📁 Directorio de salida para guardar los resultados
out_dir <- "results_jaccard"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 🔍 Patrón de búsqueda de archivos
pattern <- "_nodes_summary\\.csv$"

# 🧩 PARES de redes que quieres comparar
# (usa los nombres de archivo SIN extensión)
# Ejemplo: pares <- list(c("network_3_nodes_summary", "network_7_nodes_summary"))
pares <- list(
  c("network_3_nodes_summary", "network_7_nodes_summary"),
  c("network_4_nodes_summary", "network_8_nodes_summary")
)

# 🔢 Semilla (para reproducibilidad)
set.seed(42)

# -------------------------------
# 3. Funciones auxiliares
# -------------------------------

# Índice de Jaccard
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Matriz de Jaccard entre dos redes
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

# -------------------------------
# 4. Leer módulos desde archivos
# -------------------------------
files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2) stop("⚠️ Se necesitan al menos 2 archivos *_nodes_summary.csv para comparar.")

message("📂 Found ", length(files), " module summary files.")

# Leer los módulos
redes_modulos <- list()

for (f in files) {
  net_name <- tools::file_path_sans_ext(basename(f))
  df <- data.table::fread(f)
  
  if (!all(c("node", "membership_infomap") %in% colnames(df))) {
    stop("❌ File ", f, " must contain columns 'node' and 'membership_infomap'")
  }
  
  # Agrupar nodos por módulo
  mod_list <- split(df$node, df$membership_infomap)
  redes_modulos[[net_name]] <- mod_list
}

# -------------------------------
# 5. Detectar pares AD vs Control
# -------------------------------

# Detectar tipo de red (prefijo hasta "_counts_")
names_df <- tibble(name = names(redes_modulos)) %>%
  mutate(
    prefix = str_extract(name, "^[A-Za-z0-9_]+_counts"),
    condition = case_when(
      str_detect(name, "_AD_") ~ "AD",
      str_detect(name, "_control_") ~ "control",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(condition)

# Crear pares AD vs control por cada tipo de red
pares_auto <- names_df %>%
  group_by(prefix) %>%
  filter(all(c("AD", "control") %in% condition)) %>%
  summarise(
    ad = name[condition == "AD"][1],
    ctrl = name[condition == "control"][1],
    .groups = "drop"
  ) %>%
  mutate(pair = map2(ad, ctrl, ~c(.x, .y))) %>%
  pull(pair)

if (length(pares_auto) == 0)
  stop("⚠️ No se detectaron pares AD vs control automáticamente.")

message("🧩 Pares detectados automáticamente:")
print(pares_auto)

# -------------------------------
# 6. Comparar pares
# -------------------------------
comparar_redes <- function(lista_redes, pares) {
  resultados <- list()
  for (p in pares) {
    redA_name <- p[1]
    redB_name <- p[2]
    redA <- lista_redes[[redA_name]]
    redB <- lista_redes[[redB_name]]
    nombre <- paste(redA_name, redB_name, sep = "_vs_")
    
    message("➡️  Comparing ", nombre)
    resultados[[nombre]] <- jaccard_matrix(redA, redB)
  }
  return(resultados)
}

resultados <- comparar_redes(redes_modulos, pares_auto)

# -------------------------------
# 7. Guardar resultados
# -------------------------------
for (n in names(resultados)) {
  out_file <- file.path(out_dir, paste0("Jaccard_", n, ".csv"))
  write.csv(resultados[[n]], out_file, row.names = TRUE)
  message("💾 Saved: ", out_file)
}

message("✅ Finished computing all Jaccard matrices.")


message("✅ Finished computing all Jaccard matrices.")

library(ggplot2)
# -------------------------------
# 8. Visualizar y guardar heatmaps
# -------------------------------
message("🎨 Generando heatmaps de Jaccard...")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (n in names(resultados)) {
  mat <- resultados[[n]]
  
  df_plot <- as.data.frame(as.table(mat))
  colnames(df_plot) <- c("Modulo_A", "Modulo_B", "Jaccard")
  
  p <- ggplot(df_plot, aes(Modulo_B, Modulo_A, fill = Jaccard)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal(base_size = 12) +
    labs(
      title = str_replace(n, "_top10pct_nodes_summary", ""),
      subtitle = "Índice de similitud de Jaccard entre módulos",
      fill = "Jaccard"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 13)
    )
  
  out_png <- file.path(plot_dir, paste0("Heatmap_", n, ".png"))
  ggsave(out_png, plot = p, width = 8, height = 6, dpi = 300)
  message("🖼️  Saved plot: ", out_png)
}

message("✅ Todos los heatmaps generados y guardados en ", plot_dir)
