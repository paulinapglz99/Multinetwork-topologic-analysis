#!/usr/bin/env Rscript

#core_modules.R

#Packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
  data.table, 
  tidyverse,
  tools,
  optparse,
  ggplot2,
  igraph
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directorio con archivos *_nodes_summary.csv", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard",
              help = "Directorio de salida [default: %default]", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex para seleccionar archivos [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#Inputs
#opt$input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/"
#opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Functions

#A regex to get brain area and phenotype
extract_info <- function(filename) {
  base <- basename(filename)
  if (grepl("_counts_", base)) {
    region <- sub("_counts_.*", "", base)
    phenotype <- ifelse(grepl("_AD_", base), "AD", "Control")
    return(data.frame(filename = filename, region = region, phenotype = phenotype))
  } else {
    return(NULL)
  }
}

jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_matrix <- function(modA, modB) {
  mat <- matrix(0, nrow = length(modA), ncol = length(modB),
                dimnames = list(names(modA), names(modB)))
  for (i in seq_along(modA)) {
    for (j in seq_along(modB)) {
      mat[i, j] <- jaccard_index(modA[[i]], modB[[j]])
    }
  }
  return(mat)
}

load_modules <- function(file) {
  df <- fread(file)
  if (!all(c("node", "membership_infomap") %in% names(df))) {
    stop("File ", file, " must contain columns 'node' and 'membership_infomap'")
  }
  split(df$node, df$membership_infomap)
}

#Build pairs
build_pairs <- function(meta) {
  ad_files <- meta %>% filter(phenotype == "AD")
  control_files <- meta %>% filter(phenotype == "Control")
  
  pairs_ad_ad <- combn(ad_files$filename, 2, simplify = FALSE)
  pairs_control_control <- combn(control_files$filename, 2, simplify = FALSE)
  
  common_regions <- intersect(ad_files$region, control_files$region)
  pairs_ad_control <- lapply(common_regions, function(region) {
    ad_file <- ad_files %>% filter(region == !!region) %>% pull(filename)
    control_file <- control_files %>% filter(region == !!region) %>% pull(filename)
    list(ad_file, control_file)
  })
  
  return(list(
    AD_vs_AD = pairs_ad_ad,
    Control_vs_Control = pairs_control_control,
    AD_vs_Control = pairs_ad_control
  ))
}

#Compare networks
compare_and_export <- function(pairs, tag) {
  for (p in pairs) {
    fileA <- p[[1]]
    fileB <- p[[2]]
    nameA <- tools::file_path_sans_ext(basename(fileA))
    nameB <- tools::file_path_sans_ext(basename(fileB))
    name_out <- paste(tag, nameA, "vs", nameB, sep = "_")
    message("Comparando: ", name_out)
    
    netA <- load_modules(fileA)
    netB <- load_modules(fileB)
    
    jacmat <- jaccard_matrix(netA, netB)
    
    write.csv(jacmat, file = file.path(output_dir, paste0(name_out, "_jaccard_matrix.csv")))
  }
}

#Read files --- ---

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos dos archivos *_nodes_summary.csv.")

#Get metadata --- ---
meta <- do.call(rbind, lapply(files, extract_info))

#Get pairs!
all_pairs <- build_pairs(meta)
compare_and_export(all_pairs$AD_vs_AD, "ADvsAD")
compare_and_export(all_pairs$Control_vs_Control, "ControlvsControl")
compare_and_export(all_pairs$AD_vs_Control, "ADvsControl")

message("Comparisons DONE!. Files are now in: ", output_dir)

#Create summary of exclusive modules

generate_exclusive_module_report <- function(ad_files, control_files, region) {
  ad_file <- ad_files %>% filter(region == !!region) %>% pull(filename)
  control_file <- control_files %>% filter(region == !!region) %>% pull(filename)
  
  modAD <- load_modules(ad_file)
  modCTRL <- load_modules(control_file)
  
  jacmat <- jaccard_matrix(modAD, modCTRL)
  
  ad_exclusive <- rownames(jacmat)[apply(jacmat, 1, max) < 0.3]
  ctrl_exclusive <- colnames(jacmat)[apply(jacmat, 2, max) < 0.3]
  
  tibble(
    Region = region,
    AD_Modules_Exclusive = paste(ad_exclusive, collapse = ";"),
    Control_Modules_Exclusive = paste(ctrl_exclusive, collapse = ";")
  )
}

#Execute for each common region
summary_df <- bind_rows(
  lapply(intersect(meta$region[meta$phenotype == "AD"], meta$region[meta$phenotype == "Control"]), function(region) {
    generate_exclusive_module_report(
      ad_files = meta %>% filter(phenotype == "AD"),
      control_files = meta %>% filter(phenotype == "Control"),
      region = region
    )
  })
) 

#jeje
summary_df$Region <- gsub("^(Mayo_|ROSMAP_)", "", summary_df$Region)

write.csv(summary_df, file = file.path(output_dir, "summary_exclusive_modules.csv"), row.names = FALSE)
message("Core modules saved in 'summary_exclusive_modules.csv'")

#Visualisation

#Exclusive modules barplot per region --- ---
#How many exclusive modules has each region?

#Get long format
summary_long <- summary_df %>%
  mutate(
    AD_Count = sapply(strsplit(as.character(AD_Modules_Exclusive), ";"), function(x) sum(nzchar(x))),
    Control_Count = sapply(strsplit(as.character(Control_Modules_Exclusive), ";"), function(x) sum(nzchar(x)))
  ) %>%
  select(Region, AD_Count, Control_Count) %>%
  pivot_longer(cols = c("AD_Count", "Control_Count"), names_to = "Group", values_to = "Count") %>%
  mutate(Group = recode(Group, AD_Count = "AD", Control_Count = "Control"))

#Get cool colors
fill_colors <- c("AD" = rgb(0.3, 0.1, 0.4, 0.6),
  "Control" = rgb(0.3, 0.5, 0.4, 0.6))

#Plot
n_exclusive.p <- ggplot(summary_long, aes(x = Region, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0("n: ", Count)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),
    legend.position = "top"
  ) +
  labs(title = "Número de módulos exclusivos por región", y = NULL, x = NULL)

#Save barplot
ggsave(filename = file.path(output_dir, "barplot_exclusive_modules_ggplot.pdf"),
       plot = n_exclusive.p, width = 8, height = 5)

#Heatmap of phenotype exclusive modules --- ---

#Construcción de matriz binaria
all_modules <- list()
for (f in files) {
  name <- tools::file_path_sans_ext(basename(f))
  mods <- names(load_modules(f))
  all_modules[[name]] <- mods
}

# Unir todos los módulos
unique_modules <- unique(unlist(all_modules))
binary_matrix <- matrix(0, nrow = length(unique_modules), ncol = length(all_modules),
                        dimnames = list(unique_modules, names(all_modules)))

# Rellenar con 1 si el módulo está en esa red
for (net in names(all_modules)) {
  binary_matrix[all_modules[[net]], net] <- 1
}

# Transponer matriz
#binary_matrix_t <- t(binary_matrix)

# -----------------------------
# 🔍 Filtrar solo módulos exclusivos (AD o Control)

# Extraer listas desde summary_df
ad_exclusive_all <- unique(unlist(strsplit(summary_df$AD_Modules_Exclusive, ";")))
control_exclusive_all <- unique(unlist(strsplit(summary_df$Control_Modules_Exclusive, ";")))
exclusive_modules <- union(ad_exclusive_all, control_exclusive_all)

# Filtrar columnas de la matriz (módulos)
binary_matrix_exclusive <- binary_matrix[, colnames(binary_matrix) %in% exclusive_modules]

# -----------------------------
# 🏷 Anotar columnas por tipo de exclusividad
module_status <- ifelse(colnames(binary_matrix_exclusive) %in% ad_exclusive_all, "AD_exclusive",
                        ifelse(colnames(binary_matrix_exclusive) %in% control_exclusive_all, "Control_exclusive", "Shared"))

annotation_col_exclusive <- data.frame(Exclusivity = factor(module_status, levels = c("AD_exclusive", "Control_exclusive")))
rownames(annotation_col_exclusive) <- colnames(binary_matrix_exclusive)

# Colores
ann_colors <- list(
  Exclusivity = c(
    AD_exclusive = rgb(0.3, 0.1, 0.4, 0.8),
    Control_exclusive = rgb(0.3, 0.5, 0.4, 0.8)
  )
)

# -----------------------------
# 🔥 Heatmap solo de módulos exclusivos

pheatmap(binary_matrix_exclusive,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = c("white", "black"),
         annotation_col = annotation_col_exclusive,
         annotation_colors = ann_colors,
         main = "Módulos Exclusivos por Red",
         filename = file.path(output_dir, "heatmap_modulos_exclusivos_solo.pdf"),
         fontsize_col = 6,
         fontsize_row = 7)

message("✅ Heatmap de módulos exclusivos generado.")

# -----------------------------
# 💾 Exportar matriz binaria exclusiva como CSV

write.csv(binary_matrix_exclusive_t,
          file = file.path(output_dir, "binary_matrix_modulos_exclusivos.csv"),
          quote = FALSE)







# -----------------------------
# 🔗 Red de módulos basada en Jaccard

pacman::p_load(igraph)

message("🔗 Generando red de módulos basada en similitud Jaccard...")

# Recoger pares con Jaccard >= 0.5 de todas las comparaciones
edge_list <- list()
threshold <- 0.5

jaccard_files <- list.files(output_dir, pattern = "_jaccard_matrix.csv$", full.names = TRUE)

for (f in jaccard_files) {
  mat <- as.matrix(read.csv(f, row.names = 1))
  if (nrow(mat) > 0 && ncol(mat) > 0) {
    for (i in rownames(mat)) {
      for (j in colnames(mat)) {
        val <- mat[i, j]
        if (!is.na(val) && val >= threshold) {
          edge_list[[length(edge_list) + 1]] <- data.frame(from = paste0(i, "_A"), to = paste0(j, "_B"), weight = val)
        }
      }
    }
  }
}

edges_df <- do.call(rbind, edge_list)
if (nrow(edges_df) > 0) {
  g <- graph_from_data_frame(edges_df, directed = FALSE)
  pdf(file.path(output_dir, "modular_network_graph.pdf"), width = 10, height = 10)
  plot(g, vertex.label = NA, vertex.size = 5,
       edge.width = E(g)$weight * 3, edge.color = "gray50",
       main = "Red de Módulos con Jaccard ≥ 0.5")
  dev.off()
  message("✅ Red de módulos guardada como 'modular_network_graph.pdf'")
} else {
  message("⚠️ No se encontraron pares con Jaccard ≥ ", threshold, ". No se generó red.")
}



