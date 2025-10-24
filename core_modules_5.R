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
              help = "Regex para seleccionar archivos [default: %default]", metavar = "character"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#Inputs
#opt$input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/"
#opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules_5"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Thresholds
upp_thres <- 0.8
low_thres <- 0.5

# Functions --- #

# Extract metadata
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

#Load modules with unique names
load_modules <- function(file) {
  df <- fread(file)
  if (!all(c("node", "membership_infomap") %in% names(df))) {
    stop("File ", file, " must contain columns 'node' and 'membership_infomap'")
  }
  region <- sub("_counts_.*", "", basename(file))
  phenotype <- ifelse(grepl("_AD_", file), "AD", "Control")
  df$unique_module_id <- paste(region, phenotype, df$membership_infomap, sep = "_")
  split(df$node, df$unique_module_id)
}

# Jaccard
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#Read and Prepare Files
files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
meta <- do.call(rbind, lapply(files, extract_info))
module_list <- lapply(files, load_modules)
modules <- unlist(module_list, recursive = FALSE)

#FIRST: COMPARE NETWORKS IN EACH REGION PER PHENOTYPE TO FIND AD-exclusive and ctrl-exclusive modules in each region
results <- list()

regions <- unique(meta$region)

for (region in regions) {
  # Extraer módulos AD y Control por región
  ad_mods <- modules[grepl(paste0("^", region, "_AD_"), names(modules))]
  ctrl_mods <- modules[grepl(paste0("^", region, "_Control_"), names(modules))]
  
  for (ad_name in names(ad_mods)) {
    for (ctrl_name in names(ctrl_mods)) {
      jacc <- jaccard_index(ad_mods[[ad_name]], ctrl_mods[[ctrl_name]])
      
      # Para AD
      results[[length(results) + 1]] <- data.frame(
        module_name = ad_name,
        region = region,
        module_pheno = "AD",
        comparison_module = ctrl_name,
        jaccard = jacc,
        classification = ifelse(jacc >= upp_thres, "similar", ifelse(jacc <= low_thres, "AD-exclusive", "intermediate"))
      )
      
      # Para Control
      results[[length(results) + 1]] <- data.frame(
        module_name = ctrl_name,
        region = region,
        module_pheno = "Control",
        comparison_module = ad_name,
        jaccard = jacc,
        classification = ifelse(jacc >= upp_thres, "similar", ifelse(jacc <= low_thres, "Ctrl-exclusive", "intermediate"))
      )
    }
  }
}

# Tabla resumen
summary_df <- rbindlist(results)
summary_df <- summary_df[!duplicated(paste(module_name, module_pheno)), ]
#jeje
summary_df$region <- gsub("^(Mayo_|ROSMAP_)", "", summary_df$region)

#Factors in order
orden_regiones <- c("HCN", "PCC", "TC", "CRB","DLPFC")
summary_df$Region <- factor(summary_df$Region, levels = orden_regiones)

# Guardar resultados
fwrite(summary_df, file.path(output_dir, "modules_jaccard_summary.tsv"), sep = "\t")

# Continuar con Task_3

# Filtrar módulos exclusivos
exclusive_modules <- summary_df %>%
  filter(classification %in% c("AD-exclusive", "Ctrl-exclusive"))

# Contar por región y fenotipo
exclusive_counts <- exclusive_modules %>%
  group_by(region, classification) %>%
  summarise(n = n(), .groups = "drop")

#Get cool colors
fill_colors <- c("AD-exclusive" = rgb(0.3, 0.1, 0.4, 0.6),
                 "Ctrl-exclusive" = rgb(0.3, 0.5, 0.4, 0.6))

# Barplot
p <- ggplot(exclusive_counts, aes(x = region, y = n, fill = classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Exclusive Modules by Region and Phenotype",
       x = "Brain Region",
       y = "Number of Exclusive Modules",
       fill = "Classification") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(hjust = 1),
    legend.position = c(0.80, 0.98),  # Posición dentro del gráfico (x, y)
    legend.justification = c(0, 1),   # Justificar desde la esquina superior izquierda
    #legend.background = element_rect(fill = "transparent")  # Fondo transparente
  ) +
  #guides(fill = guide_legend(override.aes = list(shape = 21))) +  # Círculos en leyenda
  labs(title = "Exclusive modules per region", y = NULL, x = NULL)

p

# Guardar plot
ggsave(filename = file.path(output_dir, "exclusive_modules_barplot.png"), plot = p, width = 8, height = 6)

# También guardar los datos del plot por si se quieren revisar
fwrite(exclusive_counts, file.path(output_dir, "exclusive_modules_counts.tsv"), sep = "\t")
