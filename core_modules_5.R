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
opt$input_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_topos/"
opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules_5"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Functions --- ---

#Extract metadata
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

#Jaccard
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

#Comparison between phenotypes by region

#Extract unique names from regions
regions <- unique(meta$region)

#Thresholds
upp_thres <- 0.8
low_thres <- 0.25

#Loop to calculate Jaccards per region
jaccards <- list()
for (current_region in regions) {
  #Obtain module names for this region
  mods_control <- modules[grep(paste0(current_region, "_Control"), names(modules))]
  mods_ad <- modules[grep(paste0(current_region, "_AD"), names(modules))]
  
  #Compare all modules Control vs AD
  for (ctrl_name in names(mods_control)) {
    for (ad_name in names(mods_ad)) {
      jaccard_val <- jaccard_index(mods_control[[ctrl_name]], mods_ad[[ad_name]])
      
      clasificacion <- ifelse(
        jaccard_val >= upp_thres, "Similar",
        ifelse(jaccard_val <= low_thres, "Dissimilar", "Intermediate")
      )
      
      jaccards[[length(jaccards) + 1]] <- data.frame(
        Region = current_region,
        Module_Control = ctrl_name,
        Module_AD = ad_name,
        Jaccard_Index = jaccard_val,
        Classification = clasificacion,
        stringsAsFactors = FALSE
      )
    }
  }
}

#Bind results
jaccards.tb <- rbindlist(jaccards)
#jeje
jaccards.tb$Region <- gsub("^(Mayo_|ROSMAP_)", "", jaccards.tb$Region)
#Factors in order
ordered_regions <- c("HCN", "PCC", "TC", "CRB","DLPFC")
jaccards.tb$Region <- factor(jaccards.tb$Region, levels = ordered_regions)

print(table(jaccards.tb$Classification))

#Save results
fwrite(jaccards.tb, file = file.path(output_dir, "jaccard_summary_by_region.csv"))

#Histogram of Jaccard values
jaccard_hist.p <- ggplot(jaccards.tb, aes(x = Jaccard_Index)) +
  geom_histogram(binwidth = 0.005, fill = "navyblue", color = "black") +
  geom_vline(xintercept = c(low_thres, upp_thres), linetype = "dashed", color = c("red", "darkgreen"), size = 1) +
  scale_y_log10() +
  labs(
    title = "Distribution of Jaccard Index Between Modules (AD vs Control)",
    x = "Jaccard Index",
    y = "Frequency (log)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
jaccard_hist.p

#Save histogram
ggsave(filename = file.path(output_dir, "jaccard_histogram.pdf"), 
       plot = jaccard_hist.p, 
       width = 10, 
       height = 7)






#Filter exclusive modules
candidates_exclusive_modules <- jaccards.tb %>%
  filter(Classification %in% c("Dissimilar", "Similar"))

#Count by region and phenotype
candidates_exclusive_modules.c <- candidates_exclusive_modules %>%
  group_by(Region, Classification) %>%
  summarise(n = n(), .groups = "drop")

#Get cool colors
fill_colors <- c("Dissimilar" = rgb(0.3, 0.1, 0.4, 0.6),
                 "Similar" = rgb(0.3, 0.5, 0.4, 0.6))

#Barplot
barplot.p <- ggplot(candidates_exclusive_modules.c, aes(x = Region, y = n, fill = Classification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0("n: ", n)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Exclusive Modules by Region and Phenotype",
       x = "Brain Region",
       y = "Number of Exclusive Modules",
       fill = "Classification") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 0, size = 11, margin = margin(t = -5)),
    axis.text.y = element_text(size = 11),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

barplot.p

#Save barplot
ggsave(filename = file.path(output_dir, "similar_dissimilar_barplot.pdf"), 
       plot = barplot.p, 
       width = 7, height = 10)

message("First Jaccard comparison completed and table saved!")


#A module is considered unique to a phenotype (e.g., AD) if it has no similar module in the opposite phenotype in the same region.
#This means:
#If an X_AD module does not have any Y_Control modules with Jaccard_Index ≥ upp_thres → it is AD-exclusive.
#If a Y_Control module does not have any X_AD modules with Jaccard_Index ≥ upp_thres → it is Control-exclusive.

#Check if there is at least one row with Classification == ‘Similar’ → then it is NOT exclusive.
#If no row has ‘Similar’, then this Module_AD is AD-exclusive.

#Extraer módulos AD y Control por región

exclusive_AD <- jaccards.tb %>%
  group_by(Region, Module_AD) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    .groups = "drop"
  ) %>%
  filter(!has_similar) %>%
  select(Region, Module_AD) %>%
  mutate(Phenotype = "AD_exclusive")

exclusive_CTRL <- jaccards.tb %>%
  group_by(Region, Module_Control) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    .groups = "drop"
  ) %>%
  filter(!has_similar) %>%
  select(Region, Module_Control) %>%
  mutate(Phenotype = "Control_exclusive")

# 2️⃣ Unir resultados en una sola tabla
exclusive_modules <- bind_rows(
  exclusive_AD %>% rename(Module = Module_AD),
  exclusive_CTRL %>% rename(Module = Module_Control)
)

# 3️⃣ Guardar resultados
fwrite(exclusive_modules, file = file.path(output_dir, "exclusive_modules_AD_CTRL.csv"))

# 4️⃣ Resumen por región
exclusive_summary <- exclusive_modules %>%
  group_by(Region, Phenotype) %>%
  summarise(n_modules = n(), .groups = "drop")

fwrite(exclusive_summary, file = file.path(output_dir, "exclusive_modules_summary.csv"))

# 5️⃣ (Opcional) Mostrar resultados en consola
message("Resumen de módulos exclusivos por región:")
print(exclusive_summary)

# 6️⃣ (Opcional) Visualización rápida
ggplot(exclusive_summary, aes(x = Region, y = n_modules, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n_modules), vjust = -0.3, size = 3.5, position = position_dodge(width = 0.8)) +
  labs(
    title = "Número de módulos exclusivos por región y fenotipo",
    x = "Región cerebral",
    y = "Número de módulos exclusivos",
    fill = "Fenotipo"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top")

#"desacuerdo fenotípico"