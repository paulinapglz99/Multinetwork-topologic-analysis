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
opt$input_dir <- "~/Desktop/local_work/results_topos/"
opt$output_dir <- "~/Desktop/local_work/results_topos/results_core_modules"

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
#names(module_list) <- basename(files)

modules <- unlist(module_list, recursive = FALSE)
#FIRST: COMPARE NETWORKS IN EACH REGION PER PHENOTYPE TO FIND AD-exclusive and ctrl-exclusive modules in each region

#Comparison between phenotypes by region

#Extract unique names from regions
#jeje
meta$region <- gsub("^(Mayo_|ROSMAP_)", "", meta$region)
regions <- unique(meta$region)

#Count modules per network
module_counts <- tibble(
  Region = meta$region,
  Phenotype = meta$phenotype,
  N_modules = sapply(module_list, length),
  Network = paste0(Region, "_", Phenotype)
) %>%
  group_by(Region) %>%
  mutate(Total_region_modules = sum(N_modules)) %>%
  ungroup()

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
#Factors in order
ordered_regions <- c("HCN", "PCC", "TC", "CRB","DLPFC")
jaccards.tb$Region <- factor(jaccards.tb$Region, levels = ordered_regions)
print(table(jaccards.tb$Classification))

#Save results
fwrite(jaccards.tb, file = file.path(output_dir, "jaccards_all_regions.csv"))

#Histogram of Jaccard values

jaccard_hist.p <- ggplot(jaccards.tb, aes(x = Jaccard_Index)) +
  geom_histogram(binwidth = 0.005, fill = "navyblue", color = "black", alpha = 0.7) +
  geom_density(aes(y = ..count.. * 0.005), color = "orange", size = 1) +
  geom_vline(xintercept = low_thres, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = upp_thres, linetype = "dashed", color = "darkgreen", size = 1) +
  scale_y_continuous(trans = "log10", limits = c(1, NA)) +  # evita problemas con ceros
  labs(
    title = "Global",
    x = "Jaccard Index",
    y = "Frequency (log10)"
  ) +
  # facet_wrap(~Region) +
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
# 
# jaccard_region.p <- ggplot(jaccards.tb, aes(x = Jaccard_Index)) +
#   geom_histogram(binwidth = 0.005, fill = "navyblue", color = "black", alpha = 0.7) +
#   geom_density(aes(y = ..count.. * 0.005), color = "orange", size = 1) +
#   geom_vline(xintercept = low_thres, linetype = "dashed", color = "red", size = 1) +
#   geom_vline(xintercept = upp_thres, linetype = "dashed", color = "darkgreen", size = 1) +
#   scale_y_continuous(trans = "log10", limits = c(1, NA)) +  # evita problemas con ceros
#   labs(
#     title = "Distribution of Jaccard Index Between Modules (AD vs Control)",
#     x = "Jaccard Index",
#     y = "Frequency (log10)"
#   ) +
#   facet_wrap(~Region) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#     axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
#     axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
#     axis.text = element_text(size = 11),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# jaccard_region.p

library(cowplot)

# Lista de gráficos por región
jaccard_region.p <- lapply(levels(jaccards.tb$Region), function(reg) {
  ggplot(jaccards.tb %>% filter(Region == reg), aes(x = Jaccard_Index)) +
    geom_histogram(binwidth = 0.005, fill = "navyblue", color = "black", alpha = 0.7) +
    geom_density(aes(y = ..count.. * 0.005), color = "orange", size = 1) +
    geom_vline(xintercept = low_thres, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = upp_thres, linetype = "dashed", color = "darkgreen", size = 1) +
    scale_y_continuous(trans = "log10", limits = c(1, NA)) +
    labs(
      title = reg,
      x = "Jaccard Index",
      y = "Frequency (log10)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
})

# Panel de 2 filas x 3 columnas (ajustar ncol/nrow según número de regiones)
jaccard_region.p <- cowplot::plot_grid(plotlist = plots_por_region, ncol = 3)
jaccard_region.p

jaccard.p <- plot_grid(jaccard_hist.p, 
                       plots_por_region[[1]], 
                       plots_por_region[[2]],
                       plots_por_region[[3]],
                       plots_por_region[[4]],
                       plots_por_region[[5]],
                       ncol = 3)
jaccard.p

#Save histogram
ggsave(filename = file.path(output_dir, "jaccard_histogram.pdf"), 
       plot = jaccard_hist.p, 
       width = 10, 
       height = 7)

#all --- ---

#A module is considered unique to a phenotype (e.g., AD) if it has no similar module in the opposite phenotype in the same region.
#This means:
#If an X_AD module does not have any Y_Control modules with Jaccard_Index ≥ upp_thres → it is AD-exclusive.
#If a Y_Control module does not have any X_AD modules with Jaccard_Index ≥ upp_thres → it is Control-exclusive.

#Check if there is at least one row with Classification == ‘Similar’ → then it is NOT exclusive.
#If no row has ‘Similar’, then this Module_AD is AD-exclusive.

# Clasificación extendida para módulos AD
classified_AD <- jaccards.tb %>%
  group_by(Region, Module_AD) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    has_intermediate = any(Classification == "Intermediate"),
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      has_similar ~ "Similar",
      has_intermediate ~ "Intermediate",
      TRUE ~ "AD_exclusive"
    )
  ) %>%
  rename(Module = Module_AD) %>%
  mutate(Phenotype = "AD")

# Clasificación extendida para módulos Control
classified_CTRL <- jaccards.tb %>%
  group_by(Region, Module_Control) %>%
  summarise(
    has_similar = any(Classification == "Similar"),
    has_intermediate = any(Classification == "Intermediate"),
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      has_similar ~ "Similar",
      has_intermediate ~ "Intermediate",
      TRUE ~ "Control_exclusive"
    )
  ) %>%
  rename(Module = Module_Control) %>%
  mutate(Phenotype = "Control")

# Unir ambas tablas
classified_all <- bind_rows(classified_AD, classified_CTRL)

# Contar módulos por clasificación y región
summary_all <- classified_all %>%
  group_by(Region, Phenotype, Classification) %>%
  summarise(n_modules = n(), .groups = "drop") %>%
  mutate(Network = paste0(Region, "_", Phenotype)) %>% 
  select(Network, n_modules, Classification) %>%
  left_join(module_counts, by = "Network") %>%
  mutate(local_proportion = n_modules / N_modules * 100) %>% 
  mutate(global_proportion = n_modules / Total_region_modules * 100) %>%
  mutate(across(c(local_proportion, global_proportion), round, 1))

#Save result
fwrite(summary_all, file = file.path(output_dir, "modules_classification_summary.csv"))

#Visualize local proportion

#Get cool colors
fill_colors <- c(
  "Control_exclusive" = rgb(0.1, 0.4, 0.3, 0.6),   # verde oscuro
  "Intermediate"      = rgb(0.4, 0.6, 0.3, 0.6),   # verde oliva
  "Similar"           = rgb(0.4, 0.3, 0.6, 0.6),   # violeta suave
  "AD_exclusive"      = rgb(0.6, 0.2, 0.5, 0.6)    # magenta profundo
)

local.p <- ggplot(summary_all, aes(x = Region, y = local_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Phenotype) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "Local proportion of modules by classification and phenotype",
    y = "%",
    x = "Region"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
local.p

#Visualize global proportion

# Graficar proporción por región
global.p <- ggplot(summary_all, aes(x = Region, y = global_proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "Global proportion of modules by classification and phenotype",
    x = "Region",
    y = "%"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
global.p 

jaccard_heatmap.p <- ggplot(summary_all, aes(x = Region, y = Classification, fill = global_proportion)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~Phenotype) +
  theme_minimal() +
  labs(fill = "% global", title = "Proporción global por clasificación y región")
jaccard_heatmap.p

#Grid panel
library(cowplot)

plot_grid(local.p, global.p)

#Count by region and phenotype
candidates_exclusive_modules.c <- candidates_exclusive_modules %>%
  group_by(Region, Classification) %>%
  summarise(n = n(), .groups = "drop")

#Filter exclusive modules --- ---
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

#Extract AD and Control modules by region

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

#Bind results
exclusive_modules <- bind_rows(
  exclusive_AD %>% rename(Module = Module_AD),
  exclusive_CTRL %>% rename(Module = Module_Control)
)

module_gene_counts <- tibble(
  Module = names(unlist(module_list, recursive = FALSE)),
  N_genes = sapply(unlist(module_list, recursive = FALSE), length)
)

# Agregar el número de genes a exclusive_modules
exclusive_modules <- exclusive_modules %>%
  left_join(module_gene_counts, by = "Module") %>%
  arrange(Region, Phenotype, desc(N_genes))

#Save results
fwrite(exclusive_modules, file = file.path(output_dir, "exclusive_modules_AD_CTRL.csv"))

#DLPFC muestra más módulos "compartidos" entre AD y Control → podría reflejar una transición gradual o plasticidad en esta región.
#CRB y TC muestran módulos altamente específicos → pueden tener procesos moleculares altamente diferenciados entre fenotipos.

#check

exclusive_summary <- summary_all %>% filter(Classification %in% c("AD_exclusive",
                                                                  "Control_exclusive"))

#Factors in order
exclusive_summary$Region <- factor(exclusive_summary$Region, levels = ordered_regions)

#Get cool colors
fill_colors <- c("AD_exclusive" = rgb(0.3, 0.1, 0.4, 0.6),
                 "Control_exclusive" = rgb(0.3, 0.5, 0.4, 0.6))

exclusive.p <- ggplot(exclusive_summary, aes(x = Region, y = n_modules, fill = Classification))  +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0("n: ", n_modules)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Number of exclusive Modules by Region and Phenotype",
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
    legend.position = c(0.7, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank())
exclusive.p 

# exclusive_v2.p <- ggplot(exclusive_summary, aes(x = Region, fill = Classification)) +
#   geom_bar(aes(y = N_modules),
#            stat = "identity",
#            position = position_dodge(),
#            alpha = 0.3,
#            color = NA) +
#   geom_bar(aes(y = n_modules),
#            stat = "identity",
#            position = position_dodge(),
#            #color = "black",
#            size = 0.2) +
#   geom_text(aes(y = n_modules,
#                 label = sprintf("%.1f%%", Proportion)),
#             position = position_dodge(width = 0.9),
#             vjust = -0.8,
#             size = 4,
#             family = "sans") +
#   scale_fill_manual(values = fill_colors) +
#   labs(
#     title = "Exclusive Modules Across Brain Regions",
#     x = "Brain Region",
#     y = "Exclusive Module Count",
#     fill = "Phenotype"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#     axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
#     axis.text = element_text(color = "black"),
#     axis.title = element_text(face = "bold"),
#     legend.position = "top",
#     legend.title = element_text(face = "bold"),
#     legend.text = element_text(size = 12),
#     panel.grid = element_blank(),
#     #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
#   )
# exclusive_v2.p

# EL "desacuerdo fenotípico"
#¿Existen módulos AD-exclusivos que estén presentes (conservados) en varias regiones?
# ¿Hay módulos Control-exclusivos que también aparezcan en múltiples regiones?

#Tomas los módulos AD-exclusivos de cada región, por ejemplo:
#
# Mayo_CRB_AD_, Mayo_PCC_AD_3, etc.
# 
# Y luego comparas cada par entre regiones usando, por ejemplo, índice de Jaccard entre sus genes.
# 
# Si dos módulos de distintas regiones tienen un índice de Jaccard ≥ umbral (ej. 0.8), los consideras “conservados entre regiones”.

#Meta: Identificar módulos exclusivos (AD o Control) que aparecen en más de una región con alta similitud (Jaccard_Index ≥ 0.8).

# Jaccard entre módulos exclusivos de diferentes regiones (dentro de cada fenotipo)

find_conserved_modules <- function(exclusive_df, modules_list, threshold = 0.8) {
  conserved_list <- list()
  
  # Filtrar por fenotipo
  phenos <- unique(exclusive_df$Phenotype)
  
  for (pheno in phenos) {
    sub_df <- exclusive_df %>% filter(Phenotype == pheno)
    
    # Agrupar por módulos de diferentes regiones
    for (i in 1:(nrow(sub_df) - 1)) {
      for (j in (i + 1):nrow(sub_df)) {
        mod1 <- sub_df$Module[i]
        mod2 <- sub_df$Module[j]
        
        reg1 <- sub_df$Region[i]
        reg2 <- sub_df$Region[j]
        
        # Saltar si son de la misma región
        if (reg1 == reg2) next
        
        genes1 <- modules[[mod1]]
        genes2 <- modules[[mod2]]
        
        jaccard_val <- jaccard_index(genes1, genes2)
        
        if (jaccard_val >= threshold) {
          conserved_list[[length(conserved_list) + 1]] <- data.frame(
            Module_1 = mod1,
            Region_1 = reg1,
            Module_2 = mod2,
            Region_2 = reg2,
            Phenotype = pheno,
            Jaccard_Index = jaccard_val,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  bind_rows(conserved_list)
}

# Ejecutar
conserved_modules_across_regions <- find_conserved_modules(
  exclusive_modules,
  modules,
  threshold = upp_thres
)

conserved_modules_across_regions %>% 

# Guardar
fwrite(conserved_modules_across_regions,
       file = file.path(output_dir, "conserved_exclusive_modules_across_regions.csv"))

#Generar un grafo de módulos conservados entre regiones

library(igraph)

# Crear objeto grafo
g <- graph_from_data_frame(
  conserved_modules_across_regions %>%
    select(from = Module_1, to = Module_2, weight = Jaccard_Index),
  directed = FALSE
)

# Agregar atributos: región y fenotipo
# Extraer nodos únicos
all_nodes <- unique(c(conserved_modules_across_regions$Module_1,
                      conserved_modules_across_regions$Module_2))

# Tabla con info de cada nodo
node_attributes <- exclusive_modules %>%
  filter(Module %in% all_nodes) %>%
  distinct(Module, Region, Phenotype)

# Asegurar que todos los nodos estén
V(g)$name <- as.character(V(g)$name)
V(g)$Region <- node_attributes$Region[match(V(g)$name, node_attributes$Module)]
V(g)$Phenotype <- node_attributes$Phenotype[match(V(g)$name, node_attributes$Module)]

# Colores por fenotipo
phenotype_colors <- c("AD_exclusive" = "#7B3294", "Control_exclusive" = "#008837")
V(g)$color <- phenotype_colors[V(g)$Phenotype]

set.seed(42)  # reproducible layout

plot(
  g,
  vertex.label = NA,       # o usar V(g)$name para etiquetas
  vertex.size = 6,
  vertex.color = V(g)$color,
  edge.width = E(g)$weight * 5,  # grosor según Jaccard
  edge.color = "gray60",
  layout = layout_with_fr(g),  # o layout_with_kk
  main = "Grafo de módulos exclusivos conservados entre regiones"
)

