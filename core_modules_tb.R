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
  igraph, 
  cowplot, 
  pheatmap
)

#Parser
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = getwd(),
              help = "Directorio con archivos *_nodes_summary.csv", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./results_jaccard",
              help = "Directorio de salida [default: %default]", metavar = "character"),
  make_option(c("-p", "--pattern"), type = "character", default = "_nodes_summary\\.csv$",
              help = "Regex para seleccionar archivos [default: %default]", metavar = "character"), 
  make_option(c("-r", "--enrich_dir"), type = "character", default = "./results_comm",
              help = "Directorio de enriquecimientos [default: %default]", metavar = "character"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#Inputs
opt$input_dir <- "~/Desktop/local_work/results_topos/"
opt$output_dir <- "~/Desktop/local_work/results_topos/results_core_modules"
opt$enrich_dir <- "~/Desktop/local_work/results_comm/"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern
enrich_dir <- opt$enrich_dir

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

#Read enrichments
enrich_files <- list.files(enrich_dir, pattern = "_enrichment\\.csv$", full.names = TRUE)
enrich_data <- lapply(enrich_files, fread)

#Add metadata
enrich_data <- Map(function(df, filename) {
  fname <- basename(filename)
  region <- sub("^(Mayo_|ROSMAP_)", "", sub("_counts_.*", "", fname))
  phenotype <- ifelse(grepl("_AD_", fname), "AD", "Control")
  df$Region <- region
  df$Phenotype <- phenotype
  return(df)
}, enrich_data, enrich_files)

#Bind enrichment tables
enrich_all <- rbindlist(enrich_data)

enrich_all <- enrich_all %>%
  mutate(
    CommunityID = as.character(CommunityID),
    Module = paste(Region, Phenotype, CommunityID, sep = "_")
  )

modules <- unlist(module_list, recursive = FALSE)
# Filter modules with at least N genes
min_genes <- 5
modules <- modules[sapply(modules, length) >= min_genes]
#Recalculate `module_list` by filtering for minimum size
module_list <- lapply(module_list, function(mods) {
  mods[sapply(mods, length) >= min_genes]
})

#FIRST: COMPARE NETWORKS IN EACH REGION PER PHENOTYPE TO FIND AD-exclusive and ctrl-exclusive modules in each region
#Comparison between phenotypes per region

#Extract unique names from regions
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

module_gene_counts <- tibble(
  Module = names(unlist(module_list, recursive = FALSE)),
  N_genes = sapply(unlist(module_list, recursive = FALSE), length)
)

#A module is considered unique to a phenotype (e.g., AD) if it has no similar module in the opposite phenotype in the same region.
#This means:
#If an X_AD module does not have any Y_Control modules with Jaccard_Index ≥ upp_thres → it is AD-exclusive.
#If a Y_Control module does not have any X_AD modules with Jaccard_Index ≥ upp_thres → it is Control-exclusive.

#Check if there is at least one row with Classification == ‘Similar’ → then it is NOT exclusive.
#If no row has ‘Similar’, then this Module_AD is AD-exclusive.

#Thresholds
upp_thres <- 0.7
low_thres <- 0.2

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

#Save results <- RESULT #1
vroom::vroom_write(jaccards.tb, file = file.path(output_dir, "jaccards_all_regions.csv"))
message("First Jaccard comparison completed and table saved!")

#Reorder modules to have a summary table 
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

#Bind tables
classified_all <- bind_rows(classified_AD, classified_CTRL)

#Make a summary per classification and region
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
vroom::vroom_write(summary_all, file = file.path(output_dir, "modules_classification_summary.csv"))

#Filter exclusive modules --- ---
candidates_exclusive_modules <- jaccards.tb %>%
  filter(Classification %in% c("Dissimilar", "Similar"))

#Save result
vroom::vroom_write(candidates_exclusive_modules, file = file.path(output_dir, "candidates_exclusive_modules.csv"))

#Count by region and phenotype
candidates_exclusive_modules.c <- candidates_exclusive_modules %>%
  group_by(Region, Classification) %>%
  summarise(n = n(), .groups = "drop")

#Save result
vroom::vroom_write(candidates_exclusive_modules.c, file = file.path(output_dir, "n_candidates_exclusive_modules.csv"))

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

#Add number of genes per module
exclusive_modules <- exclusive_modules %>%
  left_join(module_gene_counts, by = "Module") %>%
  arrange(Region, Phenotype, desc(N_genes))

#Save results
vroom::vroom_write(exclusive_modules, file = file.path(output_dir, "exclusive_modules_AD_CTRL.csv"))

#check

exclusive_summary <- summary_all %>% filter(Classification %in% c("AD_exclusive", "Control_exclusive"))

#Save results
vroom::vroom_write(exclusive_summary, file = file.path(output_dir, "exclusive_summary.csv"))

################################## PART 2 ################################## 

#Phenotypic disagreement
#Are there AD-exclusive modules that are present (conserved) in several regions?
#Are there Control-exclusive modules that also appear in multiple regions?
#You take the AD-exclusive modules from each region, for example:
#
# Mayo_CRB_AD_, Mayo_PCC_AD_3, etc.
# 
# And then you compare each pair between regions using, for example, the Jaccard index between their genes.
# 
# If two modules from different regions have a Jaccard index ≥ threshold (e.g. 0.8), you consider them ‘conserved between regions’.

#Goal: Identify exclusive modules (AD or Control) that appear in more than one region with high similarity (Jaccard_Index ≥ 0.8).

#Jaccard between exclusive modules in different regions (within each phenotype)

find_conserved_modules <- function(exclusive_df, modules_list, threshold = 0.7) {
  conserved_list <- list()
  #Filter per phenotype
  phenos <- unique(exclusive_df$Phenotype)
  for (pheno in phenos) {
    sub_df <- exclusive_df %>% filter(Phenotype == pheno)
    
    #Group per modules in diff regions
    for (i in 1:(nrow(sub_df) - 1)) {
      for (j in (i + 1):nrow(sub_df)) {
        mod1 <- sub_df$Module[i]
        mod2 <- sub_df$Module[j]
        
        reg1 <- sub_df$Region[i]
        reg2 <- sub_df$Region[j]
        
        #Skip if it is the same region
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

#Run
conserved_modules_across_regions <- find_conserved_modules(
  exclusive_modules,
  modules,
  threshold = upp_thres
)

conserved_modules_across_regions

#Save
vroom::vroom_write(conserved_modules_across_regions, file = file.path(output_dir, "conserved_exclusive_modules_across_regions.csv"))

#Generar un grafo de módulos conservados entre regiones

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

################################## PART 2 ################################## 

#Phenotypic disagreement
#Are there AD-exclusive modules that are present (conserved) in several regions?
#Are there Control-exclusive modules that also appear in multiple regions?
#You take the AD-exclusive modules from each region, for example:
#
# Mayo_CRB_AD_, Mayo_PCC_AD_3, etc.
# 
# And then you compare each pair between regions using, for example, the Jaccard index between their genes.
# 
# If two modules from different regions have a Jaccard index ≥ threshold (e.g. 0.8), you consider them ‘conserved between regions’.

#Goal: Identify exclusive modules (AD or Control) that appear in more than one region with high similarity (Jaccard_Index ≥ 0.8).

#Jaccard between exclusive modules in different regions (within each phenotype)

# --- Jaccard entre módulos exclusivos de distintas regiones (intra-fenotipo) ---

#Extraer módulos por fenotipo
phenos <- c("AD_exclusive", "Control_exclusive")
#phenos <- unique(exclusive_df$Phenotype)

jaccard_cross_region <- list()

for (pheno in phenos) {
  # Filtrar módulos del fenotipo actual
  mods_pheno <- exclusive_modules %>%
    filter(Phenotype == pheno)
  
  # Extraer los nombres de los módulos
  mod_names <- mods_pheno$Module
  
  # Filtrar del objeto `modules` los que están en esta lista
  mod_list <- modules[names(modules) %in% mod_names]
  # 
  # # Obtener todas las combinaciones de módulos de regiones distintas
  # combs <- combn(names(mod_list), 2, simplify = FALSE)
  # 
  # for (pair in combs) {
  #   mod1 <- pair[1]
  #   mod2 <- pair[2]
  #   
  #   region1 <- strsplit(mod1, "_")[[1]][1]
  #   region2 <- strsplit(mod2, "_")[[1]][1]
  #   
  #   # Solo comparar entre regiones diferentes
  #   if (region1 != region2) {
  #     genes1 <- mod_list[[mod1]]
  #     genes2 <- mod_list[[mod2]]
  #     
  #     jaccard_val <- jaccard_index(genes1, genes2)
  #     
  #     if (jaccard_val >= upp_thres) {
  #       jaccard_cross_region[[length(jaccard_cross_region) + 1]] <- data.frame(
  #         Phenotype = pheno,
  #         Module1 = mod1,
  #         Region1 = region1,
  #         Module2 = mod2,
  #         Region2 = region2,
  #         Jaccard_Index = jaccard_val,
  #         stringsAsFactors = FALSE
  #       )
  #    }
  #  }
  #}
}

# Consolidar resultados
jaccard_cross_region.tb <- rbindlist(jaccard_cross_region)
jaccard_cross_region.tb

# Guardar resultados
vroom::vroom_write(jaccard_cross_region.tb, file = file.path(output_dir, "exclusive_modules_conserved_between_regions.csv"))



