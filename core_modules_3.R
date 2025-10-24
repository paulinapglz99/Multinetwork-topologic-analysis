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
  vroom
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

print(meta)

#Get pairs! --- ---
all_pairs <- build_pairs(meta)
compare_and_export(all_pairs$AD_vs_AD, "ADvsAD")
compare_and_export(all_pairs$Control_vs_Control, "ControlvsControl")
compare_and_export(all_pairs$AD_vs_Control, "ADvsControl")

message("Comparisons DONE!. Files are now in: ", output_dir)

#Read and classify modules ---- --- 

#Thresholds definitions 
low_thres <- 0.8 #Up to this, I consider it a similar module
upp_thres <- 0.5 #Below this, I consider it a dissimilar module

#Function to classify modules by logic
module_classifier <- function(sim_mismo, sim_opuesto) {
  if (sim_mismo > upp_thres & sim_opuesto < low_thres) return("Exclusive")
  if (sim_mismo < low_thres & sim_opuesto < low_thres) return("Region specific")
  if (sim_mismo < upp_thres & sim_opuesto > low_thres) return("Ambiguous")
  return("Shared or Unclassifiable")
}

#Analyse the ADvsControl jaccards --- ---

jaccard_files <- list.files(output_dir, pattern = "ADvsControl_.*_jaccard_matrix\\.csv$", full.names = TRUE)
tabla_modulos <- list()

for (jf in jaccard_files) {
  jac <- vroom(jf) %>% 
    column_to_rownames(var = "...1")
  mod_AD <- rownames(jac)
  mod_Control <- colnames(jac)
  
  # Para cada módulo AD: similitud máxima con Control
  #Para cada módulo de AD: ¿a qué módulo de Control se parece más?
  #Para cada módulo de AD (m), buscamos el valor máximo en su fila → o sea, el módulo de Control al que más se parece (su mayor Jaccard).
  #Eso es su sim_opuesto.
  for (m in mod_AD) {
    sim_opuesto <- max(jac[m, ])
    #Guardamos la información en la tabla
    tabla_modulos[[length(tabla_modulos) + 1]] <- data.frame(
      region = sub("ADvsControl_(.*?)_counts.*", "\\1", basename(jf)),
      modulo = m,
      fenotipo = "AD",
      sim_mismo = NA,  # se llenará luego
      sim_opuesto = sim_opuesto,
      clasificacion = NA
    )
  }
  
  # Para cada módulo Control: similitud máxima con AD
  #Cada módulo de Control ve a cuál módulo de AD se parece más (su sim_opuesto).
  for (m in mod_Control) {
    sim_opuesto <- max(jac[, m])
    tabla_modulos[[length(tabla_modulos) + 1]] <- data.frame(
      region = sub("ADvsControl_(.*?)_counts.*", "\\1", basename(jf)),
      modulo = m,
      fenotipo = "Control",
      sim_mismo = NA,
      sim_opuesto = sim_opuesto,
      clasificacion = NA
    )
  }
}

#Combine results
tabla_modulos_df <- do.call(rbind, tabla_modulos)

# Llenar sim_mismo con los valores de ADvsAD o ControlvsControl
#Esta función toma la tabla (df) y el tipo de fenotipo ("AD" o "Control"), 
#y le va a buscar para cada módulo con qué otros módulos de su mismo fenotipo se parece más, pero en otras regiones.

rellenar_sim_mismo <- function(df, tipo) {
  jaccard_files_mismo <- list.files(output_dir, pattern = paste0(tipo, "vs", tipo, "_.*_jaccard_matrix\\.csv$"), full.names = TRUE)
  #Compara módulos con otros del mismo fenotipo
  for (jf in jaccard_files_mismo) {
    jac <- vroom(jf) %>% column_to_rownames(var = "...1")
    mods <- intersect(rownames(jac), colnames(jac))
    #🔍 Para cada módulo: buscamos su valor máximo
    #Para cada módulo m que aparece en esa matriz, buscamos:
    #Si está en la tabla df.
    #Si sí, le asignamos su similitud máxima con los módulos del mismo fenotipo (su fila de la matriz).
    for (m in mods) {
      if (m %in% df$modulo) {
        idx <- which(df$modulo == m & df$fenotipo == tipo)
        if (length(idx) > 0) {
          df$sim_mismo[idx] <- max(jac[m, ], na.rm = TRUE)
        }
      }
    }
  }
  return(df)
}

tabla_modulos_df <- rellenar_sim_mismo(tabla_modulos_df, "AD")
tabla_modulos_df <- rellenar_sim_mismo(tabla_modulos_df, "Control")

#Classify
tabla_modulos_df$clasificacion <- mapply(module_classifier,
                                         tabla_modulos_df$sim_mismo,
                                         tabla_modulos_df$sim_opuesto)

#Save tables
write.csv(tabla_modulos_df, file.path(output_dir, "tabla_clasificacion_modulos.csv"), row.names = FALSE)

