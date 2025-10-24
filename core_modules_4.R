#!/usr/bin/env Rscript

# core_modules.R

# Packages
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
#opt$output_dir <- "/datos/rosmap/multiregion_networks/networks_final/networks_filtered/results_core_modules"

input_dir <- opt$input_dir
output_dir <- opt$output_dir
pattern <- opt$pattern

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Thresholds
upp_thres <- 0.8
low_thres <- 0.5

# --- Functions --- #

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

# Load modules with unique names
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

#Module pairing
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

#Compare and export matrices
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
    write.csv(jacmat, file = file.path(output_dir, paste0(name_out, "_jaccard_matrix.csv")), row.names = TRUE)
  }
}

#Logical classification of modules
module_classifier <- function(sim_mismo, sim_opuesto) {
  if (is.na(sim_mismo) | is.na(sim_opuesto)) return("Missing Data")
  if (sim_mismo > upp_thres & sim_opuesto < low_thres) return("Exclusive")
  if (sim_mismo < low_thres & sim_opuesto < low_thres) return("Region specific")
  if (sim_mismo < upp_thres & sim_opuesto > low_thres) return("Ambiguous")
  return("Shared or Unclassifiable")
}

# Llenar similitud con mismo fenotipo
# rellenar_sim_mismo <- function(df, tipo) {
#   jaccard_files_mismo <- list.files(output_dir, pattern = paste0(tipo, "vs", tipo, "_.*_jaccard_matrix\\.csv$"), full.names = TRUE)
#   
#   for (jf in jaccard_files_mismo) {
#     jac <- vroom(jf, show_col_types = FALSE)
#     first_col <- colnames(jac)[1]
#     jac <- column_to_rownames(jac, var = first_col)    filas <- rownames(jac)
#     columnas <- colnames(jac)
#     
#     for (i in which(df$fenotipo == tipo)) {
#       mod <- df$modulo[i]
#       if (mod %in% filas) {
#         df$sim_mismo[i] <- max(jac[mod, ], na.rm = TRUE)
#       } else if (mod %in% columnas) {
#         df$sim_mismo[i] <- max(jac[, mod], na.rm = TRUE)
#       }
#     }
#   }
#   
#   return(df)
# }

rellenar_sim_mismo <- function(df, tipo) {
  jaccard_files_mismo <- list.files(
    output_dir,
    pattern = paste0(tipo, "vs", tipo, "_.*_jaccard_matrix\\.csv$"),
    full.names = TRUE
  )
  
  for (jf in jaccard_files_mismo) {
    jac <- vroom(jf, show_col_types = FALSE)
    
    # Detectar si la primera columna está sin nombre o con ...1
    first_col <- colnames(jac)[1]
    if (first_col == "" || grepl("^\\.\\.\\.", first_col)) {
      message("⚠️ Primera columna sin nombre en: ", basename(jf), 
              " — usando la primera columna como rownames.")
      first_col <- colnames(jac)[1]
    }
    
    # Convertir la primera columna en rownames de forma segura
    jac <- column_to_rownames(jac, var = first_col)
    
    filas <- rownames(jac)
    columnas <- colnames(jac)
    
    # Recorrer filas de df que correspondan al fenotipo dado
    for (i in which(df$fenotipo == tipo)) {
      mod <- df$modulo[i]
      if (mod %in% filas) {
        df$sim_mismo[i] <- max(jac[mod, ], na.rm = TRUE)
      } else if (mod %in% columnas) {
        df$sim_mismo[i] <- max(jac[, mod], na.rm = TRUE)
      } else {
        message("❌ Módulo no encontrado en matriz: ", mod, " (", tipo, ") en ", basename(jf))
      }
    }
  }
  
  return(df)
}


# --- MAIN EXECUTION --- #

files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos dos archivos *_nodes_summary.csv.")

meta <- do.call(rbind, lapply(files, extract_info))
all_pairs <- build_pairs(meta)

compare_and_export(all_pairs$AD_vs_AD, "ADvsAD")
compare_and_export(all_pairs$Control_vs_Control, "ControlvsControl")
compare_and_export(all_pairs$AD_vs_Control, "ADvsControl")

message("Jaccard comparisons completed!")

# --- Clasificación lógica de módulos --- #

jaccard_files <- list.files(output_dir, pattern = "ADvsControl_.*_jaccard_matrix\\.csv$", full.names = TRUE)
tabla_modulos <- list()

for (jf in jaccard_files) {
  jac <- vroom(jf, show_col_types = FALSE)
  first_col <- colnames(jac)[1]
  jac <- column_to_rownames(jac, var = first_col)
  mod_AD <- rownames(jac)
  mod_Control <- colnames(jac)
  
  for (m in mod_AD) {
    sim_opuesto <- max(jac[m, ])
    tabla_modulos[[length(tabla_modulos) + 1]] <- data.frame(
      region = sub("ADvsControl_(.*?)_counts.*", "\\1", basename(jf)),
      modulo = m,
      fenotipo = "AD",
      sim_mismo = NA,
      sim_opuesto = sim_opuesto,
      clasificacion = NA
    )
  }
  
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

tabla_modulos_df <- do.call(rbind, tabla_modulos)
tabla_modulos_df <- rellenar_sim_mismo(tabla_modulos_df, "AD")
tabla_modulos_df <- rellenar_sim_mismo(tabla_modulos_df, "Control")

tabla_modulos_df$clasificacion <- mapply(module_classifier,
                                         tabla_modulos_df$sim_mismo,
                                         tabla_modulos_df$sim_opuesto)

write.csv(tabla_modulos_df, file.path(output_dir, "tabla_clasificacion_modulos.csv"), row.names = FALSE)

message("Clasificación terminada. Tabla guardada en: tabla_clasificacion_modulos.csv")
