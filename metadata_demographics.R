#!/usr/bin/env Rscript

#Script: metadata_demographics.R

#Libraries
pacman::p_load("tidyverse", "vroom", "stringr")

#Functions
#Read tables from path
read.f <- function(path) {
  files <- list.files(path, full.names = TRUE)
  files <- files[grepl("\\.(csv|tsv|txt)$", files)]
  setNames(
    lapply(files, vroom::vroom, delim = NULL),
    basename(files)
  )
}

#Format mean ± sd
mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf("%.1f ± %.1f", mean(x), sd(x))
}

#Read input data

#Read metadata and RNA-seq matrices
metadatas <- read.f(path = "/home/tiamat/Desktop/local_work/matrices_originales/metadata")
matrices <- read.f(path = "/home/tiamat/Desktop/local_work/matrices_originales")

#Convert gene to rownames in expression matrices
matrices <- lapply(matrices, function(x) column_to_rownames(x, var = "gene"))

#Extract sample IDs from all matrices
matrix_ids <- lapply(matrices, colnames) %>% unlist() %>% unique()

#Preprocess metadata

#Standardize column names
standardize_colnames <- function(df) {
  df %>%
    rename_with(~ gsub("ageDeath", "age_death", .x)) %>%
    rename_with(~ gsub("apoeGenotype", "apoe_genotype", .x)) %>% 
    rename_with(~ gsub("^msex$", "sex", .x)) 
}

#Apply column name standardization and convert to character
metadatas <- lapply(metadatas, function(df) {
  df <- standardize_colnames(df)
  mutate(df, across(everything(), as.character))
})

#Combine metadata tables
metadata_all <- bind_rows(metadatas, .id = "source")

#Filter metadata for samples present in matrices
metadata_all <- metadata_all %>%
  filter(specimenID %in% matrix_ids) %>%
  mutate(
    age_death = str_remove(age_death, "\\+") %>% as.numeric(),
    age_at_visit_max = str_remove(age_at_visit_max, "\\+") %>% as.numeric(),
    pmi = as.numeric(pmi),
    educ = as.numeric(educ),
    sex = recode(sex, "1" = "male", "0" = "female"),
    ceradsc = recode(ceradsc, "1" = "Alzheimer's Disease", "4" = "Control"),
    diagnosis = recode(diagnosis, "Alzheimer Disease" = "Alzheimer's Disease", "control" = "Control"),
    diagnosis_final = coalesce(diagnosis, ceradsc)
  )

#Compute demographics grouped by tissue and diagnosis

demographics_dx <- metadata_all %>%
  group_by(tissue, diagnosis_final) %>%
  summarise(
    N = n_distinct(individualID),
    Age_death = mean_sd(age_death),
    PMI = mean_sd(pmi),
    # Education = mean_sd(educ),  # Uncomment if needed
    Male = sum(sex == "male", na.rm = TRUE),
    Female = sum(sex == "female", na.rm = TRUE),
    APOE4_carriers = sum(str_detect(apoe_genotype, "4"), na.rm = TRUE),
    .groups = "drop"
  )

#Output
vroom_write(file = "demographics.txt", demographics_dx)

