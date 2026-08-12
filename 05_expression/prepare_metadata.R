#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(purrr)
  library(stringr)
})

# ===================== PATHS ===================== #
in_dir  <- "data/expression"
out_rds <- "results/metadata_sex.rds"
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "\\.txt$", full.names = TRUE)

# ===================== READ ===================== #
read_metadata <- function(file) {
  tryCatch({

    df <- vroom(
      file,
      delim = "\t",
      col_types = cols(.default = "c"),
      progress = FALSE,
      trim_ws = TRUE
    )

    df %>%
      mutate(
        source_file = basename(file),
        dataset = str_remove(basename(file), "\\.txt$")
      )

  }, error = function(e) {
    message("Error reading: ", file)
    NULL
  })
}

metadata_list <- map(files, read_metadata) %>% compact()

# ===================== HARMONIZE ===================== #
metadata_final <- map_dfr(metadata_list, function(df) {

  sex_vec <- if ("sex" %in% colnames(df)) {
    df$sex
  } else if ("msex" %in% colnames(df)) {
    df$msex
  } else {
    rep(NA_character_, nrow(df))
  }

  region_vec <- if ("tissue" %in% colnames(df)) df$tissue else rep(NA_character_, nrow(df))
  id_vec     <- if ("specimenID" %in% colnames(df)) df$specimenID else rep(NA_character_, nrow(df))

  tibble(
    region = region_vec,
    specimenID = id_vec,
    sex = sex_vec
  )
})

metadata_final <- metadata_final %>%
  mutate(
    region = case_when(
      region == "cerebellum" ~ "CRB",
      region == "temporal cortex" ~ "TC",
      region == "dorsolateral prefrontal cortex" ~ "DLPFC",
      region == "posterior cingulate cortex" ~ "PCC",
      region == "Head of caudate nucleus" ~ "HCN",
      TRUE ~ region
    ),
    sex = case_when(
      sex %in% c("0", 0) ~ "female",
      sex %in% c("1", 1) ~ "male",
      sex %in% c("female", "Female", "F") ~ "female",
      sex %in% c("male", "Male", "M") ~ "male",
      TRUE ~ NA_character_
    )
  )

print(table(metadata_final$region, metadata_final$sex))

# ===================== SAVE ===================== #
saveRDS(metadata_final, out_rds)
