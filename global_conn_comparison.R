#!/usr/bin/env Rscript
# Compute Jaccard similarity between networks (based on edge sets)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")

pacman::p_load(
  "data.table", "optparse", "tools",
  "igraph", "pheatmap", "tidyverse")

#Argument parsing
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Directory with network files", metavar = "path"),
  make_option(c("-d", "--out_dir"), type = "character", default = "jaccard_output", help = "Output directory [default: %default]", metavar = "path"),
  make_option(c("-p", "--pattern"), type = "character", default = "\\.tsv$", help = "File name pattern to match [default: %default]", metavar = "regex"),
  make_option(c("-t", "--type"), type = "character", default = "auto", help = "Network file type: edgelist, adjacency, auto [default: %default]", metavar = "type"))

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
#opt$out_dir <- "~/Desktop/fomo_networks/results_conn_comparisons"
opt$input_dir <- "~/Desktop/local_work/fomo_networks/"
opt$pattern <- "\\.tsv$"

#Read_network function
read_network <- function(path, type = opt$type) {
  message("DEBUG - received type: ", paste0(type, collapse = ","))

  type <- match.arg(type)
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "graphml") {
    g <- igraph::read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    df <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
    
    if (type == "edgelist") {
      g <- if (ncol(df) >= 3) igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
      else igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      
    } else if (type == "adjacency") {
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      
    } else if (type == "auto") {
      looks_like_adj <- all(sapply(df[,-1, drop=FALSE], is.numeric))
      if (looks_like_adj) {
        rownames(df) <- make.unique(as.character(df[, 1]))
        mat <- as.matrix(df[, -1, drop=FALSE])
        storage.mode(mat) <- "numeric"
        g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
      } else {
        g <- if (ncol(df) >= 3) igraph::graph_from_data_frame(df[, 1:3], directed = FALSE)
        else igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      }
    }
    
  } else {
    stop("Unsupported file format: ", ext)
  }
  
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  return(g)
}

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

#Jaccard function
jaccard_edges <- function(g1, g2) {
  inter <- igraph::intersection(g1, g2, byname = TRUE)
  uni <- igraph::union(g1, g2, byname = TRUE)
  return(length(E(inter)) / length(E(uni)))
}

#Load networks
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("At least 2 network files are required.")

meta <- do.call(rbind, lapply(files, extract_info))
#Extract unique names from regions
meta$region <- gsub("^(Mayo_|ROSMAP_)", "", meta$region)
meta$Network <- paste0(meta$region, "_", meta$phenotype)
networks <- lapply(files, read_network)
names(networks) <- meta$Network 
regions <- unique(meta$region)

#Compute Jaccard similarity matrix
n <- length(networks)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(networks)

for (i in 1:n) {
  for (j in i:n) {
    sim <- jaccard_edges(networks[[i]], networks[[j]])
    jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- sim
  }
}

#Save results
#csv_out <- file.path(opt$out_dir, "jaccard_edge_matrix.csv")
#fwrite(as.data.table(jaccard_matrix, keep.rownames = TRUE), csv_out)

#Make matrix long format
jaccard_matrix.l <- as.data.frame(jaccard_matrix) %>% 
  rownames_to_column("network1") %>% 
  pivot_longer(-network1, names_to = "network2", values_to = "jaccard")

#Only comparisons within the same region
region <- jaccard_matrix.l %>% 
  # unir info de cada red
  left_join(meta,  by = c("network1" = "Network"))  %>% 
  dplyr::select(-filename) %>% 
  rename(region1 = region, pheno1 = phenotype) %>% 
  left_join(meta,  by = c("network2" = "Network")) %>% 
  dplyr::select(-filename) %>% 
  rename(region2 = region, pheno2 = phenotype) %>% 
  # quedarse con pares Control vs AD
  filter(region1 == region2,
         pheno1 != pheno2)  %>% 
  transmute(
    region = region1,
    phenotype = pheno1,   # Control o AD
    jaccard = jaccard
  ) %>% 
  distinct() 

heatmap.df <- region  %>% 
  select(region, jaccard) %>% 
  distinct() %>% 
  arrange(jaccard) %>% 
  mutate(region = factor(region, levels = region))

#Plot miniheatmap

miniheat <- ggplot(heatmap.df, aes(x = "AD vs Control", 
                                   y = region,
                                fill = jaccard
                                )) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(jaccard, 3)), color = "black", size = 4) +
  scale_fill_gradient(low = "cornflowerblue", high = "firebrick") +
  labs(
    x = "",
    y = "",
    fill = " ") +
  #coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(), 
    legend.position="none"
)

#Vis
miniheat

#Save histograms
ggsave(filename = file.path(opt$out_dir, "miniheat.jpeg"), 
       plot = miniheat, 
       width = 5, 
       height = 7)

#Barplot

jaccard_global_connect <- heatmap.df %>%
  mutate(region = factor(region, levels = region[order(jaccard)])) %>%
  ggplot(aes(x = jaccard, y = region, fill = jaccard)) +
  geom_col(width = 0.95,
           color = "black",        # color del borde de la barra
           linewidth = 0.4)+         # grosor del borde) +
  geom_text(aes(label = round(jaccard, 3)),
            hjust = 1.45,
            color = "white",
            size = 4,
            fontface = "bold") +
  scale_fill_gradient(low = "steelblue", high = "firebrick") +
  labs(title = "Global connectivity overlap", 
       x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
jaccard_global_connect

#Save plot
ggsave("jaccard_global_connect.jpeg",
       plot = jaccard_global_connect,
       device = "jpeg",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300
)

#Save plot
ggsave("jaccard_global_connect.pdf",
       plot = jaccard_global_connect,
       device = "pdf",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300
)



#Get sets of edges function
get_edge_set <- function(g) {
  edgelist <- igraph::as_edgelist(g)
  apply(edgelist, 1, function(e) paste(sort(e), collapse = "|"))
}

opt$input_dir <- "~/Desktop/local_work/fomo_networks/"
opt$out_dir  <- "~/Desktop/local_work/fomo_networks/plots"

#Get files
files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) < 2) stop("Se necesitan al menos 2 archivos.")

networks <- lapply(files, read_network)
names(networks) <- basename(files)

#Get edge sets
edge_sets <- lapply(networks, get_edge_set)

#Calculate Jaccard matrix

n <- length(edge_sets)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(edge_sets)

for (i in 1:n) {
  for (j in i:n) {
    a <- edge_sets[[i]]
    b <- edge_sets[[j]]
    jaccard <- if (length(union(a, b)) > 0) {
      length(intersect(a, b)) / length(union(a, b))
    } else {
      NA
    }
    jaccard_matrix[i, j] <- jaccard_matrix[j, i] <- jaccard
  }
}

#Save matrix and heatmap
csv_out <- file.path(opt$out_dir, "jaccard_edge_matrix.csv")
pdf_out <- file.path(opt$out_dir, "jaccard_edge_heatmap.pdf")

fwrite(as.data.table(jaccard_matrix, keep.rownames = TRUE), csv_out)
pheatmap(jaccard_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, main = "Jaccard de enlaces",
         filename = pdf_out)

message("✅ Matriz guardada en: ", csv_out)
message("✅ Heatmap guardado en: ", pdf_out)