#!/usr/bin/env Rscript
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = F)) BiocManager::install("org.Hs.eg.db")

library(igraph)
library(data.table)
library(biomaRt)

set.seed(42)
outdir <- "test_data/"
dirs <- list(
  edgelist = file.path(outdir, "edgelists"),
  adjacency = file.path(outdir, "adjacency_matrices"),
  graphml = file.path(outdir, "graphmls")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Obtener 100 genes humanos de ENSEMBL
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
               mart = ensembl, 
               filter = "biotype", 
               values = "protein_coding")
# tomar 100 genes aleatorios
genes <- genes[sample(nrow(genes), 100), ]
gene_names <- genes$ensembl_gene_id  # o genes$external_gene_name

# Guardar universo de genes
writeLines(gene_names, file.path(outdir, "gene_universe.txt"))


# Generar 10 redes
for(i in 1:10) {
  if (i <= 3) {
    g <- sample_gnp(100, 0.05, directed = FALSE)
  } else {
    g <- sample_pa(100, power = 1, m = 3, directed = FALSE)
  }
  
  # asignar nombres de genes a los nodos
  V(g)$name <- gene_names
  
  # 1. Edgelist
  edgelist <- as_data_frame(g, what = "edges")
  w <- runif(nrow(edgelist), min = 0.1, max = 1)
  edgelist$weight <- w
  fwrite(edgelist, file = file.path(dirs$edgelist, paste0("network_", i, "_edgelist.csv")),
         sep = ",", col.names = TRUE)
  
  # 2. Adjacency matrix
  if (!"weight" %in% edge_attr_names(g)) E(g)$weight <- 1
  adj <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  adj_df <- data.table(node = rownames(adj), adj)
  fwrite(adj_df, file = file.path(dirs$adjacency, paste0("network_", i, "_adjacency.csv")), sep = ",")
  
  # 3. GraphML
  write_graph(g, file = file.path(dirs$graphml, paste0("network_", i, ".graphml")), format = "graphml")
}

cat("Sample files generated with ENSEMBL gene IDs:\n")
cat("Edgelists:", dirs$edgelist, "\n")
cat("Adjacency matrices:", dirs$adjacency, "\n")
cat("GraphMLs:", dirs$graphml, "\n")
