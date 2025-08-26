#!/usr/bin/env Rscript

library(igraph)
library(data.table)

set.seed(42)
outdir <- "test_data/"
# Create subfolders for each format
dirs <- list(
  edgelist = file.path(outdir, "edgelists"),
  adjacency = file.path(outdir, "adjacency_matrices"),
  graphml = file.path(outdir, "graphmls")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

#Generate 5 networks: 3 Erdos-Renyi, 7 Preferential Attachment
for(i in 1:10) {
  if (i <= 3) {
    #Erdos-Renyi G(n=100, p=0.05)
    g <- sample_gnp(100, 0.05, directed = FALSE)
  } else {
    #Preferential Attachment network
    g <- sample_pa(100, power = 1, m = 3, directed = FALSE)
  }
  
  #1.Edgelist
  edgelist <- as_data_frame(g, what = "edges")
  w <- runif(nrow(edgelist), min = 0.1, max = 1)
  edgelist$weight <- w
  fwrite(edgelist, file = file.path(dirs$edgelist, paste0("network_", i, "_edgelist.csv")),
         sep = ",", col.names = TRUE)
  
  #2. Adjacency matrix
  #Assign weight if missing
  if (!"weight" %in% edge_attr_names(g)) E(g)$weight <- 1
  adj <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  adj_df <- data.table(node = rownames(adj), adj)
  fwrite(adj_df, file = file.path(dirs$adjacency, paste0("network_", i, "_adjacency.csv")), sep = ",")
  
  #3. GraphML
  write_graph(g, file = file.path(dirs$graphml, paste0("network_", i, ".graphml")), format = "graphml")
}

cat("Sample files generated:\n")
cat("Edgelists:", dirs$edgelist, "\n")
cat("Adjacency matrices:", dirs$adjacency, "\n")
cat("GraphMLs:", dirs$graphml, "\n")

#END