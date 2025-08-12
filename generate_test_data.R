#!/usr/bin/env Rscript

library(igraph)
library(data.table)

set.seed(42)
outdir <- "test_data/edges"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#Generate 3 Erdos-Renyi networks and 2 Preferential Attachment networks
for(i in 1:5) {
  if (i <= 3) {
    #Erdos-Renyi G(n=100, p=0.05)
    g <- erdos.renyi.game(100, 0.05, directed = FALSE)
  } else {
    #Preferential Attachment Barabasi-Albert (100 nodes, m=3 edges per new node)
    g <- sample_pa(100, power = 1, m = 3, directed = FALSE)
  }
  #Extract list of unweighted edges
  edgelist <- as_data_frame(g, what = "edges")
  
  #Generate random weights between 0.1 and 1 for each edge
  w <- runif(nrow(edgelist), min = 0.1, max = 1)
  edgelist$weight <- w
  
  #Save csv file without header
  fwrite(edgelist, file = file.path(outdir, paste0("network_", i, ".csv")), sep = ",", col.names = FALSE)
}

cat("Sample files generated in", outdir, "\n")

#END