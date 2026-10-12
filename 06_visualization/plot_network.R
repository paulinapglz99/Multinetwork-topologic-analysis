#!/usr/bin/env Rscript
# Plot one co-expression network colored by Louvain community.
# Usage: Rscript plot_network.R -e <edgelist.tsv> -n <nodes_summary.csv> [-t 1000]

library(optparse)
library(tidyverse)
library(igraph)
library(ggraph)
library(pals)  # glasbey/polychrome: many maximally-distinguishable colors

# --- CLI ---------------------------------------------------------------------
option_list <- list(
  make_option(c("-e", "--edgelist"), type = "character",
              help = "Path to edgelist TSV (gene1, gene2, MI)"),
  make_option(c("-n", "--nodes"), type = "character",
              help = "Path to nodes summary CSV with a `membership` column"),
  make_option(c("-t", "--top"), type = "integer", default = 1000,
              help = "Top-N edges by MI to keep for the plot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

edgelist_file <- opt$edgelist
nodes_file    <- opt$nodes
n_top_edges   <- opt$top

# Derive output stem from the edgelist filename so the graph inherits its name.
out_stem <- edgelist_file %>%
  basename() %>%
  tools::file_path_sans_ext() %>%
  paste0("-graph_top", n_top_edges)

# --- Read data ---------------------------------------------------------------
# Edgelist: gene1, gene2, MI. MI is kept as edge weight in case we want to map
# it to alpha/width later; for now it's just carried along.
edges <- read_tsv(edgelist_file, show_col_types = FALSE) %>%
  rename(from = 1, to = 2, weight = 3) %>%
  slice_head(n = n_top_edges)

# Node metadata: `membership` (Louvain community) drives node color. Filter to
# nodes that survive the edge prune, otherwise we'd get orphan points floating
# in the layout with no connections.
surviving_nodes <- unique(c(edges$from, edges$to))

nodes <- read_csv(nodes_file, show_col_types = FALSE) %>%
  rename(name = node) %>%
  filter(name %in% surviving_nodes) %>%
  mutate(membership = as.factor(membership))

# --- Build graph -------------------------------------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Drop self-loops and duplicated edges — ARACNe shouldn't produce them but
# it's cheap insurance before layout.
g <- igraph::simplify(g)

# Keep only the giant connected component. With top-N MI edges most of the
# graph is small orphan cliques that the stress layout stacks in a grid,
# producing vertical-stripe artifacts. The GCC is where the interpretable
# community structure lives.
comps  <- components(g)
gcc_id <- which.max(comps$csize)
g      <- induced_subgraph(g, V(g)[comps$membership == gcc_id])

# --- Palette -----------------------------------------------------------------
# ggsci::aaas caps at 10 colors, and these networks can hold dozens of
# communities. glasbey gives 32 maximally-distinguishable colors; we recycle
# instead of interpolating so contrast between neighbors stays high. Rebuild
# the level count from the surviving subgraph so we don't allocate colors for
# communities that got dropped with the orphan components.
surviving_membership <- factor(V(g)$membership)
n_comm               <- nlevels(surviving_membership)
comm_colors          <- rep(pals::glasbey(32), length.out = n_comm)

# --- Plot --------------------------------------------------------------------
# 'stress' layout separates communities well. With only ~1k edges we can raise
# edge alpha and node size compared to the full-network version.
p <- ggraph(g, layout = "stress") +
  geom_edge_link(color = "black", width = 1) +
  geom_node_point(aes(color = membership), size = 7) +
  scale_color_manual(values = comm_colors) +
  theme_void() +
  theme(legend.position = "none")

# --- Save --------------------------------------------------------------------
ggsave(paste0(out_stem, ".pdf"),  p, width = 8, height = 8, device = cairo_pdf)