ok <- pacman::p_load(
  "tidyverse",
  "clusterProfiler",
  "org.Hs.eg.db"
)
opt<-NULL
opt$type <- "edgelist"
  
read_NULLread_network <- function(path, type = opt$type) {
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

PCC <- read_network("~/Desktop/local_work/fomo_networks/ROSMAP_PCC_counts_AD_topN200000.tsv", type = opt$type)
universe <- as.character(names(V(PCC)))
  
AD_conect <- c("ENSG00000119801", "ENSG00000163399", "ENSG00000165280")
ctrl_conect <- c("ENSG00000100916","ENSG00000147650","ENSG00000115561",
                 "ENSG00000144827","ENSG00000104142","ENSG00000164172","ENSG00000181852")
all_genes <- c(AD_conect, ctrl_conect)

AD_enr <- enrichGO(gene = AD_conect,
         OrgDb = org.Hs.eg.db,
         universe = universe,
         keyType = 'ENSEMBL',
         readable = TRUE,
         ont = "MF",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.10)

cnetplot(AD_enr)

ctrl_enr <- enrichGO(gene = ctrl_conect,
                   OrgDb = org.Hs.eg.db,
                   universe = universe,
                   keyType = 'ENSEMBL',
                   readable = TRUE,
                   ont = "MF",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.10)

cnetplot(ctrl_enr)


all_genes_enr <- enrichGO(gene = all_genes,
                     OrgDb = org.Hs.eg.db,
                     universe = universe,
                     keyType = 'ENSEMBL',
                     readable = TRUE,
                     ont = "MF",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.10)

cnetplot(ctrl_enr)

