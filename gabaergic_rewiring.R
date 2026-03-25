#!/usr/bin/env Rscript
# gabaergic_neighborhood_analysis.R
#
# Differential co-expression neighborhood analysis for GABAergic genes
# across AD and Control networks (multi-region).
#
# Input: GraphML files in translated_graphs/
#        Vertex attrs: name (Ensembl ID), name_trad (HGNC symbol)
#        Naming convention: {cohort}_{region}_counts_{phenotype}_topN200000_translated.graphml
#
# See: https://github.com/paulinapglz99/Multinetwork-topologic-analysis

# ── 0. Packages ──────────────────────────────────────────────────────────────

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  igraph, tidyverse, ggsci, clusterProfiler, org.Hs.eg.db, ggrepel, glue
)

# ── 1. Paths ──────────────────────────────────────────────────────────────────

network_dir <- "~/Desktop/local_work/fomo_networks/translated_graphs"
out_dir     <- "~/Desktop/local_work/fomo_networks/results/gabaergic_neighborhood"
setwd(out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Cohort prefix per region (as observed in your filenames)
cohort_map <- c(DLPFC = "ROSMAP", HCN = "ROSMAP", PCC = "ROSMAP",
                TC    = "Mayo",   CRB  = "Mayo")

regions <- names(cohort_map)

# Build file path for a region-phenotype pair
network_path <- function(region, phenotype) {
  file.path(
    network_dir,
    glue("{cohort_map[region]}_{region}_counts_{phenotype}_topN200000_translated.graphml")
  )
}

# ── 2. Read GraphML using name_trad as working vertex name ────────────────────
#
# igraph uses the `name` vertex attribute as the canonical vertex identifier.
# Our files store Ensembl IDs in `name` and HGNC symbols in `name_trad`.
# We swap `name` for `name_trad` so all neighbor lookups use gene symbols.
# Vertices without a valid translation keep their Ensembl ID as fallback.

read_network <- function(region, phenotype) {
  path <- network_path(region, phenotype)
  
  if (!file.exists(path)) {
    warning(glue("File not found: {basename(path)}"))
    return(NULL)
  }
  
  g <- igraph::read_graph(path, format = "graphml") |>
    igraph::as_undirected(mode = "collapse") |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)
  
  # Promote name_trad to the working vertex name
  hgnc    <- vertex_attr(g, "name_trad")
  ensembl <- vertex_attr(g, "name")
  
  vertex_attr(g, "name") <- if_else(
    !is.na(hgnc) & nchar(hgnc) > 0, hgnc, ensembl
  )
  
  message(glue("Loaded {region} {phenotype}: {vcount(g)} nodes, {ecount(g)} edges"))
  g
}

# ── 3. GABAergic seed gene list ───────────────────────────────────────────────
#
# Curated from:
#   Krueger-Burg (2025) Trends Neurosci — GABAergic synapse organizer proteins
#   Williams et al. (2021) Brain Commun  — synaptic DEGs in AD PFC
#   Internal results: GABRA3 = AD-exclusive hub; SYT2, SV2C, NEFH define the
#   conserved trans-regional AD module (PCC AD-25 / DLPFC AD-14)

gabaergic_genes <- list(
  
  # Ionotropic GABA-A receptor subunits
  GABAAR_subunits = c(
    "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6",
    "GABRB1", "GABRB2", "GABRB3",
    "GABRG1", "GABRG2", "GABRG3",
    "GABRD", "GABRE", "GABRP", "GABRQ"
  ),
  
  # Postsynaptic scaffolding
  Scaffolding = c(
    "GPHN",    # gephyrin: core GABAergic postsynaptic scaffold
    "ARHGEF9"  # collybistin: GEF that anchors gephyrin to membrane
  ),
  
  # Trans-synaptic adhesion / synapse organizers
  # (Krueger-Burg 2025: circuit-specific organizer complexes)
  Adhesion = c(
    "NLGN2",              # neuroligin-2: universal GABAergic organizer, PV+ perisomatic
    "NLGN3",              # neuroligin-3: CCK+ and SST+ subtypes
    "NRXN1", "NRXN2", "NRXN3",  # neurexins: massive splice diversity
    "ERBB4",              # receptor TK enriched in PV+ interneurons
    "SLITRK3",            # GABAergic synaptogenesis, PV+ and chandelier cells
    "GRID1",              # GluD1: trans-synaptic via Cbln4, SST+ specificity
    "CBLN4",              # cerebellin-4: secreted by SST+ interneurons
    "IGSF9B",             # disinhibitory (IN -> IN) synapses
    "CDH13",              # cadherin-13: PV+ and SST+; AD GWAS locus
    "CDH12",              # cadherin-12: CCK+ specificity
    "DAG1",               # dystroglycan: CCK+ basket cell synapses
    "MDGA1", "MDGA2",     # suppressors of synapse formation
    "LHFPL4",             # GARLH: GABA-AR auxiliary subunit
    "SHISA7",             # GABA-AR auxiliary subunit, benzodiazepine modulation
    "IQSEC3",             # Arf-GEF at SST+ synapses
    "SRGAP2"              # inhibitory synapse development
  ),
  
  # Presynaptic vesicle machinery
  # SYT2, SV2C, NEFH: define the conserved trans-regional AD module in our data
  Presynaptic_vesicle = c(
    "SYT1",               # Ca2+ sensor for fast release, forebrain
    "SYT2",               # Ca2+ sensor for fast GABAergic release; downregulated in AD cortex
    "SV2A", "SV2B", "SV2C",  # synaptic vesicle glycoproteins; SV2C in AD-exclusive module
    "SNAP25",             # SNARE: reduced in AD M1 module (Williams 2021)
    "STX1A",              # syntaxin-1A: SNARE complex, downregulated in AD
    "VAMP1", "VAMP2",     # synaptobrevin: vesicle-SNARE
    "NEFH",               # heavy neurofilament: cytoskeletal; in conserved AD module
    "SYN1", "SYN2"        # synapsins; SYN2 is an AD-exclusive hub in our dataset
  ),
  
  # Interneuron identity markers
  Interneuron_markers = c(
    "PVALB",  # parvalbumin: fast-spiking PV+ interneurons, severely affected in AD
    "SST",    # somatostatin: dendritic-targeting interneurons
    "VIP",    # vasoactive intestinal peptide: disinhibitory circuit
    "RELN",   # reelin: expressed by SST+ cells; AD biomarker
    "CALB1",  # calbindin
    "CALB2"   # calretinin
  )
)

all_gabaergic <- unique(unlist(gabaergic_genes))
message(glue("Total GABAergic seed genes: {length(all_gabaergic)}"))

# ── 4. Load all networks once ─────────────────────────────────────────────────

graphs <- crossing(region = regions, phenotype = c("AD", "control")) |>
  mutate(
    graph = map2(region, phenotype, read_network),
    key   = paste(region, phenotype, sep = "_")
  ) |>
  filter(!map_lgl(graph, is.null))  # drop regions where files were missing

message(glue("\nLoaded {nrow(graphs)} networks successfully."))

# ── 5. Differential neighborhood analysis ────────────────────────────────────
#
# For each seed gene g and each region:
#   N_AD(g)       = neighbors in the AD network (by HGNC symbol)
#   N_CTL(g)      = neighbors in the Control network
#   AD-exclusive  = N_AD(g) \ N_CTL(g)   — new co-expression gained in AD
#   CTL-exclusive = N_CTL(g) \ N_AD(g)   — healthy co-expression lost in AD
#   Shared        = N_AD(g) ∩ N_CTL(g)   — conserved partnerships

get_neighbors <- function(g, gene) {
  if (is.null(g) || !gene %in% V(g)$name) return(character(0))
  neighbors(g, gene)$name
}

neighborhood_results <- graphs |>
  dplyr::select(region, phenotype, graph) |>
  pivot_wider(names_from = phenotype, values_from = graph) |>
  mutate(
    region_results = map2(AD, control, function(g_ad, g_ctl) {
      map_dfr(all_gabaergic, function(gene) {
        n_ad  <- get_neighbors(g_ad,  gene)
        n_ctl <- get_neighbors(g_ctl, gene)
        
        tibble(
          seed_gene              = gene,
          present_in_AD          = gene %in% V(g_ad)$name,
          present_in_CTL         = gene %in% V(g_ctl)$name,
          n_AD_neighbors         = length(n_ad),
          n_CTL_neighbors        = length(n_ctl),
          AD_exclusive_partners  = list(setdiff(n_ad, n_ctl)),
          CTL_exclusive_partners = list(setdiff(n_ctl, n_ad)),
          shared_partners        = list(intersect(n_ad, n_ctl))
        )
      })
    })
  ) |>
  dplyr::select(region, region_results) |>
  unnest(region_results)

# ── 6. Summary metrics ────────────────────────────────────────────────────────

class_map <- imap_dfr(gabaergic_genes,
                      ~tibble(seed_gene = .x, functional_class = .y))

neighborhood_summary <- neighborhood_results |>
  mutate(
    n_AD_exclusive  = map_int(AD_exclusive_partners,  length),
    n_CTL_exclusive = map_int(CTL_exclusive_partners, length),
    n_shared        = map_int(shared_partners,        length),
    rewiring_index  = (n_AD_exclusive + n_CTL_exclusive) /
      (n_AD_neighbors + n_CTL_neighbors + 0.01)
  ) |>
  left_join(class_map, by = "seed_gene")

write_tsv(
  neighborhood_summary |> dplyr::select(-ends_with("_partners")),
  file.path(out_dir, "neighborhood_summary.tsv")
)

message("Summary table written.")

# ── 7. Visualizations ─────────────────────────────────────────────────────────

aaas6 <- ggsci::pal_aaas()(6)

# 7A. Heatmap: rewiring index per gene x region, faceted by functional class
p_heatmap <- neighborhood_summary |>
  filter(present_in_AD | present_in_CTL) |>
  ggplot(aes(x = region, y = seed_gene, fill = rewiring_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#F2B9B9", mid = "#D72E2E", high = "#5B1111",
    midpoint = 0.5, na.value = "grey85",
    name = "Rewiring\nindex"
  ) +
  facet_grid(functional_class ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "GABAergic co-expression rewiring in Alzheimer's disease",
    subtitle = "Rewiring index = changed partnerships / total neighborhood",
    x = "Brain region", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y  = element_text(size = 7),
    panel.grid   = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
    plot.title   = element_text(face = "bold")
  )

#Vis

p_heatmap

ggsave(file.path(out_dir, "heatmap_rewiring_index.pdf"),
       p_heatmap, width = 7, height = 16)

# 7B. Diverging bar: AD-exclusive (gained) vs CTL-exclusive (lost),
#     summed across all five regions per seed gene

p_diverging <- neighborhood_summary |>
  filter(present_in_AD | present_in_CTL) |>
  group_by(seed_gene, functional_class) |>
  summarise(
    gained = sum(n_AD_exclusive,  na.rm = TRUE),
    lost   = sum(n_CTL_exclusive, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(c(gained, lost), names_to = "direction", values_to = "n") |>
  mutate(
    n_plot    = if_else(direction == "lost", -n, n),
    direction = recode(direction,
                       gained = "AD-exclusive (gained)",
                       lost   = "Control-exclusive (lost)"
    )
  ) |>
  ggplot(aes(x = reorder(seed_gene, abs(n_plot)), y = n_plot, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(
    values = c("AD-exclusive (gained)"    = aaas6[2],
               "Control-exclusive (lost)" = aaas6[1]),
    name = NULL
  ) +
  facet_grid(functional_class ~ ., scales = "free_y", space = "free_y") +
  coord_flip() +
  labs(
    title    = "Differential co-expression partnerships of GABAergic genes",
    subtitle = "Summed across all five brain regions",
    x = NULL, y = "Number of exclusive co-expression partners"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text.y    = element_text(angle = 0, hjust = 0, face = "bold"),
    plot.title      = element_text(face = "bold")
  )

#Vis
p_diverging

ggsave(file.path(out_dir, "barplot_exclusive_partners.pdf"),
       p_diverging, width = 9, height = 16)

# 7C. Scatter: AD-degree vs rewiring index per region
p_scatter <- neighborhood_summary |>
  filter(present_in_AD, n_AD_neighbors > 0) |>
  ggplot(aes(
    x     = n_AD_neighbors,
    y     = rewiring_index,
    color = functional_class,
    label = seed_gene
  )) +
  geom_point(alpha = 0.7, size = 2) +
  ggrepel::geom_text_repel(
    data = neighborhood_summary |> filter(rewiring_index > 0.7 | n_AD_neighbors > 300),
    size = 2.5, max.overlaps = 20, segment.color = "grey60"
  ) +
  scale_color_manual(values = ggsci::pal_aaas()(7), name = "Functional class") +
  scale_x_log10() +
  facet_wrap(~region, ncol = 3) +
  labs(
    title    = "Connectivity vs. rewiring of GABAergic genes in AD",
    subtitle = "Upper-right: highly connected and extensively rewired nodes",
    x        = expression(paste("Degree in AD network (log"[10], ")")),
    y        = "Rewiring index"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

#Vis
p_scatter

ggsave(file.path(out_dir, "scatter_degree_rewiring.pdf"),
       p_scatter, width = 12, height = 8)

# ── 8. GO:BP enrichment on differential partners ──────────────────────────────

run_enrichment <- function(partner_genes, universe_genes, label) {
  if (length(unique(partner_genes)) < 5) {
    message(glue("Skipping {label}: fewer than 5 genes.")); return(NULL)
  }
  tryCatch(
    clusterProfiler::enrichGO(
      gene          = unique(partner_genes),
      universe      = unique(universe_genes),
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.10,
      readable      = TRUE
    ) |> as.data.frame() |> mutate(analysis_label = label),
    error = function(e) { message(glue("Enrichment failed for {label}: {e$message}")); NULL }
  )
}

# Build universe per region from loaded igraph objects
region_universe <- graphs |>
  group_by(region) |>
  summarise(universe = list(unique(unlist(map(graph, ~V(.x)$name)))),
            .groups = "drop")

enrichment_results <- neighborhood_results |>
  filter(present_in_AD | present_in_CTL) |>
  group_by(region) |>
  summarise(
    AD_partners  = list(unique(unlist(AD_exclusive_partners))),
    CTL_partners = list(unique(unlist(CTL_exclusive_partners))),
    .groups = "drop"
  ) |>
  left_join(region_universe, by = "region") |>
  pmap_dfr(function(region, AD_partners, CTL_partners, universe) {
    bind_rows(
      run_enrichment(AD_partners,  universe, glue("{region}_AD-exclusive")),
      run_enrichment(CTL_partners, universe, glue("{region}_CTL-exclusive"))
    )
  })

if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
  write_tsv(enrichment_results,
            file.path(out_dir, "enrichment_differential_partners.tsv"))
  
  p_enrich <- enrichment_results |>
    separate(analysis_label, into = c("region", "direction"),
             sep = "_", extra = "merge") |>
    group_by(region, direction) |>
    slice_min(p.adjust, n = 8) |>
    ungroup() |>
    ggplot(aes(x = reorder(Description, -p.adjust),
               y = Count, size = Count, color = direction)) +
    geom_point() +
    scale_color_manual(
      values = c("AD-exclusive"  = aaas6[2],
                 "CTL-exclusive" = aaas6[1]),
      name = "Partner type"
    ) +
    facet_wrap(~region, scales = "free_y", ncol = 2) +
    coord_flip() +
    labs(
      title    = "Functional enrichment of differential co-expression partners",
      subtitle = "Top GO:BP terms — GABAergic gene neighborhoods",
      x = NULL, y = "Gene count"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir, "dotplot_enrichment.pdf"),
         p_enrich, width = 14, height = 18)
}

# ── 9. Export long-format partner tables ──────────────────────────────────────

bind_rows(
  neighborhood_results |>
    dplyr::select(region, seed_gene, AD_exclusive_partners) |>
    unnest(AD_exclusive_partners) |>
    rename(partner_gene = AD_exclusive_partners) |>
    mutate(direction = "AD-exclusive"),
  neighborhood_results |>
    dplyr::select(region, seed_gene, CTL_exclusive_partners) |>
    unnest(CTL_exclusive_partners) |>
    rename(partner_gene = CTL_exclusive_partners) |>
    mutate(direction = "CTL-exclusive")
) |>
  left_join(class_map, by = "seed_gene") |>
  write_tsv(file.path(out_dir, "all_exclusive_partners_long.tsv"))

message("\nDone. All outputs written to: ", out_dir)

# ── 10. Ego-network star plot for a gene of interest ─────────────────────────
#
# Draws a star-layout ego network for a single seed gene in a given region,
# showing AD vs Control neighbors side by side. Nodes are colored by their
# partnership category: AD-exclusive, Control-exclusive, or Shared.
#
# Change `focal_gene` and `focal_region` to the gene and region you want.

focal_gene   <- "SYN2"   # ← change to any gene in all_gabaergic
focal_region <- "DLPFC"    # ← change to: DLPFC, HCN, PCC, TC, or CRB

# Color palette (AAAS)
col_ad      <- aaas6[2]   # AD-exclusive partners
col_ctl     <- aaas6[1]   # Control-exclusive partners
col_shared  <- aaas6[3]   # Shared partners
col_seed    <- "#333333"   # focal gene node

# Helper: build a tidy edge + node table for one ego network condition
build_ego_tbl <- function(g, seed, partner_status) {
  # partner_status: named character vector gene -> "AD-exclusive" | "CTL-exclusive" | "Shared"
  if (is.null(g) || !seed %in% V(g)$name) return(NULL)
  
  nbrs <- neighbors(g, seed)$name
  
  nodes <- tibble(gene = c(seed, nbrs)) |>
    mutate(
      is_seed = gene == seed,
      status  = case_when(
        is_seed             ~ "Seed",
        TRUE                ~ partner_status[gene]
      ),
      status = replace_na(status, "Shared")   # genes not in the map are shared
    )
  
  edges <- tibble(from = seed, to = nbrs)
  
  list(nodes = nodes, edges = edges)
}

# Helper: star layout — seed at center (0,0), neighbors on a circle
star_layout <- function(n_nodes) {
  # First node = seed (center); the rest on a unit circle
  n_nbrs <- n_nodes - 1
  angles <- seq(0, 2 * pi, length.out = n_nbrs + 1)[-(n_nbrs + 1)]
  rbind(
    c(0, 0),
    cbind(cos(angles), sin(angles))
  )
}

# Retrieve the pre-computed partner sets for this gene x region
focal_row <- neighborhood_results |>
  filter(seed_gene == focal_gene, region == focal_region)

if (nrow(focal_row) == 0) {
  warning(glue("No data for {focal_gene} in {focal_region}. Check gene name and region."))
} else {
  
  ad_excl  <- focal_row$AD_exclusive_partners[[1]]
  ctl_excl <- focal_row$CTL_exclusive_partners[[1]]
  shared   <- focal_row$shared_partners[[1]]
  
  # Named status vector for all neighbors
  status_map <- c(
    setNames(rep("AD-exclusive",      length(ad_excl)),  ad_excl),
    setNames(rep("CTL-exclusive",     length(ctl_excl)), ctl_excl),
    setNames(rep("Shared",            length(shared)),   shared)
  )
  
  # Load the two graphs for this region
  g_ad  <- graphs |> filter(region == focal_region, phenotype == "AD")     |> pull(graph) |> pluck(1)
  g_ctl <- graphs |> filter(region == focal_region, phenotype == "control") |> pull(graph) |> pluck(1)
  
  # Build ego tables
  ego_ad  <- build_ego_tbl(g_ad,  focal_gene, status_map)
  ego_ctl <- build_ego_tbl(g_ctl, focal_gene, status_map)
  
  # ggplot2 star layout: one panel per condition, faceted side by side
  make_star_df <- function(ego_tbl, condition_label) {
    if (is.null(ego_tbl)) return(NULL)
    
    nodes <- ego_tbl$nodes
    edges <- ego_tbl$edges
    n     <- nrow(nodes)
    lay   <- star_layout(n)
    
    node_df <- nodes |>
      mutate(x = lay[, 1], y = lay[, 2], condition = condition_label)
    
    # Edges: join x,y for from and to
    seed_xy <- node_df |> filter(is_seed) |> dplyr::select(x0 = x, y0 = y)
    edge_df <- edges |>
      left_join(node_df |> dplyr::select(gene, x, y), by = c("to" = "gene")) |>
      mutate(
        x0        = seed_xy$x0,
        y0        = seed_xy$y0,
        condition = condition_label
      )
    
    list(nodes = node_df, edges = edge_df)
  }
  
  star_ad  <- make_star_df(ego_ad,  "AD")
  star_ctl <- make_star_df(ego_ctl, "Control")
  
  # Combine for faceted plot
  all_nodes <- bind_rows(star_ad$nodes, star_ctl$nodes)
  all_edges <- bind_rows(star_ad$edges, star_ctl$edges)
  
  # Node color and size scale
  status_colors <- c(
    "Seed"            = col_seed,
    "AD-exclusive"    = col_ad,
    "CTL-exclusive"   = col_ctl,
    "Shared"          = col_shared
  )
  status_sizes <- c(
    "Seed"            = 6,
    "AD-exclusive"    = 2.5,
    "CTL-exclusive"   = 2.5,
    "Shared"          = 2.5
  )
  
  # Count label for subtitle
  n_ad_excl  <- length(ad_excl)
  n_ctl_excl <- length(ctl_excl)
  n_shared   <- length(shared)
  
  p_star <- ggplot() +
    # Edges (spokes of the star)
    geom_segment(
      data = all_edges,
      aes(x = x0, y = y0, xend = x, yend = y),
      color = "grey75", linewidth = 0.3, alpha = 0.7
    ) +
    # Neighbor nodes
    geom_point(
      data = all_nodes |> filter(!is_seed),
      aes(x = x, y = y, color = status),
      size = 2.2, alpha = 0.85
    ) +
    # Seed node (on top)
    geom_point(
      data = all_nodes |> filter(is_seed),
      aes(x = x, y = y),
      color = col_seed, size = 6
    ) +
    # Seed label
    geom_text(
      data = all_nodes |> filter(is_seed),
      aes(x = x, y = y, label = gene),
      color = "white", fontface = "bold", size = 3
    ) +
    # Label a subset of neighbors: AD-exclusive and CTL-exclusive only
    # (labeling all shared neighbors would be unreadable)
    ggrepel::geom_text_repel(
      data = all_nodes |> filter(!is_seed, status != "Shared"),
      aes(x = x, y = y, label = gene, color = status),
      size          = 2,
      max.overlaps  = 30,
      segment.size  = 0.2,
      segment.color = "grey60",
      show.legend   = FALSE
    ) +
    scale_color_manual(
      values = status_colors,
      name   = "Co-expression status",
      guide  = guide_legend(override.aes = list(size = 4))
    ) +
    facet_wrap(~condition, ncol = 2) +
    coord_equal() +
    labs(
      title    = glue("Ego-network of {focal_gene} in {focal_region}"),
      subtitle = glue(
        "AD-exclusive: {n_ad_excl}  |  ",
        "Control-exclusive: {n_ctl_excl}  |  ",
        "Shared: {n_shared}"
      ),
      caption  = "Labeled nodes: AD-exclusive and Control-exclusive partners only"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 9),
      plot.caption    = element_text(hjust = 0.5, color = "grey60", size = 7),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold", size = 11),
      plot.margin     = margin(10, 10, 10, 10)
    )
  
  out_name <- glue("star_ego_{focal_gene}_{focal_region}.png")
  ggsave(file.path(out_dir, out_name), p_star,
         width = 14, height = 7)
  
  message(glue("Star plot saved: {out_name}"))
}

