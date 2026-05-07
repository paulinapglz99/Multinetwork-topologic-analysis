#!/usr/bin/env Rscript
# gabaergic_neighborhood_analysis.R
#
# Differential co-expression neighborhood analysis for GABAergic genes
# across AD and Control networks (multi-region).
#
# Goal: For each seed GABAergic gene, identify which genes it co-expresses
# with in AD but not in Control (and vice versa), revealing condition-specific
# transcriptional partnerships.
#
# Usage: Adapt the input paths and run interactively or as Rscript.
# See: https://github.com/paulinapglz99/Multinetwork-topologic-analysis

# ── 0. Packages ──────────────────────────────────────────────────────────────

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  igraph,
  tidyverse,
  ggsci,       # AAAS palette
  clusterProfiler,
  org.Hs.eg.db,
  ggrepel,
  patchwork
)

# ── 1. GABAergic seed gene list ──────────────────────────────────────────────
#
# Curated from:
#   - Krueger-Burg (2025) Trends in Neurosciences: synapse organizer proteins
#   - Williams et al. (2021) Brain Communications: synaptic DEGs in AD PFC
#   - Taddei & Duff (2025) eBioMedicine: synapse vulnerability in AD
#   - Internal results: GABRA3 is an AD-exclusive hub; SYT2, SV2C, NEFH
#     define the conserved trans-regional AD synaptic module (PCC–DLPFC)
#
# Organized by functional class for interpretability:

gabaergic_genes <- list(
  
  # GABA-A receptor subunits — ionotropic, primary mediators of fast inhibition
  GABAAR_subunits = c(
    "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6",
    "GABRB1", "GABRB2", "GABRB3",
    "GABRG1", "GABRG2", "GABRG3",
    "GABRD", "GABRE", "GABRP", "GABRQ"
  ),
  
  # Postsynaptic scaffolding — organizes GABA-AR clusters at the synapse
  Scaffolding = c(
    "GPHN",      # gephyrin: core GABAergic postsynaptic scaffold
    "ARHGEF9"    # collybistin: GEF that recruits gephyrin to membrane
  ),
  
  # Trans-synaptic adhesion proteins — define synapse identity and specificity
  # (Krueger-Burg 2025: these form circuit-specific organizer complexes)
  Adhesion = c(
    "NLGN2",     # neuroligin-2: core organizer, selective for perisomatic PV+ synapses
    "NLGN3",     # neuroligin-3: CCK+ and SST+ subtypes; AD-linked
    "NRXN1", "NRXN2", "NRXN3",  # neurexins: presynaptic, massive splice diversity
    "ERBB4",     # receptor tyrosine kinase enriched in PV+ interneurons
    "SLITRK3",   # GABAergic synaptogenesis, PV+ and chandelier cells
    "GRID1",     # GluD1: trans-synaptic with Cbln4, SST+ specificity
    "CBLN4",     # cerebellin-4: secreted by SST+ interneurons
    "IGSF9B",    # disinhibitory (IN→IN) synapses
    "CDH13",     # cadherin-13: PV+ and SST+ restricted; AD GWAS hit
    "CDH12",     # cadherin-12: CCK+ specificity
    "DAG1",      # dystroglycan: CCK+ basket cell synapses
    "MDGA1", "MDGA2",  # suppressors of synapse formation via Nlgn interaction
    "LHFPL4",    # GARLH: GABA-AR auxiliary subunit, cell-type specific clustering
    "SHISA7",    # GABA-AR auxiliary subunit, benzodiazepine modulation
    "IQSEC3",    # Arf-GEF at SST+ synapses; anxiety behavior
    "SRGAP2"     # regulates inhibitory synapse development
  ),
  
  # Presynaptic vesicle machinery — directly implicated in the conserved
  # trans-regional AD module (PCC AD-25 / DLPFC AD-14) in our network analysis
  Presynaptic_vesicle = c(
    "SYT1",      # synaptotagmin-1: Ca2+ sensor for fast release in forebrain
    "SYT2",      # synaptotagmin-2: Ca2+ sensor for fast release at GABAergic
    # inhibitory synapses, esp. cerebellum; downregulated in AD cortex
    "SV2A", "SV2B", "SV2C",  # synaptic vesicle glycoproteins; SV2C in AD-exclusive module
    "SNAP25",    # SNARE component; reduced in AD (M1 module in Williams 2021)
    "STX1A",     # syntaxin-1A: SNARE complex, downregulated in AD
    "VAMP1", "VAMP2",  # synaptobrevin: vesicle-SNARE, AD-affected
    "NEFH",      # heavy neurofilament: cytoskeletal, in conserved AD module;
    # differentially expressed in AD, linked to ALS and neurodegeneration
    "SYN1", "SYN2"    # synapsins: vesicle tethering; SYN2 is an AD-exclusive hub
    # in our dataset (DLPFC, HCN, TC)
  ),
  
  # Interneuron identity markers — used to track interneuron-specific programs
  Interneuron_markers = c(
    "PVALB",     # parvalbumin: PV+ interneurons, severely affected in AD
    "SST",       # somatostatin: dendritic-targeting interneurons
    "VIP",       # vasoactive intestinal peptide: disinhibitory interneurons
    "RELN",      # reelin: expressed by SST+ and other GABAergic cells; AD biomarker
    "CALB1",     # calbindin: interneuron subtype marker
    "CALB2"      # calretinin: interneuron subtype marker
  )
)

# Flatten to a single character vector for network queries
all_gabaergic <- unique(unlist(gabaergic_genes))
message(glue::glue("Total GABAergic seed genes: {length(all_gabaergic)}"))

# ── 2. Input paths — adapt to your directory structure ────────────────────────
#
# Expects one edgelist per region-phenotype, e.g.:
#   networks/DLPFC_AD.tsv, networks/DLPFC_control.tsv, ...
# Columns: Regulator, Target, MI (or similar weight column)

network_dir <- "~/Desktop/local_work/fomo_networks/"    # ← CHANGE THIS
out_dir     <- "~/Desktop/local_work/fomo_networks/results/gabaergic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

regions    <- c("DLPFC", "HCN", "PCC", "TC", "CRB")
phenotypes <- c("AD", "control")

# ── 3. Helper: read edgelist → igraph ─────────────────────────────────────────

read_network <- function(path) {
  read_tsv(path, col_types = cols(), show_col_types = FALSE) |>
    select(1:3) |>                        # Regulator, Target, weight
    rename(from = 1, to = 2, weight = 3) |>
    igraph::graph_from_data_frame(directed = FALSE) |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)
}

# ── 4. Core analysis: differential neighborhood per seed gene ─────────────────
#
# For each seed gene g and each region:
#   N_AD(g)      = set of neighbors in the AD network
#   N_CTL(g)     = set of neighbors in the Control network
#   AD-exclusive = N_AD(g) \ N_CTL(g)
#   CTL-exclusive= N_CTL(g) \ N_AD(g)
#   Shared       = N_AD(g) ∩ N_CTL(g)
#
# This reveals rewired transcriptional partnerships under pathological conditions.

extract_neighborhood <- function(graph, seed_gene) {
  if (!seed_gene %in% V(graph)$name) return(character(0))
  neighbors(graph, seed_gene)$name
}

run_neighborhood_analysis <- function(region, network_dir, seed_genes) {
  
  # Load the two networks for this region
  path_ad  <- file.path(network_dir, paste0(region, "_AD.tsv"))
  path_ctl <- file.path(network_dir, paste0(region, "_control.tsv"))
  
  if (!file.exists(path_ad) | !file.exists(path_ctl)) {
    warning(paste("Files not found for region:", region))
    return(NULL)
  }
  
  g_ad  <- read_network(path_ad)
  g_ctl <- read_network(path_ctl)
  
  # For each seed gene, compute neighborhood sets
  map_dfr(seed_genes, function(gene) {
    
    n_ad  <- extract_neighborhood(g_ad,  gene)
    n_ctl <- extract_neighborhood(g_ctl, gene)
    
    tibble(
      region          = region,
      seed_gene       = gene,
      present_in_AD   = gene %in% V(g_ad)$name,
      present_in_CTL  = gene %in% V(g_ctl)$name,
      n_AD_neighbors  = length(n_ad),
      n_CTL_neighbors = length(n_ctl),
      # AD-exclusive partners: wiring that only exists under pathology
      AD_exclusive_partners  = list(setdiff(n_ad, n_ctl)),
      # Control-exclusive partners: healthy wiring lost in AD
      CTL_exclusive_partners = list(setdiff(n_ctl, n_ad)),
      # Shared partners: conserved across conditions
      shared_partners        = list(intersect(n_ad, n_ctl))
    )
  })
}

# Run across all regions
neighborhood_results <- map_dfr(regions, run_neighborhood_analysis,
                                network_dir = network_dir,
                                seed_genes  = all_gabaergic)

# ── 5. Summary table ──────────────────────────────────────────────────────────

neighborhood_summary <- neighborhood_results |>
  mutate(
    n_AD_exclusive  = map_int(AD_exclusive_partners,  length),
    n_CTL_exclusive = map_int(CTL_exclusive_partners, length),
    n_shared        = map_int(shared_partners,        length),
    # Rewiring index: proportion of total neighborhood that changed
    rewiring_index  = (n_AD_exclusive + n_CTL_exclusive) /
      (n_AD_neighbors + n_CTL_neighbors + 0.01)
  )

write_tsv(neighborhood_summary |> select(-ends_with("_partners")),
          file.path(out_dir, "neighborhood_summary.tsv"))

# ── 6. Visualization ──────────────────────────────────────────────────────────

aaas_pal <- ggsci::pal_aaas()(6)

# 6A. Heatmap of rewiring index per gene × region
rewiring_matrix <- neighborhood_summary |>
  filter(present_in_AD | present_in_CTL) |>  # keep genes detected in at least one network
  select(region, seed_gene, rewiring_index) |>
  pivot_wider(names_from = region, values_from = rewiring_index, values_fill = NA)

p_heatmap <- neighborhood_summary |>
  filter(present_in_AD | present_in_CTL) |>
  ggplot(aes(x = region, y = seed_gene, fill = rewiring_index)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0.5,
    na.value = "grey90",
    name     = "Rewiring\nindex"
  ) +
  labs(
    title   = "GABAergic gene co-expression rewiring in Alzheimer's disease",
    subtitle = "Rewiring index = (AD-exclusive + Control-exclusive partners) / total neighborhood",
    x = "Brain region",
    y = "GABAergic seed gene"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y     = element_text(size = 7),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "heatmap_rewiring_index.pdf"),
       p_heatmap, width = 8, height = 14)

# 6B. Bar plot: total exclusive partners per gene, stacked by direction
exclusive_long <- neighborhood_summary |>
  filter(present_in_AD | present_in_CTL) |>
  select(region, seed_gene, n_AD_exclusive, n_CTL_exclusive) |>
  pivot_longer(
    cols      = c(n_AD_exclusive, n_CTL_exclusive),
    names_to  = "direction",
    values_to = "n_partners"
  ) |>
  mutate(
    direction = recode(direction,
                       n_AD_exclusive  = "AD-exclusive (gained)",
                       n_CTL_exclusive = "Control-exclusive (lost)"
    ),
    # Flip control-exclusive to negative for a diverging bar
    n_plot = if_else(str_detect(direction, "lost"), -n_partners, n_partners)
  )

p_diverging <- exclusive_long |>
  group_by(seed_gene, direction) |>
  summarise(n_plot = sum(n_plot, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = reorder(seed_gene, n_plot), y = n_plot, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(
    values = c("AD-exclusive (gained)"      = aaas_pal[2],
               "Control-exclusive (lost)"   = aaas_pal[1]),
    name = "Partnership direction"
  ) +
  coord_flip() +
  labs(
    title    = "Differential co-expression partnerships of GABAergic genes",
    subtitle = "Summed across all five brain regions",
    x        = NULL,
    y        = "Number of exclusive co-expression partners"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "barplot_exclusive_partners.pdf"),
       p_diverging, width = 9, height = 12)

# 6C. Scatter: rewiring index vs total degree (connectivity) in AD network
#     Genes that are both highly connected and highly rewired are the most
#     pathologically reconfigured GABAergic nodes

degree_ad <- neighborhood_summary |>
  select(region, seed_gene, n_AD_neighbors, rewiring_index) |>
  filter(!is.na(rewiring_index))

p_scatter <- degree_ad |>
  ggplot(aes(x = n_AD_neighbors, y = rewiring_index,
             color = region, label = seed_gene)) +
  geom_point(alpha = 0.7, size = 2) +
  ggrepel::geom_text_repel(
    data = . |> filter(rewiring_index > 0.7 | n_AD_neighbors > 500),
    size = 2.5, max.overlaps = 20, segment.color = "grey60"
  ) +
  scale_color_manual(values = ggsci::pal_aaas()(5), name = "Region") +
  scale_x_log10() +
  labs(
    title    = "Connectivity vs. rewiring of GABAergic genes in AD",
    subtitle = "Points in the upper-right quadrant represent highly connected, extensively rewired nodes",
    x        = "Degree in AD network (log\u2081\u2080)",
    y        = "Rewiring index"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "scatter_degree_vs_rewiring.pdf"),
       p_scatter, width = 9, height = 6)

# ── 7. GO enrichment on differential partners ─────────────────────────────────
#
# For each phenotype direction (AD-exclusive / CTL-exclusive), pool the
# differential partners of all GABAergic seed genes per region and run
# GO:BP enrichment against the region's full network gene universe.

run_enrichment <- function(partner_genes, universe_genes, label) {
  
  if (length(partner_genes) < 5) {
    message("Skipping enrichment for ", label, " — fewer than 5 genes.")
    return(NULL)
  }
  
  tryCatch({
    clusterProfiler::enrichGO(
      gene          = partner_genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.10,
      readable      = TRUE
    ) |>
      as.data.frame() |>
      mutate(label = label)
  }, error = function(e) {
    message("Enrichment failed for ", label, ": ", e$message)
    NULL
  })
}

# Flatten partner lists per region × direction
enrichment_results <- neighborhood_results |>
  filter(present_in_AD | present_in_CTL) |>
  group_by(region) |>
  summarise(
    AD_partners  = list(unique(unlist(AD_exclusive_partners))),
    CTL_partners = list(unique(unlist(CTL_exclusive_partners))),
    .groups = "drop"
  ) |>
  pmap_dfr(function(region, AD_partners, CTL_partners) {
    
    # Region universe: all genes present in either network
    # (load both networks again for the universe — adapt path as needed)
    path_ad  <- file.path(network_dir, paste0(region, "_AD.tsv"))
    path_ctl <- file.path(network_dir, paste0(region, "_control.tsv"))
    universe <- union(
      read_tsv(path_ad,  col_types = cols(), show_col_types = FALSE) |> pull(1:2) |> unlist() |> unique(),
      read_tsv(path_ctl, col_types = cols(), show_col_types = FALSE) |> pull(1:2) |> unlist() |> unique()
    )
    
    bind_rows(
      run_enrichment(AD_partners,  universe, paste0(region, "_AD-exclusive")),
      run_enrichment(CTL_partners, universe, paste0(region, "_CTL-exclusive"))
    )
  })

if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
  
  write_tsv(enrichment_results,
            file.path(out_dir, "enrichment_differential_partners.tsv"))
  
  # Dot plot of top enriched terms per region × direction
  p_enrich <- enrichment_results |>
    group_by(label) |>
    slice_min(p.adjust, n = 8) |>
    ungroup() |>
    separate(label, into = c("region", "direction"), sep = "_", extra = "merge") |>
    ggplot(aes(
      x    = reorder(Description, -p.adjust),
      y    = Count,
      size = Count,
      color = direction
    )) +
    geom_point() +
    scale_color_manual(
      values = c("AD-exclusive"  = aaas_pal[2],
                 "CTL-exclusive" = aaas_pal[1]),
      name = "Partner type"
    ) +
    facet_wrap(~region, scales = "free_y", ncol = 2) +
    coord_flip() +
    labs(
      title    = "Functional enrichment of differential co-expression partners",
      subtitle = "Top GO:BP terms for AD-exclusive and Control-exclusive GABAergic partners",
      x        = NULL,
      y        = "Gene count"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir, "dotplot_enrichment.pdf"),
         p_enrich, width = 14, height = 18)
}

# ── 8. Export full partner tables for downstream use ─────────────────────────

# AD-exclusive partners (long format, one row per seed gene × partner)
ad_exclusive_long <- neighborhood_results |>
  select(region, seed_gene, AD_exclusive_partners) |>
  unnest(AD_exclusive_partners) |>
  rename(partner_gene = AD_exclusive_partners) |>
  mutate(direction = "AD-exclusive")

ctl_exclusive_long <- neighborhood_results |>
  select(region, seed_gene, CTL_exclusive_partners) |>
  unnest(CTL_exclusive_partners) |>
  rename(partner_gene = CTL_exclusive_partners) |>
  mutate(direction = "CTL-exclusive")

all_exclusive <- bind_rows(ad_exclusive_long, ctl_exclusive_long)

write_tsv(all_exclusive,
          file.path(out_dir, "all_exclusive_partners_long.tsv"))

message("Done. All outputs written to: ", out_dir)