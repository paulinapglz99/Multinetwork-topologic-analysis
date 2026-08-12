# Multiscale co-expression network analysis reveals region-specific topological reorganization uncoupled from transcript abundance in Alzheimer's disease pathology

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R ≥ 4.0](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

This repository contains the analysis pipeline accompanying the manuscript. Scripts cover preprocessing of multi-cohort RNA-seq data, co-expression network construction and filtering, multiscale topological analysis (global, local, and mesoscale), community detection, hub identification, and functional enrichment — across five brain regions in normal aging and Alzheimer's disease (AD).

---

## Datasets

Raw expression and metadata were obtained from the following AMP-AD cohorts, available through the [AD Knowledge Portal](https://adknowledgeportal.synapse.org/):

| Cohort | Brain regions | Access |
|--------|--------------|--------|
| ROSMAP | DLPFC | Synapse: [syn3388564](https://www.synapse.org/#!Synapse:syn3388564) |
| Mayo Clinic | TC, CRB | Synapse: [syn5550404](https://www.synapse.org/#!Synapse:syn5550404) |

> Raw data are not included in this repository. Access requires registration at the AD Knowledge Portal.

---

## Repository structure

```
.
├── 00_preprocessing/         # Cohort-specific expression matrix preparation
├── 01_network_construction/  # Network filtering, format conversion, and test data generation
├── 02_topology/              # Global and local topological metrics
├── 03_community_detection/   # Module identification and functional enrichment
├── 04_hub_analysis/          # Hub gene detection and GABAergic rewiring
├── 05_expression/            # Differential expression and module eigengene summaries
│   └── exploratory/          # Scratch scripts, not part of the reproducible pipeline
├── 06_visualization/         # Figures for global metrics and Jaccard similarity matrices
├── bash/                     # Shell wrappers, run from the repository root
├── config/
│   └── paths.R               # Single source of truth for INPUTS_DIR/OUTPUTS_DIR/PLOTS_DIR
├── inputs/                   # (gitignored) real pipeline data — see "Setup" below
├── outputs/                  # (gitignored) intermediate results, regenerable
├── plots/                    # Final manuscript figures (tracked) — see "Figures" below
└── test_data/                # Small synthetic networks for pipeline testing
    └── outputs/               # Outputs of a test run on test_data/ (tracked, small)
```

All scripts are run from the repository root, e.g. `Rscript 06_visualization/3.plots_global.R`,
and locate their I/O via `config/paths.R` (no absolute or machine-specific paths).

---

## Setup: populating `inputs/`

`inputs/` is gitignored — real pipeline data must be placed there locally before running the
scripts that produce the manuscript figures. Expected layout:

```
inputs/
├── networks/             # filtered co-expression edgelists (*.tsv), one per region × phenotype
├── topology/             # networks_summary.csv + *_nodes_summary.csv (from 02_topology, --per_node)
├── enrichment/           # *_enrichment.csv (from 03_community_detection/community_enrichment_MF_BP_KEGG.R)
├── DEGS/                 # limma_*_allGenes.tsv differential expression tables (one per region)
├── expression_matrices/  # per-region *_counts_{AD,control}.tsv + metadata .txt (05_expression)
└── graphs/               # translated (gene-symbol) .graphml networks, used by gabaergic_rewiring.R
```

`00_preprocessing/*.R` and `05_expression/metadata_demographics.R` read raw counts/metadata
directly from a Synapse-mounted path (`/datos/rosmap/...`) — this data is access-controlled
(AD Knowledge Portal registration required) and is never copied into `inputs/`.

---

## Analysis pipeline

Scripts are numbered and organized by stage. Run them in the following order:

### 00 — Preprocessing

Cohort-specific scripts that normalize and filter expression matrices, retain samples with complete CERAD/diagnosis metadata, and write per-region expression tables used downstream.

| Script | Cohort |
|--------|--------|
| `2.pre-pro-mRNA_ROSMAP.R` | ROSMAP (DLPFC) |
| `2.pre-pro-mRNA_Mayo.R` | Mayo (TC, CRB) |

### 01 — Network construction

| Script | Description |
|--------|-------------|
| `cut_networks.R` | Filters networks by mutual information threshold |
| `cut_networks_n.R` | Variant with n-based MI threshold |
| `translate_graphs.R` | Converts Ensembl IDs to gene symbols across network files |
| `generate_test_data.R` | Generates synthetic edge lists for pipeline testing |

Run with: `bash/run_cut.sh`, `bash/run_translate.sh`

### 02 — Topology

| Script | Description |
|--------|-------------|
| `1.network_topology.R` | Core script — computes global, local, and mesoscale metrics in parallel (CLI) |
| `degree_distributions.R` | Degree distribution fitting and plots |
| `edge_threshold_sensitivity.R` | Sensitivity analysis across MI thresholds |
| `null_distributions.R` | Null model generation for metric comparison |
| `centrality_comparison_AD_vs_control.R` | Per-gene degree, PageRank, and betweenness in paired AD/control networks; percentile-rank deltas and Spearman correlation per region |

Run with: `bash/run_topos.sh`, `bash/run_degree.sh`, `bash/run_edge_threshold_s.sh`, `bash/run_null_dis.sh`

#### CLI usage for `1.network_topology.R`

```bash
Rscript 02_topology/1.network_topology.R \
  --input_dir data/networks \
  --pattern ".*\\.csv$" \
  --out_dir results \
  --workers 4 \
  --per_node TRUE \
  --type edgelist
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i, --input_dir` | text | *(required)* | Folder with network files |
| `-p, --pattern` | text | `.*\.(txt\|tsv\|csv)$` | Regex for file selection |
| `-o, --out_dir` | text | `results` | Output folder |
| `-w, --workers` | integer | `2` | Parallel processes |
| `--per_node` | flag | `FALSE` | Save per-node metrics |
| `--make_html` | flag | `FALSE` | Generate HTML report |
| `--percol_steps` | integer | `51` | Steps for percolation simulation |
| `--seed` | integer | `42` | Random seed |
| `--type` | text | `auto` | Input format: `auto`, `graphml`, `edgelist`, `adjacency` |

### 03 — Community detection

| Script | Description |
|--------|-------------|
| `6.core_modules.R` | Louvain community detection; computes NMI and Jaccard similarity across networks |
| `community_enrichment_MF_BP_KEGG.R` | GO (MF, BP) and KEGG enrichment per module via clusterProfiler |

Run with: `bash/run_comm_enr.sh`, `bash/run_comm_enr_2net.sh`

### 04 — Hub analysis

| Script | Description |
|--------|-------------|
| `4.wide_hubs.R` | Identifies hub genes by degree and PageRank across all networks |
| `GABA_ego_analysis.R` | Ego network analysis centered on GABAergic marker genes |
| `gabaergic_rewiring.R` | Quantifies connectivity rewiring of GABAergic genes in AD vs. controls |

### 05 — Expression

| Script | Description |
|--------|-------------|
| `metadata_demographics.R` | Summarizes cohort demographics and CERAD score distributions |
| `explore_DEGS.R` | Differential expression analysis and overlap with network hubs |
| `module_expression_summary.R` | Computes module eigengenes and correlates them with clinical variables |
| `prepare_metadata.R` | Harmonizes sample metadata (region, sex) across cohorts |
| `dge_limma_AD_vs_control.R` | limma-trend differential expression (AD vs. control, sex-adjusted) per region |
| `gsea_reactome.R` | GSEA on limma t-statistics against Reactome gene sets; exports full, significant, and top terms per region |
| `dge_degree_integration.R` | Joins differential expression with degree percentile changes per region |

Run with: `bash/run_module_expression_summary.sh`

### 06 — Visualization

| Script | Description |
|--------|-------------|
| `3.plots_global.R` | Plots global topological metrics across regions and phenotypes |
| `raw_jaccard_matrices.R` | Plots Jaccard similarity heatmaps across network pairs |
| `degree_AD_vs_control_panels.R` | Gene-overlap barplot plus per-region degree scatterplots (AD vs. control) with Spearman rho |
| `volcano_by_region.R` | Volcano plots of AD vs. control differential expression, one row per region |
| `gsea_heatmap.R` | NES heatmap of top Reactome terms per region with FDR stars |
| `dge_gsea_panel.R` | Composite figure combining the volcano rows and the GSEA heatmap |
| `dge_degree_scatter.R` | log2 fold-change vs. degree percentile change for differentially expressed genes across regions |

Run with: `bash/run_jaccard.sh`

---

## Figures

Final manuscript figures are written to `plots/`, prefixed with their number in the manuscript.
Everything else these scripts produce (secondary/exploratory plots) also lands in `plots/`,
unprefixed.

| Figure | File in `plots/` | Generated by |
|--------|-------------------|--------------|
| Fig. 1 | `Figure1_final_global_figure.jpeg` | `06_visualization/3.plots_global.R` |
| Fig. 2 | `Figure2_panel_module_counts.jpeg` | `03_community_detection/6.core_modules.R` |
| Fig. 3 | `Figure3_cnet_conserved_pairs.jpeg` | `03_community_detection/6.core_modules.R` |
| Fig. 4 | `Figure4_hub_heatmaps_logfc.jpeg` | `04_hub_analysis/4.wide_hubs.R` |
| Fig. 5 | `Figure5_degree_AD_vs_control_panels.{png,pdf}` | `06_visualization/degree_AD_vs_control_panels.R` |
| Fig. 6 | `Figure6_DGE_GSEA_panel.{png,pdf}` | `06_visualization/dge_gsea_panel.R` |
| Fig. S1 | `FigureS1_contributions_pca.jpeg` | `06_visualization/3.plots_global.R` |
| Fig. S2 | `FigureS2_jaccard_histogram.jpeg` | `03_community_detection/6.core_modules.R` |
| Fig. S3 | `FigureS3_DGE_degree_scatter.{png,pdf}` | `06_visualization/dge_degree_scatter.R` |

Fig. 6 is assembled by `dge_gsea_panel.R` from the `.rds` components written by
`volcano_by_region.R` and `gsea_heatmap.R`, so those two run first.

Fig. 3 additionally depends on two hardcoded conserved module pairs
(`DLPFC_AD_14__PCC_AD_25`, `DLPFC_Control_3__PCC_Control_2`) inside `6.core_modules.R` — if the
underlying community detection changes, these keys may need to be updated to match whatever
pairs `cnet_networks` actually contains.

---

## Requirements

- R ≥ 4.0
- Key packages:

```r
install.packages(c(
  "igraph", "tidyverse", "ggplot2", "ggsci",
  "future.apply", "optparse", "data.table",
  "jsonlite", "stringr", "rmarkdown",
  "vroom", "ggrepel", "patchwork", "cowplot", "msigdbr"
))

# Bioconductor
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DESeq2",
                       "limma", "biomaRt"))
```

---

## Citation

> Perez-Gonzalez et al. *Multiscale co-expression network analysis reveals region-specific topological reorganization uncoupled from transcript abundance in Alzheimer's disease pathology.* (in preparation)

---
