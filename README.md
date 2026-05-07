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
├── 06_visualization/         # Figures for global metrics and Jaccard similarity matrices
├── bash/                     # Shell wrappers to run scripts in parallel across networks
├── test_data/                # Minimal synthetic networks for pipeline testing
└── results/                  # Output directory (generated at runtime)
```

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

Run with: `bash/run_module_expression_summary.sh`

### 06 — Visualization

| Script | Description |
|--------|-------------|
| `3.plots_global.R` | Plots global topological metrics across regions and phenotypes |
| `raw_jaccard_matrices.R` | Plots Jaccard similarity heatmaps across network pairs |

Run with: `bash/run_jaccard.sh`

---

## Requirements

- R ≥ 4.0
- Key packages:

```r
install.packages(c(
  "igraph", "tidyverse", "ggplot2", "ggsci",
  "future.apply", "optparse", "data.table",
  "jsonlite", "stringr", "rmarkdown"
))

# Bioconductor
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DESeq2"))
```

---

## Citation

> Perez-Gonzalez et al. *Multiscale co-expression network analysis reveals region-specific topological reorganization uncoupled from transcript abundance in Alzheimer's disease pathology.* (in preparation)

---
