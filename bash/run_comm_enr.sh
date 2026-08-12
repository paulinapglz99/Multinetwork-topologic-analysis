#!/usr/bin/env bash
# Run from the repository root: bash bash/run_comm_enr.sh
# NOTE: previously pointed at "community_enrichment.R", which was renamed to
# community_enrichment_MF_BP_KEGG.R. Needs a node index CSV (columns:
# file,node_table) — e.g. the networks_nodes_index.csv produced by
# 02_topology/1.network_topology.R with --per_node — placed under inputs/topology/.
Rscript 03_community_detection/community_enrichment_MF_BP_KEGG.R \
  --index_file inputs/topology/networks_nodes_index.csv \
  --out_dir outputs/enrichment \
  --workers 2
