#!/usr/bin/env bash
# Run from the repository root: bash bash/run_comm_enr_2net.sh
# Reproduces community enrichment on the small synthetic test networks.
# NOTE: previously passed "-u test_data/gene_universe.txt", but the current
# script computes the gene universe internally from the index file's node
# tables — there is no --universe/-u option anymore, so it was dropped here.
Rscript 03_community_detection/community_enrichment_MF_BP_KEGG.R \
  --index_file test_data/outputs/topology/networks_nodes_index.csv \
  --out_dir test_data/outputs/enrichment \
  --workers 2
