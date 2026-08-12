#!/usr/bin/env bash
# Run from the repository root: bash bash/run_translate.sh
# Translates Ensembl IDs to gene symbols for the real filtered networks.
Rscript 01_network_construction/translate_graphs.R \
  -i inputs/networks \
  -o outputs/translated_graphs \
  -p "\\.tsv$" \
  -t edgelist
