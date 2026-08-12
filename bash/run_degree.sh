#!/usr/bin/env bash
# Run from the repository root: bash bash/run_degree.sh
# Degree distribution fitting on the real filtered networks (inputs/networks/).
Rscript 02_topology/degree_distributions.R \
  -i inputs/networks \
  -p ".*\\.tsv$" \
  -o outputs/degree_distributions \
  -w 2 \
  --make_html TRUE
