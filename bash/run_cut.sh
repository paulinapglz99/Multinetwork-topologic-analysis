#!/usr/bin/env bash
# Run from the repository root: bash bash/run_cut.sh
# Filters synthetic test edgelists by top MI percentage.
Rscript 01_network_construction/cut_networks.R \
  --input_dir test_data/edgelists \
  --out_dir test_data/outputs/edgelists_filtered \
  --pattern "\\.csv$" \
  --top_percent 0.10 \
  --workers 6
