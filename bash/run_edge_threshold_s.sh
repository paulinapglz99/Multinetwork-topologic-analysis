#!/usr/bin/env bash
# Run from the repository root: bash bash/run_edge_threshold_s.sh
# NOTE: needs the *unfiltered* MI edgelists (pre top-N cut), not the filtered
# networks in inputs/networks/. These are not included in the repo (too
# large) — place them in inputs/raw_networks/ before running.
Rscript 02_topology/edge_threshold_sensitivity.R \
  --input_dir inputs/raw_networks \
  --out_dir outputs/sensitivity \
  --workers 8 \
  --thresholds "10000,50000,100000,200000,300000,500000"
