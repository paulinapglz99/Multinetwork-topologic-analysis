#!/usr/bin/env bash
# Run from the repository root: bash bash/run_null_dis.sh
# Null-model comparison on the small synthetic test networks.
Rscript 02_topology/null_distributions.R \
  -i test_data/graphmls/ \
  -o test_data/outputs/null_dist/ \
  -w 4 \
  --make_html
