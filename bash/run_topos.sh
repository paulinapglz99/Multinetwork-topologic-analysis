#!/usr/bin/env bash
# Run from the repository root: bash bash/run_topos.sh
# Reproduces the topology pipeline on the small synthetic test networks.
Rscript 02_topology/1.network_topology.R \
  -i test_data/graphmls \
  -p ".*\\.graphml$" \
  -o test_data/outputs/topology \
  -w 2 \
  --per_node \
  --make_html TRUE
