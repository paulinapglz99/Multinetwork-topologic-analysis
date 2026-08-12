#!/usr/bin/env bash
# Run from the repository root: bash bash/run_module_expression_summary.sh
# NOTE: --input_dir expects the normalized expression matrices produced by
# 00_preprocessing/ (not included in the repo — Synapse-gated raw data).
# Place them locally under inputs/expression_matrices/ before running.
Rscript 05_expression/module_expression_summary.R \
  --input_dir inputs/expression_matrices \
  --module_topos inputs/topology \
  --out_dir outputs/module_expression
