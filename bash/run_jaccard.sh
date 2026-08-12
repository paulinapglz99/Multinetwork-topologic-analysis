#!/usr/bin/env bash
# Run from the repository root: bash bash/run_jaccard.sh
# NOTE: previously pointed at a "compare_networks_jaccard.R" that no longer
# exists in the repo; the current script is 06_visualization/raw_jaccard_matrices.R.
Rscript 06_visualization/raw_jaccard_matrices.R \
  --input_dir inputs/topology \
  --output_dir outputs/jaccard \
  --pattern "_nodes_summary\\.csv$"
