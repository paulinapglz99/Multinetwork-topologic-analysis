Rscript ~/Multinetwork-topologic-analysis/edge_threshold_sensitivity.R \
  --input_dir /STORAGE/csbig/rosmap_networks/networks_final \
  --out_dir /STORAGE/csbig/rosmap_networks/networks_final/sensitivity \
  --workers 8 \
  --thresholds "10000,50000,100000,200000,300000,500000"
