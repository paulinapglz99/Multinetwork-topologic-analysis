Rscript ~/Multinetwork-topologic-analysis/cut_networks.R \
  --input_dir ~/Multinetwork-topologic-analysis/test_data/edgelists \
  --out_dir ~/Multinetwork-topologic-analysis/test_data/edgelists/filtered \
  --pattern "\\.csv$" \
  --top_percent 0.10 \
  --workers 6
