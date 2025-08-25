Rscript ~/Multinetwork-classic-analysis/network_topology.R \
-i input_dir \
-p *.tsv \
-o results_dir \
-w 2 \
-- per_node TRUE \
--make_html TRUE \
--type edgelist
