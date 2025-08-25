Rscript ~/Multinetwork-topologic-analysis/network_topology.R \
-i /datos/rosmap/multiregion_networks/data/PCC/PCC_AD_output/full_net \
-p *.tsv \
-o /datos/rosmap/multiregion_networks/data/PCC/PCC_AD_output/full_net/results_dir \
-w 2 \
--per_node TRUE \
--make_html TRUE \
--type edgelist
