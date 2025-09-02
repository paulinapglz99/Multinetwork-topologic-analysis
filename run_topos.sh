Rscript ~/Multinetwork-topologic-analysis/network_topology.R \
-i /datos/rosmap/multiregion_networks/all_nets \
-p ".*\\.graphml$" \
-o /datos/rosmap/multiregion_networks/all_nets/results_degree \
-w 2 \
--make_html TRUE

