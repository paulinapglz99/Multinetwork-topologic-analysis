Rscript ~/Multinetwork-topologic-analysis/network_topology.R \
-i ~/Multinetwork-topologic-analysis/test_data/graphmls \
-p ".*\\.graphml$" \
-o ~/Multinetwork-topologic-analysis/results_graphml \
-w 2 \
--per_node \
--make_html TRUE

