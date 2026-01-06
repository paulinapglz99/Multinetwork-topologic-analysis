Rscript ~/Desktop/local_work/Multinetwork-topologic-analysis/network_topology.R \
-i ~/Desktop/local_work/Multinetwork-topologic-analysis/test_data/graphmls \
-p ".*\\.graphml$" \
-o ~/Desktop/local_work/Multinetwork-topologic-analysis/results_graphml \
-w 2 \
--per_node \
--make_html TRUE
