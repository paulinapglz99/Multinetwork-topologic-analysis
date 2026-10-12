# cd into wherever you want the -graph_top1000 PDFs/JPEGs to land
cd ~/Desktop/Multinetwork-topologic-analysis/figures/networks

NETWORKS_DIR=~/Desktop/Multinetwork-topologic-analysis/inputs/networks
TOPOLOGY_DIR=~/Desktop/Multinetwork-topologic-analysis/inputs/topology

for name in \
    Mayo_CRB_counts_AD_topN200000 \
    Mayo_CRB_counts_control_topN200000 \
    Mayo_TC_counts_AD_topN200000 \
    Mayo_TC_counts_control_topN200000 \
    ROSMAP_DLPFC_counts_AD_topN200000 \
    ROSMAP_DLPFC_counts_control_topN200000 \
    ROSMAP_HCN_counts_AD_topN200000 \
    ROSMAP_HCN_counts_control_topN200000 \
    ROSMAP_PCC_counts_AD_topN200000 \
    ROSMAP_PCC_counts_control_topN200000
do
    echo ">>> Plotting ${name}"
    Rscript /home/tiamat/Desktop/Multinetwork-topologic-analysis/06_visualization/plot_network.R \
        --edgelist "${NETWORKS_DIR}/${name}.tsv" \
        --nodes    "${TOPOLOGY_DIR}/${name}_nodes_summary.csv"\
	--top 500
done
