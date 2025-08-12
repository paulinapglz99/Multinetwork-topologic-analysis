# Multinetwork-classic-analysis

Sometimes we want to analyse many networks at once. This repository performs a classic topological analysis of networks in parallel for several networks.

Here we calculate

1. Global

    Number of nodes — size of the set of vertices.

    Number of edges — size of the set of edges.

    Average path length — average distance between all pairs.

    Global density — proportion of existing links to possible links.

    Diameter — largest minimum distance between two nodes.

    Size of giant component — fraction of nodes in the largest component.

    Assortativity — degree correlation between connected nodes.

    Percolation threshold — point at which the network loses global connectivity.

2. Local

    Degree distribution — frequency of each degree kk in the network.

    PageRank — relative importance of each node according to link flow.

    K-core index — maximum level of dense core in which the node is found.

3. Mesoscale (based on Infomap patitioning)

    Number of components — how many disconnected subsets there are.

    Clustering coefficient — local cohesion of neighbourhoods (can be averaged globally but is a measure of intermediate structure).

    Q modularity — quality of partitioning into communities.
