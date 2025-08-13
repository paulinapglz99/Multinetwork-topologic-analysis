# Multinetwork-classic-analysis

Sometimes we want to analyse many networks at once. This repository performs a classic topological analysis of networks in parallel for several networks.
This script processes edge lists from networks and calculates global metrics, node metrics, and percolation simulations.
It can be run from the console with different parameters to customise its operation.

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
    
    
# Packages needed

* igraph
* data.table
* future.apply
* ggplot2
* rmarkdown
* tidyverse

# How to run

| Option            | Type       | Default                | Description                                            |
| ----------------- | ---------- | ---------------------- | ------------------------------------------------------ |
| `-i, --input_dir` | **texto**  | *(Mandatory )*        | Folder containing edge list files.                   |
| `-p, --pattern`   | **texto**  | `.*\.(txt\|tsv\|csv)$` | Regular expression for selecting files.           |
| `-o, --out_dir`   | **texto**  | `results`              | Folder where results are stored.                      |
| `-w, --workers`   | **entero** | `2`                    | Number of parallel processes requested.                      |
| `--per_node`      | **flag**   | `FALSE`                | If TRUE, save metrics per node (CSV per network). |
| `--make_html`     | **flag**   | `FALSE`                | If TRUE, generate an HTML report with `rmarkdown`. |
| `--percol_steps`  | **entero** | `51`                   | Number of steps in the percolation simulation. This greatly increases RAM usage, so proceed with caution.      |
| `--seed`          | **entero** | `42`                   | Reproducibility seed                         |


## Basic excecution

Rscript script.R -i path/to/directory/ -o results_dir/

## Example of running

1. When data is an edgelist in .csv

Rscript analyze_networks.R --input_dir data/edges --pattern "*.csv" --out_dir results --workers 4 --per_node TRUE --make_html TRUE

2.When network is in graphml format with extension .graphml

Rscript analyze_networks.R --input_dir test_data/graphmls --pattern "*.graphml" --out_dir results_graphml --workers 4 --per_node TRUE --make_html TRUE

## Help

Rscript script.R --help



