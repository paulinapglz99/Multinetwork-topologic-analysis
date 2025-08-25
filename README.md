# Multinetwork-topology-analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  
[![R](https://img.shields.io/badge/R-%3E=4.0.0-blue.svg)](https://www.r-project.org/)

A command-line tool for performing **parallel topological analysis of multiple networks**.  
It supports different input formats (GraphML, edge lists, adjacency matrices) and computes **global, local, and mesoscale metrics**, as well as **percolation simulations**.

---

## Features

- **Parallel processing** of multiple networks.
- **Automatic file type detection**:
  - GraphML
  - CSV/TSV adjacency matrices
  - CSV/TSV edge lists (2–3 columns)
- **Weighted graph support** (if weights are provided).
- **Handles loops and multiple edges**.
- **Optional HTML reporting** via R Markdown.

---

## Metrics Calculated

### Global
- Number of nodes
- Number of edges
- Average path length
- Global density
- Diameter
- Size of giant component
- Assortativity
- Percolation threshold

### Local
- Degree distribution
- PageRank
- K-core index

### Mesoscale (Infomap partitioning)
- Number of components
- Clustering coefficient
- Q modularity

---

## Requirements

- [R (≥ 4.5.0)](https://cran.r-project.org/)
- Packages:
  - `igraph (2.1.4)`
  - `data.table`
  - `future.apply`
  - `ggplot2`
  - `tidyverse`
  - `jsonlite`
  - `stringr`
  - `optparse`
  - `tools`
  - `rmarkdown`

Install the dependencies in R:

```r
install.packages(c(
  "igraph", "data.table", "future.apply",
  "ggplot2", "tidyverse", "jsonlite",
  "stringr", "optparse", "tools", "rmarkdown"
))
````

---

## Usage

Run the script from the command line:

```bash
Rscript analyze_networks.R [options]
```

### Options

| Option            | Type    | Default                | Description                                                                 |
| ----------------- | ------- | ---------------------- | --------------------------------------------------------------------------- |
| `-i, --input_dir` | text    | *(mandatory)*          | Folder containing edge list / adjacency / GraphML files.                    |
| `-p, --pattern`   | text    | `.*\.(txt\|tsv\|csv)$` | Regular expression for selecting files.                                     |
| `-o, --out_dir`   | text    | `results`              | Folder where results are stored. Script writes it automatically if it doesn't exist|
| `-w, --workers`   | integer | `2`                    | Number of parallel processes (threads to use).                                               |
| `--per_node`      | flag    | `FALSE`                | If `TRUE`, saves per-node metrics (CSV per network).                        |
| `--make_html`     | flag    | `FALSE`                | If `TRUE`, generates an HTML report with `rmarkdown`.                       |
| `--percol_steps`  | integer | `51`                   | Number of steps for percolation simulation (high RAM usage, use with care). |
| `--seed`          | integer | `42`                   | Random seed for reproducibility.                                            |
| `--type`          | text | `auto`                   | Format type of the input network ('auto' (default), 'graphml','edgelist' or 'adjacency')|

Show help:

```bash
Rscript analyze_networks.R --help
```

<<<<<<< HEAD
* Handles loops and multiple edges automatically.
    
# Packages needed 

* igraph 2.1.4
* data.table
* future.apply
* ggplot2
* tidyverse
* jsonlite
* stringr
* optparse
* rmarkdown


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
| `--type`          | **texto** | `auto`                   | Graph format (Input format: 'auto' (default), 'edgelist' or 'adjacency)                       |
=======
---

## Examples

### 1. Edge list in `.csv`

```bash
Rscript analyze_networks.R \
  --input_dir data/edges \
  --pattern "*.csv" \
  --out_dir results \
  --workers 4 \
  --per_node TRUE \
  --make_html TRUE
```
>>>>>>> 873b586536dd59adab1cf10a98b2917f6c697210

### 2. GraphML format

```bash
Rscript analyze_networks.R \
  --input_dir test_data/graphmls \
  --pattern "*.graphml" \
  --out_dir results_graphml \
  --workers 4 \
  --per_node TRUE \
  --make_html TRUE
```

### 3. Adjacency matrix in `.csv`

```bash
Rscript analyze_networks.R \
  --input_dir test_data/adjacency_matrices \
  --pattern "*.csv" \
  --out_dir results_adj_m \
  --workers 4 \
  --per_node TRUE \
  --make_html TRUE
```

---

## Notes

* **Adjacency matrices**:

  * Do **not** use a column with names.
  * Add a row with node/gene names at the top.
  * Example:

    ```
    V1,V2,V3
    0,0,0
    0,1,0
    0,0,1
    ```

* **Edge lists**:

  * Do **not** include a header row.

* **Mixed formats**:

  * Do not mix network formats in the same directory.
  * Use one directory per format if needed.

---

## Contributing

Contributions are welcome!

* Open an [issue](../../issues) for bug reports or feature requests.
* Submit a pull request with improvements.
* Please ensure code is well-documented and tested.

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Summary of Improvements

* **Added project description** at the top with badges for license and R support.
* **Reorganized sections** into Features, Metrics, Requirements, Usage, Examples, Notes, Contributing, License.
* **Added installation instructions** for required R packages.
* **Formatted metrics** into structured lists (Global, Local, Mesoscale).
* **Improved clarity** in examples (line breaks, code formatting).
* **Linked** to CRAN and GitHub issues for references.
* **Added placeholders** where details may need to be filled (e.g., LICENSE file link).

This structure follows [GitHub’s best practices](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes) and common patterns in open source README files.

```

Would you like me to also add a **"Project Structure"** section (explaining what each file/script does) or keep it minimal?
```
