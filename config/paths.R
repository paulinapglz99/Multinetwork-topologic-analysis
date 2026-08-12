# config/paths.R
# Single source of truth for repo-relative I/O locations.
# Assumes scripts are run from the repository root, e.g.:
#   Rscript 06_visualization/3.plots_global.R
#
# inputs/  — real pipeline data (gitignored). Populate locally; see README.md
#            for the expected layout (networks/, topology/, enrichment/, DEGS/, graphs/).
# outputs/ — intermediate results written by the pipeline (gitignored, regenerable).
# plots/   — final manuscript figures (tracked in git).

REPO_ROOT   <- normalizePath(".")
INPUTS_DIR  <- file.path(REPO_ROOT, "inputs")
OUTPUTS_DIR <- file.path(REPO_ROOT, "outputs")
PLOTS_DIR   <- file.path(REPO_ROOT, "plots")

dir.create(OUTPUTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)
