# ======================================================================
# main.R — Run the full analysis end-to-end.
#
# Usage:
#   1. Place the CAMELS-GB v2 files in data/ (see data/README.md).
#   2. Install the required packages once: source("install_packages.R")
#   3. source("main.R")
#
# The whole pipeline takes a few minutes depending on machine and the
# number of catchments processed.
# ======================================================================

cat("=== Bristol CAMELS-GB flooding analysis ===\n\n")

cat("[1/3] Data preparation...\n")
source(here::here("R", "01_data_preparation.R"))

cat("\n[2/3] Statistical analysis...\n")
source(here::here("R", "02_analysis.R"))

cat("\n[3/3] Visualisation...\n")
source(here::here("R", "03_visualisation.R"))

cat("\nDone. Figures are in outputs/.\n")
