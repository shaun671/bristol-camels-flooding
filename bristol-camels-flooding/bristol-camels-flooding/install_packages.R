# Install all packages required by the analysis. Run once.

required_packages <- c(
  # Core / tidyverse
  "here", "readr", "dplyr", "tidyr", "purrr", "stringr", "lubridate",
  "data.table", "tools",
  # Spatial
  "sf", "tmap", "gstat",
  # Modelling and diagnostics
  "car", "lmtest", "DescTools", "performance", "relaimpo", "proxy",
  # Visualisation
  "ggplot2", "patchwork", "viridis", "ggcorrplot",
  # Excel I/O (used by some auxiliary CAMELS files)
  "readxl"
)

new_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(new_packages) > 0) {
  message("Installing: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}
