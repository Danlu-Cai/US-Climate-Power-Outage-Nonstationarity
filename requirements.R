
# =============================================================================
# Script: requirements.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Install the R packages required for Scripts 01–04.
# =============================================================================

required_packages <- c(
  # Data processing
  "data.table",
  "lubridate",
  "TTR",

  # Data manipulation and input/output
  "dplyr",
  "readr",
  "tidyr",
  "tibble",

  # Machine learning
  "xgboost",
  "ranger",

  # Visualization
  "ggplot2",
  "patchwork"
)

# Identify missing packages
missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

# Install missing packages
if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

# Verify required packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg)
  }
}

message("All required R packages are available.")
