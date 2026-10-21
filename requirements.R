
# =============================================================================
# Script: requirements.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Install R packages required for data processing, machine learning,
# and figure generation.
# =============================================================================

required_packages <- c(
  "data.table",
  "lubridate",
  "TTR",
  "dplyr",
  "readr",
  "tidyr",
  "tibble",
  "purrr",
  "stringr",
  "ggplot2",
  "patchwork",
  "cowplot",
  "ggrepel",
  "scales",
  "xgboost",
  "ranger",
  "iml",
  "Matrix",
  "sf",
  "usmap",
  "lmtest",
  "sandwich",
  "conflicted",
  "tidyverse"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package unavailable: ", pkg)
  }
}

message("All required R packages are available.")
