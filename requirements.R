
# =============================================================================
# Script: requirements.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Install the R packages required to reproduce the data processing
# workflows for EAGLE-I power outages and CONUS404 climate data.
#
# Supported scripts:
# 01_process_eaglei_high_impact_outages.R
# 02_aggregate_conus404_degree_days.R
# 03_identify_conus404_temperature_extremes.R
# =============================================================================

# Required R packages
required_packages <- c(
  "data.table",
  "lubridate",
  "TTR"
)

# Identify packages not yet installed
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

# Confirm package availability
invisible(lapply(
  required_packages,
  function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required package unavailable: ", pkg)
    }
  }
))

message("All required R packages are available.")
