
# R package dependencies for the climate–power outage analysis

required_packages <- c(
  "tidyverse",
  "lubridate",
  "conflicted"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

message("Required R packages are installed.")
