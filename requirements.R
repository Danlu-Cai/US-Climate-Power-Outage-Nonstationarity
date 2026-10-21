
# R package dependencies
# EAGLE-I high-impact outage processing

required_packages <- c(
  "data.table",
  "lubridate"
)

missing_packages <- setdiff(
  required_packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
