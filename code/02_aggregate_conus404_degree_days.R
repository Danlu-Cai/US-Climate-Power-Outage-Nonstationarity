# =============================================================================
# CONUS404 county-level degree-day aggregation (October 1979–September 2022)
# Original analysis: Danlu Cai
# Reproducible adaptation: configurable paths, input checks, coverage diagnostics
# =============================================================================

suppressPackageStartupMessages(library(data.table))


INPUT_DIR <- Sys.getenv(
  "CONUS404_DAILY_INPUT_DIR",
  unset = "data/raw/conus404/daily_county"
)

OUTPUT_DIR <- Sys.getenv(
  "CONUS404_DD_OUTPUT_DIR",
  unset = "results/conus404_degree_days"
)

SCALE <- "county"
BASE_TEMP_C <- 18.3
START_MONTH <- "197910"
END_MONTH <- "202209"

# Each input CSV is assumed to contain one row per county, GEOID, and columns
# Day_1 ... Day_31 (or Day_01 ... Day_31) with DAILY MEAN 2-m temperature in C.
# Check these assumptions against the upstream CONUS404 extraction workflow.
if (!dir.exists(INPUT_DIR)) stop("Input directory does not exist: ", INPUT_DIR)
for (folder in c("Daily", "Monthly", "Yearly", "QA")) {
  dir.create(file.path(OUTPUT_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}

all_files <- list.files(INPUT_DIR, pattern = paste0("_", SCALE, "_merged\\.csv$"), full.names = TRUE)
file_month <- regmatches(basename(all_files), regexpr("(?<![0-9])[0-9]{6}(?![0-9])", basename(all_files), perl = TRUE))
file_month[lengths(file_month) == 0L] <- NA_character_
files <- data.table(path = all_files, Month_Ref = as.character(file_month))
files <- files[!is.na(Month_Ref) & Month_Ref >= START_MONTH & Month_Ref <= END_MONTH]
setorder(files, Month_Ref)

expected_months <- format(seq(as.Date("1979-10-01"), as.Date("2022-09-01"), by = "month"), "%Y%m")
if (anyDuplicated(files$Month_Ref)) stop("Multiple source files for the same month: ", paste(unique(files$Month_Ref[duplicated(files$Month_Ref)]), collapse = ", "))
missing_months <- setdiff(expected_months, files$Month_Ref)
if (length(missing_months)) stop("Missing monthly input files: ", paste(missing_months, collapse = ", "))
message("Processing ", nrow(files), " monthly files...")

monthly_results <- vector("list", nrow(files))
qa_results <- vector("list", nrow(files))

for (i in seq_len(nrow(files))) {
  ym <- files$Month_Ref[i]
  message("Processing ", ym)
  dt <- fread(files$path[i], colClasses = list(character = "GEOID"))
  if (!"GEOID" %in% names(dt)) stop("Missing GEOID: ", files$path[i])
  dt[, GEOID := trimws(as.character(GEOID))]
  if (anyNA(dt$GEOID) || any(dt$GEOID == "")) stop("Missing GEOID in ", ym)
  if (anyDuplicated(dt$GEOID)) stop("Duplicate GEOID in ", ym)

  day_cols <- grep("^Day_", names(dt), value = TRUE)
  if (!length(day_cols)) stop("No Day_ columns in ", ym)
  # Fail explicitly for unknown column conventions instead of assigning wrong dates.
  if (any(!grepl("^Day_0?[1-9]$|^Day_[12][0-9]$|^Day_3[01]$", day_cols))) {
    stop("Unexpected Day_ column naming in ", ym, ": ", paste(day_cols[!grepl("^Day_0?[1-9]$|^Day_[12][0-9]$|^Day_3[01]$", day_cols)], collapse = ", "))
  }
  day_numbers <- as.integer(sub("^Day_", "", day_cols))
  days_in_month <- as.integer(format(seq(as.Date(paste0(substr(ym, 1, 4), "-", substr(ym, 5, 6), "-01")), by = "month", length.out = 2)[2] - 1, "%d"))
  if (anyDuplicated(day_numbers) || !setequal(day_numbers, seq_len(days_in_month))) {
    stop("Expected exactly Day_1 ... Day_", days_in_month, " in ", ym, "; found ", paste(day_cols, collapse = ", "))
  }

  daily <- melt(dt, id.vars = "GEOID", measure.vars = day_cols,
                variable.name = "Day_ID", value.name = "T2_C", variable.factor = FALSE)
  daily[, day_number := as.integer(sub("^Day_", "", Day_ID))]
  daily[, Date := as.Date(sprintf("%s-%s-%02d", substr(ym, 1, 4), substr(ym, 5, 6), day_number))]
  daily[, T2_C := suppressWarnings(as.numeric(T2_C))]
  if (any(!is.na(daily$T2_C) & !is.finite(daily$T2_C))) stop("Non-finite temperature in ", ym)
  daily[, `:=`(
    CDD = fifelse(is.na(T2_C), NA_real_, pmax(T2_C - BASE_TEMP_C, 0)),
    HDD = fifelse(is.na(T2_C), NA_real_, pmax(BASE_TEMP_C - T2_C, 0)),
    Month_Ref = ym
  )]
  setcolorder(daily, c("GEOID", "Date", "Month_Ref", "Day_ID", "T2_C", "CDD", "HDD"))
  fwrite(daily[, .(GEOID, Date, Month_Ref, Day_ID, T2_C, CDD, HDD)],
         file.path(OUTPUT_DIR, "Daily", paste0("Daily_DD_", ym, ".csv")))

  monthly_results[[i]] <- daily[, .(
    mean_T2C = if (all(is.na(T2_C))) NA_real_ else mean(T2_C, na.rm = TRUE),
    monthly_CDD = if (all(is.na(CDD))) NA_real_ else sum(CDD, na.rm = TRUE),
    monthly_HDD = if (all(is.na(HDD))) NA_real_ else sum(HDD, na.rm = TRUE),
    valid_days = sum(!is.na(T2_C)),
    expected_days = days_in_month,
    missing_days = sum(is.na(T2_C)),
    sum_valid_T2C = sum(T2_C, na.rm = TRUE)
  ), by = .(GEOID, Month_Ref)]
  qa_results[[i]] <- data.table(Month_Ref = ym, counties = nrow(dt),
                                missing_temperature_cells = sum(is.na(daily$T2_C)),
                                total_temperature_cells = nrow(daily))
  rm(dt, daily)
}

monthly <- rbindlist(monthly_results)
monthly[, coverage_pct := round(100 * valid_days / expected_days, 2)]
setcolorder(monthly, c("GEOID", "Month_Ref", "mean_T2C", "monthly_CDD", "monthly_HDD",
                       "valid_days", "expected_days", "missing_days", "coverage_pct", "sum_valid_T2C"))
fwrite(monthly[, !"sum_valid_T2C"], file.path(OUTPUT_DIR, "Monthly", "CONUS404_Monthly_Degree_Days_1979_2022.csv"))

monthly[, Year := substr(Month_Ref, 1, 4)]
yearly <- monthly[, .(
  annual_CDD = if (all(is.na(monthly_CDD))) NA_real_ else sum(monthly_CDD, na.rm = TRUE),
  annual_HDD = if (all(is.na(monthly_HDD))) NA_real_ else sum(monthly_HDD, na.rm = TRUE),
  avg_T2C = if (sum(valid_days) == 0L) NA_real_ else sum(sum_valid_T2C) / sum(valid_days),
  months_covered = .N,
  valid_days = sum(valid_days),
  expected_days = sum(expected_days),
  missing_days = sum(missing_days)
), by = .(GEOID, Year)]
yearly[, `:=`(
  coverage_pct = round(100 * valid_days / expected_days, 2),
  complete_calendar_year = (months_covered == 12L & missing_days == 0L)
)]
fwrite(yearly, file.path(OUTPUT_DIR, "Yearly", "CONUS404_Yearly_Degree_Days_1979_2022.csv"))
fwrite(rbindlist(qa_results), file.path(OUTPUT_DIR, "QA", "monthly_input_missingness.csv"))
message("Completed. Calendar years 1979 and 2022 are partial by design; check coverage fields before using annual totals.")
