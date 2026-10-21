# =============================================================================
# Script: 03_identify_conus404_temperature_extremes.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Identify county-level extreme heat and cold days from CONUS404 daily mean
# temperature (October 1979–September 2022), with daily, monthly, calendar-
# year, and quality-assurance outputs.
#
# Method (preserved from original Extremes_D_M_Y_saving.R):
# - County-specific, calendar-month 2nd, 5th, 95th and 98th percentiles,
#   calculated using all available days in the 1979-10 to 2022-09 baseline.
# - Candidate heat (cold) runs: >=3 consecutive days strictly above the 95th
#   (below the 5th) percentile. Within each candidate run, test moving-average
#   windows of length n, n-1, ..., 3 against the 98th (2nd) percentile.
#   Mark the first longest qualifying window, choosing the hottest/coldest
#   window if several qualify. Only that window is flagged.
# - IMPORTANT: Runs are evaluated WITHIN EACH CALENDAR MONTH, matching the
#   original analysis. Cross-month events are NOT joined. Changing this
#   convention may change the resulting extreme-day counts.
# - Missing temperature observations BREAK runs; they are not treated as zero
#   or removed before evaluating continuity.
# - A calendar year at the beginning/end of the record may be incomplete.
#
# Requirements: data.table, TTR
# Run from the repository root, or set the two environment variables below.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TTR)
})

# ---- CONFIGURATION -----------------------------------------------------------
INPUT_DIR <- Sys.getenv("CONUS404_DAILY_INPUT_DIR",
                        unset = "data/raw/conus404/daily_county")
OUTPUT_DIR <- Sys.getenv("CONUS404_EXTREMES_OUTPUT_DIR",
                         unset = "results/conus404_temperature_extremes")
SCALE <- "county"
START_MONTH <- "197910"
END_MONTH <- "202209"

# ---- HELPERS -----------------------------------------------------------------
month_days <- function(ym) {
  start <- as.Date(paste0(substr(ym, 1, 4), "-", substr(ym, 5, 6), "-01"))
  as.integer(format(seq(start, by = "month", length.out = 2L)[2L] - 1, "%d"))
}

read_month <- function(path, ym) {
  dt <- fread(path, colClasses = list(character = "GEOID"))
  if (!"GEOID" %in% names(dt)) stop("Missing GEOID: ", path)
  dt[, GEOID := trimws(GEOID)]
  if (anyNA(dt$GEOID) || any(!grepl("^[0-9]{1,5}$", dt$GEOID))) {
    stop("Invalid GEOID: ", path)
  }
  dt[, GEOID := sprintf("%05d", as.integer(GEOID))]
  if (any(dt$GEOID == "00000") || anyDuplicated(dt$GEOID)) {
    stop("Zero or duplicate GEOID: ", path)
  }
  days <- grep("^Day_", names(dt), value = TRUE)
  day_numbers <- suppressWarnings(as.integer(sub("^Day_", "", days)))
  if (!length(days) || anyNA(day_numbers) || anyDuplicated(day_numbers) ||
      !setequal(day_numbers, seq_len(month_days(ym)))) {
    stop("Invalid/incomplete Day_ columns for ", ym, " in ", path)
  }
  long <- melt(dt, id.vars = "GEOID", measure.vars = days,
               variable.name = "Day_ID", value.name = "T2_C",
               variable.factor = FALSE)
  long[, Date := as.Date(sprintf("%s-%s-%02d", substr(ym, 1, 4),
                                  substr(ym, 5, 6),
                                  as.integer(sub("^Day_", "", Day_ID))))]
  # Do not silently turn malformed nonempty strings into missing values.
  original <- trimws(as.character(long$T2_C))
  temp <- suppressWarnings(as.numeric(original))
  bad <- !is.na(original) & nzchar(original) & is.na(temp)
  if (any(bad) || any(!is.na(temp) & !is.finite(temp))) {
    stop("Invalid nonnumeric/nonfinite temperatures in ", path)
  }
  long[, T2_C := temp]
  long[, Month_Ref := ym]
  setorder(long, GEOID, Date)
  long[, .(GEOID, Date, Month_Ref, T2_C)]
}

# This implements the original two-tier selection within one consecutive run.
# A single best qualifying window is flagged, not necessarily the whole run.
mark_best_window <- function(values, strict_threshold, type = c("heat", "cold")) {
  type <- match.arg(type)
  n <- length(values)
  result <- integer(n)
  if (n < 3L || anyNA(values) || is.na(strict_threshold)) return(result)
  for (k in seq.int(n, 3L, by = -1L)) {
    means <- as.numeric(TTR::SMA(values, n = k))
    good <- if (type == "heat") which(means >= strict_threshold) else
      which(means <= strict_threshold)
    if (length(good)) {
      endpoint <- if (type == "heat") which.max(means) else which.min(means)
      result[seq.int(endpoint - k + 1L, endpoint)] <- 1L
      break
    }
  }
  result
}

mark_extremes <- function(temps, candidate_threshold, strict_threshold,
                          type = c("heat", "cold")) {
  type <- match.arg(type)
  output <- rep(NA_integer_, length(temps))
  output[!is.na(temps)] <- 0L
  if (is.na(candidate_threshold) || is.na(strict_threshold)) return(output)
  candidate <- !is.na(temps) & if (type == "heat")
    temps > candidate_threshold else temps < candidate_threshold
  runs <- rle(candidate)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (j in which(runs$values & runs$lengths >= 3L)) {
    ix <- seq.int(starts[j], ends[j])
    output[ix] <- mark_best_window(temps[ix], strict_threshold, type)
  }
  output
}

# ---- FIND AND VALIDATE INPUTS ------------------------------------------------
if (!dir.exists(INPUT_DIR)) stop("Input directory does not exist: ", INPUT_DIR)
for (subdir in c("Daily", "Monthly", "Yearly", "QA")) {
  dir.create(file.path(OUTPUT_DIR, subdir), recursive = TRUE,
             showWarnings = FALSE)
}
paths <- list.files(INPUT_DIR, pattern = paste0("_", SCALE, "_merged\\.csv$"),
                    full.names = TRUE, recursive = TRUE)
extract_month <- function(path) {
  match <- regmatches(basename(path),
                     regexpr("(?<![0-9])[0-9]{6}(?![0-9])",
                             basename(path), perl = TRUE))
  if (length(match) == 0L) NA_character_ else match
}
files <- data.table(path = paths,
                    Month_Ref = vapply(paths, extract_month, character(1)))
files <- files[!is.na(Month_Ref) & Month_Ref >= START_MONTH &
                 Month_Ref <= END_MONTH]
setorder(files, Month_Ref)
expected <- format(seq(as.Date("1979-10-01"), as.Date("2022-09-01"),
                       by = "month"), "%Y%m")
if (anyDuplicated(files$Month_Ref)) stop("Duplicate monthly input files detected")
missing <- setdiff(expected, files$Month_Ref)
if (length(missing)) stop("Missing input months: ", paste(missing, collapse = ", "))
message("Input months: ", nrow(files))

# ---- PASS 1: HISTORICAL COUNTY x CALENDAR-MONTH THRESHOLDS --------------------
# Thresholds use all available daily values, not a moving/annual baseline.
baseline_parts <- vector("list", nrow(files))
for (i in seq_len(nrow(files))) {
  ym <- files$Month_Ref[i]
  if (i %% 24L == 1L || i == nrow(files)) message("Baseline: ", i, "/", nrow(files))
  d <- read_month(files$path[i], ym)
  d[, month_val := substr(Month_Ref, 5, 6)]
  baseline_parts[[i]] <- d[!is.na(T2_C), .(GEOID, month_val, T2_C)]
  rm(d)
}
baseline <- rbindlist(baseline_parts, use.names = TRUE)
rm(baseline_parts)
thresholds <- baseline[, .(
  q02 = as.numeric(quantile(T2_C, 0.02, na.rm = TRUE)),
  q05 = as.numeric(quantile(T2_C, 0.05, na.rm = TRUE)),
  q95 = as.numeric(quantile(T2_C, 0.95, na.rm = TRUE)),
  q98 = as.numeric(quantile(T2_C, 0.98, na.rm = TRUE)),
  baseline_valid_days = .N
), by = .(GEOID, month_val)]
rm(baseline)
invisible(gc())
fwrite(thresholds, file.path(OUTPUT_DIR, "QA", "county_month_baseline_thresholds.csv"))

# ---- PASS 2: DAILY EXTREMES, MONTHLY COUNTS, AND QA ---------------------------
monthly_parts <- vector("list", nrow(files))
qa_parts <- vector("list", nrow(files))
for (i in seq_len(nrow(files))) {
  ym <- files$Month_Ref[i]
  message("Extremes: ", ym, " (", i, "/", nrow(files), ")")
  d <- read_month(files$path[i], ym)
  d[, month_val := substr(Month_Ref, 5, 6)]
  d <- merge(d, thresholds[, .(GEOID, month_val, q02, q05, q95, q98)],
             by = c("GEOID", "month_val"), all.x = TRUE, sort = FALSE)
  setorder(d, GEOID, Date)
  if (anyNA(d$q02) || anyNA(d$q05) || anyNA(d$q95) || anyNA(d$q98)) {
    stop("Missing baseline thresholds for one or more counties in ", ym)
  }
  # Grouping by GEOID within a single month preserves the original month-
  # boundary behavior; NA temperatures remain in place and break sequences.
  d[, `:=`(
    heat_ext = mark_extremes(T2_C, q95[1], q98[1], "heat"),
    cold_ext = mark_extremes(T2_C, q05[1], q02[1], "cold")
  ), by = GEOID]

  fwrite(d[, .(GEOID, Date, Month_Ref, T2_C, heat_ext, cold_ext)],
         file.path(OUTPUT_DIR, "Daily", paste0("Daily_Extremes_", ym, ".csv")))
  monthly_parts[[i]] <- d[, .(
    heat_days = sum(heat_ext, na.rm = TRUE),
    cold_days = sum(cold_ext, na.rm = TRUE),
    valid_days = sum(!is.na(T2_C)),
    missing_days = sum(is.na(T2_C)),
    expected_days = .N
  ), by = .(GEOID, Month_Ref)]
  qa_parts[[i]] <- data.table(Month_Ref = ym, counties = uniqueN(d$GEOID),
                              total_temperature_cells = nrow(d),
                              missing_temperature_cells = sum(is.na(d$T2_C)))
  rm(d)
}

monthly <- rbindlist(monthly_parts)
monthly[, coverage_pct := round(100 * valid_days / expected_days, 2)]
setorder(monthly, GEOID, Month_Ref)
fwrite(monthly, file.path(OUTPUT_DIR, "Monthly",
                          "CONUS404_Monthly_Extremes_1979_2022.csv"))

# ---- CALENDAR-YEAR SUMMARY ---------------------------------------------------
monthly[, Year := as.integer(substr(Month_Ref, 1, 4))]
yearly <- monthly[, .(
  annual_heat_days = sum(heat_days),
  annual_cold_days = sum(cold_days),
  months_covered = .N,
  valid_days = sum(valid_days),
  missing_days = sum(missing_days),
  expected_days = sum(expected_days)
), by = .(GEOID, Year)]
yearly[, `:=`(
  coverage_pct = round(100 * valid_days / expected_days, 2),
  complete_calendar_year = months_covered == 12L & missing_days == 0L
)]
setorder(yearly, GEOID, Year)
fwrite(yearly, file.path(OUTPUT_DIR, "Yearly",
                         "CONUS404_Yearly_Extremes_1979_2022.csv"))
fwrite(rbindlist(qa_parts), file.path(OUTPUT_DIR, "QA",
                                      "monthly_input_missingness.csv"))
message("Completed: daily, monthly, yearly, baseline thresholds, and QA outputs.")
