# =============================================================================
# 01_process_eaglei_high_impact_outages.R
# Climate–Power Outage Relationships Across the United States
#
# Purpose: Aggregate EAGLE-I 15-minute county outage snapshots (2015–2024),
#          identify qualifying continuous county-level outage events, and
#          create hourly, daily, monthly and annual summaries.
#
# Definition: A high-impact county-day has (1) observed daily COH > 5,000
#             AND (2) overlaps an event with > 8 consecutive complete clock
#             hours of positive customer-outage observations.
#
# IMPORTANT ASSUMPTIONS (verify against original data and Methods):
#   - One observation per county per 15-minute slot; duplicates are rejected.
#   - A qualifying outage hour requires all four quarter-hour slots present
#     AND customers_out > 0 in all four slots. Missing slots break continuity.
#   - COH = sum(customers_out * 0.25) across observed 15-minute slots.
#   - Daily outage hours = number of positive 15-minute observations * 0.25;
#     this is county-level positive-snapshot time, NOT customer interruption
#     duration or a measure of service continuity for individual customers.
#   - Each day overlapping a qualifying event must separately exceed 5,000 COH.
#   - Timestamps must be expressed in one consistent timezone. Default UTC;
#     change TIMEZONE if your source timestamps use a different convention.
#   - All annual files must be available to detect events crossing year ends.
#
# Dependencies: data.table, lubridate
# Install once: install.packages(c('data.table', 'lubridate'))
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
})

# ---- USER CONFIGURATION ------------------------------------------------------
INPUT_DIR <- Sys.getenv('EAGLEI_INPUT_DIR', unset = 'data/raw/eagle-i')
OUTPUT_DIR <- Sys.getenv('EAGLEI_OUTPUT_DIR', unset = 'results/eaglei_high_impact')
YEARS <- 2015:2024
TIMEZONE <- 'UTC'       # Verify timezone of run_start_time in source files.
DURATION_THRESHOLD_HOURS <- 8L  # Strictly greater than 8 hours.
DAILY_COH_THRESHOLD <- 5000     # Strictly greater than 5,000 COH.

# ---- HELPERS -----------------------------------------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

empty_daily <- function() {
  data.table(fips_code = character(), county = character(), state = character(),
             date_stamp = as.Date(character()), daily_coh = numeric(),
             total_outage_hours = numeric(), max_continuous_streak = integer(),
             observed_quarter_hours = integer(), missing_quarter_hours = integer())
}

# Accept the formats used in the original script. Timezone is explicit.
parse_times <- function(x) {
  suppressWarnings(parse_date_time(
    as.character(x),
    orders = c('mdy HMS', 'mdy HM', 'mdy I p', 'mdy Ip',
               'ymd HMS', 'ymd HM', 'ymd'),
    truncated = 0, tz = TIMEZONE, quiet = TRUE
  ))
}

# ---- INPUT VALIDATION AND HOURLY AGGREGATION ---------------------------------
ensure_dir(OUTPUT_DIR)
for (folder in c('Hourly_Summary', 'Daily_Summary',
                 'Monthly_Summary', 'Yearly_Summary')) {
  ensure_dir(file.path(OUTPUT_DIR, folder))
}

input_files <- file.path(INPUT_DIR, sprintf('eaglei_outages_%d.csv', YEARS))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop('Missing annual input files; cannot reliably detect year-boundary events:\n',
       paste(missing_files, collapse = '\n'))
}

hourly_files <- character(length(YEARS))
for (i in seq_along(YEARS)) {
  yr <- YEARS[i]
  message('Preparing hourly summaries: ', yr)
  x <- fread(input_files[i], colClasses = list(character = 'fips_code'))
  if ('sum' %in% names(x) && !'customers_out' %in% names(x)) {
    setnames(x, 'sum', 'customers_out')
  }
  required <- c('fips_code', 'county', 'state', 'run_start_time', 'customers_out')
  absent <- setdiff(required, names(x))
  if (length(absent)) stop('Missing columns in ', input_files[i], ': ', paste(absent, collapse=', '))

  x[, `:=`(fips_code = as.character(fips_code),
           customers_out = as.numeric(customers_out),
           datetime = parse_times(run_start_time))]
  invalid <- x[is.na(datetime) | is.na(customers_out) |
                 !is.finite(customers_out) | customers_out < 0 |
                 is.na(fips_code) | fips_code == '']
  if (nrow(invalid)) {
    stop('Invalid timestamps, FIPS, or outage counts in ', input_files[i],
         ' (', nrow(invalid), ' rows). Fix source records before continuing.')
  }
  # Require exact 15-minute boundaries. Rounding could silently alter COH.
  x[, slot := floor_date(datetime, '15 minutes')]
  if (any(as.numeric(difftime(x$datetime, x$slot, units = 'secs')) != 0)) {
    stop('Non-quarter-hour timestamps found in ', input_files[i],
         '. Check timestamp conventions before aggregation.')
  }
  # Check within-year duplicates; cross-year duplicates are checked below.
  if (anyDuplicated(x, by = c('fips_code', 'slot'))) {
    stop('Duplicate county/15-minute timestamps in ', input_files[i],
         '. Resolve duplicates instead of double-counting COH.')
  }
  x[, hour_stamp := floor_date(slot, 'hour')]
  hourly <- x[, .(
    hourly_coh = sum(customers_out) * 0.25,
    avg_customers_out = mean(customers_out),
    observed_quarter_hours = .N,
    positive_quarter_hours = sum(customers_out > 0),
    # Count as a complete positive clock hour only if all four slots are
    # observed and all four have positive outage counts.
    complete_positive_hour = (.N == 4L && all(customers_out > 0))
  ), by = .(fips_code, county, state, hour_stamp)]
  if (any(hourly$observed_quarter_hours > 4L)) {
    stop('More than four observations in an hourly county group for ', yr)
  }
  setorder(hourly, fips_code, hour_stamp)
  hourly_files[i] <- file.path(OUTPUT_DIR, 'Hourly_Summary',
                               sprintf('eaglei_outages_hourly_%d.csv', yr))
  fwrite(hourly, hourly_files[i])
  rm(x, hourly)
  invisible(gc())
}

# ---- CONTINUOUS EVENTS ACROSS MIDNIGHT AND YEAR BOUNDARIES -------------------
# This loads hourly summaries for all years. For very large datasets, use a
# county-partitioned version to reduce peak memory usage.
message('Combining hourly summaries to identify cross-year events...')
h <- rbindlist(lapply(hourly_files, fread), use.names = TRUE)
h[, hour_stamp := as.POSIXct(hour_stamp, tz = TIMEZONE)]
if (anyNA(h$hour_stamp)) stop('Could not parse hourly summary timestamps.')
if (anyDuplicated(h, by = c('fips_code', 'hour_stamp'))) {
  stop('Duplicate county/hour across annual files. Check overlapping input years.')
}
setorder(h, fips_code, hour_stamp)
h[, complete_positive_hour := as.logical(complete_positive_hour)]

# An event starts at a complete positive hour when the preceding hour is not
# complete-positive OR is missing OR is not exactly one clock hour earlier.
h[, previous_time := shift(hour_stamp), by = fips_code]
h[, previous_positive := shift(complete_positive_hour), by = fips_code]
h[, starts_event := complete_positive_hour & (
  is.na(previous_time) | is.na(previous_positive) | !previous_positive |
    as.numeric(difftime(hour_stamp, previous_time, units = 'hours')) != 1
)]
h[, event_id := cumsum(starts_event), by = fips_code]

# Event IDs are county-specific; count qualifying hours only within positive runs.
events <- h[complete_positive_hour == TRUE,
            .(event_duration_hrs = .N,
              event_start = min(hour_stamp), event_end = max(hour_stamp)),
            by = .(fips_code, event_id)]
events[, qualifies_duration := event_duration_hrs > DURATION_THRESHOLD_HOURS]

h <- merge(h, events[, .(fips_code, event_id, event_duration_hrs,
                         qualifies_duration)],
           by = c('fips_code', 'event_id'), all.x = TRUE, sort = FALSE)
h[is.na(qualifies_duration), qualifies_duration := FALSE]
h[is.na(event_duration_hrs), event_duration_hrs := 0L]
h[, date_stamp := as.Date(hour_stamp, tz = TIMEZONE)]

# ---- DAILY OUTAGE METRICS AND HIGH-IMPACT FILTER -----------------------------
# Hourly COH sums observed quarter-hour snapshots. This remains an observed
# COH estimate when slots are missing; missing data are not imputed.
daily_all <- h[, .(
  daily_coh = sum(hourly_coh),
  total_outage_hours = sum(positive_quarter_hours) * 0.25,
  observed_quarter_hours = sum(observed_quarter_hours),
  missing_quarter_hours = 96L - sum(observed_quarter_hours),
  # Only event hours that overlap THIS calendar day contribute to qualification.
  overlaps_qualifying_event = any(complete_positive_hour & qualifies_duration),
  max_continuous_streak = if (any(complete_positive_hour & qualifies_duration))
    as.integer(max(event_duration_hrs[complete_positive_hour & qualifies_duration]))
    else 0L
), by = .(fips_code, county, state, date_stamp)]

if (any(daily_all$missing_quarter_hours < 0L)) {
  stop('A county-day has more than 96 observed quarter-hour snapshots.')
}

daily <- daily_all[overlaps_qualifying_event == TRUE &
                   daily_coh > DAILY_COH_THRESHOLD]
daily[, overlaps_qualifying_event := NULL]

# ---- WRITE ANNUAL FILES (INCLUDING EMPTY OUTPUTS) ----------------------------
for (yr in YEARS) {
  message('Writing high-impact summaries: ', yr)
  dy <- daily[year(date_stamp) == yr]
  daily_path <- file.path(OUTPUT_DIR, 'Daily_Summary',
                          sprintf('eaglei_outages_daily_%d.csv', yr))
  fwrite(dy, daily_path)

  if (nrow(dy)) {
    dy[, month_stamp := as.Date(floor_date(date_stamp, 'month'))]
    mo <- dy[, .(monthly_coh = sum(daily_coh),
                 monthly_outage_hours = sum(total_outage_hours)),
             by = .(fips_code, county, state, month_stamp)]
    an <- mo[, .(annual_coh = sum(monthly_coh),
                 annual_outage_hours = sum(monthly_outage_hours)),
             by = .(fips_code, county, state)]
    an[, year := yr]
  } else {
    mo <- data.table(fips_code = character(), county = character(),
                     state = character(), month_stamp = as.Date(character()),
                     monthly_coh = numeric(), monthly_outage_hours = numeric())
    an <- data.table(fips_code = character(), county = character(),
                     state = character(), annual_coh = numeric(),
                     annual_outage_hours = numeric(), year = integer())
  }
  fwrite(mo, file.path(OUTPUT_DIR, 'Monthly_Summary',
                        sprintf('eaglei_outages_monthly_%d.csv', yr)))
  fwrite(an, file.path(OUTPUT_DIR, 'Yearly_Summary',
                        sprintf('eaglei_outages_yearly_%d.csv', yr)))
  message('  Qualifying county-days: ', nrow(dy))
}

message('Completed. Check data coverage, timezone, and comparison with the original outputs.')
