# =============================================================================
# Script: Figure_1_outage_trends.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Generate individual Figure 1 panels describing county-level outage-hour
# trends, regional heterogeneity, national annual outage hours and COH,
# and monthly/seasonal comparisons from EAGLE-I daily summaries (2015–2024).
#
# Important:
# County trends are expressed in outage hours per year. The original script
# multiplied year-based regression slopes by 12; that unit error is corrected.
# Review missing-data treatment before interpreting trends: this script
# retains the original zero-filling rules for missing records/months.
# =============================================================================

library(tidyverse)
library(lubridate)
library(scales)
library(usmap)       # County & State spatial data
library(sf)
library(grid)        # Graphical object manipulations
library(lmtest)      # Coeftest for HAC robust regression
library(sandwich)    # HAC covariance matrix estimation
library(cowplot)
library(patchwork)

# ==============================================================================
# 0. GLOBAL SETUP, PATHS & TYPOGRAPHY
# ==============================================================================

daily_summary_dir <- Sys.getenv("EAGLEI_DAILY_SUMMARY_DIR", unset = "results/eaglei_high_impact/Daily_Summary")
output_dir <- Sys.getenv("FIGURE1_OUTPUT_DIR", unset = "results/figures/figure1")


if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

conus_fips <- c(
  "01", "04", "05", "06", "08", "09", "10", "11", "12", "13",
  "16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
  "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
  "36", "37", "38", "39", "40", "41", "42", "44", "45", "46",
  "47", "48", "49", "50", "51", "53", "54", "55", "56"
)

# Date Bounds
start_yr <- 2015
end_yr   <- 2024

# --- NATIVE PANEL DIMENSIONS ---
PANEL_WIDTH     <- 10.0
PANEL_HEIGHT    <- 9.5

# Typography & Graphic Settings
FONT_AXIS_TITLE   <- 24   # Axis Titles & Section Headers
FONT_LEG_TITLE    <- 24   # Legend Headers
FONT_TICK_TEXT    <- 22   # Axis Ticks & Colorbar Ticks
FONT_AXIS_TICK    <- 22   # Standardized Tick Font Size
FONT_LEG_TEXT     <- 22   # Legend Items & Text Sub-labels
GEOM_TEXT_SIZE    <- 6.5  # Geom Text Size

# Palette Definitions
COLOR_HOURS_BAR   <- "#EE5B43" # Salmon fill for outage hours
COLOR_TEXT_LEFT   <- "#D9533F" # Matching red for left axis text/ticks
COLOR_COH         <- "#0A284B" # Dark navy for COH line/points
COLOR_CHANGE_TEXT <- "#D9533F" # Color for percentage change annotations

# Standardized Minimal Theme
theme_nature_super <- function() {
  theme_minimal(base_size = 20) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.6),
      axis.line.x        = element_line(color = "black", linewidth = 1.2),
      axis.ticks.x       = element_line(color = "black", linewidth = 1.2),
      axis.text.x        = element_text(face = "bold", size = FONT_AXIS_TICK, color = "black"),
      axis.line.y        = element_line(color = "black", linewidth = 1.2),
      axis.ticks.y       = element_line(color = "black", linewidth = 1.2),
      axis.title.y       = element_text(face = "bold", size = FONT_AXIS_TITLE, color = "black", margin = margin(r = 12)),
      axis.text.y        = element_text(face = "bold", size = FONT_AXIS_TICK, color = "black")
    )
}

# ==============================================================================
# 1. DATA PREPROCESSING & REGIONAL MAPPINGS
# ==============================================================================

message("[INFO] Loading daily outage data...")
file_list <- list.files(daily_summary_dir, pattern = "eaglei_outages_daily_\\d{4}\\.csv", full.names = TRUE)

read_and_standardize <- function(file_path) {
  df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  names(df) <- tolower(trimws(names(df)))
  return(df)
}

df_raw <- file_list %>% map_df(read_and_standardize)

coh_col_name <- names(df_raw)[grep("coh", names(df_raw))[1]]
if (is.na(coh_col_name)) coh_col_name <- "daily_coh"

outage_col_name <- names(df_raw)[grep("outage.*hour|outage.*hr|tot.*outage|outage", names(df_raw))[1]]
if (is.na(outage_col_name)) outage_col_name <- "total_outage_hours"

df_all <- df_raw %>%
  mutate(
    fips_std          = str_pad(as.character(fips_code), width = 5, pad = "0"),
    state_fips        = substr(fips_std, 1, 2),
    coh_val           = suppressWarnings(readr::parse_number(as.character(.data[[coh_col_name]]))),
    coh_val           = ifelse(is.na(coh_val), 0, coh_val),
    outage_hours_val = suppressWarnings(readr::parse_number(as.character(.data[[outage_col_name]]))),
    outage_hours_val = ifelse(is.na(outage_hours_val), 0, outage_hours_val),
    date_clean        = ymd(date_stamp),
    year              = year(date_clean),
    month_num         = month(date_clean)
  ) %>%
  filter(state_fips %in% conus_fips, year >= start_yr & year <= end_yr)

counties_sf <- usmap::us_map(regions = "counties") %>%
  filter(substr(fips, 1, 2) %in% conus_fips)

states_sf      <- usmap::us_map(regions = "states")
all_conus_fips <- unique(counties_sf$fips)

ncei_colors <- c(
  "Northeast"        = "#008B8B", "Northwest"        = "#2E8B57", "Ohio Valley"      = "#A16207",
  "Northern Rockies" = "#8B8B8B", "Upper Midwest"    = "#D97706", "Southeast"        = "#F97316",
  "Southwest"        = "#E65100", "West"             = "#1D4ED8", "South"            = "#EF4444"
)

state_to_ncei <- c(
  "CT"="Northeast", "DE"="Northeast", "ME"="Northeast", "MD"="Northeast", "MA"="Northeast", "NH"="Northeast", "NJ"="Northeast", "NY"="Northeast", "PA"="Northeast", "RI"="Northeast", "VT"="Northeast",
  "ID"="Northwest", "OR"="Northwest", "WA"="Northwest",
  "IL"="Ohio Valley", "IN"="Ohio Valley", "KY"="Ohio Valley", "MO"="Ohio Valley", "OH"="Ohio Valley", "TN"="Ohio Valley", "WV"="Ohio Valley",
  "MT"="Northern Rockies", "NE"="Northern Rockies", "ND"="Northern Rockies", "SD"="Northern Rockies", "WY"="Northern Rockies",
  "IA"="Upper Midwest", "MI"="Upper Midwest", "MN"="Upper Midwest", "WI"="Upper Midwest",
  "AL"="Southeast", "FL"="Southeast", "GA"="Southeast", "NC"="Southeast", "SC"="Southeast", "VA"="Southeast", 
  "AZ"="Southwest", "CO"="Southwest", "NM"="Southwest", "UT"="Southwest",
  "CA"="West", "NV"="West",
  "AR"="South", "KS"="South", "LA"="South", "MS"="South", "OK"="South", "TX"="South"
)

fips_to_state_map <- c(
  "01"="AL", "02"="AK", "04"="AZ", "05"="AR", "06"="CA", "08"="CO", "09"="CT",
  "10"="DE", "11"="DC", "12"="FL", "13"="GA", "15"="HI", "16"="ID", "17"="IL",
  "18"="IN", "19"="IA", "20"="KS", "21"="KY", "22"="LA", "23"="ME", "24"="MD",
  "25"="MA", "26"="MI", "27"="MN", "28"="MS", "29"="MO", "30"="MT", "31"="NE",
  "32"="NV", "33"="NH", "34"="NJ", "35"="NM", "36"="NY", "37"="NC", "38"="ND",
  "39"="OH", "40"="OK", "41"="OR", "42"="PA", "44"="RI", "45"="SC", "46"="SD",
  "47"="TN", "48"="TX", "49"="UT", "50"="VT", "51"="VA", "53"="WA", "54"="WV",
  "55"="WI", "56"="WY"
)

ncei_regions_sf <- states_sf %>%
  mutate(region = state_to_ncei[abbr]) %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(.groups = "drop")

ncei_items_3x3_fixed <- data.frame(
  name  = c("Northeast", "Northwest", "Ohio Valley", "Northern Rockies", "Upper Midwest", "Southeast", "Southwest", "West", "South"),
  col   = c("#008B8B", "#2E8B57", "#A16207", "#8B8B8B", "#D97706", "#F97316", "#E65100", "#1D4ED8", "#EF4444"),
  x_box = c(0.01, 0.35, 0.69, 0.01, 0.35, 0.69, 0.01, 0.35, 0.69),
  y_pos = c(0.22, 0.22, 0.22, 0.13, 0.13, 0.13, 0.04, 0.04, 0.04)
)

# ==============================================================================
# 2. FIGURE 1A: CUMULATIVE OUTAGE HOURS HAC SPATIAL TREND
# ==============================================================================
message("[INFO] Processing and saving Figure 1a...")

calc_hac_trend_1a <- function(df_sub) {
  if (all(df_sub$monthly_hours == 0, na.rm = TRUE)) {
    return(tibble(slope = 0, p_val = 1.0))
  }
  
  mod <- tryCatch(
    lm(monthly_hours ~ year_frac + factor(month), data = df_sub),
    error = function(e) NULL
  )
  
  if (is.null(mod)) return(tibble(slope = 0, p_val = 1.0))
  
  hac_cov <- tryCatch(sandwich::vcovHAC(mod), error = function(e) NULL)
  
  if (is.null(hac_cov)) {
    slope <- coef(mod)["year_frac"] * 12
    p_val <- summary(mod)$coefficients["year_frac", "Pr(>|t|)"]
  } else {
    ct    <- lmtest::coeftest(mod, vcov. = hac_cov)
    slope <- ct["year_frac", "Estimate"] * 12
    p_val <- ct["year_frac", "Pr(>|t|)"]
  }
  
  return(tibble(slope = as.numeric(slope), p_val = as.numeric(p_val)))
}

county_trends_1a <- df_all %>%
  group_by(fips_std, year, month = month_num) %>%
  summarise(monthly_hours = sum(outage_hours_val, na.rm = TRUE), .groups = "drop") %>%
  complete(fips_std = all_conus_fips, year = 2015:2024, month = 1:12, fill = list(monthly_hours = 0)) %>%
  mutate(year_frac = year + (month - 1) / 12) %>%
  group_by(fips_std) %>%
  group_modify(~calc_hac_trend_1a(.x)) %>%
  ungroup() %>%
  rename(fips = fips_std)

map_data_1a <- counties_sf %>%
  left_join(county_trends_1a, by = "fips") %>%
  mutate(poly_alpha = ifelse(!is.na(p_val) & p_val < 0.05, 1.0, 0.25))

map_colors_1a <- c("#08519c", "#3182bd", "#6baed6", "#9ecae1", "#fcfbfd", "#fee3d6", "#fc9272", "#de2d26", "#a50f15", "#67000d")

plot_map_1a <- ggplot(map_data_1a) +
  geom_sf(aes(fill = slope, alpha = poly_alpha), color = "gray88", linewidth = 0.04) +
  scale_alpha_identity() +
  geom_sf(data = ncei_regions_sf, aes(color = region), fill = NA, linewidth = 1.2) +
  scale_fill_stepsn(
    colors   = map_colors_1a,
    breaks   = c(-100, -50, -10, 0, 10, 50, 100),
    limits   = c(-105, 105),
    oob      = scales::squish
  ) +
  scale_color_manual(values = ncei_colors) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 0, l = 10, unit = "pt"))

plot_cb_1a <- ggplot(data.frame(x = seq(-105, 105, length.out = 150), y = 1), aes(x = x, y = y, fill = x)) +
  geom_tile(height = 1.5) +
  scale_fill_gradientn(colors = map_colors_1a) +
  scale_x_continuous(breaks = c(-100, -50, 0, 50, 100), labels = c("≤ -100", "-50", "0", "50", "≥ 100")) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x      = element_text(size = FONT_TICK_TEXT, face = "bold", color = "black", margin = margin(t = 6))
  )

legend_canvas_1a <- ggdraw() +
  draw_label("Trend (outage hours per year)", x = 0.28, y = 0.96, fontface = "bold", size = FONT_LEG_TITLE) +
  draw_plot(plot_cb_1a, x = 0.04, y = 0.70, width = 0.50, height = 0.16) +
  draw_label("Decrease", x = 0.08, y = 0.54, fontface = "bold", color = "#08519c", size = FONT_LEG_TEXT) +
  draw_label("No trend",  x = 0.29, y = 0.54, fontface = "bold", color = "black",   size = FONT_LEG_TEXT) +
  draw_label("Increase", x = 0.50, y = 0.54, fontface = "bold", color = "#a50f15", size = FONT_LEG_TEXT) +
  draw_grob(grid::rectGrob(
    x      = unit(0.76, "npc"),
    y      = unit(0.76, "npc"),
    width  = unit(0.35, "npc"),
    height = unit(0.24, "npc"),
    gp     = grid::gpar(fill = "#f8f9fa", col = "black", lwd = 1.2)
  )) +
  draw_label("P < 0.05: 100% Opacity",    x = 0.59, y = 0.83, fontface = "bold", size = FONT_LEG_TEXT, hjust = 0) +
  draw_label("P \u2265 0.05: 25% Opacity", x = 0.59, y = 0.69, fontface = "bold", size = FONT_LEG_TEXT, hjust = 0) +
  draw_label("NCEI Climate Regions", x = 0.50, y = 0.35, fontface = "bold", size = FONT_LEG_TITLE)

for (i in 1:nrow(ncei_items_3x3_fixed)) {
  item <- ncei_items_3x3_fixed[i, ]
  legend_canvas_1a <- legend_canvas_1a +
    draw_grob(grid::rectGrob(x = unit(item$x_box, "npc"), y = unit(item$y_pos, "npc"), width = unit(7.0, "mm"), height = unit(5.5, "mm"), gp = gpar(fill = "transparent", col = item$col, lwd = 3.2))) +
    draw_label(item$name, x = item$x_box + 0.045, y = item$y_pos, fontface = "bold", size = FONT_LEG_TEXT, hjust = 0)
}

plot_1a_final <- plot_grid(plot_map_1a, legend_canvas_1a, ncol = 1, rel_heights = c(1, 0.52))

ggsave(
  filename = file.path(output_dir, "Figure_1a_Outage_Hours_Final.png"),
  plot     = plot_1a_final,
  width    = PANEL_WIDTH,
  height   = PANEL_HEIGHT,
  dpi      = 300,
  bg       = "white"
)

# ==============================================================================
# 3. FIGURE 1B: REGIONAL HETEROGENEITY (BROKEN AXIS) [SWAPPED WITH 1C]
# ==============================================================================
message("[INFO] Processing and saving Figure 1b (Regional Heterogeneity)...")

sf::sf_use_s2(FALSE)
county_areas <- usmap::us_map(regions = "counties") %>%
  st_as_sf() %>%
  mutate(
    fips_clean = str_pad(as.character(fips), width = 5, pad = "0"),
    area_km2   = as.numeric(sf::st_area(.)) / 1e6
  ) %>%
  st_drop_geometry() %>%
  select(fips, area_km2) %>%
  distinct(fips, .keep_all = TRUE)
sf::sf_use_s2(TRUE)

wtd_quantile <- function(x, w, probs = 0.5) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  if (length(x) == 0) return(NA_real_)
  ord <- order(x)
  x   <- x[ord]
  w   <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  sapply(probs, function(p) {
    idx <- which(cum_w >= p)[1]
    if (is.na(idx)) return(x[length(x)])
    x[idx]
  })
}

df_summary <- county_trends_1a %>%
  filter(!is.na(slope)) %>%
  mutate(
    fips_clean = str_pad(as.character(fips), width = 5, side = "left", pad = "0"),
    state_fips = substr(fips_clean, 1, 2),
    state_abb  = fips_to_state_map[state_fips],
    region     = state_to_ncei[state_abb]
  ) %>%
  filter(!is.na(region)) %>%
  left_join(county_areas, by = c("fips_clean" = "fips")) %>%
  mutate(area_km2 = if_else(is.na(area_km2), median(area_km2, na.rm = TRUE), area_km2)) %>%
  group_by(region) %>%
  summarise(
    med_wtd          = wtd_quantile(slope[!is.na(p_val) & p_val < 0.05], area_km2[!is.na(p_val) & p_val < 0.05], probs = 0.5),
    q25_wtd          = wtd_quantile(slope[!is.na(p_val) & p_val < 0.05], area_km2[!is.na(p_val) & p_val < 0.05], probs = 0.25),
    q75_wtd          = wtd_quantile(slope[!is.na(p_val) & p_val < 0.05], area_km2[!is.na(p_val) & p_val < 0.05], probs = 0.75),
    pct_area_sig_pos = (sum(area_km2[slope > 0 & (!is.na(p_val) & p_val < 0.05)], na.rm = TRUE) / sum(area_km2, na.rm = TRUE)) * 100,
    .groups          = "drop"
  ) %>%
  arrange(desc(med_wtd), desc(pct_area_sig_pos)) %>%
  mutate(
    region      = factor(region, levels = rev(region)),
    fmt_pct_sig = sprintf("%.1f%%", pct_area_sig_pos)
  )

df_west   <- df_summary %>% filter(region == "West")
df_others <- df_summary %>% filter(region != "West")

p_west <- ggplot(df_west, aes(y = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_segment(aes(x = -25, xend = 0, y = region, yend = region, color = region), alpha = 0.35, linewidth = 0.8) +
  geom_errorbarh(aes(xmin = q25_wtd, xmax = q75_wtd, color = region), width = 0.3, linewidth = 1.3) +
  geom_point(aes(x = -25, color = region), size = 6.5) +
  geom_point(aes(x = med_wtd, color = region), size = 6.2) +
  geom_text(aes(x = -35, label = as.character(region)), color = "black", hjust = 1, fontface = "bold", size = 8.5) +
  geom_text(aes(x = 280, label = fmt_pct_sig), color = "black", hjust = 1, fontface = "bold", size = 8.0) +
  annotate("text", x = 280, y = 1.85, label = "% of regional area\nsig. pos. trend", fontface = "bold", size = 7.5, hjust = 1, vjust = 0, lineheight = 0.95) +
  scale_color_manual(values = ncei_colors) +
  scale_x_continuous(breaks = c(0, 100, 200), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(-140, 290), ylim = c(0.5, 2.5), clip = "off") +
  theme_classic(base_size = 20) +
  theme(
    legend.position  = "none",
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    axis.title.x     = element_blank(),
    axis.line.x      = element_line(color = "black", linewidth = 1.0),
    axis.text.x      = element_text(face = "bold", size = FONT_TICK_TEXT, color = "black"),
    plot.margin      = margin(t = 25, r = 25, b = 4, l = 100)
  )

p_others <- ggplot(df_others, aes(y = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_segment(aes(x = -5, xend = 0, y = region, yend = region, color = region), alpha = 0.35, linewidth = 0.8) +
  geom_errorbarh(aes(xmin = q25_wtd, xmax = q75_wtd, color = region), width = 0.3, linewidth = 1.3) +
  geom_point(aes(x = -5, color = region), size = 6.5) +
  geom_point(aes(x = med_wtd, color = region), size = 6.2) +
  geom_text(aes(x = -7, label = as.character(region)), color = "black", hjust = 1, fontface = "bold", size = 8.5) +
  geom_text(aes(x = 62, label = fmt_pct_sig), color = "black", hjust = 1, fontface = "bold", size = 8.0) +
  scale_color_manual(values = ncei_colors) +
  scale_x_continuous(breaks = c(0, 20, 40), expand = c(0, 0)) +
  labs(x = "Area-weighted median trend (outage hours per year)", y = NULL) +
  coord_cartesian(xlim = c(-28, 64), clip = "off") +
  theme_classic(base_size = 20) +
  theme(
    legend.position  = "none",
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.y     = element_blank(),
    axis.title.x     = element_text(face = "bold", size = FONT_AXIS_TITLE, color = "black", margin = margin(t = 10)),
    axis.line.x      = element_line(color = "black", linewidth = 1.0),
    axis.text.x      = element_text(face = "bold", size = FONT_TICK_TEXT, color = "black"),
    plot.margin      = margin(t = 4, r = 25, b = 10, l = 100)
  )

plot_1b_final <- (p_west / p_others) + plot_layout(heights = c(1.4, 8))

ggsave(
  filename = file.path(output_dir, "Figure_1b_Regional_Outage_Hours_Trends_Final.png"),
  plot     = plot_1b_final,
  width    = PANEL_WIDTH,
  height   = PANEL_HEIGHT,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# ==============================================================================
# FIGURE C: ANNUAL CUMULATIVE OUTAGE DURATION & COH TIME SERIES (ENLARGED)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

message("[INFO] Processing and saving Figure C...")

# ------------------------------------------------------------------------------
# 1. TYPOGRAPHY & COLOR CONSTANTS (ENLARGED)
# ------------------------------------------------------------------------------
FONT_AXIS_TITLE_C <- 30   # Axis Titles (Left Y, Right Y, X)
FONT_TICK_TEXT_C  <- 26   # Axis Tick Labels
FONT_LEG_TEXT_C   <- 26   # Legend Text Labels
ANNOT_TEXT_SIZE_C <- 9.0  # Headers & Mean Text Labels (~26pt font, identical size)

COLOR_BAR  <- "#f05a3e"   # Coral/Red-orange for Outage Duration bars
COLOR_LINE <- "#0f2c4d"   # Dark Navy for COH points/lines

# ------------------------------------------------------------------------------
# 2. DATA AGGREGATION
# ------------------------------------------------------------------------------
annual_summary <- df_all %>%
  filter(year >= 2015 & year <= 2024) %>%
  group_by(year) %>%
  summarise(
    outage_duration_m = sum(outage_hours_val, na.rm = TRUE) / 1e6, # Million hours
    coh_m             = sum(coh_val, na.rm = TRUE) / 1e6,          # Million customer-hours
    .groups = "drop"
  )

# Compute historical and recent mean outage durations
mean_hist   <- mean(annual_summary$outage_duration_m[annual_summary$year %in% 2015:2018], na.rm = TRUE)
mean_recent <- mean(annual_summary$outage_duration_m[annual_summary$year %in% 2019:2022], na.rm = TRUE)

# ------------------------------------------------------------------------------
# 3. DUAL-AXIS SCALING REGIME
# ------------------------------------------------------------------------------
y1_max    <- 0.85
y2_max    <- 2000
scale_fac <- y2_max / y1_max  # Scale factor

# ------------------------------------------------------------------------------
# 4. PLOT BUILD
# ------------------------------------------------------------------------------
plot_c <- ggplot(annual_summary, aes(x = year)) +
  # Background Shaded Period Rectangles
  annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 0, ymax = y1_max, fill = "#f0f4f8", alpha = 0.7) +
  annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 0, ymax = y1_max, fill = "#fdf0ed", alpha = 0.7) +
  
  # Primary Y Geom: Outage Duration Bars
  geom_col(
    aes(y = outage_duration_m, fill = "Cumulative outage duration\n(million hours)"),
    width = 0.55, alpha = 0.9
  ) +
  
  # Secondary Y Geom: Thicker COH Line & Larger Points
  geom_line(
    aes(y = coh_m / scale_fac, color = "High-impact COH\n(million customer-hours)"),
    linewidth = 2.2
  ) +
  geom_point(
    aes(y = coh_m / scale_fac, color = "High-impact COH\n(million customer-hours)"),
    size = 7.0
  ) +
  
  # Vertical Section Divider Lines
  geom_vline(xintercept = 2018.5, linetype = "dashed", color = "gray50", linewidth = 1.0) +
  geom_vline(xintercept = 2022.5, linetype = "dashed", color = "gray50", linewidth = 1.0) +
  
  # Period Mean Horizontal Dashed Lines
  annotate("segment", x = 2014.5, xend = 2018.5, y = mean_hist, yend = mean_hist,
           color = "#bd3217", linetype = "dashed", linewidth = 1.2) +
  annotate("segment", x = 2018.5, xend = 2022.5, y = mean_recent, yend = mean_recent,
           color = "#bd3217", linetype = "dashed", linewidth = 1.2) +
  
  # Period Header Annotations
  annotate("text", x = 2016.5, y = 0.77, label = "Historical\n(2015–2018)",
           fontface = "bold", size = ANNOT_TEXT_SIZE_C, color = "#0f2c4d", lineheight = 0.9) +
  annotate("text", x = 2020.5, y = 0.77, label = "Recent\n(2019–2022)",
           fontface = "bold", size = ANNOT_TEXT_SIZE_C, color = "#a8240d", lineheight = 0.9) +
  annotate("text", x = 2023.5, y = 0.77, label = "Outage obs. only\n(2023–2024)",
           fontface = "bold", size = ANNOT_TEXT_SIZE_C, color = "black", lineheight = 0.9) +
  
  # Period Mean Text Labels (SAME SIZE AS HEADERS: ANNOT_TEXT_SIZE_C)
  annotate("text", x = 2016.5, y = mean_hist + 0.18, label = sprintf("Mean: %.2f\nmillion hours", mean_hist),
           fontface = "bold", size = ANNOT_TEXT_SIZE_C, color = "#0f2c4d", lineheight = 0.9) +
  annotate("text", x = 2020.5, y = mean_recent + 0.14, label = sprintf("Mean: %.2f\nmillion hours", mean_recent),
           fontface = "bold", size = ANNOT_TEXT_SIZE_C, color = "#a8240d", lineheight = 0.9) +
  
  # Aesthetic Mapping
  scale_fill_manual(name = NULL, values = c("Cumulative outage duration\n(million hours)" = COLOR_BAR)) +
  scale_color_manual(name = NULL, values = c("High-impact COH\n(million customer-hours)" = COLOR_LINE)) +
  
  # Axis Configuration
  scale_x_continuous(breaks = 2015:2024, expand = c(0.02, 0.02)) +
  scale_y_continuous(
    name   = "Cumulative outage duration\n(million hours)",
    limits = c(0, y1_max),
    breaks = seq(0, 0.8, by = 0.1),
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . * scale_fac,
      name   = "High-impact COH\n(million customer-hours)",
      breaks = seq(0, y2_max, by = 400),
      labels = comma
    )
  ) +
  labs(x = "Year") +
  theme_minimal(base_size = 22) +
  theme(
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_line(color = "#e5e5e5", linewidth = 0.6),
    
    # Left Axis Styling (Coral/Red-Orange)
    axis.line.y.left     = element_line(color = COLOR_BAR, linewidth = 1.6),
    axis.ticks.y.left    = element_line(color = COLOR_BAR, linewidth = 1.6),
    axis.text.y.left     = element_text(face = "bold", size = FONT_TICK_TEXT_C, color = COLOR_BAR),
    axis.title.y.left    = element_text(face = "bold", size = FONT_AXIS_TITLE_C, color = COLOR_BAR, margin = margin(r = 18)),
    
    # Right Axis Styling (Dark Navy)
    axis.line.y.right    = element_line(color = COLOR_LINE, linewidth = 1.6),
    axis.ticks.y.right   = element_line(color = COLOR_LINE, linewidth = 1.6),
    axis.text.y.right    = element_text(face = "bold", size = FONT_TICK_TEXT_C, color = COLOR_LINE),
    axis.title.y.right   = element_text(face = "bold", size = FONT_AXIS_TITLE_C, color = COLOR_LINE, angle = 90, margin = margin(l = 18)),
    
    # X Axis Styling
    axis.line.x          = element_line(color = "black", linewidth = 1.6),
    axis.ticks.x         = element_line(color = "black", linewidth = 1.6),
    axis.text.x          = element_text(face = "bold", size = FONT_TICK_TEXT_C, color = "black"),
    axis.title.x         = element_text(face = "bold", size = FONT_AXIS_TITLE_C, color = "black", margin = margin(t = 15)),
    
    # Legend Placement & Layout
    legend.position      = "top",
    legend.direction     = "horizontal",
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.text          = element_text(size = FONT_LEG_TEXT_C, face = "bold", lineheight = 0.9),
    legend.key.size      = unit(1.4, "cm"),
    legend.spacing.x     = unit(0.6, "cm"),
    legend.margin        = margin(b = 20),
    
    plot.margin          = margin(t = 25, r = 30, b = 20, l = 30, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.9)),
    color = guide_legend(order = 2, override.aes = list(linewidth = 2.2, size = 7.0))
  )

# ------------------------------------------------------------------------------
# 5. SAVE OUTPUT
# ------------------------------------------------------------------------------
ggsave(
  filename = file.path(output_dir, "Figure_C_Annual_Time_Series.png"),
  plot     = plot_c,
  width    = 16.0,
  height   = 10.0,
  dpi      = 300,
  bg       = "white"
)


# Print Figure 1b summary table to console
print(df_summary %>% select(
  Region = region,
  `Median Trend (hrs/yr)` = med_wtd,
  `25th Pct` = q25_wtd,
  `75th Pct` = q75_wtd,
  `Sig Pos Area (%)` = fmt_pct_sig
))

# Export to CSV
write_csv(df_summary, file.path(output_dir, "Table_Figure_1b_Regional_Statistics.csv"))


# Add Period labels and Period Means to the dataframe
annual_summary_table <- annual_summary %>%
  mutate(
    Period = case_when(
      year %in% 2015:2018 ~ "Historical (2015–2018)",
      year %in% 2019:2022 ~ "Recent (2019–2022)",
      year %in% 2023:2024 ~ "Outage Obs. Only (2023–2024)"
    ),
    Period_Mean_Outage_M = case_when(
      year %in% 2015:2018 ~ mean_hist,
      year %in% 2019:2022 ~ mean_recent,
      TRUE ~ NA_real_
    )
  )

# Print Figure C annual table
print(annual_summary_table)

# Export to CSV
write_csv(annual_summary_table, file.path(output_dir, "Table_Figure_C_Annual_Statistics.csv"))

# ==============================================================================
# MAXIMUM-TYPOGRAPHY SETUP FOR 18-INCH PANELS
# ==============================================================================
FONT_AXIS_TITLE_18 <- 28   # Axis Titles (Y-left, Y-right)
FONT_TICK_TEXT_18  <- 24   # Axis Ticks (X & Y)
FONT_LEG_TEXT_18   <- 24   # Legend Text Labels
GEOM_TEXT_SIZE_18  <- 8.0  # Bar Percentage Values (8.0 pt geom ~ 23pt font)

# ==============================================================================
# 5. FIGURE 1D (i): MONTHLY DUAL-AXIS OUTAGE DAYS & CHANGE (EXTRA LARGE LABELS)
# ==============================================================================
message("[INFO] Processing and saving Figure 1d (i) with enlarged labels...")

monthly_comp <- df_all %>%
  filter(year >= 2015 & year <= 2022) %>%
  mutate(
    period = case_when(
      year >= 2015 & year <= 2018 ~ "Early",
      year >= 2019 & year <= 2022 ~ "Recent",
      TRUE ~ NA_character_
    ),
    month = month(date_clean, label = TRUE, abbr = TRUE)
  ) %>%
  filter(!is.na(period)) %>%
  filter(outage_hours_val > 0) %>% 
  group_by(month, period) %>%
  summarise(total_outage = n_distinct(paste(fips_std, date_clean)), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = total_outage, values_fill = 0) %>%
  mutate(pct_change = ifelse(Early > 0, ((Recent - Early) / Early) * 100, 0))

# Dynamic scale factor with 35% extra top margin for large bar labels
y_max_data  <- max(c(monthly_comp$Early, monthly_comp$Recent), na.rm = TRUE)
pct_max     <- max(monthly_comp$pct_change, na.rm = TRUE)
y_limit_top <- ceiling(y_max_data * 1.35)
scale_fac   <- y_limit_top / (ceiling(pct_max / 50) * 50)

plot_1d_i <- ggplot(monthly_comp, aes(x = month)) +
  geom_col(aes(y = pct_change * scale_fac, fill = "Change (%)"), width = 0.45, alpha = 0.7) +
  geom_text(
    aes(y = pct_change * scale_fac + (y_limit_top * 0.035), label = sprintf("%.0f%%", pct_change)),
    vjust = 0, color = COLOR_CHANGE_TEXT, fontface = "bold", size = GEOM_TEXT_SIZE_18
  ) +
  geom_line(aes(y = Early, group = 1, color = "Historical (2015–2018)"), linewidth = 1.5) +
  geom_point(aes(y = Early, color = "Historical (2015–2018)"), shape = 16, size = 5.5) +
  geom_line(aes(y = Recent, group = 1, color = "Recent (2019–2022)"), linewidth = 1.5) +
  geom_point(aes(y = Recent, color = "Recent (2019–2022)"), shape = 16, size = 5.5) +
  scale_color_manual(
    name   = NULL,
    values = c("Historical (2015–2018)" = "#0045a5", "Recent (2019–2022)" = "#e31a1c")
  ) +
  scale_fill_manual(name = NULL, values = c("Change (%)" = "#fcae91")) +
  scale_y_continuous(
    name     = "High-impact outage days\n(county-days)",
    expand   = c(0, 0),
    labels   = comma,
    sec.axis = sec_axis(~ . / scale_fac, name = "Change (%)", labels = function(x) paste0(x, "%"))
  ) +
  coord_cartesian(ylim = c(0, y_limit_top), clip = "off") +
  labs(x = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_line(color = "#e5e5e5", linewidth = 0.6),
    axis.line.x          = element_line(color = "black", linewidth = 1.4),
    axis.ticks.x         = element_line(color = "black", linewidth = 1.4),
    axis.text.x          = element_text(face = "bold", size = FONT_TICK_TEXT_18, color = "black"),
    axis.line.y.left     = element_line(color = "black", linewidth = 1.4),
    axis.ticks.y.left    = element_line(color = "black", linewidth = 1.4),
    axis.text.y.left     = element_text(face = "bold", size = FONT_TICK_TEXT_18, color = "black"),
    axis.title.y.left    = element_text(face = "bold", size = FONT_AXIS_TITLE_18, color = "black", margin = margin(r = 16)),
    axis.line.y.right    = element_line(color = COLOR_CHANGE_TEXT, linewidth = 1.4),
    axis.ticks.y.right   = element_line(color = COLOR_CHANGE_TEXT, linewidth = 1.4),
    axis.text.y.right    = element_text(face = "bold", size = FONT_TICK_TEXT_18, color = COLOR_CHANGE_TEXT),
    axis.title.y.right   = element_text(face = "bold", size = FONT_AXIS_TITLE_18, color = COLOR_CHANGE_TEXT, angle = 90, margin = margin(l = 16)),
    legend.position      = c(0.01, 0.98),
    legend.justification = c(0, 1),
    legend.direction     = "horizontal",
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.text          = element_text(size = FONT_LEG_TEXT_18, face = "bold"),
    legend.key.size      = unit(1.0, "cm"),
    legend.spacing.x     = unit(0.4, "cm"),
    plot.margin          = margin(t = 35, r = 25, b = 15, l = 25, unit = "pt")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 1.5, size = 5.5)),
    fill  = guide_legend(order = 2)
  )

ggsave(
  filename = file.path(output_dir, "Figure_1d_i.png"),
  plot     = plot_1d_i,
  width    = 18.0,
  height   = PANEL_HEIGHT,
  dpi      = 300,
  bg       = "white"
)

# ==============================================================================
# 6. FIGURE 1D (ii): SEASONAL CONTRIBUTION (EXTRA LARGE LABELS)
# ==============================================================================
message("[INFO] Processing and saving Figure 1d (ii) with enlarged labels...")

monthly_comp_1d <- df_all %>%
  filter(year >= 2015 & year <= 2022) %>%
  mutate(
    period = case_when(
      year >= 2015 & year <= 2018 ~ "Early",
      year >= 2019 & year <= 2022 ~ "Recent",
      TRUE ~ NA_character_
    ),
    month = month(date_clean, label = TRUE, abbr = TRUE)
  ) %>%
  filter(!is.na(period)) %>%
  group_by(month, period) %>%
  summarise(total_outage = sum(outage_hours_val, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = total_outage, values_fill = 0)

season_map <- c(
  "Jan" = "Winter\n(Dec–Feb)", "Feb" = "Winter\n(Dec–Feb)", "Dec" = "Winter\n(Dec–Feb)",
  "Mar" = "Spring\n(Mar–May)", "Apr" = "Spring\n(Mar–May)", "May" = "Spring\n(Mar–May)",
  "Jun" = "Summer\n(Jun–Aug)", "Jul" = "Summer\n(Jun–Aug)", "Aug" = "Summer\n(Jun–Aug)",
  "Sep" = "Autumn\n(Sep–Nov)", "Oct" = "Autumn\n(Sep–Nov)", "Nov" = "Autumn\n(Sep–Nov)"
)

seasonal_comp <- monthly_comp_1d %>%
  mutate(season = season_map[as.character(month)]) %>%
  group_by(season) %>%
  summarise(
    Early  = sum(Early, na.rm = TRUE),
    Recent = sum(Recent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Early, Recent), names_to = "period", values_to = "outage_val") %>%
  group_by(period) %>%
  mutate(
    pct          = (outage_val / sum(outage_val)) * 100,
    period_label = ifelse(period == "Early", "Historical (2015–2018)", "Recent (2019–2022)"),
    season       = factor(season, levels = c("Winter\n(Dec–Feb)", "Spring\n(Mar–May)", "Summer\n(Jun–Aug)", "Autumn\n(Sep–Nov)"))
  ) %>%
  arrange(period, season)

plot_1d_ii <- ggplot(seasonal_comp, aes(x = season, y = pct, fill = period_label)) +
  geom_col(position = position_dodge(width = 0.55), width = 0.45) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    position = position_dodge(width = 0.55),
    vjust = -0.4,
    fontface = "bold",
    size = GEOM_TEXT_SIZE_18,
    color = "black"
  ) +
  scale_fill_manual(name = NULL, values = c("Historical (2015–2018)" = "#0045a5", "Recent (2019–2022)" = "#e31a1c")) +
  scale_y_continuous(
    name   = "Seasonal share of total outage hours (%)",
    limits = c(0, 58),
    breaks = seq(0, 50, by = 10),
    expand = c(0, 0)
  ) +
  labs(x = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.grid.major.y   = element_line(color = "#e5e5e5", linewidth = 0.6),
    axis.line.x          = element_line(color = "black", linewidth = 1.4),
    axis.ticks.x         = element_line(color = "black", linewidth = 1.4),
    axis.text.x          = element_text(face = "bold", size = FONT_TICK_TEXT_18, color = "black"),
    axis.line.y          = element_line(color = "black", linewidth = 1.4),
    axis.ticks.y         = element_line(color = "black", linewidth = 1.4),
    axis.text.y          = element_text(face = "bold", size = FONT_TICK_TEXT_18, color = "black"),
    axis.title.y         = element_text(face = "bold", size = FONT_AXIS_TITLE_18, color = "black", margin = margin(r = 16)),
    legend.position      = c(0.01, 0.98),
    legend.justification = c(0, 1),
    legend.direction     = "horizontal",
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.text          = element_text(size = FONT_LEG_TEXT_18, face = "bold"),
    legend.key.size      = unit(1.0, "cm"),
    legend.spacing.x     = unit(0.4, "cm"),
    plot.margin          = margin(t = 35, r = 25, b = 15, l = 25, unit = "pt")
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(output_dir, "Figure_1d_ii.png"),
  plot     = plot_1d_ii,
  width    = 18.0,
  height   = PANEL_HEIGHT,
  dpi      = 300,
  bg       = "white"
)

