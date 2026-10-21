# =============================================================================
# Script: Figure_4c_threshold_shift_maps.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description: Figure 4c: state-level maps of standardized threshold shifts.
# Input: Figure 3b ALE zero-crossing threshold summary CSV files.
# Dependencies: dplyr, tidyr, ggplot2, stringr, usmap, sf, patchwork, ggrepel
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(usmap)      # For US CONUS map geometry
library(sf)         # For spatial feature handling
library(patchwork)  # For multi-panel grid layout
library(ggrepel)    # For non-overlapping state callout text with arrows

input_base_dir <- Sys.getenv("FIG3_INPUT_DIR", unset = "data/processed")
output_dir     <- Sys.getenv("FIG4_OUTPUT_DIR", unset = "results/figure3")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

non_conus_abbr <- c("AK", "HI", "PR", "VI", "GU", "MP", "AS")

# ==============================================================================
# VARIABLE CONFIGURATIONS & TITLES
# ==============================================================================
var_config_4c <- list(
  "mean_T2C" = list(
    title = "Mean temperature (°C)", 
    unit = "°C",
    show_unit = TRUE
  ),
  "monthly_CDD" = list(
    title = "CDD (base 18 °C)", 
    unit = "°C-days",
    show_unit = TRUE
  ),
  "monthly_HDD" = list(
    title = "HDD (base 18 °C)", 
    unit = "°C-days",
    show_unit = TRUE
  ),
  "heat_days" = list(
    title = "Extreme heat days (days month⁻¹)", 
    unit = "",
    show_unit = FALSE
  ),
  "cold_days" = list(
    title = "Extreme cold days (days month⁻¹)", 
    unit = "",
    show_unit = FALSE
  ),
  "monthly_max_mag" = list(
    title = "Maximum magnitude (severity)", 
    unit = "",
    show_unit = FALSE
  ),
  "storm_days" = list(
    title = "Storm duration (days month⁻¹)", 
    unit = "",
    show_unit = FALSE
  ),
  "monthly_event_count" = list(
    title = "Event frequency (events month⁻¹)", 
    unit = "",
    show_unit = FALSE
  )
)

# ==============================================================================
# HELPER: BUILD A SINGLE MAP PANEL WITH TITLE & INSIDE LABELS (NO BOX)
# ==============================================================================
build_single_map <- function(df_var, var_key) {
  cfg <- var_config_4c[[var_key]]
  
  us_states <- usmap::us_map(regions = "states") %>%
    filter(!abbr %in% non_conus_abbr)
  
  if (inherits(us_states, "sf")) {
    centroids <- us_states %>%
      sf::st_centroid(of_largest_polygon = TRUE) %>%
      mutate(x = sf::st_coordinates(.)[, 1], y = sf::st_coordinates(.)[, 2]) %>%
      sf::st_drop_geometry() %>% select(abbr, x, y)
  } else {
    centroids <- us_states %>%
      group_by(abbr) %>% summarise(x = mean(x), y = mean(y), .groups = "drop")
  }
  
  map_df <- us_states %>% left_join(df_var, by = "abbr")
  
  # Calculate 25th and 75th percentiles
  q25 <- quantile(df_var$delta_threshold, 0.25, na.rm = TRUE)
  q75 <- quantile(df_var$delta_threshold, 0.75, na.rm = TRUE)
  
  # Filter within 25th–75th percentile range
  filtered_df <- df_var %>%
    filter(!is.na(delta_threshold)) %>%
    filter(delta_threshold >= q25 & delta_threshold <= q75) %>%
    arrange(delta_threshold)
  
  n_count <- nrow(filtered_df)
  
  # Ensure exactly top 2 lowest and top 2 highest are selected if n >= 4
  if (n_count >= 4) {
    callout_df <- filtered_df %>% slice(c(1, 2, n_count - 1, n_count))
  } else {
    callout_df <- filtered_df
  }
  
  callout_df <- callout_df %>%
    mutate(
      val_formatted = if_else(
        abs(delta_threshold) >= 10,
        sprintf("%+d", round(delta_threshold)),
        sprintf("%+.1f", delta_threshold)
      ),
      label_str = if (cfg$show_unit) {
        paste0(abbr, "\n", val_formatted, " ", cfg$unit)
      } else {
        paste0(abbr, "\n", val_formatted)
      }
    )
  
  callout_coords <- centroids %>% inner_join(callout_df, by = "abbr")
  
  p <- ggplot()
  
  if (inherits(us_states, "sf")) {
    p <- p + geom_sf(
      data = map_df, aes(fill = delta_std_clamped),
      color = "#ffffff", linewidth = 0.6
    ) + 
      coord_sf(default_crs = NULL, clip = "off")  # <-- ALLOW LABELS OUTSIDE MAP CANVAS
  } else {
    p <- p + geom_polygon(
      data = map_df, aes(x = x, y = y, group = group, fill = delta_std_clamped),
      color = "#ffffff", linewidth = 0.6
    ) + 
      coord_equal(clip = "off")                   # <-- ALLOW LABELS OUTSIDE MAP CANVAS
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#e33935", mid = "#ffffff", high = "#1d70b8", midpoint = 0,
      limits = c(-2, 2), na.value = "#cccccc", guide = "none"
    ) +
    geom_text_repel(
      data = callout_coords,
      aes(x = x, y = y, label = label_str),
      size = 7.0,                         # Slightly tuned size to fit dense areas
      fontface = "bold",
      color = "#111111",
      bg.color = "#ffffff",
      bg.r = 0.15,
      box.padding = 0.3,
      point.padding = 0.2,
      force = 3,                          # Pushes overlapping labels apart
      segment.color = "#222222",
      segment.size = 0.8,
      arrow = arrow(length = unit(0.012, "npc"), type = "closed", ends = "first"),
      min.segment.length = 0,
      max.overlaps = Inf                  # Forces ALL labels to render regardless of overlap
    ) +
    labs(title = cfg$title) +
    theme_void() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 30, face = "bold", color = "#111111", hjust = 0, margin = margin(b = 8, l = 2)),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 15) # Expanded margins so labels aren't clipped
    )
  
  return(p)
}

# ==============================================================================
# HELPER: VERTICAL SECTION BANNER GRAPHIC (SIDEWAYS LABEL)
# ==============================================================================
build_vertical_banner <- function(title_text, bg_color = "#dbebf9") {
  ggplot() +
    annotate(
      "rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
      fill = bg_color, color = NA
    ) +
    annotate(
      "text", x = 0.5, y = 0.5, label = title_text,
      fontface = "bold", size = 12.0, color = "#111111", angle = 90, hjust = 0.5, vjust = 0.5
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(r = 4, l = 2))
}

# ==============================================================================
# HELPER: FIXED LEGEND GRAPHIC
# ==============================================================================
build_fig4c_legend <- function() {
  legend_df <- data.frame(x = seq(-2, 2, length.out = 300), y = 0)
  
  ggplot() +
    # Continuous Gradient Bar
    geom_tile(data = legend_df, aes(x = x, y = y, fill = x), height = 0.24) +
    scale_fill_gradient2(
      low = "#e33935", mid = "#ffffff", high = "#1d70b8", midpoint = 0
    ) +
    
    # Text safely ABOVE gradient bar
    annotate("text", x = -2, y = 0.55, label = "Lower threshold (in recent period)", 
             fontface = "bold", size = 8.5, color = "#e33935", hjust = 0.5, lineheight = 0.95) +
    annotate("text", x =  2, y = 0.55, label = "Higher threshold (in recent period)", 
             fontface = "bold", size = 8.5, color = "#1d70b8", hjust = 0.5, lineheight = 0.95) +
    
    # Tick Marks & Numbers BELOW gradient bar
    annotate("segment", x = -2, xend = 2, y = -0.18, yend = -0.18, color = "#333333", linewidth = 1.0) +
    annotate("text", x = c(-2, -1, 0, 1, 2), y = -0.42, label = c("-2", "-1", "0", "1", "2"), 
             fontface = "bold", size = 9.0) +
    
    # Axis Title at bottom
    annotate("text", x = 0, y = -0.75, label = "Standardized threshold shift (Δ / historical IQR)", 
             fontface = "bold", size = 9.5, color = "#111111") +
    
    # Grey Box & Label on Right Side
    annotate("rect", xmin = 2.8, xmax = 3.15, ymin = -0.12, ymax = 0.12, fill = "#cccccc", color = "#aaaaaa", linewidth = 0.5) +
    annotate("text", x = 3.25, y = 0, label = "Not estimable\n(insufficient data or no threshold)", 
             fontface = "bold", size = 7.5, color = "#555555", hjust = 0, lineheight = 0.95) +
    
    coord_cartesian(xlim = c(-3.5, 5.2), ylim = c(-0.95, 0.95), clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(t = 10, b = 10))
}

# ==============================================================================
# MAIN BUILDER FOR FIGURE 4C (PANEL MAP GRID WITH VERTICAL BANNERS)
# ==============================================================================
build_fig4c_map_grid <- function(df_raw) {
  
  fips_col <- NULL
  for (col in c("fips", "geoid", "county_fips", "region", "abbr")) {
    if (col %in% names(df_raw)) { fips_col <- col; break }
  }
  if (is.null(fips_col)) return(NULL)
  
  df_proc <- df_raw %>%
    mutate(
      abbr = str_trim(get(fips_col)),
      abbr = gsub("^000", "", abbr),
      recent_threshold_x     = as.numeric(recent_threshold_x),
      historical_threshold_x = as.numeric(historical_threshold_x),
      delta_threshold        = recent_threshold_x - historical_threshold_x
    ) %>%
    filter(!abbr %in% non_conus_abbr, !is.na(delta_threshold)) %>%
    group_by(variable) %>%
    mutate(
      iqr_val             = IQR(delta_threshold, na.rm = TRUE),
      iqr_val             = ifelse(iqr_val == 0, 1, iqr_val),
      delta_std           = delta_threshold / iqr_val,
      delta_std_clamped   = pmin(pmax(delta_std, -2), 2)
    ) %>%
    ungroup()
  
  plots_list <- list()
  for (var_key in names(var_config_4c)) {
    df_var <- df_proc %>% filter(variable == var_key)
    plots_list[[var_key]] <- build_single_map(df_var, var_key)
  }
  
  # Vertical Section Banners
  vbanner1 <- build_vertical_banner("Thermal conditions", "#dbebf9")
  vbanner2 <- build_vertical_banner("Temperature extremes", "#fde2e1")
  vbanner3 <- build_vertical_banner("Storm characteristics", "#e2f2e5")
  
  # Row 1: Vertical Banner + 3 maps
  row1_maps <- (plots_list[["mean_T2C"]] | plots_list[["monthly_CDD"]] | plots_list[["monthly_HDD"]])
  row1 <- (vbanner1 | row1_maps) + plot_layout(widths = c(0.04, 1))
  
  # Row 2: Vertical Banner + 2 maps centered
  row2_maps <- (plot_spacer() | plots_list[["heat_days"]] | plots_list[["cold_days"]] | plot_spacer()) + 
    plot_layout(widths = c(0.5, 1, 1, 0.5))
  row2 <- (vbanner2 | row2_maps) + plot_layout(widths = c(0.04, 1))
  
  # Row 3: Vertical Banner + 3 maps
  row3_maps <- (plots_list[["monthly_max_mag"]] | plots_list[["storm_days"]] | plots_list[["monthly_event_count"]])
  row3 <- (vbanner3 | row3_maps) + plot_layout(widths = c(0.04, 1))
  
  legend_plot <- build_fig4c_legend()
  
  # Full Assembly
  full_grid <- (
    row1 / 
      row2 / 
      row3 / 
      legend_plot
  ) + 
    plot_layout(heights = c(3, 3, 3, 0.8))
  
  return(full_grid)
}

# ==============================================================================
# EXECUTE & EXPORT FIGURE 4C (29 x 20 INCHES)
# ==============================================================================
run_fig4c_export <- function(metric_name = "monthly_coh", model_type = "RF") {
  
  fig3b_folder <- file.path(input_base_dir, paste0("Figure3b_", metric_name, "_", model_type))
  csv_input_path <- file.path(fig3b_folder, "All_Regions_TwoPeriod_SafetyWindow_Thresholds.csv")
  
  if (!file.exists(csv_input_path)) {
    warning("Input CSV not found: ", csv_input_path)
    return(NULL)
  }
  
  df_raw <- read.csv(csv_input_path, stringsAsFactors = FALSE, colClasses = "character")
  names(df_raw) <- tolower(names(df_raw))
  
  p_4c <- build_fig4c_map_grid(df_raw)
  
  if (!is.null(p_4c)) {
    png_path_4c <- file.path(output_dir, paste0("Fig4c_", metric_name, "_", model_type, "_MapGrid.png"))
    ggsave(
      png_path_4c,
      plot = p_4c,
      width = 29,
      height = 19,
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    message(">> Exported updated Figure 4c to: ", png_path_4c)
  }
}

# Run execution
run_fig4c_export(metric_name = "monthly_outage_hours", model_type = "RF")
}
