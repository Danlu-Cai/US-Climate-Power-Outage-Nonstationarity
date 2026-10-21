# =============================================================================
# Script: Figure_3c_ALE_threshold_heatmap.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description: Heatmap of descriptive ALE zero-crossing shifts using Figure 3b outputs.
# Dependencies: dplyr, tidyr, ggplot2, stringr, cowplot, grid, scales
# =============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(grid)
library(scales)

# ==============================================================================
# CONFIGURATION & PATHS
# ==============================================================================
input_base_dir <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related"
output_dir     <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related_final"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

non_conus_abbr <- c("AK", "HI", "PR", "VI", "GU", "MP", "AS")

# Variable order & clean display labels
var_label_map <- c(
  "mean_T2C"             = "Mean temperature (°C)",
  "monthly_CDD"          = "CDD (base 18 °C)",
  "monthly_HDD"          = "HDD (base 18 °C)",
  "heat_days"            = "Extreme heat days (days month⁻¹)",
  "cold_days"            = "Extreme cold days (days month⁻¹)",
  "monthly_max_mag"      = "Maximum magnitude (severity)",
  "storm_days"           = "Storm days (days month⁻¹)",
  "monthly_event_count" = "Event frequency (events month⁻¹)"
)

var_order <- c(
  "Mean temperature (°C)",
  "CDD (base 18 °C)",
  "HDD (base 18 °C)",
  "Extreme heat days (days month⁻¹)",
  "Extreme cold days (days month⁻¹)",
  "Maximum magnitude (severity)",
  "Storm days (days month⁻¹)",
  "Event frequency (events month⁻¹)"
)

# Region display mappings matching Fig 3a exactly
region_name_map <- c(
  "West"                        = "West",
  "Northwest"                   = "Northwest",
  "South"                       = "South",
  "Central"                     = "Ohio\nValley",
  "Ohio Valley"                 = "Ohio\nValley",
  "East North Central"          = "Upper\nMidwest",
  "Upper Midwest"               = "Upper\nMidwest",
  "Southeast"                   = "Southeast",
  "Northeast"                   = "Northeast",
  "West North Central"          = "Northern\nRockies",
  "Northern Rockies"            = "Northern\nRockies",
  "Northern Rockies and Plains" = "Northern\nRockies",
  "Southwest"                   = "Southwest",
  "All CONUS"                   = "All\nCONUS",
  "CONUS"                       = "All\nCONUS"
)

# Region display order aligned with Fig 3a
region_order <- c(
  "All\nCONUS", "Northern\nRockies", "South", "Northeast", 
  "Northwest", "Upper\nMidwest", "Southwest", "Ohio\nValley", 
  "Southeast", "West"
)

# ==============================================================================
# HELPER: ROBUST LEGEND EXTRACTION
# ==============================================================================
get_heatmap_legend <- function(p) {
  leg <- tryCatch({
    cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(leg)) {
    leg <- tryCatch({
      cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
    }, error = function(e) NULL)
  }
  
  if (is.null(leg)) {
    leg <- cowplot::get_legend(p)
  }
  return(leg)
}

# ==============================================================================
# DATA PROCESSING FUNCTION (WITH CONUS DEDUPLICATION)
# ==============================================================================
process_fig3c_data <- function(csv_path) {
  if (!file.exists(csv_path)) {
    warning("Target CSV not found at path: ", csv_path)
    return(NULL)
  }
  
  df_raw <- read.csv(csv_path, stringsAsFactors = FALSE)
  names(df_raw) <- tolower(names(df_raw))
  
  fips_col <- NULL
  for (col in c("fips", "geoid", "county_fips", "region", "abbr")) {
    if (col %in% names(df_raw)) { fips_col <- col; break }
  }
  
  df_proc <- df_raw %>%
    mutate(
      abbr = if (!is.null(fips_col)) str_trim(get(fips_col)) else NA_character_,
      abbr = gsub("^000", "", abbr),
      recent_threshold_x     = as.numeric(recent_threshold_x),
      historical_threshold_x = as.numeric(historical_threshold_x),
      delta_threshold        = recent_threshold_x - historical_threshold_x,
      variable_clean         = var_label_map[variable]
    ) %>%
    filter(!abbr %in% non_conus_abbr, !is.na(variable_clean))
  
  # Baseline for IQR estimation
  df_state_level <- if ("region_type" %in% names(df_proc) && any(tolower(df_proc$region_type) == "state")) {
    df_proc %>% filter(tolower(region_type) == "state")
  } else {
    df_proc
  }
  
  df_iqr <- df_state_level %>%
    group_by(variable_clean) %>%
    summarise(
      iqr_hist = IQR(historical_threshold_x, na.rm = TRUE),
      iqr_hist = ifelse(is.na(iqr_hist) | iqr_hist < 0.05, 1.0, iqr_hist),
      .groups  = "drop"
    )
  
  # Regional summaries (excluding raw CONUS rows to prevent duplicate text overlay)
  df_regional <- df_proc %>%
    filter(
      (tolower(region_type) == "ncei region" | region %in% names(region_name_map)),
      !region %in% c("All CONUS", "CONUS", "National", "All")
    ) %>%
    mutate(region_disp = recode(region, !!!region_name_map)) %>%
    group_by(region_disp, variable_clean) %>%
    summarise(
      med_hist = median(historical_threshold_x, na.rm = TRUE),
      med_rec  = median(recent_threshold_x, na.rm = TRUE),
      delta    = med_rec - med_hist,
      .groups  = "drop"
    )
  
  # CONUS Summary
  df_conus <- df_proc %>%
    group_by(variable_clean) %>%
    summarise(
      med_hist = median(historical_threshold_x, na.rm = TRUE),
      med_rec  = median(recent_threshold_x, na.rm = TRUE),
      delta    = med_rec - med_hist,
      .groups  = "drop"
    ) %>%
    mutate(region_disp = "All\nCONUS")
  
  # Combine & Standardize
  df_combined <- bind_rows(df_regional, df_conus) %>%
    distinct(region_disp, variable_clean, .keep_all = TRUE) %>%
    left_join(df_iqr, by = "variable_clean") %>%
    mutate(
      std_shift = delta / iqr_hist,
      val_label = ifelse(is.na(std_shift), "", sprintf("%.1f", std_shift)),
      variable_clean = factor(variable_clean, levels = rev(var_order)),
      region_disp    = factor(region_disp, levels = region_order)
    ) %>%
    complete(variable_clean, region_disp)
  
  return(df_combined)
}

# ==============================================================================
# PLOT RENDERER (ORIGINAL COLORS WITH TRANSPARENCY + LARGE LEGEND LABELS)
# ==============================================================================
build_fig3c_heatmap <- function(df_data) {
  if (is.null(df_data) || nrow(df_data) == 0) return(NULL)
  
  # Original colors softened with 70% opacity (alpha = 0.70)
  orig_red_alpha  <- scales::alpha("#d73027", 0.70)
  orig_mid_alpha  <- "#ffffff"
  orig_blue_alpha <- scales::alpha("#1a73e8", 0.70)
  orig_gray_alpha <- scales::alpha("#b0b5b9", 0.70)
  
  # 1. Main Heatmap Plot Panel
  p_main <- ggplot(df_data, aes(x = region_disp, y = variable_clean)) +
    geom_tile(aes(fill = std_shift), color = "white", linewidth = 1.2) +
    geom_text(aes(label = val_label), color = "#111111", size = 5.8, fontface = "bold") +
    scale_fill_gradient2(
      low      = orig_red_alpha,
      mid      = orig_mid_alpha,
      high     = orig_blue_alpha,
      midpoint = 0,
      limits   = c(-2, 2),
      oob      = scales::squish,
      na.value = orig_gray_alpha
    ) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x.top   = element_text(size = 14, face = "bold", color = "#111111", lineheight = 0.95, margin = margin(b = 8)),
      axis.text.y       = element_text(size = 14, face = "bold", color = "#111111", hjust = 1),
      panel.grid        = element_blank(),
      legend.position   = "none",
      plot.margin       = margin(t = 15, r = 20, b = 10, l = 20)
    )
  
  # 2. Standalone Legend Source (Large Labels & Ticks)
  p_legend_source <- ggplot(df_data, aes(x = region_disp, y = variable_clean, fill = std_shift)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = orig_red_alpha,
      mid      = orig_mid_alpha,
      high     = orig_blue_alpha,
      midpoint = 0,
      limits   = c(-2, 2),
      oob      = scales::squish,
      name     = "Standardized threshold shift (Δ / historical IQR)"
    ) +
    theme_minimal() + 
    theme(
      legend.position   = "bottom",
      legend.key.width  = unit(4.2, "cm"),   # Wide colorbar for spacious labels
      legend.key.height = unit(0.50, "cm"),  # Clear gradient bar
      legend.title      = element_text(size = 15, face = "bold", color = "#111111", margin = margin(t = 12)), # 15 pt bold title
      legend.text       = element_text(size = 15, face = "bold", color = "#111111", margin = margin(t = 5))    # 15 pt bold ticks
    ) +
    guides(
      fill = guide_colorbar(
        title.position  = "bottom",
        title.hjust     = 0.5,
        ticks.colour    = "black",
        ticks.linewidth = 1.2,
        frame.colour    = "#333333",
        frame.linewidth = 1.0
      )
    )
  
  legend_grob <- get_heatmap_legend(p_legend_source)
  
  # 3. Canvas Layout Assembly
  p_final <- ggdraw() +
    # Heatmap Grid Plot
    draw_plot(p_main, x = 0, y = 0.19, width = 1.0, height = 0.81) +
    
    # Center Colorbar Legend
    draw_plot(legend_grob, x = 0.20, y = -0.01, width = 0.45, height = 0.21) +
    
    # Left Red Annotation (14 pt bold, original red hue)
    draw_label("Lower threshold\n(in recent period)", 
               x = 0.21, y = 0.08, 
               fontface = "bold", color = "#d73027", size = 14, hjust = 1) +
    
    # Right Blue Annotation (14 pt bold, original blue hue)
    draw_label("Higher threshold\n(in recent period)", 
               x = 0.64, y = 0.08, 
               fontface = "bold", color = "#1a73e8", size = 14, hjust = 0) +
    
    # Far-Right Grey Box & "Not estimable" Legend (12.5 pt)
    draw_grob(
      grid::rectGrob(
        x = 0.83, y = 0.08, width = 0.022, height = 0.038,
        gp = grid::gpar(fill = orig_gray_alpha, col = "white", lwd = 1)
      )
    ) +
    draw_label("Not estimable\n(insufficient data or no threshold)", 
               x = 0.85, y = 0.08, 
               fontface = "plain", color = "#222222", size = 12.5, hjust = 0)
  
  return(p_final)
}

# ==============================================================================
# EXECUTION RUNNER
# ==============================================================================
run_fig3c_export <- function(metric_name = "monthly_coh", model_type = "RF") {
  fig3b_folder <- file.path(input_base_dir, paste0("Figure3b_", metric_name, "_", model_type))
  csv_input_path <- file.path(fig3b_folder, "All_Regions_TwoPeriod_SafetyWindow_Thresholds.csv")
  
  if (!file.exists(csv_input_path)) {
    warning("Input CSV not found: ", csv_input_path)
    return(NULL)
  }
  
  df_summary <- process_fig3c_data(csv_input_path)
  
  if (!is.null(df_summary)) {
    # Export CSV summary table
    csv_out_path <- file.path(output_dir, paste0("Fig3c_", metric_name, "_", model_type, "_Heatmap_Data.csv"))
    write.csv(df_summary, csv_out_path, row.names = FALSE)
    
    # Export Heatmap Plot
    p_fig3c <- build_fig3c_heatmap(df_summary)
    png_out_path <- file.path(output_dir, paste0("Fig3c_", metric_name, "_", model_type, "_Heatmap.png"))
    
    ggsave(
      png_out_path,
      plot   = p_fig3c,
      width  = 16.0,
      height = 7.5,
      dpi    = 300,
      bg     = "white",
      limitsize = FALSE
    )
    message(">> Successfully exported Figure 3c Heatmap to: ", png_out_path)
  }
}

# Execute for all 4 model runs
run_fig3c_export(metric_name = "monthly_coh",          model_type = "XGBoost")
run_fig3c_export(metric_name = "monthly_coh",          model_type = "RF")
run_fig3c_export(metric_name = "monthly_outage_hours", model_type = "XGBoost")
run_fig3c_export(metric_name = "monthly_outage_hours", model_type = "RF")

# Execute for all 4 model runs
run_fig3c_export(metric_name = "monthly_coh",          model_type = "XGBoost")
run_fig3c_export(metric_name = "monthly_coh",          model_type = "RF")
run_fig3c_export(metric_name = "monthly_outage_hours", model_type = "XGBoost")
run_fig3c_export(metric_name = "monthly_outage_hours", model_type = "RF")
