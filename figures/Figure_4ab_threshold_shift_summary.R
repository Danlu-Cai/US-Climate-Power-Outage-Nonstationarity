# =============================================================================
# Script: Figure_4ab_threshold_shift_summary.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description: Figure 4a–4b: state-level threshold-shift proportions and standardized distributions.
# Input: Figure 3b ALE zero-crossing threshold summary CSV files.
# Dependencies: dplyr, tidyr, ggplot2, stringr
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Match the directory structure used by Figure 3b.
fig3_output_dir <- Sys.getenv("FIG3_OUTPUT_DIR", unset = "results/figure3")
output_dir <- Sys.getenv("FIG4_OUTPUT_DIR", unset = "results/figure4")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

non_conus_abbr <- c("AK", "HI", "PR", "VI", "GU", "MP", "AS")

# ==============================================================================
# UNIFIED VARIABLE MAPPINGS & PALETTES
# ==============================================================================
var_label_map_4a <- c(
  "mean_T2C"             = "Mean temperature\n(°C)",
  "monthly_CDD"          = "CDD (base 18 °C)",
  "monthly_HDD"          = "HDD (base 18 °C)",
  "heat_days"            = "Extreme heat days\n(days month⁻¹)",
  "cold_days"            = "Extreme cold days\n(days month⁻¹)",
  "monthly_max_mag"      = "Maximum magnitude\n(severity)",
  "storm_days"           = "Storm duration\n(days month⁻¹)",
  "monthly_event_count" = "Event frequency\n(events month⁻¹)"
)

var_order_top_to_bottom_4a <- c(
  "Mean temperature\n(°C)",
  "CDD (base 18 °C)",
  "HDD (base 18 °C)",
  "Extreme heat days\n(days month⁻¹)",
  "Extreme cold days\n(days month⁻¹)",
  "Maximum magnitude\n(severity)",
  "Storm duration\n(days month⁻¹)",
  "Event frequency\n(events month⁻¹)"
)

unit_map_4a <- c(
  "Mean temperature\n(°C)"             = "°C",
  "CDD (base 18 °C)"                   = "°C-days",
  "HDD (base 18 °C)"                   = "°C-days",
  "Extreme heat days\n(days month⁻¹)" = "days month⁻¹",
  "Extreme cold days\n(days month⁻¹)" = "days month⁻¹",
  "Maximum magnitude\n(severity)"     = "severity units",
  "Storm duration\n(days month⁻¹)"    = "days month⁻¹",
  "Event frequency\n(events month⁻¹)" = "events month⁻¹"
)

cat_levels_4a <- c(
  "Lower threshold",
  "No change",
  "Higher threshold"
)

cat_colors_4a <- c(
  "Lower threshold"  = "#e33935",
  "No change"        = "#cccccc",
  "Higher threshold" = "#1d70b8"
)

point_colors_4b <- c(
  "positive" = "#1d70b8",
  "negative" = "#e33935"
)

format_shift_text <- function(med, q25, q75, unit) {
  if (is.na(med)) return("")
  if (abs(med) >= 10) {
    med_str <- sprintf("%+d", round(med))
    q25_str <- sprintf("%d", round(q25))
    q75_str <- sprintf("%d", round(q75))
  } else {
    med_str <- sprintf("%+.1f", med)
    q25_str <- sprintf("%.1f", q25)
    q75_str <- sprintf("%.1f", q75)
  }
  paste0(med_str, " [", q25_str, ", ", q75_str, "] ", unit)
}

# ==============================================================================
# BUILD FIGURE 4A (STACKED BAR CHART)
# ==============================================================================
build_fig4_barchart <- function(df_raw, eps = 0.05) {
  
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
    mutate(
      category = case_when(
        delta_threshold < -eps ~ "Lower threshold",
        delta_threshold > eps  ~ "Higher threshold",
        TRUE                   ~ "No change"
      ),
      variable_clean = var_label_map_4a[variable]
    ) %>%
    filter(!is.na(variable_clean))
  
  # 1. Stacked percentages
  df_counts <- df_proc %>%
    group_by(variable_clean, category) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(
      variable_clean = unname(var_label_map_4a),
      category = cat_levels_4a,
      fill = list(count = 0)
    ) %>%
    group_by(variable_clean) %>%
    mutate(
      total = sum(count),
      pct = ifelse(total > 0, (count / total) * 100, 0),
      pct_round = round(pct),
      pct_label = ifelse(pct_round > 0, paste0(pct_round, "%"), "")
    ) %>%
    ungroup() %>%
    mutate(
      variable_clean = factor(variable_clean, levels = rev(var_order_top_to_bottom_4a)),
      category       = factor(category, levels = cat_levels_4a)
    )
  
  # 2. Right column shift metrics
  df_shift <- df_proc %>%
    group_by(variable_clean) %>%
    summarise(
      med = median(delta_threshold, na.rm = TRUE),
      q25 = quantile(delta_threshold, 0.25, na.rm = TRUE),
      q75 = quantile(delta_threshold, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      unit       = unit_map_4a[variable_clean],
      shift_text = mapply(format_shift_text, med, q25, q75, unit),
      variable_clean = factor(variable_clean, levels = rev(var_order_top_to_bottom_4a))
    )
  
  p_bar <- ggplot(df_counts, aes(x = pct, y = variable_clean, fill = category)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.72, color = "white", linewidth = 0.6) +
    
    # Increased percentage label size inside bars
    geom_text(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      color = ifelse(df_counts$category == "No change", "#222222", "#ffffff"),
      size = 8.0,
      fontface = "bold"
    ) +
    
    # Increased shift text label size on the right side
    geom_text(
      data = df_shift,
      aes(x = 104, y = variable_clean, label = shift_text),
      inherit.aes = FALSE,
      size = 7.2,
      fontface = "bold",
      color = "#111111",
      hjust = 0
    ) +
    
    # Increased shift table header annotation size
    annotate(
      "text",
      x = 125,
      y = length(var_order_top_to_bottom_4a) + 0.70,
      label = "Median shift [IQR]",
      fontface = "bold",
      size = 7.8,
      color = "#111111",
      hjust = 0.5,
      vjust = 0,
      lineheight = 0.95
    ) +
    
    scale_fill_manual(values = cat_colors_4a) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      limits = c(0, 155),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = "Percentage of states (%)", y = NULL) +
    theme_classic(base_size = 26) +
    theme(
      legend.position      = "top",
      legend.location      = "plot",
      legend.direction     = "horizontal",
      legend.title         = element_blank(),
      legend.text          = element_text(size = 23, face = "bold", color = "#111111"),
      legend.key.size      = unit(0.9, "cm"),
      legend.spacing.x     = unit(0.4, "cm"),
      legend.margin        = margin(b = 15),
      
      axis.text.y          = element_text(size = 23, face = "bold", color = "#111111", hjust = 1, lineheight = 0.85),
      axis.text.x          = element_text(size = 22, face = "bold", color = "#111111"),
      axis.title.x         = element_text(size = 25, face = "bold", color = "#111111", margin = margin(t = 14)),
      
      axis.line.y          = element_blank(),
      axis.ticks.y         = element_blank(),
      axis.line.x          = element_line(color = "#333333", linewidth = 1.4),
      axis.ticks.x         = element_line(color = "#333333", linewidth = 1.4),
      plot.margin          = margin(t = 55, r = 45, b = 20, l = 20)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  return(p_bar)
}

# ==============================================================================
# BUILD FIGURE 4B (REVISED STANDARDIZED DOT PLOT)
# ==============================================================================
build_panel_4b_dotplot <- function(df_raw) {
  
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
    mutate(variable_clean = var_label_map_4a[variable]) %>%
    filter(!is.na(variable_clean))
  
  df_summary <- df_proc %>%
    group_by(variable_clean) %>%
    summarise(
      iqr_val   = IQR(delta_threshold, na.rm = TRUE),
      iqr_val   = ifelse(iqr_val == 0, 1, iqr_val),
      median_Δ  = median(delta_threshold, na.rm = TRUE) / iqr_val,
      q25_Δ     = quantile(delta_threshold, 0.25, na.rm = TRUE) / iqr_val,
      q75_Δ     = quantile(delta_threshold, 0.75, na.rm = TRUE) / iqr_val,
      q05_Δ     = quantile(delta_threshold, 0.05, na.rm = TRUE) / iqr_val,
      q95_Δ     = quantile(delta_threshold, 0.95, na.rm = TRUE) / iqr_val,
      .groups   = "drop"
    ) %>%
    mutate(
      color_cat      = ifelse(median_Δ >= 0, "positive", "negative"),
      variable_clean = factor(variable_clean, levels = rev(var_order_top_to_bottom_4a))
    )
  
  y_top <- length(var_order_top_to_bottom_4a) + 0.95
  
  p_4b <- ggplot(df_summary, aes(y = variable_clean, color = color_cat)) +
    # Vertical reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "#a0a0a0", linewidth = 1.0) +
    
    # 1. Outer 5th-95th percentile errorbar
    geom_errorbarh(
      aes(xmin = q05_Δ, xmax = q95_Δ),
      height = 0.20,
      linewidth = 0.8,
      alpha = 0.45
    ) +
    
    # 2. Inner IQR (25th-75th) errorbar
    geom_errorbarh(
      aes(xmin = q25_Δ, xmax = q75_Δ),
      height = 0.28,
      linewidth = 1.5
    ) +
    
    # 3. Median point
    geom_point(aes(x = median_Δ), size = 6.0) +
    
    # ==========================================================================
  # TOP LEGEND ANNOTATIONS (BIGGER LABELS)
  # ==========================================================================
  # --- Item 1: Median ---
  annotate("point", x = -1.55, y = y_top, size = 6.0, color = "#1d70b8") +
    annotate("text",  x = -1.38, y = y_top, label = "Median", fontface = "bold", size = 7.2, color = "#222222", hjust = 0) +
    
    # --- Item 2: IQR (25th–75th) ---
    annotate("segment", x = -0.55, xend = -0.15, y = y_top, yend = y_top, color = "#1d70b8", linewidth = 1.5) +
    annotate("segment", x = -0.55, xend = -0.55, y = y_top - 0.12, yend = y_top + 0.12, color = "#1d70b8", linewidth = 1.5) +
    annotate("segment", x = -0.15, xend = -0.15, y = y_top - 0.12, yend = y_top + 0.12, color = "#1d70b8", linewidth = 1.5) +
    annotate("text",    x = -0.02, y = y_top, label = "IQR (25th–75th)", fontface = "bold", size = 7.2, color = "#222222", hjust = 0) +
    
    # --- Item 3: 5th–95th ---
    annotate("segment", x =  1.05, xend =  1.45, y = y_top, yend = y_top, color = "#1d70b8", linewidth = 0.8, alpha = 0.45) +
    annotate("segment", x =  1.05, xend =  1.05, y = y_top - 0.10, yend = y_top + 0.10, color = "#1d70b8", linewidth = 0.8, alpha = 0.45) +
    annotate("segment", x =  1.45, xend =  1.45, y = y_top - 0.10, yend = y_top + 0.10, color = "#1d70b8", linewidth = 0.8, alpha = 0.45) +
    annotate("text",    x =  1.55, y = y_top, label = "5th–95th", fontface = "bold", size = 7.2, color = "#222222", hjust = 0) +
    
    scale_color_manual(values = point_colors_4b) +
    scale_x_continuous(
      limits = c(-2.2, 2.2),
      breaks = seq(-2, 2, 1)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = expression(bold("Standardized threshold shift (" * Delta * " / IQR of shifts)")),
      y = NULL
    ) +
    theme_classic(base_size = 26) +
    theme(
      legend.position      = "none",
      
      panel.grid.major.x   = element_line(color = "#e5e5e5", linewidth = 0.6),
      panel.grid.minor     = element_blank(),
      panel.grid.major.y   = element_blank(),
      
      axis.line.x          = element_line(color = "#333333", linewidth = 1.4),
      axis.line.y          = element_line(color = "#333333", linewidth = 1.4),
      axis.ticks.x         = element_line(color = "#333333", linewidth = 1.4),
      axis.ticks.y         = element_line(color = "#333333", linewidth = 1.4),
      axis.text.y          = element_text(size = 23, face = "bold", color = "#111111", hjust = 1, lineheight = 0.85),
      axis.text.x          = element_text(size = 22, face = "bold", color = "#111111"),
      axis.title.x         = element_text(size = 25, face = "bold", color = "#111111", margin = margin(t = 14)),
      
      plot.margin          = margin(t = 60, r = 35, b = 20, l = 20)
    )
  
  return(p_4b)
}

# ==============================================================================
# EXECUTE & EXPORT
# ==============================================================================
run_fig4_export_all <- function(metric_name = "monthly_coh", model_type = "RF", eps = 0.05) {
  
  fig3b_folder <- file.path(fig3_output_dir, paste0("Figure3b_", metric_name, "_", model_type))
  csv_input_path <- file.path(fig3b_folder, "All_Regions_TwoPeriod_SafetyWindow_Thresholds.csv")
  
  if (!file.exists(csv_input_path)) {
    warning("Input CSV not found: ", csv_input_path)
    return(NULL)
  }
  
  df_raw <- read.csv(csv_input_path, stringsAsFactors = FALSE, colClasses = "character")
  names(df_raw) <- tolower(names(df_raw))
  if (!"region_type" %in% names(df_raw)) {
    stop("Figure 3b threshold CSV must contain Region_Type to identify state-level rows.")
  }
  df_raw <- df_raw[tolower(trimws(df_raw$region_type)) %in% c("statewise", "state"), , drop = FALSE]
  if (nrow(df_raw) == 0) {
    warning("No state-level thresholds available in: ", csv_input_path)
    return(NULL)
  }
  
  # 1. Export Figure 4a
  p_4a <- build_fig4_barchart(df_raw, eps = eps)
  if (!is.null(p_4a)) {
    png_path_4a <- file.path(output_dir, paste0("Fig4a_", metric_name, "_", model_type, "_Bar_AllVars.png"))
    ggsave(
      png_path_4a,
      plot = p_4a,
      width = 15.0,
      height = 8.5,
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    message(">> Exported Figure 4a (Bar Chart) to: ", png_path_4a)
  }
  
  # 2. Export Figure 4b
  p_4b <- build_panel_4b_dotplot(df_raw)
  if (!is.null(p_4b)) {
    png_path_4b <- file.path(output_dir, paste0("Fig4b_", metric_name, "_", model_type, "_DotPlot_AllVars.png"))
    ggsave(
      png_path_4b,
      plot = p_4b,
      width = 15.0,
      height = 8.5,
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    message(">> Exported Figure 4b (Dot Plot) to: ", png_path_4b)
  }
}

# Run execution
run_fig4_export_all(metric_name = "monthly_coh",          model_type = "XGBoost", eps = 0.05)
run_fig4_export_all(metric_name = "monthly_coh",          model_type = "RF",      eps = 0.05)
run_fig4_export_all(metric_name = "monthly_outage_hours", model_type = "XGBoost", eps = 0.05)
run_fig4_export_all(metric_name = "monthly_outage_hours", model_type = "RF",      eps = 0.05)