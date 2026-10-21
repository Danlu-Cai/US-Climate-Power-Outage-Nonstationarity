# =============================================================================
# Script: Figure_3b_climate_driver_ALE.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description: Bootstrap ALE curves and descriptive zero-crossing thresholds.
# Dependencies: sf, xgboost, ranger, dplyr, ggplot2, cowplot, stringr, tidyr, iml, Matrix, scales
# =============================================================================

# ==============================================================================
# PROJECT: DigiCARES - Publication-Grade ALE Curves & Safety Window Thresholds
# SCALE:   Monthly Analytics (CONUS + 9 NCEI Regions + Statewise)
# OUTPUT:  Figure3b_MultiVariable_Analysis (8 Climate Drivers Grid Plots & CSV)
# ==============================================================================

# --- 1. Setup & Libraries ---
library(sf)
library(xgboost)
library(ranger)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(tidyr)
library(iml)
library(Matrix)
library(scales)

# --- 2. Configuration & Metadata ---
input_path <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related/Digicare_Merged_Monthly_2015_2022_CONUS.csv" 
output_dir <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related_final"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

period_colors <- c(
  "Historical (2015–2018)" = "#00468B", 
  "Recent (2019–2022)"     = "#ED0000"
)

# 4 Key Climate Variables matching the 2x2 panel figure layout
variable_metadata <- list(
  "cold_days"           = list(title = "Extreme cold duration",   unit = "days/month"),
  "monthly_event_count" = list(title = "Extreme event frequency", unit = "events/month"),
  "monthly_max_mag"     = list(title = "Maximum event magnitude", unit = "severity index"),
  "storm_days"          = list(title = "Storm duration",            unit = "days/month")
)

climate_vars <- names(variable_metadata)

# --- 3. Helper Functions ---
get_rising_zero_crossing <- function(ale_df) {
  df_sorted <- ale_df %>% arrange(x_val) %>% filter(!is.na(x_val), !is.na(.value))
  if (nrow(df_sorted) < 2) return(NA_real_)
  
  x_vals <- df_sorted$x_val
  y_vals <- df_sorted$.value
  rising_idx <- which(y_vals[-length(y_vals)] <= 0 & y_vals[-1] > 0)
  
  if (length(rising_idx) == 0) {
    if (min(y_vals) > 0) return(min(x_vals))
    return(NA_real_)
  }
  
  idx <- tail(rising_idx, 1)
  x1 <- x_vals[idx];     x2 <- x_vals[idx + 1]
  y1 <- y_vals[idx];     y2 <- y_vals[idx + 1]
  
  return(x1 + (0 - y1) * (x2 - x1) / (y2 - y1))
}

prepare_nexus_data_log <- function(df, target_col) {
  df_clean <- df %>%
    mutate(target = log10(as.numeric(.data[[target_col]]) + 1)) %>% 
    filter(!is.infinite(target), !is.na(target)) %>%
    mutate(
      parsed_date = as.Date(month_stamp, tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%m/%d/%Y %H:%M", "%m/%d/%Y")),
      Year = as.numeric(format(parsed_date, "%Y")),
      Period = case_when(
        Year >= 2015 & Year <= 2018 ~ "Historical (2015–2018)",
        Year >= 2019 & Year <= 2022 ~ "Recent (2019–2022)",
        TRUE ~ NA_character_
      ),
      Month = factor(format(parsed_date, "%b"), levels = month.abb)
    ) %>%
    filter(!is.na(Month), !is.na(Period)) %>%
    select(
      target, heat_days, cold_days, mean_T2C, monthly_CDD, monthly_HDD, 
      monthly_event_count, monthly_max_mag, storm_days, pop_density, 
      Period, Month, month_stamp                                                  
    ) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  return(df_clean)
}

get_bootstrapped_ale <- function(p_data, var_name, current_period, current_region, model_type = "RF", n_boot = 30, grid_size = 50) {
  x_vals_orig <- p_data[[var_name]]
  x_min <- min(x_vals_orig, na.rm = TRUE)
  x_max <- max(x_vals_orig, na.rm = TRUE)
  
  if (x_min == x_max) return(NULL)
  
  eval_grid <- seq(x_min, x_max, length.out = grid_size)
  boot_matrix <- matrix(NA, nrow = grid_size, ncol = n_boot)
  
  for (b in 1:n_boot) {
    set.seed(b * 100)
    boot_idx <- sample(nrow(p_data), replace = TRUE)
    b_data   <- p_data[boot_idx, ]
    
    y_b <- b_data$target
    X_b <- b_data %>% select(-target, -month_stamp, -Period, -any_of(c("county_area", "area", "State", "NCEI_Region")))
    X_b$Month <- factor(X_b$Month, levels = month.abb)
    
    if (model_type == "XGBoost") {
      X_mat <- model.matrix(~ . - 1, data = X_b)
      mod_b <- xgboost(data = X_mat, label = y_b, nrounds = 100, verbose = 0, nthread = 1)
      pred_b <- Predictor$new(mod_b, data = as.data.frame(X_mat), y = y_b)
    } else {
      mod_b <- ranger(target ~ ., data = cbind(target = y_b, X_b), num.trees = 150)
      pred_b <- Predictor$new(mod_b, data = X_b, y = y_b)
    }
    
    fe_res <- tryCatch({
      FeatureEffect$new(pred_b, feature = var_name, method = "ale", grid.size = grid_size)$results
    }, error = function(e) NULL)
    
    if (!is.null(fe_res) && nrow(fe_res) >= 2) {
      boot_matrix[, b] <- approx(x = fe_res[[var_name]], y = fe_res$.value, xout = eval_grid, rule = 2)$y
    }
  }
  
  res_df <- data.frame(
    x_val    = eval_grid,
    .value   = rowMeans(boot_matrix, na.rm = TRUE),
    lower_ci = apply(boot_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
    upper_ci = apply(boot_matrix, 1, quantile, probs = 0.975, na.rm = TRUE),
    Period   = current_period,
    Variable = var_name,
    Region   = current_region
  )
  return(res_df)
}

# --- 4. Main Execution Runner ---
run_fig3b <- function(metric_name, model_type = "RF", n_boot = 30) {
  
  png_filename <- paste0("Fig3b_", metric_name, "_", model_type, ".png")
  csv_filename <- paste0("Fig3b_", metric_name, "_", model_type, ".csv")
  
  message("\n==========================================================")
  message(sprintf("Executing Figure 3b | Metric: %s | Model: %s", metric_name, model_type))
  message("==========================================================")
  
  heat_data <- read.csv(input_path)
  full_data_log <- prepare_nexus_data_log(heat_data, metric_name)
  
  current_region <- "CONUS"
  reg_type <- "CONUS"
  region_data <- full_data_log
  
  region_ale_accumulator <- list()
  region_thr_accumulator <- list()
  
  for (current_period in c("Historical (2015–2018)", "Recent (2019–2022)")) {
    p_data <- region_data %>% filter(Period == current_period)
    if (nrow(p_data) < 10) next
    
    for (var_name in climate_vars) {
      if (!var_name %in% colnames(p_data)) next
      
      boot_ale <- get_bootstrapped_ale(p_data, var_name, current_period, current_region, model_type = model_type, n_boot = n_boot)
      if (is.null(boot_ale)) next
      
      region_ale_accumulator[[paste(current_period, var_name, sep = "_")]] <- boot_ale
      
      exact_thr <- get_rising_zero_crossing(boot_ale)
      region_thr_accumulator[[paste(current_period, var_name, sep = "_")]] <- data.frame(
        Metric = metric_name,
        Model = model_type,
        Region = current_region,
        Region_Type = reg_type,
        Period = current_period, 
        Variable = var_name, 
        Safety_Window_Threshold_x = exact_thr
      )
    }
  }
  
  combined_region_ale <- bind_rows(region_ale_accumulator)
  combined_region_thr <- bind_rows(region_thr_accumulator)
  
  # --- Generate CONUS 2x2 Grid Plot ---
  if (nrow(combined_region_ale) > 0) {
    plot_list <- list()
    
    for (i in seq_along(climate_vars)) {
      var_name <- climate_vars[i]
      var_ale  <- combined_region_ale %>% filter(Variable == var_name)
      var_thr  <- combined_region_thr %>% filter(Variable == var_name)
      
      if (nrow(var_ale) == 0) next
      
      meta <- variable_metadata[[var_name]]
      hist_thr <- var_thr %>% filter(Period == "Historical (2015–2018)") %>% pull(Safety_Window_Threshold_x)
      rec_thr  <- var_thr %>% filter(Period == "Recent (2019–2022)") %>% pull(Safety_Window_Threshold_x)
      
      hist_str <- if (length(hist_thr) == 1 && !is.na(hist_thr)) sprintf("Historical [Thr = %.1f]", hist_thr) else "Historical"
      rec_str  <- if (length(rec_thr) == 1 && !is.na(rec_thr)) sprintf("Recent [Thr = %.1f]", rec_thr) else "Recent"
      
      legend_mapping <- c(
        "Historical (2015–2018)" = hist_str,
        "Recent (2019–2022)"     = rec_str
      )
      
      x_axis_title <- paste0(meta$title, "\n(", meta$unit, ")")
      
      # Suppress individual plot y-axis titles to use a single unified y-axis label
      p <- ggplot(var_ale, aes(x = x_val, group = Period, color = Period)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.9) +
        geom_line(aes(y = lower_ci), linetype = "dashed", linewidth = 0.9, alpha = 0.85) +
        geom_line(aes(y = upper_ci), linetype = "dashed", linewidth = 0.9, alpha = 0.85) +
        geom_line(aes(y = `.value`), linewidth = 1.8) +
        scale_color_manual(values = period_colors, labels = legend_mapping) +
        guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
        theme_classic(base_size = 16) +
        labs(
          x = x_axis_title,
          y = NULL
        ) +
        theme(
          # --- ENLARGED FONT CONFIGURATIONS ---
          axis.title.x = element_text(face = "bold", size = 16, color = "black", lineheight = 0.9, margin = margin(t = 8)),
          axis.text    = element_text(face = "bold", size = 15, color = "black"),
          
          # Top legend configuration with enlarged text
          legend.position      = "top",
          legend.justification = "left",
          legend.direction     = "vertical",
          legend.title         = element_blank(),
          legend.text          = element_text(size = 15, face = "bold"),
          legend.key.width     = unit(1.2, "cm"),
          legend.key.height    = unit(0.5, "cm"),
          legend.margin        = margin(b = 2),
          legend.box.margin    = margin(b = -4),
          legend.background    = element_rect(fill = "transparent", color = NA),
          
          plot.margin          = margin(t = 6, r = 14, b = 10, l = 12),
          axis.line            = element_line(linewidth = 1.0, color = "black"),
          axis.ticks           = element_line(linewidth = 1.0, color = "black")
        )
      
      if (length(hist_thr) == 1 && !is.na(hist_thr)) {
        p <- p + geom_vline(xintercept = hist_thr, linetype = "dashed", color = period_colors["Historical (2015–2018)"], linewidth = 0.9, alpha = 0.85)
      }
      if (length(rec_thr) == 1 && !is.na(rec_thr)) {
        p <- p + geom_vline(xintercept = rec_thr, linetype = "dashed", color = period_colors["Recent (2019–2022)"], linewidth = 0.9, alpha = 0.85)
      }
      
      plot_list[[var_name]] <- p
    }
    
    if (length(plot_list) > 0) {
      # Combine the 2x2 grid without individual y-axis titles
      grid_plots <- plot_grid(plotlist = plot_list, ncol = 2, nrow = 2, align = "vh")
      
      # Create single shared Y-axis title label on the far left
      y_title_label <- ggdraw() + 
        draw_label(
          "ALE (log cumulative outage hours)", 
          angle = 90, 
          fontface = "bold", 
          size = 18, 
          x = 0.5, 
          y = 0.5
        )
      
      # Attach shared Y-axis title to the left of the 2x2 grid
      final_figure <- plot_grid(y_title_label, grid_plots, rel_widths = c(0.04, 0.96))
      
      full_png_path <- file.path(output_dir, png_filename)
      ggsave(
        full_png_path, 
        final_figure, width = 8.0, height = 8.0, dpi = 300
      )
      message(">> Exported figure directly to: ", full_png_path)
    }
  }
  
  if (nrow(combined_region_thr) > 0) {
    final_threshold_table <- combined_region_thr %>%
      pivot_wider(names_from = Period, values_from = Safety_Window_Threshold_x) %>%
      rename(
        Historical_Threshold_x = `Historical (2015–2018)`,
        Recent_Threshold_x     = `Recent (2019–2022)`
      ) %>%
      mutate(Safety_Window_Shift = Recent_Threshold_x - Historical_Threshold_x)
    
    full_csv_path <- file.path(output_dir, csv_filename)
    write.csv(final_threshold_table, full_csv_path, row.names = FALSE)
    message(">> Exported CSV thresholds directly to: ", full_csv_path)
  }
}

# --- Execute 4 Target Metric & Model Combinations ---
run_fig3b(metric_name = "monthly_coh",          model_type = "XGBoost")
run_fig3b(metric_name = "monthly_coh",          model_type = "RF")
run_fig3b(metric_name = "monthly_outage_hours", model_type = "XGBoost")
run_fig3b(metric_name = "monthly_outage_hours", model_type = "RF")
