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
Input_ROOT <- Sys.getenv("FIG3_INPUT_DIR", unset = "data/processed")
output_root <- Sys.getenv("FIG3_OUTPUT_DIR", unset = "results/figure3")
input_path <- file.path(Input_ROOT, "Digicare_Merged_Monthly_2015_2022_CONUS.csv")

period_colors <- c(
  "Historical (2015–2018)" = "#00468B", 
  "Recent (2019–2022)"     = "#ED0000"
)

variable_metadata <- list(
  "monthly_CDD"         = list(title = "Cooling degree days",        unit = "CDD, base 18 °C"),
  "monthly_HDD"         = list(title = "Heating degree days",        unit = "HDD, base 18 °C"),
  "mean_T2C"            = list(title = "Mean 2-m temperature",       unit = "°C"),
  "heat_days"           = list(title = "Extreme heat frequency",      unit = "days/month"),
  "cold_days"           = list(title = "Extreme cold frequency",      unit = "days/month"),
  "monthly_event_count" = list(title = "Extreme event frequency",    unit = "events/month"), #could be several times per day
  "monthly_max_mag"     = list(title = "Maximum event magnitude",    unit = "severity index"),
  "storm_days"          = list(title = "Storm days",             unit = "days/month")
)

climate_vars <- names(variable_metadata)

# --- 3. Regional Lookup Mapping (CONUS + 9 NCEI Regions + 50 States) ---
ncei_lookup <- data.frame(
  Abbr = c("IL", "IN", "KY", "MO", "OH", "WV",
           "IA", "MI", "MN", "WI",
           "CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT",
           "ID", "OR", "WA", 
           "AR", "LA", "MS", "OK", "TX",
           "AL", "FL", "GA", "NC", "SC", "VA",
           "AZ", "CO", "NM", "UT",
           "CA", "NV",
           "MT", "NE", "ND", "SD", "WY"),
  NCEI_Region = c(rep("Central", 6),
                  rep("East North Central", 4),
                  rep("Northeast", 11),
                  rep("Northwest", 3),
                  rep("South", 5),
                  rep("Southeast", 6),
                  rep("Southwest", 4),
                  rep("West", 2),
                  rep("West North Central", 5))
)

# --- 4. Helper Functions ---
get_rising_zero_crossing <- function(ale_df) {
  df_sorted <- ale_df %>% arrange(x_val) %>% filter(!is.na(x_val), !is.na(.value))
  if (nrow(df_sorted) < 2) return(NA_real_)
  
  x_vals <- df_sorted$x_val
  y_vals <- df_sorted$.value
  rising_idx <- which(y_vals[-length(y_vals)] <= 0 & y_vals[-1] > 0)
  
  if (length(rising_idx) == 0) {
    if (min(y_vals) > 0) return(NA_real_)
    return(NA_real_)
  }
  
  idx <- tail(rising_idx, 1)
  x1 <- x_vals[idx];     x2 <- x_vals[idx + 1]
  y1 <- y_vals[idx];     y2 <- y_vals[idx + 1]
  
  return(x1 + (0 - y1) * (x2 - x1) / (y2 - y1))
}

prepare_nexus_data_log <- function(df, target_col) {
  state_lookup <- data.frame(
    code = sprintf("%02d", c(1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56)),
    StateName = c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming")
  )
  
  df_clean <- df %>%
    mutate(state_code_clean = str_pad(substr(str_pad(fips, 5, pad = "0"), 1, 2), 2, pad = "0")) %>%
    left_join(state_lookup, by = c("state_code_clean" = "code")) %>%
    mutate(Abbr = state.abb[match(StateName, state.name)]) %>%
    left_join(ncei_lookup, by = "Abbr") %>%
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
      State = as.character(Abbr), 
      Month = factor(format(parsed_date, "%b"), levels = month.abb)
    ) %>%
    filter(!is.na(Month), !is.na(State), !is.na(Period)) %>%
    select(
      target, heat_days, cold_days, mean_T2C, monthly_CDD, monthly_HDD, 
      monthly_event_count, monthly_max_mag, storm_days, pop_density, 
      State, NCEI_Region, Period, Month, month_stamp                                                           
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
    X_b <- b_data %>% select(-target, -month_stamp, -State, -NCEI_Region, -Period, -any_of(c("county_area", "area")))
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
  
  if (all(is.na(boot_matrix))) return(NULL)
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

# --- 5. Main Execution Runner ---
run_fig3b <- function(metric_name, model_type = "RF", n_boot = 30) {
  
  save_path <- file.path(output_root, paste0("Figure3b_", metric_name, "_", model_type))
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  message("\n==========================================================")
  message(sprintf("Executing Figure 3b | Target: %s | Model: %s", metric_name, model_type))
  message("==========================================================")
  
  if (!file.exists(input_path)) stop("Missing monthly input: ", input_path)
  heat_data <- read.csv(input_path)
  if (!metric_name %in% names(heat_data)) stop("Missing outcome column: ", metric_name)
  full_data_log <- prepare_nexus_data_log(heat_data, metric_name)
  
  all_states <- sort(unique(full_data_log$State))
  ncei_regions <- sort(unique(na.omit(full_data_log$NCEI_Region)))
  target_regions <- c("CONUS", ncei_regions, all_states)
  
  all_threshold_results <- list()
  
  for (current_region in target_regions) {
    if (current_region == "CONUS") {
      region_data <- full_data_log
      reg_type <- "CONUS"
    } else if (current_region %in% ncei_regions) {
      region_data <- full_data_log %>% filter(NCEI_Region == current_region)
      reg_type <- "NCEI Region"
    } else {
      region_data <- full_data_log %>% filter(State == current_region)
      reg_type <- "Statewise"
    }
    
    if (nrow(region_data) < 20) next
    message("--> Processing Region: ", current_region, " (", reg_type, ")")
    
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
    
    if (nrow(combined_region_thr) > 0) {
      all_threshold_results[[current_region]] <- combined_region_thr
    }
    
    # --- Generate 8 Climate Variables in 1 Plot Grid with Super Big Labels ---
    if (nrow(combined_region_ale) > 0) {
      plot_list <- list()
      
      for (var_name in climate_vars) {
        var_ale <- combined_region_ale %>% filter(Variable == var_name)
        var_thr <- combined_region_thr %>% filter(Variable == var_name)
        
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
        
        # Placed unit on a new line (\n) so labels fit without truncation
        x_axis_title <- paste0(meta$title, "\n(", meta$unit, ")")
        
        p <- ggplot(var_ale, aes(x = x_val, group = Period, color = Period)) +
          geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1.0) +
          geom_line(aes(y = lower_ci), linetype = "dashed", linewidth = 1.1, alpha = 0.85) +
          geom_line(aes(y = upper_ci), linetype = "dashed", linewidth = 1.1, alpha = 0.85) +
          geom_line(aes(y = `.value`), linewidth = 2.2) +
          scale_color_manual(values = period_colors, labels = legend_mapping) +
          guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
          theme_classic(base_size = 20) +
          labs(
            title = meta$title,
            x = x_axis_title,
            y = "ALE (log10(outcome + 1))"
          ) +
          theme(
            plot.title = element_text(face = "bold", size = 24, color = "black", hjust = 0),
            axis.title.y = element_text(face = "bold", size = 22, color = "black"),
            axis.title.x = element_text(face = "bold", size = 20, color = "black", lineheight = 0.85),
            axis.text   = element_text(face = "bold", size = 18, color = "black"),
            legend.position = "top",
            legend.justification = "left",
            legend.title = element_blank(),
            legend.text  = element_text(size = 18, face = "bold"),
            legend.key.width  = unit(1.8, "cm"),
            legend.key.height = unit(0.8, "cm"),
            legend.margin = margin(b = -2),
            plot.margin = margin(t = 10, r = 14, b = 10, l = 14),
            axis.line  = element_line(linewidth = 1.0, color = "black"),
            axis.ticks = element_line(linewidth = 1.0, color = "black")
          )
        
        if (length(hist_thr) == 1 && !is.na(hist_thr)) {
          p <- p + geom_vline(xintercept = hist_thr, linetype = "dashed", color = period_colors["Historical (2015–2018)"], linewidth = 1.1, alpha = 0.8)
        }
        if (length(rec_thr) == 1 && !is.na(rec_thr)) {
          p <- p + geom_vline(xintercept = rec_thr, linetype = "dashed", color = period_colors["Recent (2019–2022)"], linewidth = 1.1, alpha = 0.8)
        }
        
        plot_list[[var_name]] <- p
      }
      
      if (length(plot_list) > 0) {
        combined_header <- ggdraw() + 
          draw_label(
            paste0("b  Climate driver impact on outage risk — ", current_region, " (", reg_type, ")"), 
            fontface = "bold", size = 32, x = 0.02, hjust = 0
          )
        
        grid_8var <- plot_grid(plotlist = plot_list, ncol = 4, nrow = 2, align = "vh")
        final_figure <- plot_grid(combined_header, grid_8var, ncol = 1, rel_heights = c(0.08, 1))
        
        ggsave(
          file.path(save_path, paste0(current_region, "_All8Variables_TwoPeriod_ALE_Grid.png")), 
          final_figure, width = 25, height = 13.5, dpi = 300
        )
      }
    }
  }
  
  if (length(all_threshold_results) > 0) {
    final_threshold_table <- bind_rows(all_threshold_results) %>%
      pivot_wider(names_from = Period, values_from = Safety_Window_Threshold_x) %>%
      rename(
        Historical_Threshold_x = `Historical (2015–2018)`,
        Recent_Threshold_x     = `Recent (2019–2022)`
      ) %>%
      mutate(Safety_Window_Shift = Recent_Threshold_x - Historical_Threshold_x) %>%
      arrange(factor(Region_Type, levels = c("CONUS", "NCEI Region", "Statewise")), Region, Variable)
    
    csv_out_path <- file.path(save_path, "All_Regions_TwoPeriod_SafetyWindow_Thresholds.csv")
    write.csv(final_threshold_table, csv_out_path, row.names = FALSE)
    message(">> Exported complete CSV with CONUS, 9 NCEI regions, and Statewise data to: ", csv_out_path)
  }
}

# --- Figure 3b Executions: 2 Target Metrics x 2 Models ---
run_fig3b(metric_name = "monthly_coh",          model_type = "XGBoost")
run_fig3b(metric_name = "monthly_coh",          model_type = "RF")
run_fig3b(metric_name = "monthly_outage_hours", model_type = "XGBoost")
run_fig3b(metric_name = "monthly_outage_hours", model_type = "RF")