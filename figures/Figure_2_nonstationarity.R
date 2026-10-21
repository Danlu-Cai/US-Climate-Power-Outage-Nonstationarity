# =============================================================================
# Script: Figure_2_nonstationarity.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
#   Generate Figure 2a–2f panels for temporal nonstationarity, out-of-time
#   prediction, within-year signal-to-noise, predictor importance, and
#   temperature-envelope analyses. Uses county-year merged CONUS data.
#
# Input: data/processed/Digicare_Merged_Yearly_2015_2022_CONUS.csv
#        Override directory with FIGURE2_INPUT_DIR.
# Output: results/figures/figure2/ (override with FIGURE2_OUTPUT_DIR)
# Dependencies: dplyr, readr, ggplot2, tidyr, purrr, xgboost, ranger,
#               stringr, conflicted, ggrepel, tidyverse, patchwork, scales
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)
library(xgboost)
library(ranger)
library(stringr)
library(conflicted)
library(ggrepel)

conflicts_prefer(dplyr::lag)

# ==============================================================================
# --- 1. Global Setup & Core Helper Functions ---
# ==============================================================================

Input_ROOT <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related" 
save_path  <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related_final" 
csv_path   <- file.path(save_path, "Fig2_csv_results") 

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
if (!dir.exists(csv_path))  dir.create(csv_path, recursive = TRUE)

load_and_prep_data <- function(metric_name) {
  file_path <- file.path(Input_ROOT, "Digicare_Merged_Yearly_2015_2022_CONUS.csv")
  if (!file.exists(file_path)) stop("Data file not found at: ", file_path)
  
  df_full <- read_csv(file_path, show_col_types = FALSE)
  
  pop_col     <- base::intersect(c("population", "Pop", "pop"), names(df_full))[1]
  density_col <- base::intersect(c("pop_density", "Pop_Density", "density"), names(df_full))[1]
  spatial_exposure_vars <- if (!is.na(density_col)) density_col else if (!is.na(pop_col)) pop_col else character(0)
  
  unweighted_degree_vars <- c("annual_CDD", "annual_HDD", spatial_exposure_vars)
  base_hazard_vars       <- c("avg_T2C", "annual_event_count", "annual_max_mag", "annual_heat_days", "annual_cold_days", "annual_storm_days", "state")
  predictors             <- setdiff(c(base_hazard_vars, unweighted_degree_vars), c("county_area", "County_Area", "area", "Area"))
  vars_to_use            <- base::intersect(predictors, names(df_full))
  
  df_transformed <- df_full %>%
    mutate(!!metric_name := log10(as.numeric(.data[[metric_name]]) + 1)) %>%
    filter(!is.infinite(.data[[metric_name]]), !is.na(.data[[metric_name]]))
  
  county_col <- base::intersect(c("fips", "fips_code", "county"), names(df_transformed))[1]
  if (is.na(county_col)) stop("Error: County column missing.")
  
  if ("year" %in% names(df_transformed)) {
    df_transformed <- df_transformed %>% mutate(Year_Internal = as.numeric(year))
  } else if ("Year" %in% names(df_transformed)) {
    df_transformed <- df_transformed %>% mutate(Year_Internal = as.numeric(Year))
  } else if ("date" %in% names(df_transformed)) {
    df_transformed <- df_transformed %>% mutate(Year_Internal = as.numeric(format(as.Date(date), "%Y")))
  }
  
  clean_df <- df_transformed %>%
    select(county_group = !!sym(county_col), Year_Internal, target = !!sym(metric_name), all_of(vars_to_use)) %>%
    mutate(across(where(is.character) & !county_group, ~as.numeric(as.factor(.)))) %>%
    mutate(across(-county_group, ~as.numeric(as.character(.)))) %>%
    drop_na()
  
  set.seed(42)
  unique_counties <- unique(clean_df$county_group)
  county_fold_map <- data.frame(
    county_group = unique_counties,
    fold = sample(rep(1:5, length.out = length(unique_counties)))
  )
  
  clean_df %>% left_join(county_fold_map, by = "county_group")
}

fit_predict <- function(tr_df, te_df, feat_cols, model_type) {
  if (model_type == "XGBoost") {
    mod <- xgboost(data = as.matrix(tr_df[, feat_cols]), label = tr_df$target,
                   nrounds = 100, objective = "reg:squarederror", verbosity = 0, 
                   params = list(eta = 0.1, max_depth = 6))
    return(predict(mod, as.matrix(te_df[, feat_cols])))
  } else {
    mod <- ranger(target ~ ., data = tr_df[, c("target", feat_cols)], num.trees = 200)
    return(predict(mod, data = te_df[, c("target", feat_cols)])$predictions)
  }
}

evaluate_decay_data <- function(clean_df, model_type, n_folds = 5) {
  feature_cols <- setdiff(names(clean_df), c("county_group", "Year_Internal", "target", "fold"))
  target_years <- 2019:2022
  hist_windows <- list(
    "2015–2016" = 2015:2016, "2015–2017" = 2015:2017,
    "2015–2018" = 2015:2018, "2015–2019" = 2015:2019
  )
  
  decay_results <- list()
  for (yr in target_years) {
    f_act <- numeric(0); f_contemp <- numeric(0)
    for (k in 1:n_folds) {
      test_sub   <- clean_df %>% filter(Year_Internal == yr & fold == k)
      contemp_tr <- clean_df %>% filter(Year_Internal >= 2015 & Year_Internal < yr & fold != k)
      if (nrow(test_sub) == 0) next
      p_c <- fit_predict(contemp_tr, test_sub, feature_cols, model_type)
      f_act <- c(f_act, test_sub$target); f_contemp <- c(f_contemp, p_c)
    }
    tss <- sum((f_act - mean(f_act))^2)
    r2_contemp <- max(0, 1 - (sum((f_act - f_contemp)^2) / tss))
    
    for (win_name in names(hist_windows)) {
      win_years <- hist_windows[[win_name]]
      if (max(win_years) >= yr) next 
      
      f_hist <- numeric(0)
      for (k in 1:n_folds) {
        test_sub <- clean_df %>% filter(Year_Internal == yr & fold == k)
        hist_tr  <- clean_df %>% filter(Year_Internal %in% win_years & fold != k)
        if (nrow(test_sub) == 0) next
        p_h <- fit_predict(hist_tr, test_sub, feature_cols, model_type)
        f_hist <- c(f_hist, p_h)
      }
      
      r2_hist <- max(0, 1 - (sum((f_act - f_hist)^2) / tss))
      decay_results[[paste(win_name, yr, sep = "_")]] <- data.frame(
        Model_Window    = win_name, Target_Year = as.character(yr),
        R2_Out_Of_Time  = round(r2_hist, 4), R2_Contemporary = round(r2_contemp, 4),
        Delta_R2        = round(r2_hist - r2_contemp, 4)
      )
    }
  }
  bind_rows(decay_results)
}

evaluate_snr_data <- function(clean_df, model_type = "XGBoost", n_folds = 5) {
  feature_cols <- setdiff(names(clean_df), c("county_group", "Year_Internal", "target", "fold"))
  years <- 2015:2022
  snr_results <- list()
  
  for (yr in years) {
    # 1. Filter data exclusively to the current target year
    year_data <- clean_df %>% filter(Year_Internal == yr)
    if (nrow(year_data) < 50) next
    
    unique_counties <- unique(year_data$county_group)
    if (length(unique_counties) < n_folds) next
    
    # 2. Perform 5-fold spatial cross-validation WITHIN this single year
    set.seed(42 + yr)
    county_fold_map <- data.frame(
      county_group = unique_counties,
      fold_in_year = sample(rep(1:n_folds, length.out = length(unique_counties)))
    )
    
    year_data <- year_data %>% left_join(county_fold_map, by = "county_group")
    oof_preds <- numeric(nrow(year_data))
    
    for (k in 1:n_folds) {
      tr_sub <- year_data %>% filter(fold_in_year != k)
      te_sub <- year_data %>% filter(fold_in_year == k)
      te_idx <- which(year_data$fold_in_year == k)
      
      if (nrow(te_sub) == 0) next
      preds <- fit_predict(tr_sub, te_sub, feature_cols, model_type)
      oof_preds[te_idx] <- preds
    }
    
    # 3. Calculate Signal SD and Residual Noise SD
    signal_sd <- sd(oof_preds, na.rm = TRUE)
    noise_sd  <- sd(year_data$target - oof_preds, na.rm = TRUE)
    
    # 4. Standard Deviation SNR Ratio (SD_Signal / SD_Noise)
    snr_val <- if (noise_sd > 0) signal_sd / noise_sd else 0
    snr_results[[as.character(yr)]] <- data.frame(
      Year = yr,
      SNR  = round(snr_val, 2)
    )
  }
  bind_rows(snr_results)
}

# --- GLOBAL LARGE-FONT THEME FOR ALL FIGURES ---
theme_large_nature <- function() {
  theme_classic() +
    theme(
      panel.grid            = element_blank(),
      panel.background      = element_blank(),
      plot.background       = element_blank(),
      legend.background     = element_blank(),
      legend.box.background = element_blank(),
      legend.key            = element_blank(),
      plot.title            = element_text(face = "bold", size = 18, margin = margin(b = 12)),
      axis.title            = element_text(face = "bold", size = 16),
      axis.text             = element_text(face = "bold", size = 14, color = "black"),
      legend.title          = element_text(face = "bold", size = 15),
      legend.text           = element_text(face = "bold", size = 13),
      axis.line             = element_line(color = "black", linewidth = 0.8)
    )
}


# ==============================================================================
# --- 2. INDIVIDUAL MODULE: Figure 2a ---
# ==============================================================================

run_fig2a <- function(metric_name = "annual_coh", model_type = "XGBoost", n_folds = 5) {
  message(paste("[FIGURE 2a] Running for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data(metric_name)
  feature_cols <- setdiff(names(clean_data), c("county_group", "Year_Internal", "target", "fold"))
  
  res_list <- list()
  for (yr in 2017:2022) {
    period_label <- paste0("2015–", yr)
    f_act <- numeric(0); f_hist <- numeric(0); f_contemp <- numeric(0)
    
    for (k in 1:n_folds) {
      test_sub    <- clean_data %>% filter(Year_Internal == yr & fold == k)
      hist_tr     <- clean_data %>% filter(Year_Internal %in% c(2015, 2016) & fold != k)
      contemp_tr  <- clean_data %>% filter(Year_Internal >= 2015 & Year_Internal < yr & fold != k)
      if (nrow(test_sub) == 0) next
      
      p_h <- fit_predict(hist_tr, test_sub, feature_cols, model_type)
      p_c <- fit_predict(contemp_tr, test_sub, feature_cols, model_type)
      f_act <- c(f_act, test_sub$target); f_hist <- c(f_hist, p_h); f_contemp <- c(f_contemp, p_c)
    }
    
    tss <- sum((f_act - mean(f_act))^2)
    res_list[[period_label]] <- data.frame(
      Training_Period = period_label, Test_Year = as.character(yr),
      R2_Historical   = round(max(0, 1 - (sum((f_act - f_hist)^2) / tss)), 4),
      R2_Contemporary = round(max(0, 1 - (sum((f_act - f_contemp)^2) / tss)), 4)
    )
  }
  
  res_df <- bind_rows(res_list) %>% mutate(Training_Period = factor(Training_Period, levels = unique(Training_Period)))
  write_csv(res_df, file.path(csv_path, paste0("Fig2a_Data_", metric_name, "_", model_type, ".csv")))
  
  metric_label <- case_when(
    metric_name == "annual_coh" ~ "COH",
    metric_name == "annual_outage_hours" ~ "TOH",
    TRUE ~ metric_name
  )
  
  all_r2   <- c(res_df$R2_Historical, res_df$R2_Contemporary)
  data_min <- min(all_r2)
  data_max <- max(all_r2)
  
  y_min <- floor((data_min - 0.08) * 10) / 10
  y_max <- ceiling((data_max + 0.25) * 10) / 10 
  
  lbl_hist <- "Historical model (trained on 2015–2016 only)"
  lbl_cont <- "Contemporary model (rolling, up to current year)"
  
  p <- ggplot(res_df, aes(x = Training_Period)) +
    annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, fill = "#eef4fa", alpha = 0.8) +
    annotate("rect", xmin = 3, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#fdf0ed", alpha = 0.8) +
    
    geom_line(aes(y = R2_Contemporary, group = 1, color = lbl_cont), linewidth = 1.3) +
    geom_point(aes(y = R2_Contemporary, color = lbl_cont), size = 4.5) +
    geom_line(aes(y = R2_Historical, group = 1, color = lbl_hist), linewidth = 1.3) +
    geom_point(aes(y = R2_Historical, color = lbl_hist), size = 4.5) +
    
    geom_text(aes(y = R2_Contemporary, label = sprintf("%.2f", R2_Contemporary)), vjust = -1.2, color = "#d9381e", fontface = "bold", size = 5.5) +
    geom_text(aes(y = R2_Historical, label = sprintf("%.2f", R2_Historical)), vjust = 1.8, color = "#1f77b4", fontface = "bold", size = 5.5) +
    
    scale_color_manual(
      values = setNames(c("#1f77b4", "#d9381e"), c(lbl_hist, lbl_cont)),
      breaks = c(lbl_hist, lbl_cont)
    ) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = seq(max(0, y_min), y_max, by = 0.1)) +
    coord_cartesian(clip = "off") +
    labs(
      title    = NULL,
      subtitle = NULL,
      x        = "Training period (years)", 
      y        = expression(bold("Out-of-sample ") ~ bold(R^2)), 
      color    = NULL
    ) +
    theme_large_nature() +
    theme(
      legend.position      = c(0.02, 0.96),
      legend.justification = c(0, 1),
      legend.direction     = "vertical",
      plot.margin          = margin(t = 15, r = 15, b = 15, l = 15)
    )
  
  ggsave(file.path(save_path, paste0("Fig2a_", metric_name, "_", model_type, ".png")), p, width = 8, height = 4, dpi = 300, bg = "white")
  return(p)
}

# --- Figure 2a Executions ---
run_fig2a(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2a(metric_name = "annual_coh", model_type = "RF")
run_fig2a(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2a(metric_name = "annual_outage_hours", model_type = "RF")


# ==============================================================================
# --- 3. INDIVIDUAL MODULE: Figure 2b ---
# ==============================================================================
run_fig2b <- function(metric_name = "annual_outage_hours", model_type = "RF", n_folds = 5, plot_loss = FALSE) {
  message(paste("[FIGURE 2b] Running for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data(metric_name)
  feature_cols <- setdiff(names(clean_data), c("county_group", "Year_Internal", "target", "fold"))
  
  train_periods <- list(
    "2015–2016" = c(2015, 2016),
    "2015–2017" = 2015:2017,
    "2015–2018" = 2015:2018,
    "2015–2019" = 2015:2019
  )
  
  contemp_r2_map <- list()
  if (plot_loss) {
    for (t_yr in 2020:2022) {
      f_act <- numeric(0); f_pred <- numeric(0)
      for (k in 1:n_folds) {
        test_sub  <- clean_data %>% filter(Year_Internal == t_yr & fold == k)
        train_sub <- clean_data %>% filter(Year_Internal >= 2015 & Year_Internal < t_yr & fold != k)
        if (nrow(test_sub) == 0) next
        preds <- fit_predict(train_sub, test_sub, feature_cols, model_type)
        f_act  <- c(f_act, test_sub$target)
        f_pred <- c(f_pred, preds)
      }
      tss <- sum((f_act - mean(f_act))^2)
      contemp_r2_map[[as.character(t_yr)]] <- max(0, 1 - (sum((f_act - f_pred)^2) / tss))
    }
  }
  
  res_list <- list()
  for (p_name in names(train_periods)) {
    tr_years <- train_periods[[p_name]]
    eval_years <- 2020:2022 
    
    for (t_yr in eval_years) {
      f_act <- numeric(0); f_pred <- numeric(0)
      for (k in 1:n_folds) {
        test_sub  <- clean_data %>% filter(Year_Internal == t_yr & fold == k)
        train_sub <- clean_data %>% filter(Year_Internal %in% tr_years & fold != k)
        if (nrow(test_sub) == 0) next
        
        preds  <- fit_predict(train_sub, test_sub, feature_cols, model_type)
        f_act  <- c(f_act, test_sub$target)
        f_pred <- c(f_pred, preds)
      }
      
      tss <- sum((f_act - mean(f_act))^2)
      r2_val <- max(0, 1 - (sum((f_act - f_pred)^2) / tss))
      
      if (plot_loss) {
        val_to_plot <- contemp_r2_map[[as.character(t_yr)]] - r2_val
      } else {
        val_to_plot <- r2_val
      }
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Trained_Period = p_name,
        Target_Year    = as.character(t_yr),
        Value          = round(val_to_plot, 4)
      )
    }
  }
  
  res_df <- bind_rows(res_list)
  write_csv(res_df, file.path(csv_path, paste0("Fig2b_Data_", metric_name, "_", model_type, ".csv")))
  
  metric_label <- case_when(
    metric_name == "annual_coh" ~ "COH",
    metric_name == "annual_outage_hours" ~ "TOH",
    TRUE ~ metric_name
  )
  
  y_label <- if (plot_loss) expression(bold(R^2) ~ bold("Loss (vs. contemporary)")) else expression(bold("Out-of-time ") ~ bold(R^2))
  
  # Y-axis bounds set to [0.1, 0.5]
  y_min <- 0.1
  y_max <- 0.5
  
  p <- ggplot(res_df, aes(x = Target_Year, y = Value, color = Trained_Period, group = Trained_Period)) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 4) +
    
    geom_text_repel(
      aes(label = sprintf("%.2f", Value)),
      fontface          = "bold",
      size              = 5.5,
      box.padding       = 0.6,
      point.padding     = 0.5,
      force             = 8,
      force_pull        = 0.5,
      max.overlaps      = Inf,
      min.segment.length = 0.1,
      segment.color     = "grey50",
      segment.linewidth = 0.5,
      show.legend       = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "2015–2016" = "#1f77b4",
        "2015–2017" = "#ff7f0e",
        "2015–2018" = "#2ca02c",
        "2015–2019" = "#9467bd"
      )
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, by = 0.1)) +
    coord_cartesian(clip = "off") +
    labs(
      title = NULL,
      x     = "Target year (out-of-time)",
      y     = y_label,
      color = "Model trained on:"
    ) +
    guides(color = guide_legend(nrow = 1, title.position = "left")) +
    theme_large_nature() +
    theme(
      legend.position       = c(0.02, 0.96),
      legend.justification  = c(0, 1),
      legend.direction      = "horizontal",
      legend.background     = element_blank(),
      legend.box.background = element_blank(),
      plot.margin           = margin(t = 15, r = 15, b = 15, l = 15)
    )
  
  ggsave(file.path(save_path, paste0("Fig2b_", metric_name, "_", model_type, ".png")), p, width = 8, height = 4, dpi = 300, bg = "white")
  return(p)
}

# --- Figure 2b Executions ---
run_fig2b(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2b(metric_name = "annual_coh", model_type = "RF")
run_fig2b(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2b(metric_name = "annual_outage_hours", model_type = "RF")


# ==============================================================================
# --- 4. INDIVIDUAL MODULE: Figure 2c ---
# ==============================================================================
evaluate_decay_data <- function(clean_df, model_type, n_folds = 5) {
  feature_cols <- setdiff(names(clean_df), c("county_group", "Year_Internal", "target", "fold"))
  target_years <- 2019:2022
  hist_windows <- list(
    "2015–2016" = 2015:2016, "2015–2017" = 2015:2017,
    "2015–2018" = 2015:2018, "2015–2019" = 2015:2019
  )
  
  decay_results <- list()
  for (yr in target_years) {
    f_act <- numeric(0); f_contemp <- numeric(0)
    for (k in 1:n_folds) {
      test_sub   <- clean_df %>% filter(Year_Internal == yr & fold == k)
      contemp_tr <- clean_df %>% filter(Year_Internal >= 2015 & Year_Internal < yr & fold != k)
      if (nrow(test_sub) == 0) next
      p_c <- fit_predict(contemp_tr, test_sub, feature_cols, model_type)
      f_act <- c(f_act, test_sub$target); f_contemp <- c(f_contemp, p_c)
    }
    tss <- sum((f_act - mean(f_act))^2)
    r2_contemp <- max(0, 1 - (sum((f_act - f_contemp)^2) / tss))
    
    for (win_name in names(hist_windows)) {
      win_years <- hist_windows[[win_name]]
      if (max(win_years) >= yr) next 
      
      f_hist <- numeric(0)
      for (k in 1:n_folds) {
        test_sub <- clean_df %>% filter(Year_Internal == yr & fold == k)
        hist_tr  <- clean_df %>% filter(Year_Internal %in% win_years & fold != k)
        if (nrow(test_sub) == 0) next
        p_h <- fit_predict(hist_tr, test_sub, feature_cols, model_type)
        f_hist <- c(f_hist, p_h)
      }
      
      r2_hist <- max(0, 1 - (sum((f_act - f_hist)^2) / tss))
      decay_results[[paste(win_name, yr, sep = "_")]] <- data.frame(
        Model_Window    = win_name, Target_Year = as.character(yr),
        R2_Out_Of_Time  = round(r2_hist, 4), R2_Contemporary = round(r2_contemp, 4),
        Delta_R2        = round(r2_hist - r2_contemp, 4)
      )
    }
  }
  bind_rows(decay_results)
}

run_fig2c <- function(metric_name = "annual_coh", model_type = "XGBoost") {
  message(paste("[FIGURE 2c] Running for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data(metric_name)
  
  # Fetch decay data and format factor levels
  perf_df <- evaluate_decay_data(clean_data, model_type) %>% 
    filter(Target_Year >= 2020) %>% 
    mutate(
      Model_Window = factor(
        Model_Window, 
        levels = c("2015–2016", "2015–2017", "2015–2018", "2015–2019")
      )
    )
  
  write_csv(perf_df, file.path(csv_path, paste0("Fig2c_Data_", metric_name, "_", model_type, ".csv")))
  
  # Dynamic Y-limits with clean step intervals
  min_val  <- min(perf_df$Delta_R2, na.rm = TRUE)
  max_val  <- max(perf_df$Delta_R2, na.rm = TRUE)
  y_min    <- min(-0.25, floor((min_val - 0.08) * 10) / 10)
  y_max    <- max(0.10, ceiling((max_val + 0.15) * 10) / 10)
  
  # zapsmall removes floating precision artifacts (e.g. 5.55e-17 becomes 0)
  y_breaks <- zapsmall(seq(y_min, y_max, by = 0.1))
  
  p <- ggplot(perf_df, aes(x = factor(Target_Year), y = Delta_R2, fill = Model_Window)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    
    geom_text(
      aes(
        label = sprintf("%.2f", Delta_R2),
        y     = ifelse(Delta_R2 >= 0, Delta_R2 + 0.02, Delta_R2 - 0.02),
        vjust = ifelse(Delta_R2 >= 0, 0, 1)
      ),
      position = position_dodge(width = 0.8),
      fontface = "bold",
      size     = 4.8,
      color    = "black"
    ) +
    
    scale_fill_manual(
      values = c(
        "2015–2016" = "#1f77b4",
        "2015–2017" = "#ff7f0e",
        "2015–2018" = "#2ca02c",
        "2015–2019" = "#9467bd"
      ),
      name = "Model trained on:"
    ) +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    
    # Formats ticks cleanly as -0.3, -0.2, -0.1, 0.0, 0.1
    scale_y_continuous(
      breaks = y_breaks,
      labels = sprintf("%.1f", y_breaks)
    ) +
    labs(
      title = NULL,
      x     = "Target year (out-of-time)",
      y     = expression(bold(Delta * R^2 ~ "(out-of-time)"))
    ) +
    theme_large_nature() +
    theme(
      legend.position   = "top",
      legend.box.margin = margin(b = -5),
      plot.margin       = margin(t = 15, r = 15, b = 15, l = 15)
    )
  
  ggsave(
    file.path(save_path, paste0("Fig2c_", metric_name, "_", model_type, ".png")),
    p,
    width  = 8,
    height = 4,
    dpi    = 300,
    bg     = "white"
  )
  
  return(p)
}

# --- Figure 2c Executions ---
run_fig2c(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2c(metric_name = "annual_coh", model_type = "RF")
run_fig2c(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2c(metric_name = "annual_outage_hours", model_type = "RF")


# ==============================================================================
# --- 5. INDIVIDUAL MODULE: Figure 2d ---
# ==============================================================================

run_fig2d <- function(metric_name = "annual_coh", model_type = "XGBoost") {
  message(paste("[FIGURE 2d] Running for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data(metric_name)
  snr_df     <- evaluate_snr_data(clean_data, model_type)
  
  write_csv(snr_df, file.path(csv_path, paste0("Fig2d_", metric_name, "_", model_type, ".csv")))
  
  p <- ggplot(snr_df, aes(x = Year, y = SNR)) +
    # Background shaded regions
    annotate("rect", xmin = 2014.5, xmax = 2018.5, ymin = 0.8, ymax = 1.1, fill = "#eef4fa", alpha = 0.8) +
    annotate("rect", xmin = 2018.5, xmax = 2022.5, ymin = 0.8, ymax = 1.1, fill = "#fdf0ed", alpha = 0.8) +
    
    # Period header annotations
    annotate("text", x = 2016.5, y = 1.07, label = "Historical period\n(2015–2018)", fontface = "bold", color = "#1a3b68", size = 4.8) +
    annotate("text", x = 2020.5, y = 1.07, label = "Recent period\n(2019–2022)", fontface = "bold", color = "#333333", size = 4.8) +
    
    # Line and points
    geom_line(color = "#1f77b4", linewidth = 1.3) +
    geom_point(color = "#1f77b4", shape = 15, size = 4.5) +
    
    # Data labels
    geom_text(
      aes(label = sprintf("%.2f", SNR)),
      vjust    = -1.2,
      fontface = "bold",
      size     = 5.5,
      color    = "black"
    ) +
    
    scale_x_continuous(breaks = 2015:2022, limits = c(2014.5, 2022.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.8, 1.1), breaks = seq(0.8, 1.1, by = 0.1)) +
    coord_cartesian(clip = "off") +
    labs(
      title = NULL,
      x     = "Year",
      y     = "Signal-to-noise ratio"
    ) +
    theme_large_nature() +
    theme(
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )
  
  ggsave(
    file.path(save_path, paste0("Fig2d_", metric_name, "_", model_type, ".png")),
    p,
    width  = 8,
    height = 4,
    dpi    = 300,
    bg     = "white"
  )
  
  return(p)
}

# --- Figure 2d Executions ---
run_fig2d(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2d(metric_name = "annual_coh", model_type = "RF")
run_fig2d(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2d(metric_name = "annual_outage_hours", model_type = "RF")


# ==============================================================================
# --- 6. INDIVIDUAL MODULE: Figure 2e ---
# ==============================================================================

evaluate_importance_data <- function(clean_df, model_type = "XGBoost") {
  
  # ---------------------------------------------------------------------------
  # 1. Identify predictor columns
  # ---------------------------------------------------------------------------
  feature_cols <- setdiff(
    names(clean_df),
    c("county_group", "Year_Internal", "target", "fold")
  )
  
  # ---------------------------------------------------------------------------
  # 2. Publication labels
  # ---------------------------------------------------------------------------
  label_map <- c(
    "avg_T2C"            = "Temperature",
    "annual_CDD"          = "Cooling demand (CDD)",
    "annual_HDD"          = "Heating demand (HDD)",
    "annual_heat_days"   = "Extreme heat days",
    "annual_cold_days"   = "Extreme cold days",
    "annual_storm_days"  = "Precipitation extremes",
    "annual_event_count" = "Event frequency",
    "annual_max_mag"     = "Maximum event magnitude",
    "pop_density"        = "Population density",
    "population"         = "Population",
    "Pop"                = "Population",
    "state"              = "State"
  )
  
  # ---------------------------------------------------------------------------
  # 3. Define historical and recent periods
  # ---------------------------------------------------------------------------
  periods <- list(
    "Historical" = 2015:2018,
    "Recent"     = 2019:2022
  )
  
  # ---------------------------------------------------------------------------
  # 4. Fit separate models and calculate relative importance
  # ---------------------------------------------------------------------------
  imp_results <- list()
  
  for (p_name in names(periods)) {
    p_years <- periods[[p_name]]
    sub_df  <- clean_df %>% filter(Year_Internal %in% p_years)
    
    if (model_type == "XGBoost") {
      mod <- xgboost(
        data      = as.matrix(sub_df[, feature_cols]),
        label     = sub_df$target,
        nrounds   = 100,
        objective = "reg:squarederror",
        verbosity = 0,
        params    = list(eta = 0.1, max_depth = 6)
      )
      imp <- xgb.importance(model = mod)
      imp_df <- data.frame(
        Feature    = imp$Feature,
        Importance = imp$Gain
      )
    } else if (model_type == "RF") {
      mod <- ranger(
        target ~ .,
        data       = sub_df[, c("target", feature_cols)],
        num.trees  = 200,
        importance = "permutation"
      )
      imp_vals <- importance(mod)
      imp_df   <- data.frame(
        Feature    = names(imp_vals),
        Importance = as.numeric(imp_vals)
      )
    } else {
      stop("model_type must be either 'XGBoost' or 'RF'.")
    }
    
    # Normalize importance within each period
    imp_df <- imp_df %>% mutate(Importance = pmax(0, Importance))
    total_importance <- sum(imp_df$Importance, na.rm = TRUE)
    
    if (total_importance <= 0) {
      stop(paste("Total feature importance is zero for", p_name, model_type))
    }
    
    imp_df <- imp_df %>%
      mutate(
        Rel_Importance = Importance / total_importance * 100,
        Feature_Label  = ifelse(Feature %in% names(label_map), unname(label_map[Feature]), Feature)
      )
    
    imp_results[[p_name]] <- imp_df
  }
  
  # ---------------------------------------------------------------------------
  # 5. Combine historical and recent importance
  # ---------------------------------------------------------------------------
  hist_df <- imp_results[["Historical"]] %>% select(Feature_Label, Rel_Hist = Rel_Importance)
  recent_df <- imp_results[["Recent"]] %>% select(Feature_Label, Rel_Recent = Rel_Importance)
  
  diff_df <- full_join(hist_df, recent_df, by = "Feature_Label") %>%
    mutate(
      Rel_Hist   = ifelse(is.na(Rel_Hist), 0, Rel_Hist),
      Rel_Recent = ifelse(is.na(Rel_Recent), 0, Rel_Recent),
      Delta_pp   = Rel_Recent - Rel_Hist
    )
  
  # ---------------------------------------------------------------------------
  # 6. Target features & aggregation
  # ---------------------------------------------------------------------------
  target_labels <- c(
    "Cooling demand (CDD)",
    "State",
    "Population density",
    "Temperature",
    "Heating demand (HDD)",
    "Event frequency",
    "Precipitation extremes"
  )
  
  main_features  <- diff_df %>% filter(Feature_Label %in% target_labels)
  other_features <- diff_df %>% filter(!Feature_Label %in% target_labels)
  
  if (nrow(other_features) > 0) {
    other_row <- data.frame(
      Feature_Label = "Other",
      Rel_Hist      = sum(other_features$Rel_Hist, na.rm = TRUE),
      Rel_Recent    = sum(other_features$Rel_Recent, na.rm = TRUE)
    ) %>%
      mutate(Delta_pp = Rel_Recent - Rel_Hist)
    
    final_df <- bind_rows(main_features, other_row)
  } else {
    final_df <- main_features
  }
  
  final_df <- final_df %>%
    mutate(
      Change_Type = case_when(
        Delta_pp > 0 ~ "Increase",
        Delta_pp < 0 ~ "Decrease",
        TRUE         ~ "No change"
      )
    )
  
  feature_order <- c(
    "Other",
    "Precipitation extremes",
    "Event frequency",
    "Heating demand (HDD)",
    "Temperature",
    "Population density",
    "State",
    "Cooling demand (CDD)"
  )
  
  existing_order <- intersect(feature_order, final_df$Feature_Label)
  final_df$Feature_Label <- factor(final_df$Feature_Label, levels = existing_order)
  
  final_df <- final_df %>%
    mutate(
      Model = model_type,
      Importance_Metric = ifelse(model_type == "XGBoost", "Relative gain", "Relative permutation importance")
    )
  
  return(final_df)
}

# ==============================================================================
# --- Figure 2e Plotting Function ---
# ==============================================================================

run_fig2e <- function(
    metric_name = "annual_coh",
    model_type  = "XGBoost"
) {
  
  message(paste("[FIGURE 2e] Running for Metric:", metric_name, "| Model:", model_type))
  
  clean_data <- load_and_prep_data(metric_name)
  imp_df     <- evaluate_importance_data(clean_data, model_type)
  
  write_csv(
    imp_df,
    file.path(csv_path, paste0("Fig2e_Data_", metric_name, "_", model_type, ".csv"))
  )
  
  max_abs_delta <- max(abs(imp_df$Delta_pp), na.rm = TRUE)
  x_limit       <- max(6, ceiling(max_abs_delta + 1))
  n_features    <- length(unique(imp_df$Feature_Label))
  
  # ---------------------------------------------------------------------------
  # 1. Model-specific x-axis label with line break (\n)
  # ---------------------------------------------------------------------------
  if (model_type == "XGBoost") {
    x_axis_label <- "Change in relative gain\n(Recent \u2212 Historical, percentage points)"
  } else {
    x_axis_label <- "Change in relative permutation importance\n(Recent \u2212 Historical, percentage points)"
  }
  
  # ---------------------------------------------------------------------------
  # 2. Plotting
  # ---------------------------------------------------------------------------
  p <- ggplot(
    imp_df,
    aes(x = Delta_pp, y = Feature_Label, fill = Change_Type)
  ) +
    geom_col(width = 0.65) +
    
    # Zero reference
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      color      = "black",
      linewidth  = 0.6
    ) +
    
    # Direction colors
    scale_fill_manual(
      values = c(
        "Increase"  = "#d9381e",
        "Decrease"  = "#1f77b4",
        "No change" = "grey70"
      )
    ) +
    
    # Value labels
    geom_text(
      aes(
        label = sprintf(ifelse(Delta_pp > 0, "+%.1f pp", "%.1f pp"), Delta_pp),
        hjust = ifelse(Delta_pp > 0, -0.15, ifelse(Delta_pp < 0, 1.15, -0.15))
      ),
      fontface = "bold",
      size     = 5.2,
      color    = "black"
    ) +
    
    # -------------------------------------------------------------------------
  # 3. Prominent Annotations & Larger Arrows (Using Abbreviation)
  # -------------------------------------------------------------------------
  # Decrease annotation (left)
  annotate(
    "text",
    x        = -x_limit * 0.50,
    y        = n_features + 0.80,
    label    = "Decrease in rel. importance",
    fontface = "bold",
    size     = 5.0,
    color    = "grey20"
  ) +
    annotate(
      "segment",
      x         = -x_limit * 0.12,
      xend      = -x_limit * 0.88,
      y         = n_features + 1.25,
      yend      = n_features + 1.25,
      arrow     = arrow(type = "closed", length = unit(0.10, "inches")),
      linewidth = 0.9,
      color     = "grey30"
    ) +
    
    # Increase annotation (right)
    annotate(
      "text",
      x        = x_limit * 0.50,
      y        = n_features + 0.80,
      label    = "Increase in rel. importance",
      fontface = "bold",
      size     = 5.0,
      color    = "grey20"
    ) +
    annotate(
      "segment",
      x         = x_limit * 0.12,
      xend      = x_limit * 0.88,
      y         = n_features + 1.25,
      yend      = n_features + 1.25,
      arrow     = arrow(type = "closed", length = unit(0.10, "inches")),
      linewidth = 0.9,
      color     = "grey30"
    ) +
    
    # -------------------------------------------------------------------------
  # Axes & Expansion
  # -------------------------------------------------------------------------
  scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = seq(-x_limit, x_limit, by = 2)
  ) +
    scale_y_discrete(
      expand = expansion(mult = c(0.08, 0.28)) # Gives extra top space for arrows
    ) +
    coord_cartesian(clip = "off") +
    
    labs(
      title    = NULL,
      subtitle = NULL,
      x        = x_axis_label,
      y        = NULL
    ) +
    
    theme_large_nature() +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
      plot.margin  = margin(t = 20, r = 25, b = 15, l = 15)
    )
  
  # Save plot
  ggsave(
    file.path(save_path, paste0("Fig2e_", metric_name, "_", model_type, ".png")),
    p,
    width  = 8,
    height = 5,
    dpi    = 300,
    bg     = "white"
  )
  
  return(p)
}

# ==============================================================================
# --- Figure 2e Executions ---
# ==============================================================================

# Customer outage hours
run_fig2e(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2e(metric_name = "annual_coh", model_type = "RF")

# Cumulative outage hours
run_fig2e(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2e(metric_name = "annual_outage_hours", model_type = "RF")

# ==============================================================================
# --- 7. INDIVIDUAL MODULE: Figure 2f (Climate Envelope Expansion & SHAP) ---
# ==============================================================================
library(tidyverse)
library(xgboost)
library(ranger)
library(patchwork)
library(scales)

evaluate_fig2f_data <- function(clean_df, model_type = "XGBoost") {
  # 1. Detect Temperature Variable
  temp_candidates <- c("summer_max_T2C", "max_T2C", "avg_T2C")
  temp_var <- intersect(temp_candidates, names(clean_df))[1]
  if (is.na(temp_var)) stop("Error: No valid temperature column found in clean_df.")
  
  feature_cols <- setdiff(names(clean_df), c("county_group", "Year_Internal", "target", "fold"))
  
  df_proc <- clean_df %>%
    mutate(
      Period    = ifelse(Year_Internal <= 2018, "2015–2018 (Historical)", "2019–2022 (Recent)"),
      Period    = factor(Period, levels = c("2015–2018 (Historical)", "2019–2022 (Recent)")),
      Temp_Plot = .data[[temp_var]]
    )
  
  hist_df   <- df_proc %>% filter(Year_Internal <= 2018)
  recent_df <- df_proc %>% filter(Year_Internal > 2018)
  
  # 90th percentile threshold on historical period
  temp_thresh_upper <- quantile(hist_df$Temp_Plot, 0.90, na.rm = TRUE)
  pct_temp_exceed   <- mean(recent_df$Temp_Plot > temp_thresh_upper, na.rm = TRUE) * 100
  
  # --- Calculate SHAP Values ---
  X_mat <- as.matrix(df_proc[, feature_cols])
  
  if (model_type == "XGBoost") {
    dtrain <- xgb.DMatrix(data = X_mat, label = df_proc$target)
    params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
    mod    <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
    
    shap_raw <- predict(mod, dtrain, predcontrib = TRUE)
    shap_mat <- if (is.matrix(shap_raw)) shap_raw else matrix(shap_raw, ncol = length(feature_cols) + 1, byrow = TRUE)
    
    df_proc$SHAP_Temp <- if (temp_var %in% colnames(shap_mat)) {
      shap_mat[, temp_var]
    } else {
      shap_mat[, match(temp_var, feature_cols)]
    }
  } else {
    mod       <- ranger(target ~ ., data = df_proc[, c("target", feature_cols)], num.trees = 200)
    orig_data <- df_proc[, c("target", feature_cols)]
    base_data <- orig_data
    base_data[[temp_var]] <- mean(hist_df[[temp_var]], na.rm = TRUE)
    df_proc$SHAP_Temp     <- predict(mod, data = orig_data)$predictions - predict(mod, data = base_data)$predictions
  }
  
  df_proc <- df_proc %>%
    mutate(
      Envelope_Status = ifelse(Temp_Plot > temp_thresh_upper, "> Historical P90", "≤ Historical P90"),
      Envelope_Status = factor(Envelope_Status, levels = c("≤ Historical P90", "> Historical P90"))
    )
  
  list(
    data              = df_proc,
    temp_thresh_upper = temp_thresh_upper,
    pct_temp_exceed   = round(pct_temp_exceed),
    temp_name         = ifelse(grepl("max", temp_var), "Max Temperature (°C)", "Temperature (°C)")
  )
}

run_fig2f <- function(metric_name = "annual_coh", model_type = "XGBoost") {
  message(paste("[FIGURE 2f] Running Climate Envelope & SHAP for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data(metric_name)
  eval_res   <- evaluate_fig2f_data(clean_data, model_type)
  df_plot    <- eval_res$data
  
  # Exact color palette matching Figure 2
  period_colors <- c("2015–2018 (Historical)" = "#a0c4df", "2019–2022 (Recent)" = "#e8989a")
  colors_map    <- c("≤ Historical P90" = "#6ba2e6", "> Historical P90" = "#e05252")
  thresh_red    <- "#d9381e"
  
  # --- Unified Typography & Style Theme ---
  common_theme <- theme_classic() +
    theme(
      axis.title        = element_text(size = 14, face = "bold", color = "black"),
      axis.text         = element_text(size = 12, face = "bold", color = "black"),
      axis.line         = element_line(color = "black", linewidth = 0.8),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      legend.position   = "bottom",
      legend.direction  = "vertical",
      legend.justification = "left",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 12, face = "bold", color = "black"),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      plot.margin       = margin(t = 10, r = 15, b = 10, l = 10)
    )
  
  # --- Subplot 1: Histogram ---
  p_temp <- ggplot(df_plot, aes(x = Temp_Plot, fill = Period)) +
    geom_histogram(position = "identity", alpha = 0.65, bins = 35, color = NA) +
    geom_vline(xintercept = eval_res$temp_thresh_upper, linetype = "dashed", color = thresh_red, linewidth = 0.8) +
    scale_fill_manual(values = period_colors) +
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "k")) +
    coord_cartesian(clip = "off") +
    labs(x = "Temperature (°C)", y = "Count") +
    
    # Left Annotation (Matching size/font with right)
    annotate(
      "text",
      x          = eval_res$temp_thresh_upper + 0.6,
      y          = Inf,
      label      = paste0(eval_res$pct_temp_exceed, "% above\nHistorical P90"),
      fontface   = "bold",
      color      = thresh_red,
      size       = 4,
      vjust      = 1.3,
      hjust      = 0,
      lineheight = 0.95
    ) +
    common_theme +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
  
  # --- Subplot 2: SHAP Scatter ---
  p_shap <- ggplot(df_plot, aes(x = Temp_Plot, y = SHAP_Temp, color = Envelope_Status)) +
    geom_point(alpha = 0.7, size = 1.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = eval_res$temp_thresh_upper, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    scale_color_manual(values = colors_map) +
    coord_cartesian(clip = "off") +
    labs(x = "Temperature (°C)", y = "SHAP value") +
    
    # # Right Annotation (Matching size/font with left)
    # annotate(
    #   "text",
    #   x          = eval_res$temp_thresh_upper + 0.6,
    #   y          = max(df_plot$SHAP_Temp, na.rm = TRUE) * 0.95,
    #   label      = "Exceeds historical\nrange: SHAP values\nincrease",
    #   fontface   = "bold",
    #   color      = "black",
    #   size       = 4,
    #   vjust      = 1.0,
    #   hjust      = 0,
    #   lineheight = 0.95
    # ) +
    common_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  # --- Combine Subplots ---
  p_combined <- (p_temp + p_shap) +
    plot_layout(widths = c(1, 1))
  
  ggsave(
    file.path(save_path, paste0("Fig2f_Hist_", metric_name, "_", model_type, ".png")),
    p_combined,
    width  = 8.5,
    height = 5,
    dpi    = 300,
    bg     = "white"
  )
  
  return(p_combined)
}

# --- Execution Calls ---
run_fig2f(metric_name = "annual_coh", model_type = "XGBoost")
run_fig2f(metric_name = "annual_coh", model_type = "RF")
run_fig2f(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig2f(metric_name = "annual_outage_hours", model_type = "RF")
