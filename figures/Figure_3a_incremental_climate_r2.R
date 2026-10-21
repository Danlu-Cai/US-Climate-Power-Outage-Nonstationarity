# =============================================================================
# Script: Figure_3a_incremental_climate_r2.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description: Grouped county cross-validation of incremental climate predictive R-squared.
# Dependencies: dplyr, readr, ggplot2, tidyr, purrr, xgboost, ranger, stringr, conflicted
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

conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)

# ==============================================================================
# --- 1. Global Setup & Core Helper Functions ---
# ==============================================================================

Input_ROOT <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related" 
save_path  <- "//uahdata/rgroup/urbanrs/Danlu_data/2026-Digicare/0-eaglei_8hr_5000COH_Thresholds/NatureCC_related_final" 
csv_path   <- file.path(save_path, "Fig3_csv_results") 

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
if (!dir.exists(csv_path))  dir.create(csv_path, recursive = TRUE)

assign_us_region <- function(county_vec, state_vec = NULL) {
  fips_clean <- str_pad(as.character(county_vec), width = 5, side = "left", pad = "0")
  state_fips <- substr(fips_clean, 1, 2)
  
  ne_fips <- c("09", "10", "23", "24", "25", "33", "34", "36", "42", "44", "50", "11") 
  um_fips <- c("19", "26", "27", "55")                                             
  ov_fips <- c("17", "18", "21", "29", "39", "47", "54")                            
  se_fips <- c("01", "12", "13", "37", "45", "51")                                 
  nr_fips <- c("30", "31", "38", "46", "56")                                       
  so_fips <- c("05", "20", "22", "28", "40", "48")                                 
  sw_fips <- c("04", "08", "35", "49")                                             
  nw_fips <- c("16", "41", "53")                                                   
  we_fips <- c("06", "32")                                                         
  
  region <- case_when(
    state_fips %in% ne_fips ~ "Northeast",
    state_fips %in% um_fips ~ "Upper Midwest",
    state_fips %in% ov_fips ~ "Ohio Valley",
    state_fips %in% se_fips ~ "Southeast",
    state_fips %in% nr_fips ~ "Northern Rockies",
    state_fips %in% so_fips ~ "South",
    state_fips %in% sw_fips ~ "Southwest",
    state_fips %in% nw_fips ~ "Northwest",
    state_fips %in% we_fips ~ "West",
    TRUE ~ NA_character_
  )
  
  if (!is.null(state_vec)) {
    st_clean <- toupper(trimws(as.character(state_vec)))
    ne_abbr <- c("CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT", "DC")
    um_abbr <- c("IA", "MI", "MN", "WI")
    ov_abbr <- c("IL", "IN", "KY", "MO", "OH", "TN", "WV")
    se_abbr <- c("AL", "FL", "GA", "NC", "SC", "VA")
    nr_abbr <- c("MT", "NE", "ND", "SD", "WY")
    so_abbr <- c("AR", "KS", "LA", "MS", "OK", "TX")
    sw_abbr <- c("AZ", "CO", "NM", "UT")
    nw_abbr <- c("ID", "OR", "WA")
    we_abbr <- c("CA", "NV")
    
    region <- case_when(
      !is.na(region) ~ region,
      st_clean %in% ne_abbr ~ "Northeast",
      st_clean %in% um_abbr ~ "Upper Midwest",
      st_clean %in% ov_abbr ~ "Ohio Valley",
      st_clean %in% se_abbr ~ "Southeast",
      st_clean %in% nr_abbr ~ "Northern Rockies",
      st_clean %in% so_abbr ~ "South",
      st_clean %in% sw_abbr ~ "Southwest",
      st_clean %in% nw_abbr ~ "Northwest",
      st_clean %in% we_abbr ~ "West",
      TRUE ~ "South"
    )
  } else {
    region <- ifelse(is.na(region), "South", region)
  }
  
  return(region)
}

load_and_prep_data_fig3 <- function(metric_name) {
  file_path <- file.path(Input_ROOT, "Digicare_Merged_Yearly_2015_2022_CONUS.csv")
  if (!file.exists(file_path)) stop("Data file not found at: ", file_path)
  
  df_full <- read_csv(file_path, show_col_types = FALSE)
  
  pop_col     <- base::intersect(c("population", "Pop", "pop"), names(df_full))[1]
  density_col <- base::intersect(c("pop_density", "Pop_Density", "density"), names(df_full))[1]
  spatial_exposure_vars <- if (!is.na(density_col)) density_col else if (!is.na(pop_col)) pop_col else character(0)
  
  unweighted_degree_vars <- c("annual_CDD", "annual_HDD", spatial_exposure_vars)
  base_hazard_vars        <- c("avg_T2C", "annual_event_count", "annual_max_mag", "annual_heat_days", "annual_cold_days", "annual_storm_days")
  predictors              <- unique(c(base_hazard_vars, unweighted_degree_vars, "state"))
  vars_to_use             <- base::intersect(predictors, names(df_full))
  
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
  
  state_val <- if ("state" %in% names(df_transformed)) df_transformed$state else NULL
  df_transformed$Region <- assign_us_region(df_transformed[[county_col]], state_val)
  
  clean_df <- df_transformed %>%
    select(county_group = !!sym(county_col), Year_Internal, Region, target = !!sym(metric_name), all_of(vars_to_use)) %>%
    mutate(across(where(is.character) & !c(county_group, Region), ~as.numeric(as.factor(.)))) %>%
    mutate(across(-c(county_group, Region), ~as.numeric(as.character(.)))) %>%
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

theme_nature_dynamic <- function() {
  theme_classic(base_size = 16) +
    theme(
      panel.grid            = element_blank(),
      panel.background      = element_blank(),
      plot.background       = element_blank(),
      legend.background     = element_blank(),
      legend.box.background = element_blank(),
      legend.key            = element_blank(),
      plot.title.position   = "plot", 
      plot.title            = element_text(face = "bold", size = 18, color = "black", hjust = 0, margin = margin(b = 4)),
      plot.subtitle         = element_text(size = 14, color = "black", hjust = 0, margin = margin(b = 10)),
      axis.title.y          = element_text(face = "bold", size = 16, color = "black", margin = margin(r = 6)),
      axis.title.x          = element_blank(),
      axis.text.y           = element_text(face = "bold", size = 14, color = "black"),
      axis.text.x           = element_text(face = "bold", size = 12, color = "black", lineheight = 0.85),
      legend.title          = element_blank(),
      legend.text           = element_text(face = "bold", size = 14, color = "black"),
      axis.line             = element_line(color = "black", linewidth = 0.8),
      axis.ticks            = element_line(color = "black", linewidth = 0.8),
      plot.margin           = margin(t = 20, r = 24, b = 12, l = 16, unit = "pt")
    )
}

# ==============================================================================
# --- 2. Figure 3a Computation & Dynamic Plotting ---
# ==============================================================================

compute_incremental_r2 <- function(clean_df, model_type = "XGBoost", n_folds = 5) {
  all_features <- setdiff(names(clean_df), c("county_group", "Year_Internal", "Region", "target", "fold"))
  climate_vars <- c("avg_T2C", "annual_CDD", "annual_HDD", "annual_heat_days", 
                    "annual_cold_days", "annual_storm_days", "annual_event_count", "annual_max_mag")
  non_climate_features <- setdiff(all_features, climate_vars)
  
  if (length(non_climate_features) == 0) {
    non_climate_features <- intersect(c("pop_density", "population", "state"), names(clean_df))
  }
  
  periods <- list(
    "Historical (2015–2018)" = 2015:2018,
    "Recent (2019–2022)"     = 2019:2022
  )
  
  regions <- c(
    "All CONUS", "Northeast", "Upper Midwest", "Ohio Valley", "Southeast", 
    "Northern Rockies", "South", "Southwest", "Northwest", "West"
  )
  
  results <- list()
  
  calc_r2_cv <- function(sub_data, feat_cols) {
    if (nrow(sub_data) < 30 || length(unique(sub_data$county_group)) < n_folds) return(NA_real_)
    f_act <- numeric(0); f_pred <- numeric(0)
    for (k in 1:n_folds) {
      tr_sub <- sub_data %>% filter(fold != k)
      te_sub <- sub_data %>% filter(fold == k)
      if (nrow(te_sub) == 0) next
      p <- fit_predict(tr_sub, te_sub, feat_cols, model_type)
      f_act  <- c(f_act, te_sub$target)
      f_pred <- c(f_pred, p)
    }
    tss <- sum((f_act - mean(f_act))^2)
    if (tss == 0) return(0)
    # Retain raw cross-validated R² without artificial max(0, ...) lower bound
    1 - (sum((f_act - f_pred)^2) / tss)
  }
  
  for (p_name in names(periods)) {
    p_years <- periods[[p_name]]
    p_df    <- clean_df %>% filter(Year_Internal %in% p_years)
    
    for (reg in regions) {
      reg_df <- if (reg == "All CONUS") p_df else p_df %>% filter(Region == reg)
      
      r2_full        <- calc_r2_cv(reg_df, all_features)
      r2_non_climate <- calc_r2_cv(reg_df, non_climate_features)
      
      # Retain actual untruncated difference (r2_full - r2_non_climate)
      inc_r2 <- if (!is.na(r2_full) && !is.na(r2_non_climate)) (r2_full - r2_non_climate) else NA_real_
      
      results[[length(results) + 1]] <- data.frame(
        Period         = p_name,
        Region         = reg,
        Incremental_R2 = inc_r2
      )
    }
  }
  
  bind_rows(results)
}

run_fig3a <- function(metric_name = "annual_coh", model_type = "XGBoost") {
  message(paste("[FIGURE 3a] Running for Metric:", metric_name, "| Model:", model_type))
  clean_data <- load_and_prep_data_fig3(metric_name)
  
  res_df <- compute_incremental_r2(clean_data, model_type = model_type)
  
  res_df <- res_df %>%
    mutate(Incremental_R2 = round(Incremental_R2, 2))
  
  min_val <- min(res_df$Incremental_R2, na.rm = TRUE)
  max_val <- max(res_df$Incremental_R2, na.rm = TRUE)
  
  # Dynamic range handling supporting negative incremental R² values
  y_min  <- ifelse(min_val < 0, floor((min_val - 0.05) * 20) / 20, 0)
  y_max  <- ceiling((max_val + 0.08) * 20) / 20
  y_step <- 0.05
  
  delta_df <- res_df %>%
    pivot_wider(names_from = Period, values_from = Incremental_R2) %>%
    mutate(
      Hist_Val    = `Historical (2015–2018)`,
      Rec_Val     = `Recent (2019–2022)`,
      Delta_Val   = Rec_Val - Hist_Val,
      Delta_Label = sprintf("\u0394 = %+.2f", Delta_Val),
      Delta_Color = ifelse(Delta_Val >= 0, "#ED0000", "#00468B")
    )
  
  # FIX: Ensure Label_Y is always placed above y = 0 when both bars are negative
  delta_df <- delta_df %>%
    mutate(
      Label_Y = pmax(0, Hist_Val, Rec_Val) + ((y_max - y_min) * 0.075)
    ) %>%
    arrange(Region != "All CONUS", desc(Delta_Val))
  
  region_order <- delta_df$Region
  
  res_df <- res_df %>%
    mutate(
      Region = factor(Region, levels = region_order),
      Period = factor(Period, levels = c("Historical (2015–2018)", "Recent (2019–2022)"))
    )
  
  delta_df <- delta_df %>%
    mutate(Region = factor(Region, levels = region_order))
  
  write_csv(res_df, file.path(csv_path, paste0("Fig3a_Data_", metric_name, "_", model_type, ".csv")))
  
  format_region <- function(x) {
    x <- gsub("All CONUS", "All\nCONUS", x)
    x <- gsub("Upper Midwest", "Upper\nMidwest", x)
    x <- gsub("Ohio Valley", "Ohio\nValley", x)
    x <- gsub("Northern Rockies", "Northern\nRockies", x)
    return(x)
  }
  
  p <- ggplot(res_df, aes(x = Region, y = Incremental_R2, fill = Period)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    
    # Baseline line at 0 for clear zero threshold separation
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    
    # Bar value labels dynamically placed above/below bars depending on sign
    geom_text(
      aes(
        group = Period, 
        label = sprintf("%.2f", Incremental_R2),
        vjust = ifelse(Incremental_R2 >= 0, -0.4, 1.4)
      ),
      position = position_dodge(width = 0.72),
      fontface = "bold",
      size = 3.5,
      color = "black"
    ) +
    
    # Delta change labels placed clearly above the highest point / baseline
    geom_text(
      data = delta_df,
      aes(
        x = Region,
        y = Label_Y,
        label = Delta_Label,
        color = Delta_Color
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.0
    ) +
    
    scale_fill_manual(
      values = c(
        "Historical (2015–2018)" = "#00468B", 
        "Recent (2019–2022)"     = "#ED0000"
      )
    ) +
    scale_color_identity() +
    scale_x_discrete(labels = format_region) +
    scale_y_continuous(
      breaks = seq(y_min, y_max, by = y_step),
      labels = sprintf("%.2f", seq(y_min, y_max, by = y_step)),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    
    annotate(
      "text",
      x = 4.0,
      y = y_max - 0.02,
      label = "\u0394 = Change (Recent \u2013 Historical)",
      color = "#D92B1D",
      fontface = "bold",
      size = 4.6,
      hjust = 0
    ) +
    
    labs(
      title    = NULL,
      subtitle = NULL,
      x        = NULL,
      y        = expression(bold("Incremental cross-validated R"^2 ~ "(full − non-climate)"))
    ) +
    theme_nature_dynamic() +
    theme(
      axis.title.y     = element_text(face = "bold", size = 16, color = "black", margin = margin(r = 8)),
      axis.text.y      = element_text(face = "bold", size = 15, color = "black"),
      axis.text.x      = element_text(face = "bold", size = 13, color = "black", lineheight = 0.85),
      legend.text      = element_text(face = "bold", size = 15, color = "black"),
      
      legend.position  = c(0.01, 0.96),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.key.size  = unit(0.55, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    ) +
    guides(fill = guide_legend(ncol = 1, title = NULL))
  
  ggsave(
    file.path(save_path, paste0("Fig3a_", metric_name, "_", model_type, ".png")),
    p,
    width  = 11.0,
    height = 7.5,
    dpi    = 300,
    bg     = "white"
  )
  
  return(p)
}

# --- Figure 3a Executions ---
run_fig3a(metric_name = "annual_coh", model_type = "XGBoost")
run_fig3a(metric_name = "annual_coh", model_type = "RF")
run_fig3a(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig3a(metric_name = "annual_outage_hours", model_type = "RF")
run_fig3a(metric_name = "annual_coh", model_type = "RF")
run_fig3a(metric_name = "annual_outage_hours", model_type = "XGBoost")
run_fig3a(metric_name = "annual_outage_hours", model_type = "RF")
